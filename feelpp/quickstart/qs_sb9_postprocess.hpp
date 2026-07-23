/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

/**
 * \file qs_sb9_postprocess.hpp
 * \brief SB9 shell strain/stress postprocessing utilities for the quickstart applications.
 *
 * This header gathers the element-local SB9 evaluation code used by
 * `qs_sb9.cpp` to build semantic symmetric tensor fields for strain and
 * stress.  The output tensor fields can be exported normally or printed with
 * \c printMatlab() in the component order used by the MATLAB validation data:
 * `xx, yy, zz, xy, xz, yz`.
 *
 * \par Example
 * \code{.cpp}
 * auto tensorSpace = Pdhms<1>( mesh );
 * auto epsilon = tensorSpace->element( "epsilon" );
 * auto sigma = tensorSpace->element( "sigma" );
 *
 * Feel::Quickstart::SB9::fillSymmetricFields(
 *     mesh, uh, alphah, epsilon, sigma,
 *     lambda, mu, alphaScale, pinchingBpzScale, shellShearFactor );
 *
 * auto ctx = Feel::Quickstart::SB9::matlabFieldContext(
 *     tensorSpace, caseConfig.referenceChecks, fallbackPoint );
 * epsilon.printMatlab( "epsilon_matlab", ctx,
 *                      Feel::Quickstart::SB9::matlabEpsilonFormat(),
 *                      true, "epsilon" );
 * sigma.printMatlab( "sigma_matlab", ctx,
 *                    Feel::Quickstart::SB9::matlabSigmaFormat(),
 *                    true, "sigma" );
 * \endcode
 */
#ifndef FEELPP_QUICKSTART_QS_SB9_POSTPROCESS_HPP
#define FEELPP_QUICKSTART_QS_SB9_POSTPROCESS_HPP 1

#include <feel/feeldiscr/tensorformat.hpp>
#include <feel/feelvf/sb9_bending.hpp>
#include <feel/feelvf/sb9_pinching.hpp>
#include <feel/feelvf/sb9_shear.hpp>
#include <feel/feelvf/vf.hpp>

#include "qs_elasticity_checks.hpp"

#include <Eigen/Core>

#include <numbers>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

/**
 * \namespace Feel::Quickstart::SB9
 * \brief Quickstart helpers for SB9 postprocessing and validation output.
 *
 * The functions in this namespace are intentionally kept outside the core
 * variational-expression headers.  They are diagnostic/application utilities:
 * they convert solved SB9 fields into discontinuous symmetric tensor fields,
 * print those fields in MATLAB-compatible order, and compare them with JSON
 * validation references.
 */
namespace Feel::Quickstart::SB9
{
namespace qsec = Feel::Quickstart::ElasticityChecks;

/**
 * \brief Return the MATLAB validation format for SB9 strain vectors.
 *
 * The validation files store strain as `xx, yy, zz, xy, xz, yz` and use
 * engineering shear entries for the off-diagonal components.  This is the
 * format expected by \c FunctionSpace::Element::printMatlab() when printing
 * `epsilon`.
 *
 * \return Symmetric tensor format with diagonal-first order and engineering
 *         shear scaling.
 */
inline SymmetricTensorFormat
matlabEpsilonFormat()
{
    return { SymmetricTensorOrder::DiagonalFirst,
             SymmetricTensorScaling::EngineeringShear };
}

/**
 * \brief Return the MATLAB validation format for SB9 stress vectors.
 *
 * Stress uses the same diagonal-first component order as strain but keeps
 * tensor shear values, i.e. off-diagonal components are not doubled.
 *
 * \return Symmetric tensor format with diagonal-first order and tensor shear
 *         scaling.
 */
inline SymmetricTensorFormat
matlabSigmaFormat()
{
    return { SymmetricTensorOrder::DiagonalFirst,
             SymmetricTensorScaling::Tensor };
}

/**
 * \brief Convert an Eigen 3D point to a Feel++ node.
 *
 * \param p Three-dimensional physical point.
 * \return Feel++ node containing the same coordinates.
 */
inline node_type
toNode( qsec::point_type<3> const& p )
{
    node_type n( 3 );
    n( 0 ) = p( 0 );
    n( 1 ) = p( 1 );
    n( 2 ) = p( 2 );
    return n;
}

/**
 * \brief Materialize one SB9 coefficient cache as a 6-by-24 element matrix.
 *
 * The SB9 variational expressions operate through basis proxies.  For
 * postprocessing we instead need explicit local matrices so that the solved
 * displacement vector can be multiplied directly at a physical evaluation
 * point.
 *
 * \tparam Kind Compile-time coefficient family selector, such as
 *         \c SB9BendingKind::Bm0, \c SB9PinchingKind::Bpc, or
 *         \c SB9ShearKind::Bc0.
 * \tparam CacheType SB9 kernel cache type providing
 *         \c fillVectorCoefficients().
 * \param cache Element-local SB9 coefficient cache.
 * \return Matrix with six symmetric-storage rows and 24 displacement columns
 *         ordered by component block then node.
 */
template <auto Kind, typename CacheType>
Eigen::Matrix<double, 6, 24>
coefficientMatrix( CacheType const& cache )
{
    Eigen::Matrix<double, 6, 24> matrix;
    matrix.setZero();

    for ( uint16_type component = 0; component < 3; ++component )
    {
        for ( uint16_type node = 0; node < 8; ++node )
        {
            Eigen::Matrix<double, 6, 1> coeff;
            coeff.setZero();
            cache.template fillVectorCoefficients<Kind>( coeff, node, component );
            matrix.col( component*8 + node ) = coeff;
        }
    }

    return matrix;
}

/**
 * \brief Evaluate a symmetric tensor field at one point in a requested format.
 *
 * This is the postprocessed-field counterpart of the validation checker:
 * it samples a \c Pdhms field and returns the six components in the supplied
 * \p format.
 *
 * \tparam SpaceType Symmetric tensor function-space pointer type.
 * \tparam FieldElementType Element type of the symmetric tensor field.
 * \param Sh Function space used to build the evaluation context.
 * \param field Symmetric tensor element, typically `epsilon` or `sigma`.
 * \param point Physical evaluation point.
 * \param format Component order and scaling requested by the caller.
 * \return Six-component vector in \p format.
 */
template <typename SpaceType, typename FieldElementType>
qsec::voigt_type
evaluateSymmetricFieldReference( SpaceType const& Sh,
                                 FieldElementType const& field,
                                 qsec::point_type<3> const& point,
                                 SymmetricTensorFormat const& format )
{
    auto ctx = Sh->context();
    ctx.add( toNode( point ) );
    auto values = field.evaluateSymmetric( ctx, format, true );
    if ( values.rows() < 1 || values.cols() != 6 )
        throw std::runtime_error( "SB9 tensor field evaluation expected one row with six components" );

    qsec::voigt_type result;
    for ( int c = 0; c < 6; ++c )
        result( c ) = values( 0, c );
    return result;
}

/**
 * \brief Physical symmetric tensor entries in diagonal-first notation.
 *
 * The helper stores true tensor entries.  When values come from an engineering
 * shear vector, the off-diagonal components are converted back to tensor
 * entries before assignment to the semantic tensor field.
 */
struct SymmetricTensorComponents
{
    /// Normal xx component.
    double xx = 0.0;
    /// Normal yy component.
    double yy = 0.0;
    /// Normal zz component.
    double zz = 0.0;
    /// Tensor xy component.
    double xy = 0.0;
    /// Tensor xz component.
    double xz = 0.0;
    /// Tensor yz component.
    double yz = 0.0;
};

/**
 * \brief Convert a validation-order vector to physical tensor entries.
 *
 * \param values Six values ordered as `xx, yy, zz, xy, xz, yz`.
 * \param engineeringShear Whether the shear entries are engineering shear
 *        values and must therefore be halved before tensor-field assignment.
 * \return Physical symmetric tensor entries.
 */
inline SymmetricTensorComponents
toSymmetricTensorComponents( qsec::voigt_type const& values, bool engineeringShear )
{
    double const shearScale = engineeringShear ? 0.5 : 1.0;
    return { values( 0 ),
             values( 1 ),
             values( 2 ),
             shearScale * values( 3 ),
             shearScale * values( 4 ),
             shearScale * values( 5 ) };
}

/**
 * \brief Assign one constant symmetric tensor value on one mesh element.
 *
 * The assignment uses the semantic tensor-component API so that symmetric
 * aliases such as `(X,Y)` and `(Y,X)` stay consistent.
 *
 * \tparam FieldElementType Symmetric tensor field element type.
 * \tparam MeshType Mesh type.
 * \param field Tensor field to modify.
 * \param mesh Mesh owning \p elementId.
 * \param elementId Element identifier receiving the constant value.
 * \param values Validation-order values `xx, yy, zz, xy, xz, yz`.
 * \param engineeringShear Whether \p values stores engineering shear entries.
 */
template <typename FieldElementType, typename MeshType>
void
assignSymmetricTensorOnElement( FieldElementType& field,
                                std::shared_ptr<MeshType> const& mesh,
                                index_type elementId,
                                qsec::voigt_type const& values,
                                bool engineeringShear )
{
    auto const tensor = toSymmetricTensorComponents( values, engineeringShear );
    field.tensorComponent( Component::X, Component::X ).on( _range=idedelements( mesh, elementId ), _expr=cst( tensor.xx ) );
    field.tensorComponent( Component::Y, Component::Y ).on( _range=idedelements( mesh, elementId ), _expr=cst( tensor.yy ) );
    field.tensorComponent( Component::Z, Component::Z ).on( _range=idedelements( mesh, elementId ), _expr=cst( tensor.zz ) );
    field.tensorComponent( Component::X, Component::Y ).on( _range=idedelements( mesh, elementId ), _expr=cst( tensor.xy ) );
    field.tensorComponent( Component::X, Component::Z ).on( _range=idedelements( mesh, elementId ), _expr=cst( tensor.xz ) );
    field.tensorComponent( Component::Y, Component::Z ).on( _range=idedelements( mesh, elementId ), _expr=cst( tensor.yz ) );
}

/**
 * \brief Evaluate SB9 strain and stress at a geometric mapping context.
 *
 * The function reconstructs the element-local SB9 operators at the supplied
 * context, multiplies them by the solved displacement and internal scalar
 * dofs, and returns the validation-order strain/stress vectors.
 *
 * The strain vector uses engineering shear components.  The stress vector uses
 * tensor shear components and isotropic 3D Hooke law with Lamé coefficients
 * \p lambda and \p mu.
 *
 * \tparam GmcPtrType Geometric mapping context pointer type.
 * \tparam DisplacementElementType Displacement field element type.
 * \tparam AlphaElementType Internal SB9 scalar field element type.
 * \param gmc Geometric context at the physical point/reference coordinate.
 * \param uh Solved displacement field.
 * \param alphah Solved internal SB9 scalar field.
 * \param lambda First Lamé coefficient.
 * \param mu Shear modulus.
 * \param alphaScale Scale applied to the internal scalar correction.
 * \param pinchingBpzScale Scale applied to the `zeta*Bpz` pinching term.
 * \param shellShearFactor Coefficient in the transverse shear shape function.
 * \return Pair `(epsilon, sigma)` in MATLAB validation order.
 */
template <typename GmcPtrType, typename DisplacementElementType, typename AlphaElementType>
std::pair<qsec::voigt_type, qsec::voigt_type>
evaluateAtGmc( GmcPtrType const& gmc,
               DisplacementElementType const& uh,
               AlphaElementType const& alphah,
               double lambda,
               double mu,
               double alphaScale,
               double pinchingBpzScale,
               double shellShearFactor )
{
    auto const data = Feel::vf::detail::computeShellCellGeometry( gmc.get() );
    using geometry_data_type = std::decay_t<decltype( data )>;
    Feel::vf::detail::SB9BendingKernelCache<geometry_data_type> const bendingCache( data );
    Feel::vf::detail::SB9PinchingKernelCache<geometry_data_type> const pinchingCache( data );
    Feel::vf::detail::SB9ShearKernelCache<geometry_data_type> const shearCache( data );

    Eigen::Matrix<double, 24, 1> uLocal;
    uLocal.setZero();
    index_type const elementId = gmc->element().id();
    for ( uint16_type component = 0; component < 3; ++component )
        for ( uint16_type node = 0; node < 8; ++node )
            uLocal( component*8 + node ) = uh.localToGlobal( elementId, node, component );

    auto const bm0 = coefficientMatrix<Feel::vf::detail::SB9BendingKind::Bm0>( bendingCache );
    auto const bb0 = coefficientMatrix<Feel::vf::detail::SB9BendingKind::Bb0>( bendingCache );
    auto const bpc = coefficientMatrix<Feel::vf::detail::SB9PinchingKind::Bpc>( pinchingCache );
    auto const bpz = coefficientMatrix<Feel::vf::detail::SB9PinchingKind::Bpz>( pinchingCache );
    auto const bc0 = coefficientMatrix<Feel::vf::detail::SB9ShearKind::Bc0>( shearCache );

    double const zetaValue = gmc->xRefs()( 2, 0 );
    double const shearWeight = shellShearFactor * ( 1.0 - zetaValue*zetaValue );
    auto const membrane = ( bm0 + zetaValue*bb0 ) * uLocal;
    auto const pinching = ( bpc + pinchingBpzScale*zetaValue*bpz ) * uLocal;
    auto const shear = bc0 * uLocal;

    Eigen::Matrix<double, 6, 1> epsilonMandel;
    epsilonMandel.setZero();
    epsilonMandel( 0 ) = membrane( 0 );
    epsilonMandel( 1 ) = std::numbers::sqrt2_v<double> * membrane( 3 );
    epsilonMandel( 2 ) = std::numbers::sqrt2_v<double> * shearWeight * shear( 4 );
    epsilonMandel( 3 ) = membrane( 1 );
    epsilonMandel( 4 ) = std::numbers::sqrt2_v<double> * shearWeight * shear( 5 );
    epsilonMandel( 5 ) = pinching( 2 ) +
                         alphaScale * ( -4.0*zetaValue/data.thickness ) * alphah.localToGlobal( elementId, 0, 0 );

    qsec::voigt_type epsilon;
    epsilon << epsilonMandel( 0 ),
               epsilonMandel( 3 ),
               epsilonMandel( 5 ),
               std::numbers::sqrt2_v<double> * epsilonMandel( 1 ),
               std::numbers::sqrt2_v<double> * epsilonMandel( 2 ),
               std::numbers::sqrt2_v<double> * epsilonMandel( 4 );

    qsec::voigt_type sigma;
    double const traceEpsilon = epsilon( 0 ) + epsilon( 1 ) + epsilon( 2 );
    sigma << lambda*traceEpsilon + 2.0*mu*epsilon( 0 ),
             lambda*traceEpsilon + 2.0*mu*epsilon( 1 ),
             lambda*traceEpsilon + 2.0*mu*epsilon( 2 ),
             mu*epsilon( 3 ),
             mu*epsilon( 4 ),
             mu*epsilon( 5 );

    return { epsilon, sigma };
}

/**
 * \brief Fill Pdhms tensor fields with SB9 center-point strain and stress.
 *
 * Each element receives the SB9 value evaluated at its reference center.  This
 * gives a compact discontinuous tensor field that can be exported and sampled
 * with the same semantic tensor-field API used by the validation checks.
 *
 * \tparam MeshPtrType Mesh shared-pointer type.
 * \tparam DisplacementElementType Displacement field element type.
 * \tparam AlphaElementType Internal SB9 scalar field element type.
 * \tparam EpsilonElementType Symmetric tensor field element type for strain.
 * \tparam SigmaElementType Symmetric tensor field element type for stress.
 * \param mesh Shell mesh.
 * \param uh Solved displacement field.
 * \param alphah Solved internal SB9 scalar field.
 * \param epsilonh Output strain field. Values are assigned as tensor entries;
 *        MATLAB printing should use \ref matlabEpsilonFormat().
 * \param sigmah Output stress field. Values are assigned as tensor entries;
 *        MATLAB printing should use \ref matlabSigmaFormat().
 * \param lambda First Lamé coefficient.
 * \param mu Shear modulus.
 * \param alphaScale Scale applied to the internal scalar correction.
 * \param pinchingBpzScale Scale applied to the `zeta*Bpz` pinching term.
 * \param shellShearFactor Coefficient in the transverse shear shape function.
 */
template <typename MeshPtrType, typename DisplacementElementType, typename AlphaElementType,
          typename EpsilonElementType, typename SigmaElementType>
void
fillSymmetricFields( MeshPtrType const& mesh,
                     DisplacementElementType const& uh,
                     AlphaElementType const& alphah,
                     EpsilonElementType& epsilonh,
                     SigmaElementType& sigmah,
                     double lambda,
                     double mu,
                     double alphaScale,
                     double pinchingBpzScale,
                     double shellShearFactor )
{
    using mesh_type = std::remove_cvref_t<decltype( *mesh )>;
    typename mesh_type::gm_type::matrix_node_t_type referenceCenter( mesh_type::nDim, 1 );
    referenceCenter.clear();
    auto gmpc = mesh->gm()->preCompute( mesh->gm(), referenceCenter );

    for ( auto const& elementRef : elements( mesh ) )
    {
        auto const& element = boost::unwrap_ref( elementRef );
        auto gmc = mesh->gm()->template context<vm::POINT|vm::JACOBIAN|vm::KB>( element, gmpc );
        auto const [epsilonValue, sigmaValue] =
            evaluateAtGmc( gmc, uh, alphah, lambda, mu, alphaScale, pinchingBpzScale, shellShearFactor );

        assignSymmetricTensorOnElement( epsilonh, mesh, element.id(), epsilonValue, true );
        assignSymmetricTensorOnElement( sigmah, mesh, element.id(), sigmaValue, false );
    }
}

/**
 * \brief Build the evaluation context used by MATLAB tensor diagnostics.
 *
 * If the case contains field-reference points, all of them are added in JSON
 * order.  Otherwise the supplied fallback point is used.  The resulting
 * context is suitable for `epsilon.printMatlab()` and `sigma.printMatlab()`.
 *
 * \tparam SpaceType Symmetric tensor function-space pointer type.
 * \param Sh Function space used to create the context.
 * \param config Parsed elasticity reference-check configuration.
 * \param fallbackPoint Point used when no field references are configured.
 * \return A Feel++ function-space context containing the diagnostic points.
 */
template <typename SpaceType>
auto
matlabFieldContext( SpaceType const& Sh,
                    qsec::Config<3> const& config,
                    qsec::point_type<3> const& fallbackPoint )
{
    auto ctx = Sh->context();
    if ( config.fieldReferences.empty() )
    {
        ctx.add( toNode( fallbackPoint ) );
        return ctx;
    }

    for ( auto const& fieldReference : config.fieldReferences )
        ctx.add( toNode( fieldReference.point ) );
    return ctx;
}

/**
 * \brief Test whether the reference configuration contains a named target.
 *
 * Both tensor field references and cantilever references are inspected because
 * SB9 validation cases may use either reference family.
 *
 * \param config Parsed reference-check configuration.
 * \param target Target name, for example `sb9g25`.
 * \return true if at least one reference item has \p target.
 */
inline bool
hasTarget( qsec::Config<3> const& config, std::string const& target )
{
    for ( auto const& fieldReference : config.fieldReferences )
        if ( fieldReference.targets.find( target ) != fieldReference.targets.end() )
            return true;
    for ( auto const& cantilever : config.cantilevers )
        if ( cantilever.targets.find( target ) != cantilever.targets.end() )
            return true;
    return false;
}

/**
 * \brief Select the default validation target for the SB9 quickstart.
 *
 * Shared shell JSON files often default to `hexa8` because they are also used
 * by the standard elasticity quickstart.  For the SB9 quickstart this helper
 * prefers `sb9g25`, then `sb9`, unless the user explicitly set
 * `--checks.target`.
 *
 * \param config Parsed reference-check configuration.
 * \param explicitTarget Whether `--checks.target` was explicitly provided.
 * \return Target name to use for reference checks.
 */
inline std::string
preferredTarget( qsec::Config<3> const& config, bool explicitTarget )
{
    if ( explicitTarget && !config.target.empty() )
        return config.target;
    if ( hasTarget( config, "sb9g25" ) )
        return "sb9g25";
    if ( hasTarget( config, "sb9" ) )
        return "sb9";
    return config.target;
}

/**
 * \brief Compare postprocessed SB9 tensor fields against JSON references.
 *
 * The comparison samples the supplied `epsilon` and `sigma` fields at each
 * configured field-reference point, using the MATLAB validation formats from
 * \ref matlabEpsilonFormat() and \ref matlabSigmaFormat().
 *
 * \tparam SpaceType Symmetric tensor function-space pointer type.
 * \tparam EpsilonElementType Strain field element type.
 * \tparam SigmaElementType Stress field element type.
 * \param config Parsed reference-check configuration.
 * \param Sh Symmetric tensor function space.
 * \param epsilonh Postprocessed SB9 strain field.
 * \param sigmah Postprocessed SB9 stress field.
 * \param target Reference target name.
 * \param checkReference If true, throw when configured tolerances fail.  If
 *        false, comparisons are printed in report-only mode.
 * \param os Output stream for the comparison report.
 * \return zero on success.
 */
template <typename SpaceType, typename EpsilonElementType, typename SigmaElementType>
int
checkFieldReferences( qsec::Config<3> const& config,
                      SpaceType const& Sh,
                      EpsilonElementType const& epsilonh,
                      SigmaElementType const& sigmah,
                      std::string const& target,
                      bool checkReference,
                      std::ostream& os = std::cout )
{
    for ( auto const& fieldReference : config.fieldReferences )
    {
        auto const targetIt = fieldReference.targets.find( target );
        if ( targetIt == fieldReference.targets.end() )
            throw std::invalid_argument( "field reference '" + fieldReference.name + "' has no target '" + target + "'" );

        auto const computedEpsilon = evaluateSymmetricFieldReference( Sh, epsilonh, fieldReference.point, matlabEpsilonFormat() );
        auto const computedSigma = evaluateSymmetricFieldReference( Sh, sigmah, fieldReference.point, matlabSigmaFormat() );

        os << "SB9 field reference check '" << fieldReference.name << "' target '" << target << "'\n";
        if ( !fieldReference.description.empty() )
            os << "SB9 field reference description = " << fieldReference.description << "\n";

        bool ok = true;
        auto const& reference = targetIt->second;
        bool const enabled = checkReference && fieldReference.check;
        if ( reference.hasEpsilon )
            ok = qsec::reportVoigtComparison( "epsilon", computedEpsilon, reference.epsilon,
                                              fieldReference.tolerance, enabled, os ) && ok;
        if ( reference.hasSigma )
            ok = qsec::reportVoigtComparison( "sigma", computedSigma, reference.sigma,
                                              fieldReference.tolerance, enabled, os ) && ok;
        if ( !ok )
            throw std::runtime_error( "SB9 field reference '" + fieldReference.name + "' check failed" );
    }
    return 0;
}

} // namespace Feel::Quickstart::SB9

#endif
