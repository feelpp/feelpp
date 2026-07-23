/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
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

namespace Feel::Quickstart::SB9
{
namespace qsec = Feel::Quickstart::ElasticityChecks;

inline SymmetricTensorFormat
matlabEpsilonFormat()
{
    return { SymmetricTensorOrder::DiagonalFirst,
             SymmetricTensorScaling::EngineeringShear };
}

inline SymmetricTensorFormat
matlabSigmaFormat()
{
    return { SymmetricTensorOrder::DiagonalFirst,
             SymmetricTensorScaling::Tensor };
}

inline node_type
toNode( qsec::point_type<3> const& p )
{
    node_type n( 3 );
    n( 0 ) = p( 0 );
    n( 1 ) = p( 1 );
    n( 2 ) = p( 2 );
    return n;
}

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

struct SymmetricTensorComponents
{
    double xx = 0.0;
    double yy = 0.0;
    double zz = 0.0;
    double xy = 0.0;
    double xz = 0.0;
    double yz = 0.0;
};

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
