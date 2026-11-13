/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file piolatransform.hpp
 * @brief Pointwise covariant and contravariant Piola transformations.
 */
#if !defined( FEELPP_FEELPOLY_PIOLATRANSFORM_HPP )
#define FEELPP_FEELPOLY_PIOLATRANSFORM_HPP 1

#include <Eigen/Core>
#include <Eigen/QR>

#include <array>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace Feel::geometry
{
/**
 * @brief Return the orthogonal projector onto a tangent space.
 *
 * For a nonzero normal vector `n`, computes
 * `P_T = I - n n^T / (n.n)`. The normalization is performed internally so
 * callers may provide either a unit normal or a Piola-scaled normal.
 *
 * @param normal normal or conormal defining the tangent space
 * @return symmetric idempotent tangent projector
 * @throws std::invalid_argument if the normal is numerically zero
 */
template <typename VectorDerived>
[[nodiscard]] auto tangentProjector(
    Eigen::MatrixBase<VectorDerived> const& normal )
{
    using value_type = typename VectorDerived::Scalar;
    constexpr int rows = VectorDerived::RowsAtCompileTime;
    Eigen::Matrix<value_type,rows,rows> projector(
        normal.rows(), normal.rows() );
    value_type const squaredNorm = normal.squaredNorm();
    if ( squaredNorm <= Eigen::NumTraits<value_type>::dummy_precision() )
        throw std::invalid_argument(
            "A tangent projector requires a nonzero normal" );
    projector.setIdentity();
    projector.noalias() -= normal*normal.transpose()/squaredNorm;
    return projector;
}

/**
 * @brief Return the orthogonal projector onto a column-vector subspace.
 * @param orthonormalFrame orthonormal basis vectors stored by columns
 * @return `Q Q^T`
 */
template <typename MatrixDerived>
[[nodiscard]] auto subspaceProjector(
    Eigen::MatrixBase<MatrixDerived> const& orthonormalFrame )
{
    return ( orthonormalFrame*orthonormalFrame.transpose() ).eval();
}

/**
 * @brief Construct an oriented orthonormal frame with a prescribed normal.
 *
 * The first candidate is the supplied normal; subsequent candidates span the
 * tangent space. A Householder QR factorization provides a stable frame. Each
 * column sign follows its corresponding candidate after orthogonalization,
 * and a square frame is made direct by adjusting its last tangent.
 *
 * @param normal prescribed nonzero normal or conormal
 * @param tangentCandidates independent tangent candidates stored by columns
 * @return orthonormal frame whose first column follows `normal`
 */
template <typename VectorDerived, typename TangentDerived>
[[nodiscard]] auto orientedOrthonormalFrame(
    Eigen::MatrixBase<VectorDerived> const& normal,
    Eigen::MatrixBase<TangentDerived> const& tangentCandidates )
{
    using value_type = typename VectorDerived::Scalar;
    constexpr int rows = VectorDerived::RowsAtCompileTime;
    constexpr int tangentColumns = TangentDerived::ColsAtCompileTime;
    constexpr int columns = tangentColumns == Eigen::Dynamic
                                ? Eigen::Dynamic
                                : tangentColumns+1;
    Eigen::Matrix<value_type,rows,columns> candidates(
        normal.rows(), tangentCandidates.cols()+1 );
    candidates.col( 0 ) = normal;
    candidates.rightCols( tangentCandidates.cols() ) = tangentCandidates;

    Eigen::HouseholderQR<decltype( candidates )> factorization( candidates );
    Eigen::Matrix<value_type,rows,columns> frame(
        candidates.rows(), candidates.cols() );
    frame = factorization.householderQ() *
        Eigen::Matrix<value_type,rows,columns>::Identity(
            candidates.rows(), candidates.cols() );

    for ( int column = 0; column < frame.cols(); ++column )
        if ( frame.col( column ).dot( candidates.col( column ) ) < value_type{} )
            frame.col( column ) *= value_type( -1 );

    if constexpr ( rows == Eigen::Dynamic || columns == Eigen::Dynamic )
    {
        if ( frame.rows() == frame.cols() && frame.rows() > 1 &&
             frame.determinant() < value_type{} )
            frame.col( frame.cols()-1 ) *= value_type( -1 );
    }
    else if constexpr ( rows == columns && rows > 1 )
    {
        if ( frame.determinant() < value_type{} )
            frame.col( frame.cols()-1 ) *= value_type( -1 );
    }
    return frame;
}

//!
//! \brief Push an H(curl) vector from reference to physical coordinates.
//!
//! Computes `K^{-T} vHat`, where `inverseJacobianTranspose` is `K^{-T}` for
//! full-dimensional cells and the corresponding metric pseudoinverse transpose
//! for embedded cells.
//!
//! \param inverseJacobianTranspose pointwise inverse Jacobian transpose
//! \param referenceVector vector in reference coordinates
//! \return physical covariant-Piola vector
template <typename MatrixDerived, typename VectorDerived>
[[nodiscard]] auto covariantPiolaPushForward(
    Eigen::MatrixBase<MatrixDerived> const& inverseJacobianTranspose,
    Eigen::MatrixBase<VectorDerived> const& referenceVector )
{
    return ( inverseJacobianTranspose * referenceVector ).eval();
}

//!
//! \brief Pull an H(curl) vector from physical to reference coordinates.
//!
//! Computes `K^T v`, the inverse of the covariant Piola push-forward.
//!
//! \param jacobian pointwise physical Jacobian `K`
//! \param physicalVector vector in physical coordinates
//! \return vector in reference coordinates
template <typename MatrixDerived, typename VectorDerived>
[[nodiscard]] auto covariantPiolaPullback(
    Eigen::MatrixBase<MatrixDerived> const& jacobian,
    Eigen::MatrixBase<VectorDerived> const& physicalVector )
{
    return ( jacobian.transpose() * physicalVector ).eval();
}

//!
//! \brief Push an H(div) vector from reference to physical coordinates.
//!
//! Computes `K vHat / J`, with `J` the pointwise integration Jacobian.
//!
//! \param jacobian pointwise physical Jacobian `K`
//! \param integrationJacobian pointwise positive integration Jacobian
//! \param referenceVector vector in reference coordinates
//! \return physical contravariant-Piola vector
template <typename MatrixDerived, typename Scalar, typename VectorDerived>
[[nodiscard]] auto contravariantPiolaPushForward(
    Eigen::MatrixBase<MatrixDerived> const& jacobian,
    Scalar const& integrationJacobian,
    Eigen::MatrixBase<VectorDerived> const& referenceVector )
{
    return ( jacobian * referenceVector / integrationJacobian ).eval();
}

//!
//! \brief Pull an H(div) vector from physical to reference coordinates.
//!
//! Computes `J K^{-1} v = J (K^{-T})^T v`.
//!
//! \param inverseJacobianTranspose pointwise inverse Jacobian transpose
//! \param integrationJacobian pointwise positive integration Jacobian
//! \param physicalVector vector in physical coordinates
//! \return vector in reference coordinates
template <typename MatrixDerived, typename Scalar, typename VectorDerived>
[[nodiscard]] auto contravariantPiolaPullback(
    Eigen::MatrixBase<MatrixDerived> const& inverseJacobianTranspose,
    Scalar const& integrationJacobian,
    Eigen::MatrixBase<VectorDerived> const& physicalVector )
{
    return ( integrationJacobian * inverseJacobianTranspose.transpose() *
             physicalVector )
        .eval();
}

//!
//! \brief Transform a reference tangent used by an H(curl) functional.
//!
//! The dual test direction is `K tHat`; pairing it with a physical covariant
//! field reproduces `vHat . tHat` without normalizing the tangent.
template <typename MatrixDerived, typename VectorDerived>
[[nodiscard]] auto covariantPiolaTestDirection(
    Eigen::MatrixBase<MatrixDerived> const& jacobian,
    Eigen::MatrixBase<VectorDerived> const& referenceTangent )
{
    return ( jacobian * referenceTangent ).eval();
}

//!
//! \brief Transform a reference normal used by an H(div) functional.
//!
//! Computes `J K^{-T} nHat` and optionally includes a reference-facet measure.
//! It is the pointwise cofactor direction preserving normal fluxes.
template <typename MatrixDerived, typename Scalar, typename VectorDerived>
[[nodiscard]] auto contravariantPiolaTestDirection(
    Eigen::MatrixBase<MatrixDerived> const& inverseJacobianTranspose,
    Scalar const& integrationJacobian,
    Eigen::MatrixBase<VectorDerived> const& referenceNormal,
    Scalar const& referenceFacetMeasure = Scalar( 1 ) )
{
    return ( integrationJacobian * referenceFacetMeasure *
             inverseJacobianTranspose * referenceNormal )
        .eval();
}

/**
 * @brief Construct the mixed Piola frame of an oriented reference facet.
 *
 * Column zero is the contravariant normal-flux direction
 * `J K^{-T} nHat |fHat|`; remaining columns are covariant tangential
 * directions `K tHat_i`. Unlike an orthonormal convenience frame, this frame
 * preserves the reference functional parametrization and must be used for
 * H(div)/H(curl) trace pairings.
 *
 * @param jacobian physical Jacobian `K`
 * @param inverseJacobianTranspose inverse Jacobian transpose `K^{-T}`
 * @param integrationJacobian positive integration Jacobian `J`
 * @param referenceNormal oriented reference unit normal
 * @param referenceTangents oriented reference tangent vectors by columns
 * @param referenceFacetMeasure measure of the reference facet
 * @return pointwise Piola trace frame
 */
template <typename JacobianDerived, typename InverseDerived, typename Scalar,
          typename NormalDerived, typename TangentDerived>
[[nodiscard]] auto piolaTraceFrame(
    Eigen::MatrixBase<JacobianDerived> const& jacobian,
    Eigen::MatrixBase<InverseDerived> const& inverseJacobianTranspose,
    Scalar const& integrationJacobian,
    Eigen::MatrixBase<NormalDerived> const& referenceNormal,
    Eigen::MatrixBase<TangentDerived> const& referenceTangents,
    Scalar const& referenceFacetMeasure )
{
    using value_type = typename JacobianDerived::Scalar;
    constexpr int rows = JacobianDerived::RowsAtCompileTime;
    constexpr int tangentColumns = TangentDerived::ColsAtCompileTime;
    constexpr int columns = tangentColumns == Eigen::Dynamic
                                ? Eigen::Dynamic
                                : tangentColumns+1;
    Eigen::Matrix<value_type,rows,columns> frame(
        jacobian.rows(), referenceTangents.cols()+1 );
    frame.col( 0 ) = contravariantPiolaTestDirection(
        inverseJacobianTranspose, integrationJacobian, referenceNormal,
        referenceFacetMeasure );
    frame.rightCols( referenceTangents.cols() ).noalias() =
        jacobian*referenceTangents;
    return frame;
}

//!
//! \brief Pair a physical H(curl) field with a reference test vector.
//!
//! Accessors keep interpolation loops allocation-free while centralizing the
//! identity `vHat . qHat = v . (K qHat)`.
//!
//! \param jacobian pointwise physical Jacobian
//! \param physicalComponent returns one physical vector component
//! \param referenceComponent returns one reference vector component
//! \param physicalDimension number of physical components
//! \param referenceDimension number of reference components
template <typename MatrixDerived, typename PhysicalAccessor,
          typename ReferenceAccessor>
[[nodiscard]] auto covariantPiolaDualPairing(
    Eigen::MatrixBase<MatrixDerived> const& jacobian,
    PhysicalAccessor&& physicalComponent,
    ReferenceAccessor&& referenceComponent,
    int physicalDimension,
    int referenceDimension )
{
    using value_type = std::decay_t<decltype( physicalComponent( 0 ) * referenceComponent( 0 ) )>;
    value_type result = 0;
    for ( int physical = 0; physical < physicalDimension; ++physical )
    {
        value_type transformedTest = 0;
        for ( int reference = 0; reference < referenceDimension; ++reference )
            transformedTest += jacobian( physical, reference ) *
                               referenceComponent( reference );
        result += physicalComponent( physical ) * transformedTest;
    }
    return result;
}

//!
//! \brief Pair a physical H(div) field with a reference test vector.
//!
//! Centralizes `vHat . qHat = J v . (K^{-T} qHat)` without allocating a
//! temporary vector in quadrature loops.
template <typename MatrixDerived, typename Scalar, typename PhysicalAccessor,
          typename ReferenceAccessor>
[[nodiscard]] auto contravariantPiolaDualPairing(
    Eigen::MatrixBase<MatrixDerived> const& inverseJacobianTranspose,
    Scalar const& integrationJacobian,
    PhysicalAccessor&& physicalComponent,
    ReferenceAccessor&& referenceComponent,
    int physicalDimension,
    int referenceDimension )
{
    using value_type = std::decay_t<decltype( physicalComponent( 0 ) * referenceComponent( 0 ) )>;
    value_type result = 0;
    for ( int physical = 0; physical < physicalDimension; ++physical )
    {
        value_type transformedTest = 0;
        for ( int reference = 0; reference < referenceDimension; ++reference )
            transformedTest += inverseJacobianTranspose( physical, reference ) *
                               referenceComponent( reference );
        result += physicalComponent( physical ) * transformedTest;
    }
    return integrationJacobian * result;
}

//!
//! \brief Map the divergence of an H(div) field to physical coordinates.
//!
//! Implements the commuting identity `div(v) = divHat(vHat) / J`.
template <typename Scalar>
[[nodiscard]] constexpr auto mappedContravariantDivergence(
    Scalar const& referenceDivergence,
    Scalar const& integrationJacobian )
{
    return referenceDivergence / integrationJacobian;
}

//!
//! \brief Map the scalar two-dimensional curl of an H(curl) field.
//!
//! Implements `curl(v) = curlHat(vHat) / J` for a full-dimensional map.
template <typename Scalar>
[[nodiscard]] constexpr auto mappedCovariantCurl2D(
    Scalar const& referenceCurl,
    Scalar const& integrationJacobian )
{
    return referenceCurl / integrationJacobian;
}

//!
//! \brief Map the three-dimensional curl of an H(curl) field.
//!
//! Implements the commuting identity `curl(v) = K curlHat(vHat) / J`.
template <typename MatrixDerived, typename Scalar, typename VectorDerived>
[[nodiscard]] auto mappedCovariantCurl3D(
    Eigen::MatrixBase<MatrixDerived> const& jacobian,
    Scalar const& integrationJacobian,
    Eigen::MatrixBase<VectorDerived> const& referenceCurl )
{
    return contravariantPiolaPushForward(
        jacobian, integrationJacobian, referenceCurl );
}

//!
//! \brief Compute the full physical gradient of a covariant Piola field.
//!
//! For `v = K^{-T} vHat`, applies the complete chain rule, including the
//! derivative of `K^{-T}`. Matrix rows are vector components and columns are
//! physical derivative directions.
//!
//! \param jacobian physical Jacobian `K`
//! \param inverseJacobianTranspose inverse Jacobian transpose `K^{-T}`
//! \param mappingHessian returns `d K(a,A) / d xi(B)`
//! \param referenceVector reference field value
//! \param referenceGradient reference gradient `(component, derivative)`
template <typename JacobianDerived, typename InverseDerived,
          typename HessianAccessor, typename VectorDerived,
          typename GradientDerived>
[[nodiscard]] auto covariantPiolaGradient(
    Eigen::MatrixBase<JacobianDerived> const& jacobian,
    Eigen::MatrixBase<InverseDerived> const& inverseJacobianTranspose,
    HessianAccessor&& mappingHessian,
    Eigen::MatrixBase<VectorDerived> const& referenceVector,
    Eigen::MatrixBase<GradientDerived> const& referenceGradient )
{
    using value_type = typename JacobianDerived::Scalar;
    Eigen::Matrix<value_type, JacobianDerived::RowsAtCompileTime,
                  JacobianDerived::RowsAtCompileTime>
        result( jacobian.rows(), jacobian.rows() );
    result.setZero();
    for ( int physicalComponent = 0;
          physicalComponent < jacobian.rows(); ++physicalComponent )
        for ( int physicalDerivative = 0;
              physicalDerivative < jacobian.rows(); ++physicalDerivative )
            for ( int referenceDerivative = 0;
                  referenceDerivative < jacobian.cols(); ++referenceDerivative )
            {
                value_type referenceDerivativeValue = 0;
                for ( int referenceComponent = 0;
                      referenceComponent < jacobian.cols(); ++referenceComponent )
                {
                    value_type inverseDerivative = 0;
                    for ( int mappedComponent = 0;
                          mappedComponent < jacobian.rows(); ++mappedComponent )
                        for ( int mappedDerivative = 0;
                              mappedDerivative < jacobian.cols(); ++mappedDerivative )
                            inverseDerivative -=
                                inverseJacobianTranspose( mappedComponent,
                                                          referenceComponent ) *
                                mappingHessian( mappedComponent,
                                                mappedDerivative,
                                                referenceDerivative ) *
                                inverseJacobianTranspose( physicalComponent,
                                                          mappedDerivative );
                    referenceDerivativeValue +=
                        inverseDerivative * referenceVector( referenceComponent ) +
                        inverseJacobianTranspose( physicalComponent,
                                                  referenceComponent ) *
                            referenceGradient( referenceComponent,
                                               referenceDerivative );
                }
                result( physicalComponent, physicalDerivative ) +=
                    inverseJacobianTranspose( physicalDerivative,
                                              referenceDerivative ) *
                    referenceDerivativeValue;
            }
    return result;
}

//!
//! \brief Compute the full physical gradient of a contravariant Piola field.
//!
//! For `v = K vHat / J`, applies the complete chain rule, including the
//! derivatives of both `K` and the positive integration Jacobian `J`.
//!
//! \param jacobian physical Jacobian `K`
//! \param inverseJacobianTranspose inverse Jacobian transpose `K^{-T}`
//! \param integrationJacobian positive integration Jacobian `J`
//! \param mappingHessian returns `d K(a,A) / d xi(B)`
//! \param referenceVector reference field value
//! \param referenceGradient reference gradient `(component, derivative)`
template <typename JacobianDerived, typename InverseDerived, typename Scalar,
          typename HessianAccessor, typename VectorDerived,
          typename GradientDerived>
[[nodiscard]] auto contravariantPiolaGradient(
    Eigen::MatrixBase<JacobianDerived> const& jacobian,
    Eigen::MatrixBase<InverseDerived> const& inverseJacobianTranspose,
    Scalar const& integrationJacobian,
    HessianAccessor&& mappingHessian,
    Eigen::MatrixBase<VectorDerived> const& referenceVector,
    Eigen::MatrixBase<GradientDerived> const& referenceGradient )
{
    using value_type = typename JacobianDerived::Scalar;
    Eigen::Matrix<value_type, JacobianDerived::RowsAtCompileTime,
                  JacobianDerived::RowsAtCompileTime>
        result( jacobian.rows(), jacobian.rows() );
    result.setZero();
    for ( int referenceDerivative = 0;
          referenceDerivative < jacobian.cols(); ++referenceDerivative )
    {
        value_type logarithmicJacobianDerivative = 0;
        for ( int mappedComponent = 0;
              mappedComponent < jacobian.rows(); ++mappedComponent )
            for ( int mappedDerivative = 0;
                  mappedDerivative < jacobian.cols(); ++mappedDerivative )
                logarithmicJacobianDerivative +=
                    inverseJacobianTranspose( mappedComponent,
                                              mappedDerivative ) *
                    mappingHessian( mappedComponent, mappedDerivative,
                                    referenceDerivative );

        for ( int physicalComponent = 0;
              physicalComponent < jacobian.rows(); ++physicalComponent )
            for ( int physicalDerivative = 0;
                  physicalDerivative < jacobian.rows(); ++physicalDerivative )
            {
                value_type referenceDerivativeValue = 0;
                for ( int referenceComponent = 0;
                      referenceComponent < jacobian.cols(); ++referenceComponent )
                    referenceDerivativeValue +=
                        mappingHessian( physicalComponent, referenceComponent,
                                        referenceDerivative ) *
                            referenceVector( referenceComponent ) +
                        jacobian( physicalComponent, referenceComponent ) *
                            ( referenceGradient( referenceComponent,
                                                 referenceDerivative ) -
                              logarithmicJacobianDerivative *
                                  referenceVector( referenceComponent ) );
                result( physicalComponent, physicalDerivative ) +=
                    inverseJacobianTranspose( physicalDerivative,
                                              referenceDerivative ) *
                    referenceDerivativeValue / integrationJacobian;
            }
    }
    return result;
}

/**
 * @brief Complete second-order chain rule for a full-dimensional Piola field.
 *
 * Cache the value map A and its first/second reference derivatives once per
 * geometry point, then reuse them for every basis function. A is K/J for
 * H(div), or K^{-T} for H(curl). Third derivatives of F differentiate K twice.
 * The returned Hessians use (component, physical derivative, physical derivative).
 */
template <typename Scalar, int Dim, bool Contravariant>
class PiolaHessianTransform
{
    static_assert( Dim == 2 || Dim == 3 );
    using matrix_type = Eigen::Matrix<Scalar, Dim, Dim>;
    using vector_type = Eigen::Matrix<Scalar, Dim, 1>;

  public:
    template <typename JacobianDerived, typename InverseDerived,
              typename HessianAccessor, typename ThirdDerivativeAccessor>
    PiolaHessianTransform(
        Eigen::MatrixBase<JacobianDerived> const& K,
        Eigen::MatrixBase<InverseDerived> const& B, Scalar const& J,
        HessianAccessor&& mappingHessian, ThirdDerivativeAccessor&& mappingThirdDerivative )
        : M_B( B )
    {
        if ( K.rows() != Dim || K.cols() != Dim )
            throw std::invalid_argument( "Piola Hessians require a full-dimensional geometry" );
        std::array<matrix_type, Dim> inverseDerivatives;
        vector_type logJacobianGradient;
        for ( int a = 0; a < Dim; ++a )
        {
            for ( int component = 0; component < Dim; ++component )
                for ( int direction = 0; direction < Dim; ++direction )
                    M_Kderivatives[a]( component, direction ) = mappingHessian( component, direction, a );
            inverseDerivatives[a] = -B * M_Kderivatives[a].transpose() * B;
            logJacobianGradient( a ) = ( B.transpose() * M_Kderivatives[a] ).trace();
        }
        if constexpr ( Contravariant )
            M_valueTransform = K / J;
        else
            M_valueTransform = B;
        for ( int a = 0; a < Dim; ++a )
        {
            if constexpr ( Contravariant )
                M_first[a] = ( M_Kderivatives[a] - logJacobianGradient( a ) * K ) / J;
            else
                M_first[a] = inverseDerivatives[a];
            for ( int b = 0; b < Dim; ++b )
            {
                matrix_type secondK;
                for ( int component = 0; component < Dim; ++component )
                    for ( int direction = 0; direction < Dim; ++direction )
                        secondK( component, direction ) = mappingThirdDerivative( component, direction, a, b );
                if constexpr ( Contravariant )
                {
                    Scalar const secondLogJacobian =
                        ( inverseDerivatives[b].transpose() * M_Kderivatives[a] + B.transpose() * secondK ).trace();
                    M_second[a * Dim + b] =
                        ( secondK - logJacobianGradient( a ) * M_Kderivatives[b] -
                          logJacobianGradient( b ) * M_Kderivatives[a] +
                          ( logJacobianGradient( a ) * logJacobianGradient( b ) - secondLogJacobian ) * K ) /
                        J;
                }
                else
                    M_second[a * Dim + b] =
                        -inverseDerivatives[b] * M_Kderivatives[a].transpose() * B -
                        B * secondK.transpose() * B - B * M_Kderivatives[a].transpose() * inverseDerivatives[b];
            }
        }
    }

    template <typename ValueAccessor, typename GradientAccessor, typename HessianAccessor>
    [[nodiscard]] std::array<matrix_type, Dim> apply(
        ValueAccessor&& value, GradientAccessor&& gradient, HessianAccessor&& hessian ) const
    {
        vector_type referenceValue;
        matrix_type referenceGradient;
        for ( int component = 0; component < Dim; ++component )
        {
            referenceValue( component ) = value( component );
            for ( int a = 0; a < Dim; ++a )
                referenceGradient( component, a ) = gradient( component, a );
        }
        matrix_type mappedReferenceGradient;
        for ( int a = 0; a < Dim; ++a )
            mappedReferenceGradient.col( a ) = M_first[a] * referenceValue + M_valueTransform * referenceGradient.col( a );
        matrix_type const physicalGradient = mappedReferenceGradient * M_B.transpose();
        std::array<matrix_type, Dim> referenceHessians;
        for ( int a = 0; a < Dim; ++a )
            for ( int b = 0; b < Dim; ++b )
            {
                vector_type referenceSecond;
                for ( int component = 0; component < Dim; ++component )
                    referenceSecond( component ) = hessian( component, a, b );
                vector_type const mappedSecond =
                    M_second[a * Dim + b] * referenceValue +
                    M_first[a] * referenceGradient.col( b ) + M_first[b] * referenceGradient.col( a ) +
                    M_valueTransform * referenceSecond;
                for ( int component = 0; component < Dim; ++component )
                    referenceHessians[component]( a, b ) = mappedSecond( component ) -
                                                           ( physicalGradient.row( component ) * M_Kderivatives[b].col( a ) )( 0, 0 );
            }
        std::array<matrix_type, Dim> result;
        for ( int component = 0; component < Dim; ++component )
            result[component] = M_B * referenceHessians[component] * M_B.transpose();
        return result;
    }

  private:
    matrix_type M_B, M_valueTransform;
    std::array<matrix_type, Dim> M_Kderivatives, M_first;
    std::array<matrix_type, Dim * Dim> M_second;
};

/**
 * @brief Differentiate the commuting H(div) identity div(v) = divHat(vHat)/J.
 *
 * Only first derivatives of the Jacobian are required, even on curved cells.
 * The returned vector is expressed in physical derivative directions.
 */
template<typename InverseDerived, typename Scalar, typename HessianAccessor,
         typename ReferenceScalar, typename GradientDerived>
[[nodiscard]] auto mappedContravariantDivergenceGradient(
    Eigen::MatrixBase<InverseDerived> const& inverseJacobianTranspose,
    Scalar const& integrationJacobian, HessianAccessor&& mappingHessian,
    ReferenceScalar const& referenceDivergence,
    Eigen::MatrixBase<GradientDerived> const& referenceDivergenceGradient )
{
    auto referenceDerivative = referenceDivergenceGradient.eval();
    for ( int derivative = 0; derivative < inverseJacobianTranspose.cols(); ++derivative )
    {
        Scalar logarithmicJacobianDerivative = 0;
        for ( int component = 0; component < inverseJacobianTranspose.rows(); ++component )
            for ( int direction = 0; direction < inverseJacobianTranspose.cols(); ++direction )
                logarithmicJacobianDerivative +=
                    inverseJacobianTranspose( component, direction ) *
                    mappingHessian( component, direction, derivative );
        referenceDerivative( derivative ) -= referenceDivergence * logarithmicJacobianDerivative;
    }
    return ( inverseJacobianTranspose * referenceDerivative / integrationJacobian ).eval();
}

/** @brief Physical curl-curl of a covariant Piola field on a 2D cell. */
template<typename InverseDerived, typename Scalar, typename HessianAccessor,
         typename ReferenceScalar, typename GradientDerived>
[[nodiscard]] auto mappedCovariantCurlCurl2D(
    Eigen::MatrixBase<InverseDerived> const& inverseJacobianTranspose,
    Scalar const& integrationJacobian, HessianAccessor&& mappingHessian,
    ReferenceScalar const& referenceCurl,
    Eigen::MatrixBase<GradientDerived> const& referenceCurlGradient )
{
    auto const gradient = mappedContravariantDivergenceGradient(
        inverseJacobianTranspose, integrationJacobian,
        std::forward<HessianAccessor>( mappingHessian ), referenceCurl,
        referenceCurlGradient );
    Eigen::Matrix<typename InverseDerived::Scalar, 2, 1> result;
    result << gradient( 1 ), -gradient( 0 );
    return result;
}

/**
 * @brief Physical curl-curl of a covariant Piola field on a 3D cell.
 *
 * Its first curl is a contravariant Piola field. Differentiating that field
 * needs the geometry Hessian, rather than the third geometry derivatives of
 * a full Hessian of the original covariant field.
 */
template<typename JacobianDerived, typename InverseDerived, typename Scalar,
         typename HessianAccessor, typename VectorDerived, typename GradientDerived>
[[nodiscard]] auto mappedCovariantCurlCurl3D(
    Eigen::MatrixBase<JacobianDerived> const& jacobian,
    Eigen::MatrixBase<InverseDerived> const& inverseJacobianTranspose,
    Scalar const& integrationJacobian, HessianAccessor&& mappingHessian,
    Eigen::MatrixBase<VectorDerived> const& referenceCurl,
    Eigen::MatrixBase<GradientDerived> const& referenceCurlGradient )
{
    auto const gradient = contravariantPiolaGradient(
        jacobian, inverseJacobianTranspose, integrationJacobian,
        std::forward<HessianAccessor>( mappingHessian ), referenceCurl,
        referenceCurlGradient );
    Eigen::Matrix<typename JacobianDerived::Scalar, 3, 1> result;
    result << gradient( 2, 1 ) - gradient( 1, 2 ),
              gradient( 0, 2 ) - gradient( 2, 0 ),
              gradient( 1, 0 ) - gradient( 0, 1 );
    return result;
}

} // namespace Feel::geometry

#endif // FEELPP_FEELPOLY_PIOLATRANSFORM_HPP
