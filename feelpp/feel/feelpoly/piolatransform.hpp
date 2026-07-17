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

#include <type_traits>
#include <utility>

namespace Feel::geometry
{
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

} // namespace Feel::geometry

#endif // FEELPP_FEELPOLY_PIOLATRANSFORM_HPP
