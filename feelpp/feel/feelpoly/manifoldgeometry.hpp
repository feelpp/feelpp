/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file manifoldgeometry.hpp
 * @brief Pointwise differential geometry kernels for embedded mappings.
 */
#if !defined( FEELPP_FEELPOLY_MANIFOLDGEOMETRY_HPP )
#define FEELPP_FEELPOLY_MANIFOLDGEOMETRY_HPP 1

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <Eigen/QR>

#include <algorithm>
#include <cmath>
#include <concepts>
#include <stdexcept>
#include <type_traits>

namespace Feel::geometry
{
/**
 * @brief Return the first fundamental form `g = K^T K`.
 * @param jacobian physical mapping Jacobian, with reference derivatives by columns
 * @return symmetric positive definite metric tensor for a regular map
 */
template <typename MatrixDerived>
[[nodiscard]] auto metricTensor(
    Eigen::MatrixBase<MatrixDerived> const& jacobian )
{
    return ( jacobian.transpose()*jacobian ).eval();
}

/**
 * @brief Return the inverse metric from the metric pseudoinverse transpose.
 *
 * If `B = K (K^T K)^{-1}`, then `B^T B = (K^T K)^{-1}`.
 *
 * @param inverseJacobianTranspose `B`, equal to `K^{-T}` for square maps
 * @return inverse first fundamental form
 */
template <typename MatrixDerived>
[[nodiscard]] auto inverseMetricTensor(
    Eigen::MatrixBase<MatrixDerived> const& inverseJacobianTranspose )
{
    return ( inverseJacobianTranspose.transpose()*
             inverseJacobianTranspose ).eval();
}

/**
 * @brief Return the positive manifold integration Jacobian `sqrt(det(g))`.
 * @param metric first fundamental form
 * @return positive local measure scaling
 */
template <typename MatrixDerived>
[[nodiscard]] auto manifoldJacobianFromMetric(
    Eigen::MatrixBase<MatrixDerived> const& metric )
{
    using value_type = typename MatrixDerived::Scalar;
    using std::sqrt;
    return sqrt( std::max( value_type{}, metric.determinant() ) );
}

/**
 * @brief Return the positive manifold integration Jacobian of a mapping.
 * @param jacobian physical mapping Jacobian
 * @return `sqrt(det(K^T K))`
 */
template <typename MatrixDerived>
[[nodiscard]] auto manifoldJacobian(
    Eigen::MatrixBase<MatrixDerived> const& jacobian )
{
    return manifoldJacobianFromMetric( metricTensor( jacobian ) );
}

/**
 * @brief Return the orthogonal tangent projector `P = K B^T`.
 * @param jacobian physical mapping Jacobian
 * @param inverseJacobianTranspose metric pseudoinverse transpose
 * @return symmetric idempotent projector onto the mapped tangent space
 */
template <typename JacobianDerived, typename InverseDerived>
[[nodiscard]] auto tangentProjector(
    Eigen::MatrixBase<JacobianDerived> const& jacobian,
    Eigen::MatrixBase<InverseDerived> const& inverseJacobianTranspose )
{
    return ( jacobian*inverseJacobianTranspose.transpose() ).eval();
}

/**
 * @brief Return the orthogonal normal projector `I-P`.
 * @param tangentProjectorValue orthogonal tangent projector
 * @return symmetric idempotent projector onto the normal space
 */
template <typename MatrixDerived>
[[nodiscard]] auto normalProjectorFromTangent(
    Eigen::MatrixBase<MatrixDerived> const& tangentProjectorValue )
{
    using value_type = typename MatrixDerived::Scalar;
    constexpr int rows = MatrixDerived::RowsAtCompileTime;
    Eigen::Matrix<value_type, rows, rows> projector(
        tangentProjectorValue.rows(), tangentProjectorValue.cols() );
    projector.setIdentity();
    projector -= tangentProjectorValue;
    return projector;
}

/**
 * @brief Return the orthogonal projector onto the normal space.
 * @param jacobian physical mapping Jacobian
 * @param inverseJacobianTranspose metric pseudoinverse transpose
 * @return `I-K B^T`
 */
template <typename JacobianDerived, typename InverseDerived>
[[nodiscard]] auto normalProjector(
    Eigen::MatrixBase<JacobianDerived> const& jacobian,
    Eigen::MatrixBase<InverseDerived> const& inverseJacobianTranspose )
{
    return normalProjectorFromTangent(
        tangentProjector( jacobian, inverseJacobianTranspose ) );
}

/**
 * @brief Return the oriented unit normal of a codimension-one curve or surface.
 *
 * Curves in the plane use the left normal `(-t_y,t_x)`. Surfaces in three
 * dimensions use `K_0 x K_1`. Higher codimension deliberately has no overload:
 * it has a normal space, not a canonical single normal.
 *
 * @param jacobian physical mapping Jacobian
 * @return oriented unit normal
 */
template <typename MatrixDerived>
    requires( ( MatrixDerived::RowsAtCompileTime == 2 &&
                MatrixDerived::ColsAtCompileTime == 1 ) ||
              ( MatrixDerived::RowsAtCompileTime == 3 &&
                MatrixDerived::ColsAtCompileTime == 2 ) )
[[nodiscard]] auto canonicalUnitNormal(
    Eigen::MatrixBase<MatrixDerived> const& jacobian )
{
    using value_type = typename MatrixDerived::Scalar;
    constexpr int physicalDimension = MatrixDerived::RowsAtCompileTime;
    Eigen::Matrix<value_type, physicalDimension, 1> normal;
    if constexpr ( physicalDimension == 2 )
    {
        normal << -jacobian( 1, 0 ), jacobian( 0, 0 );
    }
    else
        normal = jacobian.col( 0 ).cross( jacobian.col( 1 ) );

    value_type const norm = normal.norm();
    if ( norm <= Eigen::NumTraits<value_type>::dummy_precision() )
        throw std::domain_error(
            "A unit normal requires a full-rank codimension-one Jacobian" );
    return ( normal/norm ).eval();
}

/**
 * @struct canonical_normal_frame_policy
 * @brief Normal-frame policy for oriented codimension-one mappings.
 */
struct canonical_normal_frame_policy
{
    /** @brief Return the canonical normal as a one-column frame. */
    template <typename MatrixDerived>
    [[nodiscard]] auto operator()(
        Eigen::MatrixBase<MatrixDerived> const& jacobian ) const
    {
        auto const normal = canonicalUnitNormal( jacobian );
        return normal.eval();
    }
};

/**
 * @struct qr_normal_frame_policy
 * @brief Deterministic orthonormal basis of a higher-codimension normal space.
 *
 * Individual columns are a numerical frame and are not geometrically
 * canonical. Algorithms invariant under a change of normal frame should use
 * the normal projector whenever possible.
 */
struct qr_normal_frame_policy
{
    /** @brief Complete the tangent columns to an orthonormal ambient frame. */
    template <typename MatrixDerived>
    [[nodiscard]] auto operator()(
        Eigen::MatrixBase<MatrixDerived> const& jacobian ) const
    {
        using matrix_type = typename MatrixDerived::PlainObject;
        using value_type = typename MatrixDerived::Scalar;
        constexpr int physicalDimension = MatrixDerived::RowsAtCompileTime;
        constexpr int referenceDimension = MatrixDerived::ColsAtCompileTime;
        constexpr int normalDimension =
            physicalDimension == Eigen::Dynamic ||
                    referenceDimension == Eigen::Dynamic
                ? Eigen::Dynamic
                : physicalDimension-referenceDimension;

        if ( jacobian.rows() <= jacobian.cols() )
            throw std::invalid_argument(
                "A normal frame requires an embedded full-rank Jacobian" );
        Eigen::HouseholderQR<matrix_type> factorization( jacobian );
        Eigen::Matrix<value_type, physicalDimension, physicalDimension> frame(
            jacobian.rows(), jacobian.rows() );
        frame = factorization.householderQ()*
                Eigen::Matrix<value_type, physicalDimension,
                              physicalDimension>::Identity(
                    jacobian.rows(), jacobian.rows() );
        Eigen::Matrix<value_type, physicalDimension, normalDimension> normals(
            jacobian.rows(), jacobian.rows()-jacobian.cols() );
        normals = frame.rightCols( normals.cols() );
        return normals;
    }
};

/**
 * @brief Build a normal frame using an explicit policy.
 * @param jacobian physical mapping Jacobian
 * @param policy normal-frame construction policy
 * @return normal vectors stored by columns
 */
template <typename MatrixDerived, typename Policy = qr_normal_frame_policy>
    requires std::invocable<Policy, Eigen::MatrixBase<MatrixDerived> const&>
[[nodiscard]] auto normalFrame(
    Eigen::MatrixBase<MatrixDerived> const& jacobian,
    Policy policy = {} )
{
    return policy( jacobian );
}

/**
 * @brief Compute the codimension-one second fundamental form.
 *
 * The convention is `b_ij = n . d_ij F`; consequently an outward-oriented
 * sphere has negative principal curvatures.
 *
 * @param jacobian physical mapping Jacobian
 * @param unitNormal oriented unit normal
 * @param mappingHessian mapping Hessian indexed `(ambient,i,j)`
 * @return symmetric reference-coordinate second fundamental form
 */
template <typename JacobianDerived, typename NormalDerived, typename Hessian>
[[nodiscard]] auto secondFundamentalForm(
    Eigen::MatrixBase<JacobianDerived> const& jacobian,
    Eigen::MatrixBase<NormalDerived> const& unitNormal,
    Hessian const& mappingHessian )
{
    using value_type = typename NormalDerived::Scalar;
    constexpr int referenceDimension =
        JacobianDerived::ColsAtCompileTime;
    Eigen::Matrix<value_type, referenceDimension, referenceDimension> form;
    for ( int i = 0; i < jacobian.cols(); ++i )
        for ( int j = 0; j < jacobian.cols(); ++j )
        {
            value_type value{};
            for ( int component = 0; component < unitNormal.rows();
                  ++component )
                value += unitNormal( component )*
                         mappingHessian( component, i, j );
            form( i, j ) = value;
        }
    return form;
}

/**
 * @brief Compute the Weingarten shape operator `S = g^{-1} b`.
 * @param inverseMetric inverse first fundamental form
 * @param secondFundamentalFormValue second fundamental form
 * @return reference-coordinate shape operator
 */
template <typename MetricDerived, typename FormDerived>
[[nodiscard]] auto shapeOperator(
    Eigen::MatrixBase<MetricDerived> const& inverseMetric,
    Eigen::MatrixBase<FormDerived> const& secondFundamentalFormValue )
{
    return ( inverseMetric*secondFundamentalFormValue ).eval();
}

/**
 * @brief Return the averaged signed mean curvature `trace(S)/dim(Gamma)`.
 * @param shapeOperatorValue shape operator
 * @return signed mean curvature
 */
template <typename MatrixDerived>
[[nodiscard]] auto meanCurvature(
    Eigen::MatrixBase<MatrixDerived> const& shapeOperatorValue )
{
    return shapeOperatorValue.trace()/shapeOperatorValue.rows();
}

/**
 * @brief Return Gaussian curvature of a two-dimensional surface.
 * @param shapeOperatorValue two-by-two shape operator
 * @return determinant of the shape operator
 */
template <typename MatrixDerived>
    requires( MatrixDerived::RowsAtCompileTime == 2 &&
              MatrixDerived::ColsAtCompileTime == 2 )
[[nodiscard]] auto gaussianCurvature(
    Eigen::MatrixBase<MatrixDerived> const& shapeOperatorValue )
{
    return shapeOperatorValue.determinant();
}

/**
 * @brief Project an ambient scalar gradient onto a tangent space.
 * @param tangentProjectorValue tangent projector
 * @param ambientGradient ambient column gradient
 * @return tangential gradient
 */
template <typename ProjectorDerived, typename GradientDerived>
[[nodiscard]] auto tangentialGradient(
    Eigen::MatrixBase<ProjectorDerived> const& tangentProjectorValue,
    Eigen::MatrixBase<GradientDerived> const& ambientGradient )
{
    return ( tangentProjectorValue*ambientGradient ).eval();
}

/**
 * @brief Return intrinsic divergence of an ambient vector extension.
 * @param tangentProjectorValue tangent projector
 * @param ambientGradient matrix `(component,physical derivative)`
 * @return `P : grad(v)`
 */
template <typename ProjectorDerived, typename GradientDerived>
[[nodiscard]] auto surfaceDivergence(
    Eigen::MatrixBase<ProjectorDerived> const& tangentProjectorValue,
    Eigen::MatrixBase<GradientDerived> const& ambientGradient )
{
    return ( tangentProjectorValue*ambientGradient ).trace();
}

/**
 * @brief Return the covariant derivative of a tangential vector extension.
 * @param tangentProjectorValue tangent projector
 * @param ambientGradient matrix `(component,physical derivative)`
 * @return `P grad(v) P`
 */
template <typename ProjectorDerived, typename GradientDerived>
[[nodiscard]] auto covariantDerivative(
    Eigen::MatrixBase<ProjectorDerived> const& tangentProjectorValue,
    Eigen::MatrixBase<GradientDerived> const& ambientGradient )
{
    return ( tangentProjectorValue*ambientGradient*
             tangentProjectorValue ).eval();
}

/**
 * @brief Return the oriented scalar surface curl of a vector in `R^3`.
 * @param unitNormal oriented surface normal
 * @param ambientGradient matrix `(component,physical derivative)`
 * @return `n . curl(v)`
 */
template <typename NormalDerived, typename GradientDerived>
    requires( NormalDerived::RowsAtCompileTime == 3 &&
              GradientDerived::RowsAtCompileTime == 3 &&
              GradientDerived::ColsAtCompileTime == 3 )
[[nodiscard]] auto surfaceCurl(
    Eigen::MatrixBase<NormalDerived> const& unitNormal,
    Eigen::MatrixBase<GradientDerived> const& ambientGradient )
{
    using value_type = typename NormalDerived::Scalar;
    Eigen::Matrix<value_type, 3, 1> curl;
    curl << ambientGradient( 2, 1 )-ambientGradient( 1, 2 ),
            ambientGradient( 0, 2 )-ambientGradient( 2, 0 ),
            ambientGradient( 1, 0 )-ambientGradient( 0, 1 );
    return unitNormal.dot( curl );
}

/**
 * @brief Return the oriented surface curl of a scalar in `R^3`.
 * @param unitNormal oriented surface normal
 * @param tangentialGradientValue tangential scalar gradient
 * @return `n x grad_Gamma(u)`
 */
template <typename NormalDerived, typename GradientDerived>
    requires( NormalDerived::RowsAtCompileTime == 3 &&
              GradientDerived::RowsAtCompileTime == 3 &&
              GradientDerived::ColsAtCompileTime == 1 )
[[nodiscard]] auto surfaceCurl(
    Eigen::MatrixBase<NormalDerived> const& unitNormal,
    Eigen::MatrixBase<GradientDerived> const& tangentialGradientValue )
{
    return unitNormal.cross( tangentialGradientValue ).eval();
}

/**
 * @brief Compute the scalar Laplace--Beltrami operator from reference data.
 *
 * This evaluates
 * `g^{ij}(d_ij u-Gamma^k_ij d_k u)` without materializing Christoffel
 * symbols. It is invariant under ambient extensions of the scalar field.
 *
 * @param jacobian mapping Jacobian `K`
 * @param inverseMetric inverse metric `g^{-1}`
 * @param mappingHessian mapping Hessian `(ambient,i,j)`
 * @param referenceGradient scalar reference gradient
 * @param referenceHessian scalar reference Hessian
 * @return intrinsic scalar Laplacian
 */
template <typename JacobianDerived, typename MetricDerived, typename Hessian,
          typename GradientDerived, typename ScalarHessianDerived>
[[nodiscard]] auto laplaceBeltrami(
    Eigen::MatrixBase<JacobianDerived> const& jacobian,
    Eigen::MatrixBase<MetricDerived> const& inverseMetric,
    Hessian const& mappingHessian,
    Eigen::MatrixBase<GradientDerived> const& referenceGradient,
    Eigen::MatrixBase<ScalarHessianDerived> const& referenceHessian )
{
    using value_type = typename JacobianDerived::Scalar;
    const int referenceDimension = jacobian.cols();
    value_type result{};
    for ( int i = 0; i < referenceDimension; ++i )
        for ( int j = 0; j < referenceDimension; ++j )
        {
            value_type covariantHessian = referenceHessian( i, j );
            for ( int k = 0; k < referenceDimension; ++k )
            {
                value_type christoffel{};
                for ( int l = 0; l < referenceDimension; ++l )
                {
                    value_type tangentDotHessian{};
                    for ( int component = 0; component < jacobian.rows();
                          ++component )
                        tangentDotHessian +=
                            jacobian( component, l )*
                            mappingHessian( component, i, j );
                    christoffel += inverseMetric( k, l )*
                                   tangentDotHessian;
                }
                covariantHessian -= christoffel*referenceGradient( k );
            }
            result += inverseMetric( i, j )*covariantHessian;
        }
    return result;
}

} // namespace Feel::geometry

#endif // FEELPP_FEELPOLY_MANIFOLDGEOMETRY_HPP
