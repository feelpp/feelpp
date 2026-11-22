/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file test_manifold_geometry.cpp
 * @brief Analytic embedded metric, curvature, and differential-operator tests.
 */
#define BOOST_TEST_MODULE test_manifold_geometry
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelalg/glas.hpp>
#include <feel/feelpoly/manifoldgeometry.hpp>

#include <Eigen/Core>
#include <Eigen/LU>

#include <cmath>
#include <tuple>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
using surface_hessian_type = tensor3_fixed_size_t<3, 2, 2, double>;

/** True only when the mapping has a supported canonical codimension-one normal. */
template <typename MatrixType>
concept has_canonical_unit_normal = requires( MatrixType const& jacobian ) {
    geometry::canonicalUnitNormal( jacobian );
};

static_assert( has_canonical_unit_normal<Eigen::Matrix<double, 3, 2>> );
static_assert( !has_canonical_unit_normal<Eigen::Matrix<double, 3, 1>> );

/** @brief Return `K(K^T K)^{-1}`, the metric pseudoinverse transpose. */
template <typename MatrixDerived>
auto inverseJacobianTranspose( Eigen::MatrixBase<MatrixDerived> const& jacobian )
{
    return ( jacobian*( jacobian.transpose()*jacobian ).inverse() ).eval();
}

/** @brief Check metric and orthogonal-projector identities. */
template <typename MatrixDerived>
void checkMetricAndProjectors(
    Eigen::MatrixBase<MatrixDerived> const& jacobian,
    double tolerance = 2e-13 )
{
    auto const inverse = inverseJacobianTranspose( jacobian );
    auto const metric = geometry::metricTensor( jacobian );
    auto const inverseMetric = geometry::inverseMetricTensor( inverse );
    auto const tangent = geometry::tangentProjector( jacobian, inverse );
    auto const normal = geometry::normalProjector( jacobian, inverse );
    auto const identity = Eigen::MatrixXd::Identity( jacobian.rows(),
                                                     jacobian.rows() );

    BOOST_CHECK_SMALL( ( metric*inverseMetric-
                         Eigen::MatrixXd::Identity( jacobian.cols(),
                                                    jacobian.cols() ) ).norm(),
                       tolerance );
    BOOST_CHECK_SMALL( ( tangent-tangent.transpose() ).norm(), tolerance );
    BOOST_CHECK_SMALL( ( tangent*tangent-tangent ).norm(), tolerance );
    BOOST_CHECK_SMALL( ( normal-normal.transpose() ).norm(), tolerance );
    BOOST_CHECK_SMALL( ( normal*normal-normal ).norm(), tolerance );
    BOOST_CHECK_SMALL( ( tangent+normal-identity ).norm(), tolerance );
    BOOST_CHECK_SMALL( ( tangent*normal ).norm(), tolerance );
    BOOST_CHECK_SMALL( ( tangent*jacobian-jacobian ).norm(), tolerance );
}

/** @brief Evaluate curvature data shared by the analytic surfaces. */
auto curvatureData( Eigen::Matrix<double, 3, 2> const& jacobian,
                    surface_hessian_type const& hessian )
{
    auto const inverse = inverseJacobianTranspose( jacobian );
    auto const normal = geometry::canonicalUnitNormal( jacobian );
    auto const second = geometry::secondFundamentalForm(
        jacobian, normal, hessian );
    auto const shape = geometry::shapeOperator(
        geometry::inverseMetricTensor( inverse ), second );
    return std::tuple{ normal, second, shape };
}
} // namespace

BOOST_AUTO_TEST_SUITE( manifold_geometry_suite )

/** @test Verify intrinsic metric and projector identities in codimension one and two. */
BOOST_AUTO_TEST_CASE( metric_projector_and_normal_frame_contracts )
{
    Eigen::Matrix<double, 3, 2> surfaceJacobian;
    surfaceJacobian << 1.0, 0.2,
                       0.1, 1.3,
                       0.4, -0.2;
    checkMetricAndProjectors( surfaceJacobian );
    BOOST_CHECK_CLOSE_FRACTION(
        geometry::manifoldJacobian( surfaceJacobian ),
        std::sqrt( geometry::metricTensor( surfaceJacobian ).determinant() ),
        2e-14 );

    auto const unitNormal = geometry::canonicalUnitNormal( surfaceJacobian );
    BOOST_CHECK_CLOSE_FRACTION( unitNormal.norm(), 1.0, 2e-14 );
    BOOST_CHECK_SMALL( ( surfaceJacobian.transpose()*unitNormal ).norm(),
                       2e-14 );

    Eigen::Matrix<double, 3, 1> spaceCurveJacobian;
    spaceCurveJacobian << 1.0, 2.0, -0.5;
    checkMetricAndProjectors( spaceCurveJacobian );
    auto const normalFrame = geometry::normalFrame(
        spaceCurveJacobian, geometry::qr_normal_frame_policy{} );
    BOOST_CHECK_EQUAL( normalFrame.cols(), 2 );
    BOOST_CHECK_SMALL(
        ( normalFrame.transpose()*normalFrame-
          Eigen::Matrix2d::Identity() ).norm(), 2e-14 );
    BOOST_CHECK_SMALL( ( spaceCurveJacobian.transpose()*normalFrame ).norm(),
                       2e-14 );
}

/** @test Verify signed sphere curvature and a manufactured Laplace--Beltrami value. */
BOOST_AUTO_TEST_CASE( sphere_curvature_and_laplace_beltrami )
{
    // Unit sphere at theta=pi/2, phi=0, oriented by (theta,phi).
    Eigen::Matrix<double, 3, 2> jacobian;
    jacobian << 0.0, 0.0,
                0.0, 1.0,
               -1.0, 0.0;
    surface_hessian_type hessian;
    hessian.setZero();
    hessian( 0, 0, 0 ) = -1.0;
    hessian( 0, 1, 1 ) = -1.0;

    auto const [normal, second, shape] = curvatureData( jacobian, hessian );
    BOOST_CHECK_SMALL( ( normal-Eigen::Vector3d::UnitX() ).norm(), 2e-14 );
    BOOST_CHECK_SMALL( ( second+Eigen::Matrix2d::Identity() ).norm(), 2e-14 );
    BOOST_CHECK_SMALL( ( shape+Eigen::Matrix2d::Identity() ).norm(), 2e-14 );
    BOOST_CHECK_CLOSE_FRACTION( geometry::meanCurvature( shape ), -1.0,
                                2e-14 );
    BOOST_CHECK_CLOSE_FRACTION( geometry::gaussianCurvature( shape ), 1.0,
                                2e-14 );

    // u(theta,phi)=sin(theta)cos(phi)=x, so Delta_Gamma u=-2u=-2.
    Eigen::Vector2d referenceGradient = Eigen::Vector2d::Zero();
    Eigen::Matrix2d referenceHessian = -Eigen::Matrix2d::Identity();
    auto const inverse = inverseJacobianTranspose( jacobian );
    BOOST_CHECK_CLOSE_FRACTION(
        geometry::laplaceBeltrami(
            jacobian, geometry::inverseMetricTensor( inverse ), hessian,
            referenceGradient, referenceHessian ),
        -2.0, 2e-14 );
}

/** @test Verify cylinder and torus principal-curvature conventions. */
BOOST_AUTO_TEST_CASE( cylinder_and_torus_curvature )
{
    Eigen::Matrix<double, 3, 2> cylinderJacobian;
    cylinderJacobian << 0.0, 0.0,
                        1.0, 0.0,
                        0.0, 1.0;
    surface_hessian_type cylinderHessian;
    cylinderHessian.setZero();
    cylinderHessian( 0, 0, 0 ) = -1.0;
    auto const [cylinderNormal, cylinderSecond, cylinderShape] =
        curvatureData( cylinderJacobian, cylinderHessian );
    BOOST_CHECK_SMALL( ( cylinderNormal-Eigen::Vector3d::UnitX() ).norm(),
                       2e-14 );
    BOOST_CHECK_CLOSE_FRACTION( geometry::meanCurvature( cylinderShape ),
                                -0.5, 2e-14 );
    BOOST_CHECK_SMALL( geometry::gaussianCurvature( cylinderShape ), 2e-14 );

    constexpr double majorRadius = 2.0;
    constexpr double minorRadius = 0.5;
    Eigen::Matrix<double, 3, 2> torusJacobian;
    torusJacobian << 0.0, 0.0,
                     majorRadius+minorRadius, 0.0,
                     0.0, minorRadius;
    surface_hessian_type torusHessian;
    torusHessian.setZero();
    torusHessian( 0, 0, 0 ) = -( majorRadius+minorRadius );
    torusHessian( 0, 1, 1 ) = -minorRadius;
    auto const [torusNormal, torusSecond, torusShape] =
        curvatureData( torusJacobian, torusHessian );
    BOOST_CHECK_SMALL( ( torusNormal-Eigen::Vector3d::UnitX() ).norm(),
                       2e-14 );
    BOOST_CHECK_CLOSE_FRACTION(
        geometry::meanCurvature( torusShape ),
        -0.5*( 1.0/( majorRadius+minorRadius )+1.0/minorRadius ),
        2e-14 );
    BOOST_CHECK_CLOSE_FRACTION(
        geometry::gaussianCurvature( torusShape ),
        1.0/( minorRadius*( majorRadius+minorRadius ) ), 2e-14 );
}

/** @test Verify component conventions of the pointwise surface operators. */
BOOST_AUTO_TEST_CASE( surface_differential_operator_conventions )
{
    Eigen::Matrix3d tangent = Eigen::Matrix3d::Zero();
    tangent.diagonal() << 1.0, 1.0, 0.0;
    Eigen::Vector3d ambientScalarGradient( 2.0, -3.0, 5.0 );
    auto const tangential = geometry::tangentialGradient(
        tangent, ambientScalarGradient );
    BOOST_CHECK_SMALL( ( tangential-Eigen::Vector3d( 2.0, -3.0, 0.0 ) ).norm(),
                       2e-14 );

    Eigen::Matrix3d ambientVectorGradient;
    ambientVectorGradient << 1.0, 2.0, 3.0,
                             4.0, 5.0, 6.0,
                             7.0, 8.0, 9.0;
    BOOST_CHECK_CLOSE_FRACTION(
        geometry::surfaceDivergence( tangent, ambientVectorGradient ),
        6.0, 2e-14 );
    auto const covariant = geometry::covariantDerivative(
        tangent, ambientVectorGradient );
    Eigen::Matrix3d expectedCovariant = Eigen::Matrix3d::Zero();
    expectedCovariant.topLeftCorner<2, 2>() =
        ambientVectorGradient.topLeftCorner<2, 2>();
    BOOST_CHECK_SMALL( ( covariant-expectedCovariant ).norm(), 2e-14 );

    Eigen::Vector3d const normal = Eigen::Vector3d::UnitZ();
    BOOST_CHECK_CLOSE_FRACTION(
        geometry::surfaceCurl( normal, ambientVectorGradient ), 2.0, 2e-14 );
    BOOST_CHECK_SMALL(
        ( geometry::surfaceCurl( normal, tangential )-
          Eigen::Vector3d( 3.0, 2.0, 0.0 ) ).norm(), 2e-14 );
}

BOOST_AUTO_TEST_SUITE_END()
