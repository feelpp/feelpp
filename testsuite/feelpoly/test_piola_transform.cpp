/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file test_piola_transform.cpp
 * @brief Pointwise covariant and contravariant Piola identity tests.
 */
#define BOOST_TEST_MODULE test_piola_transform
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelpoly/piolatransform.hpp>

#include <Eigen/Core>
#include <Eigen/LU>

#include <cmath>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
/** Check both Piola maps and their dual functional directions. */
template <typename JacobianType, typename ReferenceVectorType>
void checkPiolaIdentities( JacobianType const& jacobian,
                           ReferenceVectorType const& referenceVector,
                           ReferenceVectorType const& referenceTest )
{
    using value_type = typename JacobianType::Scalar;
    auto const metric = jacobian.transpose() * jacobian;
    auto const inverseJacobianTranspose =
        jacobian * metric.inverse();
    value_type const integrationJacobian = std::sqrt( metric.determinant() );

    auto const hcurlPhysical = geometry::covariantPiolaPushForward(
        inverseJacobianTranspose, referenceVector );
    auto const hcurlReference = geometry::covariantPiolaPullback(
        jacobian, hcurlPhysical );
    BOOST_CHECK_SMALL( ( hcurlReference - referenceVector ).norm(), 2e-13 );

    auto const tangentDirection = geometry::covariantPiolaTestDirection(
        jacobian, referenceTest );
    BOOST_CHECK_SMALL(
        hcurlPhysical.dot( tangentDirection ) -
            referenceVector.dot( referenceTest ),
        2e-13 );

    auto const hdivPhysical = geometry::contravariantPiolaPushForward(
        jacobian, integrationJacobian, referenceVector );
    auto const hdivReference = geometry::contravariantPiolaPullback(
        inverseJacobianTranspose, integrationJacobian, hdivPhysical );
    BOOST_CHECK_SMALL( ( hdivReference - referenceVector ).norm(), 2e-13 );

    value_type const referenceFacetMeasure = value_type( 1.7 );
    auto const normalDirection = geometry::contravariantPiolaTestDirection(
        inverseJacobianTranspose, integrationJacobian, referenceTest,
        referenceFacetMeasure );
    BOOST_CHECK_SMALL(
        hdivPhysical.dot( normalDirection ) -
            referenceFacetMeasure * referenceVector.dot( referenceTest ),
        2e-13 );

    auto const covariantPairing = geometry::covariantPiolaDualPairing(
        jacobian,
        [&]( int component )
        { return hcurlPhysical( component ); },
        [&]( int component )
        { return referenceTest( component ); },
        jacobian.rows(), jacobian.cols() );
    BOOST_CHECK_SMALL(
        covariantPairing - referenceVector.dot( referenceTest ), 2e-13 );

    auto const contravariantPairing = geometry::contravariantPiolaDualPairing(
        inverseJacobianTranspose, integrationJacobian,
        [&]( int component )
        { return hdivPhysical( component ); },
        [&]( int component )
        { return referenceTest( component ); },
        jacobian.rows(), jacobian.cols() );
    BOOST_CHECK_SMALL(
        contravariantPairing - referenceVector.dot( referenceTest ), 2e-13 );
}
} // namespace

BOOST_AUTO_TEST_SUITE( piola_transform_suite )

/** Verify full-dimensional affine Piola identities. */
BOOST_AUTO_TEST_CASE( affine_full_dimensional_identities )
{
    Eigen::Matrix2d jacobian;
    jacobian << 2.0, 0.3,
        0.1, 1.5;
    Eigen::Vector2d referenceVector( 0.7, -1.2 );
    Eigen::Vector2d referenceTest( -0.4, 0.9 );
    checkPiolaIdentities( jacobian, referenceVector, referenceTest );
}

/** Verify metric Piola identities for an embedded surface. */
BOOST_AUTO_TEST_CASE( affine_embedded_surface_identities )
{
    Eigen::Matrix<double, 3, 2> jacobian;
    jacobian << 1.2, 0.1,
        -0.2, 1.1,
        0.4, 0.3;
    Eigen::Vector2d referenceVector( 0.8, -0.6 );
    Eigen::Vector2d referenceTest( 0.25, 1.4 );
    checkPiolaIdentities( jacobian, referenceVector, referenceTest );
}

/** Verify that curved pointwise Jacobians produce varying but exact maps. */
BOOST_AUTO_TEST_CASE( nonaffine_pointwise_identities )
{
    Eigen::Vector2d referenceVector( 0.6, -0.35 );
    Eigen::Vector2d referenceTest( -0.2, 1.1 );
    Eigen::Vector2d firstTangent;
    for ( int point = 0; point < 4; ++point )
    {
        double const coordinate = -0.7 + 0.4 * point;
        Eigen::Matrix2d jacobian;
        jacobian << 1.4 + 0.2 * coordinate, 0.15 * coordinate,
            0.1 - 0.05 * coordinate, 1.1 + 0.1 * coordinate * coordinate;
        checkPiolaIdentities( jacobian, referenceVector, referenceTest );
        auto const tangent = geometry::covariantPiolaTestDirection(
            jacobian, referenceTest );
        if ( point == 0 )
            firstTangent = tangent;
        else
            BOOST_CHECK_GT( ( tangent - firstTangent ).norm(), 1e-4 );
    }
}

/** Verify the commuting divergence and two-dimensional curl identities. */
BOOST_AUTO_TEST_CASE( mapped_scalar_differential_operators )
{
    double const jacobian = 2.75;
    double const referenceDivergence = -1.4;
    double const referenceCurl = 0.85;
    BOOST_CHECK_CLOSE(
        geometry::mappedContravariantDivergence(
            referenceDivergence, jacobian ),
        referenceDivergence / jacobian, 1e-12 );
    BOOST_CHECK_CLOSE(
        geometry::mappedCovariantCurl2D( referenceCurl, jacobian ),
        referenceCurl / jacobian, 1e-12 );
}

/** Verify the three-dimensional mapped curl at varying curved-map points. */
BOOST_AUTO_TEST_CASE( mapped_three_dimensional_curl )
{
    Eigen::Vector3d referenceCurl( 0.4, -0.7, 1.2 );
    Eigen::Vector3d firstPhysicalCurl;
    for ( int point = 0; point < 3; ++point )
    {
        double const coordinate = -0.5 + 0.5 * point;
        Eigen::Matrix3d jacobian;
        jacobian << 1.3 + 0.1 * coordinate, 0.2, 0.05 * coordinate,
            0.1, 1.1 - 0.08 * coordinate, 0.15,
            0.03, 0.12 * coordinate, 1.4 + 0.05 * coordinate;
        double const determinant = std::abs( jacobian.determinant() );
        auto const physicalCurl = geometry::mappedCovariantCurl3D(
            jacobian, determinant, referenceCurl );
        BOOST_CHECK_SMALL(
            ( physicalCurl - jacobian * referenceCurl / determinant ).norm(),
            2e-13 );
        if ( point == 0 )
            firstPhysicalCurl = physicalCurl;
        else
            BOOST_CHECK_GT( ( physicalCurl - firstPhysicalCurl ).norm(), 1e-4 );
    }
}

/** Verify complete curved Piola gradients against centered differences. */
BOOST_AUTO_TEST_CASE( curved_full_gradient_chain_rule )
{
    auto jacobian = []( Eigen::Vector2d const& point )
    {
        Eigen::Matrix2d result;
        result << 1.0 + 0.3*point( 1 ) + 0.4*point( 0 ),
            0.3*point( 0 ),
            0.5*point( 0 ),
            1.0 + 0.2*point( 1 );
        return result;
    };
    auto referenceValue = []( Eigen::Vector2d const& point )
    {
        return Eigen::Vector2d( 0.7 + 0.4*point( 0 ) - 0.2*point( 1 ),
                                -0.3 + 0.1*point( 0 ) + 0.6*point( 1 ) );
    };
    Eigen::Matrix2d referenceGradient;
    referenceGradient << 0.4, -0.2, 0.1, 0.6;
    auto mappingHessian = []( int component, int derivative1,
                              int derivative2 )
    {
        static constexpr double values[2][2][2] = {
            { { 0.4, 0.3 }, { 0.3, 0.0 } },
            { { 0.5, 0.0 }, { 0.0, 0.2 } }
        };
        return values[component][derivative1][derivative2];
    };

    Eigen::Vector2d const point( 0.17, -0.23 );
    Eigen::Matrix2d const K = jacobian( point );
    Eigen::Matrix2d const B = K.inverse().transpose();
    double const J = std::abs( K.determinant() );
    Eigen::Vector2d const value = referenceValue( point );
    auto const hcurlGradient = geometry::covariantPiolaGradient(
        K, B, mappingHessian, value, referenceGradient );
    auto const hdivGradient = geometry::contravariantPiolaGradient(
        K, B, J, mappingHessian, value, referenceGradient );

    constexpr double epsilon = 2e-7;
    Eigen::Matrix2d hcurlReferenceDerivative;
    Eigen::Matrix2d hdivReferenceDerivative;
    for ( int derivative = 0; derivative < 2; ++derivative )
    {
        Eigen::Vector2d plus = point;
        Eigen::Vector2d minus = point;
        plus( derivative ) += epsilon;
        minus( derivative ) -= epsilon;
        auto const Kplus = jacobian( plus );
        auto const Kminus = jacobian( minus );
        hcurlReferenceDerivative.col( derivative ) =
            ( Kplus.inverse().transpose()*referenceValue( plus ) -
              Kminus.inverse().transpose()*referenceValue( minus ) ) /
            ( 2*epsilon );
        hdivReferenceDerivative.col( derivative ) =
            ( Kplus*referenceValue( plus )/std::abs( Kplus.determinant() ) -
              Kminus*referenceValue( minus )/std::abs( Kminus.determinant() ) ) /
            ( 2*epsilon );
    }
    Eigen::Matrix2d const expectedHcurl =
        hcurlReferenceDerivative*K.inverse();
    Eigen::Matrix2d const expectedHdiv =
        hdivReferenceDerivative*K.inverse();
    BOOST_CHECK_SMALL( ( hcurlGradient-expectedHcurl ).norm(), 2e-9 );
    BOOST_CHECK_SMALL( ( hdivGradient-expectedHdiv ).norm(), 2e-9 );
}

BOOST_AUTO_TEST_SUITE_END()
