/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE quadrature_batch
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelpoly/pointsetquadrature.hpp>
#include <feel/feelpoly/quadratureevaluation.hpp>
#include <feel/feelvf/detail/localform_tensor.hpp>

#include <array>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{

/** @brief Minimal unit-Jacobian context used to initialize quadrature weights. */
struct unit_jacobian_context
{
    /** @return unit volume Jacobian at every quadrature point. */
    double J( std::uint32_t /*q*/ ) const noexcept { return 1.; }

    /** @return unit boundary normal norm at every quadrature point. */
    double normalNorm( std::uint32_t /*q*/ ) const noexcept { return 1.; }
};

/**
 * @brief Scalar reference tensor implementing only the pointwise contract.
 */
struct scalar_tensor
{
    using value_type = double;
    using table_type = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                     Eigen::RowMajor>;

    explicit scalar_tensor( Eigen::Index points )
        : coefficient( points ), testValues( 3, points ), trialValues( 3, points )
    {
        for ( Eigen::Index q = 0; q < points; ++q )
        {
            coefficient[q] = 1. + 0.125*static_cast<double>( q );
            for ( Eigen::Index i = 0; i < 3; ++i )
            {
                testValues( i, q ) = 0.2*( i + 1 ) + 0.03*q;
                trialValues( i, q ) = 0.15*( i + 2 ) - 0.02*q;
            }
        }
    }

    /** @return bilinear tensor value at one quadrature point. */
    double evalijq( std::uint32_t i, std::uint32_t j,
                    std::uint16_t /*c1*/, std::uint16_t /*c2*/,
                    std::uint32_t q ) const
    {
        ++scalarCalls;
        return coefficient[q]*testValues( i, q )*trialValues( j, q );
    }

    /** @return test tensor value at one quadrature point. */
    double evaliq( std::uint32_t i, std::uint16_t /*c1*/,
                   std::uint16_t /*c2*/, std::uint32_t q ) const
    {
        ++scalarCalls;
        return coefficient[q]*testValues( i, q );
    }

    /** @return coefficient value at one quadrature point. */
    double evalq( std::uint16_t /*c1*/, std::uint16_t /*c2*/,
                  std::uint32_t q ) const
    {
        ++scalarCalls;
        return coefficient[q];
    }

    Eigen::VectorXd coefficient;
    table_type testValues;
    table_type trialValues;
    mutable std::size_t scalarCalls = 0;
};

/**
 * @brief Tensor adding the G.7c Eigen quadrature-block capability.
 */
struct batch_tensor : scalar_tensor
{
    using scalar_tensor::scalar_tensor;

    /**
     * @brief Evaluate one complete point block with Eigen operations.
     * @param request value, test, or bilinear indices
     * @param values caller-owned output vector
     */
    void evaluateQuadratureBatch(
        quadrature::evaluation_request const& request,
        Eigen::Ref<Eigen::VectorXd> values ) const
    {
        ++batchCalls;
        switch ( request.kind )
        {
        case quadrature::evaluation_kind::value:
            values = coefficient;
            break;
        case quadrature::evaluation_kind::test:
            values.array() = coefficient.array() *
                             testValues.row( request.test ).transpose().array();
            break;
        case quadrature::evaluation_kind::bilinear:
            values.array() = coefficient.array() *
                             testValues.row( request.test ).transpose().array() *
                             trialValues.row( request.trial ).transpose().array();
            break;
        }
    }

    mutable std::size_t batchCalls = 0;
};

/** @brief Minimal scalar FE context for complete local-matrix checks. */
struct scalar_fe_batch_context
{
    using value_type = double;
    using table_type = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                     Eigen::RowMajor>;
    static constexpr uint16_type NDim = 2;
    static constexpr uint16_type nComponents1 = 1;

    scalar_fe_batch_context( Eigen::Index dofs, Eigen::Index points )
        : values( dofs, points ), gradients{ table_type( dofs, points ),
                                             table_type( dofs, points ) }
    {
        for ( Eigen::Index i = 0; i < dofs; ++i )
            for ( Eigen::Index q = 0; q < points; ++q )
            {
                values( i, q ) = 0.3 + 0.07*i + 0.02*q;
                gradients[0]( i, q ) = -0.2 + 0.03*i + 0.01*q;
                gradients[1]( i, q ) = 0.4 - 0.02*i + 0.015*q;
            }
    }

    /** @return number of local basis functions. */
    Eigen::Index nDofs() const noexcept { return values.rows(); }

    /** @return complete scalar value table. */
    table_type const& valueBatchTable() const noexcept { return values; }

    /** @return the scalar value component for every DoF and point. */
    table_type const& valueComponentBatchTable( uint16_type = 0 ) const noexcept
    {
        return values;
    }

    /** @return one basis-value quadrature row. */
    table_type::ConstRowXpr idBatch( Eigen::Index i,
                                    uint16_type = 0 ) const
    {
        return values.row( i );
    }

    /** @return one physical-gradient quadrature row. */
    table_type::ConstRowXpr gradientBatch( Eigen::Index i, uint16_type,
                                          uint16_type direction ) const
    {
        return gradients[direction].row( i );
    }

    /** @return one gradient direction for every DoF and point. */
    table_type const& gradientComponentBatchTable(
        uint16_type, uint16_type direction ) const
    {
        return gradients[direction];
    }

    /** @return one scalar basis value. */
    double id( Eigen::Index i, uint16_type, uint16_type,
               Eigen::Index q ) const
    {
        return values( i, q );
    }

    /** @return one physical-gradient component. */
    double grad( Eigen::Index i, uint16_type, uint16_type direction,
                 Eigen::Index q ) const
    {
        return gradients[direction]( i, q );
    }

    table_type values;
    std::array<table_type, NDim> gradients;
};

static_assert( !QuadratureBatchEvaluable<scalar_tensor, double> );
static_assert( QuadratureBatchEvaluable<batch_tensor, double> );

} // namespace

BOOST_AUTO_TEST_SUITE( quadrature_batch_suite )

/** Verify pointwise/batched parity and compile-time batch selection. */
BOOST_AUTO_TEST_CASE( eigen_batch_matches_scalar_fallback )
{
    PointSetQuadrature<Simplex<2>, double> quadratureRule( 6 );
    quadratureRule.update( unit_jacobian_context{} );
    auto const pointCount = static_cast<Eigen::Index>( quadratureRule.nPoints() );

    scalar_tensor scalar( pointCount );
    batch_tensor batch( pointCount );

    auto const scalarBilinear = quadratureRule( scalar, 1u, 2u, 0u, 0u );
    auto const batchBilinear = quadratureRule( batch, 1u, 2u, 0u, 0u );
    BOOST_CHECK_SMALL( batchBilinear - scalarBilinear, 1e-13 );
    BOOST_TEST( scalar.scalarCalls == quadratureRule.nPoints() );
    BOOST_TEST( batch.scalarCalls == 0u );
    BOOST_TEST( batch.batchCalls == 1u );

    scalar.scalarCalls = 0;
    batch.batchCalls = 0;
    auto const scalarLinear = quadratureRule( scalar, 2u, 0u, 0u );
    auto const batchLinear = quadratureRule( batch, 2u, 0u, 0u );
    BOOST_CHECK_SMALL( batchLinear - scalarLinear, 1e-13 );
    BOOST_TEST( scalar.scalarCalls == quadratureRule.nPoints() );
    BOOST_TEST( batch.batchCalls == 1u );

    scalar.scalarCalls = 0;
    batch.batchCalls = 0;
    auto const scalarValue = quadratureRule( scalar, 0u, 0u );
    auto const batchValue = quadratureRule( batch, 0u, 0u );
    BOOST_CHECK_SMALL( batchValue - scalarValue, 1e-13 );
    BOOST_TEST( scalar.scalarCalls == quadratureRule.nPoints() );
    BOOST_TEST( batch.batchCalls == 1u );
}

/** Verify mass and H1 strategies form every local matrix entry correctly. */
BOOST_AUTO_TEST_CASE( complete_local_tensor_matches_explicit_sums )
{
    constexpr Eigen::Index testDofs = 4;
    constexpr Eigen::Index trialDofs = 3;
    constexpr Eigen::Index points = 7;
    scalar_fe_batch_context test( testDofs, points );
    scalar_fe_batch_context trial( trialDofs, points );
    Eigen::VectorXd weightedCoefficient( points );
    for ( Eigen::Index q = 0; q < points; ++q )
        weightedCoefficient[q] = 0.15 + 0.04*q;

    Eigen::MatrixXd mass( testDofs, trialDofs );
    Eigen::MatrixXd stiffness( testDofs, trialDofs );
    Eigen::MatrixXd testScratch;
    Eigen::MatrixXd trialScratch;
    vf::detail::scalar_mass_contraction::contract(
        test, trial, weightedCoefficient, mass, testScratch, trialScratch );
    vf::detail::scalar_h1_stiffness_contraction::contract(
        test, trial, weightedCoefficient, stiffness,
        testScratch, trialScratch );

    for ( Eigen::Index i = 0; i < testDofs; ++i )
        for ( Eigen::Index j = 0; j < trialDofs; ++j )
        {
            double expectedMass = 0.;
            double expectedStiffness = 0.;
            for ( Eigen::Index q = 0; q < points; ++q )
            {
                expectedMass += weightedCoefficient[q]*
                                test.values( i, q )*trial.values( j, q );
                for ( uint16_type direction = 0;
                      direction < scalar_fe_batch_context::NDim; ++direction )
                    expectedStiffness += weightedCoefficient[q]*
                        test.gradients[direction]( i, q )*
                        trial.gradients[direction]( j, q );
            }
            BOOST_CHECK_SMALL( mass( i, j ) - expectedMass, 1e-13 );
            BOOST_CHECK_SMALL( stiffness( i, j ) - expectedStiffness, 1e-13 );
        }
}

/** Verify fixed-shape Eigen workspaces retain their steady-state storage. */
BOOST_AUTO_TEST_CASE( reusable_eigen_workspace_keeps_allocation )
{
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        matrix;
    vf::detail::resizeEigenWorkspace( matrix, 7, 11 );
    auto const* matrixStorage = matrix.data();
    vf::detail::resizeEigenWorkspace( matrix, 7, 11 );
    BOOST_CHECK_EQUAL( matrix.data(), matrixStorage );

    Eigen::VectorXd vector;
    vf::detail::resizeEigenVectorWorkspace( vector, 13 );
    auto const* vectorStorage = vector.data();
    vf::detail::resizeEigenVectorWorkspace( vector, 13 );
    BOOST_CHECK_EQUAL( vector.data(), vectorStorage );

    vf::detail::resizeEigenWorkspace( matrix, 5, 3 );
    vf::detail::resizeEigenVectorWorkspace( vector, 4 );
    BOOST_CHECK_EQUAL( matrix.rows(), 5 );
    BOOST_CHECK_EQUAL( matrix.cols(), 3 );
    BOOST_CHECK_EQUAL( vector.size(), 4 );
}

BOOST_AUTO_TEST_SUITE_END()
