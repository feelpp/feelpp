/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/** @file test_polynomial_solve.cpp
 * @brief Checked interpolation solves, reusable factors and scalar-range behavior.
 */
#define BOOST_TEST_MODULE polynomial_solve
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelpoly/polynomialsolve.hpp>
#include <boost/mp11/list.hpp>
#include <boost/test/data/test_case.hpp>
#include <array>
#include <limits>

using namespace Feel;
FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
/** Scalar types exercise one shared solve contract. */
using Scalars = boost::mp11::mp_list<float, double>;

/** @brief Verify left and right solves using explicit expected polynomial coefficients. */
template<typename T>
void checkScaledSolve( T scale )
{
    using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    Matrix system( 2, 2 );
    system << T( 4 ), T( 1 ), T( 2 ), T( 3 );
    system *= scale;
    Matrix expected( 2, 3 );
    expected << T( .1 ), T( 0 ), T( -2 ), T( .6 ), T( 1 ), T( 3 );
    Matrix rhs = system * expected;
    ::Feel::detail::PolynomialLuSolver<T> solver( system, "test interpolation" );
    T const tolerance = T( 100 ) * Eigen::NumTraits<T>::epsilon();
    BOOST_CHECK_SMALL( ( solver.solve( rhs ) - expected ).norm(), tolerance );
    Matrix rightRhs = expected.transpose() * system;
    BOOST_CHECK_SMALL( ( solver.solveRight( rightRhs ) - expected.transpose() ).norm(), tolerance );
    BOOST_CHECK_SMALL( solver.solve( Matrix::Zero( 2, 1 ) ).norm(), tolerance );
    // The factorization owns its input; changing the source cannot change a solve.
    system.setZero();
    BOOST_CHECK_SMALL( ( solver.solve( rhs ) - expected ).norm(), tolerance );
    BOOST_CHECK_THROW( (void)solver.solve( Matrix::Zero( 3, 1 ) ), std::invalid_argument );
    BOOST_CHECK_THROW( (void)solver.solveRight( Matrix::Zero( 1, 3 ) ), std::invalid_argument );
    BOOST_CHECK_EQUAL( solver.solve( Matrix::Zero( 2, 0 ) ).cols(), 0 );
    BOOST_CHECK_EQUAL( solver.solveRight( Matrix::Zero( 0, 2 ) ).rows(), 0 );
}
} // namespace

/** @test Reused ordinary/transposed factors preserve solutions across uniform scales. */
BOOST_AUTO_TEST_CASE_TEMPLATE( scaledInterpolationSolves, T, Scalars )
{
    for ( T scale : std::array<T, 3>{ T( 1e-20 ), T( 1 ), T( 1e20 ) } )
        BOOST_TEST_CONTEXT( "scale=" << scale ) { checkScaledSolve( scale ); }
}

/** @test Finite double-precision extremes do not overflow residual norms. */
BOOST_DATA_TEST_CASE( extremeInterpolationScales,
                     boost::unit_test::data::make( std::array<double, 2>{ 1e-300, 1e300 } ), scale )
{
    checkScaledSolve( scale );
}

/** @test Invalid systems and right-hand sides fail in Release as well as debug builds. */
BOOST_AUTO_TEST_CASE( invalidInterpolationSystems )
{
    using Solver = ::Feel::detail::PolynomialLuSolver<double>;
    using Matrix = Solver::Matrix;
    BOOST_CHECK_THROW( (void)Solver( Matrix::Zero( 2, 3 ), "rectangular" ), std::invalid_argument );
    BOOST_CHECK_THROW( (void)Solver( Matrix::Zero( 2, 2 ), "zero" ), std::runtime_error );
    Matrix singular( 2, 2 );
    singular << 1., 2., 2., 4.;
    BOOST_CHECK_THROW( (void)Solver( singular, "singular" ), std::runtime_error );
    Matrix invalid = Matrix::Identity( 2, 2 );
    invalid( 0, 0 ) = std::numeric_limits<double>::infinity();
    BOOST_CHECK_THROW( (void)Solver( invalid, "infinite" ), std::invalid_argument );
    Solver solver( Matrix::Identity( 2, 2 ), "finite" );
    invalid( 0, 0 ) = std::numeric_limits<double>::quiet_NaN();
    BOOST_CHECK_THROW( (void)solver.solve( invalid ), std::invalid_argument );
    Solver empty( Matrix( 0, 0 ), "empty" );
    BOOST_CHECK_EQUAL( empty.solve( Matrix( 0, 3 ) ).cols(), 3 );
    BOOST_CHECK_EQUAL( empty.solveRight( Matrix( 3, 0 ) ).rows(), 3 );
}
