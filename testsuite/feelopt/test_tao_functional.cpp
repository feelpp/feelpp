/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

#define BOOST_TEST_MODULE tao functional testsuite

#include <cmath>
#include <vector>

#include <feel/feelcore/testsuite.hpp>
#include <feel/feelopt/solveroptimizationfield.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelvf/vf.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace Feel
{

/**
 * Check and report a deterministic TAO convergence history.
 *
 * Every MPI size executes these same checks, which verifies that global
 * objective and gradient residuals are independent of the partition.
 *
 * @tparam MonitorRecord structured TAO monitor-record type
 * @tparam Result structured optimization-result type
 * @param history records collected by the C++ TAO monitor
 * @param expectedInitialObjective analytic objective at the zero initial field
 * @param result final structured solver result
 */
template<typename MonitorRecord, typename Result>
void checkConvergenceHistory( std::vector<MonitorRecord> const& history,
                              double expectedInitialObjective,
                              Result const& result )
{
    BOOST_REQUIRE( history.size() >= 2 );
    BOOST_CHECK_EQUAL( history.front().iteration, 0 );
    BOOST_CHECK_CLOSE_FRACTION( history.front().objective,
                                expectedInitialObjective, 1e-10 );

    for ( std::size_t index = 0; index < history.size(); ++index )
    {
        auto const& record = history[index];
        BOOST_TEST_MESSAGE( "TAO iteration=" << record.iteration
                            << ", objective=" << record.objective
                            << ", gradient residual=" << record.gradientNorm
                            << ", step=" << record.stepNorm );
        BOOST_CHECK( std::isfinite( record.objective ) );
        BOOST_CHECK( std::isfinite( record.gradientNorm ) );
        if ( index > 0 )
        {
            auto const& previous = history[index - 1];
            BOOST_CHECK_GE( record.iteration, previous.iteration );
            BOOST_CHECK_LE( record.objective, previous.objective + 1e-12 );
            BOOST_CHECK_LE( record.gradientNorm,
                            previous.gradientNorm + 1e-12 );
        }
    }

    BOOST_CHECK_SMALL( std::abs( history.back().objective - result.objective ),
                       1e-12 );
    BOOST_CHECK_SMALL(
        std::abs( history.back().gradientNorm - result.gradientNorm ), 1e-12 );
}

BOOST_AUTO_TEST_SUITE( tao_functional_suite )

/** Optimize a finite-element field with Feel++ variational forms. */
BOOST_AUTO_TEST_CASE( variational_function_space_objective_and_hessian )
{
    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto testFunction = Xh->element( "v" );
    auto trialFunction = Xh->element( "u" );
    auto target = Xh->element( expr( "1+x+y:x:y" ), "target" );
    auto variable = Xh->element( "control" );

    auto optimization = optimizationField(
        _space = Xh, _backend = "petsc", _name = "tao_variational" );
    using monitor_record_type =
        typename decltype( optimization )::monitor_record_type;
    std::vector<monitor_record_type> history;
    auto const monitorId = optimization.addMonitor(
        [&history]( monitor_record_type const& record )
        {
            history.push_back( record );
        } );
    static_cast<void>( monitorId );
    optimization
        .algorithm( "nls" )
        .maxIterations( 20 )
        .gradientTolerance( 1e-10 )
        .objectiveGradient(
            [mesh, target, testFunction]( auto const& field,
                                          auto& gradient )
        {
            gradient = integrate(
                _range = elements( mesh ),
                _expr = ( idv( field ) - idv( target ) ) * id( testFunction ) );

            auto const difference = idv( field ) - idv( target );
            return 0.5 * integrate( _range = elements( mesh ),
                                    _expr = difference * difference )
                             .evaluate()( 0, 0 );
        } )
        .hessian(
            [mesh, trialFunction, testFunction]( auto const&, auto& hessian )
            {
                hessian = integrate(
                    _range = elements( mesh ),
                    _expr = idt( trialFunction ) * id( testFunction ) );
            } );

    auto const result = optimization.solve( variable );
    auto const error = normL2(
        _range = elements( mesh ),
        _expr = idv( variable ) - idv( target ) );

    checkConvergenceHistory( history, 25.0 / 12.0, result );

    BOOST_CHECK( result.converged );
    BOOST_CHECK_EQUAL( result.solverType, "nls" );
    BOOST_CHECK_SMALL( result.objective, 1e-18 );
    BOOST_CHECK_SMALL( result.gradientNorm, 1e-8 );
    BOOST_CHECK_SMALL( error, 1e-9 );
}

/** Optimize a product-space field through blockform1 and blockform2 syntax. */
BOOST_AUTO_TEST_CASE( variational_product_space_objective_and_hessian )
{
    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto productSpaceValue = product( Xh, Xh );
    auto Wh = std::make_shared<decltype( productSpaceValue )>(
        productSpaceValue );

    auto test0 = Xh->element( "v0" );
    auto test1 = Xh->element( "v1" );
    auto trial0 = Xh->element( "u0" );
    auto trial1 = Xh->element( "u1" );
    auto target0 = Xh->element( expr( "1+x+y:x:y" ), "target0" );
    auto target1 = Xh->element( expr( "2-x+y:x:y" ), "target1" );
    auto variable = Wh->element();

    auto optimization = optimizationField(
        _space = Wh, _backend = "petsc", _name = "tao_product_variational" );
    using monitor_record_type =
        typename decltype( optimization )::monitor_record_type;
    std::vector<monitor_record_type> history;
    auto const monitorId = optimization.addMonitor(
        [&history]( monitor_record_type const& record )
        {
            history.push_back( record );
        } );
    static_cast<void>( monitorId );
    optimization
        .algorithm( "nls" )
        .maxIterations( 20 )
        .gradientTolerance( 1e-10 )
        .objectiveGradient(
            [mesh, target0, target1, test0, test1]( auto const& field,
                                                    auto& gradient )
            {
                auto const difference0 = idv( field( 0_c ) ) - idv( target0 );
                auto const difference1 = idv( field( 1_c ) ) - idv( target1 );

                gradient( 0_c ) = integrate(
                    _range = elements( mesh ),
                    _expr = difference0 * id( test0 ) );
                gradient( 1_c ) += integrate(
                    _range = elements( mesh ),
                    _expr = difference1 * id( test1 ) );

                return 0.5 * integrate(
                    _range = elements( mesh ),
                    _expr = difference0 * difference0 +
                            difference1 * difference1 )
                                 .evaluate()( 0, 0 );
            } )
        .hessian(
            [mesh, trial0, trial1, test0, test1]( auto const&,
                                                  auto& hessian )
            {
                hessian( 0_c, 0_c ) = integrate(
                    _range = elements( mesh ),
                    _expr = idt( trial0 ) * id( test0 ) );
                hessian( 1_c, 1_c ) += integrate(
                    _range = elements( mesh ),
                    _expr = idt( trial1 ) * id( test1 ) );
            } );

    auto const result = optimization.solve( variable );
    auto const error0 = normL2(
        _range = elements( mesh ),
        _expr = idv( variable( 0_c ) ) - idv( target0 ) );
    auto const error1 = normL2(
        _range = elements( mesh ),
        _expr = idv( variable( 1_c ) ) - idv( target1 ) );

    checkConvergenceHistory( history, 25.0 / 6.0, result );

    BOOST_CHECK( result.converged );
    BOOST_CHECK_EQUAL( result.solverType, "nls" );
    BOOST_CHECK_SMALL( result.objective, 1e-18 );
    BOOST_CHECK_SMALL( result.gradientNorm, 1e-8 );
    BOOST_CHECK_SMALL( error0, 1e-9 );
    BOOST_CHECK_SMALL( error1, 1e-9 );
}

BOOST_AUTO_TEST_SUITE_END()

} // namespace Feel
