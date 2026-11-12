/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 * @file
 * @brief FeelSystem regression test.
 */

#define BOOST_TEST_MODULE feelsystem_lifetime
#include <feel/feelcore/testsuite.hpp>

#include <cstddef>
#include <functional>
#include <memory>
#include <stdexcept>
#include <utility>

#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelsystem/system.hpp>
#include <feel/feelvf/dirichletconstraints.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/options.hpp>

namespace
{

Feel::po::options_description
makeOptions()
{
    Feel::po::options_description options( "FeelSystem lifetime test options" );
    options.add( Feel::backend_options( "lifetime-linear-temporary" ) )
        .add( Feel::backend_options( "lifetime-boundary" ) )
        .add( Feel::backend_options( "lifetime-evolution" ) )
        .add( Feel::backend_options( "lifetime-nonlinear" ) )
        .add( Feel::backend_options( "lifetime-picard" ) )
        .add( Feel::bdf_options( "lifetime-time" ) )
        .add( Feel::ts_options( "lifetime-time" ) );
    return options;
}

template <typename ElementType>
void setConstantClosed( ElementType& element, double value )
{
    element.setConstant( value );
    element.close();
}

template <typename MatrixPtr>
void checkMatrixClose( MatrixPtr const& actual, MatrixPtr const& expected )
{
    auto delta = actual->clone();
    delta->addMatrix( -1.0, expected );
    delta->close();
    BOOST_CHECK_SMALL( delta->l1Norm(), 1e-9 );
    BOOST_CHECK_SMALL( delta->linftyNorm(), 1e-9 );
}

template <typename VectorPtr>
void checkVectorClose( VectorPtr const& actual, VectorPtr const& expected )
{
    auto delta = Feel::vf::cloneVectorWithValues( actual );
    delta->add( -1.0, *expected );
    delta->close();
    BOOST_CHECK_SMALL( delta->linftyNorm(), 1e-10 );
}

} // namespace

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAboutDefault( "feelpp_test_system_lifetime" ), makeOptions() )

BOOST_AUTO_TEST_SUITE( feelsystem_lifetime_suite )

BOOST_AUTO_TEST_CASE( nonlinear_callback_guard_restores_callbacks_on_scope_exit_and_exception )
{
    struct fake_solver
    {
        std::function<void()> residual;
        std::function<void()> jacobian;
    };

    fake_solver solver;
    int residualOld = 0;
    int jacobianOld = 0;
    int residualNew = 0;
    int jacobianNew = 0;

    solver.residual = [&residualOld]()
    { ++residualOld; };
    solver.jacobian = [&jacobianOld]()
    { ++jacobianOld; };

    {
        Feel::systems::detail::NonlinearCallbackGuard guard( solver );
        solver.residual = [&residualNew]()
        { ++residualNew; };
        solver.jacobian = [&jacobianNew]()
        { ++jacobianNew; };
    }

    solver.residual();
    solver.jacobian();
    BOOST_CHECK_EQUAL( residualOld, 1 );
    BOOST_CHECK_EQUAL( jacobianOld, 1 );
    BOOST_CHECK_EQUAL( residualNew, 0 );
    BOOST_CHECK_EQUAL( jacobianNew, 0 );

    try
    {
        Feel::systems::detail::NonlinearCallbackGuard guard( solver );
        solver.residual = [&residualNew]()
        { ++residualNew; };
        solver.jacobian = [&jacobianNew]()
        { ++jacobianNew; };
        throw std::runtime_error( "expected callback guard test exception" );
    }
    catch ( std::runtime_error const& )
    {
    }

    solver.residual();
    solver.jacobian();
    BOOST_CHECK_EQUAL( residualOld, 2 );
    BOOST_CHECK_EQUAL( jacobianOld, 2 );
    BOOST_CHECK_EQUAL( residualNew, 0 );
    BOOST_CHECK_EQUAL( jacobianNew, 0 );
}

BOOST_AUTO_TEST_CASE( temporary_forms_and_constraint_callback_solve_manufactured_reaction_diffusion )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto uh = Xh->element( "u" );
    auto u = trial( Xh );
    auto v = test( Xh );
    auto constraintCalls = std::make_shared<std::size_t>( 0 );

    auto makeProblem = [&]()
    {
        auto exact = cst( 1.0 ) + Px() + Py();
        auto lhsTerm = integrate( _range = elements( mesh ),
                                  _expr = inner( grad( u ), grad( v ) ) + u * v );
        auto rhsTerm = integrate( _range = elements( mesh ),
                                  _expr = exact * v );

        return linearProblem( _space = Xh, _name = "lifetime-linear-temporary" )
            .lhs( lhsTerm )
            .rhs( rhsTerm )
            .constraint( [mesh, uh, exact, constraintCalls]( auto& a, auto& l )
                         {
                ++*constraintCalls;
                a += on( _range = boundaryfaces( mesh ),
                         _rhs = l,
                         _element = uh,
                         _expr = exact ); } );
    };

    auto problem = makeProblem();
    auto copiedProblem = problem;
    auto system = copiedProblem.assemble();
    BOOST_CHECK_EQUAL( *constraintCalls, 1 );

    uh.zero();
    auto result = system.solve( _solution = uh );
    BOOST_CHECK( result.converged );

    auto exact = cst( 1.0 ) + Px() + Py();
    BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                               _expr = uh - exact ),
                       1e-8 );

    auto u2 = Xh->element( "u2" );
    u2.zero();
    system.assemble();
    BOOST_CHECK_EQUAL( *constraintCalls, 2 );
    result = system.solve( _solution = u2 );
    BOOST_CHECK( result.converged );
    BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                               _expr = u2 - exact ),
                       1e-8 );
}

BOOST_AUTO_TEST_CASE( boundary_sugar_and_parameter_updates_solve_manufactured_natural_problem )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto uh = Xh->element( "u" );
    auto u = trial( Xh );
    auto v = test( Xh );

    auto makeProblem = [&]()
    {
        auto exact = cst( 1.0 ) + Px() + Py();
        auto normalFlux = Nx() + Ny();
        auto alpha = expr( "alpha:alpha" );
        auto beta = expr( "beta:beta" );

        return linearProblem( _space = Xh, _name = "lifetime-boundary" )
            .setParameterValues( { { "alpha", 2.0 }, { "beta", 3.0 } } )
            .diffusion( _range = elements( mesh ),
                        _expr = inner( grad( u ), grad( v ) ) )
            .reaction( _range = elements( mesh ),
                       _expr = u * v )
            .source( _range = elements( mesh ),
                     _expr = exact * v )
            .neumann( _range = boundaryfaces( mesh ),
                      _test = v,
                      _expr = normalFlux )
            .robin( _range = boundaryfaces( mesh ),
                    _trial = u,
                    _test = v,
                    _alpha = alpha,
                    _expr = alpha * exact )
            .fourier( _range = boundaryfaces( mesh ),
                      _trial = u,
                      _test = v,
                      _alpha = beta,
                      _expr = exact );
    };

    auto problem = makeProblem();
    auto system = problem.assemble();
    BOOST_CHECK_EQUAL( system.backendName(), "lifetime-boundary" );

    auto exact = cst( 1.0 ) + Px() + Py();
    auto normalFlux = Nx() + Ny();
    for ( auto const parameters : { std::pair{ 2.0, 3.0 }, std::pair{ 5.0, 7.0 } } )
    {
        problem.setParameterValue( "alpha", parameters.first )
            .setParameterValue( "beta", parameters.second );

        uh.zero();
        system.assemble();

        auto directMatrix = system.backend()->newMatrix( _test = Xh, _trial = Xh );
        auto directRhs = system.backend()->newVector( _test = Xh );
        directMatrix->zero();
        directRhs->zero();
        form2( _test = Xh, _trial = Xh, _matrix = directMatrix ) +=
            integrate( _range = elements( mesh ),
                       _expr = inner( grad( u ), grad( v ) ) + u * v );
        form1( _test = Xh, _vector = directRhs ) +=
            integrate( _range = elements( mesh ),
                       _expr = exact * v );
        form1( _test = Xh, _vector = directRhs ) +=
            integrate( _range = boundaryfaces( mesh ),
                       _expr = normalFlux * v );
        form2( _test = Xh, _trial = Xh, _matrix = directMatrix ) +=
            integrate( _range = boundaryfaces( mesh ),
                       _expr = cst( parameters.first ) * u * v );
        form1( _test = Xh, _vector = directRhs ) +=
            integrate( _range = boundaryfaces( mesh ),
                       _expr = cst( parameters.first ) * exact * v );
        form2( _test = Xh, _trial = Xh, _matrix = directMatrix ) +=
            integrate( _range = boundaryfaces( mesh ),
                       _expr = cst( parameters.second ) * u * v );
        form1( _test = Xh, _vector = directRhs ) +=
            integrate( _range = boundaryfaces( mesh ),
                       _expr = cst( parameters.second ) * exact * v );
        directMatrix->close();
        directRhs->close();
        checkMatrixClose( system.matrixPtr(), directMatrix );
        checkVectorClose( system.rhsPtr(), directRhs );

        auto result = system.solve( _solution = uh );
        BOOST_CHECK( result.converged );
        BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                                   _expr = uh - exact ),
                           1e-8 );
    }
}

BOOST_AUTO_TEST_CASE( bdf_time_parameter_terms_solve_manufactured_linear_evolution )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto uh = Xh->element( "u" );
    auto u = trial( Xh );
    auto v = test( Xh );

    auto makeProblem = [&]()
    {
        auto source = expr( "time/time:time" );
        auto timeScheme = bdf( _space = Xh,
                               _name = "lifetime-time",
                               _prefix = "lifetime-time",
                               _order = 1,
                               _initial_time = 0.0,
                               _final_time = 0.2,
                               _time_step = 0.1,
                               _save = false );

        return parabolicProblem( _space = Xh, _name = "lifetime-evolution" )
            .timeSymbol( "time" )
            .timeScheme( timeScheme )
            .mass( _range = elements( mesh ),
                   _expr = u * v )
            .rhs( _range = elements( mesh ),
                  _expr = source * v );
    };

    auto problem = makeProblem();
    auto evolution = problem.evolve();

    auto exact = Xh->element( "exact" );
    auto u0 = Xh->element( "u0" );
    setConstantClosed( u0, 1.0 );
    evolution.start( u0 );

    for ( int step = 0; step != 2; ++step )
    {
        evolution.assemble();
        BOOST_CHECK_CLOSE( evolution.problem().parameterValues().at( "time" ), evolution.time(), 1e-10 );

        uh.zero();
        auto result = evolution.solve( _solution = uh );
        BOOST_CHECK( result.converged );

        exact.on( _range = elements( mesh ),
                  _expr = cst( 1.0 + evolution.time() ) );
        exact.close();
        BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                                   _expr = uh - exact ),
                           1e-8 );

        evolution.advance( uh );
    }
}

BOOST_AUTO_TEST_CASE( nonlinear_callbacks_solve_manufactured_constant_root_after_scope_exit )
{
    using namespace Feel;
    using namespace Feel::vf;

    constexpr double beta = 0.1;
    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto uh = Xh->element( "u" );
    auto u = trial( Xh );
    auto v = test( Xh );
    auto residualCalls = std::make_shared<std::size_t>( 0 );
    auto jacobianCalls = std::make_shared<std::size_t>( 0 );

    auto makeProblem = [&]()
    {
        return nonlinearProblem( _space = Xh, _name = "lifetime-nonlinear" )
            .reaction( _range = elements( mesh ),
                       _expr = u * v )
            .source( _range = elements( mesh ),
                     _expr = cst( 1.0 + beta ) * v )
            .residual( [mesh, v, residualCalls]( auto const& state, auto& residual )
                       {
                ++*residualCalls;
                residual += integrate( _range = elements( mesh ),
                                       _expr = beta * pow( idv( state ), 3 ) * v ); } )
            .jacobian( [mesh, u, v, jacobianCalls]( auto const& state, auto& jacobian )
                       {
                ++*jacobianCalls;
                jacobian += integrate( _range = elements( mesh ),
                                      _expr = 3.0 * beta * pow( idv( state ), 2 ) * u * v ); } );
    };

    auto problem = makeProblem();
    auto copiedProblem = problem;
    setConstantClosed( uh, 0.5 );

    auto result = copiedProblem.solve( _solution = uh,
                                       _maxit = 25,
                                       _rtolerance = 1e-12 );
    BOOST_CHECK( result.converged );
    BOOST_CHECK_GE( *residualCalls, std::size_t{ 1 } );
    BOOST_CHECK_GE( *jacobianCalls, std::size_t{ 1 } );
    BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                               _expr = uh - cst( 1.0 ) ),
                       1e-8 );
}

BOOST_AUTO_TEST_CASE( picard_callback_solve_manufactured_constant_root_after_problem_copy )
{
    using namespace Feel;
    using namespace Feel::vf;

    constexpr double beta = 0.1;
    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto uh = Xh->element( "u" );
    auto u = trial( Xh );
    auto v = test( Xh );
    auto iterationCalls = std::make_shared<std::size_t>( 0 );

    auto makeProblem = [&]()
    {
        return picardProblem( _space = Xh, _name = "lifetime-picard" )
            .iteration( [mesh, u, v, iterationCalls]( auto const& state, auto& a, auto& l )
                        {
                ++*iterationCalls;
                a += integrate( _range = elements( mesh ),
                                _expr = ( 1.0 + beta * pow( idv( state ), 2 ) ) *
                                        u * v );
                l += integrate( _range = elements( mesh ),
                                _expr = cst( 1.0 + beta ) * v ); } );
    };

    auto problem = makeProblem();
    auto copiedProblem = problem;
    setConstantClosed( uh, 0.5 );

    auto result = copiedProblem.solve( _solution = uh,
                                       _strategy = systems::picard( _maxit = 100,
                                                                    _rtolerance = 1e-12 ) );
    BOOST_CHECK( result.converged );
    BOOST_CHECK_GT( *iterationCalls, std::size_t{ 1 } );
    BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                               _expr = uh - cst( 1.0 ) ),
                       1e-8 );
}

BOOST_AUTO_TEST_SUITE_END()
