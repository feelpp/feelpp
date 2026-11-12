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

#define BOOST_TEST_MODULE feelsystem
#include <feel/feelcore/testsuite.hpp>

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
    Feel::po::options_description options( "FeelSystem test options" );
    options.add( Feel::backend_options( "linear-smoke" ) )
           .add( Feel::backend_options( "linear-smoke-alt" ) )
           .add( Feel::backend_options( "natural-boundaries" ) );
    return options;
}

}

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAboutDefault( "feelpp_test_system" ), makeOptions() )

BOOST_AUTO_TEST_SUITE( feelsystem_suite )

BOOST_AUTO_TEST_CASE( linear_problem_records_and_assembles )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto uh = Xh->element( "u" );
    auto u = trial( Xh );
    auto v = test( Xh );
    auto backend = Feel::backend( _worldcomm = Xh->worldCommPtr() );

    auto lhsTerm = integrate( _range = elements( mesh ),
                              _expr = inner( grad( u ), grad( v ) ) + u * v );
    auto rhsTerm = integrate( _range = elements( mesh ),
                              _expr = v );

    auto problem =
        linearProblem( _space = Xh, _name = "linear-smoke" )
            .lhs( lhsTerm )
            .rhs( rhsTerm );

    static_assert( systems::LinearProblemLike<decltype( problem )> );
    BOOST_CHECK_EQUAL( problem.name(), "linear-smoke" );
    BOOST_CHECK_EQUAL( problem.lhsTerms().size(), 1 );
    BOOST_CHECK_EQUAL( problem.rhsTerms().size(), 1 );

    auto system = assemble( problem, _backend = backend );
    BOOST_CHECK( system.assembled() );
    BOOST_CHECK( system.matrixPtr() );
    BOOST_CHECK( system.rhsPtr() );
    BOOST_CHECK( system.backend() == backend );

    auto defaultSystem = assemble( problem );
    BOOST_CHECK_EQUAL( defaultSystem.backendName(), problem.name() );
    BOOST_CHECK_EQUAL( defaultSystem.backend()->prefix(), problem.name() );

    auto defaultSystemAgain = assemble( problem );
    BOOST_CHECK( defaultSystemAgain.backend() == defaultSystem.backend() );

    auto memberSystem = problem.assemble();
    BOOST_CHECK( memberSystem.assembled() );
    BOOST_CHECK_EQUAL( memberSystem.backendName(), problem.name() );
    BOOST_CHECK( memberSystem.backend() == defaultSystem.backend() );

    auto namedSystem = assemble( problem, _name = "linear-smoke-alt" );
    BOOST_CHECK_EQUAL( namedSystem.backendName(), "linear-smoke-alt" );
    BOOST_CHECK_EQUAL( namedSystem.backend()->prefix(), "linear-smoke-alt" );
    BOOST_CHECK( namedSystem.backend() != defaultSystem.backend() );

    auto memberNamedSystem = problem.assemble( _name = "linear-smoke-alt" );
    BOOST_CHECK_EQUAL( memberNamedSystem.backendName(), "linear-smoke-alt" );
    BOOST_CHECK( memberNamedSystem.backend() == namedSystem.backend() );

    auto explicitNamedSystem = assemble( problem, _backend = backend, _name = "ignored-by-explicit-backend" );
    BOOST_CHECK( explicitNamedSystem.backend() == backend );

    auto memberExplicitSystem = problem.assemble( _backend = backend );
    BOOST_CHECK( memberExplicitSystem.backend() == backend );

    auto directMatrix = backend->newMatrix( _test = Xh, _trial = Xh );
    auto directRhs = backend->newVector( _test = Xh );
    directMatrix->zero();
    directRhs->zero();
    form2( _test = Xh, _trial = Xh, _matrix = directMatrix ) += lhsTerm;
    form1( _test = Xh, _vector = directRhs ) += rhsTerm;
    directMatrix->close();
    directRhs->close();

    auto matrixDelta = system.matrixPtr()->clone();
    matrixDelta->addMatrix( -1.0, directMatrix );
    matrixDelta->close();
    BOOST_CHECK_SMALL( matrixDelta->l1Norm(), 1e-9 );
    BOOST_CHECK_SMALL( matrixDelta->linftyNorm(), 1e-9 );

    auto rhsDelta = vf::cloneVectorWithValues( system.rhsPtr() );
    rhsDelta->add( -1.0, *directRhs );
    rhsDelta->close();
    BOOST_CHECK_SMALL( rhsDelta->linftyNorm(), 1e-10 );

    system.assemble();
    matrixDelta = system.matrixPtr()->clone();
    matrixDelta->addMatrix( -1.0, directMatrix );
    matrixDelta->close();
    BOOST_CHECK_SMALL( matrixDelta->l1Norm(), 1e-9 );

    auto activeSystem = system.activeSystem();
    BOOST_CHECK( activeSystem.first == system.matrixPtr() );
    BOOST_CHECK( activeSystem.second == system.rhsPtr() );

    uh.zero();
    auto result = systems::solve( system, _solution = uh );
    BOOST_CHECK( result.converged );
    BOOST_CHECK_GE( result.iterations, 0 );

    auto residual = backend->newVector( _test = Xh );
    system.matrixPtr()->multVector( uh, *residual );
    residual->add( -1.0, *system.rhsPtr() );
    residual->close();
    BOOST_CHECK_SMALL( residual->linftyNorm(), 1e-8 );

    using space_type = typename std::remove_cvref_t<decltype( Xh )>::element_type;
    systems::Problem<space_type, systems::linear_t, systems::stationary_t> explicitProblem( Xh, "explicit-linear" );
    explicitProblem.add( systems::lhs, lhsTerm )
                   .add( systems::rhs, rhsTerm );

    auto explicitSystem = assemble( explicitProblem, _backend = backend );
    BOOST_CHECK( explicitSystem.assembled() );
    BOOST_CHECK_EQUAL( explicitProblem.lhsTerms().size(), 1 );
    BOOST_CHECK_EQUAL( explicitProblem.rhsTerms().size(), 1 );

    auto u2 = Xh->element( "u2" );
    auto oneShotResult = systems::solve( problem, _solution = u2 );
    BOOST_CHECK( oneShotResult.converged );

    auto u3 = Xh->element( "u3" );
    auto memberOneShotResult = problem.solve( _solution = u3 );
    BOOST_CHECK( memberOneShotResult.converged );
}

BOOST_AUTO_TEST_CASE( natural_boundary_sugar_matches_direct_forms )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto u = trial( Xh );
    auto v = test( Xh );

    auto lhsTerm = integrate( _range = elements( mesh ),
                              _expr = inner( grad( u ), grad( v ) ) );
    auto rhsTerm = integrate( _range = elements( mesh ),
                              _expr = v );
    auto neumannFlux = cst( 2.0 );
    auto robinAlpha = cst( 3.0 );
    auto robinValue = cst( 4.0 );
    auto fourierCoeff = cst( 5.0 );
    auto exteriorValue = cst( 6.0 );

    auto problem =
        linearProblem( _space = Xh, _name = "natural-boundaries" )
            .lhs( lhsTerm )
            .rhs( rhsTerm )
            .neumann( _range = boundaryfaces( mesh ),
                      _test = v,
                      _expr = neumannFlux )
            .robin( _range = boundaryfaces( mesh ),
                    _trial = u,
                    _test = v,
                    _alpha = robinAlpha,
                    _expr = robinValue )
            .fourier( _range = boundaryfaces( mesh ),
                      _trial = u,
                      _test = v,
                      _alpha = fourierCoeff,
                      _expr = exteriorValue );

    BOOST_CHECK_EQUAL( problem.lhsTerms().size(), 3 );
    BOOST_CHECK_EQUAL( problem.rhsTerms().size(), 4 );

    auto system = assemble( problem );
    BOOST_CHECK( !system.constrained() );
    BOOST_CHECK_EQUAL( system.backendName(), problem.name() );

    auto backend = system.backend();
    auto directMatrix = backend->newMatrix( _test = Xh, _trial = Xh );
    auto directRhs = backend->newVector( _test = Xh );
    directMatrix->zero();
    directRhs->zero();
    form2( _test = Xh, _trial = Xh, _matrix = directMatrix ) += lhsTerm;
    form1( _test = Xh, _vector = directRhs ) += rhsTerm;
    form1( _test = Xh, _vector = directRhs ) +=
        integrate( _range = boundaryfaces( mesh ),
                   _expr = neumannFlux * v );
    form2( _test = Xh, _trial = Xh, _matrix = directMatrix ) +=
        integrate( _range = boundaryfaces( mesh ),
                   _expr = robinAlpha * u * v );
    form1( _test = Xh, _vector = directRhs ) +=
        integrate( _range = boundaryfaces( mesh ),
                   _expr = robinValue * v );
    form2( _test = Xh, _trial = Xh, _matrix = directMatrix ) +=
        integrate( _range = boundaryfaces( mesh ),
                   _expr = fourierCoeff * u * v );
    form1( _test = Xh, _vector = directRhs ) +=
        integrate( _range = boundaryfaces( mesh ),
                   _expr = fourierCoeff * exteriorValue * v );
    directMatrix->close();
    directRhs->close();

    auto matrixDelta = system.matrixPtr()->clone();
    matrixDelta->addMatrix( -1.0, directMatrix );
    matrixDelta->close();
    BOOST_CHECK_SMALL( matrixDelta->linftyNorm(), 1e-9 );

    auto rhsDelta = vf::cloneVectorWithValues( system.rhsPtr() );
    rhsDelta->add( -1.0, *directRhs );
    rhsDelta->close();
    BOOST_CHECK_SMALL( rhsDelta->linftyNorm(), 1e-10 );
}

BOOST_AUTO_TEST_SUITE_END()
