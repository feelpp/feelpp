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

#define BOOST_TEST_MODULE feelsystem_constraints
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
    Feel::po::options_description options( "FeelSystem constraint test options" );
    options.add( Feel::backend_options( "dirichlet-smoke" ) )
           .add( Feel::backend_options( "constraint-builder" ) );
    return options;
}

}

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAboutDefault( "feelpp_test_system_constraints" ), makeOptions() )

BOOST_AUTO_TEST_SUITE( feelsystem_constraints_suite )

BOOST_AUTO_TEST_CASE( dirichlet_constraint_matches_direct_forms )
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
        linearProblem( _space = Xh, _name = "dirichlet-smoke" )
            .lhs( lhsTerm )
            .rhs( rhsTerm )
            .dirichlet( _range = boundaryfaces( mesh ),
                        _element = uh,
                        _expr = cst( 0.0 ) );

    BOOST_CHECK_EQUAL( problem.lhsTerms().size(), 1 );
    BOOST_CHECK_EQUAL( problem.rhsTerms().size(), 1 );
    BOOST_CHECK_EQUAL( problem.constraintTerms().size(), 1 );
    BOOST_CHECK_EQUAL( problem.bcTerms().size(), 1 );

    auto system = assemble( problem, _backend = backend );
    BOOST_CHECK( system.matrixPtr() );
    BOOST_CHECK( system.rhsPtr() );
    BOOST_CHECK( system.activeMatrixPtr() );
    BOOST_CHECK( system.activeVectorPtr() );
    BOOST_CHECK( system.constrained() );
    BOOST_CHECK( system.activeMatrixPtr() == system.matrixPtr() );
    BOOST_CHECK( system.activeVectorPtr() != system.rhsPtr() );

    auto memberSystem = problem.assemble( _backend = backend );
    BOOST_CHECK( memberSystem.matrixPtr() );
    BOOST_CHECK( memberSystem.rhsPtr() );
    BOOST_CHECK( memberSystem.activeMatrixPtr() );
    BOOST_CHECK( memberSystem.activeVectorPtr() );
    BOOST_CHECK( memberSystem.constrained() );
    BOOST_CHECK( memberSystem.activeMatrixPtr() == memberSystem.matrixPtr() );
    BOOST_CHECK( memberSystem.activeVectorPtr() != memberSystem.rhsPtr() );
    BOOST_CHECK( memberSystem.backend() == backend );

    auto directMatrix = backend->newMatrix( _test = Xh, _trial = Xh );
    auto directRhs = backend->newVector( _test = Xh );
    directMatrix->zero();
    directRhs->zero();
    auto directA = form2( _test = Xh, _trial = Xh, _matrix = directMatrix );
    auto directF = form1( _test = Xh, _vector = directRhs );
    directA += lhsTerm;
    directF += rhsTerm;
    directA += on( _range = boundaryfaces( mesh ),
                   _rhs = directF,
                   _element = uh,
                   _expr = cst( 0.0 ) );
    directMatrix->close();
    directRhs->close();
    auto directActive = directA.activeSystem( directF );

    auto baseDelta = system.matrixPtr()->clone();
    baseDelta->addMatrix( -1.0, directMatrix );
    baseDelta->close();
    BOOST_CHECK_SMALL( baseDelta->l1Norm(), 1e-9 );

    auto activeDelta = system.activeMatrixPtr()->clone();
    activeDelta->addMatrix( -1.0, directActive.first );
    activeDelta->close();
    BOOST_CHECK_SMALL( activeDelta->l1Norm(), 1e-9 );
    BOOST_CHECK_SMALL( activeDelta->linftyNorm(), 1e-9 );

    auto rhsDelta = vf::cloneVectorWithValues( system.activeVectorPtr() );
    rhsDelta->add( -1.0, *directActive.second );
    rhsDelta->close();
    BOOST_CHECK_SMALL( rhsDelta->linftyNorm(), 1e-10 );

    auto result = systems::solve( system, _solution = uh );
    BOOST_CHECK( result.converged );

    auto boundaryError = Xh->element( "boundary-error" );
    boundaryError.on( _range = boundaryfaces( mesh ), _expr = abs( uh ) );
    BOOST_CHECK_SMALL( boundaryError.max(), 1e-8 );
}

BOOST_AUTO_TEST_CASE( constraint_builder_supports_existing_on_api )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto uh = Xh->element( "u" );
    auto u = trial( Xh );
    auto v = test( Xh );

    auto problem =
        linearProblem( _space = Xh, _name = "constraint-builder" )
            .lhs( integrate( _range = elements( mesh ),
                             _expr = inner( grad( u ), grad( v ) ) + u * v ) )
            .rhs( integrate( _range = elements( mesh ),
                             _expr = v ) )
            .constraint( [mesh, uh]( auto& a, auto& l )
            {
                a += on( _range = boundaryfaces( mesh ),
                         _rhs = l,
                         _element = uh,
                         _expr = cst( 2.0 ) );
            } );

    auto system = assemble( problem );
    BOOST_CHECK( system.constrained() );
    BOOST_CHECK( system.activeMatrixPtr() == system.matrixPtr() );
    BOOST_CHECK( system.activeVectorPtr() != system.rhsPtr() );

    auto result = systems::solve( system, _solution = uh );
    BOOST_CHECK( result.converged );

    auto boundaryError = Xh->element( "boundary-error" );
    boundaryError.on( _range = boundaryfaces( mesh ), _expr = abs( uh - cst( 2.0 ) ) );
    BOOST_CHECK_SMALL( boundaryError.max(), 1e-8 );
}

BOOST_AUTO_TEST_SUITE_END()
