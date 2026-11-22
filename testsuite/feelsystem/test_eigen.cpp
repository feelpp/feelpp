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

#define BOOST_TEST_MODULE feelsystem_eigen
#include <feel/feelcore/testsuite.hpp>

#include <cmath>

#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelsystem/system.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/options.hpp>

namespace
{

Feel::po::options_description
makeOptions()
{
    Feel::po::options_description options( "FeelSystem eigen test options" );
    options.add( Feel::backend_options( "eigen-smoke" ) )
           .add( Feel::backend_options( "eigen-smoke-alt" ) )
           .add( Feel::backend_options( "eigen-standard" ) )
           .add( Feel::solvereigen_options( "eigen-smoke" ) )
           .add( Feel::solvereigen_options( "eigen-smoke-alt" ) )
           .add( Feel::solvereigen_options( "eigen-standard" ) );
    return options;
}

template <typename SystemType, typename ModeType>
void
checkGeneralizedEigenResidual( SystemType const& system, double eigenvalue, ModeType const& mode )
{
    using namespace Feel;

    auto residual = system.backend()->newVector( _test = system.problem().space() );
    auto bMode = system.backend()->newVector( _test = system.problem().space() );
    residual->zero();
    bMode->zero();
    system.operatorA()->multVector( mode, *residual );
    system.operatorB()->multVector( mode, *bMode );
    residual->add( -eigenvalue, *bMode );
    residual->close();
    BOOST_CHECK_SMALL( residual->linftyNorm(), 1e-7 );
}

template <typename SystemType, typename ModeType>
void
checkStandardEigenResidual( SystemType const& system, double eigenvalue, ModeType const& mode )
{
    using namespace Feel;

    auto residual = system.backend()->newVector( _test = system.problem().space() );
    residual->zero();
    system.operatorA()->multVector( mode, *residual );
    residual->add( -eigenvalue, mode );
    residual->close();
    BOOST_CHECK_SMALL( residual->linftyNorm(), 1e-7 );
}

} // namespace

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAboutDefault( "feelpp_test_system_eigen" ), makeOptions() )

BOOST_AUTO_TEST_SUITE( feelsystem_eigen_suite )

BOOST_AUTO_TEST_CASE( generalized_eigen_problem_solves_known_minimum_reaction_diffusion_pencil )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto u = trial( Xh );
    auto v = test( Xh );
    auto stiffness = integrate( _range = elements( mesh ),
                                _expr = inner( grad( u ), grad( v ) ) );
    auto mass = integrate( _range = elements( mesh ),
                           _expr = u * v );

    auto problem =
        generalizedEigenvalueProblem( _space = Xh, _name = "eigen-smoke" )
            .operatorA( [&]( auto& a )
            {
                a += stiffness;
                a += mass;
            } )
            .operatorB( mass )
            .nev( 1 )
            .ncv( 3 )
            .tolerance( 1e-10 )
            .spectrum( systems::smallestMagnitude{} )
            .transform( systems::shiftInvert{} );

    static_assert( systems::GeneralizedEigenvalueProblemLike<decltype( problem )> );
    BOOST_CHECK_EQUAL( problem.operatorATerms().size(), 1 );
    BOOST_CHECK_EQUAL( problem.operatorBTerms().size(), 1 );
    BOOST_CHECK_EQUAL( problem.nev(), 1 );
    BOOST_CHECK_EQUAL( problem.ncv(), 3 );
    BOOST_CHECK( problem.hasSpectrum() );

    auto system = assemble( problem );
    BOOST_CHECK( system.assembled() );
    BOOST_CHECK_EQUAL( system.backendName(), problem.name() );
    BOOST_CHECK_EQUAL( system.backend()->prefix(), problem.name() );
    BOOST_CHECK( system.operatorA() );
    BOOST_CHECK( system.operatorB() );
    BOOST_CHECK( system.generalized() );

    auto memberSystem = problem.assemble();
    BOOST_CHECK( memberSystem.assembled() );
    BOOST_CHECK_EQUAL( memberSystem.backendName(), problem.name() );
    BOOST_CHECK( memberSystem.backend() == system.backend() );

    auto namedSystem = assemble( problem, _name = "eigen-smoke-alt" );
    BOOST_CHECK_EQUAL( namedSystem.backendName(), "eigen-smoke-alt" );
    BOOST_CHECK_EQUAL( namedSystem.backend()->prefix(), "eigen-smoke-alt" );

    auto memberNamedSystem = problem.assemble( _name = "eigen-smoke-alt" );
    BOOST_CHECK_EQUAL( memberNamedSystem.backendName(), "eigen-smoke-alt" );
    BOOST_CHECK( memberNamedSystem.backend() == namedSystem.backend() );

    auto result = systems::solve( system,
                                  _strategy = systems::eigen( _maxit = 1000 ) );
#if defined( FEELPP_HAS_SLEPC ) && defined( FEELPP_HAS_PETSC_H )
    BOOST_CHECK_GE( result.converged, 1 );
    BOOST_CHECK_EQUAL( system.eigenModes().size(), system.eigenvalues().size() );
    BOOST_CHECK_GE( system.eigenModes().size(), 1 );
    if ( !system.eigenvalues().empty() )
        BOOST_CHECK_SMALL( std::abs( system.eigenvalues().front() - 1.0 ), 1e-8 );
    for ( auto const& mode : system.eigenModes() )
        checkGeneralizedEigenResidual( system, mode.first, mode.second );

    auto memberResult = problem.solve( _strategy = systems::eigen( _nev = 1, _ncv = 3 ) );
    BOOST_CHECK_GE( memberResult.converged, 1 );
#else
    BOOST_CHECK_EQUAL( result.converged, 0 );
    BOOST_CHECK( system.eigenModes().empty() );
#endif
}

BOOST_AUTO_TEST_CASE( standard_eigen_problem_solves_exact_identity_operator )
{
    using namespace Feel;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    using space_type = typename std::remove_cvref_t<decltype( Xh )>::element_type;
    using problem_type = systems::Problem<space_type, systems::stationary_t, systems::generalized_eigen_t>;
    using space_ptrtype = typename problem_type::space_ptrtype;
    using matrix_ptrtype = typename problem_type::matrix_ptrtype;

    auto problem =
        eigenvalueProblem( _space = Xh, _name = "eigen-standard" )
            .operatorA( []( space_ptrtype const&, matrix_ptrtype const& matrix )
            {
                auto const nLocalRows = matrix->mapRow().nLocalDofWithoutGhost();
                for ( size_type k = 0; k != nLocalRows; ++k )
                    matrix->set( k, k, 1.0 );
            } )
            .nev( 2 )
            .spectrum( systems::largestMagnitude{} );

    auto system = problem.assemble();
    BOOST_CHECK( system.assembled() );
    BOOST_CHECK( system.operatorA() );
    BOOST_CHECK( !system.operatorB() );
    BOOST_CHECK( !system.generalized() );

    auto result = systems::solve( system,
                                  _strategy = systems::eigen( _ncv = 4,
                                                              _tolerance = 1e-10 ) );
#if defined( FEELPP_HAS_SLEPC ) && defined( FEELPP_HAS_PETSC_H )
    BOOST_CHECK_GE( result.converged, 1 );
    BOOST_CHECK_GE( system.eigenModes().size(), 1 );
    for ( auto const& mode : system.eigenModes() )
    {
        BOOST_CHECK_SMALL( std::abs( mode.first - 1.0 ), 1e-8 );
        checkStandardEigenResidual( system, mode.first, mode.second );
    }
#else
    BOOST_CHECK_EQUAL( result.converged, 0 );
#endif
}

BOOST_AUTO_TEST_SUITE_END()
