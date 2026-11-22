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

#define BOOST_TEST_MODULE feelsystem_evolution
#include <feel/feelcore/testsuite.hpp>

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
    Feel::po::options_description options( "FeelSystem evolution test options" );
    options.add( Feel::backend_options( "heat-step" ) )
        .add( Feel::backend_options( "heat-step-alt" ) )
        .add( Feel::backend_options( "wave-step" ) )
        .add( Feel::bdf_options( "heat-time" ) )
        .add( Feel::ts_options( "heat-time" ) )
        .add( Feel::ts_options( "wave-time" ) );
    return options;
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

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAboutDefault( "feelpp_test_system_evolution" ), makeOptions() )

BOOST_AUTO_TEST_SUITE( feelsystem_evolution_suite )

BOOST_AUTO_TEST_CASE( bdf_step_assembles_mass_stiffness_and_history )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto u = trial( Xh );
    auto v = test( Xh );
    auto backend = Feel::backend( _worldcomm = Xh->worldCommPtr() );

    auto massTerm = integrate( _range = elements( mesh ),
                               _expr = u * v );
    auto stiffnessTerm = integrate( _range = elements( mesh ),
                                    _expr = inner( grad( u ), grad( v ) ) );
    auto rhsTerm = integrate( _range = elements( mesh ),
                              _expr = cst( 3.0 ) * v );

    auto heatTime = bdf( _space = Xh,
                         _name = "heat-time",
                         _prefix = "heat-time",
                         _order = 1,
                         _initial_time = 0.0,
                         _final_time = 0.2,
                         _time_step = 0.1,
                         _save = false );

    auto problem =
        parabolicProblem( _space = Xh, _name = "heat-step" )
            .timeScheme( heatTime )
            .mass( massTerm )
            .stiffness( stiffnessTerm )
            .rhs( rhsTerm );

    static_assert( systems::FirstOrderEvolutionProblemLike<decltype( problem )> );
    BOOST_CHECK( problem.hasTimeScheme() );
    BOOST_CHECK( problem.hasBdfTimeScheme() );
    BOOST_CHECK_EQUAL( problem.timeSchemeTerms().size(), 1 );

    auto u0 = Xh->element( "u0" );
    u0.on( _range = elements( mesh ), _expr = cst( 2.0 ) );

    auto evolution = problem.evolve( _backend = backend );
    evolution.start( u0 );
    auto const context = evolution.context();
    BOOST_CHECK( context.isBdf() );
    BOOST_CHECK_CLOSE( context.time, 0.1, 1e-10 );
    evolution.assemble();
    BOOST_CHECK( evolution.activeMatrixPtr() );
    BOOST_CHECK( evolution.activeVectorPtr() );
    BOOST_CHECK( !evolution.constrained() );

    auto directMatrix = backend->newMatrix( _test = Xh, _trial = Xh );
    auto directRhs = backend->newVector( _test = Xh );
    auto directMass = backend->newMatrix( _test = Xh, _trial = Xh );
    auto history = backend->newVector( _test = Xh );
    directMatrix->zero();
    directRhs->zero();
    directMass->zero();

    form2( _test = Xh, _trial = Xh, _matrix = directMatrix ) +=
        integrate( _range = elements( mesh ),
                   _expr = context.derivativeCoefficient * u * v +
                           inner( grad( u ), grad( v ) ) );
    form1( _test = Xh, _vector = directRhs ) += rhsTerm;
    form2( _test = Xh, _trial = Xh, _matrix = directMass ) += massTerm;
    directMatrix->close();
    directRhs->close();
    directMass->close();
    history->zero();
    directMass->multVector( heatTime->polyDeriv(), *history );
    history->close();
    directRhs->add( 1.0, *history );
    directRhs->close();

    checkMatrixClose( evolution.matrixPtr(), directMatrix );
    checkVectorClose( evolution.rhsPtr(), directRhs );
}

BOOST_AUTO_TEST_CASE( bdf_step_updates_symbolic_time_and_backend_names )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto u = trial( Xh );
    auto v = test( Xh );

    auto heatTime = bdf( _space = Xh,
                         _name = "heat-time",
                         _prefix = "heat-time",
                         _order = 1,
                         _initial_time = 0.0,
                         _final_time = 0.1,
                         _time_step = 0.1,
                         _save = false );

    auto timeSource = expr( "t:t" );
    auto problem =
        parabolicProblem( _space = Xh, _name = "heat-step" )
            .timeScheme( heatTime )
            .mass( integrate( _range = elements( mesh ),
                              _expr = u * v ) )
            .rhs( _range = elements( mesh ),
                  _expr = timeSource * v );

    auto u0 = Xh->element( "u0" );
    u0.zero();
    auto evolution = problem.evolve( _name = "heat-step-alt" );
    evolution.start( u0 );
    evolution.assemble();

    BOOST_CHECK_EQUAL( evolution.backendName(), "heat-step-alt" );
    BOOST_CHECK_EQUAL( evolution.backend()->prefix(), "heat-step-alt" );
    BOOST_CHECK_CLOSE( evolution.problem().parameterValues().at( "t" ), evolution.time(), 1e-10 );

    auto directRhs = evolution.backend()->newVector( _test = Xh );
    timeSource.setParameterValues( { { "t", evolution.time() } } );
    directRhs->zero();
    form1( _test = Xh, _vector = directRhs ) +=
        integrate( _range = elements( mesh ),
                   _expr = timeSource * v );
    directRhs->close();

    checkVectorClose( evolution.rhsPtr(), directRhs );
}

BOOST_AUTO_TEST_CASE( newmark_step_assembles_mass_damping_and_stiffness )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto u = trial( Xh );
    auto v = test( Xh );
    auto backend = Feel::backend( _worldcomm = Xh->worldCommPtr() );

    auto waveTime = newmark( _space = Xh,
                             _name = "wave-time",
                             _prefix = "wave-time",
                             _initial_time = 0.0,
                             _final_time = 0.1,
                             _time_step = 0.1,
                             _gamma = 0.5,
                             _beta = 0.25,
                             _save = false );

    auto problem =
        hyperbolicProblem( _space = Xh, _name = "wave-step" )
            .timeScheme( waveTime )
            .mass( integrate( _range = elements( mesh ),
                              _expr = 2.0 * u * v ) )
            .damping( integrate( _range = elements( mesh ),
                                 _expr = 0.5 * u * v ) )
            .stiffness( integrate( _range = elements( mesh ),
                                   _expr = inner( grad( u ), grad( v ) ) ) )
            .rhs( integrate( _range = elements( mesh ),
                             _expr = v ) );

    static_assert( systems::SecondOrderEvolutionProblemLike<decltype( problem )> );

    auto u0 = Xh->element( "u0" );
    u0.zero();
    auto evolution = evolve( problem, _backend = backend );
    evolution.start( u0 );
    auto const context = evolution.context();
    BOOST_CHECK( context.isNewmark() );
    evolution.assemble();

    auto directMatrix = backend->newMatrix( _test = Xh, _trial = Xh );
    auto directRhs = backend->newVector( _test = Xh );
    directMatrix->zero();
    directRhs->zero();
    form2( _test = Xh, _trial = Xh, _matrix = directMatrix ) +=
        integrate( _range = elements( mesh ),
                   _expr = context.secondDerivativeCoefficient * 2.0 * u * v +
                           context.firstDerivativeCoefficient * 0.5 * u * v +
                           inner( grad( u ), grad( v ) ) );
    form1( _test = Xh, _vector = directRhs ) +=
        integrate( _range = elements( mesh ),
                   _expr = v );
    directMatrix->close();
    directRhs->close();

    checkMatrixClose( evolution.matrixPtr(), directMatrix );
    checkVectorClose( evolution.rhsPtr(), directRhs );
}

BOOST_AUTO_TEST_SUITE_END()
