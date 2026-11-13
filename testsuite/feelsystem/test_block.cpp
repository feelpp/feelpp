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

#define BOOST_TEST_MODULE feelsystem_block
#include <feel/feelcore/testsuite.hpp>

#include <boost/hana/integral_constant.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/dh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelsystem/system.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/options.hpp>

namespace
{

Feel::po::options_description
makeOptions()
{
    Feel::po::options_description options( "FeelSystem block test options" );
    options.add( Feel::backend_options( "block-stokes" ) )
        .add( Feel::backend_options( "block-darcy" ) )
        .add( Feel::backend_options( "block-static-condensation" ) )
        .add( Feel::backend_options( Feel::prefixvm( "block-static-condensation", "sc" ) ) );
    return options;
}

} // namespace

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAboutDefault( "feelpp_test_system_block" ), makeOptions() )

BOOST_AUTO_TEST_SUITE( feelsystem_block_suite )

BOOST_AUTO_TEST_CASE( stokes_manufactured_solution_block_problem )
{
    using namespace Feel;
    using namespace Feel::vf;
    using namespace boost::hana::literals;

    auto mesh = unitSquare();
    auto Vh = Pchv<2>( mesh );
    auto Qh = Pch<1>( mesh );
    auto Xh = product( Vh, Qh );
    auto solution = Xh.element();
    auto boundaryVelocity = Vh->element( "u_dirichlet" );

    auto u = trial( Vh );
    auto v = test( Vh );
    auto p = trial( Qh );
    auto q = test( Qh );

    auto uExact = vec( Px() * Px(), -2.0 * Px() * Py() );
    auto fExact = vec( cst( -2.0 ), cst( 0.0 ) );
    double const pressurePenalty = 1e-6;

    auto problem =
        mixedProblem( _space = Xh, _name = "block-stokes" )
            .block( 0_c, 0_c,
                    integrate( _range = elements( mesh ),
                               _expr = inner( grad( u ), grad( v ) ) ) )
            .block( 0_c, 1_c,
                    integrate( _range = elements( mesh ),
                               _expr = -p * div( v ) ) )
            .block( 1_c, 0_c,
                    integrate( _range = elements( mesh ),
                               _expr = -div( u ) * q ) )
            .block( 1_c, 1_c,
                    integrate( _range = elements( mesh ),
                               _expr = pressurePenalty * p * q ) )
            .row( 0_c,
                  integrate( _range = elements( mesh ),
                             _expr = inner( fExact, v ) ) )
            .constraint( [mesh, boundaryVelocity, uExact]( auto& a, auto& l )
                         { a.row( 0_c ) += on( _range = boundaryfaces( mesh ),
                                               _rhs = l( 0_c ),
                                               _element = boundaryVelocity,
                                               _expr = uExact,
                                               _type = "elimination" ); } );

    static_assert( Feel::systems::BlockProblemLike<decltype( problem )> );
    BOOST_CHECK_EQUAL( problem.lhsTerms().size(), 4 );
    BOOST_CHECK_EQUAL( problem.rhsTerms().size(), 1 );
    BOOST_CHECK_EQUAL( problem.constraintTerms().size(), 1 );

    auto system = problem.assemble();
    BOOST_CHECK( system.assembled() );
    BOOST_CHECK_EQUAL( system.backendName(), problem.name() );
    BOOST_CHECK( system.constrained() );
    BOOST_CHECK( system.matrixForm().hasDirichletConstraints() );
    BOOST_CHECK( system.matrixForm().hasMaterializedConstrainedOperator() );
    BOOST_CHECK( system.matrixForm().dirichletInPlace() );
    BOOST_CHECK( !system.matrixForm().hasUnconstrainedMatrix() );
    BOOST_CHECK( system.activeMatrixPtr() );
    BOOST_CHECK( system.activeVectorPtr() );

    auto result = system.solve( _solution = solution );
    BOOST_CHECK( result.converged );

    double const velocityError = normL2( _range = elements( mesh ),
                                         _expr = solution( 0_c ) - uExact );
    double const pressureNorm = normL2( _range = elements( mesh ),
                                        _expr = solution( 1_c ) );
    double const boundaryError = normL2( _range = boundaryfaces( mesh ),
                                         _expr = solution( 0_c ) - uExact );

    BOOST_CHECK_SMALL( velocityError, 1e-8 );
    BOOST_CHECK_SMALL( pressureNorm, 1e-8 );
    BOOST_CHECK_SMALL( boundaryError, 1e-10 );
}

BOOST_AUTO_TEST_CASE( darcy_manufactured_solution_block_problem )
{
    using namespace Feel;
    using namespace Feel::vf;
    using namespace boost::hana::literals;

    auto mesh = unitSquare();
    auto Vh = Dh<0>( mesh );
    auto Qh = Pdh<0>( mesh );
    auto Xh = product( Vh, Qh );
    auto solution = Xh.element();

    auto sigma = trial( Vh );
    auto tau = test( Vh );
    auto p = trial( Qh );
    auto q = test( Qh );

    auto sigmaExact = oneX();
    double const pressurePenalty = 1.0;

    auto problem =
        mixedProblem( _space = Xh, _name = "block-darcy" )
            .block( 0_c, 0_c,
                    integrate( _range = elements( mesh ),
                               _expr = inner( sigma, tau ) ) )
            .block( 0_c, 1_c,
                    integrate( _range = elements( mesh ),
                               _expr = -p * div( tau ) ) )
            .block( 1_c, 0_c,
                    integrate( _range = elements( mesh ),
                               _expr = div( sigma ) * q ) )
            .block( 1_c, 1_c,
                    integrate( _range = elements( mesh ),
                               _expr = pressurePenalty * p * q ) )
            .row( 0_c,
                  integrate( _range = elements( mesh ),
                             _expr = inner( sigmaExact, tau ) ) );

    auto result = problem.solve( _solution = solution );
    BOOST_CHECK( result.converged );

    double const fluxError = normL2( _range = elements( mesh ),
                                     _expr = solution( 0_c ) - sigmaExact );
    double const pressureNorm = normL2( _range = elements( mesh ),
                                        _expr = solution( 1_c ) );
    double const divergenceNorm = normL2( _range = elements( mesh ),
                                          _expr = divv( solution( 0_c ) ) );

    BOOST_CHECK_SMALL( fluxError, 1e-10 );
    BOOST_CHECK_SMALL( pressureNorm, 1e-10 );
    BOOST_CHECK_SMALL( divergenceNorm, 1e-12 );
}

BOOST_AUTO_TEST_CASE( static_condensation_block_problem_matches_monolithic )
{
    using namespace Feel;
    using namespace Feel::vf;
    using namespace boost::hana::literals;

    auto mesh = unitSquare();
    auto Uh = Pch<1>( mesh );
    auto Ah = Pdh<0>( mesh );
    auto Xh = product( Uh, Ah );

    auto W = Xh.element();
    auto T = Xh.element();
    auto u = W( 0_c );
    auto alpha = W( 1_c );
    auto v = T( 0_c );
    auto beta = T( 1_c );

    auto fillOperator = [mesh, u, alpha, v, beta]( auto& a )
    {
        a( 0_c, 0_c ) += integrate( _range = elements( mesh ),
                                    _expr = inner( gradt( u ), grad( v ) ) + idt( u ) * id( v ) );
        a( 0_c, 1_c ) += integrate( _range = elements( mesh ),
                                    _expr = idt( alpha ) * id( v ) );
        a( 1_c, 0_c ) += integrate( _range = elements( mesh ),
                                    _expr = idt( u ) * id( beta ) );
        a( 1_c, 1_c ) += integrate( _range = elements( mesh ),
                                    _expr = cst( 2.0 ) * idt( alpha ) * id( beta ) );
    };
    auto fillRhs = [mesh, v]( auto& l )
    {
        l( 0_c ) += integrate( _range = elements( mesh ),
                               _expr = cst( 1.0 ) * id( v ) );
    };
    auto applyBoundary = [mesh, u]( auto& a, auto& l )
    {
        a.row( 0_c ) += on( _range = boundaryfaces( mesh ),
                            _rhs = l( 0_c ),
                            _element = u,
                            _expr = cst( 0.0 ),
                            _type = "elimination" );
    };

    auto monolithicProblem =
        blockProblem( _space = Xh, _name = "block-static-condensation" )
            .lhs( fillOperator )
            .rhs( fillRhs )
            .constraint( applyBoundary );

    auto condensedProblem =
        blockProblem( _space = Xh, _name = "block-static-condensation" )
            .strategy( solve::strategy::static_condensation )
            .lhs( fillOperator )
            .rhs( fillRhs )
            .constraint( applyBoundary );

    auto monolithicSolution = Xh.element();
    auto condensedSolution = Xh.element();

    auto monolithicResult = monolithicProblem.solve( _solution = monolithicSolution );
    auto condensedSystem = condensedProblem.assemble();
    BOOST_CHECK( condensedSystem.activeMatrixPtr() );
    BOOST_CHECK( condensedSystem.activeVectorPtr() );
    auto condensedResult = condensedSystem.solve( _solution = condensedSolution,
                                                  _condense = true,
                                                  _condenser = condenser_sb9() );

    BOOST_CHECK( monolithicResult.converged );
    BOOST_CHECK( condensedResult.converged );

    double const uError = normL2( _range = elements( mesh ),
                                  _expr = monolithicSolution( 0_c ) - condensedSolution( 0_c ) );
    double const alphaError = normL2( _range = elements( mesh ),
                                      _expr = monolithicSolution( 1_c ) - condensedSolution( 1_c ) );
    double const boundaryError = normL2( _range = boundaryfaces( mesh ),
                                         _expr = condensedSolution( 0_c ) );

    BOOST_CHECK_SMALL( uError, 2e-10 );
    BOOST_CHECK_SMALL( alphaError, 2e-10 );
    BOOST_CHECK_SMALL( boundaryError, 1e-12 );
}

BOOST_AUTO_TEST_SUITE_END()
