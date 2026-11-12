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

#define BOOST_TEST_MODULE feelsystem_hdg
#include <feel/feelcore/testsuite.hpp>

#include <boost/hana/integral_constant.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelsystem/system.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/options.hpp>

namespace
{

Feel::po::options_description
makeOptions()
{
    Feel::po::options_description options( "FeelSystem HDG problem test options" );
    options.add( Feel::backend_options( "hdg-role" ) )
           .add( Feel::backend_options( "hdg-role-explicit" ) )
           .add( Feel::backend_options( "hdg-static-mono" ) )
           .add( Feel::backend_options( "hdg-static-condensed" ) )
           .add( Feel::backend_options( Feel::prefixvm( "hdg-static-condensed", "sc" ) ) );
    return options;
}

template <typename MatrixPtr>
void
checkMatrixClose( MatrixPtr const& actual, MatrixPtr const& expected )
{
    auto delta = actual->clone();
    delta->addMatrix( -1.0, expected );
    delta->close();
    BOOST_CHECK_SMALL( delta->l1Norm(), 1e-9 );
    BOOST_CHECK_SMALL( delta->linftyNorm(), 1e-9 );
}

template <typename VectorPtr>
void
checkVectorClose( VectorPtr const& actual, VectorPtr const& expected )
{
    auto delta = Feel::vf::cloneVectorWithValues( actual );
    delta->add( -1.0, *expected );
    delta->close();
    BOOST_CHECK_SMALL( delta->linftyNorm(), 1e-10 );
}

} // namespace

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAboutDefault( "feelpp_test_system_hdg" ), makeOptions() )

BOOST_AUTO_TEST_SUITE( feelsystem_hdg_suite )

BOOST_AUTO_TEST_CASE( hdg_problem_records_policy_and_backend_names )
{
    using namespace Feel;

    auto mesh = unitSquare();
    auto H = hdg::poissonSpaces<1>( mesh );

    auto problem =
        hdgProblem( _spaces = H, _name = "hdg-role" )
            .condensation( hdg::condensation()
                               .enabled( true )
                               .condenser( condenser_poisson() )
                               .backendName( prefixvm( "hdg-role", "sc" ) ) );

    static_assert( systems::HdgProblemLike<decltype( problem )> );
    BOOST_CHECK_EQUAL( problem.name(), "hdg-role" );
    BOOST_CHECK_EQUAL( problem.nTerms(), 0 );
    BOOST_CHECK( problem.condensationPolicy().enabled() );
    BOOST_CHECK( problem.strategy() == solve::strategy::static_condensation );
    BOOST_CHECK_EQUAL( problem.condensedBackendName(), prefixvm( "hdg-role", "sc" ) );
    BOOST_CHECK_EQUAL( problem.postBackendName(), prefixvm( "hdg-role", "post" ) );
    BOOST_CHECK_EQUAL( problem.condensationPolicy().solveBackendBaseName( "hdg-role" ), "hdg-role" );

    auto system = problem.assemble();
    BOOST_CHECK( system.assembled() );
    BOOST_CHECK_EQUAL( system.backendName(), "hdg-role" );
    BOOST_CHECK_EQUAL( system.condensedBackendName(), prefixvm( "hdg-role", "sc" ) );
    BOOST_CHECK_EQUAL( system.postBackendName(), prefixvm( "hdg-role", "post" ) );
    BOOST_CHECK( system.strategy() == solve::strategy::static_condensation );

    auto customPolicy =
        hdg::condensation().backendName( prefixvm( "custom-hdg-role", "sc" ) );
    BOOST_CHECK_EQUAL( customPolicy.solveBackendBaseName( "hdg-role" ), "custom-hdg-role" );
}

BOOST_AUTO_TEST_CASE( role_based_hdg_assembly_matches_explicit_block_indices )
{
    using namespace Feel;
    using namespace Feel::vf;
    using namespace boost::hana::literals;
    namespace role = Feel::systems::roles;

    auto mesh = unitSquare();
    auto H = hdg::poissonSpaces<1>( mesh );

    auto u = H.flux().trial();
    auto v = H.flux().test();
    auto p = H.potential().trial();
    auto q = H.potential().test();
    auto phat = H.trace().trial();

    auto fillTraceBlock = [mesh]( auto const& X, auto& a )
    {
        auto uLocal = X.flux().trial();
        auto phat = X.trace().trial();
        auto mu = X.trace().test();
        a( role::trace, role::flux ) +=
            integrate( _range = facet::interior( mesh ),
                       _expr = facet::trace( mu ) * facet::normalSum( uLocal ) )
            + integrate( _range = facet::boundary( mesh ),
                         _expr = facet::trace( mu ) * facet::normal( uLocal ) );
        a( role::trace, role::trace ) +=
            integrate( _range = boundaryfaces( mesh ),
                       _expr = phat * mu );
    };

    auto fillTraceBlockExplicit = [mesh, H]( auto& a )
    {
        auto uLocal = H.flux().trial();
        auto phat = H.trace().trial();
        auto mu = H.trace().test();
        a( 2_c, 0_c ) +=
            integrate( _range = internalfaces( mesh ),
                       _expr = mu * ( leftfacet( normal( uLocal ) ) +
                                      rightfacet( normal( uLocal ) ) ) )
            + integrate( _range = boundaryfaces( mesh ),
                         _expr = mu * normal( uLocal ) );
        a( 2_c, 2_c ) +=
            integrate( _range = boundaryfaces( mesh ),
                       _expr = phat * mu );
    };

    auto fillRhs = [mesh]( auto const& X, auto& l )
    {
        auto qLocal = X.potential().test();
        auto mu = X.trace().test();
        l( role::potential ) +=
            integrate( _range = elements( mesh ),
                       _expr = qLocal );
        l( role::trace ) +=
            integrate( _range = boundaryfaces( mesh ),
                       _expr = mu );
    };

    auto fillRhsExplicit = [mesh, H]( auto& l )
    {
        auto qLocal = H.potential().test();
        auto mu = H.trace().test();
        l( 1_c ) +=
            integrate( _range = elements( mesh ),
                       _expr = qLocal );
        l( 2_c ) +=
            integrate( _range = boundaryfaces( mesh ),
                       _expr = mu );
    };

    auto roleProblem =
        hdgProblem( _spaces = H, _name = "hdg-role" )
            .block( role::flux, role::flux,
                    integrate( _range = elements( mesh ),
                               _expr = inner( u, v ) ) )
            .block( role::potential, role::potential,
                    integrate( _range = elements( mesh ),
                               _expr = p * q ) )
            .block( role::flux, role::trace,
                    integrate( _range = facet::interior( mesh ),
                               _expr = facet::trace( phat ) * facet::normalSum( v ) )
                    + integrate( _range = facet::boundary( mesh ),
                                 _expr = facet::trace( phat ) * facet::normal( v ) ) )
            .lhs( fillTraceBlock )
            .rhs( fillRhs );

    auto explicitProblem =
        blockProblem( _space = H.product(), _name = "hdg-role-explicit" )
            .block( 0_c, 0_c,
                    integrate( _range = elements( mesh ),
                               _expr = inner( u, v ) ) )
            .block( 1_c, 1_c,
                    integrate( _range = elements( mesh ),
                               _expr = p * q ) )
            .block( 0_c, 2_c,
                    integrate( _range = internalfaces( mesh ),
                               _expr = phat * ( leftface( normal( v ) ) +
                                                rightface( normal( v ) ) ) )
                    + integrate( _range = boundaryfaces( mesh ),
                                 _expr = phat * normal( v ) ) )
            .lhs( fillTraceBlockExplicit )
            .rhs( fillRhsExplicit );

    BOOST_CHECK_EQUAL( roleProblem.lhsTerms().size(), explicitProblem.lhsTerms().size() );
    BOOST_CHECK_EQUAL( roleProblem.rhsTerms().size(), explicitProblem.rhsTerms().size() );

    auto roleSystem = assemble( roleProblem );
    auto explicitSystem = assemble( explicitProblem );

    BOOST_CHECK( roleSystem.assembled() );
    BOOST_CHECK( explicitSystem.assembled() );
    checkMatrixClose( roleSystem.activeMatrixPtr(), explicitSystem.activeMatrixPtr() );
    checkVectorClose( roleSystem.activeVectorPtr(), explicitSystem.activeVectorPtr() );
}

BOOST_AUTO_TEST_CASE( hdg_condensation_policy_drives_solve_arguments )
{
    using namespace Feel;
    using namespace Feel::vf;
    using namespace boost::hana::literals;
    namespace role = Feel::systems::roles;

    auto mesh = unitSquare();
    auto facetMesh = createSubmesh( _mesh = mesh, _range = faces( mesh ), _update = 0 );
    auto Uh = Pch<1>( mesh );
    auto Ah = Pdh<0>( mesh );
    auto H = systems::makeHdgSpaces<systems::roles::potential_t,
                                    systems::roles::constraint_t>(
        mesh, facetMesh, product( Uh, Ah ) );

    auto W = H.element();
    auto T = H.element();
    auto u = W( 0_c );
    auto alpha = W( 1_c );
    auto v = T( 0_c );
    auto beta = T( 1_c );

    auto fillOperator = [mesh, u, alpha, v, beta]( auto& a )
    {
        a( role::potential, role::potential ) +=
            integrate( _range = elements( mesh ),
                       _expr = inner( gradt( u ), grad( v ) ) + idt( u ) * id( v ) );
        a( role::potential, role::constraint ) +=
            integrate( _range = elements( mesh ),
                       _expr = idt( alpha ) * id( v ) );
        a( role::constraint, role::potential ) +=
            integrate( _range = elements( mesh ),
                       _expr = idt( u ) * id( beta ) );
        a( role::constraint, role::constraint ) +=
            integrate( _range = elements( mesh ),
                       _expr = cst( 2.0 ) * idt( alpha ) * id( beta ) );
    };
    auto fillRhs = [mesh, v]( auto& l )
    {
        l( role::potential ) +=
            integrate( _range = elements( mesh ),
                       _expr = cst( 1.0 ) * id( v ) );
    };
    auto applyBoundary = [mesh, u]( auto& a, auto& l )
    {
        a.row( role::potential ) +=
            on( _range = boundaryfaces( mesh ),
                _rhs = l( role::potential ),
                _element = u,
                _expr = cst( 0.0 ),
                _type = "elimination" );
    };

    auto monolithicProblem =
        hdgProblem( _spaces = H, _name = "hdg-static-mono" )
            .lhs( fillOperator )
            .rhs( fillRhs )
            .constraint( applyBoundary );

    auto condensedProblem =
        hdgProblem( _spaces = H, _name = "hdg-static-condensed" )
            .condensation( hdg::condensation()
                               .enabled( true )
                               .condenser( condenser_sb9() )
                               .backendName( prefixvm( "hdg-static-condensed", "sc" ) ) )
            .lhs( fillOperator )
            .rhs( fillRhs )
            .constraint( applyBoundary );

    auto monolithicSolution = H.element();
    auto condensedSolution = H.element();

    auto monolithicResult = monolithicProblem.solve( _solution = monolithicSolution );
    auto condensedSystem = condensedProblem.assemble();
    BOOST_CHECK_EQUAL( condensedSystem.condensedBackendName(), prefixvm( "hdg-static-condensed", "sc" ) );
    auto condensedResult = condensedSystem.solve( _solution = condensedSolution );

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
