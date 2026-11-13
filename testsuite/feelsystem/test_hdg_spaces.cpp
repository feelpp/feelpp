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

#define BOOST_TEST_MODULE feelsystem_hdg_spaces
#include <feel/feelcore/testsuite.hpp>

#include <string>
#include <string_view>

#include <boost/hana/integral_constant.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelsystem/system.hpp>
#include <feel/options.hpp>

struct boundary_state_role_t {};

namespace Feel::systems
{
template <>
struct block_role_traits<boundary_state_role_t>
{
    inline static constexpr std::string_view name = "boundary_state";
};
} // namespace Feel::systems

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAboutDefault( "feelpp_test_system_hdg_spaces" ),
                                 Feel::po::options_description( "FeelSystem HDG spaces test options" ) )

BOOST_AUTO_TEST_SUITE( feelsystem_hdg_spaces_suite )

BOOST_AUTO_TEST_CASE( generic_block_roles_are_extensible )
{
    namespace role = Feel::systems::roles;

    static_assert( Feel::systems::BlockRole<decltype( role::flux )> );
    static_assert( Feel::systems::BlockRole<decltype( role::potential )> );
    static_assert( Feel::systems::BlockRole<boundary_state_role_t> );

    BOOST_CHECK_EQUAL( std::string( Feel::systems::blockRoleName( role::flux ) ), "flux" );
    BOOST_CHECK_EQUAL( std::string( Feel::systems::blockRoleName( role::potential ) ), "potential" );
    BOOST_CHECK_EQUAL( std::string( Feel::systems::blockRoleName( boundary_state_role_t{} ) ), "boundary_state" );
}

BOOST_AUTO_TEST_CASE( poisson_hdg_spaces_match_quickstart_product_layout )
{
    using namespace Feel;
    using namespace boost::hana::literals;
    namespace role = Feel::systems::roles;

    auto mesh = unitSquare();
    auto H = hdg::poissonSpaces<1>( mesh );

    using hdg_spaces_type = std::remove_cvref_t<decltype( H )>;
    static_assert( hdg_spaces_type::nroles == 3 );
    static_assert( hdg_spaces_type::indexValue( role::flux ) == 0 );
    static_assert( hdg_spaces_type::indexValue( role::potential ) == 1 );
    static_assert( hdg_spaces_type::indexValue( role::trace ) == 2 );

    BOOST_CHECK( H.mesh() == mesh );
    BOOST_CHECK( H.facetMesh() );
    BOOST_CHECK_EQUAL( H.numberOfRoles(), 3 );
    BOOST_CHECK_EQUAL( H.numberOfSpaces(), 3 );
    BOOST_CHECK( H.hasRole( role::flux ) );
    BOOST_CHECK( H.hasRole( role::potential ) );
    BOOST_CHECK( H.hasRole( role::trace ) );
    BOOST_CHECK( !H.hasRole( role::pressure ) );

    auto flux = H.flux();
    auto potential = H.potential();
    auto trace = H.trace();

    BOOST_CHECK( flux.index() == 0_c );
    BOOST_CHECK( potential.index() == 1_c );
    BOOST_CHECK( trace.index() == 2_c );
    BOOST_CHECK_EQUAL( std::string( flux.name() ), "flux" );
    BOOST_CHECK_EQUAL( std::string( potential.name() ), "potential" );
    BOOST_CHECK_EQUAL( std::string( trace.name() ), "trace" );

    BOOST_CHECK( flux.space() == H.product()[0_c] );
    BOOST_CHECK( potential.space() == H.product()[1_c] );
    BOOST_CHECK( trace.space() == H.product()[2_c] );

    auto manualProduct = product( flux.space(), potential.space(), trace.space() );
    BOOST_CHECK( manualProduct == H.product() );
    BOOST_CHECK_EQUAL( H.product().nDof(),
                       flux.space()->nDof() + potential.space()->nDof() + trace.space()->nDof() );

    auto U = H.element();
    BOOST_CHECK_EQUAL( U.functionSpace().numberOfSpaces(), 3 );

    auto u = flux.trial();
    auto v = flux.test();
    auto p = potential.trial();
    auto q = potential.test();
    auto phat = trace.trial();
    auto mu = trace.test();

    static_assert( Feel::vf::is_trial_basis_v<decltype( u )> );
    static_assert( Feel::vf::is_test_basis_v<decltype( v )> );
    static_assert( Feel::vf::is_trial_basis_v<decltype( p )> );
    static_assert( Feel::vf::is_test_basis_v<decltype( q )> );
    static_assert( Feel::vf::is_trial_basis_v<decltype( phat )> );
    static_assert( Feel::vf::is_test_basis_v<decltype( mu )> );
}

BOOST_AUTO_TEST_CASE( explicit_hdg_spaces_accept_user_defined_roles )
{
    using namespace Feel;
    namespace role = Feel::systems::roles;

    auto mesh = unitSquare();
    auto facetMesh = createSubmesh( _mesh = mesh, _range = faces( mesh ), _update = 0 );
    auto Ch = Pdh<0>( facetMesh );

    auto H = systems::makeHdgSpaces<boundary_state_role_t>(
        mesh, facetMesh, product( Ch ) );

    using hdg_spaces_type = std::remove_cvref_t<decltype( H )>;
    static_assert( hdg_spaces_type::nroles == 1 );
    static_assert( hdg_spaces_type::indexValue( boundary_state_role_t{} ) == 0 );

    BOOST_CHECK( H.hasRole( boundary_state_role_t{} ) );
    BOOST_CHECK( !H.hasRole( role::trace ) );

    auto boundaryState = H.component( boundary_state_role_t{} );
    BOOST_CHECK_EQUAL( std::string( boundaryState.name() ), "boundary_state" );
    BOOST_CHECK_EQUAL( boundaryState.space()->nDof(), Ch->nDof() );
}

BOOST_AUTO_TEST_SUITE_END()
