/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file test_broken_fe.cpp
 * @brief Contract and assembly-topology tests for Broken<Family>.
 */
#define BOOST_TEST_MODULE test_broken_fe
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/brokenh.hpp>
#include <feel/feelpoly/brezzidouglasmarini.hpp>
#include <feel/feelpoly/lagrange.hpp>
#include <feel/feelpoly/nedelec.hpp>
#include <feel/feelpoly/raviartthomas.hpp>
#include <feel/feelfilters/unitsquare.hpp>

#include <algorithm>
#include <set>
#include <vector>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
template<class BrokenFE, class WrappedFE>
void checkLocalBrokenContract( BrokenFE const& broken, WrappedFE const& wrapped )
{
    static_assert( ReferenceFiniteElement<BrokenFE> );
    static_assert( is_broken_v<BrokenFE> );
    static_assert( BrokenFE::continuity_type::is_discontinuous_totally );
    static_assert( !HDivFiniteElement<BrokenFE> );
    static_assert( !HCurlFiniteElement<BrokenFE> );

    BOOST_CHECK_EQUAL( broken.order(), wrapped.order() );
    BOOST_CHECK_EQUAL( broken.polynomialDegree(), wrapped.polynomialDegree() );
    BOOST_CHECK_EQUAL( broken.localDofCount(), wrapped.localDofCount() );
    auto const brokenDualOnPrimal = broken.dual()( broken.primal() );
    auto const wrappedDualOnPrimal = wrapped.dual()( wrapped.primal() );
    BOOST_CHECK_EQUAL( brokenDualOnPrimal.size1(), wrappedDualOnPrimal.size1() );
    BOOST_CHECK_EQUAL( brokenDualOnPrimal.size2(), wrappedDualOnPrimal.size2() );

    BOOST_CHECK_EQUAL( broken.localDofCountOnFacet( 0 ), 0 );
    BOOST_CHECK_EQUAL( broken.localDofCountOnEntity( 0, 0 ), 0 );
    BOOST_CHECK_EQUAL( broken.localDofCountOnEntity( 1, 0 ), 0 );
    BOOST_CHECK_EQUAL( broken.localDofCountOnEntity( BrokenFE::nDim, 0 ), broken.localDofCount() );

    std::vector<bool> ordinals( broken.localDofPerComponent(), false );
    for ( uint16_type localDof = 0; localDof < broken.localDofCount(); ++localDof )
    {
        auto const layout = broken.localDofLayout( localDof );
        BOOST_CHECK_EQUAL( layout.attachment.entityDim, BrokenFE::nDim );
        BOOST_CHECK_EQUAL( layout.attachment.entityId, 0 );
        BOOST_REQUIRE_LT( layout.attachment.ordinal, ordinals.size() );
        if ( layout.component == 0 )
            ordinals[layout.attachment.ordinal] = true;
    }
    BOOST_CHECK( std::all_of( ordinals.begin(), ordinals.end(), []( bool value ) { return value; } ) );
}

template<class SpacePtr>
void checkGlobalBrokenCardinality( SpacePtr const& Xh )
{
    using space_type = typename SpacePtr::element_type;
    using size_type = typename space_type::size_type;
    auto const globalElements = nelements( elements( Xh->mesh() ), true );
    auto const localDof = Xh->basis()->localDofCount();
    BOOST_CHECK_EQUAL( Xh->nDof(), globalElements * localDof );

    for ( auto const& eltWrap : elements( Xh->mesh(), entity_process_t::LOCAL_ONLY ) )
    {
        auto const& elt = unwrap_ref( eltWrap );
        auto const indices = Xh->dof()->localToGlobalIndices( elt.id() );
        BOOST_CHECK_EQUAL( indices.size(), localDof );
        std::set<size_type> unique( indices.begin(), indices.end() );
        BOOST_CHECK_EQUAL( unique.size(), localDof );
    }
}
} // namespace

BOOST_AUTO_TEST_SUITE( broken_fe_suite )

BOOST_AUTO_TEST_CASE( local_static_and_runtime_wrappers_preserve_ciarlet_data )
{
    using rt_static = typename RaviartThomas<1>::template apply<2>::type;
    using rt_broken = typename Broken<RaviartThomas<1>>::template apply<2>::type;
    using rt_dynamic = typename RaviartThomas<Dynamic>::template apply<2>::type;
    using rt_dynamic_broken = typename Broken<RaviartThomas<Dynamic>>::template apply<2>::type;
    using bdm_static = typename BrezziDouglasMarini<1>::template apply<2>::type;
    using bdm_broken = typename Broken<BrezziDouglasMarini<1>>::template apply<2>::type;
    using ned_static = typename Nedelec<0>::template apply<2>::type;
    using ned_broken = typename Broken<Nedelec<0>>::template apply<2>::type;

    checkLocalBrokenContract( rt_broken{}, rt_static{} );
    checkLocalBrokenContract( rt_dynamic_broken{ RuntimeOrder{ 1 } }, rt_dynamic{ RuntimeOrder{ 1 } } );
    checkLocalBrokenContract( bdm_broken{}, bdm_static{} );
    checkLocalBrokenContract( ned_broken{}, ned_static{} );
}

BOOST_AUTO_TEST_CASE( global_dofs_are_cell_local_for_enabled_families )
{
    auto mesh = unitSquare( 0.4 );

    checkGlobalBrokenCardinality( brokenh<RaviartThomas<0>>( mesh ) );
    checkGlobalBrokenCardinality( brokenh<RaviartThomas<Dynamic>>( mesh, RuntimeOrder{ 1 } ) );
    checkGlobalBrokenCardinality( brokenh<BrezziDouglasMarini<1>>( mesh ) );
    checkGlobalBrokenCardinality( brokenh<BrezziDouglasMarini<Dynamic>>( mesh, RuntimeOrder{ 1 } ) );
    checkGlobalBrokenCardinality( brokenh<Nedelec<0>>( mesh ) );
    checkGlobalBrokenCardinality( brokenh<Nedelec<Dynamic>>( mesh, RuntimeOrder{ 0 } ) );
}

BOOST_AUTO_TEST_SUITE_END()
