/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/** @file testmixedassembly.cpp
 * @brief Actual scalar/component/two-cell assembly with independently static/runtime FE orders.
 */
#define BOOST_TEST_MODULE mixedassembly
#include <boost/mp11/list.hpp>
#include <feel/feelcore/testsuite.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feelfilters/unithypercube.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
/** @brief Test and trial order-policy dataset, independently of field components. */
template <bool RuntimeTest, bool RuntimeTrial>
struct OrderCase
{
    /** Test degree is P2 in every case. */
    static constexpr int testOrder = RuntimeTest ? Dynamic : 2;
    /** Trial degree is P1, making rectangular local extents observable. */
    static constexpr int trialOrder = RuntimeTrial ? Dynamic : 1;
};
/** All four independent order-policy combinations share one assembly check. */
using OrderCases = boost::mp11::mp_list<OrderCase<false, false>, OrderCase<false, true>, OrderCase<true, false>, OrderCase<true, true>>;

/** @brief Construct scalar or vector discontinuous spaces with the same effective degree. */
template <int Order, bool Vector, typename Mesh>
auto makeSpace( std::shared_ptr<Mesh> const& mesh, int degree )
{
    if constexpr ( Vector )
    {
        if constexpr ( Order == Dynamic )
            return Pdhv<Order>( mesh, RuntimeOrder{ degree } );
        else
            return Pdhv<Order>( mesh );
    }
    else
    {
        if constexpr ( Order == Dynamic )
            return Pdh<Order>( mesh, RuntimeOrder{ degree } );
        else
            return Pdh<Order>( mesh );
    }
}

/** @brief Compare assembled element and four-block face actions with exact constant integrals. */
template <typename Case, bool Vector>
void checkAssembly()
{
    auto mesh = unitHypercube<2>( .7 );
    auto test = makeSpace<Case::testOrder, Vector>( mesh, 2 );
    auto trial = makeSpace<Case::trialOrder, Vector>( mesh, 1 );
    auto v = test->element();
    auto u = trial->element();
    auto constant = []
    {
        if constexpr ( Vector )
            return vec( cst( 1. ), cst( 2. ) );
        else
            return cst( 1. );
    }();
    v.on( _range = elements( mesh ), _expr = constant );
    u.on( _range = elements( mesh ), _expr = constant );
    double const scale = Vector ? 5. : 1.;
    double const volume = integrate( _range = elements( mesh ), _expr = cst( 1. ) ).evaluate()( 0, 0 );
    double const interface = integrate( _range = internalfaces( mesh ), _expr = cst( 1. ) ).evaluate()( 0, 0 );

    // Default component assembly exercises rectangular P2-by-P1 element matrices.
    auto mass = form2( _test = test, _trial = trial );
    mass = integrate( _range = elements( mesh ), _expr = inner( id( v ), idt( u ) ) );
    mass.close();
    BOOST_CHECK_SMALL( mass( v, u ) - scale * volume, 2e-10 );

    // Both traces contribute, forcing all four adjacent-element blocks to be assembled.
    auto face = form2( _test = test, _trial = trial, _pattern = size_type( Pattern::EXTENDED ) );
    face = integrate( _range = internalfaces( mesh ),
                      _expr = inner( leftface( id( v ) ) + rightface( id( v ) ), leftfacet( idt( u ) ) + rightfacet( idt( u ) ) ) );
    face.close();
    BOOST_CHECK_SMALL( face( v, u ) - 4 * scale * interface, 2e-9 );
}
} // namespace

BOOST_AUTO_TEST_CASE_TEMPLATE( scalarMixedAssembly, Case, OrderCases )
{
    checkAssembly<Case, false>();
}
BOOST_AUTO_TEST_CASE_TEMPLATE( componentMixedAssembly, Case, OrderCases )
{
    checkAssembly<Case, true>();
}
