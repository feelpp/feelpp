/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2026-05-09

  Copyright (C) 2026 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
 * @file test_fe_runtime_contract.cpp
 * @brief FE-family compile-time/runtime-order contract checks.
 */
#define BOOST_TEST_MODULE test_fe_runtime_contract
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelpoly/brezzidouglasmarini.hpp>
#include <feel/feelpoly/crouzeixraviart.hpp>
#include <feel/feelpoly/lagrange.hpp>
#include <feel/feelpoly/nedelec.hpp>
#include <feel/feelpoly/raviartthomas.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( fe_runtime_contract_suite )

BOOST_AUTO_TEST_CASE( public_fe_families_expose_static_and_dynamic_order_traits )
{
    using lag_static = Lagrange<1, Scalar>;
    using lag_dynamic = Lagrange<Dynamic, Scalar>;
    using rt_static = RaviartThomas<0>;
    using rt_dynamic = RaviartThomas<Dynamic>;
    using ned_static = Nedelec<0, NedelecKind::NED1>;
    using ned_dynamic = Nedelec<Dynamic, NedelecKind::NED1>;
    using bdm_static = BrezziDouglasMarini<0>;
    using bdm_dynamic = BrezziDouglasMarini<Dynamic>;

    static_assert( lag_static::is_order_static );
    static_assert( !lag_static::is_order_dynamic );
    static_assert( lag_dynamic::is_order_dynamic );
    static_assert( !lag_dynamic::is_order_static );

    static_assert( rt_static::is_order_static );
    static_assert( !rt_static::is_order_dynamic );
    static_assert( rt_dynamic::is_order_dynamic );
    static_assert( !rt_dynamic::is_order_static );

    static_assert( ned_static::is_order_static );
    static_assert( !ned_static::is_order_dynamic );
    static_assert( ned_dynamic::is_order_dynamic );
    static_assert( !ned_dynamic::is_order_static );

    static_assert( bdm_static::is_order_static );
    static_assert( !bdm_static::is_order_dynamic );
    static_assert( bdm_dynamic::is_order_dynamic );
    static_assert( !bdm_dynamic::is_order_static );

    static_assert( !orderIsDynamic<lag_static> );
    static_assert( orderIsDynamic<lag_dynamic> );
    static_assert( !orderIsDynamic<rt_static> );
    static_assert( orderIsDynamic<rt_dynamic> );
    static_assert( !orderIsDynamic<ned_static> );
    static_assert( orderIsDynamic<ned_dynamic> );
    static_assert( !orderIsDynamic<bdm_static> );
    static_assert( orderIsDynamic<bdm_dynamic> );

    BOOST_CHECK_EQUAL( lag_static::nOrder, 1 );
    BOOST_CHECK_EQUAL( lag_dynamic::nOrder, 0 );
    BOOST_CHECK_EQUAL( rt_static::nOrder, 0 );
    BOOST_CHECK_EQUAL( rt_dynamic::nOrder, 0 );
    BOOST_CHECK_EQUAL( ned_static::nOrder, 0 );
    BOOST_CHECK_EQUAL( ned_dynamic::nOrder, 0 );
    BOOST_CHECK_EQUAL( bdm_static::nOrder, 0 );
    BOOST_CHECK_EQUAL( bdm_dynamic::nOrder, 0 );
}

BOOST_AUTO_TEST_CASE( dynamic_family_component_bases_preserve_dynamic_order )
{
    using lag_dynamic_component = typename Lagrange<Dynamic, Vectorial>::component_basis_type;
    using rt_dynamic_component = typename RaviartThomas<Dynamic>::component_basis_type;
    using ned_dynamic_component = typename Nedelec<Dynamic, NedelecKind::NED1>::component_basis_type;
    using bdm_dynamic_component = typename BrezziDouglasMarini<Dynamic>::component_basis_type;

    static_assert( lag_dynamic_component::is_order_dynamic );
    static_assert( rt_dynamic_component::is_order_dynamic );
    static_assert( ned_dynamic_component::is_order_dynamic );
    static_assert( bdm_dynamic_component::is_order_dynamic );
}

BOOST_AUTO_TEST_CASE( static_low_order_simplex_factories_still_apply )
{
    using lag_p1 = typename Lagrange<1, Scalar>::template apply<2, 2, double, Simplex<2>>::type;
    using rt0 = typename RaviartThomas<0>::template apply<2, 2, double, Simplex<2>>::type;
    using ned0 = typename Nedelec<0, NedelecKind::NED1>::template apply<2, 2, double, Simplex<2>>::type;
    using bdm0 = typename BrezziDouglasMarini<0>::template apply<2, 2, double, Simplex<2>>::type;
    using cr1 = typename CrouzeixRaviart<1>::template apply<2, 2, double, Simplex<2>>::type;

    static_assert( lag_p1::nDim == 2 );
    static_assert( lag_p1::is_scalar );
    static_assert( rt0::nDim == 2 );
    static_assert( rt0::is_vectorial );
    static_assert( ned0::nDim == 2 );
    static_assert( ned0::is_vectorial );
    static_assert( bdm0::nDim == 2 );
    static_assert( bdm0::is_vectorial );
    static_assert( cr1::nDim == 2 );
    static_assert( cr1::is_scalar );

    lag_p1 lag;
    rt0 rt;
    ned0 ned;
    bdm0 bdm;
    cr1 cr;

    BOOST_CHECK_EQUAL( lag.familyName(), "lagrange" );
    BOOST_CHECK_EQUAL( rt.familyName(), "raviartthomas" );
    BOOST_CHECK_EQUAL( ned.familyName(), "nedelec" );
    BOOST_CHECK_EQUAL( bdm.familyName(), "brezzidouglasmarini" );
    BOOST_CHECK_EQUAL( cr.familyName(), "CrouzeixRaviart" );
    BOOST_CHECK_EQUAL( cr.localDofCount(), cr.localDofPerComponent() );
    BOOST_CHECK_EQUAL( cr.localDofCountOnFacet( 0, true ), 1 );
}

BOOST_AUTO_TEST_SUITE_END()
