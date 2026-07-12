/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file test_fe_runtime_contract.cpp
 * @brief FE-family compile-time/runtime-order contract checks.
 */
#define BOOST_TEST_MODULE test_fe_runtime_contract
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/doflayout.hpp>
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

/** @test Classify representative concrete families with the layered FE concepts. */
BOOST_AUTO_TEST_CASE( concrete_families_satisfy_one_layered_fe_contract )
{
    using h1_p1 = typename Lagrange<1, Scalar>::template apply<2, 2, double, Simplex<2>>::type;
    using l2_p1 = typename Lagrange<1, Scalar, Discontinuous>::template apply<2, 2, double, Simplex<2>>::type;
    using rt0 = typename RaviartThomas<0>::template apply<2, 2, double, Simplex<2>>::type;
    using bdm0 = typename BrezziDouglasMarini<0>::template apply<2, 2, double, Simplex<2>>::type;
    using ned0 = typename Nedelec<0, NedelecKind::NED1>::template apply<2, 2, double, Simplex<2>>::type;
    using cr1 = typename CrouzeixRaviart<1>::template apply<2, 2, double, Simplex<2>>::type;

    static_assert( CiarletFiniteElement<h1_p1> );
    static_assert( CiarletFiniteElement<l2_p1> );
    static_assert( CiarletFiniteElement<rt0> );
    static_assert( CiarletFiniteElement<bdm0> );
    static_assert( CiarletFiniteElement<ned0> );
    static_assert( CiarletFiniteElement<cr1> );

    static_assert( H1FiniteElement<h1_p1> );
    static_assert( L2FiniteElement<l2_p1> );
    static_assert( HDivFiniteElement<rt0> );
    static_assert( HDivFiniteElement<bdm0> );
    static_assert( HCurlFiniteElement<ned0> );
    static_assert( NonconformingH1FiniteElement<cr1> );

    static_assert( !H1FiniteElement<cr1> );
    static_assert( !HDivFiniteElement<ned0> );
    static_assert( !HCurlFiniteElement<rt0> );

    h1_p1 h1;
    rt0 rt;
    bdm0 bdm;
    ned0 ned;

    BOOST_CHECK_EQUAL( h1.polynomialDegree(), 1 );
    BOOST_CHECK_EQUAL( h1.localDof(), h1.localDofCount( true ) );
    BOOST_CHECK_EQUAL( rt.localDof(), rt.localDofCount( true ) );
    BOOST_CHECK_EQUAL( bdm.localDof(), bdm.localDofCount( true ) );
    BOOST_CHECK_EQUAL( ned.localDof(), ned.localDofCount( true ) );
}

/** @test Check zero-copy Eigen representations while retaining mathematical functional objects. */
BOOST_AUTO_TEST_CASE( functional_objects_expose_eigen_representations_without_losing_math_api )
{
    using fe_type = typename Lagrange<1, Scalar>::template apply<2, 2, double, Simplex<2>>::type;
    using primal_type = typename fe_type::primal_space_type;
    using functional_type = Functional<primal_type>;

    static_assert( !std::is_polymorphic_v<functional_type> );

    fe_type fe;
    auto functionals = functional::makePointEvaluationFunctionals( fe.primal(), fe.points() );
    FunctionalSet<primal_type> dual( fe.primal(), functionals );

    BOOST_REQUIRE_EQUAL( dual.size(), functionals.size() );
    BOOST_CHECK_EQUAL( dual.dualMatrix().size1(), functionals.size() );
    BOOST_CHECK_EQUAL( dual.dualMatrix().size2(), fe.primal().polynomialDimension() );

    auto const eigenDual = dual.eigenDualMatrix();
    BOOST_CHECK_EQUAL( eigenDual.rows(), dual.dualMatrix().size1() );
    BOOST_CHECK_EQUAL( eigenDual.cols(), dual.dualMatrix().size2() );
    BOOST_CHECK_EQUAL( eigenDual( 0, 0 ), dual.dualMatrix()( 0, 0 ) );

    auto const riesz = functionals.front().rieszRepresentation();
    BOOST_CHECK_EQUAL( riesz.rows(), 1 );
    BOOST_CHECK_EQUAL( riesz.cols(), fe.primal().polynomialDimension() );
}

BOOST_AUTO_TEST_SUITE_END()
