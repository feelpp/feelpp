/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
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

namespace
{
template<uint16_type Dim, uint16_type Order>
void checkNedelec1PrimalStaticRuntimeParity()
{
    NedelecPolynomialSet<Dim, Order, NedelecKind::NED1> staticPrimal;
    Nedelec1DynamicPolynomialSet<Dim> runtimePrimal( RuntimeOrder{ Order } );
    auto const& staticCoefficients = staticPrimal.coeff();
    auto const& runtimeCoefficients = runtimePrimal.coeff();
    BOOST_REQUIRE_EQUAL( staticCoefficients.size1(), runtimeCoefficients.size1() );
    BOOST_REQUIRE_EQUAL( staticCoefficients.size2(), runtimeCoefficients.size2() );
    for ( std::size_t i = 0; i < staticCoefficients.size1(); ++i )
        for ( std::size_t j = 0; j < staticCoefficients.size2(); ++j )
            BOOST_CHECK_SMALL( staticCoefficients( i, j )-runtimeCoefficients( i, j ), 1e-12 );
}

template<uint16_type Dim, uint16_type Order>
void checkNedelec1DualStaticRuntimeParity()
{
    using static_primal_type = NedelecPolynomialSet<Dim, Order, NedelecKind::NED1>;
    using static_dual_type = fem::detail::NedelecDualFirstKind<
        static_primal_type, PointSetEquiSpaced>;
    using runtime_primal_type = Nedelec1DynamicPolynomialSet<Dim>;
    using runtime_dual_type = fem::detail::Nedelec1DynamicDual<
        runtime_primal_type, PointSetEquiSpaced>;

    static_primal_type staticPrimal;
    static_dual_type staticDual( staticPrimal );
    runtime_primal_type runtimePrimal( RuntimeOrder{ Order } );
    runtime_dual_type runtimeDual( runtimePrimal );

    BOOST_CHECK_EQUAL( runtimeDual.publicOrder(), Order );
    BOOST_CHECK_EQUAL( runtimeDual.runtimeDofPerEdge(),
                       Feel::detail::nedelec1EdgeDof( Order ) );
    BOOST_CHECK_EQUAL( runtimeDual.runtimeDofPerFace(),
                       Feel::detail::nedelec1FaceDof( Dim, Order ) );
    BOOST_CHECK_EQUAL( runtimeDual.runtimeDofPerVolume(),
                       Dim == 3 ? Feel::detail::nedelec1VolumeDof( Dim, Order ) : 0 );
    BOOST_CHECK_EQUAL( runtimeDual.runtimeLocalDof(),
                       Feel::detail::nedelec1TotalDof( Dim, Order ) );

    auto const staticMatrix = staticDual( staticPrimal );
    auto const runtimeMatrix = runtimeDual( runtimePrimal );
    BOOST_REQUIRE_EQUAL( staticMatrix.size1(), runtimeMatrix.size1() );
    BOOST_REQUIRE_EQUAL( staticMatrix.size2(), runtimeMatrix.size2() );
    for ( std::size_t i = 0; i < staticMatrix.size1(); ++i )
        for ( std::size_t j = 0; j < staticMatrix.size2(); ++j )
            BOOST_CHECK_SMALL( staticMatrix( i, j )-runtimeMatrix( i, j ), 1e-11 );
}

template<uint16_type Dim, uint16_type Order,
         NedelecKind Kind = NedelecKind::NED1>
void checkNedelecFiniteElementStaticRuntimeParity()
{
    using static_fe_type = fem::NedelecElement<Dim, Order, Kind>;
    using runtime_fe_type = fem::NedelecElement<Dim, Dynamic, Kind>;
    static_fe_type staticFe;
    runtime_fe_type runtimeFe( RuntimeOrder{ Order } );

    BOOST_CHECK_EQUAL( runtimeFe.order(), Order );
    BOOST_CHECK_EQUAL( runtimeFe.polynomialDegree(), Order+1 );
    BOOST_CHECK_EQUAL( runtimeFe.localDofCount(), staticFe.localDofCount() );
    BOOST_CHECK_EQUAL( runtimeFe.localDofPerComponent(), staticFe.localDofPerComponent() );
    BOOST_CHECK_EQUAL( runtimeFe.dofPerEdge(), staticFe.nDofPerEdge );
    BOOST_CHECK_EQUAL( runtimeFe.dofPerFace(), staticFe.nDofPerFace );
    BOOST_CHECK_EQUAL( runtimeFe.dofPerVolume(), staticFe.nDofPerVolume );
    BOOST_CHECK_EQUAL( runtimeFe.points().size2(), staticFe.points().size2() );

    auto const& staticSampling = staticFe.dual().samplingLayout();
    auto const& runtimeSampling = runtimeFe.dual().samplingLayout();
    BOOST_CHECK_EQUAL( runtimeSampling.edgePointCount, staticSampling.edgePointCount );
    BOOST_CHECK_EQUAL( runtimeSampling.faces.size(), staticSampling.faces.size() );
    BOOST_CHECK_EQUAL( runtimeSampling.cellFirstDof, staticSampling.cellFirstDof );
    BOOST_CHECK_EQUAL( runtimeSampling.cellMomentCount, staticSampling.cellMomentCount );
    BOOST_CHECK_EQUAL( runtimeSampling.cellPointOffset, staticSampling.cellPointOffset );
    BOOST_CHECK_EQUAL( runtimeSampling.cellPointCount, staticSampling.cellPointCount );
    BOOST_CHECK_EQUAL( runtimeSampling.totalPointCount, staticSampling.totalPointCount );
    for ( std::size_t face = 0; face < staticSampling.faces.size(); ++face )
    {
        auto const& staticBlock = staticSampling.faces[face];
        auto const& runtimeBlock = runtimeSampling.faces[face];
        BOOST_CHECK_EQUAL( runtimeBlock.face, staticBlock.face );
        BOOST_CHECK_EQUAL( runtimeBlock.firstDof, staticBlock.firstDof );
        BOOST_CHECK_EQUAL( runtimeBlock.momentCount, staticBlock.momentCount );
        BOOST_CHECK_EQUAL( runtimeBlock.elementPointOffset, staticBlock.elementPointOffset );
        BOOST_CHECK_EQUAL( runtimeBlock.facetPointOffset, staticBlock.facetPointOffset );
        BOOST_CHECK_EQUAL( runtimeBlock.pointCount, staticBlock.pointCount );
    }

    auto const& staticCoefficients = staticFe.coeff();
    auto const& runtimeCoefficients = runtimeFe.coeff();
    BOOST_REQUIRE_EQUAL( staticCoefficients.size1(), runtimeCoefficients.size1() );
    BOOST_REQUIRE_EQUAL( staticCoefficients.size2(), runtimeCoefficients.size2() );
    for ( std::size_t i = 0; i < staticCoefficients.size1(); ++i )
        for ( std::size_t j = 0; j < staticCoefficients.size2(); ++j )
            BOOST_CHECK_SMALL( staticCoefficients( i, j )-runtimeCoefficients( i, j ), 1e-10 );

    for ( uint16_type localDof = 0; localDof < runtimeFe.localDofCount(); ++localDof )
    {
        auto const staticAttachment = staticFe.dofAttachment( localDof );
        auto const runtimeAttachment = runtimeFe.dofAttachment( localDof );
        BOOST_CHECK_EQUAL( runtimeAttachment.entityDim, staticAttachment.entityDim );
        BOOST_CHECK_EQUAL( runtimeAttachment.entityId, staticAttachment.entityId );
        BOOST_CHECK_EQUAL( runtimeAttachment.ordinal, staticAttachment.ordinal );
        BOOST_CHECK_EQUAL( runtimeFe.dofFunctionalKind( localDof ),
                           staticFe.dofFunctionalKind( localDof ) );
    }

    typename static_fe_type::points_type points( Dim, 3 );
    for ( uint16_type component = 0; component < Dim; ++component )
    {
        points( component, 0 ) = -0.4 + 0.1*component;
        points( component, 1 ) = -0.2 + 0.05*component;
        points( component, 2 ) = -0.1 - 0.05*component;
    }
    auto const staticDerivatives = staticFe.derivate( points );
    auto const runtimeDerivatives = runtimeFe.derivate( points );
    BOOST_REQUIRE_EQUAL( staticDerivatives.size(), runtimeDerivatives.size() );
    for ( std::size_t direction = 0; direction < staticDerivatives.size(); ++direction )
    {
        BOOST_REQUIRE_EQUAL( staticDerivatives[direction].size1(),
                             runtimeDerivatives[direction].size1() );
        BOOST_REQUIRE_EQUAL( staticDerivatives[direction].size2(),
                             runtimeDerivatives[direction].size2() );
        for ( std::size_t i = 0; i < staticDerivatives[direction].size1(); ++i )
            for ( std::size_t q = 0; q < staticDerivatives[direction].size2(); ++q )
                BOOST_CHECK_SMALL( staticDerivatives[direction]( i, q )-
                                   runtimeDerivatives[direction]( i, q ), 1e-10 );
    }
}
} // namespace

BOOST_AUTO_TEST_SUITE( fe_runtime_contract_suite )

BOOST_AUTO_TEST_CASE( nedelec1_dynamic_primal_matches_static_dimensions )
{
    for ( uint16_type order = 0; order <= 3; ++order )
    {
        Nedelec1DynamicPolynomialSet<2> primal2d( RuntimeOrder{ order } );
        Nedelec1DynamicPolynomialSet<3> primal3d( RuntimeOrder{ order } );
        BOOST_CHECK_EQUAL( primal2d.publicOrder(), order );
        BOOST_CHECK_EQUAL( primal3d.publicOrder(), order );
        BOOST_CHECK_EQUAL( primal2d.polynomialDimension(),
                           Feel::detail::nedelec1TotalDof( 2, order ) );
        BOOST_CHECK_EQUAL( primal3d.polynomialDimension(),
                           Feel::detail::nedelec1TotalDof( 3, order ) );
    }

    checkNedelec1PrimalStaticRuntimeParity<2,0>();
    checkNedelec1PrimalStaticRuntimeParity<2,1>();
    checkNedelec1PrimalStaticRuntimeParity<2,2>();
    checkNedelec1PrimalStaticRuntimeParity<2,3>();
    checkNedelec1PrimalStaticRuntimeParity<3,0>();
    checkNedelec1PrimalStaticRuntimeParity<3,1>();
    checkNedelec1PrimalStaticRuntimeParity<3,2>();
    checkNedelec1PrimalStaticRuntimeParity<3,3>();
}

/** @test Runtime NED1 duals reuse the static mathematical functional kernel. */
BOOST_AUTO_TEST_CASE( nedelec1_dynamic_dual_matches_static_functionals )
{
    checkNedelec1DualStaticRuntimeParity<2,0>();
    checkNedelec1DualStaticRuntimeParity<2,1>();
    checkNedelec1DualStaticRuntimeParity<2,2>();
    checkNedelec1DualStaticRuntimeParity<2,3>();
    checkNedelec1DualStaticRuntimeParity<3,0>();
    checkNedelec1DualStaticRuntimeParity<3,1>();
    checkNedelec1DualStaticRuntimeParity<3,2>();
    checkNedelec1DualStaticRuntimeParity<3,3>();
}

/** @test The public runtime FE shell is identical to static NED1 through order three. */
BOOST_AUTO_TEST_CASE( nedelec1_dynamic_finite_element_matches_static_shell )
{
    checkNedelecFiniteElementStaticRuntimeParity<2,0>();
    checkNedelecFiniteElementStaticRuntimeParity<2,1>();
    checkNedelecFiniteElementStaticRuntimeParity<2,2>();
    checkNedelecFiniteElementStaticRuntimeParity<2,3>();
    checkNedelecFiniteElementStaticRuntimeParity<3,0>();
    checkNedelecFiniteElementStaticRuntimeParity<3,1>();
    checkNedelecFiniteElementStaticRuntimeParity<3,2>();
    checkNedelecFiniteElementStaticRuntimeParity<3,3>();
}

/** @test NED2 uses the same true runtime FE shell as its static instances. */
BOOST_AUTO_TEST_CASE( nedelec2_dynamic_finite_element_matches_static_shell )
{
    checkNedelecFiniteElementStaticRuntimeParity<2,0,NedelecKind::NED2>();
    checkNedelecFiniteElementStaticRuntimeParity<2,1,NedelecKind::NED2>();
    checkNedelecFiniteElementStaticRuntimeParity<2,2,NedelecKind::NED2>();
    checkNedelecFiniteElementStaticRuntimeParity<2,3,NedelecKind::NED2>();
    checkNedelecFiniteElementStaticRuntimeParity<3,0,NedelecKind::NED2>();
    checkNedelecFiniteElementStaticRuntimeParity<3,1,NedelecKind::NED2>();
    checkNedelecFiniteElementStaticRuntimeParity<3,2,NedelecKind::NED2>();
    checkNedelecFiniteElementStaticRuntimeParity<3,3,NedelecKind::NED2>();
}

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

/** @test Bind static and runtime descriptors through the same public family shells. */
BOOST_AUTO_TEST_CASE( family_descriptors_bind_to_unified_shells )
{
    using rt_static = typename RaviartThomas<0>::template apply<2, 2, double, Simplex<2>>::type;
    using rt_dynamic = typename RaviartThomas<Dynamic>::template apply<2, 2, double, Simplex<2>>::type;
    using bdm_static = typename BrezziDouglasMarini<0>::template apply<2, 2, double, Simplex<2>>::type;
    using bdm_dynamic = typename BrezziDouglasMarini<Dynamic>::template apply<2, 2, double, Simplex<2>>::type;
    using ned_static = typename Nedelec<0, NedelecKind::NED1>::template apply<2, 2, double, Simplex<2>>::type;
    using ned_dynamic = typename Nedelec<Dynamic, NedelecKind::NED1>::template apply<2, 2, double, Simplex<2>>::type;

    static_assert( std::same_as<rt_static, fem::RaviartThomasElement<2, 0, double, Simplex, 0>> );
    static_assert( std::same_as<rt_dynamic, fem::RaviartThomasElement<2, Dynamic, double, Simplex, 0>> );
    static_assert( std::same_as<bdm_static, fem::BrezziDouglasMariniElement<2, 0, double, Simplex, 0>> );
    static_assert( std::same_as<bdm_dynamic, fem::BrezziDouglasMariniElement<2, Dynamic, double, Simplex, 0>> );
    static_assert( std::same_as<ned_static, fem::NedelecElement<2, 0, NedelecKind::NED1, double, 0>> );
    static_assert( std::same_as<ned_dynamic, fem::NedelecElement<2, Dynamic, NedelecKind::NED1, double, 0>> );

    static_assert( HDivFiniteElement<rt_static> && HDivFiniteElement<rt_dynamic> );
    static_assert( HDivFiniteElement<bdm_static> && HDivFiniteElement<bdm_dynamic> );
    static_assert( HCurlFiniteElement<ned_static> && HCurlFiniteElement<ned_dynamic> );

    ned_dynamic ned0( RuntimeOrder{ 0 } );
    ned_dynamic ned1( RuntimeOrder{ 1 } );
    ned_dynamic ned3( RuntimeOrder{ 3 } );
    BOOST_CHECK_EQUAL( ned0.order(), 0 );
    BOOST_CHECK_EQUAL( ned0.polynomialDegree(), 1 );
    BOOST_CHECK_EQUAL( ned1.order(), 1 );
    BOOST_CHECK_EQUAL( ned1.polynomialDegree(), 2 );
    BOOST_CHECK_EQUAL( ned3.order(), 3 );
    BOOST_CHECK_EQUAL( ned3.polynomialDegree(), 4 );
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
    static_assert( ReferenceFiniteElement<h1_p1> );
    static_assert( ReferenceFiniteElement<l2_p1> );
    static_assert( ReferenceFiniteElement<rt0> );
    static_assert( ReferenceFiniteElement<bdm0> );
    static_assert( ReferenceFiniteElement<ned0> );
    static_assert( ReferenceFiniteElement<cr1> );

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
