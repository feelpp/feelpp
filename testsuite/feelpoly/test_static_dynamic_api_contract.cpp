/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 * @file test_static_dynamic_api_contract.cpp
 * @brief Contract checks for a seamless static/dynamic order API.
 */
#define BOOST_TEST_MODULE test_static_dynamic_api_contract
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelmesh/geond.hpp>
#include <feel/feelpoly/brezzidouglasmarini.hpp>
#include <feel/feelpoly/broken.hpp>
#include <feel/feelpoly/concepts.hpp>
#include <feel/feelpoly/equispaced.hpp>
#include <feel/feelpoly/lagrange.hpp>
#include <feel/feelpoly/nedelec.hpp>
#include <feel/feelpoly/orthonormalpolynomialset.hpp>
#include <feel/feelpoly/raviartthomas.hpp>

#include <boost/mp11/list.hpp>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
/** @brief Check public modal metadata and its applied implementation together.
 * @tparam Order Static expansion degree or `Dynamic`.
 * @tparam Cell Simplex or hypercube reference cell.
 */
template<int Order, typename Cell>
void checkModalDescriptor()
{
    using descriptor_type = OrthonormalPolynomialSet<Order, Scalar, 7>;
    using modal_type = typename descriptor_type::template apply<2, 2, double, Cell>::type;
    using retagged_type = typename descriptor_type::template ChangeTag<9>::type;
    static_assert( descriptor_type::order_spec == modal_type::nOrder );
    static_assert( descriptor_type::is_order_dynamic == modal_type::is_order_dynamic );
    static_assert( descriptor_type::is_order_static == modal_type::is_order_static );
    static_assert( descriptor_type::TAG == modal_type::TAG );
    static_assert( descriptor_type::TAG == 7 && retagged_type::TAG == 9 );
    static_assert( retagged_type::order_spec == Order );
    static_assert( descriptor_type::component_basis_type::order_spec == Order );
    static_assert( std::is_base_of_v<FiniteElementFamilyDescriptor<Order, 7>, descriptor_type> );
    if constexpr ( Order != Dynamic )
        static_assert( staticOrderValue<descriptor_type>() == Order );

    for ( int degree : { 0, 1, 3 } )
    {
        if constexpr ( Order != Dynamic )
            if ( degree != Order )
                continue;
        auto const order = validatedOrder<descriptor_type>( RuntimeOrder{ degree } );
        auto modal = makeOrderAware<modal_type>( order );
        BOOST_TEST_CONTEXT( "degree=" << degree << " tensorProduct=" << Cell::is_hypercube )
        {
            BOOST_CHECK_EQUAL( modal.order(), degree );
            BOOST_CHECK_EQUAL( modal.localDof(), Cell::polyDims( degree ) );
            BOOST_CHECK_EQUAL( modal.coeff().size1(), modal.localDof() );
            BOOST_CHECK_EQUAL( modal.coeff().size2(), modal.localDof() );
        }
    }
}

template<int Order>
using lagrange_fe_type =
    typename Lagrange<Order, Scalar>::template apply<
        2, 2, double, Simplex<2>>::type;

template<int Order>
using rt_fe_type =
    typename RaviartThomas<Order>::template apply<
        2, 2, double, Simplex<2>>::type;

template<int Order>
using bdm_fe_type =
    typename BrezziDouglasMarini<Order>::template apply<
        2, 2, double, Simplex<2>>::type;

template<int Order>
using nedelec_fe_type =
    typename Nedelec<Order, NedelecKind::NED1>::template apply<
        2, 2, double, Simplex<2>>::type;

template<int Order>
using broken_lagrange_fe_type =
    typename Broken<Lagrange<Order, Scalar>>::template apply<
        2, 2, double, Simplex<2>>::type;

template<int Order>
using simplex_orthonormal_type =
    Feel::detail::OrthonormalPolynomialSet<
        2, Order, 2, Scalar, double, 0, Simplex>;

template<typename StaticFE, typename DynamicFE>
void checkCommonFeTopology( StaticFE const& staticFe,
                            DynamicFE const& dynamicFe )
{
    static_assert( EffectiveOrderAware<StaticFE> );
    static_assert( EffectiveOrderAware<DynamicFE> );

    BOOST_CHECK_EQUAL( staticFe.localDof(), dynamicFe.localDof() );
    for ( uint16_type entityDim = 0; entityDim <= 2; ++entityDim )
        BOOST_CHECK_EQUAL(
            staticFe.dofPerEntity( entityDim, 0 ),
            dynamicFe.dofPerEntity( entityDim, 0 ) );
}

template<HasDofTopology FE>
void checkNamedDofTopology2D( FE const& fe )
{
    static_assert( FE::nDim == 2 );

    BOOST_CHECK_EQUAL( fe.dofPerVertex(), fe.dofPerEntity( 0, 0 ) );
    BOOST_CHECK_EQUAL( fe.dofPerEdge(), fe.dofPerEntity( 1, 0 ) );
    BOOST_CHECK_EQUAL( fe.dofPerFace(), fe.dofPerEntity( 2, 0 ) );
    BOOST_CHECK_EQUAL( fe.dofPerVolume(), 0 );
}

template<typename FE, typename Factory>
void checkFiniteElementValueSemantics( Factory&& makeFe )
{
    std::optional<FE> survivingCopy;
    uint16_type expectedOrder = 0;
    uint16_type expectedLocalDof = 0;

    {
        FE source = makeFe();
        expectedOrder = source.order();
        expectedLocalDof = source.localDof();
        survivingCopy.emplace( source );
        BOOST_CHECK( &survivingCopy->primal() ==
                     &survivingCopy->dual().primalSpace() );

        FE assigned = makeFe();
        assigned = source;
        BOOST_CHECK( &assigned.primal() ==
                     &assigned.dual().primalSpace() );
        BOOST_CHECK_EQUAL( assigned.order(), expectedOrder );
        BOOST_CHECK_EQUAL( assigned.localDof(), expectedLocalDof );
    }

    BOOST_REQUIRE( survivingCopy.has_value() );
    BOOST_CHECK( &survivingCopy->primal() ==
                 &survivingCopy->dual().primalSpace() );
    BOOST_CHECK_EQUAL( survivingCopy->order(), expectedOrder );
    BOOST_CHECK_EQUAL( survivingCopy->localDof(), expectedLocalDof );

    std::vector<FE> copies;
    copies.reserve( 1 );
    copies.push_back( *survivingCopy );
    copies.push_back( *survivingCopy );
    for ( auto const& copy : copies )
    {
        BOOST_CHECK( &copy.primal() == &copy.dual().primalSpace() );
        BOOST_CHECK_EQUAL( copy.order(), expectedOrder );
        BOOST_CHECK_EQUAL( copy.localDof(), expectedLocalDof );
    }

    FE moved = std::move( copies.back() );
    BOOST_CHECK( &moved.primal() == &moved.dual().primalSpace() );
    BOOST_CHECK_EQUAL( moved.order(), expectedOrder );
    BOOST_CHECK_EQUAL( moved.localDof(), expectedLocalDof );
}
} // namespace

BOOST_AUTO_TEST_SUITE( static_dynamic_api_contract_suite )

/** @test Effective order and point cardinality use the same object API. */
BOOST_AUTO_TEST_CASE( geometry_and_pointsets_expose_common_accessors )
{
    using static_simplex = Simplex<2, 2>;
    using dynamic_simplex = Simplex<2, Dynamic>;
    using static_hypercube = Hypercube<2, 2>;
    using dynamic_hypercube = Hypercube<2, Dynamic>;

    static_assert( EffectiveOrderAware<static_simplex> );
    static_assert( EffectiveOrderAware<dynamic_simplex> );
    static_assert( EffectiveOrderAware<static_hypercube> );
    static_assert( EffectiveOrderAware<dynamic_hypercube> );
    static_assert( HasStaticOrderAccessor<static_simplex> );
    static_assert( !HasStaticOrderAccessor<dynamic_simplex> );
    static_assert( HasStaticOrderAccessor<static_hypercube> );
    static_assert( !HasStaticOrderAccessor<dynamic_hypercube> );

    static_simplex staticSimplex;
    dynamic_simplex dynamicSimplex( RuntimeOrder{ 2 } );
    static_hypercube staticHypercube;
    dynamic_hypercube dynamicHypercube( RuntimeOrder{ 2 } );
    BOOST_CHECK_EQUAL( staticSimplex.order(), dynamicSimplex.order() );
    BOOST_CHECK_EQUAL( staticHypercube.order(), dynamicHypercube.order() );

    using static_simplex_entity = GeoND<2, static_simplex>;
    using dynamic_simplex_entity = GeoND<2, dynamic_simplex>;
    using static_hypercube_entity = GeoND<2, static_hypercube>;
    using dynamic_hypercube_entity = GeoND<2, dynamic_hypercube>;

    static_assert( HasPointCardinality<static_simplex_entity> );
    static_assert( HasPointCardinality<dynamic_simplex_entity> );
    static_assert( HasPointCardinality<static_hypercube_entity> );
    static_assert( HasPointCardinality<dynamic_hypercube_entity> );

    static_simplex_entity staticSimplexEntity;
    dynamic_simplex_entity dynamicSimplexEntity( RuntimeOrder{ 2 } );
    static_hypercube_entity staticHypercubeEntity;
    dynamic_hypercube_entity dynamicHypercubeEntity( RuntimeOrder{ 2 } );
    BOOST_CHECK_EQUAL(
        staticSimplexEntity.nPoints(), dynamicSimplexEntity.nPoints() );
    BOOST_CHECK_EQUAL(
        staticHypercubeEntity.nPoints(), dynamicHypercubeEntity.nPoints() );

    using static_pointset = PointSetEquiSpaced<Simplex<2>, 2, double>;
    using dynamic_pointset = PointSetEquiSpaced<Simplex<2>, Dynamic, double>;
    static_assert( EffectiveOrderAware<static_pointset> );
    static_assert( EffectiveOrderAware<dynamic_pointset> );
    static_assert( HasPointCardinality<static_pointset> );
    static_assert( HasPointCardinality<dynamic_pointset> );
    static_assert( HasStaticOrderAccessor<static_pointset> );
    static_assert( !HasStaticOrderAccessor<dynamic_pointset> );

    static_pointset staticPoints;
    dynamic_pointset dynamicPoints( RuntimeOrder{ 2 } );
    BOOST_CHECK_EQUAL( staticPoints.order(), dynamicPoints.order() );
    BOOST_CHECK_EQUAL( staticPoints.nPoints(), dynamicPoints.nPoints() );
}

/** @test Record the current common-topology contract across FE families. */
BOOST_AUTO_TEST_CASE( finite_elements_share_effective_order_and_layout_accessors )
{
    using lagrange_static = lagrange_fe_type<2>;
    using lagrange_dynamic = lagrange_fe_type<Dynamic>;
    using rt_static = rt_fe_type<0>;
    using rt_dynamic = rt_fe_type<Dynamic>;
    using bdm_static = bdm_fe_type<0>;
    using bdm_dynamic = bdm_fe_type<Dynamic>;
    using nedelec_static = nedelec_fe_type<0>;
    using nedelec_dynamic = nedelec_fe_type<Dynamic>;

    static_assert( EffectiveOrderAware<lagrange_static> );
    static_assert( EffectiveOrderAware<lagrange_dynamic> );
    static_assert( EffectiveOrderAware<rt_static> );
    static_assert( EffectiveOrderAware<rt_dynamic> );
    static_assert( EffectiveOrderAware<bdm_static> );
    static_assert( EffectiveOrderAware<bdm_dynamic> );
    static_assert( EffectiveOrderAware<nedelec_static> );
    static_assert( EffectiveOrderAware<nedelec_dynamic> );

    static_assert( HasDofTopology<lagrange_static> );
    static_assert( HasDofTopology<lagrange_dynamic> );
    static_assert( HasDofTopology<rt_static> );
    static_assert( HasDofTopology<rt_dynamic> );
    static_assert( HasDofTopology<bdm_static> );
    static_assert( HasDofTopology<bdm_dynamic> );
    static_assert( HasDofTopology<nedelec_static> );
    static_assert( HasDofTopology<nedelec_dynamic> );

    lagrange_static staticLagrange;
    lagrange_dynamic dynamicLagrange( RuntimeOrder{ 2 } );
    BOOST_CHECK_EQUAL( staticLagrange.order(), dynamicLagrange.order() );
    checkCommonFeTopology( staticLagrange, dynamicLagrange );
    checkNamedDofTopology2D( staticLagrange );
    checkNamedDofTopology2D( dynamicLagrange );

    rt_static staticRt;
    rt_dynamic dynamicRt( RuntimeOrder{ 0 } );
    BOOST_CHECK_EQUAL( staticRt.order(), dynamicRt.order() );
    checkCommonFeTopology( staticRt, dynamicRt );
    checkNamedDofTopology2D( staticRt );
    checkNamedDofTopology2D( dynamicRt );

    bdm_static staticBdm;
    bdm_dynamic dynamicBdm( RuntimeOrder{ 0 } );
    BOOST_CHECK_EQUAL( staticBdm.order(), dynamicBdm.order() );
    checkCommonFeTopology( staticBdm, dynamicBdm );
    checkNamedDofTopology2D( staticBdm );
    checkNamedDofTopology2D( dynamicBdm );

    nedelec_static staticNedelec;
    nedelec_dynamic dynamicNedelec( RuntimeOrder{ 0 } );
    BOOST_CHECK_EQUAL( staticNedelec.order(), dynamicNedelec.order() );
    checkCommonFeTopology( staticNedelec, dynamicNedelec );
    checkNamedDofTopology2D( staticNedelec );
    checkNamedDofTopology2D( dynamicNedelec );
}

/** Static degrees and the runtime-order policy applied to both reference cells. */
using modal_order_specs = boost::mp11::mp_list_c<int, 0, 1, 3, Dynamic>;

/** @test Public modal descriptors preserve their order policy through apply/rebind. */
BOOST_AUTO_TEST_CASE_TEMPLATE( modal_descriptors_preserve_order_metadata, Order, modal_order_specs )
{
    checkModalDescriptor<Order::value, Simplex<2>>();
    checkModalDescriptor<Order::value, Hypercube<2>>();
}

/** @test Modal and broken producers expose the same topology vocabulary. */
BOOST_AUTO_TEST_CASE( orthonormal_and_broken_producers_share_the_contract )
{
    using orthonormal_static = simplex_orthonormal_type<2>;
    using orthonormal_dynamic = simplex_orthonormal_type<Dynamic>;
    using broken_static = broken_lagrange_fe_type<2>;
    using broken_dynamic = broken_lagrange_fe_type<Dynamic>;

    static_assert( HasDofTopology<orthonormal_static> );
    static_assert( HasDofTopology<orthonormal_dynamic> );
    static_assert( HasDofTopology<broken_static> );
    static_assert( HasDofTopology<broken_dynamic> );

    orthonormal_static staticModal;
    orthonormal_dynamic dynamicModal( RuntimeOrder{ 2 } );
    BOOST_CHECK_EQUAL( staticModal.order(), dynamicModal.order() );
    checkCommonFeTopology( staticModal, dynamicModal );
    checkNamedDofTopology2D( staticModal );
    checkNamedDofTopology2D( dynamicModal );

    broken_static staticBroken;
    broken_dynamic dynamicBroken( RuntimeOrder{ 2 } );
    BOOST_CHECK_EQUAL( staticBroken.order(), dynamicBroken.order() );
    checkCommonFeTopology( staticBroken, dynamicBroken );
    checkNamedDofTopology2D( staticBroken );
    checkNamedDofTopology2D( dynamicBroken );
    BOOST_CHECK_EQUAL( staticBroken.dofPerFace(), staticBroken.localDof() );
    BOOST_CHECK_EQUAL( dynamicBroken.dofPerFace(), dynamicBroken.localDof() );
}

/** @test Enforce zero-overhead static order and bounded dynamic storage. */
BOOST_AUTO_TEST_CASE( order_storage_layout_contract )
{
    static_assert( std::is_empty_v<OrderBase<2>> );
    static_assert( sizeof( OrderBase<Dynamic> ) == sizeof( RuntimeOrder ) );
    BOOST_CHECK( std::is_empty_v<OrderBase<2>> );
    BOOST_CHECK_EQUAL(
        sizeof( OrderBase<Dynamic> ), sizeof( RuntimeOrder ) );

    using static_shape = Simplex<2, 2>;
    using dynamic_shape = Simplex<2, Dynamic>;
    using static_entity = GeoND<2, static_shape>;
    using dynamic_entity = GeoND<2, dynamic_shape>;
    using static_hypercube_shape = Hypercube<2, 2>;
    using dynamic_hypercube_shape = Hypercube<2, Dynamic>;
    using static_hypercube_entity = GeoND<2, static_hypercube_shape>;
    using dynamic_hypercube_entity = GeoND<2, dynamic_hypercube_shape>;

    static_assert( static_shape::numPoints_static == 6 );
    static_assert( static_hypercube_shape::numPoints_static == 9 );
    static_assert( sizeof( static_shape ) <= sizeof( dynamic_shape ) );
    static_assert(
        sizeof( static_hypercube_shape ) <= sizeof( dynamic_hypercube_shape ) );
    static_assert( sizeof( static_entity ) == sizeof( dynamic_entity ) );
    static_assert(
        sizeof( static_hypercube_entity ) == sizeof( dynamic_hypercube_entity ) );
    using basis_type = typename lagrange_fe_type<2>::basis_type;
    static_assert( std::is_empty_v<Feel::detail::PolynomialBasisState<basis_type,2>> );
    static_assert( !std::is_empty_v<Feel::detail::PolynomialBasisState<basis_type,Dynamic>> );
    static_assert( orderIsStatic<typename lagrange_fe_type<2>::pointset_type> );
    static_assert( orderIsDynamic<typename lagrange_fe_type<Dynamic>::pointset_type> );
    // Whole FE size includes the point backend's reference geometry. It does
    // not measure order-policy overhead or total heap allocation footprint.

    BOOST_TEST_MESSAGE(
        "Static/dynamic API layout: OrderBase<2>=" << sizeof( OrderBase<2> )
        << ", OrderBase<Dynamic>=" << sizeof( OrderBase<Dynamic> )
        << ", Simplex<2,2>=" << sizeof( static_shape )
        << ", Simplex<2,Dynamic>=" << sizeof( dynamic_shape )
        << ", GeoND static=" << sizeof( static_entity )
        << ", GeoND dynamic=" << sizeof( dynamic_entity )
        << ", Hypercube GeoND static=" << sizeof( static_hypercube_entity )
        << ", Hypercube GeoND dynamic=" << sizeof( dynamic_hypercube_entity )
        << ", Lagrange static=" << sizeof( lagrange_fe_type<2> )
        << ", Lagrange dynamic=" << sizeof( lagrange_fe_type<Dynamic> ) );
}

BOOST_AUTO_TEST_CASE( finite_element_copies_own_their_primal_space )
{
    checkFiniteElementValueSemantics<lagrange_fe_type<2>>(
        [] { return lagrange_fe_type<2>{}; } );
    checkFiniteElementValueSemantics<lagrange_fe_type<Dynamic>>(
        [] { return lagrange_fe_type<Dynamic>{ RuntimeOrder{ 2 } }; } );
    checkFiniteElementValueSemantics<rt_fe_type<0>>(
        [] { return rt_fe_type<0>{}; } );
    checkFiniteElementValueSemantics<rt_fe_type<Dynamic>>(
        [] { return rt_fe_type<Dynamic>{ RuntimeOrder{ 0 } }; } );
    checkFiniteElementValueSemantics<bdm_fe_type<0>>(
        [] { return bdm_fe_type<0>{}; } );
    checkFiniteElementValueSemantics<bdm_fe_type<Dynamic>>(
        [] { return bdm_fe_type<Dynamic>{ RuntimeOrder{ 0 } }; } );
    checkFiniteElementValueSemantics<nedelec_fe_type<0>>(
        [] { return nedelec_fe_type<0>{}; } );
    checkFiniteElementValueSemantics<nedelec_fe_type<Dynamic>>(
        [] { return nedelec_fe_type<Dynamic>{ RuntimeOrder{ 0 } }; } );
}

BOOST_AUTO_TEST_CASE( runtime_fe_orders_fail_before_local_index_overflow )
{
    using lagrange_dynamic = lagrange_fe_type<Dynamic>;
    using rt_dynamic = rt_fe_type<Dynamic>;
    using bdm_dynamic = bdm_fe_type<Dynamic>;
    using nedelec_dynamic = nedelec_fe_type<Dynamic>;

    BOOST_CHECK_THROW( lagrange_dynamic( RuntimeOrder{ 65535 } ),
                       std::length_error );
    BOOST_CHECK_THROW( rt_dynamic( RuntimeOrder{ 65535 } ),
                       std::length_error );
    BOOST_CHECK_THROW( bdm_dynamic( RuntimeOrder{ 65535 } ),
                       std::length_error );
    BOOST_CHECK_THROW( nedelec_dynamic( RuntimeOrder{ 65535 } ),
                       std::length_error );
}

BOOST_AUTO_TEST_SUITE_END()
