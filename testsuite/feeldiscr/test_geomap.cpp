/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-02-07

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007-2010 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2010-2014 Feel++ Consortium

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
   @file test_geomap.cpp
   @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   @date 2005-02-07
   @brief GeoMap tests including dynamic order support (C++20/23 modernized)
 */
#define BOOST_TEST_MODULE geomap testsuite
#include <feel/feelcore/testsuite.hpp>

#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/debug.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/importergmsh.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelfilters/gmshsimplexdomain.hpp>
#include <feel/feelpoly/geomap.hpp>
#include <feel/feelpoly/order.hpp>
#include <feel/feelmesh/regiontree.hpp>

#include <concepts>
#include <type_traits>
#include <array>
#include <cmath>
#include <vector>

using namespace Feel;
namespace bdata = boost::unit_test::data;

//=============================================================================
// Modern C++20 helper concepts
//=============================================================================

/// @brief Concept for GeoMap-like types
template <typename T>
concept GeoMapLike = requires( T gm ) {
    { gm.dim() } -> std::convertible_to<uint16_type>;
    { gm.realDim() } -> std::convertible_to<uint16_type>;
    { gm.order() } -> std::convertible_to<uint16_type>;
    { T::is_order_static } -> std::convertible_to<bool>;
    { T::is_order_dynamic } -> std::convertible_to<bool>;
};

//=============================================================================
// Helper functions
//=============================================================================

inline auto f( node<double>::type const& __n ) -> double
{
    return ublas::sum( __n );
}

inline auto fx( node<double>::type const& __n ) -> double
{
    return __n[0];
}

template<int Order>
void compareGeomapContextStaticDynamicSimplex2D()
{
    using gm_static_type = GeoMap<2, Order, 2, double, Simplex>;
    using gm_dynamic_type = GeoMap<2, Dynamic, 2, double, Simplex>;
    using shape_type = Simplex<2, Order, 2>;
    using element_type = GeoND<2, shape_type>;
    using point_type = typename element_type::point_type;
    using node_type = typename point_type::node_type;

    constexpr double tol = 1e-12;
    constexpr size_type gmc_context_v = vm::POINT | vm::JACOBIAN | vm::KB;

    auto gm_static = std::make_shared<gm_static_type>();
    auto gm_dynamic = std::make_shared<gm_dynamic_type>( RuntimeOrder{ Order } );

    BOOST_REQUIRE_EQUAL( gm_static->order(), Order );
    BOOST_REQUIRE_EQUAL( gm_dynamic->order(), Order );

    // Runtime and compile-time behavior should agree on linear/nonlinear mapping.
    BOOST_CHECK_EQUAL( gm_static->isLinear(), Order == 1 );
    BOOST_CHECK_EQUAL( gm_dynamic->isLinear(), Order == 1 );

    element_type element;
    std::vector<point_type> point_storage;
    point_storage.reserve( shape_type::numPoints );

    for ( uint16_type i = 0; i < shape_type::numPoints; ++i )
    {
        auto ref = gm_static->refNode( i );
        const double x = ref( 0 );
        const double y = ref( 1 );

        node_type real_pt( 2 );
        // Intentionally curved map so high-order geometric nodes matter.
        real_pt( 0 ) = x + 0.17 * x * y + 0.11 * y * y;
        real_pt( 1 ) = y + 0.13 * x * x - 0.07 * x * y;

        point_storage.emplace_back( static_cast<uint32_type>( i ),
                                    real_pt,
                                    false,
                                    i < shape_type::numVertices );
        element.setPoint( i, point_storage.back() );
    }

    auto ref_pts = gm_static->points();
    auto pc_static = gm_static->preCompute( gm_static, ref_pts );
    auto pc_dynamic = gm_dynamic->preCompute( gm_dynamic, ref_pts );

    auto gmc_static = gm_static->template context<gmc_context_v>( element, pc_static );
    auto gmc_dynamic = gm_dynamic->template context<gmc_context_v>( element, pc_dynamic );

    BOOST_REQUIRE( gmc_static );
    BOOST_REQUIRE( gmc_dynamic );

    const uint16_type npts = ref_pts.size2();
    double max_err_x = 0.0;
    double max_err_J = 0.0;
    double max_err_K = 0.0;
    double max_err_B = 0.0;
    for ( uint16_type q = 0; q < npts; ++q )
    {
        auto x_static = gmc_static->xReal( q );
        auto x_dynamic = gmc_dynamic->xReal( q );
        for ( int c = 0; c < 2; ++c )
            max_err_x = std::max( max_err_x, std::abs( x_static( c ) - x_dynamic( c ) ) );

        max_err_J = std::max( max_err_J, std::abs( gmc_static->J( q ) - gmc_dynamic->J( q ) ) );

        auto const& K_static = gmc_static->K( q );
        auto const& K_dynamic = gmc_dynamic->K( q );
        auto const& B_static = gmc_static->B( q );
        auto const& B_dynamic = gmc_dynamic->B( q );

        for ( int i = 0; i < 2; ++i )
        {
            for ( int j = 0; j < 2; ++j )
            {
                max_err_K = std::max( max_err_K, std::abs( K_static( i, j ) - K_dynamic( i, j ) ) );
                max_err_B = std::max( max_err_B, std::abs( B_static( i, j ) - B_dynamic( i, j ) ) );
            }
        }
    }

    BOOST_CHECK_SMALL( max_err_x, tol );
    BOOST_CHECK_SMALL( max_err_J, tol );
    BOOST_CHECK_SMALL( max_err_K, tol );
    BOOST_CHECK_SMALL( max_err_B, tol );
}

template<int Order>
void compareGeomapHessianStaticDynamicSimplex2D()
{
    using gm_static_type = GeoMap<2, Order, 2, double, Simplex>;
    using gm_dynamic_type = GeoMap<2, Dynamic, 2, double, Simplex>;

    constexpr double tol = 1e-12;

    auto gm_static = std::make_shared<gm_static_type>();
    auto gm_dynamic = std::make_shared<gm_dynamic_type>( RuntimeOrder{ Order } );

    auto ref_pts = gm_static->points();
    auto pc_static = gm_static->preCompute( gm_static, ref_pts );
    auto pc_dynamic = gm_dynamic->preCompute( gm_dynamic, ref_pts );

    const uint16_type n_nodes = gm_static->points().size2();
    typename gm_static_type::matrix_node_t_type G( 2, n_nodes );
    for ( uint16_type i = 0; i < n_nodes; ++i )
    {
        auto ref = gm_static->refNode( i );
        const double x = ref( 0 );
        const double y = ref( 1 );
        G( 0, i ) = x + 0.17 * x * y + 0.11 * y * y;
        G( 1, i ) = y + 0.13 * x * x - 0.07 * x * y;
    }

    for ( uint16_type q = 0; q < ref_pts.size2(); ++q )
    {
        typename gm_static_type::hessian_basis_type h_basis_static( n_nodes, 2, 2 );
        typename gm_dynamic_type::hessian_basis_type h_basis_dynamic( n_nodes, 2, 2 );
        h_basis_static.setZero();
        h_basis_dynamic.setZero();

        gm_static->hessianBasisAtPoint( q, h_basis_static, pc_static.get() );
        gm_dynamic->hessianBasisAtPoint( q, h_basis_dynamic, pc_dynamic.get() );

        double max_err_h = 0.0;
        for ( int a = 0; a < 2; ++a )
        {
            for ( int b = 0; b < 2; ++b )
            {
                for ( int c = 0; c < 2; ++c )
                {
                    double h_static = 0.0;
                    double h_dynamic = 0.0;
                    for ( uint16_type i = 0; i < n_nodes; ++i )
                    {
                        h_static += G( a, i ) * h_basis_static( i, b, c );
                        h_dynamic += G( a, i ) * h_basis_dynamic( i, b, c );
                    }
                    max_err_h = std::max( max_err_h, std::abs( h_static - h_dynamic ) );
                }
            }
        }
        BOOST_CHECK_SMALL( max_err_h, tol );
    }
}

//=============================================================================
// Unified TestInterp struct for both static and dynamic order
//=============================================================================

/// @brief Test interpolation on mesh elements - unified for static and dynamic order
/// @tparam Dim Spatial dimension
/// @tparam Order Polynomial order (compile-time, or Dynamic for runtime)
/// @tparam Entity Entity template (Simplex or Hypercube)
template<int Dim, int Order, template<int, int, int> class Entity>
struct TestInterp
{
    static constexpr bool is_order_static = (Order != Dynamic);
    static constexpr bool is_order_dynamic = !is_order_static;

    // Entity type - same for both static and dynamic
    using entity_type = Entity<Dim, Order, Dim>;
    using geomap_type = GeoMap<Dim, Order, Dim, double, Entity>;

    // Mesh types only available for static order (mesh infra doesn't support dynamic yet)
    // For static order, we use these; for dynamic, we test GeoMap directly
    using mesh_type = std::conditional_t<is_order_static,
                                          Mesh<Entity<Dim, Order, Dim>>,
                                          void>;
    using mesh_ptr_type = std::conditional_t<is_order_static,
                                              std::shared_ptr<mesh_type>,
                                              void>;

    TestInterp() = default;

    /// @brief Test with compile-time order (full mesh test)
    void test( double hsize, std::string version = FEELPP_GMSH_FORMAT_VERSION ) requires( is_order_static )
    {
        testStaticOrder( hsize, version );
    }

    /// @brief Test with runtime order
    void test( uint16_type runtime_order, double hsize, std::string version = FEELPP_GMSH_FORMAT_VERSION ) requires( is_order_dynamic )
    {
        testDynamicOrder( runtime_order, hsize, version );
    }

private:
    //-------------------------------------------------------------------------
    // Static order implementation - full mesh test
    //-------------------------------------------------------------------------
    void testStaticOrder( double hsize, std::string version ) requires( is_order_static )
    {
        using ref_entity_type = Reference<entity_type, Dim, Order, Dim>;
        using gm_type = typename mesh_type::gm_type;
        using size_type = typename mesh_type::size_type;
        static constexpr size_type gmc_context_v = vm::POINT | vm::JACOBIAN | vm::HESSIAN;
        using gmc_type = typename gm_type::template Context<typename mesh_type::element_type>;
        using gmc_ptrtype = std::shared_ptr<gmc_type>;
        using gic_type = typename gm_type::Inverse;

        auto M_mesh = std::make_shared<mesh_type>();
        VLOG( 1 ) << "testing TestInterp<" << Dim << "," << Order << "> (static) with file format version " << version << "\n";

        GmshSimplexDomain td( entity_type::nDim, entity_type::nOrder );
        td.setVersion( version );
        td.setCharacteristicLength( hsize );
        auto fname = td.generate( entity_type::name().c_str() );
        ImporterGmsh<mesh_type> import( fname );
        import.setVersion( version );
        M_mesh->accept( import );

        // Verify GeoMap is_order_static
        static_assert( gm_type::is_order_static, "Static test should use static order GeoMap" );
        BOOST_CHECK_EQUAL( static_cast<int>( gm_type::nOrder ), Order );

        auto rangeElement = elements( *M_mesh );
        auto el_it = rangeElement.begin();
        auto el_en = rangeElement.end();

        ref_entity_type refelem;
        auto __geopc = std::make_shared<typename gm_type::precompute_type>(
            M_mesh->gm(), refelem.points() );

        MeshInverse<mesh_type> meshinv( M_mesh );
        meshinv.addPoints( M_mesh->points() );
        meshinv.distribute();

        std::vector<boost::tuple<size_type, uint16_type>> itab;

        std::cout << "refelem = " << refelem.points() << "\n";
        gmc_ptrtype gmc;

        for ( ; el_it != el_en; ++el_it )
        {
            auto const& meshElt = unwrap_ref( *el_it );
            if ( !gmc )
                gmc = M_mesh->gm()->template context<gmc_context_v>( meshElt, __geopc );
            else
                gmc->template update<gmc_context_v>( meshElt );
            gic_type gic( M_mesh->gm(), meshElt );

            meshinv.pointsInConvex( meshElt.id(), itab );

            for ( auto q = 0uz; q < itab.size(); ++q )
            {
                std::cout << "xref = " << meshinv.referenceCoords().find( boost::get<0>( itab[q] ) )->second << "\n";
            }

            for ( auto q = 0uz; q < refelem.points().size2(); ++q )
            {
                std::cout << "gmc xref " << q << " = " << gmc->xRef( q ) << "\n";
                std::cout << "is in gmc? = " << gmc->geometricMapping()->isIn( gmc->xRef( q ) ) << "\n";

                gic.setXReal( gmc->xReal( q ) );

                typename ref_entity_type::points_type pts( Dim, 1 );
                ublas::column( pts, 0 ) = gic.xRef();
                std::cout << "gic xref " << q << " = " << gic.xRef() << "\n";
                std::cout << "is in gic? = " << gic.geometricMapping()->isIn( gic.xRef() ) << "\n";
            }

            FEELPP_ASSERT( gic.isIn() )
            ( refelem.points() )( gmc->xReal() )
            ( meshElt.id() ).error( "invalid geometric transformation inversion" );
        }

        VLOG( 1 ) << "testing TestInterp (static) with file format version " << version << " done\n";
    }

    //-------------------------------------------------------------------------
    // Dynamic order implementation - GeoMap and entity test
    //-------------------------------------------------------------------------
    void testDynamicOrder( uint16_type runtime_order, double /*hsize*/, std::string /*version*/ ) requires( is_order_dynamic )
    {
        VLOG( 1 ) << "testing TestInterp<" << Dim << ", Dynamic> with runtime order " << runtime_order << "\n";

        // Verify compile-time dynamic detection
        static_assert( geomap_type::is_order_dynamic, "Dynamic test should use dynamic order GeoMap" );
        static_assert( entity_type::is_order_dynamic, "Dynamic test should use dynamic order Entity" );

        // Create GeoMap with runtime order
        geomap_type gm{ RuntimeOrder{ runtime_order } };

        // Verify runtime order
        BOOST_CHECK_EQUAL( gm.order(), runtime_order );
        BOOST_CHECK_EQUAL( gm.dim(), Dim );
        BOOST_CHECK( gm.is_order_dynamic );

        // Create entity type with dynamic order
        entity_type entity{ RuntimeOrder{ runtime_order } };
        BOOST_CHECK_EQUAL( entity.order(), runtime_order );
        BOOST_CHECK( entity.is_order_dynamic );

        // Verify total points match expected values
        if constexpr ( Dim == 1 )
        {
            BOOST_CHECK_EQUAL( entity.nPointsTotal(), runtime_order + 1 );
        }
        else
        {
            BOOST_CHECK( entity.nPointsTotal() > 0 );
        }

        VLOG( 1 ) << "testing TestInterp (dynamic) with runtime order " << runtime_order << " done\n";
    }
};


FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( test_geomap_suite )

//=============================================================================
// SECTION 1: Static Order Tests (full mesh)
//=============================================================================

BOOST_AUTO_TEST_CASE( test_geomap_static_order1 )
{
    TestInterp<2, 1, Simplex> test_interp;
    test_interp.test( doption( _name = "gmsh.hsize" ), FEELPP_GMSH_FORMAT_VERSION );
}

BOOST_AUTO_TEST_CASE( test_geomap_static_order2 )
{
    TestInterp<2, 2, Simplex> test_interp;
    test_interp.test( doption( _name = "gmsh.hsize" ), FEELPP_GMSH_FORMAT_VERSION );
}

//=============================================================================
// SECTION 2: Dynamic Order Tests
//=============================================================================

constexpr std::array test_orders = { 1, 2, 3 };

BOOST_DATA_TEST_CASE( test_geomap_dynamic_simplex_2d, bdata::make( test_orders ), order )
{
    TestInterp<2, Dynamic, Simplex> test;
    test.test( static_cast<uint16_type>( order ), doption( _name = "gmsh.hsize" ), FEELPP_GMSH_FORMAT_VERSION );
}

BOOST_DATA_TEST_CASE( test_geomap_dynamic_simplex_3d, bdata::make( test_orders ), order )
{
    TestInterp<3, Dynamic, Simplex> test;
    test.test( static_cast<uint16_type>( order ), doption( _name = "gmsh.hsize" ), FEELPP_GMSH_FORMAT_VERSION );
}

BOOST_DATA_TEST_CASE( test_geomap_dynamic_hypercube_2d, bdata::make( test_orders ), order )
{
    TestInterp<2, Dynamic, Hypercube> test;
    test.test( static_cast<uint16_type>( order ), doption( _name = "gmsh.hsize" ), FEELPP_GMSH_FORMAT_VERSION );
}

BOOST_DATA_TEST_CASE( test_geomap_dynamic_hypercube_3d, bdata::make( test_orders ), order )
{
    TestInterp<3, Dynamic, Hypercube> test;
    test.test( static_cast<uint16_type>( order ), doption( _name = "gmsh.hsize" ), FEELPP_GMSH_FORMAT_VERSION );
}

//=============================================================================
// SECTION 3: GeoMap Type Traits Tests
//=============================================================================

BOOST_AUTO_TEST_CASE( test_geomap_static_order_flags )
{
    using geomap_static = GeoMap<2, 1, 2, double, Simplex>;

    static_assert( GeoMapLike<geomap_static>, "GeoMap should satisfy GeoMapLike concept" );
    static_assert( geomap_static::is_order_static, "Order 1 should be static" );
    static_assert( !geomap_static::is_order_dynamic, "Order 1 should not be dynamic" );
    static_assert( geomap_static::nOrder == 1, "nOrder should be 1" );
    static_assert( geomap_static::nOrder_v == 1, "nOrder_v should be 1" );

    BOOST_CHECK( geomap_static::is_order_static );
}

BOOST_AUTO_TEST_CASE( test_geomap_dynamic_order_flags )
{
    using geomap_dynamic = GeoMap<2, Dynamic, 2, double, Simplex>;

    static_assert( GeoMapLike<geomap_dynamic>, "Dynamic GeoMap should satisfy GeoMapLike concept" );
    static_assert( !geomap_dynamic::is_order_static, "Dynamic should not be static" );
    static_assert( geomap_dynamic::is_order_dynamic, "Dynamic should be dynamic" );
    static_assert( geomap_dynamic::nOrder_v == Dynamic, "nOrder_v should be Dynamic" );
    static_assert( geomap_dynamic::nOrder == 1, "nOrder placeholder should be 1 for dynamic" );

    BOOST_CHECK( geomap_dynamic::is_order_dynamic );
}

BOOST_AUTO_TEST_CASE( test_geomap_type_distinctness )
{
    using geomap_p1 = GeoMap<2, 1, 2, double, Simplex>;
    using geomap_p2 = GeoMap<2, 2, 2, double, Simplex>;
    using geomap_dyn = GeoMap<2, Dynamic, 2, double, Simplex>;

    static_assert( !std::is_same_v<geomap_p1, geomap_p2>, "P1 and P2 should be different types" );
    static_assert( !std::is_same_v<geomap_p1, geomap_dyn>, "P1 and Dynamic should be different types" );
    static_assert( !std::is_same_v<geomap_p2, geomap_dyn>, "P2 and Dynamic should be different types" );

    BOOST_CHECK( true );
}

//=============================================================================
// SECTION 4: GeoMap Construction Tests
//=============================================================================

BOOST_AUTO_TEST_CASE( test_geomap_static_construction )
{
    using geomap_p1 = GeoMap<2, 1, 2, double, Simplex>;
    using geomap_p2 = GeoMap<2, 2, 2, double, Simplex>;

    geomap_p1 gm1;
    geomap_p2 gm2;

    BOOST_CHECK_EQUAL( gm1.order(), 1 );
    BOOST_CHECK_EQUAL( gm2.order(), 2 );
    BOOST_CHECK_EQUAL( gm1.dim(), 2 );
    BOOST_CHECK_EQUAL( gm2.dim(), 2 );
}

BOOST_AUTO_TEST_CASE( test_geomap_dynamic_construction )
{
    using geomap_dyn = GeoMap<2, Dynamic, 2, double, Simplex>;

    for ( uint16_type order : { 1, 2, 3 } )
    {
        geomap_dyn gm{ RuntimeOrder{ order } };

        BOOST_CHECK_EQUAL( gm.order(), order );
        BOOST_CHECK_EQUAL( gm.dim(), 2 );
        BOOST_CHECK( gm.is_order_dynamic );
    }
}

//=============================================================================
// SECTION 5: Hypercube GeoMap Tests
//=============================================================================

BOOST_AUTO_TEST_CASE( test_geomap_hypercube_static )
{
    using geomap_q1 = GeoMap<2, 1, 2, double, Hypercube>;
    using geomap_q2 = GeoMap<2, 2, 2, double, Hypercube>;

    static_assert( GeoMapLike<geomap_q1>, "Hypercube GeoMap should satisfy GeoMapLike" );
    static_assert( GeoMapLike<geomap_q2>, "Hypercube GeoMap should satisfy GeoMapLike" );

    geomap_q1 gm1;
    geomap_q2 gm2;

    BOOST_CHECK_EQUAL( gm1.order(), 1 );
    BOOST_CHECK_EQUAL( gm2.order(), 2 );
    BOOST_CHECK( geomap_q1::is_order_static );
    BOOST_CHECK( geomap_q2::is_order_static );
}

BOOST_AUTO_TEST_CASE( test_geomap_hypercube_dynamic )
{
    using geomap_dyn = GeoMap<2, Dynamic, 2, double, Hypercube>;

    static_assert( GeoMapLike<geomap_dyn>, "Dynamic Hypercube GeoMap should satisfy GeoMapLike" );

    for ( uint16_type order : { 1, 2, 3 } )
    {
        geomap_dyn gm{ RuntimeOrder{ order } };

        BOOST_CHECK_EQUAL( gm.order(), order );
        BOOST_CHECK( geomap_dyn::is_order_dynamic );
    }
}

//=============================================================================
// SECTION 6: Static vs Dynamic Value Comparison
//=============================================================================

BOOST_AUTO_TEST_CASE( test_static_vs_dynamic_p1 )
{
    using static_t = GeoMap<2, 1, 2, double, Simplex>;
    using dynamic_t = GeoMap<2, Dynamic, 2, double, Simplex>;

    static_t gm_static;
    dynamic_t gm_dynamic{ RuntimeOrder{ 1 } };

    BOOST_CHECK_EQUAL( gm_static.order(), gm_dynamic.order() );
    BOOST_CHECK_EQUAL( gm_static.dim(), gm_dynamic.dim() );
    BOOST_CHECK_EQUAL( gm_static.realDim(), gm_dynamic.realDim() );
}

BOOST_AUTO_TEST_CASE( test_static_vs_dynamic_p2 )
{
    using static_t = GeoMap<2, 2, 2, double, Simplex>;
    using dynamic_t = GeoMap<2, Dynamic, 2, double, Simplex>;

    static_t gm_static;
    dynamic_t gm_dynamic{ RuntimeOrder{ 2 } };

    BOOST_CHECK_EQUAL( gm_static.order(), gm_dynamic.order() );
    BOOST_CHECK_EQUAL( gm_static.dim(), gm_dynamic.dim() );
}

//=============================================================================
// SECTION 7: Entity Dynamic Order Tests
//=============================================================================

BOOST_AUTO_TEST_CASE( test_simplex_dynamic_entity )
{
    using entity_static = Simplex<2, 1, 2>;
    using entity_dynamic = Simplex<2, Dynamic, 2>;

    static_assert( entity_static::is_order_static, "Static simplex should be static" );
    static_assert( entity_dynamic::is_order_dynamic, "Dynamic simplex should be dynamic" );

    for ( uint16_type order : { 1, 2, 3 } )
    {
        entity_dynamic entity{ RuntimeOrder{ order } };
        BOOST_CHECK_EQUAL( entity.order(), order );
        BOOST_CHECK( entity.is_order_dynamic );

        uint16_type expected_points = ( order + 1 ) * ( order + 2 ) / 2;
        BOOST_CHECK_EQUAL( entity.nPointsTotal(), expected_points );
    }
}

BOOST_AUTO_TEST_CASE( test_hypercube_dynamic_entity )
{
    using entity_static = Hypercube<2, 1, 2>;
    using entity_dynamic = Hypercube<2, Dynamic, 2>;

    static_assert( entity_static::is_order_static, "Static hypercube should be static" );
    static_assert( entity_dynamic::is_order_dynamic, "Dynamic hypercube should be dynamic" );

    for ( uint16_type order : { 1, 2, 3 } )
    {
        entity_dynamic entity{ RuntimeOrder{ order } };
        BOOST_CHECK_EQUAL( entity.order(), order );
        BOOST_CHECK( entity.is_order_dynamic );

        uint16_type expected_points = ( order + 1 ) * ( order + 1 );
        BOOST_CHECK_EQUAL( entity.nPointsTotal(), expected_points );
    }
}

//=============================================================================
// SECTION 8: Concept Verification Tests
//=============================================================================

BOOST_AUTO_TEST_CASE( test_geomap_concepts )
{
    static_assert( GeoMapLike<GeoMap<1, 1, 1, double, Simplex>>, "1D P1 Simplex" );
    static_assert( GeoMapLike<GeoMap<2, 1, 2, double, Simplex>>, "2D P1 Simplex" );
    static_assert( GeoMapLike<GeoMap<3, 1, 3, double, Simplex>>, "3D P1 Simplex" );
    static_assert( GeoMapLike<GeoMap<2, Dynamic, 2, double, Simplex>>, "2D Dynamic Simplex" );
    static_assert( GeoMapLike<GeoMap<2, 1, 2, double, Hypercube>>, "2D Q1 Hypercube" );
    static_assert( GeoMapLike<GeoMap<2, Dynamic, 2, double, Hypercube>>, "2D Dynamic Hypercube" );

    BOOST_CHECK( true );
}

//=============================================================================
// SECTION 9: Static vs Dynamic Numerical Equivalence (Context + Hessian)
//=============================================================================

BOOST_AUTO_TEST_CASE( test_geomap_context_equivalence_simplex_2d_p1 )
{
    compareGeomapContextStaticDynamicSimplex2D<1>();
}

BOOST_AUTO_TEST_CASE( test_geomap_context_equivalence_simplex_2d_p2 )
{
    compareGeomapContextStaticDynamicSimplex2D<2>();
}

BOOST_AUTO_TEST_CASE( test_geomap_context_equivalence_simplex_2d_p3 )
{
    compareGeomapContextStaticDynamicSimplex2D<3>();
}

BOOST_AUTO_TEST_CASE( test_geomap_hessian_equivalence_simplex_2d_p1 )
{
    compareGeomapHessianStaticDynamicSimplex2D<1>();
}

BOOST_AUTO_TEST_CASE( test_geomap_hessian_equivalence_simplex_2d_p2 )
{
    compareGeomapHessianStaticDynamicSimplex2D<2>();
}

BOOST_AUTO_TEST_CASE( test_geomap_hessian_equivalence_simplex_2d_p3 )
{
    compareGeomapHessianStaticDynamicSimplex2D<3>();
}

BOOST_AUTO_TEST_SUITE_END()
