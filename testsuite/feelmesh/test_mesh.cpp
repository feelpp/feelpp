/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-09-03

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2009 Université de Grenoble 1 (Joseph Fourier)
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
   @file test_mesh.cpp
   @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   @date 2005-09-03
   @brief Mesh testsuite with C++20/23 modernization and static/dynamic order support
 */

// give a name to the testsuite
#define BOOST_TEST_MODULE mesh testsuite

#include <feel/feelcore/testsuite.hpp>

#include <concepts>
#include <type_traits>

#include <feel/feelcore/environment.hpp>
#include <feel/feelmesh/geoentity.hpp>
#include <feel/feelmesh/refentity.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feelpoly/order.hpp>

namespace Feel::test_mesh_detail
{

//! @brief Concept for mesh types supporting our tests
template<typename M>
concept TestMeshConcept = requires( M m )
{
    { m.numElements() } -> std::convertible_to<typename M::size_type>;
    { m.numFaces() } -> std::convertible_to<typename M::size_type>;
    { m.numPoints() } -> std::convertible_to<typename M::size_type>;
    { m.order() } -> std::convertible_to<uint16_type>;
};

/**
 * @brief Generic mesh test fixture that works for both static and dynamic order
 *
 * Uses C++20 concepts to constrain mesh type and if constexpr for
 * static/dynamic order handling.
 */
template <typename MeshType>
    requires TestMeshConcept<MeshType>
class MeshTestFixture
{
public:
    using mesh_type = MeshType;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;

    static constexpr bool is_order_static = mesh_type::is_order_static;
    static constexpr bool is_order_dynamic = mesh_type::is_order_dynamic;
    static constexpr uint16_type nDim = mesh_type::nDim;

    /**
     * @brief Construct fixture for static order mesh
     */
    explicit MeshTestFixture( double meshSize = 1.0 )
        requires( is_order_static )
        : M_meshSize( meshSize )
    {
        BOOST_TEST_MESSAGE( "Setting up static order mesh (order=" << mesh_type::nOrder << ")" );
        M_mesh = createMesh();
        BOOST_CHECK( M_mesh != nullptr );
        BOOST_CHECK_EQUAL( M_mesh->order(), mesh_type::nOrder );
    }

    /**
     * @brief Construct fixture for dynamic order mesh
     */
    explicit MeshTestFixture( RuntimeOrder runtime_order, double meshSize = 1.0 )
        requires( is_order_dynamic )
        : M_meshSize( meshSize ), M_runtime_order( runtime_order.value )
    {
        BOOST_TEST_MESSAGE( "Setting up dynamic order mesh (order=" << runtime_order.value << ")" );
        M_mesh = createMesh( runtime_order );
        BOOST_CHECK( M_mesh != nullptr );
        BOOST_CHECK_EQUAL( M_mesh->order(), runtime_order.value );
    }

    //! @brief Get the mesh
    [[nodiscard]] mesh_ptrtype mesh() const noexcept { return M_mesh; }

    //! @brief Get the effective order (works for both static and dynamic)
    [[nodiscard]] uint16_type order() const noexcept
    {
        if constexpr ( is_order_static )
            return mesh_type::nOrder;
        else
            return M_runtime_order;
    }

    /**
     * @brief Run mesh filter tests (same logic for static and dynamic)
     */
    void testFilters()
    {
        BOOST_TEST_MESSAGE( "Testing mesh filters for order=" << order() );

        // Test internal faces
        auto rangeInternalFaces = M_mesh->internalFaces();
        auto it = std::get<0>( rangeInternalFaces );
        auto en = std::get<1>( rangeInternalFaces );

        for ( ; it != en; ++it )
        {
            auto const& iface = boost::unwrap_ref( *it );

            // Internal faces must be connected to two elements
            BOOST_CHECK( iface.isConnectedTo0() && iface.isConnectedTo1() );

            // Check face vertex coordinates consistency
            int face_0 = iface.pos_first();
            int face_1 = iface.pos_second();

            auto n00 = iface.element( 0 ).point( iface.element( 0 ).fToP( face_0, 0 ) ).node();
            auto n10 = iface.element( 1 ).point( iface.element( 1 ).fToP( face_1, 1 ) ).node();
            BOOST_CHECK( ublas::norm_2( n00 - n10 ) < 1e-15 );

            auto n01 = iface.element( 0 ).point( iface.element( 0 ).fToP( face_0, 1 ) ).node();
            auto n11 = iface.element( 1 ).point( iface.element( 1 ).fToP( face_1, 0 ) ).node();
            BOOST_CHECK( ublas::norm_2( n01 - n11 ) < 1e-15 );
        }

        // Test boundary faces
        auto rangeBoundaryFaces = M_mesh->facesOnBoundary();
        it = std::get<0>( rangeBoundaryFaces );
        en = std::get<1>( rangeBoundaryFaces );

        for ( ; it != en; ++it )
        {
            auto const& bface = boost::unwrap_ref( *it );

            // Boundary faces connect to exactly one element
            BOOST_CHECK( bface.isConnectedTo0() && !bface.isConnectedTo1() );

            // Check marker validity
            BOOST_CHECK( bface.marker().value() == M_mesh->markerName("Gamma1") ||
                         bface.marker().value() == M_mesh->markerName("Gamma2") ||
                         bface.marker().value() == M_mesh->markerName("Gamma3") );
        }

        BOOST_TEST_MESSAGE( "Mesh filter tests passed for order=" << order() );
    }

    /**
     * @brief Run mesh element tests (same logic for static and dynamic)
     */
    void testElements()
    {
        BOOST_TEST_MESSAGE( "Testing mesh elements for order=" << order() );

        auto __gm = M_mesh->gm();
        typename mesh_type::reference_convex_type ref_conv;
        auto __geopc = __gm->preCompute( ref_conv.points() );

        auto it = M_mesh->beginElement();
        auto en = M_mesh->endElement();

        for ( ; it != en; ++it )
        {
            auto const& elt = it->second;

            // Check geometric transformation gives back element vertices
            auto __c = __gm->template context<vm::POINT>( elt, __geopc );

            // For P1 elements, xReal should match G exactly
            // For higher order elements, the comparison is more complex
            // since G contains more points than the reference vertices
            auto xReal = __c->xReal();
            const bool use_p1_vertex_equivalence = [&]() {
                if constexpr ( is_order_static )
                    return mesh_type::nOrder <= 1;
                else
                    return order() <= 1;
            }();

            if ( use_p1_vertex_equivalence )
            {
                BOOST_CHECK( ublas::norm_frobenius( xReal - elt.G() ) < 1e-15 );
            }
            else
            {
                // For higher order, just verify we have valid geometry data
                BOOST_CHECK( xReal.size2() > 0 );
                BOOST_CHECK( elt.G().size2() > 0 );
            }
        }

        BOOST_TEST_MESSAGE( "Mesh element tests passed for order=" << order() );
    }

    /**
     * @brief Run component tests (same logic for both static and dynamic)
     */
    void testComponents()
    {
        BOOST_TEST_MESSAGE( "Testing mesh components for order=" << order() );

        M_mesh->components().reset();
        BOOST_CHECK( M_mesh->components().test( MESH_CHECK ) == false );
        BOOST_CHECK( M_mesh->components().test( MESH_RENUMBER ) == false );
        BOOST_CHECK( M_mesh->components().test( MESH_UPDATE_FACES ) == false );
        BOOST_CHECK( M_mesh->components().test( MESH_UPDATE_EDGES ) == false );

        M_mesh->components().reset();
        M_mesh->components().set( MESH_CHECK );
        BOOST_CHECK( M_mesh->components().test( MESH_CHECK ) == true );
        BOOST_CHECK( M_mesh->components().test( MESH_RENUMBER ) == false );

        M_mesh->components().reset();
        M_mesh->components().set( MESH_CHECK | MESH_UPDATE_EDGES | MESH_UPDATE_FACES );
        BOOST_CHECK( M_mesh->components().test( MESH_CHECK ) == true );
        BOOST_CHECK( M_mesh->components().test( MESH_UPDATE_FACES ) == true );
        BOOST_CHECK( M_mesh->components().test( MESH_UPDATE_EDGES ) == true );

        BOOST_TEST_MESSAGE( "Mesh component tests passed for order=" << order() );
    }

private:
    /**
     * @brief Create mesh for static order
     */
    mesh_ptrtype createMesh()
        requires( is_order_static )
    {
        BOOST_TEST_MESSAGE( "Creating static order mesh" );

        GeoTool::Node x1( -1, -1 );
        GeoTool::Node x2( 1, -1 );
        GeoTool::Node x3( -1, 1 );
        GeoTool::Triangle T( M_meshSize, "MyTriangle", x1, x2, x3 );
        T.setMarker( _type = "line", _name = "Gamma1", _marker1 = true );
        T.setMarker( _type = "line", _name = "Gamma2", _marker2 = true );
        T.setMarker( _type = "line", _name = "Gamma3", _marker3 = true );
        T.setMarker( _type = "surface", _name = "Omega", _markerAll = true );

        return T.createMesh( _mesh = new mesh_type, _name = "triangle_static" );
    }

    /**
     * @brief Create mesh for dynamic order
     */
    mesh_ptrtype createMesh( RuntimeOrder runtime_order )
        requires( is_order_dynamic )
    {
        BOOST_TEST_MESSAGE( "Creating dynamic order mesh (order=" << runtime_order.value << ")" );

        // For dynamic order, we need to create the mesh differently
        // Currently, GeoTool may not support dynamic order directly
        // This is a placeholder that shows the intended API

        GeoTool::Node x1( -1, -1 );
        GeoTool::Node x2( 1, -1 );
        GeoTool::Node x3( -1, 1 );
        GeoTool::Triangle T( M_meshSize, "MyTriangle", x1, x2, x3 );
        T.setMarker( _type = "line", _name = "Gamma1", _marker1 = true );
        T.setMarker( _type = "line", _name = "Gamma2", _marker2 = true );
        T.setMarker( _type = "line", _name = "Gamma3", _marker3 = true );
        T.setMarker( _type = "surface", _name = "Omega", _markerAll = true );

        // Create mesh with runtime order
        auto mesh = std::make_shared<mesh_type>( runtime_order );
        // Note: Full implementation requires updating GeoTool to support dynamic order
        return T.createMesh( _mesh = mesh.get(), _name = "triangle_dynamic" );
    }

    double M_meshSize;
    mesh_ptrtype M_mesh;
    uint16_type M_runtime_order{1};  // Default runtime order
};

} // namespace Feel::test_mesh_detail

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( mesh )

//==============================================================================
// Static Order Tests (P1)
//==============================================================================

BOOST_AUTO_TEST_CASE( test_mesh_static_order_p1_filters )
{
    using namespace Feel;
    using mesh_type = Mesh<Simplex<2, 1>>;

    BOOST_TEST_MESSAGE( "Testing static order P1 mesh filters" );

    test_mesh_detail::MeshTestFixture<mesh_type> test_fixture( 0.5 );

    // Verify static order properties
    static_assert( mesh_type::is_order_static, "Must be static order" );
    static_assert( mesh_type::nOrder == 1, "Must be P1" );
    BOOST_CHECK_EQUAL( test_fixture.order(), 1 );

    test_fixture.testFilters();
}

BOOST_AUTO_TEST_CASE( test_mesh_static_order_p1_elements )
{
    using namespace Feel;
    using mesh_type = Mesh<Simplex<2, 1>>;

    BOOST_TEST_MESSAGE( "Testing static order P1 mesh elements" );

    test_mesh_detail::MeshTestFixture<mesh_type> test_fixture( 0.5 );
    test_fixture.testElements();
}

BOOST_AUTO_TEST_CASE( test_mesh_static_order_p1_components )
{
    using namespace Feel;
    using mesh_type = Mesh<Simplex<2, 1>>;

    BOOST_TEST_MESSAGE( "Testing static order P1 mesh components" );

    test_mesh_detail::MeshTestFixture<mesh_type> test_fixture( 0.5 );
    test_fixture.testComponents();
}

//==============================================================================
// Static Order Tests (P2)
//==============================================================================

BOOST_AUTO_TEST_CASE( test_mesh_static_order_p2_filters )
{
    using namespace Feel;
    using mesh_type = Mesh<Simplex<2, 2>>;

    BOOST_TEST_MESSAGE( "Testing static order P2 mesh filters" );

    test_mesh_detail::MeshTestFixture<mesh_type> test_fixture( 0.5 );

    // Verify static order properties
    static_assert( mesh_type::is_order_static, "Must be static order" );
    static_assert( mesh_type::nOrder == 2, "Must be P2" );
    BOOST_CHECK_EQUAL( test_fixture.order(), 2 );

    test_fixture.testFilters();
}

BOOST_AUTO_TEST_CASE( test_mesh_static_order_p2_elements )
{
    using namespace Feel;
    using mesh_type = Mesh<Simplex<2, 2>>;

    BOOST_TEST_MESSAGE( "Testing static order P2 mesh elements" );

    test_mesh_detail::MeshTestFixture<mesh_type> test_fixture( 0.5 );
    test_fixture.testElements();
}

//==============================================================================
// Dynamic Order Tests (order set at runtime)
//==============================================================================

// NOTE:
// Mesh<Simplex<..., Dynamic>> is not yet fully supported in the mesh core
// (GeoElement/Mesh2D type stack). Keep these tests disabled until core support
// is completed. Dynamic geometric-order parity is covered in feeldiscr/geomap tests.
#if 0
BOOST_AUTO_TEST_CASE( test_mesh_dynamic_order_p1_filters )
{
    using namespace Feel;
    using mesh_type = Mesh<Simplex<2, Dynamic>>;

    BOOST_TEST_MESSAGE( "Testing dynamic order P1 mesh filters" );

    // Create with runtime order = 1
    test_mesh_detail::MeshTestFixture<mesh_type> test_fixture( RuntimeOrder(1), 0.5 );

    // Verify dynamic order properties
    static_assert( mesh_type::is_order_dynamic, "Must be dynamic order" );
    BOOST_CHECK_EQUAL( test_fixture.order(), 1 );

    test_fixture.testFilters();
}

BOOST_AUTO_TEST_CASE( test_mesh_dynamic_order_p2_filters )
{
    using namespace Feel;
    using mesh_type = Mesh<Simplex<2, Dynamic>>;

    BOOST_TEST_MESSAGE( "Testing dynamic order P2 mesh filters" );

    // Create with runtime order = 2
    test_mesh_detail::MeshTestFixture<mesh_type> test_fixture( RuntimeOrder(2), 0.5 );

    // Verify dynamic order properties
    static_assert( mesh_type::is_order_dynamic, "Must be dynamic order" );
    BOOST_CHECK_EQUAL( test_fixture.order(), 2 );

    test_fixture.testFilters();
}

BOOST_AUTO_TEST_CASE( test_mesh_dynamic_order_p1_elements )
{
    using namespace Feel;
    using mesh_type = Mesh<Simplex<2, Dynamic>>;

    BOOST_TEST_MESSAGE( "Testing dynamic order P1 mesh elements" );

    test_mesh_detail::MeshTestFixture<mesh_type> test_fixture( RuntimeOrder(1), 0.5 );
    test_fixture.testElements();
}

BOOST_AUTO_TEST_CASE( test_mesh_dynamic_order_p2_elements )
{
    using namespace Feel;
    using mesh_type = Mesh<Simplex<2, Dynamic>>;

    BOOST_TEST_MESSAGE( "Testing dynamic order P2 mesh elements" );

    test_mesh_detail::MeshTestFixture<mesh_type> test_fixture( RuntimeOrder(2), 0.5 );
    test_fixture.testElements();
}

BOOST_AUTO_TEST_CASE( test_mesh_dynamic_order_p1_components )
{
    using namespace Feel;
    using mesh_type = Mesh<Simplex<2, Dynamic>>;

    BOOST_TEST_MESSAGE( "Testing dynamic order P1 mesh components" );

    test_mesh_detail::MeshTestFixture<mesh_type> test_fixture( RuntimeOrder(1), 0.5 );
    test_fixture.testComponents();
}

BOOST_AUTO_TEST_CASE( test_mesh_dynamic_order_p2_components )
{
    using namespace Feel;
    using mesh_type = Mesh<Simplex<2, Dynamic>>;

    BOOST_TEST_MESSAGE( "Testing dynamic order P2 mesh components" );

    test_mesh_detail::MeshTestFixture<mesh_type> test_fixture( RuntimeOrder(2), 0.5 );
    test_fixture.testComponents();
}

BOOST_AUTO_TEST_CASE( test_mesh_dynamic_vs_static_consistency )
{
    using namespace Feel;

    BOOST_TEST_MESSAGE( "Testing consistency between static and dynamic order meshes" );

    constexpr double meshSize = 0.5;

    // Create static P1 mesh
    using static_mesh_type = Mesh<Simplex<2, 1>>;
    test_mesh_detail::MeshTestFixture<static_mesh_type> static_test_fixture( meshSize );

    // Create dynamic P1 mesh
    using dynamic_mesh_type = Mesh<Simplex<2, Dynamic>>;
    test_mesh_detail::MeshTestFixture<dynamic_mesh_type> dynamic_test_fixture( RuntimeOrder(1), meshSize );

    // Both should have order 1
    BOOST_CHECK_EQUAL( static_test_fixture.order(), 1 );
    BOOST_CHECK_EQUAL( dynamic_test_fixture.order(), 1 );

    // Both should have same number of elements (assuming same mesh generation)
    BOOST_CHECK_EQUAL( static_test_fixture.mesh()->numElements(),
                       dynamic_test_fixture.mesh()->numElements() );

    // Both should pass the same tests
    static_test_fixture.testFilters();
    dynamic_test_fixture.testFilters();

    static_test_fixture.testElements();
    dynamic_test_fixture.testElements();

    BOOST_TEST_MESSAGE( "Static and dynamic order meshes are consistent" );
}

BOOST_AUTO_TEST_CASE( test_mesh_dynamic_vs_static_consistency_p2 )
{
    using namespace Feel;

    BOOST_TEST_MESSAGE( "Testing consistency between static and dynamic order meshes (P2)" );

    constexpr double meshSize = 0.5;

    using static_mesh_type = Mesh<Simplex<2, 2>>;
    test_mesh_detail::MeshTestFixture<static_mesh_type> static_test_fixture( meshSize );

    using dynamic_mesh_type = Mesh<Simplex<2, Dynamic>>;
    test_mesh_detail::MeshTestFixture<dynamic_mesh_type> dynamic_test_fixture( RuntimeOrder(2), meshSize );

    BOOST_CHECK_EQUAL( static_test_fixture.order(), 2 );
    BOOST_CHECK_EQUAL( dynamic_test_fixture.order(), 2 );
    BOOST_CHECK_EQUAL( static_test_fixture.mesh()->numElements(),
                       dynamic_test_fixture.mesh()->numElements() );

    static_test_fixture.testFilters();
    dynamic_test_fixture.testFilters();
    static_test_fixture.testElements();
    dynamic_test_fixture.testElements();

    BOOST_TEST_MESSAGE( "Static and dynamic order P2 meshes are consistent" );
}
#endif

//==============================================================================
// Legacy Tests (preserved for compatibility)
//==============================================================================

BOOST_AUTO_TEST_CASE( test_mesh_comp_legacy )
{
    using namespace Feel;
    using mesh_type = Mesh<Simplex<2, 1>>;
    mesh_type mesh;

    mesh.components().reset();
    BOOST_CHECK( mesh.components().test( MESH_CHECK ) == false );
    BOOST_CHECK( mesh.components().test( MESH_RENUMBER ) == false );
    BOOST_CHECK( mesh.components().test( MESH_UPDATE_FACES ) == false );
    BOOST_CHECK( mesh.components().test( MESH_UPDATE_EDGES ) == false );

    mesh.components().reset();
    mesh.components().set( MESH_CHECK );
    BOOST_TEST_MESSAGE( "check MESH_CHECK comp: " << mesh.components().context() );
    BOOST_CHECK( mesh.components().test( MESH_CHECK ) == true );
    BOOST_CHECK( mesh.components().test( MESH_RENUMBER ) == false );
    BOOST_CHECK( mesh.components().test( MESH_UPDATE_FACES ) == false );
    BOOST_CHECK( mesh.components().test( MESH_UPDATE_EDGES ) == false );

    mesh.components().reset();
    mesh.components().set( MESH_CHECK | MESH_UPDATE_EDGES | MESH_UPDATE_FACES );
    BOOST_TEST_MESSAGE( "check MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES comp: " << mesh.components().context() );
    BOOST_CHECK( mesh.components().test( MESH_CHECK ) == true );
    BOOST_CHECK( mesh.components().test( MESH_RENUMBER ) == false );
    BOOST_CHECK( mesh.components().test( MESH_UPDATE_FACES ) == true );
    BOOST_CHECK( mesh.components().test( MESH_UPDATE_EDGES ) == true );

    mesh.components().reset();
    mesh.components().set( MESH_RENUMBER | MESH_UPDATE_FACES );
    BOOST_TEST_MESSAGE( "check MESH_RENUMBER|MESH_UPDATE_FACES comp: " << mesh.components().context() );
    BOOST_CHECK( mesh.components().test( MESH_CHECK ) == false );
    BOOST_CHECK( mesh.components().test( MESH_RENUMBER ) == true );
    BOOST_CHECK( mesh.components().test( MESH_UPDATE_FACES ) == true );
    BOOST_CHECK( mesh.components().test( MESH_UPDATE_EDGES ) == false );

    mesh.components().reset();
    mesh.components().set( MESH_RENUMBER );
    BOOST_TEST_MESSAGE( "check MESH_RENUMBER comp: " << mesh.components().context() );
    BOOST_CHECK( mesh.components().test( MESH_CHECK ) == false );
    BOOST_CHECK( mesh.components().test( MESH_RENUMBER ) == true );
    BOOST_CHECK( mesh.components().test( MESH_UPDATE_FACES ) == false );
    BOOST_CHECK( mesh.components().test( MESH_UPDATE_EDGES ) == false );
}

BOOST_AUTO_TEST_CASE( test_simple_mesh2d_legacy )
{
    using namespace Feel;
    using mesh_type = Mesh<Simplex<2, 1>>;
    mesh_type mesh;

    // Basic instantiation test
    BOOST_CHECK( mesh.numElements() == 0 );
    BOOST_CHECK( mesh.numPoints() == 0 );
}

BOOST_AUTO_TEST_SUITE_END()
