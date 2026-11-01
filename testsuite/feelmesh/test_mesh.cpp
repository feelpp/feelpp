/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2009 Université de Grenoble 1 (Joseph Fourier)
 * SPDX-FileCopyrightText: 2005-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
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
#include <optional>
#include <sstream>
#include <type_traits>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelmesh/geoentity.hpp>
#include <feel/feelmesh/refentity.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/detail/meshcontiguousnumberingmapping.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/partitionio.hpp>
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
 * @brief Check the point count of every locally owned mesh element.
 *
 * Empty local partitions are valid in distributed tests and must not
 * dereference `beginElement()`.
 *
 * @tparam MeshPtrType mesh shared-pointer type
 * @param mesh distributed mesh
 * @param expectedPointCount expected geometry-node count per element
 */
template <typename MeshPtrType>
void checkLocalElementPointCount(
    MeshPtrType const& mesh,
    uint16_type expectedPointCount )
{
    for ( auto const& elementReference : mesh->elements() )
    {
        BOOST_CHECK_EQUAL(
            elementReference.second.nPoints(),
            expectedPointCount );
    }
}

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
        return T.createMesh(
            _mesh = mesh,
            _name = fmt::format(
                "triangle_dynamic_g{}",
                runtime_order.value ),
            _force_rebuild = true,
            _straighten = false );
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

/**
 * @test Import and evaluate two-dimensional runtime G3--G5 simplex meshes.
 */
BOOST_AUTO_TEST_CASE( test_mesh_dynamic_order_g3_to_g5_elements )
{
    using namespace Feel;
    using mesh_type = Mesh<Simplex<2, Dynamic>>;

    for ( uint16_type order = 3; order <= 5; ++order )
    {
        test_mesh_detail::MeshTestFixture<mesh_type> testFixture(
            RuntimeOrder{ order }, 0.8 );
        BOOST_REQUIRE( testFixture.mesh() );
        BOOST_CHECK_EQUAL( testFixture.mesh()->order(), order );
        testFixture.testElements();

        for ( auto const& elementReference :
              testFixture.mesh()->elements() )
        {
            auto const& element = elementReference.second;
            BOOST_CHECK_EQUAL(
                element.nPoints(),
                Feel::detail::simplexTotal( 2, order ) );
            BOOST_CHECK_EQUAL(
                element.G().size2(),
                Feel::detail::simplexTotal( 2, order ) );
        }

        if ( order == 3 )
        {
            auto mergedMesh = merge(
                testFixture.mesh(),
                testFixture.mesh() );
            BOOST_REQUIRE( mergedMesh );
            BOOST_CHECK_EQUAL(
                mergedMesh->order(), order );
            BOOST_CHECK_EQUAL(
                mergedMesh->numGlobalElements(),
                2 * testFixture.mesh()
                        ->numGlobalElements() );
            Feel::test_mesh_detail::
                checkLocalElementPointCount(
                    mergedMesh,
                    Feel::detail::simplexTotal(
                        2, order ) );
        }
    }
}

/**
 * @test Propagate runtime order through tetrahedral faces, edges, and Context.
 */
BOOST_AUTO_TEST_CASE( test_mesh_dynamic_3d_g2_and_g5_topology )
{
    using namespace Feel;
    using mesh_type = Mesh<Simplex<3, Dynamic, 3>>;

    for ( uint16_type order : {uint16_type{ 2 },
                               uint16_type{ 5 }} )
    {
        auto mesh = mesh_type::New( RuntimeOrder{ order } );
        GeoTool::Node first( 0, 0, 0 );
        GeoTool::Node last( 1, 1, 1 );
        GeoTool::Cube cube( 0.9, "Omega", first, last );
        mesh = cube.createMesh(
            _mesh = mesh,
            _name = fmt::format(
                "cube_tetra_dynamic_g{}", order ),
            _force_rebuild = true,
            _straighten = false );

        BOOST_REQUIRE( mesh );
        BOOST_CHECK_EQUAL( mesh->order(), order );
        BOOST_CHECK_GT( mesh->numGlobalElements(), 0 );

        auto const cellPointCount =
            Feel::detail::simplexTotal( 3, order );
        auto const facePointCount =
            Feel::detail::simplexTotal( 2, order );
        auto const edgePointCount =
            Feel::detail::simplexTotal( 1, order );

        for ( auto const& elementReference : mesh->elements() )
        {
            auto const& element = elementReference.second;
            BOOST_CHECK_EQUAL(
                element.nPoints(), cellPointCount );
            BOOST_CHECK_EQUAL(
                element.G().size2(), cellPointCount );
        }
        for ( auto const& faceReference : mesh->faces() )
            BOOST_CHECK_EQUAL(
                faceReference.second.nPoints(),
                facePointCount );
        for ( auto const& edgeReference : mesh->edges() )
            BOOST_CHECK_EQUAL(
                edgeReference.second.nPoints(),
                edgePointCount );

        std::optional<typename mesh_type::geometry_coordinate_matrix_type>
            nodalCoordinates;
        if ( order == 5 )
        {
            nodalCoordinates = mesh->nodalGeometryCoordinates();
            mesh->setGeometryRepresentation(
                geometry::GeometryRepresentationKind::BernsteinBezier );
            BOOST_CHECK(
                mesh->geometryRepresentation() ==
                geometry::GeometryRepresentationKind::BernsteinBezier );
            BOOST_CHECK_SMALL(
                ( mesh->nodalGeometryCoordinates() -
                  *nodalCoordinates )
                    .norm(),
                2e-9 );
        }

        if ( mesh->numElements() > 0 )
        {
            auto const& element =
                mesh->beginElement()->second;
            auto gm = mesh->gm();
            typename mesh_type::reference_convex_type
                reference;
            auto precompute =
                gm->preCompute( gm, reference.points() );
            auto context =
                gm->template context<
                    vm::POINT | vm::JACOBIAN>(
                        element, precompute );
            BOOST_CHECK( context->isValid() );
        }

        if ( order == 5 )
        {
            auto volumeSubmesh = createSubmesh(
                _mesh = mesh, _range = elements( mesh ) );
            BOOST_REQUIRE( volumeSubmesh );
            BOOST_CHECK_EQUAL(
                volumeSubmesh->order(), order );
            BOOST_CHECK(
                volumeSubmesh->geometryRepresentation() ==
                geometry::GeometryRepresentationKind::BernsteinBezier );
            BOOST_REQUIRE_GT(
                volumeSubmesh->numGlobalElements(), 0 );
            Feel::test_mesh_detail::
                checkLocalElementPointCount(
                    volumeSubmesh, cellPointCount );

            auto boundaryMesh = createSubmesh(
                _mesh = mesh,
                _range = boundaryfaces( mesh ) );
            BOOST_REQUIRE( boundaryMesh );
            BOOST_CHECK_EQUAL(
                boundaryMesh->order(), order );
            BOOST_CHECK(
                boundaryMesh->geometryRepresentation() ==
                geometry::GeometryRepresentationKind::BernsteinBezier );
            BOOST_REQUIRE_GT(
                boundaryMesh->numGlobalElements(), 0 );
            Feel::test_mesh_detail::
                checkLocalElementPointCount(
                    boundaryMesh, facePointCount );

            using export_mapping_type =
                Feel::detail::MeshContiguousNumberingMapping<
                    mesh_type, float>;
            export_mapping_type exportMapping{
                mesh.get(), true };
            BOOST_CHECK_EQUAL(
                exportMapping.pointsPerElement(),
                cellPointCount );
            auto const exportedNodalCoordinates =
                mesh->nodalGeometryCoordinates();
            auto const localRank =
                mesh->worldComm().localRank();
            for ( auto const& [part, nameAndRange] :
                  exportMapping.partIdToRangeElement() )
            {
                auto const& connectivity =
                    exportMapping.pointIdsInElements( part );
                BOOST_CHECK_EQUAL(
                    connectivity.size(),
                    exportMapping.numberOfElement(
                        part, localRank )*cellPointCount );
                auto const& nodes = exportMapping.nodes( part );
                auto const startPoint =
                    exportMapping.startPointIds(
                        part, localRank );
                for ( auto const& [pointId, pointData] :
                      exportMapping.pointIdToContiguous( part ) )
                {
                    auto const exportedPoint =
                        pointData.first-startPoint;
                    auto const pointIndex =
                        mesh->pointLocalIndex( pointId );
                    for ( uint16_type dimension = 0;
                          dimension < mesh_type::nRealDim;
                          ++dimension )
                        BOOST_CHECK_SMALL(
                            nodes[3*exportedPoint+dimension]-
                                exportedNodalCoordinates(
                                    dimension, pointIndex ),
                            2e-6 );
                }
            }

            if ( Environment::worldComm()
                     .globalSize() == 1 )
            {
                std::stringstream archiveBuffer;
                {
                    boost::archive::text_oarchive archive{
                        archiveBuffer };
                    archive << *mesh;
                }
                auto archivedMesh =
                    mesh_type::New( RuntimeOrder{ 1 } );
                {
                    boost::archive::text_iarchive archive{
                        archiveBuffer };
                    archive >> *archivedMesh;
                }
                BOOST_CHECK_EQUAL(
                    archivedMesh->order(), order );
                BOOST_CHECK(
                    archivedMesh->geometryRepresentation() ==
                    geometry::GeometryRepresentationKind::BernsteinBezier );
                BOOST_CHECK_EQUAL(
                    archivedMesh->numGlobalElements(),
                    mesh->numGlobalElements() );
                BOOST_CHECK_EQUAL(
                    archivedMesh->beginElement()
                        ->second.nPoints(),
                    cellPointCount );
                BOOST_CHECK_SMALL(
                    ( archivedMesh->nodalGeometryCoordinates() -
                      *nodalCoordinates )
                        .norm(),
                    2e-9 );

                auto const metadataPath =
                    Environment::repository().directory() /
                    "runtime_geometry_g5_partition.json";
                PartitionIO<mesh_type> writer{
                    metadataPath.string() };
                writer.write( mesh );

                auto partitionMesh =
                    mesh_type::New( RuntimeOrder{ 1 } );
                PartitionIO<mesh_type> reader{
                    metadataPath.string() };
                reader.read( partitionMesh );
                BOOST_CHECK_EQUAL(
                    partitionMesh->order(), order );
                BOOST_CHECK(
                    partitionMesh->geometryRepresentation() ==
                    geometry::GeometryRepresentationKind::BernsteinBezier );
                BOOST_CHECK_EQUAL(
                    partitionMesh->numGlobalElements(),
                    mesh->numGlobalElements() );
                BOOST_CHECK_EQUAL(
                    partitionMesh->beginElement()
                        ->second.nPoints(),
                    cellPointCount );
                BOOST_CHECK_SMALL(
                    ( partitionMesh->nodalGeometryCoordinates() -
                      *nodalCoordinates )
                        .norm(),
                    2e-9 );
            }
        }
    }
}

/**
 * @test Import complete runtime-order quadrilateral and hexahedral G3/G5
 *       meshes without falling back to statically instantiated maps.
 */
BOOST_AUTO_TEST_CASE( test_mesh_dynamic_hypercube_g3_and_g5 )
{
    using namespace Feel;
    using quadrangle_mesh_type =
        Mesh<Hypercube<2, Dynamic>>;
    using hexahedron_mesh_type =
        Mesh<Hypercube<3, Dynamic, 3>>;

    for ( uint16_type order :
          {uint16_type{ 3 }, uint16_type{ 5 }} )
    {
        auto quadrangleMesh =
            quadrangle_mesh_type::New(
                RuntimeOrder{ order } );
        quadrangleMesh = createGMSHMesh(
            _mesh = quadrangleMesh,
            _desc = domain(
                _name = fmt::format(
                    "quadrangle_dynamic_g{}", order ),
                _shape = "hypercube",
                _convex = "Hypercube",
                _dim = 2,
                _order = order,
                _h = 0.5 ),
            _structured = 1,
            _force_rebuild = true,
            _straighten = false );
        BOOST_REQUIRE( quadrangleMesh );
        BOOST_CHECK_EQUAL(
            quadrangleMesh->order(), order );
        BOOST_REQUIRE_GT(
            quadrangleMesh->numGlobalElements(), 0 );
        Feel::test_mesh_detail::
            checkLocalElementPointCount(
                quadrangleMesh,
                Feel::detail::hypercubeTotal(
                    2, order ) );

        auto hexahedronMesh =
            hexahedron_mesh_type::New(
                RuntimeOrder{ order } );
        hexahedronMesh = createGMSHMesh(
            _mesh = hexahedronMesh,
            _desc = domain(
                _name = fmt::format(
                    "hexahedron_dynamic_g{}", order ),
                _shape = "hypercube",
                _convex = "Hypercube",
                _dim = 3,
                _order = order,
                _h = 0.5 ),
            _structured = 1,
            _force_rebuild = true,
            _straighten = false );
        BOOST_REQUIRE( hexahedronMesh );
        BOOST_CHECK_EQUAL(
            hexahedronMesh->order(), order );
        BOOST_REQUIRE_GT(
            hexahedronMesh->numGlobalElements(), 0 );
        Feel::test_mesh_detail::
            checkLocalElementPointCount(
                hexahedronMesh,
                Feel::detail::hypercubeTotal(
                    3, order ) );
    }
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

    // Independent Gmsh runs may partition cells differently across ranks.
    // Mathematical parity therefore compares the global mesh cardinality.
    BOOST_CHECK_EQUAL(
        static_test_fixture.mesh()->numGlobalElements(),
        dynamic_test_fixture.mesh()->numGlobalElements() );

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
    BOOST_CHECK_EQUAL(
        static_test_fixture.mesh()->numGlobalElements(),
        dynamic_test_fixture.mesh()->numGlobalElements() );

    static_test_fixture.testFilters();
    dynamic_test_fixture.testFilters();
    static_test_fixture.testElements();
    dynamic_test_fixture.testElements();

    BOOST_TEST_MESSAGE( "Static and dynamic order P2 meshes are consistent" );
}
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
