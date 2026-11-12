/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

#define BOOST_TEST_MODULE periodic mesh testsuite
#include <feel/feelcore/testsuite.hpp>

#include <sstream>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/loadgmshmesh.hpp>
#include <feel/feelfilters/straightenmesh.hpp>
#include <feel/feelfilters/straightenmesh_impl.hpp>
#include <feel/feelmesh/filters.hpp>

using namespace Feel;

namespace
{
using mesh_type = Mesh<Simplex<2, 1>>;
using mesh_ptrtype = std::shared_ptr<mesh_type>;
using mesh3d_type = Mesh<Simplex<3, 1>>;
using mesh3d_ptrtype = std::shared_ptr<mesh3d_type>;

bool gmshDefaultIsV4()
{
    std::string version = FEELPP_GMSH_FORMAT_VERSION;
    return !version.empty() && version.front() == '4';
}

std::string squareGeoDescription( bool periodic, double h = 0.2 )
{
    std::ostringstream ostr;
    ostr << "Mesh.MshFileVersion = " << FEELPP_GMSH_FORMAT_VERSION << ";\n"
         << "h=" << h << ";\n"
         << "Point(1) = {0,0,0,h};\n"
         << "Point(2) = {1,0,0,h};\n"
         << "Point(3) = {1,1,0,h};\n"
         << "Point(4) = {0,1,0,h};\n"
         << "Line(1) = {1,2};\n"
         << "Line(2) = {2,3};\n"
         << "Line(3) = {3,4};\n"
         << "Line(4) = {4,1};\n"
         << "Line Loop(5) = {1,2,3,4};\n"
         << "Plane Surface(6) = {5};\n"
         << "Physical Surface(\"Omega\") = {6};\n"
         << "Physical Line(\"Bottom\") = {1};\n"
         << "Physical Line(\"Right\") = {2};\n"
         << "Physical Line(\"Top\") = {3};\n"
         << "Physical Line(\"Left\") = {4};\n";
    if ( periodic )
    {
        ostr << "Periodic Curve {2} = {4} Translate {1,0,0};\n"
             << "Periodic Curve {3} = {1} Translate {0,1,0};\n";
    }
    return ostr.str();
}

std::string cubeGeoDescription( bool periodic, double h = 0.3 )
{
    std::ostringstream ostr;
    ostr << "Mesh.MshFileVersion = " << FEELPP_GMSH_FORMAT_VERSION << ";\n"
         << "SetFactory(\"OpenCASCADE\");\n"
         << "h=" << h << ";\n"
         << "Box(1) = {0,0,0,1,1,1};\n"
         << "Mesh.CharacteristicLengthMin = h;\n"
         << "Mesh.CharacteristicLengthMax = h;\n"
         << "eps=1e-6;\n"
         << "sxm[] = Surface In BoundingBox{-eps,-eps,-eps,eps,1+eps,1+eps};\n"
         << "sxp[] = Surface In BoundingBox{1-eps,-eps,-eps,1+eps,1+eps,1+eps};\n"
         << "sym[] = Surface In BoundingBox{-eps,-eps,-eps,1+eps,eps,1+eps};\n"
         << "syp[] = Surface In BoundingBox{-eps,1-eps,-eps,1+eps,1+eps,1+eps};\n"
         << "szm[] = Surface In BoundingBox{-eps,-eps,-eps,1+eps,1+eps,eps};\n"
         << "szp[] = Surface In BoundingBox{-eps,-eps,1-eps,1+eps,1+eps,1+eps};\n"
         << "Physical Volume(\"Omega\") = {1};\n"
         << "Physical Surface(\"XMin\") = {sxm[]};\n"
         << "Physical Surface(\"XMax\") = {sxp[]};\n"
         << "Physical Surface(\"YMin\") = {sym[]};\n"
         << "Physical Surface(\"YMax\") = {syp[]};\n"
         << "Physical Surface(\"ZMin\") = {szm[]};\n"
         << "Physical Surface(\"ZMax\") = {szp[]};\n";
    if ( periodic )
    {
        ostr << "Periodic Surface {sxp[]} = {sxm[]} Translate {1,0,0};\n"
             << "Periodic Surface {syp[]} = {sym[]} Translate {0,1,0};\n"
             << "Periodic Surface {szp[]} = {szm[]} Translate {0,0,1};\n";
    }
    return ostr.str();
}

mesh_ptrtype createMeshFromGmsh( bool periodic, GMSH_FORMAT format, std::string const& prefix )
{
    Gmsh gmsh;
    gmsh.setDimension( 2 );
    gmsh.setOrder( 1 );
    gmsh.setVersion( FEELPP_GMSH_FORMAT_VERSION, format );
    gmsh.setPrefix( prefix );

    std::string fname;
    bool generated_or_modified = false;
    boost::tie( fname, generated_or_modified ) = gmsh.generate( prefix, squareGeoDescription( periodic ), true );
    Feel::detail::ignore_unused_variable_warning( generated_or_modified );

    return loadGMSHMesh( _mesh = new mesh_type,
                         _filename = fname,
                         _update = MESH_CHECK | MESH_UPDATE_FACES | MESH_UPDATE_EDGES );
}

mesh3d_ptrtype createMesh3dFromGmsh( bool periodic, std::string const& prefix )
{
    Gmsh gmsh;
    gmsh.setDimension( 3 );
    gmsh.setOrder( 1 );
    gmsh.setVersion( FEELPP_GMSH_FORMAT_VERSION, GMSH_FORMAT_ASCII );
    gmsh.setPrefix( prefix );

    std::string fname;
    bool generated_or_modified = false;
    boost::tie( fname, generated_or_modified ) = gmsh.generate( prefix, cubeGeoDescription( periodic ), true );
    Feel::detail::ignore_unused_variable_warning( generated_or_modified );

    return loadGMSHMesh( _mesh = new mesh3d_type,
                         _filename = fname,
                         _update = MESH_CHECK | MESH_UPDATE_FACES | MESH_UPDATE_EDGES );
}

std::set<size_type> canonicalEdgesOnMarker( mesh_ptrtype const& mesh, std::string const& marker )
{
    std::set<size_type> canonicalEdges;
    auto range = markedfaces( mesh, marker );
    for ( auto it = range.begin(), en = range.end(); it != en; ++it )
        canonicalEdges.insert( mesh->canonicalEdgeId( unwrap_ref( *it ).id() ) );
    return canonicalEdges;
}

bool markersShareCanonicalEdgeSet( mesh_ptrtype const& mesh, std::string const& markerA, std::string const& markerB )
{
    auto const canonicalA = canonicalEdgesOnMarker( mesh, markerA );
    auto const canonicalB = canonicalEdgesOnMarker( mesh, markerB );
    return !canonicalA.empty() && canonicalA == canonicalB;
}

std::unordered_set<size_type> edgeIdsOnMarker( mesh_ptrtype const& mesh, std::string const& marker )
{
    std::unordered_set<size_type> edgeIds;
    auto range = markedfaces( mesh, marker );
    for ( auto it = range.begin(), en = range.end(); it != en; ++it )
        edgeIds.insert( unwrap_ref( *it ).id() );
    return edgeIds;
}

std::unordered_set<size_type> faceIdsOnMarker( mesh3d_ptrtype const& mesh, std::string const& marker )
{
    std::unordered_set<size_type> faceIds;
    auto range = markedfaces( mesh, marker );
    for ( auto it = range.begin(), en = range.end(); it != en; ++it )
        faceIds.insert( unwrap_ref( *it ).id() );
    return faceIds;
}
}

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( periodicmeshsuite )

BOOST_AUTO_TEST_CASE( periodic_mesh_canonical_ids_v4_ascii )
{
    if ( Environment::isParallel() || !gmshDefaultIsV4() )
        return;

    auto mesh = createMeshFromGmsh( true, GMSH_FORMAT_ASCII, "periodic-mesh-canonical-v4-ascii" );
    BOOST_REQUIRE( mesh );
    BOOST_REQUIRE( mesh->isPeriodic() );

    auto const& periodicEntities = mesh->periodicEntities();
    BOOST_REQUIRE( !periodicEntities.empty() );

    int checkedPairs = 0;
    for ( auto const& e : periodicEntities )
    {
        if ( e.dim != 1 )
            continue;
        for ( auto const& [slaveId, masterId] : e.correspondingVertices )
        {
            if ( !mesh->hasPoint( slaveId ) || !mesh->hasPoint( masterId ) )
                continue;
            BOOST_CHECK_EQUAL( mesh->canonicalPointId( slaveId ), mesh->canonicalPointId( masterId ) );
            if ( slaveId != masterId )
                ++checkedPairs;
        }
    }
    BOOST_CHECK_GT( checkedPairs, 0 );

    int canonicalizedPoints = 0;
    for ( auto it = mesh->beginPoint(); it != mesh->endPoint(); ++it )
    {
        auto const pid = it->second.id();
        if ( mesh->canonicalPointId( pid ) != pid )
            ++canonicalizedPoints;
    }
    BOOST_CHECK_GT( canonicalizedPoints, 0 );

    BOOST_CHECK( markersShareCanonicalEdgeSet( mesh, "Right", "Left" ) );
    BOOST_CHECK( markersShareCanonicalEdgeSet( mesh, "Top", "Bottom" ) );
}

BOOST_AUTO_TEST_CASE( periodic_mesh_master_relations_v4_ascii_2d )
{
    if ( Environment::isParallel() || !gmshDefaultIsV4() )
        return;

    auto mesh = createMeshFromGmsh( true, GMSH_FORMAT_ASCII, "periodic-mesh-relations-v4-ascii" );
    BOOST_REQUIRE( mesh );
    BOOST_REQUIRE( mesh->isPeriodic() );

    int checkedPointMasters = 0;
    for ( auto const& periodicEntity : mesh->periodicEntities() )
    {
        if ( periodicEntity.dim != 1 )
            continue;

        for ( auto const& [slaveId, masterId] : periodicEntity.correspondingVertices )
        {
            if ( !mesh->hasPoint( slaveId ) || slaveId == masterId )
                continue;
            BOOST_CHECK( mesh->hasPeriodicPointMaster( slaveId ) );
            BOOST_CHECK_EQUAL( mesh->canonicalPointId( mesh->periodicPointMasterId( slaveId ) ),
                               mesh->canonicalPointId( masterId ) );
            if ( mesh->hasPoint( masterId ) )
                BOOST_CHECK_EQUAL( mesh->canonicalPointId( slaveId ), mesh->canonicalPointId( masterId ) );
            ++checkedPointMasters;
        }
    }
    BOOST_CHECK_GT( checkedPointMasters, 0 );

    auto const rightEdges = edgeIdsOnMarker( mesh, "Right" );
    auto const leftEdges = edgeIdsOnMarker( mesh, "Left" );
    BOOST_REQUIRE( !rightEdges.empty() );
    BOOST_REQUIRE( !leftEdges.empty() );

    int checkedEdgeMasters = 0;
    for ( auto const edgeId : rightEdges )
    {
        if ( !mesh->hasPeriodicEdgeMaster( edgeId ) )
            continue;
        auto const masterEdgeId = mesh->periodicEdgeMasterId( edgeId );
        BOOST_CHECK( leftEdges.contains( masterEdgeId ) );
        BOOST_CHECK_EQUAL( mesh->canonicalEdgeId( edgeId ), mesh->canonicalEdgeId( masterEdgeId ) );
        auto const orientation = mesh->periodicEdgeMasterOrientation( edgeId );
        BOOST_CHECK( orientation == 1 || orientation == -1 );
        ++checkedEdgeMasters;
    }
    BOOST_CHECK_GT( checkedEdgeMasters, 0 );
}

BOOST_AUTO_TEST_CASE( periodic_mesh_face_relations_v4_ascii_3d )
{
    if ( Environment::isParallel() || !gmshDefaultIsV4() )
        return;

    auto mesh = createMesh3dFromGmsh( true, "periodic-mesh-face-relations-v4-ascii" );
    BOOST_REQUIRE( mesh );
    BOOST_REQUIRE( mesh->isPeriodic() );

    auto checkFacePair = [&]( std::string const& slaveMarker, std::string const& masterMarker )
    {
        auto const slaveFaces = faceIdsOnMarker( mesh, slaveMarker );
        auto const masterFaces = faceIdsOnMarker( mesh, masterMarker );
        BOOST_REQUIRE( !slaveFaces.empty() );
        BOOST_REQUIRE( !masterFaces.empty() );

        int checkedFaceMasters = 0;
        for ( auto const faceId : slaveFaces )
        {
            if ( !mesh->hasPeriodicFaceMaster( faceId ) )
                continue;
            auto const masterFaceId = mesh->periodicFaceMasterId( faceId );
            BOOST_CHECK( masterFaces.contains( masterFaceId ) );
            BOOST_CHECK_EQUAL( mesh->canonicalFaceId( faceId ), mesh->canonicalFaceId( masterFaceId ) );

            auto const permutation = mesh->periodicFaceMasterPermutation( faceId );
            BOOST_CHECK_EQUAL( permutation.size(), mesh3d_type::face_type::numVertices );
            std::set<uint16_type> permutationValues( permutation.begin(), permutation.end() );
            BOOST_CHECK_EQUAL( permutationValues.size(), mesh3d_type::face_type::numVertices );
            for ( auto const p : permutation )
                BOOST_CHECK_LT( p, mesh3d_type::face_type::numVertices );
            ++checkedFaceMasters;
        }
        BOOST_CHECK_GT( checkedFaceMasters, 0 );
    };

    checkFacePair( "XMax", "XMin" );
    checkFacePair( "YMax", "YMin" );
    checkFacePair( "ZMax", "ZMin" );
}

BOOST_AUTO_TEST_SUITE_END()
