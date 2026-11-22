/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

#define BOOST_TEST_MODULE periodic mesh functionspace testsuite
#include <feel/feelcore/testsuite.hpp>

#include <sstream>
#include <set>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include <feel/feeldiscr/bdmh.hpp>
#include <feel/feeldiscr/dh.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
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

mesh_ptrtype createMeshFromGmsh( bool periodic, std::string const& prefix )
{
    Gmsh gmsh;
    gmsh.setDimension( 2 );
    gmsh.setOrder( 1 );
    gmsh.setVersion( FEELPP_GMSH_FORMAT_VERSION, GMSH_FORMAT_ASCII );
    gmsh.setPrefix( prefix );

    std::string fname;
    bool generated_or_modified = false;
    boost::tie( fname, generated_or_modified ) = gmsh.generate( prefix, squareGeoDescription( periodic ), true );
    Feel::detail::ignore_unused_variable_warning( generated_or_modified );

    return loadGMSHMesh( _mesh = new mesh_type,
                         _filename = fname,
                         _update = MESH_CHECK | MESH_UPDATE_FACES | MESH_UPDATE_EDGES );
}

size_type countCanonicalPoints( mesh_ptrtype const& mesh )
{
    std::set<size_type> canonicalPoints;
    for ( auto it = mesh->beginPoint(), en = mesh->endPoint(); it != en; ++it )
        canonicalPoints.insert( mesh->canonicalPointId( it->second.id() ) );
    return canonicalPoints.size();
}

size_type countCanonicalEdges( mesh_ptrtype const& mesh )
{
    std::set<size_type> canonicalEdges;
    for ( auto it = mesh->beginEdge(), en = mesh->endEdge(); it != en; ++it )
        canonicalEdges.insert( mesh->canonicalEdgeId( it->second.id() ) );
    return canonicalEdges.size();
}

std::vector<size_type> boundaryEdgeIds( mesh_ptrtype const& mesh, std::string const& marker )
{
    std::vector<size_type> edgeIds;
    auto range = markedfaces( mesh, marker );
    for ( auto it = range.begin(), en = range.end(); it != en; ++it )
        edgeIds.push_back( unwrap_ref( *it ).id() );
    return edgeIds;
}

template <typename SpacePtrType>
std::unordered_map<size_type, std::set<size_type>>
collectEdgeDofSetsByEdge( SpacePtrType const& space )
{
    std::unordered_map<size_type, std::set<size_type>> edgeDofs;
    auto const& mesh = space->mesh();
    auto const dof = space->dof();

    auto rangeElements = mesh->elementsWithProcessId( mesh->worldComm().localRank() );
    for ( auto itElement = std::get<0>( rangeElements ), enElement = std::get<1>( rangeElements );
          itElement != enElement; ++itElement )
    {
        auto const& elt = unwrap_ref( *itElement );
        for ( uint16_type i = 0; i < mesh_type::element_type::numEdges; ++i )
        {
            auto const edgeId = elt.edge( i ).id();
            auto edgeLocalDofs = dof->edgeLocalDof( elt.id(), i );
            auto& edgeSet = edgeDofs[edgeId];
            for ( auto const& edgeDof : edgeLocalDofs )
                edgeSet.insert( edgeDof.index() );
        }
    }

    return edgeDofs;
}

size_type edgeDofCardinality( std::unordered_map<size_type, std::set<size_type>> const& edgeDofs )
{
    BOOST_REQUIRE( !edgeDofs.empty() );
    size_type card = 0;
    for ( auto const& [edgeId, edgeSet] : edgeDofs )
    {
        Feel::detail::ignore_unused_variable_warning( edgeId );
        BOOST_REQUIRE( !edgeSet.empty() );
        if ( card == 0 )
            card = edgeSet.size();
        else
            BOOST_CHECK_EQUAL( edgeSet.size(), card );
    }
    return card;
}

template <typename SpacePtrType>
void checkPeriodicBoundaryEdgesShareDofs( SpacePtrType const& space,
                                          mesh_ptrtype const& mesh,
                                          std::string const& slaveMarker,
                                          std::string const& masterMarker )
{
    auto edgeDofs = collectEdgeDofSetsByEdge( space );
    auto slaveEdges = boundaryEdgeIds( mesh, slaveMarker );
    auto masterEdges = boundaryEdgeIds( mesh, masterMarker );

    BOOST_REQUIRE( !slaveEdges.empty() );
    BOOST_REQUIRE_EQUAL( slaveEdges.size(), masterEdges.size() );

    std::unordered_map<size_type, size_type> masterByCanonicalEdgeId;
    for ( auto const edgeId : masterEdges )
    {
        auto const canonicalId = mesh->canonicalEdgeId( edgeId );
        auto [itMaster, inserted] = masterByCanonicalEdgeId.emplace( canonicalId, edgeId );
        BOOST_CHECK( inserted || itMaster->second == edgeId );
    }

    for ( auto const slaveEdgeId : slaveEdges )
    {
        auto const canonicalId = mesh->canonicalEdgeId( slaveEdgeId );
        auto itMaster = masterByCanonicalEdgeId.find( canonicalId );
        BOOST_REQUIRE( itMaster != masterByCanonicalEdgeId.end() );

        auto itSlaveDofs = edgeDofs.find( slaveEdgeId );
        auto itMasterDofs = edgeDofs.find( itMaster->second );
        BOOST_REQUIRE( itSlaveDofs != edgeDofs.end() );
        BOOST_REQUIRE( itMasterDofs != edgeDofs.end() );
        BOOST_CHECK( itSlaveDofs->second == itMasterDofs->second );
    }
}

template <typename SpacePtrType>
void checkPeriodicVertexPairsShareDofIds( SpacePtrType const& space, mesh_ptrtype const& mesh,
                                          uint16_type nComponents = 1 )
{
    auto const relation = space->dof()->pointIdToDofRelation( "", false, true );
    auto const& pointToDof = relation.second;

    int checkedPairs = 0;
    for ( auto const& periodicEntity : mesh->periodicEntities() )
    {
        if ( periodicEntity.dim != 1 )
            continue;

        for ( auto const& [slaveId, masterId] : periodicEntity.correspondingVertices )
        {
            if ( !mesh->hasPoint( slaveId ) || !mesh->hasPoint( masterId ) )
                continue;

            for ( uint16_type c = 0; c < nComponents; ++c )
            {
                auto itSlave = pointToDof.find( nComponents * slaveId + c );
                auto itMaster = pointToDof.find( nComponents * masterId + c );
                BOOST_REQUIRE( itSlave != pointToDof.end() );
                BOOST_REQUIRE( itMaster != pointToDof.end() );
                BOOST_CHECK_EQUAL( itSlave->second, itMaster->second );
            }
            if ( slaveId != masterId )
                ++checkedPairs;
        }
    }
    BOOST_CHECK_GT( checkedPairs, 0 );
}
}

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( periodicspacesuite )

BOOST_AUTO_TEST_CASE( periodic_mesh_drives_dof_numbering )
{
    if ( Environment::isParallel() || !gmshDefaultIsV4() )
        return;

    auto meshPeriodic = createMeshFromGmsh( true, "periodic-space-v4" );
    auto meshNonPeriodic = createMeshFromGmsh( false, "nonperiodic-space-v4" );

    BOOST_REQUIRE( meshPeriodic );
    BOOST_REQUIRE( meshNonPeriodic );
    BOOST_REQUIRE( meshPeriodic->isPeriodic() );
    auto const periodicCanonicalPoints = countCanonicalPoints( meshPeriodic );
    auto const nonPeriodicCanonicalPoints = countCanonicalPoints( meshNonPeriodic );
    auto const periodicCanonicalEdges = countCanonicalEdges( meshPeriodic );
    auto const nonPeriodicCanonicalEdges = countCanonicalEdges( meshNonPeriodic );

    auto XhPeriodic = Pch<1>( meshPeriodic );
    auto XhNonPeriodic = Pch<1>( meshNonPeriodic );
    BOOST_CHECK( XhPeriodic->meshHasPeriodicity() );
    BOOST_CHECK( !XhNonPeriodic->meshHasPeriodicity() );
    BOOST_CHECK( !std::remove_reference_t<decltype( *XhPeriodic )>::has_type_level_periodicity );
    BOOST_CHECK_EQUAL( XhPeriodic->nDof(), periodicCanonicalPoints );
    BOOST_CHECK_EQUAL( XhNonPeriodic->nDof(), nonPeriodicCanonicalPoints );
    BOOST_CHECK_LT( XhPeriodic->nDof(), XhNonPeriodic->nDof() );
    checkPeriodicVertexPairsShareDofIds( XhPeriodic, meshPeriodic );

    auto XhvPeriodic = Pchv<1>( meshPeriodic );
    auto XhvNonPeriodic = Pchv<1>( meshNonPeriodic );
    constexpr uint16_type nVectorComponents = Mesh<Simplex<2, 1>>::nRealDim;
    BOOST_CHECK_EQUAL( XhvPeriodic->nDof(), nVectorComponents * periodicCanonicalPoints );
    BOOST_CHECK_EQUAL( XhvNonPeriodic->nDof(), nVectorComponents * nonPeriodicCanonicalPoints );
    BOOST_CHECK_LT( XhvPeriodic->nDof(), XhvNonPeriodic->nDof() );
    checkPeriodicVertexPairsShareDofIds( XhvPeriodic, meshPeriodic, nVectorComponents );

    auto RThPeriodic = Dh<0>( meshPeriodic );
    auto RThNonPeriodic = Dh<0>( meshNonPeriodic );
    auto const rtEdgeDofsPeriodic = collectEdgeDofSetsByEdge( RThPeriodic );
    auto const rtEdgeDofsNonPeriodic = collectEdgeDofSetsByEdge( RThNonPeriodic );
    auto const rtDofPerEdgePeriodic = edgeDofCardinality( rtEdgeDofsPeriodic );
    auto const rtDofPerEdgeNonPeriodic = edgeDofCardinality( rtEdgeDofsNonPeriodic );
    BOOST_CHECK_EQUAL( rtDofPerEdgePeriodic, rtDofPerEdgeNonPeriodic );
    BOOST_REQUIRE_GE( RThNonPeriodic->nDof(), nonPeriodicCanonicalEdges * rtDofPerEdgeNonPeriodic );
    auto const rtEdgeInvariantOffset = RThNonPeriodic->nDof() - nonPeriodicCanonicalEdges * rtDofPerEdgeNonPeriodic;
    BOOST_CHECK_EQUAL( RThPeriodic->nDof(), rtEdgeInvariantOffset + periodicCanonicalEdges * rtDofPerEdgePeriodic );
    BOOST_CHECK_LT( RThPeriodic->nDof(), RThNonPeriodic->nDof() );
    checkPeriodicBoundaryEdgesShareDofs( RThPeriodic, meshPeriodic, "Right", "Left" );
    checkPeriodicBoundaryEdgesShareDofs( RThPeriodic, meshPeriodic, "Top", "Bottom" );

    auto BDMhPeriodic = BDMh<0>( meshPeriodic );
    auto BDMhNonPeriodic = BDMh<0>( meshNonPeriodic );
    auto const bdmEdgeDofsPeriodic = collectEdgeDofSetsByEdge( BDMhPeriodic );
    auto const bdmEdgeDofsNonPeriodic = collectEdgeDofSetsByEdge( BDMhNonPeriodic );
    auto const bdmDofPerEdgePeriodic = edgeDofCardinality( bdmEdgeDofsPeriodic );
    auto const bdmDofPerEdgeNonPeriodic = edgeDofCardinality( bdmEdgeDofsNonPeriodic );
    BOOST_CHECK_EQUAL( bdmDofPerEdgePeriodic, bdmDofPerEdgeNonPeriodic );
    BOOST_REQUIRE_GE( BDMhNonPeriodic->nDof(), nonPeriodicCanonicalEdges * bdmDofPerEdgeNonPeriodic );
    auto const bdmEdgeInvariantOffset = BDMhNonPeriodic->nDof() - nonPeriodicCanonicalEdges * bdmDofPerEdgeNonPeriodic;
    BOOST_CHECK_EQUAL( BDMhPeriodic->nDof(), bdmEdgeInvariantOffset + periodicCanonicalEdges * bdmDofPerEdgePeriodic );
    BOOST_CHECK_LT( BDMhPeriodic->nDof(), BDMhNonPeriodic->nDof() );
    checkPeriodicBoundaryEdgesShareDofs( BDMhPeriodic, meshPeriodic, "Right", "Left" );
    checkPeriodicBoundaryEdgesShareDofs( BDMhPeriodic, meshPeriodic, "Top", "Bottom" );

    auto NedPeriodic = Ned1h<0>( meshPeriodic );
    auto NedNonPeriodic = Ned1h<0>( meshNonPeriodic );
    auto const nedEdgeDofsPeriodic = collectEdgeDofSetsByEdge( NedPeriodic );
    auto const nedEdgeDofsNonPeriodic = collectEdgeDofSetsByEdge( NedNonPeriodic );
    auto const nedDofPerEdgePeriodic = edgeDofCardinality( nedEdgeDofsPeriodic );
    auto const nedDofPerEdgeNonPeriodic = edgeDofCardinality( nedEdgeDofsNonPeriodic );
    BOOST_CHECK_EQUAL( nedDofPerEdgePeriodic, nedDofPerEdgeNonPeriodic );
    BOOST_REQUIRE_GE( NedNonPeriodic->nDof(), nonPeriodicCanonicalEdges * nedDofPerEdgeNonPeriodic );
    auto const nedEdgeInvariantOffset = NedNonPeriodic->nDof() - nonPeriodicCanonicalEdges * nedDofPerEdgeNonPeriodic;
    BOOST_CHECK_EQUAL( NedPeriodic->nDof(), nedEdgeInvariantOffset + periodicCanonicalEdges * nedDofPerEdgePeriodic );
    BOOST_CHECK_LT( NedPeriodic->nDof(), NedNonPeriodic->nDof() );
    checkPeriodicBoundaryEdgesShareDofs( NedPeriodic, meshPeriodic, "Right", "Left" );
    checkPeriodicBoundaryEdgesShareDofs( NedPeriodic, meshPeriodic, "Top", "Bottom" );
}

BOOST_AUTO_TEST_SUITE_END()
