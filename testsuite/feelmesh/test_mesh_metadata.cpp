/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
#define BOOST_TEST_MODULE mesh_metadata
#include <feel/feelcore/testsuite.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/partitionio.hpp>
#include <feel/feelfilters/partitionmetadata.hpp>
#include <feel/feelmesh/compactmarkerfragmentation.hpp>

#include <boost/mpi/collectives/all_gather.hpp>
#include <sstream>

using namespace Feel;
FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
template <typename Value>
std::vector<Value> copyValues( std::span<Value const> values )
{
    return { values.begin(), values.end() };
}

Marker<> marker( std::initializer_list<flag_type> values )
{
    Marker<> result;
    result.assign( values );
    return result;
}
}

BOOST_AUTO_TEST_CASE( streamingMetadataAndSparseFragments )
{
    std::map<std::string, std::pair<int, int>> physicals;
    Feel::detail::PartitionMetadataReader<flag_type> reader(
        [&]( auto const& name, int id, int dimension ) { physicals[name] = { id, dimension }; } );
    std::istringstream input( R"({
        "extension": {"ignored": [1, {"mesh": 42}]},
        "mesh": {
            "fragmentation": {"0": {"1000000000": [7, 1, 7], "0": [], "10": [9999999999]}, "1": {}},
            "h5": "mesh.h5", "partition": {"n": "3"},
            "physicals": {"wall/with~escape": ["7", 1], "roof": [8, "2"], "old": "9 0"}
        }})" );
    reader.read( input );
    BOOST_CHECK_EQUAL( *reader.M_h5Filename, "mesh.h5" );
    BOOST_CHECK_EQUAL( *reader.M_partitions, 3 );
    BOOST_CHECK( reader.M_hasFragmentation );
    BOOST_CHECK( physicals.at( "wall/with~escape" ) == std::make_pair( 7, 1 ) );
    BOOST_CHECK( physicals.at( "roof" ) == std::make_pair( 8, 2 ) );
    BOOST_CHECK( physicals.at( "old" ) == std::make_pair( 9, 0 ) );
    auto const& table = reader.M_fragments.at( 0 );
    BOOST_CHECK_EQUAL( table.size(), 3 );
    BOOST_CHECK( table.at( 0 ).empty() );
    BOOST_CHECK( copyValues( table.at( 1000000000 ) ) == ( std::vector<flag_type>{7, 1, 7} ) );
    BOOST_CHECK_EQUAL( table.at( 10 )[0], 9999999999LL );
    BOOST_CHECK_THROW( table.at( 11 ), std::out_of_range );
    BOOST_CHECK_LT( table.storageBytes(), 1024 );
}

BOOST_AUTO_TEST_CASE( legacyAndInvalidMetadata )
{
    Feel::detail::PartitionMetadataReader<flag_type> legacy( []( auto const&, int, int ) {} );
    std::istringstream input( R"({"mesh":{"physicals":{"alias":"12 2"},"partition":{"n":4}}})" );
    legacy.read( input );
    BOOST_CHECK( !legacy.M_hasFragmentation );
    BOOST_CHECK_EQUAL( *legacy.M_partitions, 4 );
    for ( std::string invalid : {
            R"({"mesh":{"physicals":{"a":[1]}}})",
            R"({"mesh":{"physicals":{"a":[1,2,3]}}})",
            R"({"mesh":{"fragmentation":{"0":{"1":["2"]}}}})",
            R"({"mesh":{"fragmentation":{"0":{"1":[{}]}}}})",
            R"({"mesh":{"partition":{"n":[]}}})",
            R"({"mesh":[])",
            R"({})" } )
    {
        Feel::detail::PartitionMetadataReader<flag_type> reader( []( auto const&, int, int ) {} );
        std::istringstream malformed( invalid );
        BOOST_CHECK_THROW( reader.read( malformed ), std::exception );
    }
    Feel::detail::PartitionMarkerTable<flag_type> aliases;
    aliases.append( 12, {12} );
    aliases.append( 12, {12} );
    aliases.finalize();
    BOOST_CHECK_EQUAL( aliases.size(), 1 );
}

BOOST_AUTO_TEST_CASE( compactMergePreservesLexicographicOrder )
{
    Feel::detail::MergeMarkerBuffers<flag_type> merge;
    // Empty set, [1,2], [2] versus [1], [1,2], [3]: length must not sort first.
    std::vector<flag_type> a{0, 2, 1, 2, 1, 2};
    std::vector<flag_type> b{1, 1, 2, 1, 2, 1, 3};
    std::vector<flag_type> expected{0, 1, 1, 2, 1, 2, 1, 2, 1, 3};
    BOOST_CHECK( merge( a, b ) == expected );
    BOOST_CHECK( merge( b, a ) == expected );
    BOOST_CHECK( merge( a, a ) == a );
    BOOST_CHECK( merge( {}, b ) == b );
    BOOST_CHECK( merge( a, {} ) == a );
    BOOST_CHECK( merge( merge( a, b ), {1, -4} ) == merge( a, merge( b, {1, -4} ) ) );
}

BOOST_AUTO_TEST_CASE( mpiUnionMatchesTreeSetReference )
{
    auto const& comm = Environment::worldComm().localComm();
    for ( int pattern = 0; pattern < 3; ++pattern )
    {
        std::set<Marker<>> local;
        if ( pattern == 0 || ( pattern == 1 && comm.rank() == comm.size() - 1 ) )
        {
            local.insert( marker( {} ) );
            local.insert( marker( {1} ) );
            local.insert( marker( {1, 2} ) );
            local.insert( marker( {2, flag_type( 100 + comm.rank() )} ) );
            local.insert( marker( {-9, 9999999999LL} ) );
        }
        std::vector<std::set<Marker<>>> gathered;
        boost::mpi::all_gather( comm, local, gathered );
        std::set<Marker<>> reference;
        for ( auto const& contribution : gathered )
            reference.insert( contribution.begin(), contribution.end() );
        auto actual = Feel::detail::globalMarkerFragments( local, comm );
        BOOST_CHECK( local.empty() );
        BOOST_REQUIRE_EQUAL( actual.size(), reference.size() );
        int id = 0;
        for ( auto const& values : reference )
            BOOST_CHECK( actual.at( id++ ) == values );
    }
}

BOOST_AUTO_TEST_CASE( mpiUnionCrossesBroadcastChunkBoundary )
{
    auto const& comm = Environment::worldComm().localComm();
    constexpr int count = 530000; // Two primitive broadcast chunks.
    std::set<Marker<>> local;
    if ( comm.rank() == 0 )
        for ( int id = 0; id < count; ++id )
            local.emplace_hint( local.end(), marker( {id} ) );
    auto actual = Feel::detail::globalMarkerFragments( local, comm );
    BOOST_REQUIRE_EQUAL( actual.size(), count );
    BOOST_CHECK( local.empty() );
    BOOST_CHECK( actual.at( 0 ) == marker( {0} ) );
    BOOST_CHECK( actual.at( count - 1 ) == marker( {count - 1} ) );
}

#if defined( FEELPP_HAS_HDF5 )
BOOST_AUTO_TEST_CASE( partitionRoundTripPreservesEntityMarkers )
{
    using MeshType = Mesh<Simplex<2, 1, 3>>;
    auto mesh = std::make_shared<MeshType>();
    int rank = mesh->worldComm().localRank();
    for ( int i = 0; i < 3; ++i )
    {
        MeshType::point_type point( 3 * rank + i );
        point( 0 ) = 3 * rank + ( i == 1 );
        point( 1 ) = ( i == 2 );
        point( 2 ) = 0;
        point.setProcessId( rank );
        point.setProcessIdInPartition( rank );
        if ( i == 0 )
            point.setMarker( 99 );
        mesh->addPoint( std::move( point ) );
    }
    MeshType::element_type element;
    element.setId( rank );
    element.setProcessId( rank );
    element.setProcessIdInPartition( rank );
    element.setMarker( marker( {1, flag_type( 100 + rank )} ) );
    for ( int i = 0; i < 3; ++i )
        element.setPoint( i, mesh->point( 3 * rank + i ) );
    mesh->addElement( std::move( element ) );
    MeshType::face_type face;
    face.setId( rank );
    face.setProcessIdInPartition( rank );
    face.setProcessId( rank );
    face.setMarker( marker( {4, 5} ) );
    face.setPoint( 0, mesh->point( 3 * rank ) );
    face.setPoint( 1, mesh->point( 3 * rank + 1 ) );
    mesh->addFace( std::move( face ) );
    mesh->addMarkerName( std::make_pair( "surface", std::vector<MeshType::size_type>{1, 2} ) );
    mesh->components().set( MESH_UPDATE_FACES | MESH_UPDATE_EDGES );
    mesh->updateForUse();

    for ( bool transpose : { false, true } )
    {
        auto path = ( fs::current_path() / ( transpose ? "metadata_transposed.json" : "metadata.json" ) ).string();
        PartitionIO<MeshType> writer( path, transpose );
        writer.write( mesh );
        PartitionIO<MeshType> reader( path, transpose );
        auto reloaded = std::make_shared<MeshType>();
        reader.read( reloaded );
        BOOST_CHECK_EQUAL( reloaded->numGlobalElements(), mesh->numGlobalElements() );
        BOOST_CHECK( reloaded->markerNames() == mesh->markerNames() );
        BOOST_CHECK( reloaded->meshFragmentationByMarkerByEntity() == mesh->meshFragmentationByMarkerByEntity() );
        BOOST_REQUIRE_EQUAL( reloaded->elements().size(), mesh->elements().size() );
        BOOST_CHECK( reloaded->elements().begin()->second.marker() == marker( {1, flag_type( 100 + rank )} ) );
    }
}
#endif
