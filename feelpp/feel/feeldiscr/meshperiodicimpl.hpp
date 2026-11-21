/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

   This file is part of the Feel++ library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2026-02-19

   Copyright (C) 2026 Feel++ Consortium

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#ifndef FEELPP_MESHPERIODICIMPL_HPP
#define FEELPP_MESHPERIODICIMPL_HPP 1

namespace Feel
{

template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void
Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::invalidatePeriodicCanonicalMaps()
{
    M_hasPeriodicCanonicalMaps = false;
    M_periodicCanonicalPointId.clear();
    M_periodicCanonicalEdgeId.clear();
    M_periodicCanonicalFaceId.clear();
    M_periodicPointMasterId.clear();
    M_periodicEdgeMasterId.clear();
    M_periodicEdgeMasterOrientation.clear();
    M_periodicFaceMasterId.clear();
    M_periodicFaceMasterPermutation.clear();
}

template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void
Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::buildPeriodicCanonicalMaps()
{
    this->invalidatePeriodicCanonicalMaps();

    for ( auto const& pit : this->points() )
        M_periodicCanonicalPointId.emplace( pit.first, pit.first );

    if ( !this->isPeriodic() )
    {
        M_hasPeriodicCanonicalMaps = true;
        return;
    }

    std::unordered_map<size_type, size_type> directMaster;
    for ( auto const& periodicEntity : M_periodic_entities )
    {
        for ( auto const& [slaveId, masterId] : periodicEntity.correspondingVertices )
        {
            if ( slaveId == masterId )
                continue;
            if ( this->hasPoint( slaveId ) && this->hasPoint( masterId ) )
            {
                directMaster[slaveId] = masterId;
                M_periodicPointMasterId[slaveId] = masterId;
            }
        }
    }

    for ( auto const& pit : this->points() )
    {
        size_type cur = pit.first;
        std::vector<size_type> path;
        while ( true )
        {
            auto it = directMaster.find( cur );
            if ( it == directMaster.end() || it->second == cur )
                break;
            path.push_back( cur );
            cur = it->second;
            if ( path.size() > directMaster.size() )
            {
                cur = *std::min_element( path.begin(), path.end() );
                break;
            }
        }
        for ( size_type id : path )
            M_periodicCanonicalPointId[id] = cur;
        M_periodicCanonicalPointId[pit.first] = cur;
    }

    auto canonicalPointFromCache = [this]( size_type pointId ) -> size_type
    {
        auto it = M_periodicCanonicalPointId.find( pointId );
        return ( it != M_periodicCanonicalPointId.end() ) ? it->second : pointId;
    };

    if constexpr ( nDim >= 2 )
    {
        std::map<std::pair<size_type, size_type>, size_type> edgeByVertexPair;
        std::unordered_map<size_type, size_type> edgeParent;
        std::unordered_map<size_type, std::pair<size_type, size_type>> edgeVertices;
        std::unordered_map<size_type, std::array<size_type, 2>> edgeVerticesOrdered;
        std::unordered_map<size_type, std::vector<size_type>> incidentEdgesByVertex;
        for ( auto it = this->beginEdge(), en = this->endEdge(); it != en; ++it )
        {
            auto const& edge = it->second;
            constexpr uint16_type nEdgeVertices = edge_type::numVertices;
            if constexpr ( nEdgeVertices >= 2 )
            {
                size_type const vo0 = edge.point( 0 ).id();
                size_type const vo1 = edge.point( 1 ).id();
                size_type v0 = vo0;
                size_type v1 = vo1;
                if ( v1 < v0 )
                    std::swap( v0, v1 );
                edgeByVertexPair.emplace( std::make_pair( v0, v1 ), edge.id() );
                edgeParent.emplace( edge.id(), edge.id() );
                edgeVertices.emplace( edge.id(), std::make_pair( v0, v1 ) );
                edgeVerticesOrdered.emplace( edge.id(), std::array<size_type, 2>{ vo0, vo1 } );
                incidentEdgesByVertex[v0].push_back( edge.id() );
                incidentEdgesByVertex[v1].push_back( edge.id() );
            }
        }

        auto findEdgeRoot = [&edgeParent]( size_type edgeId ) -> size_type
        {
            std::vector<size_type> path;
            size_type cur = edgeId;
            while ( true )
            {
                auto it = edgeParent.find( cur );
                if ( it == edgeParent.end() || it->second == cur )
                    break;
                if ( std::find( path.begin(), path.end(), cur ) != path.end() )
                {
                    cur = *std::min_element( path.begin(), path.end() );
                    break;
                }
                path.push_back( cur );
                cur = it->second;
            }
            for ( size_type id : path )
                edgeParent[id] = cur;
            return cur;
        };
        auto unifyEdgesDirected = [&edgeParent, &findEdgeRoot]( size_type slaveEdge, size_type masterEdge )
        {
            size_type const slaveRoot = findEdgeRoot( slaveEdge );
            size_type const masterRoot = findEdgeRoot( masterEdge );
            if ( slaveRoot == masterRoot )
                return;
            edgeParent[slaveRoot] = masterRoot;
        };

        for ( auto const& periodicEntity : M_periodic_entities )
        {
            if ( periodicEntity.dim != 1 || periodicEntity.correspondingVertices.empty() )
                continue;

            auto const& correspondingVertices = periodicEntity.correspondingVertices;
            std::unordered_set<size_type> visitedSlaveEdges;
            for ( auto const& slaveToMaster : correspondingVertices )
            {
                auto const slaveVertexId = slaveToMaster.first;
                auto itIncidentEdges = incidentEdgesByVertex.find( slaveVertexId );
                if ( itIncidentEdges == incidentEdgesByVertex.end() )
                    continue;

                for ( auto const slaveEdgeId : itIncidentEdges->second )
                {
                    if ( !visitedSlaveEdges.insert( slaveEdgeId ).second )
                        continue;

                    auto itEdgeVertices = edgeVertices.find( slaveEdgeId );
                    if ( itEdgeVertices == edgeVertices.end() )
                        continue;

                    auto const slave0 = itEdgeVertices->second.first;
                    auto const slave1 = itEdgeVertices->second.second;

                    auto itSlave0 = correspondingVertices.find( slave0 );
                    auto itSlave1 = correspondingVertices.find( slave1 );
                    if ( itSlave0 == correspondingVertices.end() || itSlave1 == correspondingVertices.end() )
                        continue;

                    size_type master0 = itSlave0->second;
                    size_type master1 = itSlave1->second;
                    if ( master1 < master0 )
                        std::swap( master0, master1 );

                    auto itMasterEdge = edgeByVertexPair.find( std::make_pair( master0, master1 ) );
                    if ( itMasterEdge == edgeByVertexPair.end() )
                        continue;

                    size_type const masterEdgeId = itMasterEdge->second;
                    if ( slaveEdgeId != masterEdgeId )
                    {
                        M_periodicEdgeMasterId[slaveEdgeId] = masterEdgeId;

                        auto const& slaveOrdered = edgeVerticesOrdered.at( slaveEdgeId );
                        auto const& masterOrdered = edgeVerticesOrdered.at( masterEdgeId );
                        auto itMapped0 = correspondingVertices.find( slaveOrdered[0] );
                        auto itMapped1 = correspondingVertices.find( slaveOrdered[1] );
                        if ( itMapped0 != correspondingVertices.end() &&
                             itMapped1 != correspondingVertices.end() )
                        {
                            size_type const mapped0 = itMapped0->second;
                            size_type const mapped1 = itMapped1->second;
                            int orientation = 1;
                            if ( mapped0 == masterOrdered[1] && mapped1 == masterOrdered[0] )
                                orientation = -1;
                            M_periodicEdgeMasterOrientation[slaveEdgeId] = orientation;
                        }
                    }

                    unifyEdgesDirected( slaveEdgeId, masterEdgeId );
                }
            }
        }

        for ( auto it = this->beginEdge(), en = this->endEdge(); it != en; ++it )
        {
            auto const& edge = it->second;
            constexpr uint16_type nEdgeVertices = edge_type::numVertices;
            if constexpr ( nEdgeVertices >= 2 )
                M_periodicCanonicalEdgeId[edge.id()] = findEdgeRoot( edge.id() );
            else
                M_periodicCanonicalEdgeId[edge.id()] = edge.id();
        }
    }

    if constexpr ( nDim == 3 )
    {
        std::map<std::vector<size_type>, size_type> faceByVertexSet;
        std::unordered_map<size_type, std::vector<size_type>> faceVertices;
        std::unordered_map<size_type, std::vector<size_type>> incidentFacesByVertex;
        std::unordered_map<size_type, size_type> faceParent;
        for ( auto it = this->beginFace(), en = this->endFace(); it != en; ++it )
        {
            auto const& face = it->second;
            std::vector<size_type> vertices;
            vertices.reserve( face_type::numVertices );
            for ( uint16_type i = 0; i < face_type::numVertices; ++i )
                vertices.push_back( face.point( i ).id() );
            faceVertices.emplace( face.id(), vertices );
            faceParent.emplace( face.id(), face.id() );
            for ( auto const vertexId : vertices )
                incidentFacesByVertex[vertexId].push_back( face.id() );

            auto key = vertices;
            std::sort( key.begin(), key.end() );
            faceByVertexSet.emplace( key, face.id() );
        }

        auto findFaceRoot = [&faceParent]( size_type faceId ) -> size_type
        {
            std::vector<size_type> path;
            size_type cur = faceId;
            while ( true )
            {
                auto it = faceParent.find( cur );
                if ( it == faceParent.end() || it->second == cur )
                    break;
                if ( std::find( path.begin(), path.end(), cur ) != path.end() )
                {
                    cur = *std::min_element( path.begin(), path.end() );
                    break;
                }
                path.push_back( cur );
                cur = it->second;
            }
            for ( size_type id : path )
                faceParent[id] = cur;
            return cur;
        };
        auto unifyFacesDirected = [&faceParent, &findFaceRoot]( size_type slaveFace, size_type masterFace )
        {
            size_type const slaveRoot = findFaceRoot( slaveFace );
            size_type const masterRoot = findFaceRoot( masterFace );
            if ( slaveRoot == masterRoot )
                return;
            faceParent[slaveRoot] = masterRoot;
        };

        std::unordered_set<size_type> explicitMasterFaces;
        for ( auto const& periodicEntity : M_periodic_entities )
        {
            if ( periodicEntity.dim != 2 || periodicEntity.correspondingVertices.empty() )
                continue;

            auto const& correspondingVertices = periodicEntity.correspondingVertices;
            std::unordered_set<size_type> visitedSlaveFaces;
            for ( auto const& [slaveVertexId, masterVertexId] : correspondingVertices )
            {
                Feel::detail::ignore_unused_variable_warning( masterVertexId );
                auto itIncidentFaces = incidentFacesByVertex.find( slaveVertexId );
                if ( itIncidentFaces == incidentFacesByVertex.end() )
                    continue;

                for ( auto const slaveFaceId : itIncidentFaces->second )
                {
                    if ( !visitedSlaveFaces.insert( slaveFaceId ).second )
                        continue;

                    auto const& slaveVertices = faceVertices.at( slaveFaceId );
                    std::vector<size_type> mappedVertices;
                    mappedVertices.reserve( slaveVertices.size() );
                    bool mappedAllVertices = true;
                    for ( auto const slaveFaceVertexId : slaveVertices )
                    {
                        auto itMappedVertex = correspondingVertices.find( slaveFaceVertexId );
                        if ( itMappedVertex == correspondingVertices.end() )
                        {
                            mappedAllVertices = false;
                            break;
                        }
                        mappedVertices.push_back( itMappedVertex->second );
                    }
                    if ( !mappedAllVertices )
                        continue;

                    auto masterKey = mappedVertices;
                    std::sort( masterKey.begin(), masterKey.end() );
                    auto itMasterFace = faceByVertexSet.find( masterKey );
                    if ( itMasterFace == faceByVertexSet.end() )
                        continue;

                    size_type const masterFaceId = itMasterFace->second;
                    if ( slaveFaceId == masterFaceId )
                        continue;

                    M_periodicFaceMasterId[slaveFaceId] = masterFaceId;
                    explicitMasterFaces.insert( masterFaceId );
                    unifyFacesDirected( slaveFaceId, masterFaceId );

                    auto const& masterVertices = faceVertices.at( masterFaceId );
                    std::vector<uint16_type> permutation( slaveVertices.size(), uint16_type( 0 ) );
                    std::vector<bool> usedMasterVertex( masterVertices.size(), false );
                    bool validPermutation = true;
                    for ( size_type i = 0; i < mappedVertices.size(); ++i )
                    {
                        auto itPos = std::find( masterVertices.begin(), masterVertices.end(), mappedVertices[i] );
                        if ( itPos == masterVertices.end() )
                        {
                            validPermutation = false;
                            break;
                        }
                        auto const pos = static_cast<size_type>( std::distance( masterVertices.begin(), itPos ) );
                        if ( pos >= usedMasterVertex.size() || usedMasterVertex[pos] )
                        {
                            validPermutation = false;
                            break;
                        }
                        usedMasterVertex[pos] = true;
                        permutation[i] = static_cast<uint16_type>( pos );
                    }
                    if ( validPermutation )
                        M_periodicFaceMasterPermutation[slaveFaceId] = std::move( permutation );
                }
            }
        }

        std::map<std::vector<size_type>, size_type> canonicalFaceKeys;
        for ( auto it = this->beginFace(), en = this->endFace(); it != en; ++it )
        {
            auto const& face = it->second;
            std::vector<size_type> key;
            key.reserve( face_type::numVertices );
            for ( uint16_type i = 0; i < face_type::numVertices; ++i )
                key.push_back( canonicalPointFromCache( face.point( i ).id() ) );

            std::sort( key.begin(), key.end() );
            bool const hasDuplicateCanonicalVertex =
                std::adjacent_find( key.begin(), key.end() ) != key.end();
            CHECK( !hasDuplicateCanonicalVertex )
                << "Degenerate periodic canonical face key for face id " << face.id();

            size_type const rootFaceId = findFaceRoot( face.id() );
            auto [itKey, inserted] = canonicalFaceKeys.emplace( key, rootFaceId );
            if ( !inserted && itKey->second != rootFaceId )
            {
                bool const existingIsMaster = explicitMasterFaces.contains( itKey->second );
                bool const currentIsMaster = explicitMasterFaces.contains( rootFaceId );
                if ( !existingIsMaster && currentIsMaster )
                    itKey->second = rootFaceId;
                else if ( !existingIsMaster && !currentIsMaster )
                    itKey->second = std::min( itKey->second, rootFaceId );
            }
            M_periodicCanonicalFaceId[face.id()] = itKey->second;
        }
    }

    M_hasPeriodicCanonicalMaps = true;
}

} // namespace Feel

#endif
