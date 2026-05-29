/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
  This file is part of the Feel library

  Copyright (C) 2012 Université de Grenoble 1

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
 * \file doftablempi.hpp
 * \author Vincent Chabannes
 */
#if !defined(FEELPP_DOFTABLE_MPI_HPP)
#define FEELPP_DOFTABLE_MPI_HPP 1

namespace Feel
{

template<typename MeshType, typename FEType, typename MortarType>
void
DofTable<MeshType, FEType, MortarType>::buildGhostDofMap( mesh_type& mesh )
{
    wc(mesh)->print(fmt::format("[DofTable::buildGhostDofMap rank={}] starts. hasMeshSupport: {}", rank(mesh), this->hasMeshSupport()), Environment::logVerbosityLevel() > 1, Environment::logVerbosityLevel() > 0, Environment::logVerbosityLevel() > 1 );

    this->buildGlobalProcessToGlobalClusterDofMapOthersMesh( mesh );

    DVLOG(2) << "[buildGhostDofMap] finish () with rank "<< this->worldComm().rank();
}


//--------------------------------------------------------------------------------------------------------//

template<typename MeshType, typename FEType, typename MortarType>
void
DofTable<MeshType, FEType, MortarType>::generateDofPoints( Range<mesh_type,MESH_ELEMENTS> const& myrange ) const
{
    if ( fe_type::is_modal )
        return;

    DVLOG(2) << "[Dof::generateDofPoints] generating dof coordinates\n";
    typedef typename gm_type::template Context<element_type> gm_context_type;
    typedef std::shared_ptr<gm_context_type> gm_context_ptrtype;

    typedef typename fe_type::template Context<vm::POINT, fe_type, gm_type, element_type> fecontext_type;

    gm_ptrtype gm( new gm_type );
    auto const& fe = this->fe();

    // Precompute some data in the reference element for
    // geometric mapping and reference finite element
    typename gm_type::precompute_ptrtype __geopc( new typename gm_type::precompute_type( gm, fe.points() ) );

    gm_context_ptrtype __c;

    std::vector<bool> dof_done( this->nLocalDofWithGhost(), false );

    for (auto const& eltWrap : myrange )
    {
        auto const& elt = unwrap_ref( eltWrap );

        if ( __c )
            __c->template update<vm::POINT>( elt );
        else
            __c = gm->template context<vm::POINT>( elt, __geopc );

        for ( auto const& ldof : this->localDof( elt.id() ) )
        {
            size_type thedof = ldof.second.index();
            uint16_type ldofId = ldof.first.localDof();
            if ( !this->localDofHasRepresentativePoint( ldofId ) )
                continue;

            if ( dof_done[thedof] )
                continue;
            dof_done[thedof] = true;

            uint16_type pointId = this->localDofRepresentativePointIndex( ldofId );
            FEELPP_ASSERT( pointId < static_cast<uint16_type>( fe.points().size2() ) )
                ( int( pointId ) )( int( fe.points().size2() ) )( int( ldofId ) )
                .error( "invalid FE representative point index" );

            if ( ( thedof >= this->firstDof() ) && ( thedof <= this->lastDof() ) )
            {
                DCHECK( thedof < this->nLocalDofWithGhost() )
                    << "invalid local dof index "
                    <<  thedof << ", " << this->nLocalDofWithGhost() << "," << this->firstDof()  << ","
                    <<  this->lastDof() << "," << elt.id() << "," << ldofId << "," << pointId;

                uint16_type comp = fe.component( ldofId );
                M_dof_points[thedof] = boost::make_tuple( __c->xReal( pointId ), thedof, comp );
            }
        }
    }

}

template<typename MeshType, typename FEType, typename MortarType>
bool
DofTable<MeshType, FEType, MortarType>::buildGlobalProcessToGlobalClusterDofMapDescriptorKeys( mesh_type& mesh )
{
    if constexpr ( !FiniteElementDofLayoutProvider<fe_type> || fe_type::is_modal ||
                   !is_continuous || is_tensor2symm || is_product )
        return false;
    else
    {
        if ( mesh.isPeriodic() )
            return false;

        if ( this->hasMeshSupport() && this->meshSupport()->isPartialSupport() )
            return false;

        if ( this->worldComm().localSize() <= 1 )
            return false;

        auto const& fe = this->fe();
        const uint16_type nLocalDof = runtimeNLocalDof();
        for ( uint16_type parentLid = 0; parentLid < nLocalDof; ++parentLid )
            if ( !fe.localDofLayout( parentLid ).attachment.isValid() )
                return false;

        const uint16_type nDofPerVertex = runtimeDofPerVertex();
        const uint16_type nDofPerEdge = runtimeDofPerEdge();
        const uint16_type nDofPerFace = runtimeDofPerFace();
        const uint16_type nDofPerVolume = runtimeDofPerVolume();
        const rank_type myRank = this->worldComm().localRank();
        const rank_type nProc = this->worldComm().localSize();
        const size_type nLocalDofWithGhost = this->M_n_localWithGhost_df[myRank];

        using descriptor_key_type = DofKey<size_type>;

        auto familyTag = []() constexpr -> uint16_type
        {
            if constexpr ( requires { fe_type::TAG; } )
                return fe_type::TAG;
            else
                return 0;
        };

        auto makeKey = [&familyTag]( uint16_type topologicalDim,
                                     rank_type canonicalPartition,
                                     size_type canonicalEntityId,
                                     uint16_type ordinal,
                                     uint16_type component,
                                     uint16_type functional )
        {
            return descriptor_key_type{
                .topologicalDim = static_cast<uint8_type>( topologicalDim ),
                .canonicalPartition = canonicalPartition,
                .canonicalEntityId = canonicalEntityId,
                .ordinal = ordinal,
                .component = component,
                .familyTag = familyTag(),
                .variant = functional
            };
        };

        auto packKey = []( descriptor_key_type const& key, std::vector<size_type>& payload )
        {
            payload.push_back( static_cast<size_type>( key.topologicalDim ) );
            payload.push_back( static_cast<size_type>( key.canonicalPartition ) );
            payload.push_back( key.canonicalEntityId );
            payload.push_back( static_cast<size_type>( key.ordinal ) );
            payload.push_back( static_cast<size_type>( key.component ) );
            payload.push_back( static_cast<size_type>( key.familyTag ) );
            payload.push_back( static_cast<size_type>( key.variant ) );
        };

        auto unpackKey = []( std::vector<size_type> const& payload, size_type offset )
        {
            return descriptor_key_type{
                .topologicalDim = static_cast<uint8_type>( payload[offset] ),
                .canonicalPartition = static_cast<rank_type>( payload[offset+1] ),
                .canonicalEntityId = payload[offset+2],
                .ordinal = static_cast<uint16_type>( payload[offset+3] ),
                .component = static_cast<uint16_type>( payload[offset+4] ),
                .familyTag = static_cast<uint16_type>( payload[offset+5] ),
                .variant = static_cast<uint16_type>( payload[offset+6] )
            };
        };

        auto entityHasOtherPartitions = []( auto const& entity )
        {
            return !entity.idInOthersPartitions().empty();
        };

        auto partitionCanonicalEntity = []( auto const& entity )
        {
            rank_type canonicalPartition = entity.pidInPartition();
            if ( canonicalPartition == invalid_rank_type_value )
                canonicalPartition = entity.processId();
            size_type canonicalId = entity.id();
            if ( canonicalPartition != entity.pidInPartition() )
            {
                auto const itId = entity.idInOthersPartitions().find( canonicalPartition );
                if ( itId != entity.idInOthersPartitions().end() && itId->second != invalid_v<size_type> )
                    canonicalId = itId->second;
            }
            for ( auto const& [pid,entityId] : entity.idInOthersPartitions() )
            {
                if ( entityId != invalid_v<size_type> && pid < canonicalPartition )
                {
                    canonicalPartition = pid;
                    canonicalId = entityId;
                }
            }
            return std::make_pair( canonicalPartition, canonicalId );
        };

        auto makeElementDofKey = [&]( element_type const& elt,
                                      uint16_type parentLid,
                                      descriptor_key_type& key,
                                      bool& isShared ) -> bool
        {
            auto const layout = fe.localDofLayout( parentLid );
            auto const& attachment = layout.attachment;
            if ( !attachment.isValid() )
                return false;

            isShared = false;
            switch ( attachment.entityDim )
            {
            case 0:
            {
                if ( nDofPerVertex == 0 || attachment.entityId >= element_type::numVertices ||
                     attachment.ordinal >= nDofPerVertex )
                    return false;

                auto const& point = elt.point( attachment.entityId );
                isShared = mesh.isInterprocessPoints( point.id() ) || entityHasOtherPartitions( point );
                auto const [entityPartition,entityId] = partitionCanonicalEntity( point );
                key = makeKey( 0, entityPartition, entityId,
                               attachment.ordinal, layout.component, attachment.kind );
                return true;
            }
            case 1:
            {
                if ( nDofPerEdge == 0 || attachment.ordinal >= nDofPerEdge )
                    return false;

                if constexpr ( nDim == 1 )
                {
                    key = makeKey( 1, myRank, elt.id(), attachment.ordinal, layout.component, attachment.kind );
                    return true;
                }
                else
                {
                    if ( attachment.entityId >= element_type::numEdges )
                        return false;

                    uint16_type ordinal = attachment.ordinal;
                    auto const edgePermutation = elt.edgePermutation( attachment.entityId );
                    if ( edgePermutation.value() == edge_permutation_type::REVERSE_PERMUTATION )
                        ordinal = static_cast<uint16_type>( nDofPerEdge - 1 - attachment.ordinal );
                    else if ( edgePermutation.value() != edge_permutation_type::IDENTITY )
                        return false;

                    if constexpr ( nDim == 2 )
                    {
                        auto const facePtr = elt.facePtr( attachment.entityId );
                        if ( !facePtr )
                            return false;

                        auto const& edge = *facePtr;
                        isShared = edge.isInterProcessDomain() || entityHasOtherPartitions( edge );
                        auto const [entityPartition,entityId] = partitionCanonicalEntity( edge );
                        key = makeKey( 1, entityPartition, entityId,
                                       ordinal, layout.component, attachment.kind );
                        return true;
                    }
                    else
                    {
                        auto const edgePtr = elt.edgePtr( attachment.entityId );
                        if ( !edgePtr )
                            return false;

                        auto const& edge = *edgePtr;
                        isShared = mesh.isInterprocessEdges( edge.id() ) || entityHasOtherPartitions( edge );
                        auto const [entityPartition,entityId] = partitionCanonicalEntity( edge );
                        key = makeKey( 1, entityPartition, entityId,
                                       ordinal, layout.component, attachment.kind );
                        return true;
                    }
                }
            }
            case 2:
            {
                if ( nDofPerFace == 0 || attachment.ordinal >= nDofPerFace )
                    return false;

                if constexpr ( nDim == 2 )
                {
                    key = makeKey( 2, myRank, elt.id(), attachment.ordinal, layout.component, attachment.kind );
                    return true;
                }
                else if constexpr ( nDim == 3 )
                {
                    if ( attachment.entityId >= element_type::numFaces )
                        return false;

                    auto const facePtr = elt.facePtr( attachment.entityId );
                    if ( !facePtr )
                        return false;

                    auto const& face = *facePtr;
                    uint16_type ordinal = attachment.ordinal;
                    auto const facePermutation = elt.facePermutation( attachment.entityId );
                    if ( facePermutation == face_permutation_type( 0 ) )
                        return false;
                    if ( nDofPerFace != 1 && facePermutation != face_permutation_type( face_permutation_type::IDENTITY ) )
                    {
                        if ( !this->hasValidFacePermutation( facePermutation, nDofPerFace ) )
                            return false;
                        auto const& permutation = this->facePermutationVector( facePermutation, nDofPerFace );
                        ordinal = static_cast<uint16_type>( permutation( attachment.ordinal ) );
                    }

                    isShared = face.isInterProcessDomain() || entityHasOtherPartitions( face );
                    auto const [entityPartition,entityId] = partitionCanonicalEntity( face );
                    key = makeKey( 2, entityPartition, entityId,
                                   ordinal, layout.component, attachment.kind );
                    return true;
                }
                else
                    return false;
            }
            case 3:
            {
                if ( nDofPerVolume == 0 || attachment.ordinal >= nDofPerVolume )
                    return false;

                key = makeKey( 3, myRank, elt.id(), attachment.ordinal, layout.component, attachment.kind );
                return true;
            }
            default:
                return false;
            }
        };

        std::map<descriptor_key_type,size_type> localSharedKeyToDof;
        auto rangeElements = elements( mesh, entity_process_t::LOCAL_ONLY );
        for ( auto const& eltWrap : rangeElements )
        {
            auto const& elt = unwrap_ref( eltWrap );
            if ( !this->isElementDone( elt.id() ) )
                continue;

            for ( uint16_type parentLid = 0; parentLid < nLocalDof; ++parentLid )
            {
                descriptor_key_type key;
                bool isShared = false;
                if ( !makeElementDofKey( elt, parentLid, key, isShared ) || !isShared )
                    continue;

                auto const layout = fe.localDofLayout( parentLid );
                auto const& gdof = this->localToGlobal( elt.id(), layout.parentLocalDofId, layout.component );
                size_type const localDofIndex = gdof.index();
                CHECK( localDofIndex < nLocalDofWithGhost )
                    << fmt::format( "[dof-mpi-key] local dof index {} out of range {}", localDofIndex, nLocalDofWithGhost );

                auto [it, inserted] = localSharedKeyToDof.emplace( key, localDofIndex );
                CHECK( inserted || it->second == localDofIndex )
                    << "[dof-mpi-key] same descriptor key resolved to different local dofs";
            }
        }

        std::vector<size_type> localKeyPayload;
        localKeyPayload.reserve( localSharedKeyToDof.size()*7 );
        for ( auto const& [key,localDofIndex] : localSharedKeyToDof )
        {
            Feel::detail::ignore_unused_variable_warning( localDofIndex );
            packKey( key, localKeyPayload );
        }

        std::map<rank_type,std::vector<size_type>> keysToSend, keysToRecv;
        for ( rank_type neighborRank : mesh.neighborSubdomains() )
            keysToSend[neighborRank] = localKeyPayload;

        std::vector<mpi::request> reqs( 2*mesh.neighborSubdomains().size() );
        int countRequest = 0;
        for ( rank_type neighborRank : mesh.neighborSubdomains() )
        {
            reqs[countRequest++] = this->worldComm().localComm().irecv( neighborRank, 0, keysToRecv[neighborRank] );
            reqs[countRequest++] = this->worldComm().localComm().isend( neighborRank, 0, keysToSend[neighborRank] );
        }
        mpi::wait_all( std::begin( reqs ), std::begin( reqs ) + countRequest );

        std::map<descriptor_key_type,std::set<rank_type>> ranksByKey;
        for ( auto const& [key,localDofIndex] : localSharedKeyToDof )
        {
            Feel::detail::ignore_unused_variable_warning( localDofIndex );
            ranksByKey[key].insert( myRank );
        }

        for ( auto const& [rankRecv,payload] : keysToRecv )
        {
            CHECK( payload.size() % 7 == 0 )
                << fmt::format( "[dof-mpi-key] invalid key payload size {} from rank {}", payload.size(), rankRecv );
            for ( size_type k = 0; k < payload.size(); k += 7 )
            {
                auto key = unpackKey( payload, k );
                if ( localSharedKeyToDof.find( key ) != localSharedKeyToDof.end() )
                    ranksByKey[key].insert( rankRecv );
            }
        }

        std::map<descriptor_key_type,rank_type> ownerByKey;
        std::vector<bool> dofIsGhost( nLocalDofWithGhost, false );
        size_type nDofNotPresent = 0;
        for ( auto const& [key,ranks] : ranksByKey )
        {
            if ( ranks.size() <= 1 )
                continue;
            rank_type const ownerRank = *ranks.begin();
            ownerByKey.emplace( key, ownerRank );
            if ( ownerRank == myRank )
                continue;

            auto const itLocal = localSharedKeyToDof.find( key );
            CHECK( itLocal != localSharedKeyToDof.end() ) << "[dof-mpi-key] missing local dof for shared key";
            if ( !dofIsGhost[itLocal->second] )
            {
                dofIsGhost[itLocal->second] = true;
                ++nDofNotPresent;
            }
        }

        CHECK( this->M_n_localWithGhost_df[myRank] >= nDofNotPresent ) << "invalid descriptor-key ghost count";
        this->M_n_localWithoutGhost_df[myRank] = this->M_n_localWithGhost_df[myRank] - nDofNotPresent;

        std::vector<std::tuple<size_type,size_type>> dataRecvFromGather;
        auto dataSendToGather = std::make_tuple( this->M_n_localWithGhost_df[myRank],
                                                 this->M_n_localWithoutGhost_df[myRank] );
        mpi::all_gather( this->worldComm(), dataSendToGather, dataRecvFromGather );

        for ( rank_type p = 0; p < nProc; ++p )
        {
            this->M_n_localWithGhost_df[p] = std::get<0>( dataRecvFromGather[p] );
            this->M_n_localWithoutGhost_df[p] = std::get<1>( dataRecvFromGather[p] );
        }

        this->M_n_dofs = 0;
        for ( rank_type proc = 0; proc < nProc; ++proc )
            this->M_n_dofs += this->M_n_localWithoutGhost_df[proc];

        this->M_first_df_globalcluster[0] = 0;
        if ( this->M_n_localWithoutGhost_df[0] > 0 )
            this->M_last_df_globalcluster[0] = this->M_first_df_globalcluster[0] + this->M_n_localWithoutGhost_df[0] - 1;
        else
            this->M_last_df_globalcluster[0] = this->M_first_df_globalcluster[0];

        for ( rank_type p = 1; p < nProc; ++p )
        {
            if ( this->M_n_localWithoutGhost_df[p-1] > 0 )
                this->M_first_df_globalcluster[p] = this->M_last_df_globalcluster[p-1] + 1;
            else
                this->M_first_df_globalcluster[p] = this->M_last_df_globalcluster[p-1];

            if ( this->M_n_localWithoutGhost_df[p] > 0 )
                this->M_last_df_globalcluster[p] = this->M_first_df_globalcluster[p] + this->M_n_localWithoutGhost_df[p] - 1;
            else
                this->M_last_df_globalcluster[p] = this->M_first_df_globalcluster[p];
        }

        this->M_mapGlobalProcessToGlobalCluster.resize( this->M_n_localWithGhost_df[myRank], invalid_v<size_type> );
        size_type nextGlobIndex = this->M_first_df_globalcluster[myRank];
        for ( size_type i = 0; i < this->M_n_localWithGhost_df[myRank]; ++i )
        {
            if ( !dofIsGhost[i] )
                this->M_mapGlobalProcessToGlobalCluster[i] = nextGlobIndex++;
        }

        std::map<rank_type,std::vector<size_type>> ownerDataToSend, ownerDataToRecv;
        for ( auto const& [key,ranks] : ranksByKey )
        {
            auto const itOwner = ownerByKey.find( key );
            if ( itOwner == ownerByKey.end() || itOwner->second != myRank )
                continue;

            auto const itLocal = localSharedKeyToDof.find( key );
            CHECK( itLocal != localSharedKeyToDof.end() ) << "[dof-mpi-key] owner key without local dof";
            size_type const localDofIndex = itLocal->second;
            size_type const gcId = this->M_mapGlobalProcessToGlobalCluster[localDofIndex];
            CHECK( gcId != invalid_v<size_type> ) << "[dof-mpi-key] owner key without global-cluster id";

            for ( rank_type sharedRank : ranks )
            {
                if ( sharedRank == myRank )
                    continue;
                packKey( key, ownerDataToSend[sharedRank] );
                ownerDataToSend[sharedRank].push_back( gcId );
                this->M_activeDofSharedOnCluster[localDofIndex].insert( sharedRank );
                this->addNeighborSubdomain( sharedRank );
            }
        }

        reqs.assign( 2*mesh.neighborSubdomains().size(), mpi::request{} );
        countRequest = 0;
        for ( rank_type neighborRank : mesh.neighborSubdomains() )
        {
            reqs[countRequest++] = this->worldComm().localComm().irecv( neighborRank, 0, ownerDataToRecv[neighborRank] );
            reqs[countRequest++] = this->worldComm().localComm().isend( neighborRank, 0, ownerDataToSend[neighborRank] );
        }
        mpi::wait_all( std::begin( reqs ), std::begin( reqs ) + countRequest );

        for ( auto const& [rankRecv,payload] : ownerDataToRecv )
        {
            CHECK( payload.size() % 8 == 0 )
                << fmt::format( "[dof-mpi-key] invalid owner payload size {} from rank {}", payload.size(), rankRecv );
            for ( size_type k = 0; k < payload.size(); k += 8 )
            {
                auto key = unpackKey( payload, k );
                size_type const gcId = payload[k+7];
                auto const itLocal = localSharedKeyToDof.find( key );
                if ( itLocal == localSharedKeyToDof.end() )
                    continue;

                size_type const localDofIndex = itLocal->second;
                if ( this->M_mapGlobalProcessToGlobalCluster[localDofIndex] == invalid_v<size_type> )
                    this->M_mapGlobalProcessToGlobalCluster[localDofIndex] = gcId;
                else
                    CHECK( this->M_mapGlobalProcessToGlobalCluster[localDofIndex] == gcId )
                        << "[dof-mpi-key] inconsistent owner global-cluster id";

                rank_type const activeProcId = this->procOnGlobalCluster( gcId );
                if ( activeProcId != myRank )
                    this->addNeighborSubdomain( activeProcId );
            }
        }

        for ( auto const& [key,ownerRank] : ownerByKey )
        {
            auto const itLocal = localSharedKeyToDof.find( key );
            CHECK( itLocal != localSharedKeyToDof.end() ) << "[dof-mpi-key] missing local key after ownership exchange";
            size_type const localDofIndex = itLocal->second;
            if ( ownerRank != myRank )
                CHECK( this->M_mapGlobalProcessToGlobalCluster[localDofIndex] != invalid_v<size_type> )
                    << "[dof-mpi-key] unresolved ghost global-cluster id";
        }

        if ( this->hasDofTableExtended() )
        {
            if ( this->hasMeshSupport() && this->meshSupport()->isPartialSupport() )
                this->buildGhostDofMapExtended( mesh, elements( this->meshSupport(), entity_process_t::GHOST_ONLY ) );
            else
                this->buildGhostDofMapExtended( mesh, elements( mesh, entity_process_t::GHOST_ONLY ) );
        }

        VLOG(1) << "[dof-mpi-key] descriptor-key MPI ownership path"
                << " local_shared_keys=" << localSharedKeyToDof.size()
                << " shared_keys=" << ownerByKey.size()
                << " ghosts=" << nDofNotPresent;
        this->M_hasDescriptorKeyClusterDofMap = true;
        return true;
    }
}


template<typename MeshType, typename FEType, typename MortarType>
void
DofTable<MeshType, FEType, MortarType>::buildGlobalProcessToGlobalClusterDofMapOthersMesh( mesh_type& mesh )
{
    DVLOG(2) << "[buildGlobalProcessToGlobalClusterDofMapOthersMesh] start\n";

    this->M_hasDescriptorKeyClusterDofMap = false;
    if ( this->buildGlobalProcessToGlobalClusterDofMapDescriptorKeys( mesh ) )
        return;

    // Runtime FE layout (supports both compile-time and dynamic polynomial orders).
    const uint16_type nDofPerVertexRt = runtimeDofPerVertex();
    const uint16_type nDofPerEdgeRt = runtimeDofPerEdge();
    const uint16_type nDofPerFaceRt = runtimeDofPerFace();

    const uint16_type nLocalDofUpToVertices = element_type::numVertices * nDofPerVertexRt;
    const uint16_type nLocalDofUpToEdges = nLocalDofUpToVertices + element_type::numEdges * nDofPerEdgeRt;
    const uint16_type nLocalDofUpToFaces =
        (nDim==1)? nLocalDofUpToVertices :
        (nDim==2)? nLocalDofUpToVertices + element_type::numEdges * nDofPerEdgeRt :
                   nLocalDofUpToVertices + element_type::numEdges * nDofPerEdgeRt + element_type::numGeometricFaces * nDofPerFaceRt;

    const uint16_type nDofPerVertexForDivision = ( nDofPerVertexRt > 0 ) ? nDofPerVertexRt : uint16_type( 1 );
    const uint16_type nDofPerEdgeForDivision = ( nDofPerEdgeRt > 0 ) ? nDofPerEdgeRt : uint16_type( 1 );
    const uint16_type nLocalDofBeforeDofsOfTopologicalFaces =
        (nDim==1)? 0 :
        (nDim==2)? nLocalDofUpToVertices :
                   nLocalDofUpToVertices + element_type::numEdges * nDofPerEdgeRt;

    const uint16_type nDofPerTopologicalFace = (nDim==1)? nDofPerVertexRt : (nDim==2)? nDofPerEdgeRt : nDofPerFaceRt;
    const uint16_type nDofPerTopologicalFaceForDivision = ( nDofPerTopologicalFace > 0 ) ? nDofPerTopologicalFace : uint16_type( 1 );

    const rank_type myRank = this->worldComm().localRank();
    const rank_type nProc = this->worldComm().localSize();

    std::set<uint16_type> localDofUsedForTensor2symm;
    if ( is_tensor2symm )
    {
        std::map<uint16_type,std::vector<uint16_type> > symm2unsymm;
        for ( uint16_type k=0;k<this->nLocalDof();++k )
            symm2unsymm[this->fe().unsymmToSymm(k)].push_back( k );
        for ( auto const& symmdof : symm2unsymm )
        {
            if ( symmdof.second.empty() )
                continue;
            localDofUsedForTensor2symm.insert( symmdof.second.front() );
        }
    }


    // create inverse mapping of vector_permutation (vector_permutation map elt ordering to face ordering)
    // thus, this mapping will map face ordering to elt ordering with respect to permutation
    std::map<face_permutation_type, permutation_vector_type> mapLocalDofFaceToFaceInElement;
    if constexpr ( nDim == 3 )
    {
        if ( nDofPerTopologicalFace > 1 )
        {
            for ( auto const& [perm,dofsMapping] : this->vector_permutation )
            {
                mapLocalDofFaceToFaceInElement[perm].resize(dofsMapping.size());
                for (int k=0;k<dofsMapping.size();++k)
                    mapLocalDofFaceToFaceInElement[perm][dofsMapping[k]] = k;
            }
        }
    }


    size_type nLocalDofWithGhost = this->M_n_localWithGhost_df[myRank];
    std::vector<bool> dofdone( nLocalDofWithGhost,false);
    std::vector<bool> dofIsGhost( nLocalDofWithGhost,false);
    size_type nDofNotPresent=0;


    // maybe a vector instead of of map??
    std::map<rank_type, std::map<size_type,std::vector<uint16_type> > > dataToSend;

    //typename MeshTraits<mesh_type>::elements_reference_wrapper_ptrtype myActiveEltsTouchInterProcess( new typename MeshTraits<mesh_type>::elements_reference_wrapper_type );

    bool hasMeshSupportPartial = this->hasMeshSupport() && this->meshSupport()->isPartialSupport();

    std::map<rank_type, std::map<size_type,std::vector<size_type> > > dataMemory;
    using periodic_dof_key_type = std::tuple<uint16_type,size_type,uint16_type>; // (entity dim, global dof key, component)
    std::map<periodic_dof_key_type,size_type> periodicKeyToLocalDof;
    std::map<periodic_dof_key_type,rank_type> periodicOwnerByKey;
    std::map<periodic_dof_key_type,std::set<rank_type>> periodicRanksByKey;
    std::map<periodic_dof_key_type,size_type> periodicMultiplicity;
    bool hasPeriodicSharedDofs = false;
    if ( is_continuous )
    {
        std::map<size_type, std::tuple<rank_type,size_type,uint16_type>> mapPtIdToDofOwnerElt; //( ptId ->( rank, eltId, ptIdInElt  ))
        std::map<size_type, std::tuple<rank_type,size_type,uint16_type,edge_permutation_type>> mapEdgeIdToDofOwnerElt; //( ptId ->( rank, eltId, edgeIdInElt, edge permutation ))
        std::map<size_type, std::tuple<rank_type,size_type,uint16_type,face_permutation_type>> mapFaceIdToDofOwnerElt; //( faceId ->( rank, eltId, faceIdInElt, face permutation))
        std::vector< std::reference_wrapper<const typename mesh_type::element_type> > activeEltTouchInterprocess;
        auto rangeElements = hasMeshSupportPartial? elements( this->meshSupport(), entity_process_t::ALL ) : elements( mesh, entity_process_t::ALL );
        for ( auto const& eltWrap : rangeElements )
        {
            auto const& elt = unwrap_ref( eltWrap );
            rank_type eltPid = elt.processId();
            size_type eltIdOtherPart = elt.idInOthersPartitions( eltPid );
            bool currentEltTouchInterprocess = false;
            for ( uint16_type n=0; n < elt.nVertices(); n++ )
            {
                auto const& point = elt.point(n);
                bool isInterprocessPoint = hasMeshSupportPartial? this->meshSupport()->isInterprocessPoints( point.id() ) : mesh.isInterprocessPoints( point.id() );
                if ( isInterprocessPoint )
                {
                    currentEltTouchInterprocess = true;
                    auto itFind = mapPtIdToDofOwnerElt.find( point.id() );
                    if ( itFind == mapPtIdToDofOwnerElt.end() )
                        mapPtIdToDofOwnerElt.emplace( point.id(), std::make_tuple( eltPid, eltIdOtherPart, n ) );
                    else
                    {
                        auto & currentData = itFind->second;
                        if ( elt.processId() < std::get<0>( currentData ) )
                            currentData = std::make_tuple( eltPid, eltIdOtherPart, n );
                    }
                }
            }
            if constexpr ( nDim == 3 )
            {
                if ( nDofPerEdgeRt > 0 && currentEltTouchInterprocess )
                {
                    for ( size_type j = 0; j < elt.nEdges(); j++ )
                    {
                        if ( !elt.edgePtr( j ) )
                            continue;
                        auto & edge = elt.edge( j );

                        bool isInterprocessEdge = hasMeshSupportPartial? this->meshSupport()->isInterprocessEdges( edge.id() ) : mesh.isInterprocessEdges( edge.id() );
                        if ( isInterprocessEdge )
                        {
                            auto itFind = mapEdgeIdToDofOwnerElt.find( edge.id() );
                            if ( itFind == mapEdgeIdToDofOwnerElt.end() )
                                mapEdgeIdToDofOwnerElt.emplace( edge.id(), std::make_tuple( eltPid, eltIdOtherPart, j, elt.edgePermutation( j ).value() ) );
                            else
                            {
                                auto & currentData = itFind->second;
                                if ( elt.processId() < std::get<0>( currentData ) )
                                    currentData = std::make_tuple( eltPid, eltIdOtherPart, j, elt.edgePermutation( j ).value() );
                            }
                        }
                    }
                }
            }

            if ( nDofPerTopologicalFace > 0 && currentEltTouchInterprocess )
            {
                // face
                for ( uint16_type j = 0; j < elt.nTopologicalFaces(); j++ )
                {
                    if ( !elt.facePtr( j ) )
                        continue;
                    auto const& face = elt.face( j );

                    // TODO: maybe ignore face that don't touch interprocess??

                    auto itFind = mapFaceIdToDofOwnerElt.find( face.id() );
                    if ( itFind == mapFaceIdToDofOwnerElt.end() )
                        mapFaceIdToDofOwnerElt.emplace( face.id(), std::make_tuple( eltPid, eltIdOtherPart, j, elt.facePermutation( j ).value() ) );
                    else
                    {
                        auto & currentData = itFind->second;
                        if ( elt.processId() < std::get<0>( currentData ) )
                            currentData = std::make_tuple( eltPid, eltIdOtherPart, j, elt.facePermutation( j ).value() );
                    }
                }
            }

            if ( !elt.isGhostCell() && currentEltTouchInterprocess )
                activeEltTouchInterprocess.push_back( std::cref( elt ) );
        }

        for ( auto const& eltWrap : activeEltTouchInterprocess )
        {
            auto const& elt = eltWrap.get();
            // prepare local dofs to analyse
            std::map<uint16_type,std::map<uint16_type,std::tuple<uint16_type,size_type> > > mapLocalDofByCompToLocalDofAllComp;
            for ( auto const& ldof : this->localDof( elt.id() ) )
            {
                uint16_type ldofId = ldof.first.localDof();
                if ( is_tensor2symm )
                    if ( localDofUsedForTensor2symm.find( ldofId ) == localDofUsedForTensor2symm.end() )
                        continue;
                uint16_type ldofParentId = this->fe().dofParent( ldofId );
                size_type gdofId = ldof.second.index();
                uint16_type comp = this->fe().component(ldofId);
                mapLocalDofByCompToLocalDofAllComp[ldofParentId][comp] = std::make_tuple( ldofId,gdofId );
            }
            // loop over local dof for the detection of ghost dofs
            for (auto const& localDofDatas : mapLocalDofByCompToLocalDofAllComp )
            {
                CHECK( !localDofDatas.second.empty() ) << "no localdof data is empty";
                auto const& localDofDataFirstComponent = *localDofDatas.second.begin();
                const uint16_type locDof = localDofDatas.first;
                //const uint16_type locDof = std::get<0>( localDofDataFirstComponent.second );
                const size_type theglobdoftest = std::get<1>( localDofDataFirstComponent.second );
                DCHECK( theglobdoftest < nLocalDofWithGhost ) << "invalid globdof " << theglobdoftest << "\n";

                if ( dofdone[theglobdoftest] )
                    continue;
                dofdone[theglobdoftest]=true;

                rank_type eltPidOwnerDof = invalid_v<rank_type>;
                size_type eltIdOwnerDof = invalid_v<size_type>;
                uint16_type localDofOwnerDof = invalid_v<uint16_type>;
                if ( locDof < nLocalDofUpToVertices )
                {
                    uint16_type pointIdInElt = locDof / nDofPerVertexForDivision;
                    auto const& point = elt.point( pointIdInElt );
                    auto itFind = mapPtIdToDofOwnerElt.find( point.id() );
                    if ( itFind == mapPtIdToDofOwnerElt.end() )
                        continue;
                    std::tie( eltPidOwnerDof, eltIdOwnerDof, localDofOwnerDof ) = itFind->second;
                }
                else if ( nDim == 3 && locDof < nLocalDofUpToEdges )
                {
                    if constexpr (nDim == 3)
                    {
                        uint16_type locDofInEgdes = locDof - nLocalDofUpToVertices;
                        uint16_type edgeIdInElt = locDofInEgdes / nDofPerEdgeForDivision;
                        uint16_type locDofInEgde = locDofInEgdes % nDofPerEdgeForDivision;

                        auto edgePtr = elt.edgePtr( edgeIdInElt );
                        if ( !edgePtr )
                            continue;
                        auto const& edge = *edgePtr;
                        auto itFind = mapEdgeIdToDofOwnerElt.find( edge.id() );
                        if ( itFind == mapEdgeIdToDofOwnerElt.end() )
                            continue;

                        edge_permutation_type edgePermutationOwnerDof;
                        uint16_type edgeIdInEltOwnerDof = invalid_v<uint16_type>;
                        std::tie( eltPidOwnerDof, eltIdOwnerDof, edgeIdInEltOwnerDof, edgePermutationOwnerDof ) = itFind->second;
                        // do nothing if owner is on curent process id
                        if ( eltPidOwnerDof == myRank )
                            continue;

                        auto edgePermutation = elt.edgePermutation( edgeIdInElt );
                        // shift local dof to relative dof edge numbering
                        localDofOwnerDof = nLocalDofUpToVertices + edgeIdInEltOwnerDof * nDofPerEdgeRt;
                        localDofOwnerDof += ( edgePermutation == edgePermutationOwnerDof )? locDofInEgde : nDofPerEdgeRt - 1 - locDofInEgde;
                    }
                }
                else if ( locDof < nLocalDofUpToFaces )
                {
                    if ( nDim < 2 && nDim != nRealDim )
                        continue;
                    uint16_type locDofInTopologicalFaces = locDof - nLocalDofBeforeDofsOfTopologicalFaces;
                    uint16_type faceIdInElt = locDofInTopologicalFaces / nDofPerTopologicalFaceForDivision;
                    uint16_type locDofInTopologicalFace = locDofInTopologicalFaces % nDofPerTopologicalFaceForDivision;

                    auto facePtr = elt.facePtr( faceIdInElt );
                    if ( !facePtr )
                        continue;
                    auto const& face = *facePtr;

                    auto itFind = mapFaceIdToDofOwnerElt.find( face.id() );
                    if ( itFind == mapFaceIdToDofOwnerElt.end() )
                        continue;

                    face_permutation_type facePermutationOwnerDof;
                    uint16_type faceIdInEltOwnerDof = invalid_v<uint16_type>;
                    std::tie( eltPidOwnerDof, eltIdOwnerDof, faceIdInEltOwnerDof, facePermutationOwnerDof ) = itFind->second;

                    // do nothing if owner is on curent process id
                    if ( eltPidOwnerDof == myRank )
                        continue;

                    // shift local dof to relative dof face numbering
                    localDofOwnerDof = nLocalDofBeforeDofsOfTopologicalFaces + faceIdInEltOwnerDof*nDofPerTopologicalFace;

                    if ( nDofPerTopologicalFace > 1 )
                    {
                        auto facePermutation = elt.facePermutation( faceIdInElt );
                        if constexpr ( nDim <= 2)
                        {
                            localDofOwnerDof += ( facePermutation == facePermutationOwnerDof )? locDofInTopologicalFace : nDofPerTopologicalFace-1-locDofInTopologicalFace;
                        }
                        else
                        {
                            if ( facePermutation == facePermutationOwnerDof )
                                localDofOwnerDof += locDofInTopologicalFace;
                            else
                            {
                                uint16_type locFaceDof = invalid_v<uint16_type>;
                                if ( facePermutation.value() == face_permutation_type::IDENTITY )
                                    locFaceDof = locDofInTopologicalFace;
                                else
                                    locFaceDof = this->vector_permutation.at(facePermutation)[locDofInTopologicalFace];

                                if ( facePermutationOwnerDof.value() == face_permutation_type::IDENTITY )
                                    localDofOwnerDof += locFaceDof;
                                else
                                    localDofOwnerDof += mapLocalDofFaceToFaceInElement.at(facePermutationOwnerDof)[locFaceDof];
                            }
                        }

                        // TODO : add debug check that verify dof point are identical!
                    }
                }

                if ( eltPidOwnerDof == invalid_v<rank_type> )
                    continue;

                if ( eltPidOwnerDof != myRank )// invalid_rank_type_value && pidDofActive < myRank )
                {
                    std::vector<std::pair<uint16_type,size_type> > compglobdofs;
                    for ( auto const& localDofAllComp : localDofDatas.second )
                    {
                        uint16_type comp = localDofAllComp.first;
                        const size_type theglobdof = std::get<1>( localDofAllComp.second );
                        const uint16_type thelocdof = std::get<0>( localDofAllComp.second );
                        dofIsGhost[theglobdof] = true;
                        //compglobdofs.push_back( std::make_pair(thelocdof,theglobdof) );
                        ++nDofNotPresent;
#if 1
                        dataToSend[eltPidOwnerDof][eltIdOwnerDof].push_back( this->localDofId( localDofOwnerDof, comp ) );
                        dataMemory[eltPidOwnerDof][eltIdOwnerDof].push_back( theglobdof );
#endif
                    }
                }

            }

        }
    }

    // Periodic MPI ownership synchronization:
    // detect dofs duplicated across ranks by FE key and force a single owner (lowest rank).
    if ( is_continuous && mesh.isPeriodic() && nProc > 1 )
    {
        this->forEachGlobalDofKeyEntry( [&]( auto const& entry )
        {
            if ( entry.localDofIndex >= nLocalDofWithGhost )
                return;
            periodicKeyToLocalDof.emplace( std::make_tuple( std::get<0>( entry.key ),
                                                            std::get<1>( entry.key ),
                                                            entry.component ),
                                           entry.localDofIndex );
        } );

        std::vector<size_type> periodicLocalKeysFlat;
        periodicLocalKeysFlat.reserve( periodicKeyToLocalDof.size()*3 );
        for ( auto const& [key,localDofIndex] : periodicKeyToLocalDof )
        {
            Feel::detail::ignore_unused_variable_warning( localDofIndex );
            periodicLocalKeysFlat.push_back( static_cast<size_type>( std::get<0>( key ) ) );
            periodicLocalKeysFlat.push_back( std::get<1>( key ) );
            periodicLocalKeysFlat.push_back( static_cast<size_type>( std::get<2>( key ) ) );
        }

        std::vector<std::vector<size_type>> periodicKeysFlatPerRank;
        mpi::all_gather( this->worldComm().localComm(), periodicLocalKeysFlat, periodicKeysFlatPerRank );

        CHECK( periodicKeysFlatPerRank.size() == static_cast<size_type>( nProc ) )
            << fmt::format( "invalid periodic all_gather size {} vs {}", periodicKeysFlatPerRank.size(), nProc );

        for ( rank_type rank = 0; rank < nProc; ++rank )
        {
            auto const& flat = periodicKeysFlatPerRank[rank];
            CHECK( flat.size() % 3 == 0 )
                << fmt::format( "invalid periodic key payload size {} on rank {}", flat.size(), rank );
            for ( size_type k = 0; k < flat.size(); k += 3 )
            {
                periodic_dof_key_type key{
                    static_cast<uint16_type>( flat[k] ),
                    flat[k+1],
                    static_cast<uint16_type>( flat[k+2] )
                };

                auto [itOwner, insertedOwner] = periodicOwnerByKey.emplace( key, rank );
                if ( !insertedOwner && rank < itOwner->second )
                    itOwner->second = rank;
                ++periodicMultiplicity[key];
                periodicRanksByKey[key].insert( rank );
            }
        }

        size_type periodicAdditionalGhostDofs = 0;
        for ( auto const& [key,ownerRank] : periodicOwnerByKey )
        {
            auto const itMultiplicity = periodicMultiplicity.find( key );
            if ( itMultiplicity == periodicMultiplicity.end() || itMultiplicity->second <= 1 )
                continue;

            hasPeriodicSharedDofs = true;
            auto const itLocal = periodicKeyToLocalDof.find( key );
            if ( itLocal == periodicKeyToLocalDof.end() )
                continue;

            size_type const localDofIndex = itLocal->second;
            if ( ownerRank != myRank && !dofIsGhost[localDofIndex] )
            {
                dofIsGhost[localDofIndex] = true;
                ++nDofNotPresent;
                ++periodicAdditionalGhostDofs;
            }
        }

        VLOG(1) << "[periodic][mpi] periodic shared dof ownership resolution"
                << " local_keys=" << periodicKeyToLocalDof.size()
                << " additional_ghosts=" << periodicAdditionalGhostDofs;
    }

    //------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------//

    // update datamap info
    CHECK( this->M_n_localWithGhost_df[myRank] >= nDofNotPresent ) << "invalid data";
    this->M_n_localWithoutGhost_df[myRank] = this->M_n_localWithGhost_df[myRank] - nDofNotPresent;

    // std::vector<boost::tuple<size_type,size_type,size_type> > dataRecvFromGather;
    // auto dataSendToGather = boost::make_tuple(this->M_first_df[myRank],this->M_n_localWithGhost_df[myRank],this->M_n_localWithoutGhost_df[myRank]);
    std::vector<std::tuple<size_type,size_type> > dataRecvFromGather;
    auto dataSendToGather = std::make_tuple(this->M_n_localWithGhost_df[myRank],this->M_n_localWithoutGhost_df[myRank]);
    mpi::all_gather( this->worldComm(),
                     dataSendToGather,
                     dataRecvFromGather );

    for (int p=0;p<this->worldComm().localSize();++p)
    {
        this->M_n_localWithGhost_df[p] = std::get<0>( dataRecvFromGather[p] );
        this->M_n_localWithoutGhost_df[p] = std::get<1>( dataRecvFromGather[p] );
    }
    // update global nDof
    this->M_n_dofs=0;
    for ( int proc=0; proc<this->worldComm().size(); ++proc )
    {
        this->M_n_dofs+=this->M_n_localWithoutGhost_df[proc];
    }

    this->M_first_df_globalcluster[0]=0;//this->M_first_df[0];
    if ( this->M_n_localWithoutGhost_df[0] > 0 )
        this->M_last_df_globalcluster[0] = this->M_first_df_globalcluster[0]+this->M_n_localWithoutGhost_df[0]-1;
    else
        this->M_last_df_globalcluster[0] = this->M_first_df_globalcluster[0];

    for ( int i=1; i<this->worldComm().size(); ++i )
    {
        if ( this->M_n_localWithoutGhost_df[i-1] >0 )
            this->M_first_df_globalcluster[i]=this->M_last_df_globalcluster[i-1]+1;
        else
            this->M_first_df_globalcluster[i]=this->M_last_df_globalcluster[i-1];

        if ( this->M_n_localWithoutGhost_df[i] >0 )
            this->M_last_df_globalcluster[i]=this->M_first_df_globalcluster[i]+this->M_n_localWithoutGhost_df[i]-1;
        else
            this->M_last_df_globalcluster[i]=this->M_first_df_globalcluster[i];
    }
    //------------------------------------------------------------------------------//
    // init map
    this->M_mapGlobalProcessToGlobalCluster.resize( this->M_n_localWithGhost_df[myRank],invalid_v<size_type> );
    //------------------------------------------------------------------------------//
    // add in map the dofs presents
    size_type firstGlobIndex = this->M_first_df_globalcluster[myRank];
    size_type nextGlobIndex = firstGlobIndex;
    for ( size_type i=0; i< this->M_n_localWithGhost_df[myRank]; ++i )
    {
        if ( !dofIsGhost[i] )
        {
            this->M_mapGlobalProcessToGlobalCluster[i]=nextGlobIndex;
            ++nextGlobIndex;
        }
    }
   //------------------------------------------------------------------------------//

    // update parallel mapping
#if 0
    this->buildGlobalProcessToGlobalClusterDofMapOthersMeshNonBlockingComm( mesh,listToSend );
#else
    this->buildGlobalProcessToGlobalClusterInterprocessDofs( mesh, dataToSend, dataMemory );
#endif

    if ( is_continuous && mesh.isPeriodic() && nProc > 1 && hasPeriodicSharedDofs )
    {
        // Broadcast owner-side global-cluster ids for duplicated periodic keys.
        std::vector<size_type> periodicOwnedClustersFlat;
        periodicOwnedClustersFlat.reserve( periodicOwnerByKey.size()*4 );
        for ( auto const& [key,ownerRank] : periodicOwnerByKey )
        {
            auto const itMultiplicity = periodicMultiplicity.find( key );
            if ( itMultiplicity == periodicMultiplicity.end() || itMultiplicity->second <= 1 )
                continue;
            if ( ownerRank != myRank )
                continue;

            auto const itLocal = periodicKeyToLocalDof.find( key );
            if ( itLocal == periodicKeyToLocalDof.end() )
                continue;

            size_type const localDofIndex = itLocal->second;
            CHECK( localDofIndex < this->M_mapGlobalProcessToGlobalCluster.size() )
                << fmt::format( "invalid periodic local dof index {} vs {}", localDofIndex, this->M_mapGlobalProcessToGlobalCluster.size() );
            size_type const localGcDof = this->M_mapGlobalProcessToGlobalCluster[localDofIndex];
            CHECK( localGcDof != invalid_v<size_type> )
                << fmt::format( "invalid global cluster id for periodic owner dof {}", localDofIndex );

            periodicOwnedClustersFlat.push_back( static_cast<size_type>( std::get<0>( key ) ) );
            periodicOwnedClustersFlat.push_back( std::get<1>( key ) );
            periodicOwnedClustersFlat.push_back( static_cast<size_type>( std::get<2>( key ) ) );
            periodicOwnedClustersFlat.push_back( localGcDof );
        }

        std::vector<std::vector<size_type>> periodicOwnedClustersFlatPerRank;
        mpi::all_gather( this->worldComm().localComm(), periodicOwnedClustersFlat, periodicOwnedClustersFlatPerRank );
        CHECK( periodicOwnedClustersFlatPerRank.size() == static_cast<size_type>( nProc ) )
            << fmt::format( "invalid periodic owner all_gather size {} vs {}", periodicOwnedClustersFlatPerRank.size(), nProc );

        std::map<periodic_dof_key_type,size_type> periodicGcByKey;
        for ( rank_type rank = 0; rank < nProc; ++rank )
        {
            auto const& flat = periodicOwnedClustersFlatPerRank[rank];
            CHECK( flat.size() % 4 == 0 )
                << fmt::format( "invalid periodic cluster payload size {} on rank {}", flat.size(), rank );
            for ( size_type k = 0; k < flat.size(); k += 4 )
            {
                periodic_dof_key_type key{
                    static_cast<uint16_type>( flat[k] ),
                    flat[k+1],
                    static_cast<uint16_type>( flat[k+2] )
                };
                size_type const gcId = flat[k+3];

                auto [itGc, insertedGc] = periodicGcByKey.emplace( key, gcId );
                if ( !insertedGc )
                    CHECK( itGc->second == gcId )
                        << fmt::format( "inconsistent periodic owner global-cluster id {} vs {}", itGc->second, gcId );
            }
        }

        size_type periodicSyncedGhostDofs = 0;
        size_type periodicSharedOwnerDofs = 0;
        for ( auto const& [key,ownerRank] : periodicOwnerByKey )
        {
            auto const itMultiplicity = periodicMultiplicity.find( key );
            if ( itMultiplicity == periodicMultiplicity.end() || itMultiplicity->second <= 1 )
                continue;

            auto const itLocal = periodicKeyToLocalDof.find( key );
            if ( itLocal == periodicKeyToLocalDof.end() )
                continue;

            size_type const localDofIndex = itLocal->second;
            if ( ownerRank != myRank )
            {
                auto const itGc = periodicGcByKey.find( key );
                CHECK( itGc != periodicGcByKey.end() )
                    << "[periodic][mpi] missing owner global-cluster id for periodic key";
                size_type const gcId = itGc->second;
                CHECK( gcId != invalid_v<size_type> ) << "[periodic][mpi] invalid periodic global-cluster id";
                CHECK( localDofIndex < this->M_mapGlobalProcessToGlobalCluster.size() )
                    << fmt::format( "invalid periodic local dof index {} vs {}", localDofIndex, this->M_mapGlobalProcessToGlobalCluster.size() );
                if ( this->M_mapGlobalProcessToGlobalCluster[localDofIndex] == invalid_v<size_type> )
                    this->M_mapGlobalProcessToGlobalCluster[localDofIndex] = gcId;
                else
                    CHECK( this->M_mapGlobalProcessToGlobalCluster[localDofIndex] == gcId )
                        << fmt::format( "inconsistent periodic global-cluster id {} vs {}",
                                        this->M_mapGlobalProcessToGlobalCluster[localDofIndex], gcId );

                rank_type const activeProcId = this->procOnGlobalCluster( gcId );
                if ( activeProcId != myRank )
                    this->addNeighborSubdomain( activeProcId );
                ++periodicSyncedGhostDofs;
            }
            else
            {
                auto itRanks = periodicRanksByKey.find( key );
                if ( itRanks == periodicRanksByKey.end() )
                    continue;
                for ( rank_type sharedRank : itRanks->second )
                {
                    if ( sharedRank == myRank )
                        continue;
                    this->M_activeDofSharedOnCluster[localDofIndex].insert( sharedRank );
                    this->addNeighborSubdomain( sharedRank );
                    ++periodicSharedOwnerDofs;
                }
            }
        }

        VLOG(1) << "[periodic][mpi] periodic global-cluster synchronization"
                << " synced_ghost_dofs=" << periodicSyncedGhostDofs
                << " owner_shared_dofs=" << periodicSharedOwnerDofs;
    }

   //------------------------------------------------------------------------------//
   //------------------------------------------------------------------------------//
   //------------------------------------------------------------------------------//

    // extended dof table
    if ( this->hasDofTableExtended() )
    {
        if ( this->hasMeshSupport() && this->meshSupport()->isPartialSupport() )
            this->buildGhostDofMapExtended( mesh, elements(this->meshSupport(),entity_process_t::GHOST_ONLY ) );
        else
            this->buildGhostDofMapExtended( mesh, elements(mesh,entity_process_t::GHOST_ONLY ) );
    }
}



template<typename MeshType, typename FEType, typename MortarType>
void
DofTable<MeshType, FEType, MortarType>::buildGlobalProcessToGlobalClusterInterprocessDofs( mesh_type& mesh,
                                                                                                           std::map<rank_type, std::map<size_type,std::vector<uint16_type> > > & dataToSend,
                                                                                                           std::map<rank_type, std::map<size_type,std::vector<size_type> > > & dataMemory )
{
    std::map<rank_type, std::map<size_type,std::vector<uint16_type> > > dataToRecv;


    int nbMaxRequest = 2*mesh.neighborSubdomains().size();
    std::vector<mpi::request> reqs( nbMaxRequest );
    int countRequest = 0;

    // step 1 :send/recv of data
    for ( rank_type neighborRank : mesh.neighborSubdomains() )
    {
        reqs[countRequest++] = this->worldComm().localComm().irecv( neighborRank , 0, dataToRecv[neighborRank] );
        reqs[countRequest++] = this->worldComm().localComm().isend( neighborRank , 0, dataToSend[neighborRank] );
    }
    // step 1 :wait all requests
    mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
    countRequest = 0;

    //------------------------------------------------------------------------------//
    // step 2 : treat recv and prepare data for response
    std::vector<size_type> tmpIndicesSet;
    std::map<rank_type, std::map<size_type,std::vector<size_type> > > dataToSendStep2, dataToRecvStep2;
    for ( auto const& [rankRecv,eltIdToLocalDofs] : dataToRecv )
    {
        auto & dataToSendStep2AtRank = dataToSendStep2[rankRecv];
        for ( auto const& [eltId,localDofIds] : eltIdToLocalDofs )
        {
            tmpIndicesSet.resize( this->getIndicesSize( eltId ) );
            this->getIndicesSet( eltId, tmpIndicesSet );

            //this->getIndicesSetOnGlobalCluster( eltId, tmpIndicesSetOnGlobalCluster );
            auto & dataToSendStep2AtElt = dataToSendStep2AtRank[eltId];
            dataToSendStep2AtElt.resize( localDofIds.size(), invalid_v<size_type> );
            //for ( size_type localDofIds : localDofIds )
            for (int k=0;k<localDofIds.size();++k)
            {
                size_type dofIdGlobalProcess = tmpIndicesSet.at( localDofIds[k] );
                dataToSendStep2AtElt[k] = this->mapGlobalProcessToGlobalCluster()[ dofIdGlobalProcess ];
#if 0
                auto ababab = this->localDof( eltId );
                std::cout << fmt::format( "SEND M_mapGlobalProcessToGlobalCluster {} = {} with ldof {} and look size:{} eltId:{}",
                                          dofIdGlobalProcess,  dataToSendStep2AtElt[k], localDofIds[k], std::distance(ababab.first,ababab.second),eltId ) << std::endl;
#endif
                // this->M_activeDofSharedOnCluster[dofIdGlobalProcess].insert(rankRecv);
            }
            // dataToSendStep2AtRank.push_back( this->getIndicesOnGlobalCluster( eltId ) );
            // for ( size_type dofIndex : this->getIndices( eltId ) )
            //     if ( !this->dofGlobalProcessIsGhost( dofIndex ) )
            //         this->M_activeDofSharedOnCluster[dofIndex].insert(rankRecv);
        }
    }
    // step 2 :send/recv of data
    for ( rank_type neighborRank : mesh.neighborSubdomains() )
    {
        reqs[countRequest++] = this->worldComm().localComm().irecv( neighborRank , 0, dataToRecvStep2[neighborRank] );
        reqs[countRequest++] = this->worldComm().localComm().isend( neighborRank , 0, dataToSendStep2[neighborRank] );
    }
    // step 2 :wait all requests
    mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
    countRequest = 0;

    //------------------------------------------------------------------------------//
    // step 3 : treat final recv
    for ( auto const& [rankRecv,eltIdToDofGlobalClusters] : dataToRecvStep2 )
    {
        auto & dataMemoryAtRank = dataMemory[rankRecv];
        CHECK( dataMemoryAtRank.size() == eltIdToDofGlobalClusters.size() ) << "sizes should be equal";
        for ( auto const& [eltId,dofGlobalClusterIds] : eltIdToDofGlobalClusters )
        {
            auto & dataMemoryAtElt = dataMemoryAtRank.at(eltId);
            CHECK( dataMemoryAtElt.size() == dofGlobalClusterIds.size() ) << "sizes should be equal";
            for (int k=0;k<dataMemoryAtElt.size();++k)
            {
                size_type dofGlobalProcessId = dataMemoryAtElt[k];
                size_type dofGlobalClusterId = dofGlobalClusterIds[k];
                CHECK( dofGlobalClusterId != invalid_v<size_type> ) << "invalid global cluster dof id";
                //update data map
                //std::cout << fmt::format( "M_mapGlobalProcessToGlobalCluster {} = {}", dofGlobalProcessId, dofGlobalClusterId ) << std::endl;
                this->M_mapGlobalProcessToGlobalCluster[dofGlobalProcessId] = dofGlobalClusterId;
            }
        }
    }
}





//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//

template<typename MeshType, typename FEType, typename MortarType>
void
DofTable<MeshType, FEType, MortarType>::buildGhostDofMapExtended( mesh_type& mesh,
                                                                                   Range<mesh_type,MESH_ELEMENTS> const& ghostEltRange )
{
    DVLOG(2) << "[buildGhostDofMap] call buildGhostDofMapExtended on rank "<<  this->worldComm().rank();

    const rank_type myRank = this->worldComm().localRank();
    const rank_type nProc = this->worldComm().localSize();

    size_type start_next_free_dof = this->M_n_localWithGhost_df[myRank];
    //------------------------------------------------------------------------------//
    // build extended dof table
    size_type next_free_dof = start_next_free_dof;
    DofFromElement<self_type,fe_type> dfe( this, *M_fe );
    for ( auto const& ghostEltWrap : ghostEltRange )
    {
        auto const& ghostElt = boost::unwrap_ref( ghostEltWrap );
        // elements doftable
        if ( !this->isElementDone( ghostElt.id() ) )
            dfe.add( ghostEltWrap, next_free_dof, myRank );
    }
    //------------------------------------------------------------------------------//
    // update local datamap
    this->M_nGhostDofAddedInExtendedDofTable = next_free_dof-start_next_free_dof;
    std::vector<size_type> dataRecvFromGather;
    mpi::all_gather( this->worldComm().localComm(),
                     this->M_nGhostDofAddedInExtendedDofTable,
                     dataRecvFromGather );
    for (rank_type p=0;p<nProc;++p)
    {
        this->M_n_localWithGhost_df[p] += dataRecvFromGather[p];
    }
    this->M_mapGlobalProcessToGlobalCluster.resize( this->M_n_localWithGhost_df[myRank],invalid_v<size_type> );

    if constexpr ( FiniteElementDofLayoutProvider<fe_type> && is_continuous &&
                   !fe_type::is_modal && !is_tensor2symm && !is_product )
    {
        if ( !mesh.isPeriodic() && !( this->hasMeshSupport() && this->meshSupport()->isPartialSupport() ) )
        {
            using descriptor_key_type = DofKey<size_type>;

            auto familyTag = []() constexpr -> uint16_type
            {
                if constexpr ( requires { fe_type::TAG; } )
                    return fe_type::TAG;
                else
                    return 0;
            };

            auto makeKey = [&familyTag]( uint16_type topologicalDim,
                                         rank_type canonicalPartition,
                                         size_type canonicalEntityId,
                                         uint16_type ordinal,
                                         uint16_type component,
                                         uint16_type functional )
            {
                return descriptor_key_type{
                    .topologicalDim = static_cast<uint8_type>( topologicalDim ),
                    .canonicalPartition = canonicalPartition,
                    .canonicalEntityId = canonicalEntityId,
                    .ordinal = ordinal,
                    .component = component,
                    .familyTag = familyTag(),
                    .variant = functional
                };
            };

            auto packKey = []( descriptor_key_type const& key, std::vector<size_type>& payload )
            {
                payload.push_back( static_cast<size_type>( key.topologicalDim ) );
                payload.push_back( static_cast<size_type>( key.canonicalPartition ) );
                payload.push_back( key.canonicalEntityId );
                payload.push_back( static_cast<size_type>( key.ordinal ) );
                payload.push_back( static_cast<size_type>( key.component ) );
                payload.push_back( static_cast<size_type>( key.familyTag ) );
                payload.push_back( static_cast<size_type>( key.variant ) );
            };

            auto unpackKey = []( std::vector<size_type> const& payload, size_type offset )
            {
                return descriptor_key_type{
                    .topologicalDim = static_cast<uint8_type>( payload[offset] ),
                    .canonicalPartition = static_cast<rank_type>( payload[offset+1] ),
                    .canonicalEntityId = payload[offset+2],
                    .ordinal = static_cast<uint16_type>( payload[offset+3] ),
                    .component = static_cast<uint16_type>( payload[offset+4] ),
                    .familyTag = static_cast<uint16_type>( payload[offset+5] ),
                    .variant = static_cast<uint16_type>( payload[offset+6] )
                };
            };

            auto elementClusterId = []( element_type const& elt )
            {
                if ( elt.isGhostCell() )
                {
                    auto const ownerRank = elt.processId();
                    auto const idInOwner = elt.idInOthersPartitions( ownerRank );
                    if ( idInOwner != invalid_v<size_type> )
                        return idInOwner;
                }
                return elt.id();
            };

            auto elementClusterPartition = []( element_type const& elt )
            {
                return elt.isGhostCell() ? elt.processId() : elt.pidInPartition();
            };

            auto partitionCanonicalEntity = []( auto const& entity )
            {
                rank_type canonicalPartition = entity.pidInPartition();
                if ( canonicalPartition == invalid_rank_type_value )
                    canonicalPartition = entity.processId();
                size_type canonicalId = entity.id();
                if ( canonicalPartition != entity.pidInPartition() )
                {
                    auto const itId = entity.idInOthersPartitions().find( canonicalPartition );
                    if ( itId != entity.idInOthersPartitions().end() && itId->second != invalid_v<size_type> )
                        canonicalId = itId->second;
                }
                for ( auto const& [pid,entityId] : entity.idInOthersPartitions() )
                {
                    if ( entityId != invalid_v<size_type> && pid < canonicalPartition )
                    {
                        canonicalPartition = pid;
                        canonicalId = entityId;
                    }
                }
                return std::make_pair( canonicalPartition, canonicalId );
            };

            auto makeElementDofKey = [&]( element_type const& elt,
                                          uint16_type parentLid,
                                          descriptor_key_type& key ) -> bool
            {
                auto const layout = this->fe().localDofLayout( parentLid );
                auto const& attachment = layout.attachment;
                if ( !attachment.isValid() )
                    return false;

                switch ( attachment.entityDim )
                {
                case 0:
                {
                    const uint16_type nDofPerVertex = runtimeDofPerVertex();
                    if ( nDofPerVertex == 0 || attachment.entityId >= element_type::numVertices ||
                         attachment.ordinal >= nDofPerVertex )
                        return false;

                    auto const& point = elt.point( attachment.entityId );
                    auto const [entityPartition,entityId] = partitionCanonicalEntity( point );
                    key = makeKey( 0, entityPartition, entityId,
                                   attachment.ordinal, layout.component, attachment.kind );
                    return true;
                }
                case 1:
                {
                    const uint16_type nDofPerEdge = runtimeDofPerEdge();
                    if ( nDofPerEdge == 0 || attachment.ordinal >= nDofPerEdge )
                        return false;

                    if constexpr ( nDim == 1 )
                    {
                        key = makeKey( 1, elementClusterPartition( elt ), elementClusterId( elt ),
                                       attachment.ordinal, layout.component, attachment.kind );
                        return true;
                    }
                    else
                    {
                        if ( attachment.entityId >= element_type::numEdges )
                            return false;

                        uint16_type ordinal = attachment.ordinal;
                        auto const edgePermutation = elt.edgePermutation( attachment.entityId );
                        if ( edgePermutation.value() == edge_permutation_type::REVERSE_PERMUTATION )
                            ordinal = static_cast<uint16_type>( nDofPerEdge - 1 - attachment.ordinal );
                        else if ( edgePermutation.value() != edge_permutation_type::IDENTITY )
                            return false;

                        if constexpr ( nDim == 2 )
                        {
                            auto const facePtr = elt.facePtr( attachment.entityId );
                            if ( !facePtr )
                                return false;
                            auto const& edge = *facePtr;
                            auto const [entityPartition,entityId] = partitionCanonicalEntity( edge );
                            key = makeKey( 1, entityPartition, entityId,
                                           ordinal, layout.component, attachment.kind );
                            return true;
                        }
                        else
                        {
                            auto const edgePtr = elt.edgePtr( attachment.entityId );
                            if ( !edgePtr )
                                return false;
                            auto const& edge = *edgePtr;
                            auto const [entityPartition,entityId] = partitionCanonicalEntity( edge );
                            key = makeKey( 1, entityPartition, entityId,
                                           ordinal, layout.component, attachment.kind );
                            return true;
                        }
                    }
                }
                case 2:
                {
                    const uint16_type nDofPerFace = runtimeDofPerFace();
                    if ( nDofPerFace == 0 || attachment.ordinal >= nDofPerFace )
                        return false;

                    if constexpr ( nDim == 2 )
                    {
                        key = makeKey( 2, elementClusterPartition( elt ), elementClusterId( elt ),
                                       attachment.ordinal, layout.component, attachment.kind );
                        return true;
                    }
                    else if constexpr ( nDim == 3 )
                    {
                        if ( attachment.entityId >= element_type::numFaces )
                            return false;

                        auto const facePtr = elt.facePtr( attachment.entityId );
                        if ( !facePtr )
                            return false;

                        auto const& face = *facePtr;
                        uint16_type ordinal = attachment.ordinal;
                        auto const facePermutation = elt.facePermutation( attachment.entityId );
                        if ( facePermutation == face_permutation_type( 0 ) )
                            return false;
                        if ( nDofPerFace != 1 && facePermutation != face_permutation_type( face_permutation_type::IDENTITY ) )
                        {
                            if ( !this->hasValidFacePermutation( facePermutation, nDofPerFace ) )
                                return false;
                            auto const& permutation = this->facePermutationVector( facePermutation, nDofPerFace );
                            ordinal = static_cast<uint16_type>( permutation( attachment.ordinal ) );
                        }

                        auto const [entityPartition,entityId] = partitionCanonicalEntity( face );
                        key = makeKey( 2, entityPartition, entityId,
                                       ordinal, layout.component, attachment.kind );
                        return true;
                    }
                    else
                        return false;
                }
                case 3:
                {
                    const uint16_type nDofPerVolume = runtimeDofPerVolume();
                    if ( nDofPerVolume == 0 || attachment.ordinal >= nDofPerVolume )
                        return false;

                    key = makeKey( 3, elementClusterPartition( elt ), elementClusterId( elt ),
                                   attachment.ordinal, layout.component, attachment.kind );
                    return true;
                }
                default:
                    return false;
                }
            };

            std::map<descriptor_key_type,size_type> knownKeyToGc;
            std::map<descriptor_key_type,std::set<size_type>> unresolvedKeyToLocalDofs;
            std::map<descriptor_key_type,std::set<size_type>> localKeyToDofs;
            bool descriptorExtendedUsable = true;

            const uint16_type nLocalDof = runtimeNLocalDof();
            auto rangeElements = elements( mesh, entity_process_t::ALL );
            for ( auto const& eltWrap : rangeElements )
            {
                auto const& elt = unwrap_ref( eltWrap );
                if ( !this->isElementDone( elt.id() ) )
                    continue;

                for ( uint16_type parentLid = 0; parentLid < nLocalDof; ++parentLid )
                {
                    descriptor_key_type key;
                    if ( !makeElementDofKey( elt, parentLid, key ) )
                        continue;

                    auto const layout = this->fe().localDofLayout( parentLid );
                    auto const& gdof = this->localToGlobal( elt.id(), layout.parentLocalDofId, layout.component );
                    size_type const localDofIndex = gdof.index();
                    if ( localDofIndex >= this->M_mapGlobalProcessToGlobalCluster.size() )
                        continue;

                    localKeyToDofs[key].insert( localDofIndex );
                    auto const gcId = this->M_mapGlobalProcessToGlobalCluster[localDofIndex];
                    if ( gcId == invalid_v<size_type> )
                        unresolvedKeyToLocalDofs[key].insert( localDofIndex );
                    else
                    {
                        auto [itKnown, insertedKnown] = knownKeyToGc.emplace( key, gcId );
                        if ( !insertedKnown && itKnown->second != gcId )
                            descriptorExtendedUsable = false;
                    }
                }
            }

            std::vector<size_type> localPayload;
            localPayload.reserve( ( knownKeyToGc.size() + unresolvedKeyToLocalDofs.size() )*8 );
            for ( auto const& [key,gcId] : knownKeyToGc )
            {
                packKey( key, localPayload );
                localPayload.push_back( gcId );
            }
            for ( auto const& [key,localDofs] : unresolvedKeyToLocalDofs )
            {
                Feel::detail::ignore_unused_variable_warning( localDofs );
                if ( knownKeyToGc.find( key ) != knownKeyToGc.end() )
                    continue;
                packKey( key, localPayload );
                localPayload.push_back( invalid_v<size_type> );
            }

            std::map<rank_type,std::vector<size_type>> payloadToSend, payloadToRecv;
            for ( rank_type neighborRank : mesh.neighborSubdomains() )
                payloadToSend[neighborRank] = localPayload;

            std::vector<mpi::request> reqs( 2*mesh.neighborSubdomains().size() );
            int countRequest = 0;
            for ( rank_type neighborRank : mesh.neighborSubdomains() )
            {
                reqs[countRequest++] = this->worldComm().localComm().irecv( neighborRank, 0, payloadToRecv[neighborRank] );
                reqs[countRequest++] = this->worldComm().localComm().isend( neighborRank, 0, payloadToSend[neighborRank] );
            }
            mpi::wait_all( std::begin( reqs ), std::begin( reqs ) + countRequest );

            std::map<descriptor_key_type,std::map<rank_type,size_type>> receivedGcByKey;
            for ( auto const& [rankRecv,payload] : payloadToRecv )
            {
                CHECK( payload.size() % 8 == 0 )
                    << fmt::format( "[dof-mpi-key][extended] invalid key payload size {} from rank {}", payload.size(), rankRecv );
                for ( size_type k = 0; k < payload.size(); k += 8 )
                {
                    auto key = unpackKey( payload, k );
                    size_type const gcId = payload[k+7];
                    if ( unresolvedKeyToLocalDofs.find( key ) != unresolvedKeyToLocalDofs.end() ||
                         localKeyToDofs.find( key ) != localKeyToDofs.end() )
                        receivedGcByKey[key][rankRecv] = gcId;
                }
            }

            size_type syncedDofs = 0;
            size_type unresolvedDofs = 0;
            for ( auto const& [key,localDofs] : unresolvedKeyToLocalDofs )
            {
                size_type gcId = invalid_v<size_type>;
                auto const itKnown = knownKeyToGc.find( key );
                if ( itKnown != knownKeyToGc.end() )
                    gcId = itKnown->second;
                else
                {
                    auto const itRecv = receivedGcByKey.find( key );
                    if ( itRecv != receivedGcByKey.end() && !itRecv->second.empty() )
                    {
                        for ( auto const& [rankRecv,receivedGcId] : itRecv->second )
                        {
                            Feel::detail::ignore_unused_variable_warning( rankRecv );
                            if ( receivedGcId != invalid_v<size_type> )
                            {
                                gcId = receivedGcId;
                                break;
                            }
                        }
                    }
                }

                if ( gcId == invalid_v<size_type> )
                {
                    unresolvedDofs += localDofs.size();
                    continue;
                }

                for ( size_type localDofIndex : localDofs )
                {
                    this->M_mapGlobalProcessToGlobalCluster[localDofIndex] = gcId;
                    rank_type const activeProcId = this->procOnGlobalCluster( gcId );
                    if ( activeProcId != myRank )
                        this->addNeighborSubdomain( activeProcId );
                    ++syncedDofs;
                }
            }

            size_type ownerSharedDofs = 0;
            for ( auto const& [key,localDofs] : localKeyToDofs )
            {
                auto const itRecv = receivedGcByKey.find( key );
                if ( itRecv == receivedGcByKey.end() )
                    continue;

                for ( size_type localDofIndex : localDofs )
                {
                    if ( localDofIndex >= this->M_mapGlobalProcessToGlobalCluster.size() )
                        continue;
                    size_type const gcId = this->M_mapGlobalProcessToGlobalCluster[localDofIndex];
                    if ( gcId == invalid_v<size_type> || gcId >= this->nDof() )
                        continue;
                    if ( this->procOnGlobalCluster( gcId ) != myRank )
                        continue;

                    for ( auto const& [sharedRank,receivedGcId] : itRecv->second )
                    {
                        if ( sharedRank == myRank )
                            continue;
                        if ( receivedGcId != invalid_v<size_type> && receivedGcId != gcId )
                            continue;
                        this->M_activeDofSharedOnCluster[localDofIndex].insert( sharedRank );
                        this->addNeighborSubdomain( sharedRank );
                        ++ownerSharedDofs;
                    }
                }
            }

            int const localDescriptorExtendedSuccess = ( descriptorExtendedUsable && unresolvedDofs == 0 ) ? 1 : 0;
            int globalDescriptorExtendedSuccess = 0;
            mpi::all_reduce( this->worldComm().localComm(), localDescriptorExtendedSuccess,
                             globalDescriptorExtendedSuccess, mpi::minimum<int>() );

            if ( globalDescriptorExtendedSuccess )
            {
                VLOG(1) << "[dof-mpi-key][extended] descriptor-key extended synchronization"
                        << " known_keys=" << knownKeyToGc.size()
                        << " unresolved_keys=" << unresolvedKeyToLocalDofs.size()
                        << " synced_dofs=" << syncedDofs
                        << " owner_shared_dofs=" << ownerSharedDofs;
                return;
            }

            LOG(WARNING) << "[dof-mpi-key][extended] fallback to legacy extended synchronization"
                         << " descriptorExtendedUsable=" << descriptorExtendedUsable
                         << " unresolved_dofs=" << unresolvedDofs
                         << " global_success=" << globalDescriptorExtendedSuccess;
        }
    }

    if ( is_continuous && mesh.isPeriodic() && nProc > 1 )
    {
        using periodic_dof_key_type = std::tuple<uint16_type,size_type,uint16_type>; // (entity dim, global dof key, component)

        std::map<periodic_dof_key_type,size_type> periodicAllLocalKeyToDof;
        std::map<periodic_dof_key_type,size_type> periodicExtendedKeyToDof;
        this->forEachGlobalDofKeyEntry( [&]( auto const& entry )
        {
            if ( entry.localDofIndex >= this->M_n_localWithGhost_df[myRank] )
                return;

            periodic_dof_key_type key{ std::get<0>( entry.key ),
                                       std::get<1>( entry.key ),
                                       entry.component };
            periodicAllLocalKeyToDof.emplace( key, entry.localDofIndex );
            if ( entry.localDofIndex >= start_next_free_dof && entry.localDofIndex < next_free_dof )
                periodicExtendedKeyToDof.emplace( key, entry.localDofIndex );
        } );

        if ( periodicExtendedKeyToDof.empty() )
        {
            VLOG(1) << "[periodic][mpi][extended] no extended periodic dof key to synchronize";
            return;
        }

        std::vector<size_type> periodicLocalPayload;
        periodicLocalPayload.reserve( periodicAllLocalKeyToDof.size()*5 );
        for ( auto const& [key,localDofIndex] : periodicAllLocalKeyToDof )
        {
            size_type gcId = invalid_v<size_type>;
            if ( localDofIndex < this->M_mapGlobalProcessToGlobalCluster.size() )
                gcId = this->M_mapGlobalProcessToGlobalCluster[localDofIndex];
            bool const isGhostLocal = this->dofGlobalProcessIsGhost( localDofIndex );

            periodicLocalPayload.push_back( static_cast<size_type>( std::get<0>( key ) ) );
            periodicLocalPayload.push_back( std::get<1>( key ) );
            periodicLocalPayload.push_back( static_cast<size_type>( std::get<2>( key ) ) );
            periodicLocalPayload.push_back( gcId );
            periodicLocalPayload.push_back( isGhostLocal ? size_type( 1 ) : size_type( 0 ) );
        }

        std::vector<std::vector<size_type>> periodicPayloadPerRank;
        mpi::all_gather( this->worldComm().localComm(), periodicLocalPayload, periodicPayloadPerRank );
        bool periodicPayloadValid = ( periodicPayloadPerRank.size() == static_cast<size_type>( nProc ) );
        if ( !periodicPayloadValid )
            LOG(WARNING) << fmt::format( "[periodic][mpi][extended] invalid all_gather size {} vs {}, fallback to legacy extended exchange",
                                         periodicPayloadPerRank.size(), nProc );

        std::map<periodic_dof_key_type,size_type> periodicGcByKey;
        std::map<periodic_dof_key_type,std::set<rank_type>> periodicRanksByKey;
        using periodic_entry_type = std::tuple<rank_type,size_type,bool>; // rank, gc id, is ghost
        std::map<periodic_dof_key_type,std::vector<periodic_entry_type>> periodicEntriesByKey;
        size_type periodicInconsistentGc = 0;
        if ( periodicPayloadValid )
        {
            for ( rank_type rank = 0; rank < nProc; ++rank )
            {
                auto const& flat = periodicPayloadPerRank[rank];
                if ( flat.size() % 5 != 0 )
                {
                    LOG(WARNING) << fmt::format( "[periodic][mpi][extended] invalid payload size {} on rank {}, fallback to legacy extended exchange",
                                                 flat.size(), rank );
                    periodicPayloadValid = false;
                    break;
                }

                for ( size_type k = 0; k < flat.size(); k += 5 )
                {
                    periodic_dof_key_type key{
                        static_cast<uint16_type>( flat[k] ),
                        flat[k+1],
                        static_cast<uint16_type>( flat[k+2] )
                    };
                    size_type const gcId = flat[k+3];
                    bool const isGhost = ( flat[k+4] != 0 );
                    periodicRanksByKey[key].insert( rank );
                    periodicEntriesByKey[key].emplace_back( rank, gcId, isGhost );
                }
            }

            if ( periodicPayloadValid )
            {
                for ( auto const& [key,entries] : periodicEntriesByKey )
                {
                    rank_type ownerRank = invalid_v<rank_type>;
                    for ( auto const& [entryRank,entryGcId,entryIsGhost] : entries )
                    {
                        Feel::detail::ignore_unused_variable_warning( entryGcId );
                        if ( entryIsGhost )
                            continue;
                        if ( ownerRank == invalid_v<rank_type> || entryRank < ownerRank )
                            ownerRank = entryRank;
                    }
                    if ( ownerRank == invalid_v<rank_type> )
                    {
                        for ( auto const& [entryRank,entryGcId,entryIsGhost] : entries )
                        {
                            Feel::detail::ignore_unused_variable_warning( entryIsGhost );
                            if ( entryGcId == invalid_v<size_type> )
                                continue;
                            if ( ownerRank == invalid_v<rank_type> || entryRank < ownerRank )
                                ownerRank = entryRank;
                        }
                    }

                    size_type selectedGc = invalid_v<size_type>;
                    if ( ownerRank != invalid_v<rank_type> )
                    {
                        for ( auto const& [entryRank,entryGcId,entryIsGhost] : entries )
                        {
                            Feel::detail::ignore_unused_variable_warning( entryIsGhost );
                            if ( entryRank == ownerRank && entryGcId != invalid_v<size_type> )
                            {
                                selectedGc = entryGcId;
                                break;
                            }
                        }
                    }
                    if ( selectedGc == invalid_v<size_type> )
                    {
                        for ( auto const& [entryRank,entryGcId,entryIsGhost] : entries )
                        {
                            Feel::detail::ignore_unused_variable_warning( entryRank );
                            Feel::detail::ignore_unused_variable_warning( entryIsGhost );
                            if ( entryGcId != invalid_v<size_type> )
                            {
                                selectedGc = entryGcId;
                                break;
                            }
                        }
                    }
                    if ( selectedGc == invalid_v<size_type> )
                        continue;

                    periodicGcByKey.emplace( key, selectedGc );
                    for ( auto const& [entryRank,entryGcId,entryIsGhost] : entries )
                    {
                        Feel::detail::ignore_unused_variable_warning( entryRank );
                        if ( entryIsGhost || entryGcId == invalid_v<size_type> )
                            continue;
                        if ( entryGcId != selectedGc )
                            ++periodicInconsistentGc;
                    }
                }
            }
        }

        size_type periodicExtendedSyncedDofs = 0;
        size_type periodicExtendedUnresolvedDofs = 0;
        if ( !periodicPayloadValid )
        {
            periodicExtendedUnresolvedDofs = periodicExtendedKeyToDof.size();
        }
        else
        {
            for ( auto const& [key,localDofIndex] : periodicExtendedKeyToDof )
            {
                auto const itGc = periodicGcByKey.find( key );
                if ( itGc == periodicGcByKey.end() || itGc->second == invalid_v<size_type> )
                {
                    ++periodicExtendedUnresolvedDofs;
                    continue;
                }

                size_type const gcId = itGc->second;
                if ( localDofIndex >= this->M_mapGlobalProcessToGlobalCluster.size() )
                {
                    ++periodicExtendedUnresolvedDofs;
                    continue;
                }
                if ( this->M_mapGlobalProcessToGlobalCluster[localDofIndex] == invalid_v<size_type> )
                    this->M_mapGlobalProcessToGlobalCluster[localDofIndex] = gcId;
                else if ( this->M_mapGlobalProcessToGlobalCluster[localDofIndex] != gcId )
                {
                    ++periodicInconsistentGc;
                    this->M_mapGlobalProcessToGlobalCluster[localDofIndex] = gcId;
                }

                if ( gcId >= this->nDof() )
                {
                    ++periodicExtendedUnresolvedDofs;
                    continue;
                }

                rank_type const activeProcId = this->procOnGlobalCluster( gcId );
                if ( activeProcId != myRank )
                    this->addNeighborSubdomain( activeProcId );
                ++periodicExtendedSyncedDofs;
            }
        }

        if ( periodicExtendedUnresolvedDofs == 0 )
        {
            size_type periodicExtendedOwnerSharedDofs = 0;
            for ( auto const& [key,localDofIndex] : periodicAllLocalKeyToDof )
            {
                if ( localDofIndex >= this->M_mapGlobalProcessToGlobalCluster.size() )
                    continue;
                size_type const gcId = this->M_mapGlobalProcessToGlobalCluster[localDofIndex];
                if ( gcId == invalid_v<size_type> || gcId >= this->nDof() )
                    continue;

                rank_type const activeProcId = this->procOnGlobalCluster( gcId );
                if ( activeProcId != myRank )
                    continue;

                auto itRanks = periodicRanksByKey.find( key );
                if ( itRanks == periodicRanksByKey.end() || itRanks->second.size() <= 1 )
                    continue;

                for ( rank_type sharedRank : itRanks->second )
                {
                    if ( sharedRank == myRank )
                        continue;
                    this->M_activeDofSharedOnCluster[localDofIndex].insert( sharedRank );
                    this->addNeighborSubdomain( sharedRank );
                    ++periodicExtendedOwnerSharedDofs;
                }
            }

            VLOG(1) << "[periodic][mpi][extended] periodic global-cluster synchronization"
                    << " all_local_keys=" << periodicAllLocalKeyToDof.size()
                    << " extended_keys=" << periodicExtendedKeyToDof.size()
                    << " synced_extended_dofs=" << periodicExtendedSyncedDofs
                    << " owner_shared_dofs=" << periodicExtendedOwnerSharedDofs
                    << " inconsistent_gc=" << periodicInconsistentGc;
            return;
        }

        LOG(WARNING) << "[periodic][mpi][extended] unresolved periodic dofs in extended sync, fallback to legacy exchange"
                     << " unresolved=" << periodicExtendedUnresolvedDofs
                     << " synced=" << periodicExtendedSyncedDofs
                     << " inconsistent_gc=" << periodicInconsistentGc;
    }

    //------------------------------------------------------------------------------//
    // TODO : maybe we can only apply send/recv from active elt (with idInOtherspartitions.size()>0) to all ghost elements

    std::map< rank_type, std::vector<index_type> > dataToSend, dataToRecv, dataMemory;

    // get dofs in extended part
    for ( auto const& ghostEltWrap : ghostEltRange )
    {
        auto const& ghostElt = boost::unwrap_ref( ghostEltWrap );
        const rank_type processIdOfGhost = ghostElt.processId();
        const size_type eltIdOfGhost = ghostElt.id();
        const size_type eltIdInOtherPartOfGhost = ghostElt.idInOthersPartitions(processIdOfGhost);

        dataToSend[processIdOfGhost].push_back( eltIdInOtherPartOfGhost );
        dataMemory[processIdOfGhost].push_back( ghostElt.id() );
    }

    //------------------------------------------------------------------------------//

    int nbMaxRequest = 2*mesh.neighborSubdomains().size();
    std::vector<mpi::request> reqs( nbMaxRequest );
    int countRequest = 0;
    // send/recv of data size
    std::map<rank_type,std::size_t> sizeRecv, sizeSended;
    for ( rank_type neighborRank : mesh.neighborSubdomains() )
    {
        sizeSended[neighborRank] = dataToSend[neighborRank].size();
        reqs[countRequest++] = this->worldComm().localComm().isend( neighborRank, 0, sizeSended[neighborRank] );
        reqs[countRequest++] = this->worldComm().localComm().irecv( neighborRank, 0, sizeRecv[neighborRank] );
    }
    // wait all requests
    mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
    countRequest = 0;
    // step 1 :send/recv of data
    for ( rank_type neighborRank : mesh.neighborSubdomains() )
    {
        std::size_t nRecvData = sizeRecv[neighborRank];
        dataToRecv[neighborRank].resize( nRecvData );
        if ( nRecvData > 0 )
            reqs[countRequest++] = this->worldComm().localComm().irecv( neighborRank , 0, dataToRecv[neighborRank].data(), nRecvData );
        std::size_t nSendData = sizeSended[neighborRank];
        if ( nSendData > 0 )
            reqs[countRequest++] = this->worldComm().localComm().isend( neighborRank , 0, dataToSend[neighborRank].data(), nSendData );
    }
    // step 1 :wait all requests
    mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
    countRequest = 0;
    //------------------------------------------------------------------------------//

    // step 2 : treat recv and prepare data for response
    std::map< rank_type, std::vector<std::vector<size_type>> > dataToSendStep2, dataToRecvStep2;
    for ( auto const& [rankRecv,eltIds] : dataToRecv )
    {
        auto & dataToSendStep2AtRank = dataToSendStep2[rankRecv];
        dataToSendStep2AtRank.reserve( eltIds.size() );
        for ( index_type eltId : eltIds )
        {
            dataToSendStep2AtRank.push_back( this->getIndicesOnGlobalCluster( eltId ) );
            DCHECK( dataToSendStep2AtRank.back().size() > 0 ) << "no dof found";
            // for ( size_type dofIndex : this->getIndices( eltId ) )
            //     if ( !this->dofGlobalProcessIsGhost( dofIndex ) )
            //         this->M_activeDofSharedOnCluster[dofIndex].insert(rankRecv);
        }
    }

    // step 2 :send/recv of data
    for ( rank_type neighborRank : mesh.neighborSubdomains() )
    {
        std::size_t nRecvData = sizeSended[neighborRank]; // use inverse size
        dataToRecvStep2[neighborRank].resize( nRecvData );
        if ( nRecvData > 0 )
            reqs[countRequest++] = this->worldComm().localComm().irecv( neighborRank , 0, dataToRecvStep2[neighborRank].data(), nRecvData );
        std::size_t nSendData = sizeRecv[neighborRank]; // use inverse size
        DCHECK( dataToSendStep2[neighborRank].size() == nSendData ) << fmt::format("incompatible data size {} vs {}",dataToSendStep2[neighborRank].size(), nSendData);
        if ( nSendData > 0 )
            reqs[countRequest++] = this->worldComm().localComm().isend( neighborRank , 0, dataToSendStep2[neighborRank].data(), nSendData );
    }
    // step 2 :wait all requests
    mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
    countRequest = 0;

    for ( auto const& [rankRecv,dofIdInElt] : dataToRecvStep2 )
    {
        if ( dofIdInElt.empty() )
            continue;
        auto & dataMemoryAtRank = dataMemory.at( rankRecv );
        DCHECK( dataMemoryAtRank.size() == dofIdInElt.size() )<< fmt::format("incompatible data size {} vs {}",dataMemoryAtRank.size(), dofIdInElt.size() );
        for (int k=0;k<dofIdInElt.size();++k)
        {
            size_type eltId = dataMemoryAtRank[k];
            auto const& dofGlobalClusterIds = dofIdInElt[k];
            auto dofIndices = this->getIndices( eltId );
            DCHECK( dofGlobalClusterIds.size() == dofIndices.size() ) << fmt::format("incompatible data size {} vs {}",dofGlobalClusterIds.size(), dofIndices.size());
            for (int l=0;l<dofGlobalClusterIds.size();++l)
            {
                size_type dofIndex = dofIndices[l];
                size_type dofGlobalClusterIndex = dofGlobalClusterIds[l];
                DCHECK( dofIndex < this->M_mapGlobalProcessToGlobalCluster.size() )<< fmt::format("incompatible dofIndex {} vs {}",dofIndex, this->M_mapGlobalProcessToGlobalCluster.size() );
                if ( this->M_mapGlobalProcessToGlobalCluster[dofIndex] == invalid_v<size_type> )
                    this->M_mapGlobalProcessToGlobalCluster[dofIndex] = dofGlobalClusterIndex;
                else
                    CHECK( dofGlobalClusterIndex == this->M_mapGlobalProcessToGlobalCluster[dofIndex] ) << fmt::format("dof index should be identical {} vs {} ", dofGlobalClusterIndex, this->M_mapGlobalProcessToGlobalCluster[dofIndex] );

                rank_type activeProcId = this->procOnGlobalCluster( dofGlobalClusterIndex );
                if ( activeProcId != myRank )
                    this->addNeighborSubdomain( activeProcId ); // TODO try outside of loop
                else
                    this->M_activeDofSharedOnCluster[dofIndex].insert(rankRecv);
            }
        }
    }
}

} // namespace Feel

#endif /* FEELPP_DOFTABLE_MPI_HPP */
