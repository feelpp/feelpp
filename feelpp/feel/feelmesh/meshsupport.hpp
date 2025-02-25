/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date:6 June 2017

 Copyright (C) 2017 Feel++ Consortium

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

#if !defined(FEELPP_MESH_SUPPORT_HPP)
#define FEELPP_MESH_SUPPORT_HPP 1

#include <feel/feelmesh/meshsupportbase.hpp>
#include <feel/feelmesh/traits.hpp>
#include <feel/feelmesh/ranges.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feeldiscr/localization.hpp>

namespace Feel
{

/**
 * \brief Description of a mesh support.
 * allows to build a function space on a range of elements
 */
template<typename MeshType>
class MeshSupport : public MeshSupportBase, public std::enable_shared_from_this<MeshSupport<MeshType>>
{
    enum class _face_attributes{ on_boundary=0, intraprocess, interprocess };
public :
    using super_type = MeshSupportBase;
    using mesh_type = typename MeshTraits<MeshType>::mesh_type;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;
    using range_elements_type = Range<mesh_type,MESH_ELEMENTS>;
    using range_faces_type = Range<mesh_type,MESH_FACES>;
    using element_type = typename mesh_type::element_type;
    using face_type = typename mesh_type::face_type;
    using point_interprocess_map_type = typename mesh_type::point_interprocess_map_type;

    static constexpr int nDim = mesh_type::nDim;

    MeshSupport() = default;
    MeshSupport( mesh_ptrtype const& mesh ) : MeshSupport( mesh, elements(mesh), true ) {}

    MeshSupport( mesh_ptrtype const& mesh, range_elements_type const& rangeElements, bool fullsupport = false )
        :
        M_mesh( mesh ),
        M_rangeElements( rangeElements ),
        M_isFullSupport( fullsupport )
        {
            if ( !M_isFullSupport )
            {
                M_localizationToolPartialSupport = std::make_shared<Localization<mesh_type>>();
                M_localizationToolPartialSupport->setMesh( M_mesh, M_rangeElements, false );
                M_mesh->attachMeshSupport( this );
            }
            this->updateForUse();
        }

    ~MeshSupport() override
        {
            if ( M_mesh )
                M_mesh->detachMeshSupport( this );
        }

    mesh_ptrtype const& mesh() const { return M_mesh; }
    worldcomm_ptr_t const& worldCommPtr() const { return M_mesh->worldCommPtr(); }

    bool isFullSupport() const override { return M_isFullSupport; }
    bool isPartialSupport() const override { return !M_isFullSupport; }

    std::shared_ptr<Localization<mesh_type>> tool_localization() const { return this->isPartialSupport()? M_localizationToolPartialSupport : M_mesh->tool_localization(); }

    template <entity_filter_t FF, entity_process_t EPT, typename ... Ts>
    range_elements_type elementsFilterImpl( Ts&&... ts ) const;

    template <entity_filter_t FF, typename ... Ts>
    range_elements_type elementsFilter( entity_process_t ept, Ts&&... ts ) const;

    template <entity_process_t EPT = entity_process_t::LOCAL_ONLY>
    range_elements_type rangeElementsProcessId( rank_type part ) const;

    template <entity_process_t EPT = entity_process_t::LOCAL_ONLY>
    range_elements_type rangeElementsMarkerByType( uint16_type markerType, std::set<flag_type> const& markerFlags, rank_type part ) const;


    template <entity_filter_t FF, entity_process_t EPT, typename ... Ts>
    range_faces_type facesFilterImpl( Ts&&... ts ) const;

    template <entity_filter_t FF, typename ... Ts>
    range_faces_type facesFilter( entity_process_t ept, Ts&&... ts ) const;

    template <entity_process_t EPT = entity_process_t::LOCAL_ONLY>
    range_faces_type rangeFacesProcessId( rank_type part ) const;

    template <entity_process_t EPT = entity_process_t::LOCAL_ONLY>
    range_faces_type rangeBoundaryFaces( rank_type part ) const;

    template <entity_process_t EPT = entity_process_t::LOCAL_ONLY>
    range_faces_type rangeInternalFaces( rank_type part ) const;

    template <entity_process_t EPT = entity_process_t::LOCAL_ONLY>
    range_faces_type rangeFacesMarkerByType( uint16_type markerType, std::set<flag_type> const& markerFlags, rank_type part ) const;


    range_faces_type rangeInterProcessFaces( rank_type part, rank_type neighbor_pid ) const;

    size_type numElements() const override
        {
            if ( M_isFullSupport )
                return M_mesh->numElements();
            else
                return M_rangeMeshElementsIdsPartialSupport.size();
        }
    bool hasElement( size_type eltId ) const override
        {
            if ( M_isFullSupport )
                return M_mesh->hasElement( eltId );
            else
                return M_rangeMeshElementsIdsPartialSupport.find( eltId ) != M_rangeMeshElementsIdsPartialSupport.end();
        }
    bool hasGhostElement( size_type eltId ) const
        {
            if ( !this->hasElement( eltId ) )
                return false;
            return M_mesh->element( eltId ).isGhostCell();
        }

    template <typename FaceType>
    bool isGhostFace( FaceType const& face ) const
        {
            if constexpr ( !std::is_same_v<FaceType,face_type> )
                return false;
            else
            {
                if ( M_isFullSupport )
                    return face.isGhostFace();
                else
                {
                    if ( !face.isInterProcessDomain() )
                        return false;
                    auto const& elt0 = face.element(0);
                    auto const& elt1 = face.element(1);
                    bool hasElt0 = this->hasElement( elt0.id() );
                    bool hasElt1 = this->hasElement( elt1.id() );
                    if ( hasElt0 && hasElt1 )
                        return face.isGhostFace();
                    else if ( hasElt0 )
                        return elt0.isGhostCell();
                    else if ( hasElt1 )
                        return elt1.isGhostCell();
                    else
                        return true;
                }
            }
        }

    //! return true if the point id is interprocess of current partition
    bool isInterprocessPoints( index_type pointId ) const
        {
            return this->findInterprocessPoints( pointId ).first;
        }
    //! try to find data of interprocess of current partition point id and return pair(bool,iterator)
    std::pair<bool,typename point_interprocess_map_type::const_iterator> findInterprocessPoints( index_type pointId ) const
        {
            auto itFind = M_interprocessPoints.find( pointId );
            return std::make_pair( itFind != M_interprocessPoints.end(), itFind );
        }

private :
    void updateForUse();
    void updateParallelDataPartialSupport();
    void resetLocalizationTool() override
        {
            if ( M_localizationToolPartialSupport )
                M_localizationToolPartialSupport->reset();
        }

private :
    mesh_ptrtype M_mesh;
    range_elements_type M_rangeElements;
    std::shared_ptr<Localization<mesh_type>> M_localizationToolPartialSupport;

    std::vector< std::reference_wrapper<const element_type> > M_orderedElements;
    std::vector< std::tuple<std::reference_wrapper<const face_type>, _face_attributes > > M_orderedFaces;
    std::unordered_set<size_type> M_rangeMeshElementsIdsPartialSupport;
    point_interprocess_map_type M_interprocessPoints;

    bool M_isFullSupport;
};


template <typename MeshType>
void
MeshSupport<MeshType>::updateForUse()
{
    if ( M_isFullSupport )
        return;

    this->updateParallelDataPartialSupport();


    // update subentities (faces) on active elements
    std::unordered_map<size_type,std::pair<const face_type*,_face_attributes /*uint8_type*/>> faceInRange;
    for ( auto const& eltWrap : M_orderedElements ) // TODO get end index of active elements
    {
        auto const& elt = eltWrap.get();
        //auto const& elt = unwrap_ref( eltWrap );
        if ( elt.isGhostCell() )
            continue;
        for ( uint16_type i = 0; i < mesh_type::element_type::numTopologicalFaces; ++i )
        {
            if ( !elt.facePtr(i) )
                continue;
            const face_type* facePtr = elt.facePtr(i);
            size_type faceId = facePtr->id();
            auto const& face = elt.face(i);

            if ( face.isInterProcessDomain() )
            {
                auto const& elt0 = face.element0();
                auto const& elt1 = face.element1();
                if ( this->hasElement( elt0.id() ) && this->hasElement( elt1.id() ) )
                {
                    faceInRange[faceId] = std::make_pair(facePtr,_face_attributes::interprocess);
                    continue;
                }
            }
            if ( faceInRange.find( faceId ) != faceInRange.end() )
                faceInRange[faceId].second = _face_attributes::intraprocess;
            else
                faceInRange[faceId] = std::make_pair(facePtr,_face_attributes::on_boundary);
        }
    }


    std::map< rank_type, std::vector<size_type> > dataToSend, dataToRecv;
    std::map< rank_type, std::vector< std::reference_wrapper<const element_type> > > dataMemory;

    for ( auto const& eltWrap : M_orderedElements ) // TODO get start index of ghost elements
    {
        auto const& elt = eltWrap.get();
        if ( !elt.isGhostCell() )
            continue;
        dataToSend[elt.processId()].push_back( elt.idInOthersPartitions( elt.processId() ) );
        dataMemory[elt.processId()].push_back( std::cref(elt) );
    }

    // mpi comm
    int neighborSubdomains = M_mesh->neighborSubdomains().size();
    int nbMaxRequest = 2*neighborSubdomains;
    std::vector<mpi::request> reqs( nbMaxRequest );
    int countRequest = 0;
    std::map<rank_type,std::size_t> sizeRecv;
    std::map<rank_type,std::size_t> sizeSend;

    // get size of data to transfer
    for ( rank_type neighborRank : M_mesh->neighborSubdomains() )
    {
        sizeSend[neighborRank] = dataToSend[neighborRank].size();
        reqs[countRequest++] = M_mesh->worldComm().localComm().isend( neighborRank, 0, sizeSend[neighborRank] );
        reqs[countRequest++] = M_mesh->worldComm().localComm().irecv( neighborRank, 0, sizeRecv[neighborRank] );
    }
    // wait all requests
    mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
    countRequest = 0;

    // send/recv data
    for ( rank_type neighborRank : M_mesh->neighborSubdomains() )
    {
        std::size_t nSendData = dataToSend[neighborRank].size();
        if ( nSendData > 0 )
            reqs[countRequest++] = M_mesh->worldComm().localComm().isend( neighborRank , 0, dataToSend[neighborRank].data(), nSendData );
        std::size_t nRecvData = sizeRecv[neighborRank];
        dataToRecv[neighborRank].resize( nRecvData );
        if ( nRecvData > 0 )
            reqs[countRequest++] = M_mesh->worldComm().localComm().irecv( neighborRank , 0, dataToRecv[neighborRank].data(), nRecvData );
    }
    // wait all requests
    mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
    countRequest = 0;

    // step2 : from active elts, prepare mpi data of subentities required (from ghost elts)
    std::map< rank_type, std::vector<std::tuple<std::vector<_face_attributes>>> > dataToSendStep2, dataToRecvStep2;
    for ( auto const& [rankRecv,eltIds] : dataToRecv )
    {
        auto & dataToSendStep2OnRank = dataToSendStep2[rankRecv];
        dataToSendStep2OnRank.resize( eltIds.size() );
        for ( int k=0; k<eltIds.size(); ++k )
        {
            size_type eltId = eltIds[k];
            auto & [dataToSendStep2OnFacesOnElt] = dataToSendStep2OnRank[k];
            auto const& elt = M_mesh->element( eltId );
            dataToSendStep2OnFacesOnElt.resize( elt.nTopologicalFaces() );
            for ( uint16_type i = 0; i < elt.nTopologicalFaces(); ++i )
            {
                if ( !elt.facePtr(i) )
                    continue;
                auto const& face = elt.face(i);
                auto itFindFace = faceInRange.find( face.id() );
                CHECK( itFindFace != faceInRange.end() ) << "face not registered, something wrong";
                dataToSendStep2OnFacesOnElt[i] = std::get<1>( itFindFace->second );
            }
        }
    }
    // step2 : send/recv data
    for ( rank_type neighborRank : M_mesh->neighborSubdomains() )
    {
        std::size_t nSendData = dataToSendStep2[neighborRank].size();
        if ( nSendData > 0 )
            reqs[countRequest++] = M_mesh->worldComm().localComm().isend( neighborRank , 0, dataToSendStep2[neighborRank].data(), nSendData );
        std::size_t nRecvData = sizeSend[neighborRank]; // use size from send of step1
        dataToRecvStep2[neighborRank].resize( nRecvData );
        if ( nRecvData > 0 )
            reqs[countRequest++] = M_mesh->worldComm().localComm().irecv( neighborRank , 0, dataToRecvStep2[neighborRank].data(), nRecvData );
    }
    // step2 : wait all requests
    mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
    countRequest = 0;

    for ( auto const& [rankRecv,dataSubentitiesByElt] : dataToRecvStep2 )
    {
        auto const& dataMemoryOnRank = dataMemory[rankRecv];
        CHECK( dataMemoryOnRank.size() == dataSubentitiesByElt.size() ) << fmt::format( "incompatible size: {} vs {}", dataMemoryOnRank.size(), dataSubentitiesByElt.size() );
        for (int k=0;k<dataSubentitiesByElt.size();++k)
        {
            auto dataFaces = std::get<0>( dataSubentitiesByElt[k] );
            auto const& ghostElt = dataMemoryOnRank[k].get();
            CHECK( dataFaces.size() == ghostElt.nTopologicalFaces() ) << fmt::format( "incompatible size: {} vs {}", dataFaces.size(), ghostElt.nTopologicalFaces() );
            for ( uint16_type i = 0; i < ghostElt.nTopologicalFaces(); ++i )
            {
                if ( !ghostElt.facePtr(i) )
                    continue;
                auto const& face = ghostElt.face(i);
                size_type faceId = face.id();
                if ( faceInRange.find( faceId ) == faceInRange.end() )
                    faceInRange[faceId] = std::make_pair(std::addressof(face),dataFaces[i]);
            }
        }
    }





    M_orderedFaces.reserve( faceInRange.size() );
    for ( auto const& faceWrap : M_mesh->orderedFaces() )
    {
        auto const& face = unwrap_ref( faceWrap );
        auto itFind = faceInRange.find( face.id() );
        if ( itFind == faceInRange.end() )
            continue;
        M_orderedFaces.push_back( std::make_tuple(std::cref(face), std::get<1>( itFind->second ) ) );
    }







    // update interprocess entties
    // TODO optimisation if we have interprocessfaces

    std::unordered_map<index_type,std::tuple<bool,std::set<rank_type>>> pointsInterprocessDetection; // ( pt id -> ( isOnActiveElt, isOnGhostEltRanks ) )
    std::unordered_map<index_type,std::tuple<bool,std::set<rank_type>>> edgesInterprocessDetection; // ( edge id -> ( isOnActiveElt, isOnGhostEltRanks ) )
#if 0
    pointsInterprocessDetection.reserve( std::distance( this->beginOrderedPoint(),
                                                        this->endOrderedPoint() ) );
    if constexpr ( nDim == 3 )
        edgesInterprocessDetection.reserve( std::distance( this->beginOrderedEdge(),
                                                           this->endOrderedEdge() ) );
#endif
    auto itPointIpDetect = pointsInterprocessDetection.begin();
    auto itEdgeIpDetect = edgesInterprocessDetection.begin();
    for ( auto const& eltWrap : M_orderedElements ) // TODO get end index of active elements
    {
        auto const& elt = eltWrap.get();

        // nothing to do if no neighbor process
        if ( elt.idInOthersPartitions().empty() )
            continue;

        for ( uint16_type n=0; n < elt.nPoints(); n++ )
        {
            auto const& point = elt.point( n );
            std::tie( itPointIpDetect,std::ignore ) = pointsInterprocessDetection.try_emplace( point.id(), false, std::set<rank_type>{} );
            if ( elt.isGhostCell() )
                std::get<1>( itPointIpDetect->second ).insert( elt.processId() );
            else
                std::get<0>( itPointIpDetect->second ) = true;
        }
#if 0
        if constexpr ( nDim == 3 )
        {
            for ( size_type j = 0; j < elt.nEdges(); j++ )
            {
                if ( !elt.edgePtr( j ) )
                    continue;
                auto const& edge = elt.edge( j );
                std::tie( itEdgeIpDetect,std::ignore ) = edgesInterprocessDetection.try_emplace( edge.id(), false, std::set<rank_type>{} );
                if ( elt.isGhostCell() )
                    std::get<1>( itEdgeIpDetect->second ).insert( elt.processId() );
                else
                    std::get<0>( itEdgeIpDetect->second ) = true;
            }
        }
#endif
    }

    M_interprocessPoints.clear();
    for ( auto const& [pointId,ipData] : pointsInterprocessDetection )
    {
        if ( !std::get<0>( ipData ) ) // not on current process
            continue;
        if ( std::get<1>( ipData ).empty() ) // not on neighbor process
            continue;
        //M_interprocessPoints.try_emplace( pointId, std::move( std::get<1>( ipData ) ) );
        M_interprocessPoints.try_emplace( pointId, std::get<1>( ipData ) );
    }





}


template <typename MeshType>
void
MeshSupport<MeshType>::updateParallelDataPartialSupport()
{
    wc(this)->print( fmt::format( "[updateParallelDataPartialSupport] starts..." ), FLAGS_v > 1, FLAGS_v > 0, FLAGS_v > 1 );
    if ( M_mesh->worldComm().localSize() == 1 )
    {
        for ( auto const& eltWrap : M_rangeElements )//this->rangeElements() )
        {
            auto const& elt = unwrap_ref( eltWrap );
            if ( elt.isGhostCell() )
                continue;
            M_rangeMeshElementsIdsPartialSupport.insert( elt.id() );
            M_orderedElements.push_back( std::cref(elt) );
        }
        return;
    }
    // prepare data to send with mpi
    std::map< rank_type, std::vector<size_type> > dataToSend;
    std::map< rank_type, std::vector<size_type> > dataToRecv;
    for ( auto const& eltWrap : M_rangeElements )//this->rangeElements() )
    {
        auto const& elt = unwrap_ref( eltWrap );
        if ( elt.isGhostCell() )
            continue;
        M_rangeMeshElementsIdsPartialSupport.insert( elt.id() );
        M_orderedElements.push_back( std::cref(elt) );

        auto const& idInOtherPart = elt.idInOthersPartitions();
        for ( auto const& idData : idInOtherPart )
            dataToSend[idData.first].push_back(idData.second);
    }
    // mpi comm
    int neighborSubdomains = M_mesh->neighborSubdomains().size();
    int nbMaxRequest = 2*neighborSubdomains;
    wc(this)->print( fmt::format( "[updateParallelDataPartialSupport - {}] nbMaxRequest={}, neighborSubdomains={}", rank(M_mesh), nbMaxRequest, neighborSubdomains ), FLAGS_v > 1, FLAGS_v > 0, FLAGS_v >1  );
    std::vector<mpi::request> reqs( nbMaxRequest );
    int countRequest = 0;
    std::map<rank_type,std::size_t> sizeRecv;
    std::map<rank_type,std::size_t> sizeSend;

    // get size of data to transfer
    for ( rank_type neighborRank : M_mesh->neighborSubdomains() )
    {
        sizeSend[neighborRank] = dataToSend[neighborRank].size();
        reqs[countRequest++] = M_mesh->worldComm().localComm().isend( neighborRank, 0, sizeSend[neighborRank] );
        reqs[countRequest++] = M_mesh->worldComm().localComm().irecv( neighborRank, 0, sizeRecv[neighborRank] );
    }
    // wait all requests
    mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
    countRequest = 0;

    // send/recv data
    for ( rank_type neighborRank : M_mesh->neighborSubdomains() )
    {
        std::size_t nSendData = dataToSend[neighborRank].size();
        if ( nSendData > 0 )
            reqs[countRequest++] = M_mesh->worldComm().localComm().isend( neighborRank , 0, dataToSend[neighborRank].data(), nSendData );
        std::size_t nRecvData = sizeRecv[neighborRank];
        dataToRecv[neighborRank].resize( nRecvData );
        if ( nRecvData > 0 )
            reqs[countRequest++] = M_mesh->worldComm().localComm().irecv( neighborRank , 0, dataToRecv[neighborRank].data(), nRecvData );
    }
    // wait all requests
    mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
    countRequest = 0;

    // update ghost elements
    for ( auto const& [rankRecv,eltIds] : dataToRecv )
    {
        for ( size_type ghostEltId : eltIds )
        {
            auto [itEltId,isInsert] = M_rangeMeshElementsIdsPartialSupport.insert( ghostEltId );
            if ( !isInsert )
                continue;
            auto const& ghostElt = M_mesh->element( ghostEltId );
            M_rangeMeshElementsIdsPartialSupport.insert( ghostElt.id() );
            M_orderedElements.push_back( std::cref(ghostElt) );
        }
    }
}


template <typename MeshType>
template <entity_process_t EPT>
typename MeshSupport<MeshType>::range_elements_type
MeshSupport<MeshType>::rangeElementsProcessId( rank_type part ) const
{
    if ( M_isFullSupport )
        return elements( M_mesh, part, EPT );

    range_elements_type rangeElts( M_mesh );
    for ( auto const& eltWrap : M_orderedElements )
    {
        auto const& elt = eltWrap.get();
        if constexpr ( EPT == entity_process_t::LOCAL_ONLY || EPT == entity_process_t::GHOST_ONLY || EPT == entity_process_t::LOCAL_AND_INTERPROCESS_ONLY )
        {
            if ( !Feel::detail::checkPartitionPredicate<EPT>( elt, part ) )
                continue;
        }
        rangeElts.push_back( boost::cref( elt ) );
    }
    rangeElts.shrink_to_fit();
    rangeElts.setMeshSupport( const_cast<MeshSupport*>( this )->shared_from_this() );
    return rangeElts;
}



template <typename MeshType>
template <entity_process_t EPT>
typename MeshSupport<MeshType>::range_elements_type
MeshSupport<MeshType>::rangeElementsMarkerByType( uint16_type markerType, std::set<flag_type> const& markerFlags, rank_type part ) const
{
    if ( M_isFullSupport )
        return markedelementsByType( M_mesh, markerType, markerFlags, part, EPT );

    range_elements_type rangeElts( M_mesh );
    for ( auto const& eltWrap : M_orderedElements )
    {
        auto const& elt = eltWrap.get();
        if ( !elt.hasMarkerType( markerType ) )
            continue;
        if ( !elt.marker( markerType ).hasOneOf( markerFlags ) )
            continue;
        if constexpr ( EPT == entity_process_t::LOCAL_ONLY || EPT == entity_process_t::GHOST_ONLY || EPT == entity_process_t::LOCAL_AND_INTERPROCESS_ONLY )
        {
            if ( !Feel::detail::checkPartitionPredicate<EPT>( elt, part ) )
                continue;
        }
        rangeElts.push_back( boost::cref( elt ) );
    }
    rangeElts.shrink_to_fit();
    rangeElts.setMeshSupport( const_cast<MeshSupport*>( this )->shared_from_this() );
    return rangeElts;
}



template <typename MeshType>
template <entity_process_t EPT>
typename MeshSupport<MeshType>::range_faces_type
MeshSupport<MeshType>::rangeFacesProcessId( rank_type part ) const
{
    if ( M_isFullSupport )
        return faces( M_mesh, part, EPT );

    range_faces_type rangeFaces( M_mesh );
    for ( auto const& [faceWrap,faceAttribute] : M_orderedFaces )
    {
        auto const& face = faceWrap.get();
        if constexpr ( EPT == entity_process_t::LOCAL_ONLY || EPT == entity_process_t::GHOST_ONLY || EPT == entity_process_t::LOCAL_AND_INTERPROCESS_ONLY )
        {
            if ( !Feel::detail::checkPartitionPredicate<EPT>( face, part ) )
                continue;
        }
        rangeFaces.push_back( boost::cref( face ) );
    }
    rangeFaces.shrink_to_fit();
    rangeFaces.setMeshSupport( const_cast<MeshSupport*>( this )->shared_from_this() );
    return rangeFaces;
}

template <typename MeshType>
template <entity_process_t EPT>
typename MeshSupport<MeshType>::range_faces_type
MeshSupport<MeshType>::rangeBoundaryFaces( rank_type part ) const
{
    if ( M_isFullSupport )
        return boundaryfaces( M_mesh, part, EPT );

    range_faces_type rangeFaces( M_mesh );
    for ( auto const& [faceWrap,faceAttribute] : M_orderedFaces )
    {
        auto const& face = faceWrap.get();
        if ( faceAttribute != _face_attributes::on_boundary )
            continue;
        if constexpr ( EPT == entity_process_t::LOCAL_ONLY || EPT == entity_process_t::GHOST_ONLY || EPT == entity_process_t::LOCAL_AND_INTERPROCESS_ONLY )
        {
            if ( !Feel::detail::checkPartitionPredicate<EPT>( face, part ) )
                continue;
        }
        rangeFaces.push_back( boost::cref( face ) );
    }
    rangeFaces.shrink_to_fit();
    rangeFaces.setMeshSupport( const_cast<MeshSupport*>( this )->shared_from_this() );
    return rangeFaces;
}

template <typename MeshType>
template <entity_process_t EPT>
typename MeshSupport<MeshType>::range_faces_type
MeshSupport<MeshType>::rangeInternalFaces( rank_type part ) const
{
    if ( M_isFullSupport )
        return internalfaces( M_mesh, part, EPT );

    range_faces_type rangeFaces( M_mesh );
    for ( auto const& [faceWrap,faceAttribute] : M_orderedFaces )
    {
        auto const& face = faceWrap.get();
        if ( faceAttribute != _face_attributes::intraprocess || faceAttribute != _face_attributes::interprocess  )
            continue;
        if constexpr ( EPT == entity_process_t::LOCAL_ONLY || EPT == entity_process_t::GHOST_ONLY || EPT == entity_process_t::LOCAL_AND_INTERPROCESS_ONLY )
        {
            if ( !Feel::detail::checkPartitionPredicate<EPT>( face, part ) )
                continue;
        }
        rangeFaces.push_back( boost::cref( face ) );
    }
    rangeFaces.shrink_to_fit();
    rangeFaces.setMeshSupport( const_cast<MeshSupport*>( this )->shared_from_this() );
    return rangeFaces;
}

template <typename MeshType>
template <entity_process_t EPT>
typename MeshSupport<MeshType>::range_faces_type
MeshSupport<MeshType>::rangeFacesMarkerByType( uint16_type markerType, std::set<flag_type> const& markerFlags, rank_type part ) const
{
    if ( M_isFullSupport )
        return markedfacesByType( M_mesh, markerType, markerFlags, part, EPT );

    range_faces_type rangeFaces( M_mesh );
    for ( auto const& [faceWrap,faceAttribute] : M_orderedFaces )
    {
        auto const& face = faceWrap.get();
        if ( !face.hasMarkerType( markerType ) )
            continue;
        if ( !face.marker( markerType ).hasOneOf( markerFlags ) )
            continue;
        if constexpr ( EPT == entity_process_t::LOCAL_ONLY || EPT == entity_process_t::GHOST_ONLY || EPT == entity_process_t::LOCAL_AND_INTERPROCESS_ONLY )
        {
            if ( !Feel::detail::checkPartitionPredicate<EPT>( face, part ) )
                continue;
        }
        rangeFaces.push_back( boost::cref( face ) );
    }
    rangeFaces.shrink_to_fit();
    rangeFaces.setMeshSupport( const_cast<MeshSupport*>( this )->shared_from_this() );
    return rangeFaces;
}


template <typename MeshType>
typename MeshSupport<MeshType>::range_faces_type
MeshSupport<MeshType>::rangeInterProcessFaces( rank_type part, rank_type neighbor_pid ) const
{
    CHECK( part == rank( M_mesh ) ) << "TODO ; currently we support only interprocess in current rank";

    if ( M_isFullSupport )
        return interprocessfaces( M_mesh, neighbor_pid );

    bool allNeighbor = ( neighbor_pid == invalid_v<rank_type> );
    range_faces_type rangeFaces( M_mesh );
    for ( auto const& [faceWrap,faceAttribute] : M_orderedFaces )
    {
        auto const& face = faceWrap.get();
        if ( faceAttribute != _face_attributes::interprocess  )
            continue;
        if ( face.partition1() != part )
            continue;
        if ( !allNeighbor && face.partition2() != neighbor_pid )
            continue;
        rangeFaces.push_back( boost::cref( face ) );
    }
    rangeFaces.shrink_to_fit();
    rangeFaces.setMeshSupport( const_cast<MeshSupport*>( this )->shared_from_this() );
    return rangeFaces;

}


template<typename MeshType>
template <entity_filter_t FF, entity_process_t EPT, typename ... Ts>
typename MeshSupport<MeshType>::range_elements_type
MeshSupport<MeshType>::elementsFilterImpl( Ts&&... ts ) const
{
    if constexpr ( FF == entity_filter_t::PROCESS_ID )
        return this->rangeElementsProcessId<EPT>( std::forward<Ts>( ts )... );
    else if constexpr ( FF == entity_filter_t::MARKER )
        return this->rangeElementsMarkerByType<EPT>( std::forward<Ts>( ts )... );
    // else if constexpr ( FF == entity_filter_t::ON_BOUNDARY )
    //     return this->rangeBoundaryFaces<EPT>( std::forward<Ts>( ts )... );
    // else if constexpr ( FF == entity_filter_t::INTERNAL )
    //     return this->internalFaces<EPT>( std::forward<Ts>( ts )... );
    CHECK( false ) << "TODO";
    return {};
}

template<typename MeshType>
template <entity_filter_t FF, typename ... Ts>
typename MeshSupport<MeshType>::range_elements_type
MeshSupport<MeshType>::elementsFilter( entity_process_t ept, Ts&&... ts ) const
{
    return std::invoke(
        [this,&ept](auto&& ... args)
            {
                switch ( ept )
                {
                default:
                case entity_process_t::LOCAL_ONLY:
                    return this->elementsFilterImpl<FF, entity_process_t::LOCAL_ONLY>( std::forward<decltype(args)>(args) ... );
                case entity_process_t::LOCAL_AND_INTERPROCESS_ONLY:
                    return this->elementsFilterImpl<FF, entity_process_t::LOCAL_AND_INTERPROCESS_ONLY>( std::forward<decltype(args)>(args) ... );
                case entity_process_t::GHOST_ONLY:
                    return this->elementsFilterImpl<FF, entity_process_t::GHOST_ONLY>( std::forward<decltype(args)>(args) ... );
                case entity_process_t::ALL:
                    return this->elementsFilterImpl<FF, entity_process_t::ALL>( std::forward<decltype(args)>(args) ... );
                }
            },
        std::forward<Ts>( ts )... );
}

template<typename MeshType>
template <entity_filter_t FF, entity_process_t EPT, typename ... Ts>
typename MeshSupport<MeshType>::range_faces_type
MeshSupport<MeshType>::facesFilterImpl( Ts&&... ts ) const
{
    if constexpr ( FF == entity_filter_t::PROCESS_ID )
        return this->rangeFacesProcessId<EPT>( std::forward<Ts>( ts )... );
    else if constexpr ( FF == entity_filter_t::MARKER )
        return this->rangeFacesMarkerByType<EPT>( std::forward<Ts>( ts )... );
    else if constexpr ( FF == entity_filter_t::ON_BOUNDARY )
        return this->rangeBoundaryFaces<EPT>( std::forward<Ts>( ts )... );
    else if constexpr ( FF == entity_filter_t::INTERNAL )
        return this->rangeInternalFaces<EPT>( std::forward<Ts>( ts )... );
    return {};
}

template<typename MeshType>
template <entity_filter_t FF, typename ... Ts>
typename MeshSupport<MeshType>::range_faces_type
MeshSupport<MeshType>::facesFilter( entity_process_t ept, Ts&&... ts ) const
{
    return std::invoke(
        [this,&ept](auto&& ... args)
            {
                switch ( ept )
                {
                default:
                case entity_process_t::LOCAL_ONLY:
                    return this->facesFilterImpl<FF, entity_process_t::LOCAL_ONLY>( std::forward<decltype(args)>(args) ... );
                case entity_process_t::LOCAL_AND_INTERPROCESS_ONLY:
                    return this->facesFilterImpl<FF, entity_process_t::LOCAL_AND_INTERPROCESS_ONLY>( std::forward<decltype(args)>(args) ... );
                case entity_process_t::GHOST_ONLY:
                    return this->facesFilterImpl<FF, entity_process_t::GHOST_ONLY>( std::forward<decltype(args)>(args) ... );
                case entity_process_t::ALL:
                    return this->facesFilterImpl<FF, entity_process_t::ALL>( std::forward<decltype(args)>(args) ... );
                }
            },
        std::forward<Ts>( ts )... );
}




template<typename MeshSupportType, std::enable_if_t<std::is_base_of_v<MeshSupportBase,unwrap_ptr_t<MeshSupportType>>,int> = 0>
using support_mesh_t = typename unwrap_ptr_t<MeshSupportType>::mesh_type;

template<typename MeshSupportType, std::enable_if_t<std::is_base_of_v<MeshSupportBase,unwrap_ptr_t<MeshSupportType>>,int> = 0>
auto
elements( MeshSupportType const& imesh, entity_process_t ept = entity_process_t::LOCAL_ONLY )
{
    return imesh->template elementsFilter<entity_filter_t::PROCESS_ID>( ept, rank( imesh->mesh() ) );
}
template<typename MeshSupportType, std::enable_if_t<std::is_base_of_v<MeshSupportBase,unwrap_ptr_t<MeshSupportType>>,int> = 0>
auto
markedelements( MeshSupportType const& imesh, boost::any markersFlag, entity_process_t ept = entity_process_t::LOCAL_ONLY )
{
    std::set<flag_type> markerFlagSet = imesh->mesh()->markersId( markersFlag );
    return imesh->template elementsFilter<entity_filter_t::MARKER>( ept, 1, markerFlagSet, rank( imesh->mesh() ) );
}
template<typename MeshSupportType, std::enable_if_t<std::is_base_of_v<MeshSupportBase,unwrap_ptr_t<MeshSupportType>>,int> = 0>
auto
faces( MeshSupportType const& imesh, entity_process_t ept = entity_process_t::LOCAL_ONLY  )
{
    return imesh->template facesFilter<entity_filter_t::PROCESS_ID>( ept, rank( imesh->mesh() ) );
}
template<typename MeshSupportType, std::enable_if_t<std::is_base_of_v<MeshSupportBase,unwrap_ptr_t<MeshSupportType>>,int> = 0>
auto
boundaryfaces( MeshSupportType const& imesh, entity_process_t ept = entity_process_t::LOCAL_ONLY )
{
    return imesh->template facesFilter<entity_filter_t::ON_BOUNDARY>( ept, rank( imesh->mesh() ) );
}
template<typename MeshSupportType, std::enable_if_t<std::is_base_of_v<MeshSupportBase,unwrap_ptr_t<MeshSupportType>>,int> = 0>
auto
internalfaces( MeshSupportType const& imesh, entity_process_t ept = entity_process_t::LOCAL_ONLY )
{
    return imesh->template facesFilter<entity_filter_t::INTERNAL>( ept, rank( imesh->mesh() ) );
}

template<typename MeshSupportType, std::enable_if_t<std::is_base_of_v<MeshSupportBase,unwrap_ptr_t<MeshSupportType>>,int> = 0>
auto
markedfaces( MeshSupportType const& imesh, boost::any markersFlag, entity_process_t ept = entity_process_t::LOCAL_ONLY )
{
    std::set<flag_type> markerFlagSet = imesh->mesh()->markersId( markersFlag );
    return imesh->template facesFilter<entity_filter_t::MARKER>( ept, 1, markerFlagSet, rank( imesh->mesh() ) );
}

template<typename MeshSupportType, std::enable_if_t<std::is_base_of_v<MeshSupportBase,unwrap_ptr_t<MeshSupportType>>,int> = 0>
auto
interprocessfaces( MeshSupportType const& imesh, rank_type neighbor_pid = invalid_v<rank_type> )
{
    return imesh->rangeInterProcessFaces( rank( imesh->mesh() ), neighbor_pid );
}


} // namespace Feel

#endif
