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
class MeshSupport : public MeshSupportBase
{
public :
    using super_type = MeshSupportBase;
    using mesh_type = typename MeshTraits<MeshType>::mesh_type;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;
    using range_elements_type = Range<mesh_type,MESH_ELEMENTS>;
    using range_faces_type = Range<mesh_type,MESH_FACES>;
    using element_type = typename mesh_type::element_type;
    using face_type = typename mesh_type::face_type;
    static constexpr int nDim = mesh_type::nDim;

    MeshSupport() = default;
    MeshSupport( mesh_ptrtype const& mesh )
        :
        M_mesh( mesh ),
        M_rangeElements( elements(mesh) ),
        M_isFullSupport( true ),
        M_hasUpdatedParallelData( false ),
        M_hasUpdatedBoundaryInternalFaces( false )
        {
        }

    MeshSupport( mesh_ptrtype const& mesh, range_elements_type const& rangeElements )
        :
        M_mesh( mesh ),
        M_rangeElements( rangeElements ),
        M_isFullSupport( false ),
        M_hasUpdatedParallelData( false ),
        M_hasUpdatedBoundaryInternalFaces( false )
        {
            for (auto const& eltWrap : M_rangeElements )
                M_rangeMeshElementsIdsPartialSupport.insert( unwrap_ref(eltWrap).id() );

            M_localizationToolPartialSupport = std::make_shared<Localization<mesh_type>>();
            M_localizationToolPartialSupport->setMesh( M_mesh, M_rangeElements, false );
            //M_localizationToolPartialSupport->/*init*/reset( M_rangeElements );

            M_mesh->attachMeshSupport( this );
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

    range_elements_type const& rangeElements() const { return M_rangeElements; }
    range_faces_type const& rangeInterProcessFaces() const { this->updateParallelData(); return M_rangeInterProcessFaces; }
    range_faces_type const& rangeBoundaryFaces() const { this->updateBoundaryInternalFaces();return M_rangeBoundaryFaces; }
    range_faces_type const& rangeInternalFaces() const { this->updateBoundaryInternalFaces();return M_rangeInternalFaces; }

    range_elements_type rangeElements( EntityProcessType entity ) const;
    //!
    //! return the set of elements of marker type marker_t with marker flag 
    //!
    range_elements_type rangeMarkedElements( uint16_type marker_t, boost::any flag );
    //!
    //! return the set of faces of marker type marker_t with marker flag 
    //!
    range_faces_type rangeMarkedFaces( uint16_type marker_t, boost::any flag );

    //!
    //! @return true if some markers in \p l are present in the mesh data structure, false otherwise
    //!
    bool hasAnyMarker( std::initializer_list<std::string> l )
        {
            if ( M_isFullSupport )
                return M_mesh->hasAnyMarker( l );
            for (auto n : l )
            {
                if ( nelements(rangeMarkedFaces(1,n), true ) )
                    return true;
            }
            return false;
        }
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
    bool hasGhostElement( size_type eltId ) const override
        {
            if ( M_isFullSupport )
            {
                if ( M_mesh->hasElement( eltId ) )
                    return M_mesh->element( eltId ).isGhostCell();
                else
                    return false;
            }
            else
                return M_rangeMeshElementsGhostIdsPartialSupport.find( eltId ) != M_rangeMeshElementsGhostIdsPartialSupport.end();
        }

    template <typename FaceType>
    bool isGhostFace( FaceType const& face,
                      typename std::enable_if_t<std::is_same_v<FaceType,face_type> >* = nullptr ) const
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

    template <typename FaceType>
    bool isGhostFace( FaceType const& face,
                      typename std::enable_if_t<!std::is_same_v<FaceType,face_type> >* = nullptr ) const
        {
            return false;
        }

    std::unordered_set<size_type> const& rangeMeshElementsIdsPartialSupport() const override { return M_rangeMeshElementsIdsPartialSupport; }
    std::unordered_set<size_type> const& rangeMeshElementsGhostIdsPartialSupport() const override { return M_rangeMeshElementsGhostIdsPartialSupport; }

    void updateParallelData() const
        {
            if ( M_hasUpdatedParallelData )
                return;

            if ( M_isFullSupport )
                this->updateParallelDataFullSupport();
            else
                this->updateParallelDataPartialSupport();

            M_hasUpdatedParallelData = true;
        }
    void updateBoundaryInternalFaces() const
        {
            if ( M_hasUpdatedBoundaryInternalFaces )
                return;

            if ( M_isFullSupport )
                this->updateBoundaryInternalFacesFullSupport();
            else
                this->updateBoundaryInternalFacesPartialSupport();

            M_hasUpdatedBoundaryInternalFaces = true;
        }
private :
    void updateParallelDataFullSupport() const
        {
            M_rangeInterProcessFaces = interprocessfaces(M_mesh);
        }
    void updateParallelDataPartialSupport() const
        {
            typename MeshTraits<mesh_type>::faces_reference_wrapper_ptrtype myipfaces( new typename MeshTraits<mesh_type>::faces_reference_wrapper_type );
            if ( M_mesh->worldComm().localSize() == 1 )
            {
                M_rangeInterProcessFaces = range( _range=boost::make_tuple( mpl::size_t<MESH_FACES>(), myipfaces->begin(),myipfaces->end(),myipfaces ), _mesh=M_mesh );
                return;
            }
            // prepare data to send with mpi
            std::map< rank_type, std::vector<size_type> > dataToSend;
            std::map< rank_type, std::vector<size_type> > dataToRecv;
            for ( auto const& eltWrap : this->rangeElements() )
            {
                auto const& elt = unwrap_ref( eltWrap );
                if ( elt.isGhostCell() )
                    continue;
                auto const& idInOtherPart = elt.idInOthersPartitions();
                for ( auto const& idData : idInOtherPart )
                    dataToSend[idData.first].push_back(idData.second);
            }
            // mpi comm
            int neighborSubdomains = M_mesh->neighborSubdomains().size();
            int nbRequest = 2*neighborSubdomains;
            mpi::request * reqs = new mpi::request[nbRequest];
            int cptRequest=0;
            std::map<rank_type,size_type> sizeRecv;

            // get size of data to transfer
            for ( rank_type neighborRank : M_mesh->neighborSubdomains() )
            {
                reqs[cptRequest++] = M_mesh->worldComm().localComm().isend( neighborRank, 0, (size_type)dataToSend[neighborRank].size() );
                reqs[cptRequest++] = M_mesh->worldComm().localComm().irecv( neighborRank, 0, sizeRecv[neighborRank] );
            }
            // wait all requests
            mpi::wait_all(reqs, reqs + cptRequest);

            // send/recv data
            cptRequest=0;
            for ( rank_type neighborRank : M_mesh->neighborSubdomains() )
            {
                int nSendData = dataToSend[neighborRank].size();
                if ( nSendData > 0 )
                    reqs[cptRequest++] = M_mesh->worldComm().localComm().isend( neighborRank , 0, &(dataToSend[neighborRank][0]), nSendData );
                int nRecvData = sizeRecv[neighborRank];
                dataToRecv[neighborRank].resize( nRecvData );
                if ( nRecvData > 0 )
                    reqs[cptRequest++] = M_mesh->worldComm().localComm().irecv( neighborRank , 0, &(dataToRecv[neighborRank][0]), nRecvData );
            }
            mpi::wait_all(reqs, reqs + cptRequest);
            delete [] reqs;
            // get elt ids of ghost in mesh
            std::unordered_set<size_type> ghostEltIdsInMesh;
            auto rangeGhostElt = M_mesh->ghostElements();
            auto itghost = std::get<0>( rangeGhostElt );
            auto enghost = std::get<1>( rangeGhostElt );
            for ( ; itghost != enghost ; ++itghost )
                ghostEltIdsInMesh.insert( unwrap_ref( *itghost ).id() );
            // get elt ids of ghost in range of mesh element
            for ( auto const& dataToRecvByProc : dataToRecv )
            {
                for ( size_type eltId : dataToRecvByProc.second )
                {
                    if ( ghostEltIdsInMesh.find( eltId ) != ghostEltIdsInMesh.end() )
                        M_rangeMeshElementsGhostIdsPartialSupport.insert( eltId );
                }
            }

            for ( size_type eltId : M_rangeMeshElementsGhostIdsPartialSupport )
                M_rangeMeshElementsIdsPartialSupport.insert( eltId );

            for ( auto const& eltWrap : this->rangeElements() )
            {
                auto const& elt = unwrap_ref( eltWrap );
                if ( elt.isGhostCell() )
                    continue;
                for ( uint16_type f = 0; f < element_type::numTopologicalFaces; ++f )
                {
                    if ( !elt.facePtr(f) )
                        continue;
                    auto const& face = elt.face( f );
                    if ( !face.isInterProcessDomain() )
                        continue;
                    auto const& elt0 = face.element0();
                    auto const& elt1 = face.element1();
                    const bool elt0isGhost = elt0.isGhostCell();
                    auto const& eltOnProc = (elt0isGhost)?elt1:elt0;
                    auto const& eltOffProc = (elt0isGhost)?elt0:elt1;
                    if ( M_rangeMeshElementsGhostIdsPartialSupport.find( eltOffProc.id() ) == M_rangeMeshElementsGhostIdsPartialSupport.end() )
                        continue;
                    myipfaces->push_back( boost::cref( face ) );
                }
            }
            M_rangeInterProcessFaces = range( _range=boost::make_tuple( mpl::size_t<MESH_FACES>(), myipfaces->begin(),myipfaces->end(),myipfaces ), _mesh=M_mesh );
        }

    void updateBoundaryInternalFacesFullSupport() const
        {
            M_rangeBoundaryFaces = boundaryfaces(M_mesh);
            M_rangeInternalFaces = internalfaces(M_mesh);
        }
    void updateBoundaryInternalFacesPartialSupport() const
        {
            this->updateParallelData();

            std::unordered_map<size_type,std::pair<const face_type*,uint8_type>> faceInRange;
            for ( auto const& eltWrap : this->rangeElements() )
            {
                auto const& elt = unwrap_ref( eltWrap );
                for ( uint16_type i = 0; i < mesh_type::element_type::numTopologicalFaces; ++i )
                {
                    if ( !elt.facePtr(i) )
                        continue;
                    const face_type* facePtr = elt.facePtr(i);
                    size_type faceId = facePtr->id();
                    auto const& face = elt.face(i);
                    if ( faceInRange.find( faceId ) != faceInRange.end() )
                        faceInRange[faceId].second = 2;
                    else
                        faceInRange[faceId] = std::make_pair(facePtr,1);
                }
            }
            for ( auto const& faceWrap : this->rangeInterProcessFaces() )
            {
                size_type faceId = unwrap_ref(faceWrap).id();
                DCHECK( faceInRange.find( faceId ) != faceInRange.end() ) << "something wrong";
                faceInRange[faceId].second = 3;
            }

            std::map<rank_type,std::vector<size_type> > dataToSend;
            std::map<rank_type,std::vector<size_type> > dataToRecv;
            typename MeshTraits<mesh_type>::faces_reference_wrapper_ptrtype mybfaces( new typename MeshTraits<mesh_type>::faces_reference_wrapper_type );
            typename MeshTraits<mesh_type>::faces_reference_wrapper_ptrtype myifaces( new typename MeshTraits<mesh_type>::faces_reference_wrapper_type );
            for ( auto const& faceDataPair : faceInRange )
            {
                auto const& faceData = faceDataPair.second;
                if ( faceData.second == 1 )
                {
                    auto const& theface = *faceData.first;
                    if ( theface.isInterProcessDomain() )
                    {
                        rank_type neighborPid = theface.partition2();
                        dataToSend[neighborPid].push_back( theface.idInOthersPartitions(neighborPid) );
                    }
                    mybfaces->push_back( boost::cref( theface ) );
                }
                else
                    myifaces->push_back( boost::cref( *faceData.first ) );
            }

            // maybe some boundary faces on interprocess faces are not detected
            // on neighbor part (because not connected to an element of partial support)
            // but should be added on range : required mpi comm
            int neighborSubdomains = M_mesh->neighborSubdomains().size();
            int nbRequest = 2*neighborSubdomains;
            mpi::request * reqs = new mpi::request[nbRequest];
            int cptRequest=0;
            std::map<rank_type,size_type> sizeRecv;

            // get size of data to transfer
            for ( rank_type neighborRank : M_mesh->neighborSubdomains() )
            {
                reqs[cptRequest++] = M_mesh->worldComm().localComm().isend( neighborRank, 0, (size_type)dataToSend[neighborRank].size() );
                reqs[cptRequest++] = M_mesh->worldComm().localComm().irecv( neighborRank, 0, sizeRecv[neighborRank] );
            }
            // wait all requests
            mpi::wait_all(reqs, reqs + cptRequest);

            // send/recv data
            cptRequest=0;
            for ( rank_type neighborRank : M_mesh->neighborSubdomains() )
            {
                int nSendData = dataToSend[neighborRank].size();
                if ( nSendData > 0 )
                    reqs[cptRequest++] = M_mesh->worldComm().localComm().isend( neighborRank , 0, &(dataToSend[neighborRank][0]), nSendData );
                int nRecvData = sizeRecv[neighborRank];
                dataToRecv[neighborRank].resize( nRecvData );
                if ( nRecvData > 0 )
                    reqs[cptRequest++] = M_mesh->worldComm().localComm().irecv( neighborRank , 0, &(dataToRecv[neighborRank][0]), nRecvData );
            }
            mpi::wait_all(reqs, reqs + cptRequest);
            delete [] reqs;

            for ( auto const& dataRecvByProc : dataToRecv )
            {
                for ( size_type faceId : dataRecvByProc.second )
                {
                    if ( faceInRange.find( faceId ) == faceInRange.end() )
                        mybfaces->push_back( boost::cref( M_mesh->face( faceId ) ) );
                }
            }


            M_rangeBoundaryFaces = range( _range=boost::make_tuple( mpl::size_t<MESH_FACES>(),mybfaces->begin(),mybfaces->end(),mybfaces ), _mesh=M_mesh );
            M_rangeInternalFaces = range( _range=boost::make_tuple( mpl::size_t<MESH_FACES>(),myifaces->begin(),myifaces->end(),myifaces ), _mesh=M_mesh );
        }

    void resetLocalizationTool() override
        {
            if ( M_localizationToolPartialSupport )
                M_localizationToolPartialSupport->reset();
        }

private :
    mesh_ptrtype M_mesh;
    range_elements_type M_rangeElements;
    std::shared_ptr<Localization<mesh_type>> M_localizationToolPartialSupport;

    mutable range_faces_type M_rangeInterProcessFaces;
    mutable range_faces_type M_rangeBoundaryFaces;
    mutable range_faces_type M_rangeInternalFaces;
    mutable std::unordered_set<size_type> M_rangeMeshElementsIdsPartialSupport;
    mutable std::unordered_set<size_type> M_rangeMeshElementsGhostIdsPartialSupport;

    bool M_isFullSupport;
    mutable bool M_hasUpdatedParallelData;
    mutable bool M_hasUpdatedBoundaryInternalFaces;

};

template <typename MeshType>
typename MeshSupport<MeshType>::range_elements_type
MeshSupport<MeshType>::rangeElements( EntityProcessType entity ) const
{
    if ( M_isFullSupport )
        return elements( M_mesh, entity );

    if ( entity == EntityProcessType::LOCAL_ONLY )
        return M_rangeElements;

    typename MeshTraits<mesh_type>::elements_reference_wrapper_ptrtype myExtendedElements( new typename MeshTraits<mesh_type>::elements_reference_wrapper_type );

    if ( entity == EntityProcessType::ALL )
    {
        for ( auto const& eltWrap : M_rangeElements )
            myExtendedElements->push_back( eltWrap );
    }

    if ( ( ( entity == EntityProcessType::GHOST_ONLY ) || ( entity == EntityProcessType::ALL ) ) && ( M_mesh->worldComm().localSize() > 1 ) )
    {
        CHECK( M_hasUpdatedParallelData ) << "parallel data must be updated";

        std::unordered_set<size_type> eltGhostDone;
        for ( auto const& faceWrap : M_rangeInterProcessFaces )
        {
            auto const& faceip = boost::unwrap_ref( faceWrap ); //*face_it );
            auto const& elt0 = faceip.element0();
            auto const& elt1 = faceip.element1();
            const bool elt0isGhost = elt0.isGhostCell();
            auto const& eltOffProc = ( elt0isGhost ) ? elt0 : elt1;
            auto const& eltOnProc = ( elt0isGhost ) ? elt1 : elt0;
            if ( eltGhostDone.find( eltOffProc.id() ) != eltGhostDone.end() )
                continue;
            myExtendedElements->push_back( boost::cref( eltOffProc ) );
            eltGhostDone.insert( eltOffProc.id() );
        }
    }
    range_elements_type rangeExtendedElements = range( _range=boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(), myExtendedElements->begin(), myExtendedElements->end(), myExtendedElements ), _mesh=M_mesh );
    return rangeExtendedElements;
}

template<typename MeshType>
typename MeshSupport<MeshType>::range_elements_type 
MeshSupport<MeshType>::rangeMarkedElements( uint16_type marker_t, boost::any flag )
{
    std::set<flag_type> markerFlagSet = Feel::unwrap_ptr( M_mesh ).markersId( flag );
    flag_type m = *markerFlagSet.begin();
    if ( M_isFullSupport )
    {
        return markedelementsByType( M_mesh, marker_t, flag );
    }

    typename MeshTraits<mesh_type>::elements_reference_wrapper_ptrtype myelements( new typename MeshTraits<mesh_type>::elements_reference_wrapper_type );
    auto insertMarkedElements = [&myelements, &marker_t,&markerFlagSet]( auto const& eltWrap)
                             {
                                 auto const& element = unwrap_ref( eltWrap );
                                 if ( !element.hasMarker( marker_t ) )
                                     return false;
                                 if ( element.marker( marker_t ).isOff() )
                                     return false;
                                 if ( markerFlagSet.find( element.marker( marker_t ).value() ) == markerFlagSet.end() )
                                     return false;

                                 myelements->push_back( boost::cref( element ) );
                                 return true;
                             };
    for ( auto const& eltWrap : this->rangeElements() )
    {
        insertMarkedElements( eltWrap );
    }
    return range( _range=boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),myelements->begin(),myelements->end(),myelements ), _mesh=M_mesh );
}

template<typename MeshType>
typename MeshSupport<MeshType>::range_faces_type 
MeshSupport<MeshType>::rangeMarkedFaces( uint16_type marker_t, boost::any flag )
{
    std::set<flag_type> markerFlagSet = Feel::unwrap_ptr( M_mesh ).markersId( flag );
    flag_type m = *markerFlagSet.begin();
    if ( M_isFullSupport )
    {
        return markedfacesByType( M_mesh, marker_t, flag );
    }
        
    typename MeshTraits<mesh_type>::faces_reference_wrapper_ptrtype myfaces( new typename MeshTraits<mesh_type>::faces_reference_wrapper_type );
    auto insertMarkedFace = [&myfaces, &marker_t,&markerFlagSet]( auto const& eltWrap)
                             {
                                 auto const& face = unwrap_ref( eltWrap );
                                 if ( !face.hasMarker( marker_t ) )
                                     return false;
                                 if ( face.marker( marker_t ).isOff() )
                                     return false;
                                 if ( markerFlagSet.find( face.marker( marker_t ).value() ) == markerFlagSet.end() )
                                     return false;
        
                                 myfaces->push_back( boost::cref( face ) );
                                 return true;
                             };
    for ( auto const& eltWrap : this->rangeBoundaryFaces() )
    {
        insertMarkedFace( eltWrap );
    }
    for ( auto const& eltWrap : this->rangeInternalFaces() )
    {
        insertMarkedFace( eltWrap );
    }
    return range( _range=boost::make_tuple( mpl::size_t<MESH_FACES>(),myfaces->begin(),myfaces->end(),myfaces ), _mesh=M_mesh );
}

template<typename MeshSupportType, std::enable_if_t<std::is_base_of_v<MeshSupportBase,unwrap_ptr_t<MeshSupportType>>,int> = 0>
using support_mesh_t = typename unwrap_ptr_t<MeshSupportType>::mesh_type;

template<typename MeshSupportType, std::enable_if_t<std::is_base_of_v<MeshSupportBase,unwrap_ptr_t<MeshSupportType>>,int> = 0>
auto
elements( MeshSupportType const& imesh )
{
    return imesh->rangeElements();
}
template<typename MeshSupportType, std::enable_if_t<std::is_base_of_v<MeshSupportBase,unwrap_ptr_t<MeshSupportType>>,int> = 0>
auto
markedelements( MeshSupportType const& imesh, boost::any flag )
{
    return imesh->rangeMarkedElements( 1, flag );
}
template<typename MeshSupportType, std::enable_if_t<std::is_base_of_v<MeshSupportBase,unwrap_ptr_t<MeshSupportType>>,int> = 0>
auto
faces( MeshSupportType const& imesh )
{
    using mesh_type = typename unwrap_ptr_t<MeshSupportType>::mesh_type;
    typename MeshTraits<mesh_type>::faces_reference_wrapper_ptrtype myfaces( new typename MeshTraits<mesh_type>::faces_reference_wrapper_type );
    for ( auto const& eltWrap : imesh->rangeBoundaryFaces() )
    {
        myfaces->push_back( eltWrap );
    }
    for ( auto const& eltWrap : imesh->rangeInternalFaces() )
    {
        myfaces->push_back( eltWrap );
    }
    return range( _range=boost::make_tuple( mpl::size_t<MESH_FACES>(),myfaces->begin(),myfaces->end(),myfaces ), _mesh=imesh->mesh() );
}
template<typename MeshSupportType, std::enable_if_t<std::is_base_of_v<MeshSupportBase,unwrap_ptr_t<MeshSupportType>>,int> = 0>
auto const&
boundaryfaces( MeshSupportType const& imesh )
{
    return imesh->rangeBoundaryFaces();
}
template<typename MeshSupportType, std::enable_if_t<std::is_base_of_v<MeshSupportBase,unwrap_ptr_t<MeshSupportType>>,int> = 0>
auto const&
internalfaces( MeshSupportType const& imesh )
{
    return imesh->rangeInternalFaces();
}


template<typename MeshSupportType, std::enable_if_t<std::is_base_of_v<MeshSupportBase,unwrap_ptr_t<MeshSupportType>>,int> = 0>
auto const&
markedfaces( MeshSupportType const& imesh, boost::any flag )
{
    return imesh->rangeMarkedFaces( 1, flag );
}

} // namespace Feel

#endif
