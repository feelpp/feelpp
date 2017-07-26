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

#include <feel/feelmesh/traits.hpp>
#include <feel/feelmesh/filters.hpp>

namespace Feel
{

/**
 * \brief Description of a mesh support.
 * A function space can use this object for built a space in a part of the mesh
 */
template<typename MeshType>
class MeshSupport
{
public :
    using mesh_type = typename MeshTraits<MeshType>::mesh_type;
    using mesh_ptrtype = boost::shared_ptr<mesh_type>;
    using range_elements_type = elements_reference_wrapper_t<mesh_type>;
    using range_faces_type = faces_reference_wrapper_t<mesh_type>;
    using element_type = typename mesh_type::element_type;
    using face_type = typename mesh_type::face_type;

    MeshSupport( mesh_ptrtype const& mesh )
        :
        M_mesh( mesh ),
        M_rangeElements( elements(mesh) ),
        M_isFullSupport( true ),
        M_hasUpdatedParallelData( false ),
        M_hasUpdatedBoundaryInternalFaces( false )
        {}

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
        }

    bool isFullSupport() const { return M_isFullSupport; }
    bool isPartialSupport() const { return !M_isFullSupport; }

    range_elements_type const& rangeElements() const { return M_rangeElements; }
    range_faces_type const& rangeInterProcessFaces() const { return M_rangeInterProcessFaces; }
    range_faces_type const& rangeBoundaryFaces() const { return M_rangeBoundaryFaces; }
    range_faces_type const& rangeInternalFaces() const { return M_rangeInternalFaces; }

    range_elements_type rangeElements( EntityProcessType entity ) const
        {
            if ( M_isFullSupport )
                return elements(M_mesh,entity );

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
                    auto const& faceip = boost::unwrap_ref( faceWrap );//*face_it );
                    auto const& elt0 = faceip.element0();
                    auto const& elt1 = faceip.element1();
                    const bool elt0isGhost = elt0.isGhostCell();
                    auto const& eltOffProc = (elt0isGhost)?elt0:elt1;
                    auto const& eltOnProc = (elt0isGhost)?elt1:elt0;
                    if ( eltGhostDone.find( eltOffProc.id() ) != eltGhostDone.end() )
                        continue;
                    myExtendedElements->push_back( boost::cref( eltOffProc ) );
                    eltGhostDone.insert( eltOffProc.id() );
                }
            }
            range_elements_type rangeExtendedElements = boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(), myExtendedElements->begin(),myExtendedElements->end(),myExtendedElements );
            return rangeExtendedElements;
        }

    size_type numElements() const
        {
            if ( M_isFullSupport )
                return M_mesh->numElements();
            else
                return M_rangeMeshElementsIdsPartialSupport.size();
        }
    bool hasElement( size_type eltId ) const
        {
            if ( M_isFullSupport )
                return M_mesh->hasElement( eltId );
            else
                return M_rangeMeshElementsIdsPartialSupport.find( eltId ) != M_rangeMeshElementsIdsPartialSupport.end();
        }
    bool hasGhostElement( size_type eltId ) const
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

    bool isGhostFace( face_type const& face ) const
        {
            if ( !face.isGhostFace() )
                return false;

            if ( M_isFullSupport )
                return true;
            else
            {
                bool hasElt0 = this->hasElement( face.element(0).id() );
                bool hasElt1 = this->hasElement( face.element(1).id() );
                return ( hasElt0 && hasElt1 );
            }
        }

    std::unordered_set<size_type> const& rangeMeshElementsIdsPartialSupport() const { return M_rangeMeshElementsIdsPartialSupport; }
    std::unordered_set<size_type> const& rangeMeshElementsGhostIdsPartialSupport() const { return M_rangeMeshElementsGhostIdsPartialSupport; }

    void updateParallelData()
        {
            if ( M_hasUpdatedParallelData )
                return;

            if ( M_isFullSupport )
                this->updateParallelDataFullSupport();
            else
                this->updateParallelDataPartialSupport();

            M_hasUpdatedParallelData = true;
        }
    void updateBoundaryInternalFaces()
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
    void updateParallelDataFullSupport()
        {
            M_rangeInterProcessFaces = interprocessfaces(M_mesh);
        }
    void updateParallelDataPartialSupport()
        {
            typename MeshTraits<mesh_type>::faces_reference_wrapper_ptrtype myipfaces( new typename MeshTraits<mesh_type>::faces_reference_wrapper_type );
            if ( M_mesh->worldComm().localSize() == 1 )
            {
                M_rangeInterProcessFaces = boost::make_tuple( mpl::size_t<MESH_FACES>(), myipfaces->begin(),myipfaces->end(),myipfaces );
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
            for ( rank_type neighborRank : M_mesh->neighborSubdomains() )
            {
                reqs[cptRequest++] = M_mesh->worldComm().localComm().isend( neighborRank , 0, dataToSend[neighborRank] );
                reqs[cptRequest++] = M_mesh->worldComm().localComm().irecv( neighborRank , 0, dataToRecv[neighborRank] );
            }
            mpi::wait_all(reqs, reqs + nbRequest);
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
            M_rangeInterProcessFaces = boost::make_tuple( mpl::size_t<MESH_FACES>(), myipfaces->begin(),myipfaces->end(),myipfaces );
        }

    void updateBoundaryInternalFacesFullSupport()
        {
            M_rangeBoundaryFaces = boundaryfaces(M_mesh);
            M_rangeInternalFaces = internalfaces(M_mesh);
        }
    void updateBoundaryInternalFacesPartialSupport()
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

            typename MeshTraits<mesh_type>::faces_reference_wrapper_ptrtype mybfaces( new typename MeshTraits<mesh_type>::faces_reference_wrapper_type );
            typename MeshTraits<mesh_type>::faces_reference_wrapper_ptrtype myifaces( new typename MeshTraits<mesh_type>::faces_reference_wrapper_type );
            for ( auto const& faceDataPair : faceInRange )
            {
                auto const& faceData = faceDataPair.second;
                if ( faceData.second == 1 )
                    mybfaces->push_back( boost::cref( *faceData.first ) );
                else
                    myifaces->push_back( boost::cref( *faceData.first ) );
            }
            M_rangeBoundaryFaces = boost::make_tuple( mpl::size_t<MESH_FACES>(),mybfaces->begin(),mybfaces->end(),mybfaces );
            M_rangeInternalFaces = boost::make_tuple( mpl::size_t<MESH_FACES>(),myifaces->begin(),myifaces->end(),myifaces );
        }
private :
    mesh_ptrtype M_mesh;
    range_elements_type M_rangeElements;
    range_faces_type M_rangeInterProcessFaces;
    range_faces_type M_rangeBoundaryFaces;
    range_faces_type M_rangeInternalFaces;
    std::unordered_set<size_type> M_rangeMeshElementsIdsPartialSupport;
    std::unordered_set<size_type> M_rangeMeshElementsGhostIdsPartialSupport;

    bool M_isFullSupport;
    bool M_hasUpdatedParallelData;
    bool M_hasUpdatedBoundaryInternalFaces;

};

} // namespace Feel

#endif
