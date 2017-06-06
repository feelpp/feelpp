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


    MeshSupport( mesh_ptrtype const& mesh )
        :
        M_mesh( mesh ),
        M_rangeElements( elements(mesh) ),
        M_isFullSupport( true ),
        M_hasUpdatedParallelData( false ),
        M_hasUpdatedBoundaryFaces( false )
        {}

    MeshSupport( mesh_ptrtype const& mesh, range_elements_type const& rangeElements )
        :
        M_mesh( mesh ),
        M_rangeElements( rangeElements ),
        M_isFullSupport( false ),
        M_hasUpdatedParallelData( false ),
        M_hasUpdatedBoundaryFaces( false )
        {}

    bool isFullSupport() const { return M_isFullSupport; }
    bool isPartialSupport() const { return !M_isFullSupport; }

    range_elements_type const& rangeElements() const { return M_rangeElements; }
    range_faces_type const& rangeInterProcessFaces() const { return M_rangeInterProcessFaces; }
    range_faces_type const& rangeBoundaryFaces() const { return M_rangeBoundaryFaces; }

    bool hasGhostElement( size_type eltId ) const { return M_rangeMeshElementsGhostIds.find( eltId ) != M_rangeMeshElementsGhostIds.end(); }

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
    void updateBoundaryFaces()
        {
            if ( M_hasUpdatedBoundaryFaces )
                return;

            if ( M_isFullSupport )
                this->updateBoundaryFacesFullSupport();
            else
                this->updateBoundaryFacesPartialSupport();

            M_hasUpdatedBoundaryFaces = true;
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
                        M_rangeMeshElementsGhostIds.insert( eltId );
                }
            }

            for ( auto const& eltWrap : this->rangeElements() )
            {
                auto const& elt = unwrap_ref( eltWrap );
                CHECK( !elt.isGhostCell() ) << "aiziaziaiea";
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
                    if ( M_rangeMeshElementsGhostIds.find( eltOffProc.id() ) == M_rangeMeshElementsGhostIds.end() )
                        continue;
                    myipfaces->push_back( boost::cref( face ) );
                }
            }
            M_rangeInterProcessFaces = boost::make_tuple( mpl::size_t<MESH_FACES>(), myipfaces->begin(),myipfaces->end(),myipfaces );

        }

    void updateBoundaryFacesFullSupport()
        {
            M_rangeBoundaryFaces = boundaryfaces(M_mesh);
        }
    void updateBoundaryFacesPartialSupport()
        {
            this->updateParallelData();

            std::unordered_map<size_type,uint8_type> faceInRange;
            for ( auto const& eltWrap : this->rangeElements() )
            {
                auto const& elt = unwrap_ref( eltWrap );
                for ( uint16_type i = 0; i < mesh_type::element_type::numTopologicalFaces; ++i )
                {
                    auto const& face = elt.face(i);
                    if ( faceInRange.find( face.id() ) != faceInRange.end() )
                        faceInRange[face.id()] = 2;
                    else
                        faceInRange[face.id()] = 1;
                }
            }
            for ( auto const& faceWrap : this->rangeInterProcessFaces() )
                faceInRange[unwrap_ref(faceWrap).id()] = 3;

            typename MeshTraits<mesh_type>::faces_reference_wrapper_ptrtype mybfaces( new typename MeshTraits<mesh_type>::faces_reference_wrapper_type );
            for ( auto const& faceData : faceInRange )
            {
                if ( faceData.second != 1 )
                    continue;
                mybfaces->push_back( boost::cref( M_mesh->face( faceData.first ) ) );
            }
            M_rangeBoundaryFaces = boost::make_tuple( mpl::size_t<MESH_FACES>(),mybfaces->begin(),mybfaces->end(),mybfaces );
        }
private :
    mesh_ptrtype M_mesh;
    range_elements_type M_rangeElements;
    range_faces_type M_rangeInterProcessFaces;
    range_faces_type M_rangeBoundaryFaces;
    std::unordered_set<size_type> M_rangeMeshElementsGhostIds;

    bool M_isFullSupport;
    bool M_hasUpdatedParallelData;
    bool M_hasUpdatedBoundaryFaces;

};

} // namespace Feel

#endif
