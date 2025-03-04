/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@imag.fr>
       Date: 2011-07-21

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2011-2016 Feel++ Consortium

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
   \file createsubmesh.hpp
   \author Vincent Chabannes <vincent.chabannes@imag.fr>
   \date 2011-07-21
 */
#ifndef FEELPP_DISCR_CREATESUBMESH_HPP
#define FEELPP_DISCR_CREATESUBMESH_HPP 1

#include <feel/feelmesh/submeshdata.hpp>
#include <feel/feelmesh/enums.hpp>

namespace Feel
{
template <typename MeshType,typename IteratorRange>
class CreateSubmeshTool : public CommObject
{
    using super = CommObject;
public :

    using mesh_type = MeshType;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;
    using value_type = typename mesh_type::value_type;
    using index_type = typename mesh_type::index_type;
    using size_type = typename mesh_type::size_type;

    using range_type = IteratorRange;
    using idim_type = typename range_type::idim_t;
    using iterator_type = typename range_type::iterator_t;
    using range_mesh_type = typename range_type::mesh_t;

    typedef SubMeshData<> smd_type;
    typedef std::shared_ptr<smd_type> smd_ptrtype;

    CreateSubmeshTool( CreateSubmeshTool const& t ) = default;
    CreateSubmeshTool( CreateSubmeshTool     && t ) = default;

    CreateSubmeshTool( IteratorRange const& range,
                       worldcomm_ptr_t const& wc,
                       size_type updateComponentsMesh  )
        :
        super( wc ),
        M_listRange(),
        M_updateComponentsMesh( updateComponentsMesh ),
        M_subMeshIsOnBoundaryFaces( false ),
        M_isView( false )
        {
            M_listRange.push_back( range );
        }

    CreateSubmeshTool( std::list<IteratorRange> const& range,
                       worldcomm_ptr_t const& wc,
                       size_type updateComponentsMesh  )
        :
        super( wc ),
        M_listRange( range ),
        M_updateComponentsMesh( updateComponentsMesh ),
        M_subMeshIsOnBoundaryFaces( false ),
        M_isView( false )
        {}


    CreateSubmeshTool & operator=( CreateSubmeshTool const& t ) = default;
    CreateSubmeshTool & operator=( CreateSubmeshTool     && t ) = default;

    /**
     * build mesh using Context ctx
     *
     * if ctx has the bit EXTRACTION_KEEP_MESH_RELATION set then sub mesh data
     * is added to the mesh
     */
    mesh_ptrtype
    build( size_type ctx )
        {
            DVLOG(2) << "[createSubmeshTool] extracting mesh with context "<<  ctx;
            if ( M_listRange.empty() )
                return {};

            auto meshRangePtr = M_listRange.front().mesh()->shared_from_this();
            M_smd = std::make_shared<smd_type>( meshRangePtr );

            auto newMesh = std::make_shared<mesh_type>( this->worldCommPtr() );
            // inherit the table of markersName
            for( auto itMark : meshRangePtr->markerNames( ) )
            {
                if ( itMark.second[1] > mesh_type::nDim )
                    continue;
                newMesh->addMarkerName( itMark.first,itMark.second[0],itMark.second[1] );
            }
            // subStructuring
            if ( idim_type::value == MESH_FACES )
                newMesh->setSubStructuring(meshRangePtr->subStructuring());

            // build mesh entities
            build( *newMesh, *meshRangePtr,mpl::int_<idim_type::value>() );

            // update submesh for use
            VLOG(2) << "[createSubmeshTool] mesh update for use components:"<< this->updateComponentsMesh();
            Context ctxMeshUpdate( this->updateComponentsMesh() );
            const bool renumberPoint = ctxMeshUpdate.test( MESH_RENUMBER );
            if ( !renumberPoint )
                newMesh->updateOrderedPoints();
            newMesh->components().reset();
            newMesh->components().set ( this->updateComponentsMesh()/*MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK*/ );
            newMesh->updateForUse();

            // update submesh context
            Context c( ctx );
            if ( c.test( EXTRACTION_KEEP_MESH_RELATION ) )
                newMesh->setSubMeshData( this->subMeshData() );
            if ( c.test( EXTRACTION_KEEP_MARKERNAMES_ONLY_PRESENT ) )
                newMesh->removeMarkerNameWithoutEntity();
            if ( M_isView )
                newMesh->addMeshWithNodesShared( meshRangePtr );
            return newMesh;
        }

    /**
     * @return submesh data such as the parent element to which the generated
     * mesh elements are associated.
     *
     * this allows to improve tremendously the performance of interpolation
     * operators between meshes with relation
     */
    smd_ptrtype subMeshData() { return M_smd; }

    size_type updateComponentsMesh() const { return M_updateComponentsMesh; }

    bool subMeshIsOnBoundaryFaces() const { return M_subMeshIsOnBoundaryFaces; }
    void subMeshIsOnBoundaryFaces( bool b ) { M_subMeshIsOnBoundaryFaces=b; }

    void setIsView( bool b ) { M_isView = b; }

private:

    void build( mesh_type & newMesh, range_mesh_type const& meshRange, mpl::int_<MESH_ELEMENTS> /**/ );
    void build( mesh_type & newMesh, range_mesh_type const& meshRange, mpl::int_<MESH_FACES> /**/ );
    void build( mesh_type & newMesh, range_mesh_type const& meshRange, mpl::int_<MESH_EDGES> /**/ );
    using mapping_requireghostcells_type = std::map<rank_type,std::vector<std::tuple<size_type,size_type,std::vector<index_type>>>>;
    template <int RangeType>
    void updateParallelSubMeshGhost( mesh_type & newMesh, range_mesh_type const& meshRange,
                                     std::map<size_type,size_type> & new_node_numbers,
                                     std::map<size_type,size_type> & new_element_id,
                                     mapping_requireghostcells_type const& requireGhostCells,
                                     bool renumberPoint );

    template <int EntityType>
    auto const&
    entityExtracted( range_mesh_type const& meshRange, index_type id ) const
        {
            if constexpr ( EntityType == MESH_ELEMENTS )
                return meshRange.element( id );
            else if constexpr ( EntityType == MESH_FACES )
                return meshRange.face( id );
            else //MESH_EDGES
                return meshRange.edge( id );
        }

    void addMarkedEdgesInSubMesh( typename range_mesh_type::element_type const& oldElt,
                                  std::map<index_type,index_type> const& new_node_numbers,
                                  mesh_type & newMesh, std::set<index_type> & oldEdgeIdsDone );
    void addMarkedEdgesInSubMesh( typename range_mesh_type::face_type const& oldFace,
                                  std::map<index_type, index_type> const& new_node_numbers,
                                  mesh_type & newMesh, std::set<index_type>& oldEdgeIdsDone );

    std::list<range_type> M_listRange;
    smd_ptrtype M_smd;
    size_type M_updateComponentsMesh;
    bool M_subMeshIsOnBoundaryFaces;
    bool M_isView;
};

template <typename MeshType,typename IteratorRange>
void
CreateSubmeshTool<MeshType,IteratorRange>::addMarkedEdgesInSubMesh( typename range_mesh_type::element_type const& oldElt,
                                                                           std::map<index_type,index_type> const& new_node_numbers,
                                                                           mesh_type & newMesh, std::set<index_type> & oldEdgeIdsDone )
{
    if constexpr ( range_mesh_type::nDim == 3 )
    {
        const rank_type proc_id = newMesh.worldComm().localRank();
        for ( uint16_type s = 0; s < range_mesh_type::element_type::numLocalEdges; s++ )
        {
            if ( !oldElt.edgePtr( s ) ) continue;
            // get the corresponding edge
            auto const& oldEdge = oldElt.edge( s );
            // ignore edge if no marker assigned
            if ( !oldEdge.hasMarker() ) continue;
            size_type oldEdgeId = oldEdge.id();
            // ignore edge if already done
            if ( oldEdgeIdsDone.find( oldEdgeId ) != oldEdgeIdsDone.end() )
                continue;

            typename mesh_type::edge_type newEdge;
            newEdge.setId( oldEdgeIdsDone.size() );
            newEdge.setMarkers( oldEdge.markers() );
            newEdge.setProcessIdInPartition( proc_id );
            newEdge.setProcessId( proc_id );
            for ( uint16_type p = 0; p < newEdge.nPoints(); ++p )
                newEdge.setPoint( p, newMesh.point( new_node_numbers.find( oldEdge.point( p ).id() )->second ) );
            // add it to the list of edges
            newMesh.addEdge( newEdge );
            oldEdgeIdsDone.insert( oldEdgeId );
        }
    }
}

template <typename MeshType,typename IteratorRange>
void
CreateSubmeshTool<MeshType,IteratorRange>::addMarkedEdgesInSubMesh( typename range_mesh_type::face_type const& oldFace,
                                                                           std::map<index_type, index_type> const& new_node_numbers,
                                                                           mesh_type & newMesh, std::set<index_type>& oldEdgeIdsDone )
{
    if constexpr ( range_mesh_type::nDim == 3 )
    {
        const rank_type proc_id = newMesh.worldComm().localRank();
        for ( uint16_type s = 0; s < range_mesh_type::face_type::numLocalEdges; s++ )
        {
            //if ( !oldFace.edgePtr( s ) ) continue;
            if ( !oldFace.facePtr( s ) ) continue;
            // get the corresponding edge
            auto const& oldEdge = oldFace.edge( s );
            // ignore edge if no marker assigned
            if ( !oldEdge.hasMarker() ) continue;
            size_type oldEdgeId = oldEdge.id();
            // ignore edge if already done
            if ( oldEdgeIdsDone.find( oldEdgeId ) != oldEdgeIdsDone.end() )
                continue;

            typename mesh_type::face_type newFace;
            newFace.setId( oldEdgeIdsDone.size() );
            newFace.setMarkers( oldEdge.markers() );
            newFace.setProcessIdInPartition( proc_id );
            newFace.setProcessId( proc_id );
            // very important! updateForUse put false for internalfaces after
            newFace.setOnBoundary( true );
            for ( uint16_type p = 0; p < newFace.nPoints(); ++p )
                newFace.setPoint( p, newMesh.point( new_node_numbers.find( oldEdge.point( p ).id() )->second ) );
            // add it to the list of edges
            newMesh.addFace( newFace );
            oldEdgeIdsDone.insert( oldEdgeId );
        }
    }
}

template <typename MeshType,typename IteratorRange>
void
CreateSubmeshTool<MeshType,IteratorRange>::build( mesh_type & newMesh, range_mesh_type const& meshRange, mpl::int_<MESH_ELEMENTS> /**/ )
{
    if constexpr ( IteratorRange::entities() == MESH_ELEMENTS )
    {
        // How the nodes on this mesh will be renumbered to nodes
        // on the new_mesh.
        std::map<size_type,size_type> new_node_numbers;
        std::map<size_type,size_type> new_element_id;

        Context c( this->updateComponentsMesh() );
        const bool renumberPoint = c.test( MESH_RENUMBER );

        // the number of nodes on the new mesh, will be incremented
        size_type n_new_nodes = 0;
        size_type n_new_faces = 0;
        size_type n_new_edges = 0;
        std::set<size_type> oldFaceIdsDone, oldEdgeIdsDone;

        const rank_type proc_id = this->worldComm().localRank();
        const rank_type nProc = this->worldComm().localSize();
        mapping_requireghostcells_type requireGhostCells;

        for (auto& itList : M_listRange)
        {
            auto it = itList.template get<1>();
            auto const en = itList.template get<2>();
            for ( ; it != en; ++ it )
            {
                auto const& oldElem = boost::unwrap_ref( *it );
    #if !defined(NDEBUG)
                VLOG(2) << "create sub mesh element from "  << oldElem.id() << "\n";google::FlushLogFiles(google::GLOG_INFO);
    #endif

                // check elt to extract
                if ( nProc > 1 && oldElem.isGhostCell() )
                    continue;

                // create new active element with a copy of marker
                typename mesh_type::element_type newElem;
                newElem.setMarkers( oldElem.markers() );
                newElem.setProcessIdInPartition( proc_id );
                newElem.setProcessId( proc_id );

                // Loop over the nodes on this element.
                // We guess newElem.nPoints <= oldElem.nPoint (createP1mesh for example)
                for ( uint16_type n=0; n < newElem.nPoints(); n++ )
                {
                    auto const& oldPoint = oldElem.point( n );
                    size_type oldPointId = oldPoint.id();
                    size_type newPtId = invalid_v<size_type>;
                    auto itFindPoint = new_node_numbers.find( oldPointId );
                    if ( itFindPoint != new_node_numbers.end() )
                    {
                        newPtId = itFindPoint->second;
                    }
                    else
                    {
                        DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] insert point " << oldPoint << "\n";
                        newPtId = (renumberPoint)? n_new_nodes++ : oldPointId;
                        new_node_numbers[oldPointId] = newPtId;
                        typename mesh_type::point_type pt( newPtId, oldPoint, false, M_isView );
                        pt.setProcessIdInPartition( proc_id );
                        pt.setProcessId( proc_id );
                        pt.setMarkers( oldPoint.markers() );
                        // Add this node to the new mesh
                        newMesh.addPoint ( pt );
                        DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] number of  points " << newMesh.numPoints() << "\n";
                    }

                    // Define this element's connectivity on the new mesh
                    if ( renumberPoint )
                        CHECK ( newPtId < newMesh.numPoints() ) <<  "invalid connectivity";

                    DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] adding point old(" << oldPointId
                            << ") as point new(" << newPtId
                            << ") in element " << newElem.id() << "\n";

                    newElem.setPoint( n, newMesh.point( newPtId ) );

                } // for (unsigned int n=0 ... )

                // update neighbor partitions (TODO : not necessary internally, we can just use idInOtherpartitions map)
                newElem.setNeighborPartitionIds( oldElem.neighborPartitionIds() );
                // init process connection, indices will be set in updateForUse
                for ( auto const&[pid,eltIdInPartition] : oldElem.idInOthersPartitions() )
                    newElem.setIdInOtherPartitions(pid, invalid_v<index_type> );

                // Add an equivalent element type to the new_mesh
                auto [eit,inserted] = newMesh.addElement( newElem,true );
                auto const& [eid,e] = *eit;
                new_element_id[oldElem.id()] = eid;
                M_smd->bm.insert( typename smd_type::bm_type::value_type( eid, oldElem.id() ) );

                // add marked faces for this element
                for ( uint16_type s=0; s<oldElem.numTopologicalFaces; s++ )
                {
                    if ( !oldElem.facePtr( s ) ) continue;
                    // get the corresponding face
                    auto const& oldFace = oldElem.face( s );
                    // ignore face if no marker assigned
                    if ( !oldFace.hasMarker() ) continue;
                    size_type oldFaceId = oldFace.id();
                    // ignore face if already done
                    if( oldFaceIdsDone.find( oldFaceId ) != oldFaceIdsDone.end() )
                        continue;

                    typename mesh_type::face_type newFace;
                    newFace.setId( n_new_faces++ );
                    newFace.setMarkers( oldFace.markers() );
                    newFace.setProcessIdInPartition( proc_id );
                    newFace.setProcessId( proc_id );
                    // very important! updateForUse put false for internalfaces after
                    newFace.setOnBoundary( true );
                    for ( uint16_type p = 0; p < newFace.nPoints(); ++p )
                        newFace.setPoint( p, newMesh.point( new_node_numbers[ oldFace.point(p).id()] ) );
                    // add it to the list of faces
                    auto addFaceRes = newMesh.addFace( newFace );
                    oldFaceIdsDone.insert( oldFaceId );
                } // for (unsigned int s=0 ... )

                // add marked edges in 3d for this element
                if constexpr ( range_mesh_type::nDim == 3 )
                    this->addMarkedEdgesInSubMesh( oldElem, new_node_numbers, newMesh, oldEdgeIdsDone );

                // update ghost requirements
                if ( !oldElem.idInOthersPartitions().empty() )
                {
                    // get elt ordering for ghost cells (identity for ELEMENTS case because ghosts are supposed to be identical)
                    std::vector<index_type> oldElementOrdering;
                    for ( auto const&[neighborPid,neighborEltId] : oldElem.idInOthersPartitions() )
                        requireGhostCells[neighborPid].push_back( std::make_tuple( eid, neighborEltId, oldElementOrdering ) );
                }
            } //  for( ; it != en; ++ it )
        } // for (auto& itList : M_listRange)

        if ( nProc > 1 )
        {
            this->updateParallelSubMeshGhost<MESH_ELEMENTS>( newMesh, meshRange, new_node_numbers, new_element_id, requireGhostCells, renumberPoint );
        }

        VLOG(2) << "CreateSubmesh(MESH_ELEMENTS) done";
    }
}


/**
 * create subMesh from a range<MESH_FACES>
 */
template <typename MeshType,typename IteratorRange>
void
CreateSubmeshTool<MeshType,IteratorRange>::build( mesh_type & newMesh, range_mesh_type const& meshRange, mpl::int_<MESH_FACES> )
{
    if constexpr ( IteratorRange::entities() == MESH_FACES )
    {
        DVLOG(2) << "CreateSubmesh(MESH_FACES) start";
        std::map<size_type,size_type> new_node_numbers;
        std::map<size_type,size_type> new_element_id;

        Context c( this->updateComponentsMesh() );
        const bool renumberPoint = c.test( MESH_RENUMBER );

        // the number of nodes on the new mesh, will be incremented
        size_type n_new_nodes = 0;
        size_type n_new_faces = 0;
        std::set<size_type> oldEdgeIdsDone;

        const rank_type proc_id = newMesh.worldComm().localRank();
        const rank_type nProc = newMesh.worldComm().localSize();
        mapping_requireghostcells_type requireGhostCells;

        //-----------------------------------------------------------//

        for (auto& itList : M_listRange)
        {
            auto it = itList.template get<1>();
            auto const en = itList.template get<2>();

            DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] extracting " << std::distance(it,en)  << " faces " << "\n";
            for ( ; it != en; ++ it )
            {
                // create a new element
                auto const& oldElem = boost::unwrap_ref( *it );
                DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh]   + face : " << oldElem.id() << "\n";

                // check face to extract
                if ( nProc > 1 && oldElem.isGhostCell() )
                    continue;

                if ( this->subMeshIsOnBoundaryFaces() )
                    CHECK( oldElem.isOnBoundary() ) << "error : use mpi optimzation subMeshIsOnBoundaryFaces but an internal face is added";

                // create new active element with a copy of marker
                typename mesh_type::element_type newElem;
                newElem.setMarkers( oldElem.markers() );
                newElem.setProcessIdInPartition( proc_id );
                newElem.setProcessId( proc_id );
                // loop over the nodes on this element.
                for ( unsigned int n=0; n < newElem.nPoints(); n++ )
                {
                    auto const& oldPoint = oldElem.point( n );
                    size_type oldPointId = oldPoint.id();
                    size_type newPtId = invalid_v<size_type>;
                    auto itFindPoint = new_node_numbers.find( oldPointId );
                    if ( itFindPoint != new_node_numbers.end() )
                    {
                        newPtId = itFindPoint->second;
                    }
                    else
                    {
                        DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] insert point " << oldPoint << "\n";
                        newPtId = (renumberPoint)? n_new_nodes++ : oldPointId;
                        new_node_numbers[oldPointId] = newPtId;
                        typename mesh_type::point_type pt( newPtId, oldPoint, false, M_isView );
                        pt.setProcessIdInPartition( proc_id );
                        pt.setProcessId( proc_id );
                        pt.setMarkers( oldPoint.markers() );
                        // Add this node to the new mesh
                        newMesh.addPoint( pt );
                        DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] number of  points " << newMesh.numPoints() << "\n";
                    }

                    newElem.setPoint( n, newMesh.point( newPtId ) );

                } // end for n

                // update neighbor partitions (TODO : not necessary internally, we can just use idInOtherpartitions map)
                newElem.setNeighborPartitionIds( oldElem.neighborPartitionIds() );
                // init process connection, indices will be set in updateForUse
                for ( auto const&[pid,eltIdInPartition] : oldElem.idInOthersPartitions() )
                    newElem.setIdInOtherPartitions(pid, invalid_v<index_type> );

                // Add an equivalent element type to the new_mesh
                auto [eit,inserted] = newMesh.addElement( newElem, true );
                auto const& [eid,e] = *eit;
                // update mesh relation
                new_element_id[oldElem.id()]= eid;
                M_smd->bm.insert( typename smd_type::bm_type::value_type( eid, oldElem.id() ) );
                DVLOG(2) << "connecting new face to " << e.id() << " face " << oldElem.id();
                // add marked edges in 3d as marked faces for this element
                if constexpr ( range_mesh_type::nDim == 3 )
                    this->addMarkedEdgesInSubMesh( oldElem, new_node_numbers, newMesh, oldEdgeIdsDone );
                // update ghost requirements
                if ( !oldElem.idInOthersPartitions().empty() )
                {
                    // get elt ordering for ghost cells
                    std::vector<index_type> oldElementOrdering( oldElem.nPoints(), invalid_v<index_type> );
                     for ( uint16_type n = 0; n < oldElem.nPoints(); n++ )
                    {
                        auto const& oldPoint = oldElem.point( n );
                        oldElementOrdering[n] = oldPoint.id();
                    }
                    for ( auto const&[neighborPid,neighborEltId] : oldElem.idInOthersPartitions() )
                        requireGhostCells[neighborPid].push_back( std::make_tuple( eid, neighborEltId, oldElementOrdering ) );
                }
            } // end for it
        } // for (auto& itList : M_listRange)

        if ( nProc > 1 )
        {
            this->updateParallelSubMeshGhost<MESH_FACES>( newMesh, meshRange, new_node_numbers, new_element_id, requireGhostCells, renumberPoint );
        }

        DVLOG(2) << "CreateSubmesh(MESH_FACES) done";
    }
}

template <typename MeshType,typename IteratorRange>
void
CreateSubmeshTool<MeshType,IteratorRange>::build( mesh_type & newMesh, range_mesh_type const& meshRange, mpl::int_<MESH_EDGES> /**/ )
{
    if constexpr ( IteratorRange::entities() == MESH_EDGES )
    {
        DVLOG(2) << "CreateSubmesh(MESH_EDGES) start";
        Context c( this->updateComponentsMesh() );
        const bool renumberPoint = c.test( MESH_RENUMBER );

        // the number of nodes on the new mesh, will be incremented
        size_type n_new_nodes = 0;
        size_type n_new_edges = 0;

        const int proc_id = this->worldComm().localRank();
        const int nProc = this->worldComm().localSize();
        std::map<size_type,size_type> new_node_numbers;
        std::map<size_type,size_type> new_element_id;
        mapping_requireghostcells_type requireGhostCells;
        //-----------------------------------------------------------//

        auto itListRange = M_listRange.begin();
        auto const enListRange = M_listRange.end();
        for ( ; itListRange!=enListRange ; ++itListRange)
        {
            auto it = itListRange->template get<1>();
            auto const en = itListRange->template get<2>();

            DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] extracting " << std::distance(it,en)  << " edges " << "\n";
            for ( ; it != en; ++ it )
            {
                auto const& oldElem = boost::unwrap_ref( *it );
                DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh]   + face : " << oldElem.id() << "\n";

                // check elt to extract
                if ( nProc > 1 && oldElem.isGhostCell() )
                    continue;

                // create new active element with a copy of marker
                typename mesh_type::element_type newElem;
                newElem.setMarkers( oldElem.markers() );
                newElem.setProcessIdInPartition( proc_id );
                newElem.setProcessId( proc_id );

                DVLOG(2) << "\n oldElem.nPoints " << oldElem.nPoints() << "\n";
                // Loop over the nodes on this element.
                for ( unsigned int n=0; n < newElem.nPoints(); n++ )
                {
                    auto const& oldPoint = oldElem.point( n );
                    size_type oldPointId = oldPoint.id();
                    size_type newPtId = invalid_v<size_type>;
                    auto itFindPoint = new_node_numbers.find( oldPointId );
                    if ( itFindPoint != new_node_numbers.end() )
                    {
                        newPtId = itFindPoint->second;
                    }
                    else
                    {
                        DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] insert point " << oldPoint << "\n";
                        newPtId = (renumberPoint)? n_new_nodes++ : oldPointId;
                        new_node_numbers[oldPointId] = newPtId;
                        typename mesh_type::point_type pt( newPtId, oldPoint, false, M_isView );
                        pt.setProcessIdInPartition( proc_id );
                        pt.setProcessId( proc_id );
                        pt.setMarkers( oldPoint.markers() );
                        // Add this node to the new mesh
                        newMesh.addPoint( pt );
                    }

                    newElem.setPoint( n, newMesh.point( newPtId ) );
                } // end for n
                DCHECK( newElem.pointPtr(0) ) << "invalid point 0 in edge";
                DCHECK( newElem.pointPtr(1) ) << "invalid point 1 in edge";
                // update neighbor partitions (TODO : not necessary internally, we can just use idInOtherpartitions map)
                newElem.setNeighborPartitionIds( oldElem.neighborPartitionIds() );
                // init process connection, indices will be set in updateForUse
                for ( auto const&[pid,eltIdInPartition] : oldElem.idInOthersPartitions() )
                    newElem.setIdInOtherPartitions(pid, invalid_v<index_type> );

                // Add an equivalent element type to the new_mesh
                auto [eit,inserted] = newMesh.addElement( newElem, true );
                auto const& [eid,e] = *eit;
                // update mesh relation
                new_element_id[oldElem.id()]= eid;
                M_smd->bm.insert( typename smd_type::bm_type::value_type( eid, oldElem.id() ) );

                // update ghost requirements
                if ( !oldElem.idInOthersPartitions().empty() )
                {
                    // get elt ordering for ghost cells
                    std::vector<index_type> oldElementOrdering( oldElem.nPoints(), invalid_v<index_type> );
                    for ( uint16_type n = 0; n < oldElem.nPoints(); n++ )
                    {
                        auto const& oldPoint = oldElem.point( n );
                        oldElementOrdering[n] = oldPoint.id();
                    }
                    for ( auto const&[neighborPid,neighborEltId] : oldElem.idInOthersPartitions() )
                        requireGhostCells[neighborPid].push_back( std::make_tuple( eid, neighborEltId, oldElementOrdering ) );
                }
            } // end for it
        } // for ( ; itListRange!=enListRange ; ++itListRange)

        if ( nProc > 1 )
        {
            this->updateParallelSubMeshGhost<MESH_EDGES>( newMesh, meshRange, new_node_numbers, new_element_id, requireGhostCells, renumberPoint );
        }

        DVLOG(2) << "CreateSubmesh(MESH_EDGES) done";
    }
}



template <typename MeshType,typename IteratorRange>
template <int RangeType>
void
CreateSubmeshTool<MeshType,IteratorRange>::updateParallelSubMeshGhost( mesh_type & newMesh, range_mesh_type const& meshRange,
                                                                       std::map<size_type,size_type> & new_node_numbers,
                                                                       std::map<size_type,size_type> & new_element_id,
                                                                       mapping_requireghostcells_type const& requireGhostCells,
                                                                       bool renumberPoint )
{
    using element_type = typename mesh_type::element_type;
    using point_type = typename mesh_type::point_type;

    const rank_type proc_id = newMesh.worldComm().localRank();
    const rank_type nProc = newMesh.worldComm().localSize();

    auto const& neighborSubdomains = meshRange.neighborSubdomains();
    //int neighborSubdomains = neighborSubdomains.size();
    int nbRequest=2*neighborSubdomains.size();

    mapping_requireghostcells_type dataToSend, dataToRecv;

    size_type n_new_nodes = new_node_numbers.size();

    // prepare mpi comm
    mpi::request * reqs = new mpi::request[nbRequest];
    int cptRequest=0;

    for ( auto const& [pid,data] : requireGhostCells )
        CHECK( neighborSubdomains.find( pid ) != neighborSubdomains.end() ) << fmt::format("pid {} not in mesh neighborSubdomains",pid);

    // get size of data to transfer
    std::map<rank_type,std::size_t> sizeRecv;
    std::map<rank_type,std::size_t> sizeSend;
    for ( rank_type neighborRank : neighborSubdomains )
    {
        auto itFind = requireGhostCells.find( neighborRank );
        if ( itFind != requireGhostCells.end() )
            dataToSend[neighborRank] = itFind->second;
        sizeSend[neighborRank] = dataToSend[neighborRank].size();
        reqs[cptRequest++] = this->worldComm().localComm().isend( neighborRank , 0, sizeSend[neighborRank] );
        reqs[cptRequest++] = this->worldComm().localComm().irecv( neighborRank , 0, sizeRecv[neighborRank] );
    }
    // wait all requests
    mpi::wait_all(reqs, reqs + cptRequest);

    // first send
    cptRequest=0;
    for ( rank_type neighborRank : neighborSubdomains )
    {
        std::size_t nSendData = dataToSend.at(neighborRank).size();
        reqs[cptRequest++] = newMesh.worldComm().localComm().isend( neighborRank, 0, dataToSend.at(neighborRank).data(), nSendData );

        std::size_t nRecvData = sizeRecv[neighborRank];
        dataToRecv[neighborRank].resize( nRecvData );
        reqs[cptRequest++] = newMesh.worldComm().localComm().irecv( neighborRank, 0, dataToRecv[neighborRank].data(), nRecvData );
    }
    // wait all requests
    mpi::wait_all(reqs, reqs + cptRequest);

    for ( auto const& [rankRecv,dataToRecvOnProc] : dataToRecv )
    {
        CHECK( proc_id != rankRecv ) << fmt::format("should be a ghost element : process rank: {} process active: {}",proc_id,rankRecv);
        for ( auto const& [activeEltId,currentEltId,currentEltOrdering] : dataToRecvOnProc )
        {
            auto itFindSubmeshElt = new_element_id.find( currentEltId );
            if ( itFindSubmeshElt == new_element_id.end() )
            {
                auto const& oldElem = this->entityExtracted<range_type::entities()>( meshRange, currentEltId );
                // create new active element with a copy of marker
                element_type newElem;
                newElem.setMarkers( oldElem.markers() );
                newElem.setProcessIdInPartition( proc_id );
                newElem.setProcessId( rankRecv );
                newElem.addNeighborPartitionId( rankRecv );
                // update neighbor partitions (TODO : not necessary internally, we can just use idInOtherpartitions map)
                // newElem.setNeighborPartitionIds( oldElem.neighborPartitionIds() );

                // Loop over the nodes on this element.
                for ( uint16_type n=0; n < newElem.nPoints(); n++ )
                {
                    auto const& oldPoint = !currentEltOrdering.empty() ? meshRange.point( currentEltOrdering.at(n) ) : oldElem.point( n );
                    size_type oldPointId = oldPoint.id();
                    size_type newPtId = invalid_v<size_type>;
                    auto itFindPoint = new_node_numbers.find( oldPointId );
                    if ( itFindPoint != new_node_numbers.end() )
                    {
                        newPtId = itFindPoint->second;
                    }
                    else
                    {
                        DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] insert point " << oldPoint << "\n";
                        newPtId = (renumberPoint)? n_new_nodes++ : oldPointId;
                        new_node_numbers[oldPointId] = newPtId;
                        point_type pt( newPtId, oldPoint, false, M_isView );
                        pt.setProcessIdInPartition( proc_id );
                        pt.setProcessId( proc_id );
                        pt.setMarkers( oldPoint.markers() );
                        // Add this node to the new mesh
                        newMesh.addPoint ( pt );
                        DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] number of  points " << newMesh.numPoints() << "\n";
                    }

                    newElem.setPoint( n, newMesh.point( newPtId ) );
                    DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] adding point old(" << oldPointId
                             << ") as point new(" << newPtId << ") in element " << newElem.id() << "\n";
                }

                // update id in other part
                newElem.setIdInOtherPartitions( rankRecv, activeEltId );

                // Add an equivalent element type to the new_mesh
                auto [eit,inserted] = newMesh.addElement( newElem,true );
                auto const& [eid,e] = *eit;
                new_element_id[oldElem.id()] = eid;
                M_smd->bm.insert( typename smd_type::bm_type::value_type( eid, oldElem.id() ) );
            } // if ( itFindSubmeshElt == new_element_id.end() )
            else
            {
                // update for use
                auto & newElem = newMesh.elementIterator( itFindSubmeshElt->second )->second;
                CHECK( newElem.processId() == rankRecv ) << fmt::format("should be a ghost element : process id: {} process active: {}",newElem.processId(),rankRecv);
                newElem.addNeighborPartitionId( rankRecv );
                newElem.setIdInOtherPartitions( rankRecv, activeEltId );
            }

        }
    }
}



namespace detail
{
template <int Mode,typename ... Ts>
auto createSubmesh( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && range = args.get(_range);

    auto && worldcomm = args.get_else_invocable( _worldcomm, [&range](){ return range.mesh()->worldCommPtr(); } );
    size_type context = args.get_else(_context, EXTRACTION_KEEP_MESH_RELATION|EXTRACTION_KEEP_MARKERNAMES_ONLY_PRESENT );
    size_type update = args.get_else(_update, MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES);
    bool only_on_boundary_faces = args.get_else(_only_on_boundary_faces, false);
    bool view = args.get_else(_view,false);

    //auto && mesh = args0.get(_mesh);
    //using mesh_type = decltype(mesh);

    using range_type = submeshrange_t<std::decay_t<decltype(range)>>;
    using iterator_type = typename range_type::iterator_t;
    using range_mesh_type = typename range_type::mesh_t;

    using mesh_type = mp11::mp_if_c<
        range_type::entities() == MESH_ELEMENTS,
        elements_mesh_t<range_mesh_type>,
        mp11::mp_if_c<
            range_type::entities() == MESH_FACES,
            facets_mesh_t<range_mesh_type>,
            edges_mesh_t<range_mesh_type>
            >
        >;

    using submeshtool_type =  mp11::mp_if_c<
        Mode == 0,
        CreateSubmeshTool<mesh_type,range_type>,
        CreateSubmeshTool<typename mesh_type::P1_mesh_type,range_type>
        >;
    submeshtool_type t( range,worldcomm,update );
    t.subMeshIsOnBoundaryFaces( only_on_boundary_faces );
    t.setIsView( view );
    return t.build(context);
#if 0
    if constexpr ( Mode == 0 )
    {
        CreateSubmeshTool<mesh_type,range_type> t( range,worldcomm,update );
        t.subMeshIsOnBoundaryFaces( only_on_boundary_faces );
        t.setIsView( view );
        return t.build(context);
    }
    else
    {
        CreateSubmeshTool<typename mesh_type::P1_mesh_type,range_type> t( range,worldcomm,update );
        t.subMeshIsOnBoundaryFaces( only_on_boundary_faces );
        t.setIsView( view );
        return t.build(context);
    }
#endif
}
}

template <typename ... Ts>
auto createSubmesh( Ts && ... v )
{
    return Feel::detail::createSubmesh<0>( std::forward<Ts>(v)... );
}
template <typename ... Ts>
auto createSubmeshP1( Ts && ... v )
{
    return Feel::detail::createSubmesh<1>( std::forward<Ts>(v)... );
}

/**
 * @}
 */

} // namespace Feel

#endif // createsubmesh
