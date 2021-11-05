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

#include <boost/mpl/if.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION >= 105600
#include <boost/phoenix/stl/algorithm.hpp>
#else
#include <boost/spirit/home/phoenix/stl/algorithm.hpp>
#endif
#include <feel/feelmesh/submeshdata.hpp>
#include <feel/feelmesh/meshsupport.hpp>


namespace Feel
{
template <typename C, typename V, int T, typename IndexT> class Mesh;

template <typename MeshType,typename IteratorRange, int TheTag=MeshType::tag>
class CreateSubmeshTool : public CommObject
{
public :
    using super = CommObject;
    
    typedef IteratorRange range_type;
    typedef typename boost::tuples::template element<0, range_type>::type idim_type;
    typedef typename boost::tuples::template element<1, range_type>::type iterator_type;

    static const uint16_type tag = TheTag;
    typedef MeshType mesh_type;
    typedef typename mesh_type::value_type value_type;
    using index_type = typename mesh_type::index_type;
    using size_type = typename mesh_type::size_type;
    
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    typedef typename mpl::if_<mpl::bool_<mesh_type::shape_type::is_simplex>,
                              mpl::identity< Mesh< Simplex< mesh_type::nDim,mesh_type::nOrder,mesh_type::nRealDim>, value_type, tag, index_type > >,
                              mpl::identity< Mesh< Hypercube<mesh_type::nDim,mesh_type::nOrder,mesh_type::nRealDim>, value_type, tag, index_type > > >::type::type mesh_elements_type;
    typedef std::shared_ptr<mesh_elements_type> mesh_elements_ptrtype;

    typedef typename mpl::if_<mpl::bool_<mesh_type::shape_type::is_simplex>,
                              mpl::identity< Mesh< Simplex< mesh_type::nDim-1,mesh_type::nOrder,mesh_type::nRealDim>, value_type, tag, index_type > >,
                              mpl::identity< Mesh< Hypercube<mesh_type::nDim-1,mesh_type::nOrder,mesh_type::nRealDim>, value_type, tag, index_type > > >::type::type mesh_faces_type;
    typedef std::shared_ptr<mesh_faces_type> mesh_faces_ptrtype;

    typedef typename mpl::if_<mpl::bool_<mesh_type::shape_type::is_simplex>,
                              mpl::identity< Mesh< Simplex< (mesh_type::nDim==3)?mesh_type::nDim-2:mesh_type::nDim-1,mesh_type::nOrder,mesh_type::nRealDim>, value_type, tag, index_type > >,
                              mpl::identity< Mesh< Hypercube<(mesh_type::nDim==3)?mesh_type::nDim-2:mesh_type::nDim-1,mesh_type::nOrder,mesh_type::nRealDim>, value_type, tag, index_type > > >::type::type mesh_edges_type;
    typedef std::shared_ptr<mesh_edges_type> mesh_edges_ptrtype;

    typedef typename mpl::if_< mpl::equal_to< idim_type ,mpl::size_t<MESH_ELEMENTS> >,
                               mpl::identity<mesh_elements_type>,
                               typename mpl::if_< mpl::equal_to< idim_type ,mpl::size_t<MESH_FACES> >,
                                                  mpl::identity<mesh_faces_type>,
                                                  mpl::identity<mesh_edges_type> >::type>::type::type mesh_build_type;

    typedef std::shared_ptr<mesh_build_type> mesh_build_ptrtype;
    typedef SubMeshData<> smd_type;
    typedef std::shared_ptr<smd_type> smd_ptrtype;

    CreateSubmeshTool( CreateSubmeshTool const& t ) = default;
    CreateSubmeshTool( CreateSubmeshTool     && t ) = default;

    CreateSubmeshTool( std::shared_ptr<MeshType> inputMesh,
                       IteratorRange const& range,
                       worldcomm_ptr_t const& wc,
                       size_type updateComponentsMesh  )
        :
        CreateSubmeshTool( std::make_shared<MeshSupport<MeshType>>( inputMesh ), range, wc, updateComponentsMesh )
        {}

    CreateSubmeshTool( std::shared_ptr<MeshType> inputMesh,
                       std::list<IteratorRange> const& range,
                       worldcomm_ptr_t const& wc,
                       size_type updateComponentsMesh  )
        :
        CreateSubmeshTool( std::make_shared<MeshSupport<MeshType>>( inputMesh ), range, wc, updateComponentsMesh )
        {}

    CreateSubmeshTool( std::shared_ptr<MeshSupport<MeshType>> inputMesh,
                       IteratorRange const& range,
                       worldcomm_ptr_t const& wc,
                       size_type updateComponentsMesh  )
        :
        super( wc ),
        M_meshSupport( inputMesh ),
        M_mesh( inputMesh->mesh() ),
        M_listRange(),
        M_smd( new smd_type( M_mesh ) ),
        M_updateComponentsMesh( updateComponentsMesh ),
        M_subMeshIsOnBoundaryFaces( false ),
        M_isView( false )
        {
            M_listRange.push_back( range );
        }

    CreateSubmeshTool( std::shared_ptr<MeshSupport<MeshType>> inputMesh,
                       std::list<IteratorRange> const& range,
                       worldcomm_ptr_t const& wc,
                       size_type updateComponentsMesh  )
        :
        super( wc ),
        M_meshSupport( inputMesh ),
        M_mesh( inputMesh->mesh() ),
        M_listRange( range ),
        M_smd( new smd_type( M_mesh ) ),
        M_updateComponentsMesh( updateComponentsMesh ),
        M_subMeshIsOnBoundaryFaces( false ),
        M_isView( false )
        {}


    CreateSubmeshTool & operator=( CreateSubmeshTool const& t ) = default;
    CreateSubmeshTool & operator=( CreateSubmeshTool     && t ) = default;

    /**
     * build mesh using Context ctx
     */
    mesh_build_ptrtype
    build()
        {
            DVLOG(2) << "[createSubmeshTool] extracting mesh\n";
            return build( mpl::int_<idim_type::value>() );
        }
    /**
     * build mesh using Context ctx
     *
     * if ctx has the bit EXTRACTION_KEEP_MESH_RELATION set then sub mesh data
     * is added to the mesh
     */
    mesh_build_ptrtype
    build( size_type ctx )
        {
            DVLOG(2) << "[createSubmeshTool] extracting mesh with context "<<  ctx;
            auto m = build( mpl::int_<idim_type::value>() );
            Context c( ctx );
            if ( c.test( EXTRACTION_KEEP_MESH_RELATION ) )
                m->setSubMeshData( this->subMeshData() );
            if ( c.test( EXTRACTION_KEEP_MARKERNAMES_ONLY_PRESENT ) )
                 m->removeMarkerNameWithoutEntity();
            if ( M_isView )
                m->addMeshWithNodesShared( M_mesh );
            return m;
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

    mesh_elements_ptrtype build( mpl::int_<MESH_ELEMENTS> /**/ );
    mesh_faces_ptrtype build( mpl::int_<MESH_FACES> /**/ );
    mesh_edges_ptrtype build( mpl::int_<MESH_EDGES> /**/ );

    void updateParallelInputRange( mesh_build_ptrtype const& newMesh, std::set<size_type> const& faceIdsInRange,
                                   std::map<rank_type,std::set<std::pair<size_type,size_type> > > const& faceIdsToCheck,
                                   std::map<size_type, rank_type> & faceIdsInRangeToBeGhost );

    template <int RangeType,typename SubMeshType>
    void updateParallelSubMesh( std::shared_ptr<SubMeshType> & newMesh,
                                std::map<size_type,size_type> & new_node_numbers,
                                std::map<size_type,size_type> const& new_element_id,
                                std::map<rank_type,std::set<boost::tuple<size_type,size_type> > > const& ghostCellsFind,
                                bool renumberPoint, size_type n_new_nodes );
    typename MeshType::element_type const&
    entityExtracted( size_type id, rank_type pid, mpl::int_<MESH_ELEMENTS> /**/ ) const;
    typename MeshType::face_type const&
    entityExtracted( size_type id, rank_type pid, mpl::int_<MESH_FACES> /**/ ) const;
    typename MeshType::edge_type const&
    entityExtracted( size_type id, rank_type pid, mpl::int_<MESH_EDGES> /**/ ) const;

    std::shared_ptr<MeshSupport<MeshType>> M_meshSupport;
    mesh_ptrtype M_mesh;
    std::list<range_type> M_listRange;
    smd_ptrtype M_smd;
    size_type M_updateComponentsMesh;
    bool M_subMeshIsOnBoundaryFaces;
    bool M_isView;
};

namespace detail
{

template<typename MeshType,typename SubMeshType, typename IndexT=typename MeshType::index_type>
void
addMarkedEdgesInSubMesh( std::shared_ptr<MeshType> const& mesh, typename MeshType::element_type const& oldElt,
                         std::map<IndexT,IndexT> const& new_node_numbers, IndexT & n_new_edges,
                         std::shared_ptr<SubMeshType> & newMesh, std::set<IndexT> & oldEdgeIdsDone )
{
    if constexpr ( MeshType::nDim == 3 )
    {
        typedef typename MeshType::edge_type edge_type;
        const int proc_id = newMesh->worldComm().localRank();

        for ( uint16_type s = 0; s < MeshType::element_type::numLocalEdges; s++ )
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

            if ( mesh->hasEdge( oldEdgeId ) )
            {
                edge_type newEdge;
                newEdge.setId( n_new_edges++ );
                newEdge.setMarkers( oldEdge.markers() );
                newEdge.setProcessIdInPartition( proc_id );
                newEdge.setProcessId( proc_id );
                for ( uint16_type p = 0; p < newEdge.nPoints(); ++p )
                    newEdge.setPoint( p, newMesh->point( new_node_numbers.find( oldEdge.point( p ).id() )->second ) );
                // add it to the list of edges
                newMesh->addEdge( newEdge );
                oldEdgeIdsDone.insert( oldEdgeId );
            }
        }
    }
}

template <typename MeshType, typename SubMeshType, typename IndexT = typename MeshType::index_type>
void addMarkedEdgesInSubMesh( std::shared_ptr<MeshType> const& mesh, typename MeshType::face_type const& oldFace,
                              std::map<IndexT, IndexT> const& new_node_numbers, IndexT& n_new_faces,
                              std::shared_ptr<SubMeshType>& newMesh, std::set<IndexT>& oldEdgeIdsDone )
{
    if constexpr ( MeshType::nDim == 3 )
    {
        typedef typename SubMeshType::face_type new_face_type;
        const int proc_id = newMesh->worldComm().localRank();

        for ( uint16_type s = 0; s < MeshType::face_type::numLocalEdges; s++ )
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

            if ( mesh->hasEdge( oldEdgeId ) )
            {
                new_face_type newFace;
                newFace.setId( n_new_faces++ );
                newFace.setMarkers( oldEdge.markers() );
                newFace.setProcessIdInPartition( proc_id );
                newFace.setProcessId( proc_id );
                // very important! updateForUse put false for internalfaces after
                newFace.setOnBoundary( true );
                for ( uint16_type p = 0; p < newFace.nPoints(); ++p )
                    newFace.setPoint( p, newMesh->point( new_node_numbers.find( oldEdge.point( p ).id() )->second ) );
                // add it to the list of edges
                newMesh->addFace( newFace );
                oldEdgeIdsDone.insert( oldEdgeId );
            }
        }
    }
}

} // namespace detail

template <typename MeshType,typename IteratorRange,int TheTag>
typename CreateSubmeshTool<MeshType,IteratorRange,TheTag>::mesh_elements_ptrtype
CreateSubmeshTool<MeshType,IteratorRange,TheTag>::build( mpl::int_<MESH_ELEMENTS> /**/ )
{
    typedef typename mesh_type::element_type element_type;
    typedef typename mesh_type::point_type point_type;
    typedef typename mesh_type::face_type face_type;

    mesh_elements_ptrtype newMesh( new mesh_elements_type(this->worldCommPtr()));

    //-----------------------------------------------------------//

    // inherit the table of markersName
    BOOST_FOREACH( auto itMark, M_mesh->markerNames() )
    {
        newMesh->addMarkerName( itMark.first,itMark.second[0],itMark.second[1] );
    }

    //-----------------------------------------------------------//

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

    const int proc_id = this->worldComm().localRank();
    const int nProc = this->worldComm().localSize();
    std::map<rank_type,std::set<boost::tuple<size_type,size_type> > > ghostCellsFind;

    for (auto& itList : M_listRange)
    {
        auto it = itList.template get<1>();
        auto const en = itList.template get<2>();
        for ( ; it != en; ++ it )
        {
            element_type const& oldElem = boost::unwrap_ref( *it );
#if !defined(NDEBUG)
            VLOG(2) << "create sub mesh element from "  << oldElem.id() << "\n";google::FlushLogFiles(google::GLOG_INFO);
#endif

            // check elt to extract
            if ( nProc > 1 && oldElem.isGhostCell() )
                continue;

            // create new active element with a copy of marker
            element_type newElem;
            newElem.setMarkers( oldElem.markers() );
            newElem.setProcessIdInPartition( proc_id );
            newElem.setProcessId( proc_id );

            // Loop over the nodes on this element.
            for ( uint16_type n=0; n < oldElem.nPoints(); n++ )
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
                    point_type pt( newPtId, oldPoint, false, M_isView );
                    pt.setProcessIdInPartition( proc_id );
                    pt.setProcessId( proc_id );
                    pt.setMarkers( oldPoint.markers() );
                    // Add this node to the new mesh
                    newMesh->addPoint ( pt );
                    DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] number of  points " << newMesh->numPoints() << "\n";

                    // save info necessary for the build of parallel mesh
                    if ( nProc > 1 && oldPoint.numberOfProcGhost() > 0  )
                    {
                        for (auto& itProcGhost : oldPoint.elementsGhost())
                        {
                            const rank_type procIdGhost = itProcGhost.first;
                            for ( size_type eltIdGhost : itProcGhost.second)
                            {
                                auto const& ghostElt = M_mesh->element( eltIdGhost );
                                ghostCellsFind[procIdGhost].insert( boost::make_tuple( ghostElt.id(),
                                                                                       ghostElt.idInOthersPartitions( ghostElt.processId() ) ) );
                            }
                        }
                    }
                }

                // Define this element's connectivity on the new mesh
                if ( renumberPoint )
                    CHECK ( newPtId < newMesh->numPoints() ) <<  "invalid connectivity";

                DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] adding point old(" << oldPointId
                         << ") as point new(" << newPtId
                         << ") in element " << newElem.id() << "\n";

                newElem.setPoint( n, newMesh->point( newPtId ) );

            } // for (unsigned int n=0 ... )

            // Add an equivalent element type to the new_mesh
            auto [eit,inserted] = newMesh->addElement( newElem,true );
            auto const& [eid,e] = *eit;
            new_element_id[oldElem.id()] = eid;
            M_smd->bm.insert( typename smd_type::bm_type::value_type( eid, oldElem.id() ) );

            // add marked faces for this element
            for ( uint16_type s=0; s<oldElem.numTopologicalFaces; s++ )
            {
                if ( !oldElem.facePtr( s ) ) continue;
                // get the corresponding face
                face_type const& oldFace = oldElem.face( s );
                // ignore face if no marker assigned
                if ( !oldFace.hasMarker() ) continue;
                size_type oldFaceId = oldFace.id();
                // ignore face if already done
                if( oldFaceIdsDone.find( oldFaceId ) != oldFaceIdsDone.end() )
                    continue;
                if ( M_mesh->hasFace( oldFaceId ) )
                {
                    face_type newFace;
                    newFace.setId( n_new_faces++ );
                    newFace.setMarkers( oldFace.markers() );
                    newFace.setProcessIdInPartition( proc_id );
                    newFace.setProcessId( proc_id );
                    // very important! updateForUse put false for internalfaces after
                    newFace.setOnBoundary( true );
                    for ( uint16_type p = 0; p < newFace.nPoints(); ++p )
                        newFace.setPoint( p, newMesh->point( new_node_numbers[ oldFace.point(p).id()] ) );
                    // add it to the list of faces
                    auto addFaceRes = newMesh->addFace( newFace );
                    oldFaceIdsDone.insert( oldFaceId );
                }
            } // for (unsigned int s=0 ... )

            // add marked edges in 3d for this element
            Feel::detail::addMarkedEdgesInSubMesh( M_mesh, oldElem, new_node_numbers, n_new_edges,
                                                   newMesh, oldEdgeIdsDone );
        } //  for( ; it != en; ++ it )
    } // for (auto& itList : M_listRange)

    if ( nProc > 1 )
    {
        this->updateParallelSubMesh<MESH_ELEMENTS>( newMesh, new_node_numbers, new_element_id, ghostCellsFind, renumberPoint, n_new_nodes );
    }

    VLOG(2) << "submesh created\n";
    // newMesh->setNumVertices( newMesh->numPoints() );
    if (!renumberPoint )
        newMesh->updateOrderedPoints();
    VLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] update face/edge info if necessary\n";
    // Prepare the new_mesh for use
    newMesh->components().reset();
    newMesh->components().set ( this->updateComponentsMesh()/*MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK*/ );
    newMesh->updateForUse();
    VLOG(2) << "CreateSubmesh(MESH_ELEMENTS) done\n";

    return newMesh;
}


/**
 * create subMesh from a range<MESH_FACES>
 */
template <typename MeshType,typename IteratorRange,int TheTag>
typename CreateSubmeshTool<MeshType,IteratorRange,TheTag>::mesh_faces_ptrtype
CreateSubmeshTool<MeshType,IteratorRange,TheTag>::build( mpl::int_<MESH_FACES> /**/ )
{
    DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] creating new mesh" << "\n";
    mesh_faces_ptrtype newMesh( new mesh_faces_type( this->worldCommPtr()) );

    newMesh->setSubStructuring(M_mesh->subStructuring());
    //-----------------------------------------------------------//
    DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] extraction mesh faces" << "\n";
    // inherit the table of markersName
    BOOST_FOREACH( auto itMark, M_mesh->markerNames() )
    {
        if ( itMark.second[1] > mesh_faces_type::nDim ) continue;
        DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] adding marker " << itMark.first <<"\n";
        newMesh->addMarkerName( itMark.first,itMark.second[0],itMark.second[1] );
    }
    //-----------------------------------------------------------//

    typedef typename mesh_type::face_type face_type;
    typedef typename mesh_faces_type::element_type new_element_type;
    typedef typename mesh_faces_type::face_type new_face_type;
    typedef typename mesh_faces_type::point_type new_point_type;

    std::map<size_type,size_type> new_node_numbers;
    std::map<size_type,size_type> new_element_id;

    Context c( this->updateComponentsMesh() );
    const bool renumberPoint = c.test( MESH_RENUMBER );

    // the number of nodes on the new mesh, will be incremented
    size_type n_new_nodes = 0;
    size_type n_new_faces = 0;
    std::set<size_type> oldEdgeIdsDone;

    const rank_type proc_id = newMesh->worldComm().localRank();
    const rank_type nProc = newMesh->worldComm().localSize();
    std::map<rank_type,std::set<boost::tuple<size_type,size_type> > > ghostCellsFind;

    //-----------------------------------------------------------//

    bool cleanInputRange = true;
    std::map<size_type, rank_type> faceIdsInRangeMustBeGhost;
    if ( cleanInputRange && nProc > 1 )
    {
        std::set<size_type> faceIdsInRange;
        std::map<rank_type,std::set<std::pair<size_type,size_type> > > faceIdsToCheck;
        for (auto& itList : M_listRange)
        {
            auto it = itList.template get<1>();
            auto const en = itList.template get<2>();
            for ( ; it != en; ++ it )
            {
                face_type const& theface = boost::unwrap_ref( *it );

                if ( theface.isConnectedTo1() )
                {
                    if ( theface.element0().isGhostCell() && theface.element1().isGhostCell() )
                        continue;
                }
                else
                {
                    if ( theface.element0().isGhostCell() )
                        continue;
                }
                faceIdsInRange.insert( theface.id() );

                if ( theface.isConnectedTo1() && (theface.element0().isGhostCell() || theface.element1().isGhostCell()) )
                {
                    rank_type otherPid = theface.partition2();
                    DCHECK( theface.idInOthersPartitions().find( otherPid ) != theface.idInOthersPartitions().end() ) << "no id stored for other partition " << otherPid;
                    faceIdsToCheck[otherPid].insert( std::make_pair(theface.id(), theface.idInOthersPartitions( otherPid )) );

                    // special cases with partial mesh support (a face is not connected to an active element in the mesh support context)
                    if ( M_meshSupport->isPartialSupport() )
                    {
                        if ( ( theface.element0().isGhostCell() && !M_meshSupport->hasElement(theface.element1().id() ) ) ||
                             ( theface.element1().isGhostCell() && !M_meshSupport->hasElement(theface.element0().id() ) ) )
                            faceIdsInRangeMustBeGhost[ theface.id() ] = invalid_v<rank_type>;
                    }
                }
            }
        }

        this->updateParallelInputRange( newMesh, faceIdsInRange, faceIdsToCheck, faceIdsInRangeMustBeGhost );
    } // nProc > 1

    //-----------------------------------------------------------//

    for (auto& itList : M_listRange)
    {
        auto it = itList.template get<1>();
        auto const en = itList.template get<2>();

        DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] extracting " << std::distance(it,en)  << " faces " << "\n";
        for ( ; it != en; ++ it )
        {
            // create a new element
            face_type const& oldElem = boost::unwrap_ref( *it );
            DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh]   + face : " << oldElem.id() << "\n";

            // check face to extract
            if ( nProc > 1 )
            {
                if ( oldElem.isGhostCell() )
                    continue;

                if ( this->subMeshIsOnBoundaryFaces() )
                    CHECK( oldElem.isOnBoundary() ) << "error : use mpi optimzation subMeshIsOnBoundaryFaces but an internal face is added";

                auto findFaceIdGhost = faceIdsInRangeMustBeGhost.find( oldElem.id() );
                if ( findFaceIdGhost != faceIdsInRangeMustBeGhost.end() )
                {
                    ghostCellsFind[findFaceIdGhost->second].insert(boost::make_tuple( findFaceIdGhost->first,
                                                                                      oldElem.idInOthersPartitions(findFaceIdGhost->second)) );
                    continue;
                }
            }

            // create new active element with a copy of marker
            new_element_type newElem;
            newElem.setMarkers( oldElem.markers() );
            newElem.setProcessIdInPartition( proc_id );
            newElem.setProcessId( proc_id );
            // loop over the nodes on this element.
            for ( unsigned int n=0; n < oldElem.nPoints(); n++ )
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
                    typename mesh_faces_type::point_type pt( newPtId, oldPoint, false, M_isView );
                    pt.setProcessIdInPartition( proc_id );
                    pt.setProcessId( proc_id );
                    pt.setMarkers( oldPoint.markers() );
                    // Add this node to the new mesh
                    newMesh->addPoint( pt );
                    DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] number of  points " << newMesh->numPoints() << "\n";

                    // save info necessary for the build of parallel mesh
                    if ( nProc > 1 && oldPoint.numberOfProcGhost() > 0 )
                    {
                        for (auto& itProcGhost : oldPoint.elementsGhost())
                        {
                            const rank_type procIdGhost = itProcGhost.first;
                            for (size_type eltIdGhost : itProcGhost.second)
                            {
                                if ( M_meshSupport->isPartialSupport() && !M_meshSupport->hasElement( eltIdGhost ) )
                                    continue;

                                auto const& ghostElt = M_mesh->element(eltIdGhost);
                                for ( uint16_type s=0; s<ghostElt.numTopologicalFaces; s++ )
                                {
                                    if ( !ghostElt.facePtr( s ) )
                                        continue;
                                    auto const& ghostFace = ghostElt.face( s );
                                     // no interprocess faces
                                    if ( ghostFace.processId() == proc_id )
                                        continue;
                                    // reduce mpi comm if possible
                                    if ( this->subMeshIsOnBoundaryFaces() && !ghostFace.isOnBoundary() )
                                        continue;

                                    bool isConnectedToActivePartition = false;
                                    for ( int pf = 0 ; pf < ghostFace.nVertices() ; ++pf )
                                        if ( ghostFace.point(pf).processId() == proc_id )
                                        {
                                            isConnectedToActivePartition = true;
                                            break;
                                        }
                                    if ( !isConnectedToActivePartition )
                                        continue;

                                    // store ghost faces to find in other process
                                    CHECK( procIdGhost == ghostElt.processId() ) << "not allow";
                                    //continue;
                                    ghostCellsFind[procIdGhost].insert(boost::make_tuple( ghostFace.id(),
                                                                                          ghostFace.idInOthersPartitions(ghostElt.processId())) );
                                }
                            }
                        }
                    } // if (nProc > 1  && oldPoint.numberOfProcGhost()>0)
                }

                newElem.setPoint( n, newMesh->point( newPtId ) );

            } // end for n

            // Add an equivalent element type to the new_mesh
            auto [eit,inserted] = newMesh->addElement( newElem, true );
            auto const& [eid,e] = *eit;
            // update mesh relation
            new_element_id[oldElem.id()]= eid;
            M_smd->bm.insert( typename smd_type::bm_type::value_type( eid, oldElem.id() ) );
            DVLOG(2) << "connecting new face to " << e.id() << " face " << oldElem.id();
            // add marked edges in 3d as marked faces for this element
            Feel::detail::addMarkedEdgesInSubMesh( M_mesh, oldElem, new_node_numbers, n_new_faces,
                                                   newMesh, oldEdgeIdsDone );
        } // end for it
    } // for (auto& itList : M_listRange)

    if ( nProc > 1 )
    {
        this->updateParallelSubMesh<MESH_FACES>( newMesh, new_node_numbers, new_element_id, ghostCellsFind, renumberPoint, n_new_nodes );
    }

    // newMesh->setNumVertices( newMesh->numPoints() );
    DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] update face/edge info if necessary\n";
    if (!renumberPoint )
        newMesh->updateOrderedPoints();
    // Prepare the new_mesh for use
    newMesh->components().reset();
    newMesh->components().set ( this->updateComponentsMesh()/*MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK*/ );
    newMesh->updateForUse();

    return newMesh;
}

template <typename MeshType,typename IteratorRange,int TheTag>
typename CreateSubmeshTool<MeshType,IteratorRange,TheTag>::mesh_edges_ptrtype
CreateSubmeshTool<MeshType,IteratorRange,TheTag>::build( mpl::int_<MESH_EDGES> /**/ )
{
    DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] creating new mesh" << "\n";
    mesh_edges_ptrtype newMesh( new mesh_edges_type( this->worldCommPtr()) );

    //-----------------------------------------------------------//
    DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] extraction mesh edges" << "\n";
    // inherit the table of markersName
    BOOST_FOREACH( auto itMark, M_mesh->markerNames() )
    {
        if ( itMark.second[1] > mesh_edges_type::nDim ) continue;
        DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] adding marker " << itMark.first <<"\n";
        newMesh->addMarkerName( itMark.first,itMark.second[0],itMark.second[1] );
    }

    //-----------------------------------------------------------//

    typedef typename mesh_edges_type::element_type new_element_type;

    Context c( this->updateComponentsMesh() );
    const bool renumberPoint = c.test( MESH_RENUMBER );

    // the number of nodes on the new mesh, will be incremented
    size_type n_new_nodes = 0;
    size_type n_new_edges = 0;

    const int proc_id = this->worldComm().localRank();
    const int nProc = this->worldComm().localSize();
    std::map<size_type,size_type> new_node_numbers;
    std::map<size_type,size_type> new_element_id;
    std::map<rank_type,std::set<boost::tuple<size_type,size_type> > > ghostCellsFind;

    bool cleanInputRange = true;
    std::map<size_type, rank_type> edgeIdsInRangeMustBeGhost;
    if ( cleanInputRange && nProc > 1 )
    {
        std::set<size_type> edgeIdsInRange;
        std::map<rank_type,std::set<std::pair<size_type,size_type> > > edgeIdsToCheck;
        for (auto& itList : M_listRange)
        {
            auto it = itList.template get<1>();
            auto const en = itList.template get<2>();
            for ( ; it != en; ++ it )
            {
                auto const& theedge = boost::unwrap_ref( *it );
                if ( theedge.isGhostCell() )
                    continue;
                edgeIdsInRange.insert( theedge.id() );

                for ( auto const& eltIdGhostBase : theedge.elementsGhost() )
                {
                    rank_type otherPid = eltIdGhostBase.first;
                    DCHECK( theedge.idInOthersPartitions().find( otherPid ) != theedge.idInOthersPartitions().end() ) << "no id stored for other partition " << otherPid;
                    edgeIdsToCheck[otherPid].insert( std::make_pair(theedge.id(), theedge.idInOthersPartitions( otherPid )) );
                }
            }
        }

        this->updateParallelInputRange( newMesh, edgeIdsInRange, edgeIdsToCheck, edgeIdsInRangeMustBeGhost );
    }

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
            if ( nProc > 1 )
            {
                if ( oldElem.isGhostCell() )
                    continue;

                auto findEdgeIdGhost = edgeIdsInRangeMustBeGhost.find( oldElem.id() );
                if ( findEdgeIdGhost != edgeIdsInRangeMustBeGhost.end() )
                {
                    ghostCellsFind[findEdgeIdGhost->second].insert(boost::make_tuple( findEdgeIdGhost->first,
                                                                                      oldElem.idInOthersPartitions(findEdgeIdGhost->second)) );
                    continue;
                }
            }

            // create new active element with a copy of marker
            typename mesh_edges_type::element_type newElem;
            newElem.setMarkers( oldElem.markers() );
            newElem.setProcessIdInPartition( proc_id );
            newElem.setProcessId( proc_id );

            DVLOG(2) << "\n oldElem.nPoints " << oldElem.nPoints() << "\n";
            // Loop over the nodes on this element.
            for ( unsigned int n=0; n < oldElem.nPoints(); n++ )
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
                    typename mesh_edges_type::point_type pt( newPtId, oldPoint, false, M_isView );
                    pt.setProcessIdInPartition( proc_id );
                    pt.setProcessId( proc_id );
                    pt.setMarkers( oldPoint.markers() );
                    // Add this node to the new mesh
                    newMesh->addPoint( pt );

                    // save info necessary for the build of parallel mesh
                    if ( nProc > 1 && oldPoint.numberOfProcGhost() > 0 )
                    {
                        for (auto& itProcGhost : oldPoint.elementsGhost())
                        {
                            const rank_type procIdGhost = itProcGhost.first;
                            for (size_type eltIdGhost : itProcGhost.second)
                            {
                                auto const& ghostElt = M_mesh->element(eltIdGhost);
                                for ( uint16_type s=0; s<ghostElt.numEdges/*numTopologicalFaces*/; s++ )
                                {
                                    if ( !ghostElt.edgePtr( s ) )
                                        continue;
                                    auto const& ghostEdge = ghostElt.edge( s );
                                     // no interprocess faces
                                    if ( ghostEdge.processId() == proc_id )
                                        continue;
                                    // reduce mpi comm if possible
                                    if ( this->subMeshIsOnBoundaryFaces() && !ghostEdge.isOnBoundary() )
                                        continue;

                                    bool isConnectedToActivePartition = false;
                                    for ( int pf = 0 ; pf < ghostEdge.nVertices() ; ++pf )
                                        if ( ghostEdge.point(pf).processId() == proc_id )
                                        {
                                            isConnectedToActivePartition = true;
                                            break;
                                        }
                                    if ( !isConnectedToActivePartition )
                                        continue;

                                    // store ghost faces to find in other process
                                    CHECK( procIdGhost == ghostElt.processId() ) << "not allow";
                                    ghostCellsFind[procIdGhost].insert(boost::make_tuple( ghostEdge.id(),
                                                                                          ghostEdge.idInOthersPartitions(ghostElt.processId())) );
                                }
                            }
                        }
                    } // if (nProc > 1  && oldPoint.numberOfProcGhost()>0)

                }

                newElem.setPoint( n, newMesh->point( newPtId ) );
#if 0
                newElem.setFace( n, newMesh->point( newPtId ) );
#endif
            } // end for n
            CHECK( newElem.pointPtr(0) ) << "invalid point 0 in edge";
            CHECK( newElem.pointPtr(1) ) << "invalid point 1 in edge";
#if 0
            CHECK( newElem.facePtr(0) ) << "invalid face 0 in edge";
            CHECK( newElem.facePtr(1) ) << "invalid face 1 in edge";
#endif
            // Add an equivalent element type to the new_mesh
            auto [eit,inserted] = newMesh->addElement( newElem, true );
            auto const& [eid,e] = *eit;
            // update mesh relation
            new_element_id[oldElem.id()]= eid;
            M_smd->bm.insert( typename smd_type::bm_type::value_type( eid, oldElem.id() ) );
        } // end for it
    } // for ( ; itListRange!=enListRange ; ++itListRange)

    if ( nProc > 1 )
    {
        this->updateParallelSubMesh<MESH_EDGES>( newMesh, new_node_numbers, new_element_id, ghostCellsFind, renumberPoint, n_new_nodes );
    }

    // newMesh->setNumVertices( newMesh->numPoints() );
    DVLOG(2) << "[CreateSubmesh] update face/edge info if necessary\n";
    if ( !renumberPoint )
        newMesh->updateOrderedPoints();
    // Prepare the new_mesh for use
    newMesh->components().reset();
    newMesh->components().set ( this->updateComponentsMesh()/*MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK*/ );
    newMesh->updateForUse();

    return newMesh;
}


template <typename MeshType,typename IteratorRange,int TheTag>
typename MeshType::element_type const&
CreateSubmeshTool<MeshType,IteratorRange,TheTag>::entityExtracted( size_type id, rank_type pid, mpl::int_<MESH_ELEMENTS> /**/ ) const
{
    CHECK( M_mesh->hasElement( id ) ) << "no element with id " << id;
    return M_mesh->element( id );
}
template <typename MeshType,typename IteratorRange,int TheTag>
typename MeshType::face_type const&
CreateSubmeshTool<MeshType,IteratorRange,TheTag>::entityExtracted( size_type id, rank_type pid, mpl::int_<MESH_FACES> /**/ ) const
{
    CHECK( M_mesh->hasFace( id ) ) << "no face with id " << id;
    return M_mesh->face( id );
}
template <typename MeshType,typename IteratorRange,int TheTag>
typename MeshType::edge_type const&
CreateSubmeshTool<MeshType,IteratorRange,TheTag>::entityExtracted( size_type id, rank_type pid, mpl::int_<MESH_EDGES> /**/ ) const
{
    CHECK( M_mesh->hasEdge( id ) ) << "no face with id " << id;
    return M_mesh->edge( id );
}

template <typename MeshType,typename IteratorRange,int TheTag>
void
CreateSubmeshTool<MeshType,IteratorRange,TheTag>::updateParallelInputRange( mesh_build_ptrtype const& newMesh, std::set<size_type> const& entityIdsInRange,
                                                                            std::map<rank_type,std::set<std::pair<size_type,size_type> > > const& entityIdsToCheck,
                                                                            std::map<size_type, rank_type> & entityIdsInRangeMustBeGhost )
{
    const rank_type proc_id = newMesh->worldComm().localRank();
    // init container
    std::map<rank_type,std::vector<size_type> > dataToSend,dataToRecv;
    std::map<rank_type,std::vector<std::tuple<size_type,bool>> > memoryDataMpi;
    for ( auto const& faceIdToCheck : entityIdsToCheck )
    {
        rank_type otherPid = faceIdToCheck.first;
        for ( auto const& faceIdPair : faceIdToCheck.second)
        {
            dataToSend[otherPid].push_back( faceIdPair.second );

            if ( entityIdsInRangeMustBeGhost.find( faceIdPair.first ) != entityIdsInRangeMustBeGhost.end() )
                memoryDataMpi[otherPid].push_back( std::make_tuple(faceIdPair.first,true ) );
            else
                memoryDataMpi[otherPid].push_back( std::make_tuple(faceIdPair.first,false ) );
        }
    }
    // prepare mpi comm
    int neighborSubdomains = M_mesh->neighborSubdomains().size();
    int nbRequest=2*neighborSubdomains;
    mpi::request * reqs = new mpi::request[nbRequest];
    int cptRequest=0;

    // get size of data to transfer
    std::map<rank_type,size_type> sizeRecv;
    for ( rank_type neighborRank : M_mesh->neighborSubdomains() )
    {
        reqs[cptRequest++] = newMesh->worldComm().localComm().isend( neighborRank , 0, (size_type)dataToSend[neighborRank].size() );
        reqs[cptRequest++] = newMesh->worldComm().localComm().irecv( neighborRank , 0, sizeRecv[neighborRank] );
    }
    // wait all requests
    mpi::wait_all(reqs, reqs + cptRequest);

    // first send/recv
    cptRequest=0;
    for ( rank_type neighborRank : M_mesh->neighborSubdomains() )
    {
        int nSendData = dataToSend[neighborRank].size();
        if ( nSendData > 0 )
            reqs[cptRequest++] = newMesh->worldComm().localComm().isend( neighborRank , 0, &(dataToSend[neighborRank][0]), nSendData );

        int nRecvData = sizeRecv[neighborRank];
        dataToRecv[neighborRank].resize( nRecvData );
        if ( nRecvData > 0 )
            reqs[cptRequest++] = newMesh->worldComm().localComm().irecv( neighborRank , 0, &(dataToRecv[neighborRank][0]), nRecvData );
    }
    // wait all requests
    mpi::wait_all(reqs, reqs + cptRequest);
    // treat recv and prepare resend data
    // reponse=0 -> entity is not in the range
    // reponse=1 -> entity is in the range and no constraint
    // reponse=2 -> entity is in the range but considered ghost
    std::map<rank_type,std::vector<int> > dataToReSend,dataToReRecv;
    for ( auto const& dataToRecvBase : dataToRecv )
    {
        rank_type procRecv = dataToRecvBase.first;
        auto const& idsRecv = dataToRecvBase.second;
        dataToReSend[procRecv].resize( idsRecv.size(), 0 );
        for (int k=0;k<idsRecv.size();++k)
        {
            size_type idRecv = idsRecv[k];
            if ( entityIdsInRange.find( idRecv ) == entityIdsInRange.end() )
                continue;

            auto itFindEntityeMustBeGhost = entityIdsInRangeMustBeGhost.find( idRecv );
            if ( itFindEntityeMustBeGhost != entityIdsInRangeMustBeGhost.end() )
            {
                CHECK( itFindEntityeMustBeGhost->second == invalid_v<rank_type> ) << "something wrong";
                dataToReSend[procRecv][k] = 2;
            }
            else
                dataToReSend[procRecv][k] = 1;
        }
    }
    // second send/recv
    cptRequest=0;
    for ( rank_type neighborRank : M_mesh->neighborSubdomains() )
    {
        int nSendData = dataToReSend[neighborRank].size();
        if ( nSendData > 0 )
            reqs[cptRequest++] = newMesh->worldComm().localComm().isend( neighborRank, 1, &(dataToReSend[neighborRank][0]), nSendData );

        int nRecvData = dataToSend[neighborRank].size();
        dataToReRecv[neighborRank].resize( nRecvData );
        if ( nRecvData > 0 )
            reqs[cptRequest++] = newMesh->worldComm().localComm().irecv( neighborRank, 1, &(dataToReRecv[neighborRank][0]), nRecvData );
    }
    // wait all requests
    mpi::wait_all(reqs, reqs + cptRequest);
    // delete reqs because finish comm
    delete [] reqs;

    for ( auto const& dataToReRecvBase : dataToReRecv )
    {
        rank_type procRecv = dataToReRecvBase.first;
        auto const& idsPresentRecv = dataToReRecvBase.second;
        for (int k=0;k<idsPresentRecv.size();++k)
        {
            if ( idsPresentRecv[k] == 0 || idsPresentRecv[k] == 2 )
                continue;

            size_type idRemove = std::get<0>( memoryDataMpi[procRecv][k] );
            bool isForcedToBeGhost =  std::get<1>( memoryDataMpi[procRecv][k] );

            if ( !isForcedToBeGhost && procRecv >= proc_id )
                continue;

            if ( idsPresentRecv[k] == 1 )
            {
                if ( entityIdsInRangeMustBeGhost.find( idRemove ) == entityIdsInRangeMustBeGhost.end() )
                    entityIdsInRangeMustBeGhost[idRemove] = procRecv;
                else
                    entityIdsInRangeMustBeGhost[idRemove] = std::min(procRecv, entityIdsInRangeMustBeGhost[idRemove] );
            }
        }
    }

}

template <typename MeshType,typename IteratorRange,int TheTag>
template <int RangeType,typename SubMeshType>
void
CreateSubmeshTool<MeshType,IteratorRange,TheTag>::updateParallelSubMesh( std::shared_ptr<SubMeshType> & newMesh,
                                                                         std::map<size_type,size_type> & new_node_numbers,
                                                                         std::map<size_type,size_type> const& new_element_id,
                                                                         std::map<rank_type,std::set<boost::tuple<size_type,size_type> > > const& ghostCellsFind,
                                                                         bool renumberPoint,size_type n_new_nodes )
{
    typedef typename SubMeshType::element_type element_type;
    typedef typename SubMeshType::point_type point_type;

    const int proc_id = newMesh->worldComm().localRank();
    const int nProc = newMesh->worldComm().localSize();

    int neighborSubdomains = M_mesh->neighborSubdomains().size();
    int nbRequest=2*neighborSubdomains;

    std::map<rank_type,std::vector<std::pair<size_type,size_type> > > memoryDataToSend;
    std::map<rank_type,std::vector<size_type> > dataToSend,dataToRecv,dataToReSend,dataToRecv2;
    // init container
    for ( rank_type neighborRank : M_mesh->neighborSubdomains() )
    {
        memoryDataToSend[neighborRank].clear();
        dataToSend[neighborRank].clear();
        dataToRecv[neighborRank].clear();
        dataToReSend[neighborRank].clear();
        dataToRecv2[neighborRank].clear();
    }
    // update container for first send
    for ( rank_type neighborRank : M_mesh->neighborSubdomains() )
    {
        auto itGhostCellsFind = ghostCellsFind.find( neighborRank );
        if ( itGhostCellsFind != ghostCellsFind.end() )
        {
            for ( auto const& pairId : itGhostCellsFind->second )
            {
                const size_type idEltInMyProc = pairId.template get<0>();
                const size_type idEltInOtherProc = pairId.template get<1>();
                dataToSend[neighborRank].push_back( idEltInOtherProc );
                memoryDataToSend[neighborRank].push_back( std::make_pair(idEltInMyProc,idEltInOtherProc) );
            }
        }
    }
    // prepare mpi comm
    mpi::request * reqs = new mpi::request[nbRequest];
    int cptRequest=0;

    // get size of data to transfer
    std::map<rank_type,size_type> sizeRecv;
    for ( rank_type neighborRank : M_mesh->neighborSubdomains() )
    {
        reqs[cptRequest++] = this->worldComm().localComm().isend( neighborRank , 0, (size_type)dataToSend[neighborRank].size() );
        reqs[cptRequest++] = this->worldComm().localComm().irecv( neighborRank , 0, sizeRecv[neighborRank] );
    }
    // wait all requests
    mpi::wait_all(reqs, reqs + cptRequest);

    // first send
    cptRequest=0;
    for ( rank_type neighborRank : M_mesh->neighborSubdomains() )
    {
        int nSendData = dataToSend[neighborRank].size();
        if ( nSendData > 0 )
            reqs[cptRequest++] = newMesh->worldComm().localComm().isend( neighborRank, 0, &(dataToSend[neighborRank][0]), nSendData );

        int nRecvData = sizeRecv[neighborRank];
        dataToRecv[neighborRank].resize( nRecvData );
        if ( nRecvData > 0 )
            reqs[cptRequest++] = newMesh->worldComm().localComm().irecv( neighborRank, 0, &(dataToRecv[neighborRank][0]), nRecvData );
    }
    // wait all requests
    mpi::wait_all(reqs, reqs + cptRequest);

    // treat first recv and resend answer
    cptRequest=0;
    for ( auto const& [rankRecv,dataToRecvOnProc] : dataToRecv )
    {
        int nData = dataToRecvOnProc.size();
        dataToReSend[rankRecv].resize( nData );
        if ( nData == 0 )
            continue;
        std::fill( dataToReSend[rankRecv].begin(), dataToReSend[rankRecv].end(), invalid_v<size_type> );
        for ( int k=0;k<nData;++k )
        {
            size_type idEltRecv = dataToRecvOnProc[k];
            // search id
            auto const itFindId = new_element_id.find(idEltRecv);
            if ( itFindId != new_element_id.end() )
            {
                size_type idEltInNewMesh = itFindId->second;
                // update NeighborPartition for this active elt
                CHECK( newMesh->hasElement( idEltInNewMesh) ) << "mesh has not elt whit id " << idEltInNewMesh << "\n";
                auto & eltToUpdate = newMesh->elementIterator( idEltInNewMesh )->second;
                eltToUpdate.addNeighborPartitionId( rankRecv );
                dataToReSend[rankRecv][k] = idEltInNewMesh;
            }
        }
        // second send
        reqs[cptRequest++] = newMesh->worldComm().localComm().isend( rankRecv, 1, &(dataToReSend[rankRecv][0]), nData );
    }

    // second recv
    for ( rank_type neighborRank : M_mesh->neighborSubdomains() )
    {
        int nRecvData = dataToSend[neighborRank].size();
        dataToRecv2[neighborRank].resize( nRecvData );
        if ( nRecvData > 0 )
            reqs[cptRequest++] = newMesh->worldComm().localComm().irecv( neighborRank, 1, &(dataToRecv2[neighborRank][0]), nRecvData );
    }
    // wait all requests
    mpi::wait_all(reqs, reqs + cptRequest/*nbRequest*/);
    // delete reqs because finish comm
    delete [] reqs;

    std::map<size_type,std::pair<size_type,rank_type> > ghostOldEltDone;
    std::map<size_type,std::vector<std::pair<rank_type,size_type> > > mapActiveEltDuplicatedInWorld;
    // treat second recv and build ghost elements
    for ( auto const& dataToRecvPair : dataToRecv2 )
    {
        rank_type rankRecv = dataToRecvPair.first;
        auto itFindRankInMemory = memoryDataToSend.find(rankRecv);
        CHECK( itFindRankInMemory != memoryDataToSend.end() ) << "missing info in memory";
        for ( size_type k=0;k<dataToRecvPair.second.size();++k )
        {
            size_type idEltActiveInOtherProc = dataToRecvPair.second[k];
            if ( idEltActiveInOtherProc == invalid_v<size_type> )
                continue;

            auto const& memoryPair = memoryDataToSend[rankRecv][k];
            const size_type oldEltId = memoryPair.first;

            auto const& oldElem = this->entityExtracted( oldEltId, rankRecv, mpl::int_<RangeType>() );
            CHECK( oldElem.id() == oldEltId ) << "invalid id";

            auto itFindNewActiveElt = new_element_id.find(oldEltId);
            bool is_not_stored = itFindNewActiveElt == new_element_id.end();
#if 0
            LOG_IF( WARNING, is_not_stored == false )
                << "this element is already present on the new mesh\n";
            if ( is_not_stored == false ) continue;
#endif
            if ( is_not_stored )
            {
                // keep only ghost which are connected to active part by a point at least
                bool hasPointConnectionToActivePart = false;
                for ( uint16_type n=0; n < oldElem.nVertices(); n++ )
                {
                    auto const& oldPoint = oldElem.point( n );
                    size_type oldPointId = oldPoint.id();
                    auto itFindPoint = new_node_numbers.find( oldPointId );
                    if ( itFindPoint == new_node_numbers.end() )
                        continue;
                    size_type newPtId = itFindPoint->second;
                    if ( newMesh->point( newPtId ).processId() == proc_id )
                    {
                        hasPointConnectionToActivePart = true;
                        break;
                    }
                }
                if ( !hasPointConnectionToActivePart )
                {
                    //std::cout << "ignore elt " << oldElem.G() << "\n";
                    continue;
                }

                // if ghost element already build, update only new neigboring data process
                auto itFindGhostOld = ghostOldEltDone.find( oldElem.id() );
                if ( itFindGhostOld != ghostOldEltDone.end() )
                {
                    auto & eltModified = newMesh->elementIterator( itFindGhostOld->second.first )->second;
                    eltModified.setIdInOtherPartitions( rankRecv, idEltActiveInOtherProc );
                    if ( rankRecv < eltModified.processId() )
                    {
                        eltModified.setProcessId( rankRecv );
                        ghostOldEltDone[oldElem.id()]= std::make_pair( eltModified.id(),rankRecv);
                    }
                    continue;
                }

                // create a new elem with partitioning infos
                CHECK( rankRecv != oldElem.pidInPartition() && proc_id == oldElem.pidInPartition() ) << "invalid rank id";
                element_type newElem;
                newElem.setMarkers( oldElem.markers() );
                newElem.setProcessIdInPartition( proc_id );
                newElem.setProcessId( rankRecv );
                newElem.addNeighborPartitionId( rankRecv );

                // Loop over the nodes on this element.
                for ( uint16_type n=0; n < oldElem.nPoints(); n++ )
                {
                    auto const& oldPoint = oldElem.point( n );
                    size_type oldPointId = oldPoint.id();
                    size_type newPtId = invalid_v<size_type>;
                    auto itFindPoint = new_node_numbers.find( oldPointId );
                    if ( itFindPoint != new_node_numbers.end() )
                        newPtId = itFindPoint->second;
                    else//if ( itFindPoint == new_node_numbers.end() )
                    {
                        DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] insert point " << oldPoint << "\n";
                        newPtId = (renumberPoint)? n_new_nodes++ : oldPointId;
                        new_node_numbers[oldPointId] = newPtId;
                        // create point and keep default process id because point not in partition ( invalid_rank_value_type )
                        point_type pt( newPtId, oldPoint, false, M_isView );
                        pt.setProcessIdInPartition( proc_id );
                        pt.setProcessId( invalid_rank_type_value );
                        pt.setMarkers( oldPoint.markers() );
                        // Add this node to the new mesh
                        newMesh->addPoint ( pt );
                        DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] number of  points " << newMesh->numPoints() << "\n";
                        // Increment the new node counter
                        //n_new_nodes++;
                    }
#if 0
                    else if ( newMesh->point( newPtId ).processId() != invalid_rank_type_value )
                    {
                        // update NeighborPartition for this points
                        auto ptToUpdate = newMesh->pointIterator( newPtId );
                        newMesh->points().modify( ptToUpdate, Feel::detail::UpdateNeighborPartition( rankRecv ) );
                    }
#endif

                    // Define this element's connectivity on the new mesh
                    if ( renumberPoint )
                        CHECK ( newPtId < newMesh->numPoints() ) << "invalid connectivity";

                    newElem.setPoint( n, newMesh->point( newPtId ) );
                    // if ( RangeType == MESH_EDGES )
                    //    newElem.setFace( n, newMesh->point( newPtId ) );

                    DVLOG(2) << "[Mesh<Shape,T>::CreateSubmesh] adding point old(" << oldPointId
                             << ") as point new(" << newPtId << ") in element " << newElem.id() << "\n";

                } // for (uint16_type n=0 ... )

                // update id in other part
                newElem.setIdInOtherPartitions( rankRecv, idEltActiveInOtherProc );

                // Add an equivalent element type to the new_mesh
                auto [eit,inserted] = newMesh->addElement( newElem, true );
                auto const& [newEltId,e] = *eit;

                ghostOldEltDone[oldElem.id()]= std::make_pair(newEltId,e.processId());
                // update mesh relation
                M_smd->bm.insert( typename smd_type::bm_type::value_type( newEltId, oldEltId ) );
            }
            else // already stored as active element
            {
                CHECK(false ) << "a duplicated active element is not allow";
                size_type newEltId = itFindNewActiveElt->second;
                mapActiveEltDuplicatedInWorld[newEltId].push_back( std::make_pair( rankRecv, idEltActiveInOtherProc ) );
            }

        } // for ( size_type k )
    } // for ( auto dataRecv2 )

#if 0
    // maybe some active elements are duplicated in parallel mesh at interprocess zone,
    // we consider an unique active element with the minimal pid
    for ( auto const& dataEltDuplicated : mapActiveEltDuplicatedInWorld )
    {
        size_type newId = dataEltDuplicated.first;
        auto eltIt = newMesh->elementIterator( newId );

        rank_type minPid = proc_id;
        std::set<rank_type> allpid;
        allpid.insert( proc_id );
        for ( auto const& dataOtherProc : dataEltDuplicated.second )
        {
            rank_type otherPid = dataOtherProc.first;
            if ( otherPid < minPid )
                minPid = otherPid;
            allpid.insert( otherPid );
        }
        allpid.erase( minPid );

        newMesh->elements().modify( eltIt, [&minPid,&allpid,&dataEltDuplicated] (element_type & e)
                                    {
                                        for ( rank_type opid : allpid )
                                            e.addNeighborPartitionId( opid );
                                        for ( auto const& dataOtherProc : dataEltDuplicated.second )
                                            e.setIdInOtherPartitions( dataOtherProc.first,dataOtherProc.second );
                                        e.setProcessId( minPid );
                                    } );
    }
#endif
}


namespace detail
{
template<typename MeshArgType,typename RangeArgType>
struct CreateSubmeshTrait
{
    using mesh_or_support_type = typename Feel::remove_shared_ptr< typename Feel::meta::remove_all< MeshArgType >::type >::type;
    using mesh_type = typename mpl::if_c< is_mesh_v<mesh_or_support_type>, mesh_or_support_type, typename mesh_or_support_type::mesh_type>::type;
    using range_type = typename Feel::meta::remove_all< RangeArgType >::type;
    using builder_type = CreateSubmeshTool<mesh_type,typename Feel::detail::submeshrangetype<range_type>::type,mesh_type::tag>;
    using ptrtype = typename builder_type::mesh_build_ptrtype;
};
}

template <typename MeshArgType,typename RangeArgType>
typename Feel::detail::CreateSubmeshTrait<MeshArgType,RangeArgType>::ptrtype
createSubmesh( NA::arguments<
               typename na::mesh::template required_as_t<MeshArgType>,
               typename na::range::template required_as_t<RangeArgType>,
               typename na::worldcomm::template required_as_t<worldcomm_ptr_t>,
               typename na::context::template required_as_t<size_type>,
               typename na::update::template required_as_t<size_type>,
               typename na::only_on_boundary_faces::template required_as_t<bool>,
               typename na::view::template required_as_t<bool>
               > && args )
{
    auto && mesh = args.get(_mesh);
    auto && range = args.get(_range);
    auto && worldcomm = args.get(_worldcomm);
    size_type context = args.get(_context);
    size_type update = args.get(_update);
    bool only_on_boundary_faces = args.get(_only_on_boundary_faces);
    bool view =  args.get(_view);

    typename Feel::detail::CreateSubmeshTrait<MeshArgType,RangeArgType>::builder_type t( mesh,range,worldcomm,update );
    t.subMeshIsOnBoundaryFaces( only_on_boundary_faces );
    t.setIsView( view );
    return t.build(context);
}

template <typename ... Ts>
auto createSubmesh( Ts && ... v )
{
    auto args0 = NA::make_arguments( std::forward<Ts>(v)... )
        .add_default_arguments( NA::make_default_argument( _context, EXTRACTION_KEEP_MESH_RELATION|EXTRACTION_KEEP_MARKERNAMES_ONLY_PRESENT ),
                                NA::make_default_argument( _update, MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES),
                                NA::make_default_argument( _only_on_boundary_faces, false ),
                                NA::make_default_argument( _view, false )
                                );
        ;
    auto && mesh = args0.get(_mesh);
    auto && range = args0.get(_range);
    auto args = std::move( args0 ).add_default_arguments( NA::make_default_argument( _worldcomm, mesh->worldCommPtr() ) );

    using arg_mesh_type = decltype(mesh);
    using arg_range_type = decltype(range);
    return createSubmesh<arg_mesh_type,arg_range_type>( std::move( args ) );
}

/**
 * @}
 */

} // namespace Feel

#endif // createsubmesh
