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
#ifndef FEELPP_CREATESUBMESH_HPP
#define FEELPP_CREATESUBMESH_HPP 1

#include <boost/mpl/identity.hpp>
#include <boost/mpl/if.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION >= 105600
#include <boost/phoenix/stl/algorithm.hpp>
#else
#include <boost/spirit/home/phoenix/stl/algorithm.hpp>
#endif
#include <feel/feelmesh/submeshdata.hpp>

namespace Feel
{
template <typename C, typename V, int T>
class Mesh;

template <typename MeshType, typename IteratorRange, int TheTag = MeshType::tag>
class CreateSubmeshTool
{
  public:
    typedef IteratorRange range_type;
    typedef typename boost::tuples::template element<0, range_type>::type idim_type;
    typedef typename boost::tuples::template element<1, range_type>::type iterator_type;

    static const uint16_type tag = TheTag;
    typedef MeshType mesh_type;
    typedef typename mesh_type::value_type value_type;

    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef typename mpl::if_<mpl::bool_<mesh_type::shape_type::is_simplex>,
                              mpl::identity<Mesh<Simplex<mesh_type::nDim, mesh_type::nOrder, mesh_type::nRealDim>, value_type, tag>>,
                              mpl::identity<Mesh<Hypercube<mesh_type::nDim, mesh_type::nOrder, mesh_type::nRealDim>, value_type, tag>>>::type::type mesh_elements_type;
    typedef boost::shared_ptr<mesh_elements_type> mesh_elements_ptrtype;

    typedef typename mpl::if_<mpl::bool_<mesh_type::shape_type::is_simplex>,
                              mpl::identity<Mesh<Simplex<mesh_type::nDim - 1, mesh_type::nOrder, mesh_type::nRealDim>, value_type, tag>>,
                              mpl::identity<Mesh<Hypercube<mesh_type::nDim - 1, mesh_type::nOrder, mesh_type::nRealDim>, value_type, tag>>>::type::type mesh_faces_type;
    typedef boost::shared_ptr<mesh_faces_type> mesh_faces_ptrtype;

    typedef typename mpl::if_<mpl::bool_<mesh_type::shape_type::is_simplex>,
                              mpl::identity<Mesh<Simplex<( mesh_type::nDim == 3 ) ? mesh_type::nDim - 2 : mesh_type::nDim - 1, mesh_type::nOrder, mesh_type::nRealDim>, value_type, tag>>,
                              mpl::identity<Mesh<Hypercube<( mesh_type::nDim == 3 ) ? mesh_type::nDim - 2 : mesh_type::nDim - 1, mesh_type::nOrder, mesh_type::nRealDim>, value_type, tag>>>::type::type mesh_edges_type;
    typedef boost::shared_ptr<mesh_edges_type> mesh_edges_ptrtype;

    typedef typename mpl::if_<mpl::equal_to<idim_type, mpl::size_t<MESH_ELEMENTS>>,
                              mpl::identity<mesh_elements_type>,
                              typename mpl::if_<mpl::equal_to<idim_type, mpl::size_t<MESH_FACES>>,
                                                mpl::identity<mesh_faces_type>,
                                                mpl::identity<mesh_edges_type>>::type>::type::type mesh_build_type;

    typedef boost::shared_ptr<mesh_build_type> mesh_build_ptrtype;
    typedef SubMeshData smd_type;
    typedef boost::shared_ptr<smd_type> smd_ptrtype;

    CreateSubmeshTool( CreateSubmeshTool const& t ) = default;
    CreateSubmeshTool( CreateSubmeshTool&& t ) = default;

    CreateSubmeshTool( boost::shared_ptr<MeshType> inputMesh,
                       IteratorRange const& range,
                       WorldComm const& wc,
                       size_type updateComponentsMesh )
        : M_mesh( inputMesh ),
          M_listRange(),
          M_smd( new smd_type( inputMesh ) ),
          M_worldComm( wc ),
          M_updateComponentsMesh( updateComponentsMesh ),
          M_subMeshIsOnBoundaryFaces( false )
    {
        M_listRange.push_back( range );
    }

    CreateSubmeshTool( boost::shared_ptr<MeshType> inputMesh,
                       std::list<IteratorRange> const& range,
                       WorldComm const& wc,
                       size_type updateComponentsMesh )
        : M_mesh( inputMesh ),
          M_listRange( range ),
          M_smd( new smd_type( inputMesh ) ),
          M_worldComm( wc ),
          M_updateComponentsMesh( updateComponentsMesh ),
          M_subMeshIsOnBoundaryFaces( false )
    {
    }

    CreateSubmeshTool& operator=( CreateSubmeshTool const& t ) = default;
    CreateSubmeshTool& operator=( CreateSubmeshTool&& t ) = default;

    /**
     * build mesh using Context ctx 
     */
    mesh_build_ptrtype
    build()
    {
        DVLOG( 2 ) << "[createSubmeshTool] extracting mesh\n";
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
        DVLOG( 2 ) << "[createSubmeshTool] extracting mesh with context " << ctx;
        auto m = build( mpl::int_<idim_type::value>() );
        Context c( ctx );
        if ( c.test( EXTRACTION_KEEP_MESH_RELATION ) )
            m->setSubMeshData( this->subMeshData() );
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
    void subMeshIsOnBoundaryFaces( bool b ) { M_subMeshIsOnBoundaryFaces = b; }

  private:
    mesh_elements_ptrtype build( mpl::int_<MESH_ELEMENTS> /**/ );
    mesh_faces_ptrtype build( mpl::int_<MESH_FACES> /**/ );
    mesh_edges_ptrtype build( mpl::int_<MESH_EDGES> /**/ );
    template <int RangeType, typename SubMeshType>
    void updateParallelSubMesh( boost::shared_ptr<SubMeshType>& newMesh,
                                std::map<size_type, size_type>& new_node_numbers,
                                std::map<size_type, size_type> const& new_element_id,
                                std::map<rank_type, std::set<boost::tuple<size_type, size_type>>> const& ghostCellsFind,
                                bool renumberPoint, size_type n_new_nodes );
    typename MeshType::element_type const&
    entityExtracted( size_type id, rank_type pid, mpl::int_<MESH_ELEMENTS> /**/ ) const;
    typename MeshType::face_type const&
    entityExtracted( size_type id, rank_type pid, mpl::int_<MESH_FACES> /**/ ) const;
    typename MeshType::edge_type const&
    entityExtracted( size_type id, rank_type pid, mpl::int_<MESH_EDGES> /**/ ) const;

    mesh_ptrtype M_mesh;
    std::list<range_type> M_listRange;
    smd_ptrtype M_smd;
    WorldComm M_worldComm;
    size_type M_updateComponentsMesh;
    bool M_subMeshIsOnBoundaryFaces;
};

namespace detail
{

template <typename MeshType, typename SubMeshType>
void addMarkedEdgesInSubMesh( boost::shared_ptr<MeshType> const& mesh, typename MeshType::element_type const& oldElt,
                              std::map<size_type, size_type> const& new_node_numbers, size_type& n_new_edges,
                              boost::shared_ptr<SubMeshType>& submesh, std::set<size_type>& oldEdgeIdsDone, mpl::false_ /**/ )
{
}
template <typename MeshType, typename SubMeshType>
void addMarkedEdgesInSubMesh( boost::shared_ptr<MeshType> const& mesh, typename MeshType::element_type const& oldElt,
                              std::map<size_type, size_type> const& new_node_numbers, size_type& n_new_edges,
                              boost::shared_ptr<SubMeshType>& newMesh, std::set<size_type>& oldEdgeIdsDone, mpl::true_ /**/ )
{
    typedef typename MeshType::edge_type edge_type;
    const int proc_id = newMesh->worldComm().localRank();

    for ( uint16_type s = 0; s < MeshType::element_type::numLocalEdges; s++ )
    {
        if ( !oldElt.edgePtr( s ) ) continue;
        // get the corresponding edge
        auto const& oldEdge = oldElt.edge( s );
        // ignore edge if no marker assigned
        if ( oldEdge.marker().isOff() ) continue;
        size_type oldEdgeId = oldEdge.id();
        // ignore edge if already done
        if ( oldEdgeIdsDone.find( oldEdgeId ) != oldEdgeIdsDone.end() )
            continue;

        if ( mesh->hasEdge( oldEdgeId ) )
        {
            edge_type newEdge;
            newEdge.setId( n_new_edges++ );
            newEdge.setMarker( oldEdge.marker().value() );
            newEdge.setMarker2( oldEdge.marker2().value() );
            newEdge.setMarker3( oldEdge.marker3().value() );
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

template <typename MeshType, typename SubMeshType>
void addMarkedEdgesInSubMesh( boost::shared_ptr<MeshType> const& mesh, typename MeshType::face_type const& oldFace,
                              std::map<size_type, size_type> const& new_node_numbers, size_type& n_new_faces,
                              boost::shared_ptr<SubMeshType>& newMesh, std::set<size_type>& oldEdgeIdsDone, mpl::false_ /**/ )
{
}
template <typename MeshType, typename SubMeshType>
void addMarkedEdgesInSubMesh( boost::shared_ptr<MeshType> const& mesh, typename MeshType::face_type const& oldFace,
                              std::map<size_type, size_type> const& new_node_numbers, size_type& n_new_faces,
                              boost::shared_ptr<SubMeshType>& newMesh, std::set<size_type>& oldEdgeIdsDone, mpl::true_ /**/ )
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
        if ( oldEdge.marker().isOff() ) continue;
        size_type oldEdgeId = oldEdge.id();
        // ignore edge if already done
        if ( oldEdgeIdsDone.find( oldEdgeId ) != oldEdgeIdsDone.end() )
            continue;

        if ( mesh->hasEdge( oldEdgeId ) )
        {
            new_face_type newFace;
            newFace.setId( n_new_faces++ );
            newFace.setMarker( oldEdge.marker().value() );
            newFace.setMarker2( oldEdge.marker2().value() );
            newFace.setMarker3( oldEdge.marker3().value() );
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

} // namespace detail

template <typename MeshType, typename IteratorRange, int TheTag>
typename CreateSubmeshTool<MeshType, IteratorRange, TheTag>::mesh_elements_ptrtype
    CreateSubmeshTool<MeshType, IteratorRange, TheTag>::build( mpl::int_<MESH_ELEMENTS> /**/ )
{
    typedef typename mesh_type::element_type element_type;
    typedef typename mesh_type::point_type point_type;
    typedef typename mesh_type::face_type face_type;

    mesh_elements_ptrtype newMesh( new mesh_elements_type( M_worldComm ) );

    //-----------------------------------------------------------//

    // inherit the table of markersName
    BOOST_FOREACH ( auto itMark, M_mesh->markerNames() )
    {
        newMesh->addMarkerName( itMark.first, itMark.second[0], itMark.second[1] );
    }

    //-----------------------------------------------------------//

    // How the nodes on this mesh will be renumbered to nodes
    // on the new_mesh.
    std::map<size_type, size_type> new_node_numbers;
    std::map<size_type, size_type> new_element_id;

    Context c( this->updateComponentsMesh() );
    const bool renumberPoint = c.test( MESH_RENUMBER );

    // the number of nodes on the new mesh, will be incremented
    size_type n_new_nodes = 0;
    size_type n_new_faces = 0;
    size_type n_new_edges = 0;
    std::set<size_type> oldFaceIdsDone, oldEdgeIdsDone;

    const int proc_id = M_worldComm.localRank();
    const int nProc = M_worldComm.localSize();
    std::map<rank_type, std::set<boost::tuple<size_type, size_type>>> ghostCellsFind;

    for ( auto& itList : M_listRange )
    {
        auto it = itList.template get<1>();
        auto const en = itList.template get<2>();
        for ( ; it != en; ++it )
        {
            element_type const& oldElem = *it;
#if !defined( NDEBUG )
            VLOG( 2 ) << "create sub mesh element from " << oldElem.id() << "\n";
            google::FlushLogFiles( google::GLOG_INFO );
#endif

            // check elt to extract
            if ( nProc > 1 && oldElem.isGhostCell() )
                continue;

            // create new active element with a copy of marker
            element_type newElem;
            newElem.setMarker( oldElem.marker().value() );
            newElem.setMarker2( oldElem.marker2().value() );
            newElem.setMarker3( oldElem.marker3().value() );
            newElem.setProcessIdInPartition( proc_id );
            newElem.setProcessId( proc_id );

            // Loop over the nodes on this element.
            for ( uint16_type n = 0; n < oldElem.nPoints(); n++ )
            {
                auto const& oldPoint = oldElem.point( n );
                size_type oldPointId = oldPoint.id();
                size_type newPtId = invalid_size_type_value;
                auto itFindPoint = new_node_numbers.find( oldPointId );
                if ( itFindPoint != new_node_numbers.end() )
                {
                    newPtId = itFindPoint->second;
                }
                else
                {
                    DVLOG( 2 ) << "[Mesh<Shape,T>::CreateSubmesh] insert point " << oldPoint << "\n";
                    newPtId = ( renumberPoint ) ? n_new_nodes++ : oldPointId;
                    new_node_numbers[oldPointId] = newPtId;
                    point_type pt( newPtId, oldPoint );
                    pt.setProcessIdInPartition( proc_id );
                    pt.setProcessId( proc_id );
                    pt.setMarker( oldPoint.marker().value() );
                    pt.setMarker2( oldPoint.marker2().value() );
                    pt.setMarker3( oldPoint.marker3().value() );
                    // Add this node to the new mesh
                    newMesh->addPoint( pt );
                    DVLOG( 2 ) << "[Mesh<Shape,T>::CreateSubmesh] number of  points " << newMesh->numPoints() << "\n";

                    // save info necessary for the build of parallel mesh
                    if ( nProc > 1 && oldPoint.numberOfProcGhost() > 0 )
                    {
                        for ( auto& itProcGhost : oldPoint.elementsGhost() )
                        {
                            const rank_type procIdGhost = itProcGhost.first;
                            for ( size_type eltIdGhost : itProcGhost.second )
                            {
                                auto const& ghostElt = M_mesh->element( eltIdGhost, procIdGhost );
                                ghostCellsFind[procIdGhost].insert( boost::make_tuple( ghostElt.id(),
                                                                                       ghostElt.idInOthersPartitions( ghostElt.processId() ) ) );
                            }
                        }
                    }
                }

                // Define this element's connectivity on the new mesh
                if ( renumberPoint )
                    CHECK( newPtId < newMesh->numPoints() ) << "invalid connectivity";

                DVLOG( 2 ) << "[Mesh<Shape,T>::CreateSubmesh] adding point old(" << oldPointId
                           << ") as point new(" << newPtId
                           << ") in element " << newElem.id() << "\n";

                newElem.setPoint( n, newMesh->point( newPtId ) );

            } // for (unsigned int n=0 ... )

            // Add an equivalent element type to the new_mesh
            auto const& e = newMesh->addElement( newElem, true );
            new_element_id[oldElem.id()] = e.id();
            M_smd->bm.insert( typename smd_type::bm_type::value_type( e.id(), oldElem.id() ) );

            // add marked faces for this element
            for ( uint16_type s = 0; s < oldElem.numTopologicalFaces; s++ )
            {
                if ( !oldElem.facePtr( s ) ) continue;
                // get the corresponding face
                face_type const& oldFace = oldElem.face( s );
                // ignore face if no marker assigned
                if ( oldFace.marker().isOff() ) continue;
                size_type oldFaceId = oldFace.id();
                // ignore face if already done
                if ( oldFaceIdsDone.find( oldFaceId ) != oldFaceIdsDone.end() )
                    continue;
                if ( M_mesh->hasFace( oldFaceId ) )
                {
                    face_type newFace;
                    newFace.setId( n_new_faces++ );
                    newFace.setMarker( oldFace.marker().value() );
                    newFace.setMarker2( oldFace.marker2().value() );
                    newFace.setMarker3( oldFace.marker3().value() );
                    newFace.setProcessIdInPartition( proc_id );
                    newFace.setProcessId( proc_id );
                    // very important! updateForUse put false for internalfaces after
                    newFace.setOnBoundary( true );
                    for ( uint16_type p = 0; p < newFace.nPoints(); ++p )
                        newFace.setPoint( p, newMesh->point( new_node_numbers[oldFace.point( p ).id()] ) );
                    // add it to the list of faces
                    auto addFaceRes = newMesh->addFace( newFace );
                    oldFaceIdsDone.insert( oldFaceId );
                }
            } // for (unsigned int s=0 ... )

            // add marked edges in 3d for this element
            Feel::detail::addMarkedEdgesInSubMesh( M_mesh, oldElem, new_node_numbers, n_new_edges,
                                                   newMesh, oldEdgeIdsDone, mpl::bool_<mesh_type::nDim == 3>() );
        } //  for( ; it != en; ++ it )
    }     // for (auto& itList : M_listRange)

    if ( nProc > 1 )
    {
        this->updateParallelSubMesh<MESH_ELEMENTS>( newMesh, new_node_numbers, new_element_id, ghostCellsFind, renumberPoint, n_new_nodes );
    }

    VLOG( 2 ) << "submesh created\n";
    newMesh->setNumVertices( newMesh->numPoints() );

    VLOG( 2 ) << "[Mesh<Shape,T>::CreateSubmesh] update face/edge info if necessary\n";
    // Prepare the new_mesh for use
    newMesh->components().reset();
    newMesh->components().set( this->updateComponentsMesh() /*MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK*/ );
    newMesh->updateForUse();
    VLOG( 2 ) << "CreateSubmesh(MESH_ELEMENTS) done\n";

    return newMesh;
}

/**
 * create subMesh from a range<MESH_FACES>
 */
template <typename MeshType, typename IteratorRange, int TheTag>
typename CreateSubmeshTool<MeshType, IteratorRange, TheTag>::mesh_faces_ptrtype
    CreateSubmeshTool<MeshType, IteratorRange, TheTag>::build( mpl::int_<MESH_FACES> /**/ )
{
    DVLOG( 2 ) << "[Mesh<Shape,T>::CreateSubmesh] creating new mesh"
               << "\n";
    mesh_faces_ptrtype newMesh( new mesh_faces_type( M_worldComm ) );

    newMesh->setSubStructuring( M_mesh->subStructuring() );
    //-----------------------------------------------------------//
    DVLOG( 2 ) << "[Mesh<Shape,T>::CreateSubmesh] extraction mesh faces"
               << "\n";
    // inherit the table of markersName
    BOOST_FOREACH ( auto itMark, M_mesh->markerNames() )
    {
        if ( itMark.second[1] > mesh_faces_type::nDim ) continue;
        DVLOG( 2 ) << "[Mesh<Shape,T>::CreateSubmesh] adding marker " << itMark.first << "\n";
        newMesh->addMarkerName( itMark.first, itMark.second[0], itMark.second[1] );
    }
    //-----------------------------------------------------------//

    typedef typename mesh_type::face_type face_type;
    typedef typename mesh_faces_type::element_type new_element_type;
    typedef typename mesh_faces_type::face_type new_face_type;
    typedef typename mesh_faces_type::point_type new_point_type;

    std::map<size_type, size_type> new_node_numbers;
    std::map<size_type, size_type> new_element_id;

    Context c( this->updateComponentsMesh() );
    const bool renumberPoint = c.test( MESH_RENUMBER );

    // the number of nodes on the new mesh, will be incremented
    size_type n_new_nodes = 0;
    size_type n_new_faces = 0;
    std::set<size_type> oldEdgeIdsDone;

    const rank_type proc_id = newMesh->worldComm().localRank();
    const rank_type nProc = newMesh->worldComm().localSize();
    std::map<rank_type, std::set<boost::tuple<size_type, size_type>>> ghostCellsFind;

    //-----------------------------------------------------------//

    for ( auto& itList : M_listRange )
    {
        auto it = itList.template get<1>();
        auto const en = itList.template get<2>();

        DVLOG( 2 ) << "[Mesh<Shape,T>::CreateSubmesh] extracting " << std::distance( it, en ) << " faces "
                   << "\n";
        for ( ; it != en; ++it )
        {
            // create a new element
            face_type const& oldElem = *it;
            DVLOG( 2 ) << "[Mesh<Shape,T>::CreateSubmesh]   + face : " << oldElem.id() << "\n";

            // check face to extract
            if ( nProc > 1 )
            {
                if ( oldElem.isGhostCell() )
                    continue;
                if ( this->subMeshIsOnBoundaryFaces() )
                    CHECK( oldElem.isOnBoundary() ) << "error : use mpi optimzation subMeshIsOnBoundaryFaces but an internal face is added";
            }

            // create new active element with a copy of marker
            new_element_type newElem;
            newElem.setMarker( oldElem.marker().value() );
            newElem.setMarker2( oldElem.marker2().value() );
            newElem.setMarker3( oldElem.marker3().value() );
            newElem.setProcessIdInPartition( proc_id );
            newElem.setProcessId( proc_id );
            // loop over the nodes on this element.
            for ( unsigned int n = 0; n < oldElem.nPoints(); n++ )
            {
                auto const& oldPoint = oldElem.point( n );
                size_type oldPointId = oldPoint.id();
                size_type newPtId = invalid_size_type_value;
                auto itFindPoint = new_node_numbers.find( oldPointId );
                if ( itFindPoint != new_node_numbers.end() )
                {
                    newPtId = itFindPoint->second;
                }
                else
                {
                    DVLOG( 2 ) << "[Mesh<Shape,T>::CreateSubmesh] insert point " << oldPoint << "\n";
                    newPtId = ( renumberPoint ) ? n_new_nodes++ : oldPointId;
                    new_node_numbers[oldPointId] = newPtId;
                    typename mesh_faces_type::point_type pt( newPtId, oldPoint );
                    pt.setProcessIdInPartition( proc_id );
                    pt.setProcessId( proc_id );
                    pt.setMarker( oldPoint.marker().value() );
                    pt.setMarker2( oldPoint.marker2().value() );
                    pt.setMarker3( oldPoint.marker3().value() );
                    // Add this node to the new mesh
                    newMesh->addPoint( pt );
                    DVLOG( 2 ) << "[Mesh<Shape,T>::CreateSubmesh] number of  points " << newMesh->numPoints() << "\n";

                    // save info necessary for the build of parallel mesh
                    if ( nProc > 1 && oldPoint.numberOfProcGhost() > 0 )
                    {
                        for ( auto& itProcGhost : oldPoint.elementsGhost() )
                        {
                            const rank_type procIdGhost = itProcGhost.first;
                            for ( size_type eltIdGhost : itProcGhost.second )
                            {
                                auto const& ghostElt = M_mesh->element( eltIdGhost, procIdGhost );
                                for ( uint16_type s = 0; s < ghostElt.numTopologicalFaces; s++ )
                                {
                                    auto const& ghostFace = ghostElt.face( s );
                                    // reduce mpi comm if possible
                                    if ( this->subMeshIsOnBoundaryFaces() && !ghostFace.isOnBoundary() )
                                        continue;
                                    // store ghost faces to find in other process
                                    CHECK( procIdGhost == ghostElt.processId() ) << "not allow";
                                    ghostCellsFind[procIdGhost].insert( boost::make_tuple( ghostFace.id(),
                                                                                           ghostFace.idInOthersPartitions( ghostElt.processId() ) ) );
                                }
                            }
                        }
                    } // if (nProc > 1  && oldPoint.numberOfProcGhost()>0)
                }

                newElem.setPoint( n, newMesh->point( newPtId ) );

            } // end for n

            // Add an equivalent element type to the new_mesh
            auto const& e = newMesh->addElement( newElem, true );
            // update mesh relation
            new_element_id[oldElem.id()] = e.id();
            M_smd->bm.insert( typename smd_type::bm_type::value_type( e.id(), oldElem.id() ) );

            // add marked edges in 3d as marked faces for this element
            Feel::detail::addMarkedEdgesInSubMesh( M_mesh, oldElem, new_node_numbers, n_new_faces,
                                                   newMesh, oldEdgeIdsDone, mpl::bool_<mesh_type::nDim == 3>() );
        } // end for it
    }     // for (auto& itList : M_listRange)

    if ( nProc > 1 )
    {
        this->updateParallelSubMesh<MESH_FACES>( newMesh, new_node_numbers, new_element_id, ghostCellsFind, renumberPoint, n_new_nodes );
    }

    newMesh->setNumVertices( newMesh->numPoints() );
    DVLOG( 2 ) << "[Mesh<Shape,T>::CreateSubmesh] update face/edge info if necessary\n";
    // Prepare the new_mesh for use
    newMesh->components().reset();
    newMesh->components().set( this->updateComponentsMesh() /*MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK*/ );
    newMesh->updateForUse();

    return newMesh;
}

template <typename MeshType, typename IteratorRange, int TheTag>
typename CreateSubmeshTool<MeshType, IteratorRange, TheTag>::mesh_edges_ptrtype
    CreateSubmeshTool<MeshType, IteratorRange, TheTag>::build( mpl::int_<MESH_EDGES> /**/ )
{
    // we don't deal with this situation yet
    M_smd.reset();

    DVLOG( 2 ) << "[Mesh<Shape,T>::CreateSubmesh] creating new mesh"
               << "\n";
    mesh_edges_ptrtype newMesh( new mesh_edges_type( M_worldComm ) );

    //-----------------------------------------------------------//
    DVLOG( 2 ) << "[Mesh<Shape,T>::CreateSubmesh] extraction mesh edges"
               << "\n";
    // inherit the table of markersName
    BOOST_FOREACH ( auto itMark, M_mesh->markerNames() )
    {
        if ( itMark.second[1] > mesh_edges_type::nDim ) continue;
        DVLOG( 2 ) << "[Mesh<Shape,T>::CreateSubmesh] adding marker " << itMark.first << "\n";
        newMesh->addMarkerName( itMark.first, itMark.second[0], itMark.second[1] );
    }

    //-----------------------------------------------------------//

    typedef typename mesh_edges_type::element_type new_element_type;

    Context c( this->updateComponentsMesh() );
    const bool renumberPoint = c.test( MESH_RENUMBER );

    // the number of nodes on the new mesh, will be incremented
    size_type n_new_nodes = 0;
    size_type n_new_edges = 0;

    const int proc_id = M_worldComm.localRank();
    const int nProc = M_worldComm.localSize();
    std::map<size_type, size_type> new_node_numbers;
    std::map<size_type, size_type> new_element_id;
    std::map<rank_type, std::set<boost::tuple<size_type, size_type>>> ghostCellsFind;

    //-----------------------------------------------------------//

    auto itListRange = M_listRange.begin();
    auto const enListRange = M_listRange.end();
    for ( ; itListRange != enListRange; ++itListRange )
    {
        auto it = itListRange->template get<1>();
        auto const en = itListRange->template get<2>();

        DVLOG( 2 ) << "[Mesh<Shape,T>::CreateSubmesh] extracting " << std::distance( it, en ) << " edges "
                   << "\n";
        for ( ; it != en; ++it )
        {
            auto const& oldElem = boost::unwrap_ref( *it );
            DVLOG( 2 ) << "[Mesh<Shape,T>::CreateSubmesh]   + face : " << oldElem.id() << "\n";

            // check elt to extract
            if ( nProc > 1 && oldElem.isGhostCell() )
                continue;

            // create new active element with a copy of marker
            typename mesh_edges_type::element_type newElem;
            newElem.setMarker( oldElem.marker().value() );
            newElem.setMarker2( oldElem.marker2().value() );
            newElem.setMarker3( oldElem.marker3().value() );
            newElem.setProcessIdInPartition( proc_id );
            newElem.setProcessId( proc_id );

            DVLOG( 2 ) << "\n oldElem.nPoints " << oldElem.nPoints() << "\n";
            // Loop over the nodes on this element.
            for ( unsigned int n = 0; n < oldElem.nPoints(); n++ )
            {
                auto const& oldPoint = oldElem.point( n );
                size_type oldPointId = oldPoint.id();
                size_type newPtId = invalid_size_type_value;
                auto itFindPoint = new_node_numbers.find( oldPointId );
                if ( itFindPoint != new_node_numbers.end() )
                {
                    newPtId = itFindPoint->second;
                }
                else
                {
                    DVLOG( 2 ) << "[Mesh<Shape,T>::CreateSubmesh] insert point " << oldPoint << "\n";
                    newPtId = ( renumberPoint ) ? n_new_nodes++ : oldPointId;
                    new_node_numbers[oldPointId] = newPtId;
                    typename mesh_edges_type::point_type pt( newPtId, oldPoint );
                    pt.setProcessIdInPartition( proc_id );
                    pt.setProcessId( proc_id );
                    pt.setMarker( oldPoint.marker().value() );
                    pt.setMarker2( oldPoint.marker2().value() );
                    pt.setMarker3( oldPoint.marker3().value() );
                    // Add this node to the new mesh
                    newMesh->addPoint( pt );

                    // save info necessary for the build of parallel mesh
                    if ( nProc > 1 && oldPoint.numberOfProcGhost() > 0 )
                    {
                        for ( auto& itProcGhost : oldPoint.elementsGhost() )
                        {
                            const rank_type procIdGhost = itProcGhost.first;
                            for ( size_type eltIdGhost : itProcGhost.second )
                            {
                                auto const& ghostElt = M_mesh->element( eltIdGhost, procIdGhost );
                                for ( uint16_type s = 0; s < ghostElt.numEdges /*numTopologicalFaces*/; s++ )
                                {
                                    auto const& ghostEdge = ghostElt.edge( s );
                                    // reduce mpi comm if possible
                                    if ( this->subMeshIsOnBoundaryFaces() && !ghostEdge.isOnBoundary() )
                                        continue;
                                    // store ghost faces to find in other process
                                    CHECK( procIdGhost == ghostElt.processId() ) << "not allow";
                                    ghostCellsFind[procIdGhost].insert( boost::make_tuple( ghostEdge.id(),
                                                                                           ghostEdge.idInOthersPartitions( ghostElt.processId() ) ) );
                                }
                            }
                        }
                    } // if (nProc > 1  && oldPoint.numberOfProcGhost()>0)
                }

                newElem.setPoint( n, newMesh->point( newPtId ) );
                newElem.setFace( n, newMesh->point( newPtId ) );
            } // end for n
            CHECK( newElem.pointPtr( 0 ) ) << "invalid point 0 in edge";
            CHECK( newElem.pointPtr( 1 ) ) << "invalid point 1 in edge";
            CHECK( newElem.facePtr( 0 ) ) << "invalid face 0 in edge";
            CHECK( newElem.facePtr( 1 ) ) << "invalid face 1 in edge";

            // Add an equivalent element type to the new_mesh
            auto const& e = newMesh->addElement( newElem, true );
            // update mesh relation
            new_element_id[oldElem.id()] = e.id();
            //M_smd->bm.insert( typename smd_type::bm_type::value_type( e.id(), oldElem.id() ) );
        } // end for it
    }     // for ( ; itListRange!=enListRange ; ++itListRange)

    if ( nProc > 1 )
    {
        this->updateParallelSubMesh<MESH_EDGES>( newMesh, new_node_numbers, new_element_id, ghostCellsFind, renumberPoint, n_new_nodes );
    }

    newMesh->setNumVertices( newMesh->numPoints() );
    DVLOG( 2 ) << "[CreateSubmesh] update face/edge info if necessary\n";
    // Prepare the new_mesh for use
    newMesh->components().reset();
    newMesh->components().set( this->updateComponentsMesh() /*MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK*/ );
    newMesh->updateForUse();

    return newMesh;
}

template <typename MeshType, typename IteratorRange, int TheTag>
typename MeshType::element_type const&
CreateSubmeshTool<MeshType, IteratorRange, TheTag>::entityExtracted( size_type id, rank_type pid, mpl::int_<MESH_ELEMENTS> /**/ ) const
{
    CHECK( M_mesh->hasElement( id, pid ) ) << "no element with id " << id << " on proc " << pid;
    return M_mesh->element( id, pid );
}
template <typename MeshType, typename IteratorRange, int TheTag>
typename MeshType::face_type const&
CreateSubmeshTool<MeshType, IteratorRange, TheTag>::entityExtracted( size_type id, rank_type pid, mpl::int_<MESH_FACES> /**/ ) const
{
    CHECK( M_mesh->hasFace( id ) ) << "no face with id " << id;
    return M_mesh->face( id );
}
template <typename MeshType, typename IteratorRange, int TheTag>
typename MeshType::edge_type const&
CreateSubmeshTool<MeshType, IteratorRange, TheTag>::entityExtracted( size_type id, rank_type pid, mpl::int_<MESH_EDGES> /**/ ) const
{
    CHECK( M_mesh->hasEdge( id ) ) << "no face with id " << id;
    return M_mesh->edge( id );
}

template <typename MeshType, typename IteratorRange, int TheTag>
template <int RangeType, typename SubMeshType>
void CreateSubmeshTool<MeshType, IteratorRange, TheTag>::updateParallelSubMesh( boost::shared_ptr<SubMeshType>& newMesh,
                                                                                std::map<size_type, size_type>& new_node_numbers,
                                                                                std::map<size_type, size_type> const& new_element_id,
                                                                                std::map<rank_type, std::set<boost::tuple<size_type, size_type>>> const& ghostCellsFind,
                                                                                bool renumberPoint, size_type n_new_nodes )
{
    typedef typename SubMeshType::element_type element_type;
    typedef typename SubMeshType::point_type point_type;

    const int proc_id = newMesh->worldComm().localRank();
    const int nProc = newMesh->worldComm().localSize();

    int neighborSubdomains = M_mesh->neighborSubdomains().size();
    int nbRequest = 2 * neighborSubdomains;

    std::map<rank_type, std::vector<std::pair<size_type, size_type>>> memoryDataToSend;
    std::map<rank_type, std::vector<size_type>> dataToSend, dataToRecv, dataToRecv2;
    // init container
    for ( rank_type neighborRank : M_mesh->neighborSubdomains() )
    {
        memoryDataToSend[neighborRank].clear();
        dataToSend[neighborRank].clear();
        dataToRecv[neighborRank].clear();
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
                memoryDataToSend[neighborRank].push_back( std::make_pair( idEltInMyProc, idEltInOtherProc ) );
            }
        }
    }
    // prepare mpi comm
    mpi::request* reqs = new mpi::request[nbRequest];
    int cptRequest = 0;
    // first send
    auto itDataToSend = dataToSend.begin();
    auto const enDataToSend = dataToSend.end();
    for ( ; itDataToSend != enDataToSend; ++itDataToSend )
    {
        reqs[cptRequest++] = newMesh->worldComm().localComm().isend( itDataToSend->first, 0, itDataToSend->second );
    }
    // first recv
    for ( auto& dataToRecvPair : dataToRecv )
    {
        reqs[cptRequest++] = newMesh->worldComm().localComm().irecv( dataToRecvPair.first, 0, dataToRecvPair.second );
    }
    // wait all requests
    mpi::wait_all( reqs, reqs + nbRequest );

    // treat first recv and resend answer
    cptRequest = 0;
    for ( auto const& dataToRecvPair : dataToRecv )
    {
        rank_type rankRecv = dataToRecvPair.first;
        std::vector<size_type> dataToReSend;
        for ( size_type idEltRecv : dataToRecvPair.second )
        {
            // search id
            size_type idEltInNewMesh = invalid_size_type_value;
            auto const itFindId = new_element_id.find( idEltRecv );
            if ( itFindId != new_element_id.end() )
            {
                idEltInNewMesh = itFindId->second;
                // update NeighborPartition for this active elt
                CHECK( newMesh->hasElement( idEltInNewMesh ) ) << "mesh has not elt whit id " << idEltInNewMesh << "\n";
                auto eltToUpdate = newMesh->elementIterator( idEltInNewMesh );
                newMesh->elements().modify( eltToUpdate, Feel::detail::UpdateNeighborPartition( rankRecv ) );
            }
            dataToReSend.push_back( idEltInNewMesh );
        }
        // second send
        reqs[cptRequest++] = newMesh->worldComm().localComm().isend( rankRecv, 0, dataToReSend );
    }
    // second recv
    for ( auto& dataToRecvPair : dataToRecv2 )
    {
        reqs[cptRequest++] = newMesh->worldComm().localComm().irecv( dataToRecvPair.first, 0, dataToRecvPair.second );
    }
    // wait all requests
    mpi::wait_all( reqs, reqs + nbRequest );
    // delete reqs because finish comm
    delete[] reqs;

    std::map<size_type, std::vector<std::pair<rank_type, size_type>>> mapActiveEltDuplicatedInWorld;
    // treat second recv and build ghost elements
    for ( auto const& dataToRecvPair : dataToRecv2 )
    {
        rank_type rankRecv = dataToRecvPair.first;
        auto itFindRankInMemory = memoryDataToSend.find( rankRecv );
        CHECK( itFindRankInMemory != memoryDataToSend.end() ) << "missing info in memory";
        for ( size_type k = 0; k < dataToRecvPair.second.size(); ++k )
        {
            size_type idEltActiveInOtherProc = dataToRecvPair.second[k];
            if ( idEltActiveInOtherProc == invalid_size_type_value )
                continue;

            auto const& memoryPair = memoryDataToSend[rankRecv][k];
            const size_type oldEltId = memoryPair.first;

            auto const& oldElem = this->entityExtracted( oldEltId, rankRecv, mpl::int_<RangeType>() );
            CHECK( oldElem.id() == oldEltId ) << "invalid id";

            auto itFindNewActiveElt = new_element_id.find( oldEltId );
            bool is_not_stored = itFindNewActiveElt == new_element_id.end();
#if 0
            LOG_IF( WARNING, is_not_stored == false )
                << "this element is already present on the new mesh\n";
            if ( is_not_stored == false ) continue;
#endif
            if ( is_not_stored )
            {
                // create a new elem with partitioning infos
                CHECK( rankRecv != oldElem.pidInPartition() && proc_id == oldElem.pidInPartition() ) << "invalid rank id";
                element_type newElem;
                newElem.setMarker( oldElem.marker().value() );
                newElem.setMarker2( oldElem.marker2().value() );
                newElem.setMarker3( oldElem.marker3().value() );
                newElem.setProcessIdInPartition( proc_id );
                newElem.setProcessId( rankRecv );
                newElem.addNeighborPartitionId( rankRecv );

                // Loop over the nodes on this element.
                for ( uint16_type n = 0; n < oldElem.nPoints(); n++ )
                {
                    auto const& oldPoint = oldElem.point( n );
                    size_type oldPointId = oldPoint.id();
                    size_type newPtId = invalid_size_type_value;
                    auto itFindPoint = new_node_numbers.find( oldPointId );
                    if ( itFindPoint != new_node_numbers.end() )
                        newPtId = itFindPoint->second;
                    else //if ( itFindPoint == new_node_numbers.end() )
                    {
                        DVLOG( 2 ) << "[Mesh<Shape,T>::CreateSubmesh] insert point " << oldPoint << "\n";
                        newPtId = ( renumberPoint ) ? n_new_nodes++ : oldPointId;
                        new_node_numbers[oldPointId] = newPtId;
                        // create point and keep default process id because point not in partition ( invalid_rank_value_type )
                        point_type pt( newPtId, oldPoint );
                        pt.setProcessIdInPartition( proc_id );
                        pt.setMarker( oldPoint.marker().value() );
                        pt.setMarker2( oldPoint.marker2().value() );
                        pt.setMarker3( oldPoint.marker3().value() );
                        // Add this node to the new mesh
                        newMesh->addPoint( pt );
                        DVLOG( 2 ) << "[Mesh<Shape,T>::CreateSubmesh] number of  points " << newMesh->numPoints() << "\n";
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
                        CHECK( newPtId < newMesh->numPoints() ) << "invalid connectivity";

                    newElem.setPoint( n, newMesh->point( newPtId ) );
                    // if ( RangeType == MESH_EDGES )
                    //    newElem.setFace( n, newMesh->point( newPtId ) );

                    DVLOG( 2 ) << "[Mesh<Shape,T>::CreateSubmesh] adding point old(" << oldPointId
                               << ") as point new(" << newPtId << ") in element " << newElem.id() << "\n";

                } // for (uint16_type n=0 ... )

                // update id in other part
                newElem.setIdInOtherPartitions( rankRecv, idEltActiveInOtherProc );

                // Add an equivalent element type to the new_mesh
                auto const& e = newMesh->addElement( newElem, true );

                const size_type newEltId = e.id();

                // mesh relation not yet activated for edges
                if ( RangeType != MESH_EDGES )
                    M_smd->bm.insert( typename smd_type::bm_type::value_type( newEltId, oldEltId ) );
            }
            else // already stored as active element
            {
                size_type newEltId = itFindNewActiveElt->second;
                mapActiveEltDuplicatedInWorld[newEltId].push_back( std::make_pair( rankRecv, idEltActiveInOtherProc ) );
            }

        } // for ( size_type k )
    }     // for ( auto dataRecv2 )

    // maybe some active elements are duplicated in parallel mesh at interprocess zone,
    // we consider an unique active element with the minimal pid
    for ( auto const& dataEltDuplicated : mapActiveEltDuplicatedInWorld )
    {
        size_type newId = dataEltDuplicated.first;
        auto eltIt = newMesh->elementIterator( newId, proc_id );

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

        newMesh->elements().modify( eltIt, [&minPid, &allpid, &dataEltDuplicated]( element_type& e ) {
            for ( rank_type opid : allpid )
                e.addNeighborPartitionId( opid );
            for ( auto const& dataOtherProc : dataEltDuplicated.second )
                e.setIdInOtherPartitions( dataOtherProc.first, dataOtherProc.second );
            e.setProcessId( minPid );
        } );
    }
}

/**
 * @addtogroup FreeFunctions
 * @{
 */
template <typename MeshType, typename IteratorRange, int TheTag = MeshType::tag>
CreateSubmeshTool<MeshType, typename Feel::detail::submeshrangetype<IteratorRange>::type, TheTag>
createSubmeshTool( boost::shared_ptr<MeshType> inputMesh,
                   IteratorRange const& range,
                   WorldComm const& wc,
                   size_type updateComponentsMesh = MESH_CHECK | MESH_UPDATE_FACES | MESH_UPDATE_EDGES,
                   bool subMeshIsOnBoundaryFaces = false /*allow to reduce mpi comm*/ )
{
    CreateSubmeshTool<MeshType, typename Feel::detail::submeshrangetype<IteratorRange>::type, TheTag>
        t( inputMesh, range, wc, updateComponentsMesh );
    t.subMeshIsOnBoundaryFaces( subMeshIsOnBoundaryFaces );
    return t;
}

/**
 * @brief create a submesh 
 * @code
 * auto mesh = unitSquare();
 * auto m = CreateSubmesh( mesh, elements(mesh) ); 
 * @endcode
 */
template <typename MeshType, typename IteratorRange, int TheTag = MeshType::tag>
typename CreateSubmeshTool<MeshType, typename Feel::detail::submeshrangetype<IteratorRange>::type, TheTag>::mesh_build_ptrtype
createSubmesh( boost::shared_ptr<MeshType> inputMesh,
               IteratorRange const& range,
               size_type ctx = EXTRACTION_KEEP_MESH_RELATION,
               size_type updateComponentsMesh = MESH_CHECK | MESH_UPDATE_FACES | MESH_UPDATE_EDGES,
               bool subMeshIsOnBoundaryFaces = false /*allow to reduce mpi comm*/ )
{
    auto t = createSubmeshTool<MeshType, IteratorRange, TheTag>( inputMesh, range, inputMesh->worldComm(), updateComponentsMesh );
    t.subMeshIsOnBoundaryFaces( subMeshIsOnBoundaryFaces );
    return t.build( ctx );
}

/**
 * @brief create a submesh from a range of (sub)elements 
 * @code
 * auto mesh = unitSquare();
 * auto m = CreateSubmesh( mesh, elements(mesh), mesh->worldComm() ); 
 * @endcode
 */
template <typename MeshType, typename IteratorRange, int TheTag = MeshType::tag>
typename CreateSubmeshTool<MeshType, typename Feel::detail::submeshrangetype<IteratorRange>::type, TheTag>::mesh_build_ptrtype
createSubmesh( boost::shared_ptr<MeshType> inputMesh,
               IteratorRange const& range,
               WorldComm const& wc,
               size_type ctx = EXTRACTION_KEEP_MESH_RELATION,
               size_type updateComponentsMesh = MESH_CHECK | MESH_UPDATE_FACES | MESH_UPDATE_EDGES,
               bool subMeshIsOnBoundaryFaces = false /*allow to reduce mpi comm*/ )

{
    auto t = createSubmeshTool<MeshType, IteratorRange, TheTag>( inputMesh, range, wc, updateComponentsMesh );
    t.subMeshIsOnBoundaryFaces( subMeshIsOnBoundaryFaces );
    return t.build( ctx );
}

/**
 * @}
 */

} // namespace Feel

#endif // createsubmesh
