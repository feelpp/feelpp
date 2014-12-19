/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@imag.fr>
       Date: 2011-07-21

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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

#ifndef __createsubmesh_H
#define __createsubmesh_H 1

#include <boost/mpl/if.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION >= 105600
#include <boost/phoenix/stl/algorithm.hpp>
#else
#include <boost/spirit/home/phoenix/stl/algorithm.hpp>
#endif
#include <feel/feelmesh/submeshdata.hpp>


namespace Feel
{
template <typename C, typename V, int T> class Mesh;

template <typename MeshType,typename IteratorRange, int TheTag=MeshType::tag>
class createSubmeshTool
{
public :

    typedef IteratorRange range_type;
    typedef typename boost::tuples::template element<0, range_type>::type idim_type;
    typedef typename boost::tuples::template element<1, range_type>::type iterator_type;

    static const uint16_type tag = TheTag;
    typedef MeshType mesh_type;
    typedef typename mesh_type::value_type value_type;

    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef typename mpl::if_<mpl::bool_<mesh_type::shape_type::is_simplex>,
                              mpl::identity< Mesh< Simplex< mesh_type::nDim,mesh_type::nOrder,mesh_type::nRealDim>, value_type, tag > >,
                              mpl::identity< Mesh< Hypercube<mesh_type::nDim,mesh_type::nOrder,mesh_type::nRealDim>, value_type, tag > > >::type::type mesh_elements_type;
    typedef boost::shared_ptr<mesh_elements_type> mesh_elements_ptrtype;

    typedef typename mpl::if_<mpl::bool_<mesh_type::shape_type::is_simplex>,
                              mpl::identity< Mesh< Simplex< mesh_type::nDim-1,mesh_type::nOrder,mesh_type::nRealDim>, value_type, tag > >,
                              mpl::identity< Mesh< Hypercube<mesh_type::nDim-1,mesh_type::nOrder,mesh_type::nRealDim>, value_type, tag > > >::type::type mesh_faces_type;
    typedef boost::shared_ptr<mesh_faces_type> mesh_faces_ptrtype;

    typedef typename mpl::if_<mpl::bool_<mesh_type::shape_type::is_simplex>,
                              mpl::identity< Mesh< Simplex< (mesh_type::nDim==3)?mesh_type::nDim-2:mesh_type::nDim-1,mesh_type::nOrder,mesh_type::nRealDim>, value_type, tag > >,
                              mpl::identity< Mesh< Hypercube<(mesh_type::nDim==3)?mesh_type::nDim-2:mesh_type::nDim-1,mesh_type::nOrder,mesh_type::nRealDim>, value_type, tag > > >::type::type mesh_edges_type;
    typedef boost::shared_ptr<mesh_edges_type> mesh_edges_ptrtype;

    typedef typename mpl::if_< mpl::equal_to< idim_type ,mpl::size_t<MESH_ELEMENTS> >,
                               mpl::identity<mesh_elements_type>,
                               typename mpl::if_< mpl::equal_to< idim_type ,mpl::size_t<MESH_FACES> >,
                                                  mpl::identity<mesh_faces_type>,
                                                  mpl::identity<mesh_edges_type> >::type>::type::type mesh_build_type;

    typedef boost::shared_ptr<mesh_build_type> mesh_build_ptrtype;
    typedef SubMeshData smd_type;
    typedef boost::shared_ptr<smd_type> smd_ptrtype;

    createSubmeshTool( boost::shared_ptr<MeshType> inputMesh,IteratorRange const& range, WorldComm const& wc, size_type updateComponentsMesh  )
        :
        M_mesh( inputMesh ),
        M_listRange(),
        M_smd( new smd_type( inputMesh ) ),
        M_worldComm( wc ),
        M_updateComponentsMesh( updateComponentsMesh ),
        M_subMeshIsOnBoundaryFaces( true )
        {
            M_listRange.push_back( range );
        }

    createSubmeshTool( boost::shared_ptr<MeshType> inputMesh,std::list<IteratorRange> const& range, WorldComm const& wc, size_type updateComponentsMesh  )
        :
        M_mesh( inputMesh ),
        M_listRange( range ),
        M_smd( new smd_type( inputMesh ) ),
        M_worldComm( wc ),
        M_updateComponentsMesh( updateComponentsMesh ),
        M_subMeshIsOnBoundaryFaces( true )
        {}


    mesh_build_ptrtype
    build()
        {
            DVLOG(2) << "[createSubmeshTool] extracting mesh\n";
            return build( mpl::int_<idim_type::value>() );
        }

    smd_ptrtype subMeshData() { return M_smd; }

    size_type updateComponentsMesh() const { return M_updateComponentsMesh; }

    bool subMeshIsOnBoundaryFaces() const { return M_subMeshIsOnBoundaryFaces; }
    void subMeshIsOnBoundaryFaces( bool b ) { M_subMeshIsOnBoundaryFaces=b; }

private:

    mesh_elements_ptrtype build( mpl::int_<MESH_ELEMENTS> /**/ );
    mesh_faces_ptrtype build( mpl::int_<MESH_FACES> /**/ );
    mesh_edges_ptrtype build( mpl::int_<MESH_EDGES> /**/ );

    mesh_ptrtype M_mesh;
    std::list<range_type> M_listRange;
    smd_ptrtype M_smd;
    WorldComm M_worldComm;
    size_type M_updateComponentsMesh;
    bool M_subMeshIsOnBoundaryFaces;
};

template <typename MeshType,typename IteratorRange,int TheTag>
typename createSubmeshTool<MeshType,IteratorRange,TheTag>::mesh_elements_ptrtype
createSubmeshTool<MeshType,IteratorRange,TheTag>::build( mpl::int_<MESH_ELEMENTS> /**/ )
{
    typedef typename mesh_type::element_type element_type;
    typedef typename mesh_type::point_type point_type;
    typedef typename mesh_type::face_type face_type;

    mesh_elements_ptrtype newMesh( new mesh_elements_type(M_worldComm));

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
    std::map<size_type, size_type> new_vertex;

    Context c( this->updateComponentsMesh() );
    const bool renumberPoint = c.test( MESH_RENUMBER );

    // the number of nodes on the new mesh, will be incremented
    unsigned int n_new_nodes = 0;
    unsigned int n_new_elem  = 0;
    size_type n_new_faces = 0;

    const int proc_id = M_worldComm.localRank();
    const int nProc = M_worldComm.localSize();

    //-----------------------------------------------------------//

    std::map<size_type,size_type> new_element_id;

    std::map<int,std::set<boost::tuple<size_type,size_type> > > ghostCellsFind;


    for (auto& itList : M_listRange)
    {
        auto it = itList.template get<1>();
        auto const en = itList.template get<2>();
        for ( ; it != en; ++ it )
        {
            element_type const& oldElem = *it;
            VLOG(2) << "create sub mesh element from "  << oldElem.id() << "\n";google::FlushLogFiles(google::GLOG_INFO);
            // copy element so that we can modify it
            element_type newElem = oldElem;

            // reset partitioning data

            newElem.setNumberOfPartitions( 1 );
            newElem.setProcessIdInPartition( proc_id );
            newElem.setProcessId( proc_id );
            newElem.clearIdInOthersPartitions();
            newElem.clearNeighborPartitionIds();

            // Loop over the nodes on this element.
            for ( uint16_type n=0; n < oldElem.nPoints(); n++ )
            {
                //FEELPP_ASSERT (oldElem.point( n ).id() < new_node_numbers.size()).error( "invalid point id()" );
                auto const& old_point = oldElem.point( n );

                if ( new_node_numbers.find( old_point.id() ) == new_node_numbers.end() )
                {
                    const size_type newPtId = (renumberPoint)? n_new_nodes : old_point.id();
                    new_node_numbers[old_point.id()] = newPtId;

                    DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] insert point " << oldElem.point( n ) << "\n";

                    point_type pt( old_point );
                    pt.setId( newPtId );
                    pt.setProcessIdInPartition( proc_id );
                    pt.setProcessId( proc_id );
                    pt.clearElementsGhost();
                    pt.clearNeighborPartitionIds();
                    pt.setOnBoundary(false);

                    // Add this node to the new mesh
                    newMesh->addPoint ( pt );

                    DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] number of  points " << newMesh->numPoints() << "\n";

                    // Increment the new node counter
                    n_new_nodes++;

                    if ( n < element_type::numVertices )
                    {
                        CHECK( new_vertex.find(old_point.id()) == new_vertex.end() ) << "already seen this point?";
                        new_vertex[old_point.id()]=1;
                    }

                    if ( old_point.numberOfProcGhost() > 0 && nProc > 1 )
                    {
                        auto const enprocghost=old_point.elementsGhost().end();
                        for (auto& itProcGhost : old_point.elementsGhost())
                        {
                            const int procIdGhost=itProcGhost.first;
                            for (auto& itEltGhost : itProcGhost.second)
                            {
                                auto const& ghostElt = M_mesh->element( itEltGhost,procIdGhost );
                                ghostCellsFind[procIdGhost].insert( boost::make_tuple( ghostElt.id(),
                                                                                       ghostElt.idInOthersPartitions( ghostElt.processId() ) ) );
                            }
                        }
                    }

                } // if ( new_node_numbers.find( old_point.id() ) == new_node_numbers.end() )

                // Define this element's connectivity on the new mesh
                if ( renumberPoint )
                    CHECK ( new_node_numbers[old_point.id()] < newMesh->numPoints() ) <<  "invalid connectivity";

                DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] adding point old(" << old_point.id()
                         << ") as point new(" << new_node_numbers[old_point.id()]
                         << ") in element " << newElem.id() << "\n";

                newElem.setPoint( n, newMesh->point( new_node_numbers[old_point.id()] ) );

            } // for (unsigned int n=0 ... )

            // set id of element
            newElem.setId ( n_new_elem );

            // increment the new element counter
            n_new_elem++;

            // Add an equivalent element type to the new_mesh
            auto const& e = newMesh->addElement( newElem );
            new_element_id[oldElem.id()]= e.id();
            M_smd->bm.insert( typename smd_type::bm_type::value_type( e.id(), oldElem.id() ) );

            // Maybe add faces for this element
            for ( unsigned int s=0; s<oldElem.numTopologicalFaces; s++ )
            {
                if ( !oldElem.facePtr( s ) ) continue;
                // only add face on the boundary: they have some data
                // (boundary ids) which cannot be retrieved otherwise
                //if ( oldElem.neighbor(s) == invalid_size_type_value )
                size_type global_face_id = oldElem.face( s ).id();

                if ( M_mesh->hasFace( global_face_id ) )
                {

                    // get the corresponding face
                    face_type const& old_face = oldElem.face( s );
                    face_type new_face = old_face;

                    // disconnect from elements of old mesh,
                    // the connection will be redone in
                    // \c updateForUse()
                    new_face.disconnect();

                    // very important! updateForUse put false for internalfaces after
                    new_face.setOnBoundary( true );

                    // update points info
                    for ( uint16_type p = 0; p < new_face.nPoints(); ++p )
                    {
                        //new_face.setPoint( p, newMesh->point( new_node_numbers[oldElem.point( oldElem.fToP( s,p ) ).id()] ) );
                        new_face.setPoint( p, newMesh->point( new_node_numbers[ old_face.point(p).id()] ) );
                    }

                    new_face.setId( n_new_faces++ );

                    // reset partitioning data
                    new_face.setNumberOfPartitions( 1 );
                    new_face.setProcessIdInPartition( proc_id );
                    new_face.setProcessId( proc_id );
                    new_face.clearIdInOthersPartitions();
                    new_face.clearNeighborPartitionIds();

                    // add it to the list of faces
                    auto addFaceRes = newMesh->addFace( new_face );
                }

            } // for (unsigned int s=0 ... )

        } //  for( ; it != en; ++ it )
    } // for (auto& itList : M_listRange)

    if ( nProc > 1 )
    {
        DVLOG(1) << "build parallel submesh...";
        auto const theWorldCommSize = newMesh->worldComm().localComm().size();
        std::vector<int> nbMsgToSend( theWorldCommSize , 0 );
        std::vector< std::map<int,size_type> > mapMsg( theWorldCommSize );

        auto itGhostFind = ghostCellsFind.begin();
        auto const enGhostFind = ghostCellsFind.end();
        for ( ; itGhostFind!=enGhostFind ; ++itGhostFind )
        {
            auto const realProcId = itGhostFind->first;
            auto itIdElt = itGhostFind->second.begin();
            auto const enIdElt = itGhostFind->second.end();
            for ( ; itIdElt!=enIdElt ; ++itIdElt)
            {
                auto const idEltInMyProc =itIdElt->template get<0>();
                auto const idEltInRealProc =itIdElt->template get<1>();
                newMesh->worldComm().localComm().send(realProcId, nbMsgToSend[realProcId], idEltInRealProc);
                mapMsg[realProcId].insert( std::make_pair( nbMsgToSend[realProcId],idEltInMyProc ) );
                ++nbMsgToSend[realProcId];
            }
        }

        // counter of msg received for each process
        std::vector<int> nbMsgToRecv;
        mpi::all_to_all( newMesh->worldComm().localComm(),
                         nbMsgToSend,
                         nbMsgToRecv );

        // recv dof asked and re-send dof in this proc
        for ( int proc=0; proc<theWorldCommSize; ++proc )
        {
            for ( int cpt=0; cpt<nbMsgToRecv[proc]; ++cpt )
            {
                //recv
                size_type idEltRecv;
                newMesh->worldComm().localComm().recv( proc, cpt, idEltRecv );
                // search id
                size_type idEltInNewMesh = invalid_size_type_value;
                auto const itFindId = new_element_id.find(idEltRecv);
                if ( itFindId != new_element_id.end() )
                {
                    idEltInNewMesh = itFindId->second;
                    auto eltToUpdate = newMesh->elementIterator( idEltInNewMesh );
                    newMesh->elements().modify( eltToUpdate, Feel::detail::UpdateNeighborPartition( proc ) );
                    //newMesh->elements().modify( eltToUpdate, Feel::detail::updateIdInOthersPartitions( proc, idEltAsked ) );
                }
                // send response
                newMesh->worldComm().localComm().send( proc, cpt, idEltInNewMesh );
            }
        }

        // get response to initial request and update Feel::Mesh::Faces data
        for ( int proc=0; proc<theWorldCommSize; ++proc )
        {
            for ( int cpt=0; cpt<nbMsgToSend[proc]; ++cpt )
            {
                size_type idEltAsked;
                newMesh->worldComm().localComm().recv( proc, cpt, idEltAsked );

                if (idEltAsked != invalid_size_type_value)
                {
                    element_type const& oldElem = M_mesh->element( mapMsg[proc][cpt], proc );

                    //if (new_element_id.find(oldElem.id())!=new_element_id.end() ) continue;
                    CHECK( new_element_id.find(oldElem.id())==new_element_id.end() ) << "this element is present on the new mesh\n";
                    //if (new_element_id.find(oldElem.id())==new_element_id.end() )
                    //{
                    // copy element so that we can modify it
                    element_type newElem = oldElem;

                    // partitioning update
                    newElem.setProcessIdInPartition( oldElem.pidInPartition() );
                    newElem.setNumberOfPartitions(2);
                    CHECK( proc==oldElem.processId() ) << "invalid process id\n";
                    newElem.setProcessId(oldElem.processId());
                    newElem.clearIdInOthersPartitions();

                    std::vector<rank_type> newNeighborPartitionIds(1);
                    newNeighborPartitionIds[0]=newMesh->worldComm().localRank();
                    newElem.setNeighborPartitionIds( newNeighborPartitionIds );

                    // Loop over the nodes on this element.
                    for ( unsigned int n=0; n < oldElem.nPoints(); n++ )
                    {
                        auto const& oldPoint = oldElem.point( n );
                        if ( new_node_numbers.find( oldPoint.id() ) == new_node_numbers.end() )
                        {
                            const size_type newPtId = (renumberPoint)? n_new_nodes : oldPoint.id();
                            new_node_numbers[oldPoint.id()] = newPtId;

                            DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] insert point " << oldElem.point( n ) << "\n";

                            point_type pt( oldElem.point( n ) );
                            pt.setId( newPtId );
                            pt.clearElementsGhost();
                            pt.setProcessIdInPartition( proc_id );
                            pt.setProcessId( invalid_uint16_type_value );
                            pt.clearNeighborPartitionIds();
                            pt.setOnBoundary(false);

                            // Add this node to the new mesh
                            newMesh->addPoint ( pt );

                            DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] number of  points " << newMesh->numPoints() << "\n";

                            // Increment the new node counter
                            n_new_nodes++;

                            if ( n < element_type::numVertices )
                            {
                                CHECK( new_vertex.find(oldPoint.id()) == new_vertex.end() ) << "already seen this point?";
                                new_vertex[oldPoint.id()]=1;
                            }
                        }
                        else if ( newMesh->point( new_node_numbers[oldPoint.id()] ).processId() != invalid_uint16_type_value  )
                        {
                            // update NeighborPartition for this points
                            auto ptToUpdate = newMesh->pointIterator( new_node_numbers[oldPoint.id()] );
                            newMesh->points().modify( ptToUpdate, Feel::detail::UpdateNeighborPartition( proc ) );
                        }

                        // Define this element's connectivity on the new mesh
                        if ( renumberPoint )
                            CHECK ( new_node_numbers[oldPoint.id()] < newMesh->numPoints() ) << "invalid connectivity";

                        DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] adding point old(" << oldPoint.id()
                                 << ") as point new(" << new_node_numbers[oldPoint.id()]
                                 << ") in element " << newElem.id() << "\n";

                        newElem.setPoint( n, newMesh->point( new_node_numbers[oldPoint.id()] ) );

                    } // for (unsigned int n=0 ... )

                    // set id of element
                    newElem.setId ( n_new_elem );

                    // increment the new element counter
                    n_new_elem++;
                    // Add an equivalent element type to the new_mesh
                    auto const& e = newMesh->addElement( newElem );
                    new_element_id[oldElem.id()]= e.id();

                    M_smd->bm.insert( typename smd_type::bm_type::value_type( e.id(), oldElem.id() ) );

                    // update id elt in other partition
                    auto eltToUpdate = newMesh->elementIterator( e.id(),  proc );
                    newMesh->elements().modify( eltToUpdate, Feel::detail::updateIdInOthersPartitions( proc, idEltAsked ) );
                    //} // if (new_element_id.find(oldElem.id())==new_element_id.end() )

                } //  if (idEltAsked != invalid_size_type_value)

            } // for ( int cpt=0; cpt<nbMsgToSend[proc]; ++cpt )

        } // for ( int proc=0; proc<theWorldCommSize; ++proc )

    } // if ( nProc > 1 )


    VLOG(2) << "submesh created\n";google::FlushLogFiles(google::GLOG_INFO);
    newMesh->setNumVertices( std::accumulate( new_vertex.begin(), new_vertex.end(), 0,
                                              []( int lhs, std::pair<int,int> const& rhs )
                                              {
                                                  return lhs+rhs.second;
                                              } ) );

    VLOG(2) << "[Mesh<Shape,T>::createSubmesh] update face/edge info if necessary\n";google::FlushLogFiles(google::GLOG_INFO);
    // Prepare the new_mesh for use
    newMesh->components().reset();
    newMesh->components().set ( this->updateComponentsMesh()/*MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK*/ );
    newMesh->updateForUse();
    VLOG(2) << "createSubmesh(MESH_ELEMENTS) done\n";google::FlushLogFiles(google::GLOG_INFO);

    return newMesh;
}


/**
 * create subMesh from a range<MESH_FACES>
 */
template <typename MeshType,typename IteratorRange,int TheTag>
typename createSubmeshTool<MeshType,IteratorRange,TheTag>::mesh_faces_ptrtype
createSubmeshTool<MeshType,IteratorRange,TheTag>::build( mpl::int_<MESH_FACES> /**/ )
{
    DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] creating new mesh" << "\n";
    mesh_faces_ptrtype newMesh( new mesh_faces_type( M_worldComm) );

    newMesh->setSubStructuring(M_mesh->subStructuring());

    //-----------------------------------------------------------//
    DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] extraction mesh faces" << "\n";
    // inherit the table of markersName
    BOOST_FOREACH( auto itMark, M_mesh->markerNames() )
    {
        DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] adding marker " << itMark.first <<"\n";
        newMesh->addMarkerName( itMark.first,itMark.second[0],itMark.second[1] );
    }

    //-----------------------------------------------------------//

    typedef typename mesh_type::face_type face_type;
    typedef typename mesh_faces_type::element_type new_element_type;
    typedef typename mesh_faces_type::face_type new_face_type;

    std::map<size_type,size_type> new_node_numbers;
    std::map<size_type,size_type> new_vertex;

    Context c( this->updateComponentsMesh() );
    const bool renumberPoint = c.test( MESH_RENUMBER );

    // the number of nodes on the new mesh, will be incremented
    unsigned int n_new_nodes = 0;
    unsigned int n_new_elem  = 0;
    size_type n_new_faces = 0;

    const int proc_id = newMesh->worldComm().localRank();
    const int nProc = newMesh->worldComm().localSize();
    std::map<size_type,size_type> new_element_id;
    std::map<int,std::set<boost::tuple<size_type,size_type> > > ghostCellsFind;

    //-----------------------------------------------------------//

    for (auto& itList : M_listRange)
    {
        auto it = itList.template get<1>();
        auto const en = itList.template get<2>();

        DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] extracting " << std::distance(it,en)  << " faces " << "\n";
        for ( ; it != en; ++ it )
        {
            // create a new element
            face_type const& oldElem = *it;
            DVLOG(2) << "[Mesh<Shape,T>::createSubmesh]   + face : " << oldElem.id() << "\n";

            if ( nProc > 1 && this->subMeshIsOnBoundaryFaces() )
                CHECK( oldElem.isOnBoundary() ) << "error : use mpi optimzation subMeshIsOnBoundaryFaces but an internal face is added";

            // copy element so that we can modify it
            new_element_type newElem;// = oldElem;

            // get element markers
            newElem.setMarker( oldElem.marker().value() );
            newElem.setMarker2( oldElem.marker2().value() );
            newElem.setMarker3( oldElem.marker3().value() );

            //CHECK( !oldElem.isGhostCell() ) << "only actif elt\n";
            // reset partitioning data
            newElem.setProcessIdInPartition( proc_id );
            newElem.setNumberOfPartitions( 1 );
            newElem.setProcessId( proc_id );
            newElem.clearIdInOthersPartitions();
            newElem.clearNeighborPartitionIds();

            // Loop over the nodes on this element.
            for ( unsigned int n=0; n < oldElem.nPoints(); n++ )
            {
                auto const& oldPoint = oldElem.point( n );

                if ( new_node_numbers.find(oldPoint.id()) == new_node_numbers.end() )
                {
                    const size_type newPtId = (renumberPoint)? n_new_nodes : oldPoint.id();
                    new_node_numbers[oldPoint.id()] = newPtId;

                    DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] insert point " << oldPoint << "\n";

                    typename mesh_faces_type::point_type pt( oldPoint );
                    pt.setId( newPtId );
                    pt.clearElementsGhost();
                    pt.setProcessIdInPartition( proc_id );
                    pt.setProcessId( proc_id );
                    pt.clearNeighborPartitionIds();
                    pt.setOnBoundary(false);

                    // Add this node to the new mesh
                    newMesh->addPoint( pt );

                    DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] number of  points " << newMesh->numPoints() << "\n";

                    // Increment the new node counter
                    n_new_nodes++;

                    if ( n < new_element_type::numVertices )
                    {
                        CHECK( new_vertex.find(oldPoint.id()) == new_vertex.end() ) << "already seen this point?";
                        new_vertex[oldPoint.id()]=1;
                    }

                    if ( nProc > 1 && oldPoint.numberOfProcGhost() > 0 )
                    {
                        auto itprocghost = oldPoint.elementsGhost().begin();
                        auto const enprocghost = oldPoint.elementsGhost().end();
                        for ( ; itprocghost!=enprocghost ; ++itprocghost )
                        {
                            const int procIdGhost=itprocghost->first;
                            auto iteltghost = itprocghost->second.begin();
                            auto const eneltghost = itprocghost->second.end();
                            for ( ; iteltghost!=eneltghost ; ++iteltghost )
                            {
                                auto const eltIdGhost = *iteltghost;
                                auto const& ghostElt = M_mesh->element(eltIdGhost,procIdGhost);
                                for ( unsigned int s=0; s<ghostElt.numTopologicalFaces; s++ )
                                {
                                    auto const& ghostFace = ghostElt.face( s );
                                    // reduce mpi comm if possible
                                    if ( this->subMeshIsOnBoundaryFaces() && !ghostFace.isOnBoundary() )
                                        continue;
                                    // store ghost faces to find in other process
                                    ghostCellsFind[procIdGhost].insert(boost::make_tuple( ghostFace.id(),
                                                                                          ghostFace.idInOthersPartitions(ghostElt.processId())) );
                                }
                            }
                        }
                    } // if (oldPoint.numberOfProcGhost()>0)

                }

                newElem.setPoint( n, newMesh->point( new_node_numbers[oldPoint.id()] ) );

            } // end for n

            // set id of element
            newElem.setId ( n_new_elem );

            // increment the new element counter
            n_new_elem++;

            // Add an equivalent element type to the new_mesh
            auto const& e = newMesh->addElement( newElem );
            new_element_id[oldElem.id()]= e.id();
            M_smd->bm.insert( typename smd_type::bm_type::value_type( e.id(), oldElem.id() ) );


        } // end for it
    } // for (auto& itList : M_listRange)


    if ( nProc > 1 )
    {
        auto const theWorldCommSize = newMesh->worldComm().localComm().size();
        std::vector<int> nbMsgToSend( theWorldCommSize , 0 );
        std::vector< std::map<int,size_type> > mapMsg( theWorldCommSize );

        for (auto& itGhostFind : ghostCellsFind)
        {
            auto const realProcId = itGhostFind.first;
            for (auto& itIdElt : itGhostFind.second)
            {
                auto const idEltInMyProc =itIdElt.template get<0>();
                auto const idEltInRealProc =itIdElt.template get<1>();
                newMesh->worldComm().localComm().send(realProcId, nbMsgToSend[realProcId], idEltInRealProc);
                mapMsg[realProcId].insert( std::make_pair( nbMsgToSend[realProcId],idEltInMyProc ) );
                ++nbMsgToSend[realProcId];
            }
        }

        // counter of msg received for each process
        std::vector<int> nbMsgToRecv;
        mpi::all_to_all( newMesh->worldComm().localComm(),
                         nbMsgToSend,
                         nbMsgToRecv );


        // recv dof asked and re-send dof in this proc
        for ( int proc=0; proc<theWorldCommSize; ++proc )
        {
            for ( int cpt=0; cpt<nbMsgToRecv[proc]; ++cpt )
            {
                //recv
                size_type idEltRecv;
                newMesh->worldComm().localComm().recv( proc, cpt, idEltRecv );
                // search id
                size_type idEltInNewMesh = invalid_size_type_value;
                auto const itFindId = new_element_id.find(idEltRecv);
                if ( itFindId != new_element_id.end() )
                {
                    idEltInNewMesh = itFindId->second;
                    // update NeighborPartition for this elt
                    auto eltToUpdate = newMesh->elementIterator( idEltInNewMesh );
                    newMesh->elements().modify( eltToUpdate, Feel::detail::UpdateNeighborPartition( proc ) );
                }
                // send response
                newMesh->worldComm().localComm().send( proc, cpt, idEltInNewMesh );
            }
        }

        // get response to initial request and update Feel::Mesh::Faces data
        for ( int proc=0; proc<theWorldCommSize; ++proc )
        {
            for ( int cpt=0; cpt<nbMsgToSend[proc]; ++cpt )
            {
                size_type idEltAsked;
                newMesh->worldComm().localComm().recv( proc, cpt, idEltAsked );

                if (idEltAsked != invalid_size_type_value)
                {
                    auto const& old_elem = M_mesh->face( mapMsg[proc][cpt] );

                    //CHECK( new_element_id.find(old_elem.id())==new_element_id.end() ) << "this element is already present on the new mesh\n";

                    // create a new elem
                    new_element_type newElem;
                    // partitioning update
                    newElem.setProcessIdInPartition( old_elem.pidInPartition() );
                    newElem.setNumberOfPartitions( 2 );
                    //CHECK( proc==old_elem.processId() ) << "invalid process id\n";
                    newElem.setProcessId( proc );
                    //newElem.clearIdInOthersPartitions();
                    std::vector<rank_type> newNeighborPartitionIds(1);
                    newNeighborPartitionIds[0]=newMesh->worldComm().localRank();
                    newElem.setNeighborPartitionIds( newNeighborPartitionIds );

                    // Loop over the nodes on this element.
                    for ( unsigned int n=0; n < old_elem.nPoints(); n++ )
                    {
                        auto const& oldPoint = old_elem.point( n );
                        if ( new_node_numbers.find(oldPoint.id()) == new_node_numbers.end() )
                        {
                            const size_type newPtId = (renumberPoint)? n_new_nodes : oldPoint.id();
                            new_node_numbers[oldPoint.id()] = newPtId;

                            DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] insert point " << oldPoint << "\n";

                            typename mesh_faces_type::point_type pt( oldPoint );
                            pt.setId( newPtId );
                            pt.setProcessIdInPartition( proc_id );
                            pt.setProcessId( invalid_uint16_type_value );
                            pt.clearElementsGhost();
                            pt.clearNeighborPartitionIds();
                            pt.setOnBoundary(false);

                            // Add this node to the new mesh
                            newMesh->addPoint ( pt );

                            DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] number of  points " << newMesh->numPoints() << "\n";

                            // Increment the new node counter
                            n_new_nodes++;

                            if ( n < new_element_type::numVertices )
                            {
                                CHECK( new_vertex.find(oldPoint.id()) == new_vertex.end() ) << "already seen this point?";
                                new_vertex[oldPoint.id()]=1;
                            }
                        }
                        else if ( newMesh->point( new_node_numbers[oldPoint.id()] ).processId() != invalid_uint16_type_value  )
                        {
                            // update NeighborPartition for this points
                            auto ptToUpdate = newMesh->pointIterator( new_node_numbers[oldPoint.id()] );
                            newMesh->points().modify( ptToUpdate, Feel::detail::UpdateNeighborPartition( proc ) );
                        }

                        // Define this element's connectivity on the new mesh
                        if ( renumberPoint )
                            CHECK ( new_node_numbers[old_elem.point( n ).id()] < newMesh->numPoints() ) << "invalid connectivity";

                        DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] adding point old(" << old_elem.point( n ).id()
                                 << ") as point new(" << new_node_numbers[old_elem.point( n ).id()]
                                 << ") in element " << newElem.id() << "\n";

                        newElem.setPoint( n, newMesh->point( new_node_numbers[oldPoint.id()] ) );

                    } // for (unsigned int n=0 ... )

                    // set id of element
                    newElem.setId ( n_new_elem );

                    // increment the new element counter
                    n_new_elem++;
                    // Add an equivalent element type to the new_mesh
                    auto const& e = newMesh->addElement( newElem );
                    new_element_id[old_elem.id()]= e.id();
                    M_smd->bm.insert( typename smd_type::bm_type::value_type( e.id(), old_elem.id() ) );

                    // save idEltAsked;
                    auto eltToUpdate = newMesh->elementIterator( e.id(),  proc );
                    newMesh->elements().modify( eltToUpdate, Feel::detail::updateIdInOthersPartitions( proc, idEltAsked ) );
                } //  if (idEltAsked != invalid_size_type_value)

            }
        }

    } // if ( nProc > 1 )

    newMesh->setNumVertices( std::accumulate( new_vertex.begin(), new_vertex.end(), 0,
                                              []( int lhs, std::pair<int,int> const& rhs )
                                              {
                                                  return lhs+rhs.second;
                                              } ) );

    DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] update face/edge info if necessary\n";

    // Prepare the new_mesh for use
    newMesh->components().reset();
    newMesh->components().set ( this->updateComponentsMesh()/*MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK*/ );
    newMesh->updateForUse();

    return newMesh;
}

template <typename MeshType,typename IteratorRange,int TheTag>
typename createSubmeshTool<MeshType,IteratorRange,TheTag>::mesh_edges_ptrtype
createSubmeshTool<MeshType,IteratorRange,TheTag>::build( mpl::int_<MESH_EDGES> /**/ )
{
    // we don't deal with this situation yet
    M_smd.reset();

    DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] creating new mesh" << "\n";
    mesh_edges_ptrtype newMesh( new mesh_edges_type( M_worldComm) );
    //mesh_edges_ptrtype newMesh( new mesh_edges_type );

    //-----------------------------------------------------------//
    DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] extraction mesh edges" << "\n";
    // inherit the table of markersName
    BOOST_FOREACH( auto itMark, M_mesh->markerNames() )
    {
        DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] adding marker " << itMark.first <<"\n";
        newMesh->addMarkerName( itMark.first,itMark.second[0],itMark.second[1] );
    }

    //-----------------------------------------------------------//

    typedef typename mesh_edges_type::element_type new_element_type;

    std::map<size_type,size_type> new_node_numbers;
    std::map<size_type,size_type> new_vertex;

    // the number of nodes on the new mesh, will be incremented
    unsigned int n_new_nodes = 0;
    unsigned int n_new_elem  = 0;
    size_type n_new_edges = 0;

    //-----------------------------------------------------------//

    auto itListRange = M_listRange.begin();
    auto const enListRange = M_listRange.end();
    for ( ; itListRange!=enListRange ; ++itListRange)
    {
        auto it = itListRange->template get<1>();
        auto const en = itListRange->template get<2>();

        DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] extracting " << std::distance(it,en)  << " edges " << "\n";
        for ( ; it != en; ++ it )
        {
            auto const& oldElem = boost::unwrap_ref( *it );
            DVLOG(2) << "[Mesh<Shape,T>::createSubmesh]   + face : " << oldElem.id() << "\n";

            // create new element
            new_element_type newElem;// = oldElem;

            // get element markers
            newElem.setMarker( oldElem.marker().value() );
            newElem.setMarker2( oldElem.marker2().value() );
            newElem.setMarker3( oldElem.marker3().value() );

            // partitioning update
            newElem.setProcessIdInPartition( oldElem.pidInPartition() );
            newElem.setNumberOfPartitions(oldElem.numberOfPartitions());
            newElem.setProcessId(oldElem.processId());
            //newElem.setIdInPartition( oldElem.pidInPartition(), n_new_elem );
            newElem.setNeighborPartitionIds(oldElem.neighborPartitionIds());// TODO


            DVLOG(2) << "\n oldElem.nPoints " << oldElem.nPoints() << "\n";
            // Loop over the nodes on this element.
            for ( unsigned int n=0; n < oldElem.nPoints(); n++ )
            {

                if ( new_node_numbers.find(oldElem.point( n ).id()) == new_node_numbers.end() )
                {
                    new_node_numbers[oldElem.point( n ).id()] = n_new_nodes;

                    DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] insert point " << oldElem.point(n) << "\n";

                    typename mesh_edges_type::point_type pt( oldElem.point( n ) );
                    pt.setId( n_new_nodes );

                    // Add this node to the new mesh
                    newMesh->addPoint( pt );

                    DVLOG(2) << "[Mesh<Shape,T>::createSubmesh] number of  points " << newMesh->numPoints() << "\n";

                    // Increment the new node counter
                    n_new_nodes++;

                    if ( n < new_element_type::numVertices )
                    {
                        CHECK( new_vertex.find(oldElem.point( n ).id()) == new_vertex.end() ) << "already seen this point?";
                        new_vertex[oldElem.point( n ).id()]=1;
                    }

                }

                newElem.setPoint( n, newMesh->point( new_node_numbers[oldElem.point( n ).id()] ) );
                newElem.setFace( n, newMesh->point( new_node_numbers[oldElem.point( n ).id()] ) );
            } // end for n
            CHECK( newElem.pointPtr(0) ) << "invalid point 0 in edge";
            CHECK( newElem.pointPtr(1) ) << "invalid point 1 in edge";
            CHECK( newElem.facePtr(0) ) << "invalid face 0 in edge";
            CHECK( newElem.facePtr(1) ) << "invalid face 1 in edge";

            // set id of element
            newElem.setId ( n_new_elem );
            newElem.setProcessId ( oldElem.processId() );

            // increment the new element counter
            n_new_elem++;

            // Add an equivalent element type to the new_mesh
            newMesh->addElement( newElem );
        } // end for it
    } // for ( ; itListRange!=enListRange ; ++itListRange)


    newMesh->setNumVertices( std::accumulate( new_vertex.begin(), new_vertex.end(), 0,
                                              []( int lhs, std::pair<int,int> const& rhs )
                                              {
                                                  return lhs+rhs.second;
                                              } ) );

    DVLOG(2) << "[createSubmesh] update face/edge info if necessary\n";
    // Prepare the new_mesh for use
    newMesh->components().reset();
    newMesh->components().set ( this->updateComponentsMesh()/*MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK*/ );
    newMesh->updateForUse();



    return newMesh;
}

namespace detail
{
template <typename RangeType>
struct submeshrangetype
{
    typedef typename mpl::if_< boost::is_std_list<RangeType>,
                               mpl::identity<RangeType>,
                               mpl::identity<std::list<RangeType> > >::type::type::value_type type;
};
}
template <typename MeshType,typename IteratorRange, int TheTag = MeshType::tag>
typename createSubmeshTool<MeshType,typename Feel::detail::submeshrangetype<IteratorRange>::type,TheTag>::mesh_build_ptrtype
createSubmesh( boost::shared_ptr<MeshType> inputMesh,
               IteratorRange const& range,
               size_type ctx = EXTRACTION_KEEP_MESH_RELATION,
               size_type updateComponentsMesh = MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES,
               bool subMeshIsOnBoundaryFaces = true )
{
    //DVLOG(2) << "[createSubmesh] extracting " << range.template get<0>() << " nb elements :"
    //<< std::distance(range.template get<1>(),range.template get<2>()) << "\n";

    createSubmeshTool<MeshType,typename Feel::detail::submeshrangetype<IteratorRange>::type,TheTag> cSmT( inputMesh,range,inputMesh->worldComm(),updateComponentsMesh );
    cSmT.subMeshIsOnBoundaryFaces( subMeshIsOnBoundaryFaces );
    auto m = cSmT.build();
    Context c( ctx );
    if ( c.test( EXTRACTION_KEEP_MESH_RELATION ) )
        m->setSubMeshData( cSmT.subMeshData() );
    return m;
}

template <typename MeshType,typename IteratorRange, int TheTag = MeshType::tag>
typename createSubmeshTool<MeshType,typename Feel::detail::submeshrangetype<IteratorRange>::type,TheTag>::mesh_build_ptrtype
createSubmesh( boost::shared_ptr<MeshType> inputMesh,
               IteratorRange const& range,
               WorldComm wc,
               size_type ctx = EXTRACTION_KEEP_MESH_RELATION,
               size_type updateComponentsMesh = MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES,
               bool subMeshIsOnBoundaryFaces = true )

{
    //DVLOG(2) << "[createSubmesh] extracting " << range.template get<0>() << " nb elements :"
    //<< std::distance(range.template get<1>(),range.template get<2>()) << "\n";

    createSubmeshTool<MeshType,typename Feel::detail::submeshrangetype<IteratorRange>::type,TheTag> cSmT( inputMesh,range, wc,updateComponentsMesh );
    cSmT.subMeshIsOnBoundaryFaces( subMeshIsOnBoundaryFaces );
    auto m = cSmT.build();
    Context c( ctx );
    if ( c.test( EXTRACTION_KEEP_MESH_RELATION ) )
        m->setSubMeshData( cSmT.subMeshData() );
    return m;
}


} // namespace Feel

#endif // createsubmesh
