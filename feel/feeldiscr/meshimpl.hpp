/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-08-07

  Copyright (C) 2007-2010 Université Joseph Fourier (Grenoble I)

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
   \file meshimpl.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-08-07
 */
#ifndef __MESHIMPL_HPP
#define __MESHIMPL_HPP 1

#include <boost/preprocessor/comparison/greater_equal.hpp>
#include <feel/feelmesh/geoentity.hpp>
#include <feel/feelmesh/simplex.hpp>
#include <feel/feelmesh/hypercube.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelalg/solvernonlinearpetsc.hpp>
#include <feel/feelfilters/gmshenums.hpp>

namespace Feel
{
template<typename Shape, typename T, int Tag>
Mesh<Shape, T, Tag>::Mesh( WorldComm const& worldComm )
    :
    super(worldComm),
    _M_gm( new gm_type ),
    _M_gm1( new gm1_type ),
    M_meas( 0 ),
    M_measbdy( 0 ),
    //M_part(),
    M_tool_localization( new Localization() )
{
    Debug( 4015 ) << "[Mesh] constructor called\n";
}
template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::partition ( const uint16_type n_parts )
{
    //M_part->partition( *this, n_parts );
}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::updateForUse()
{
    Debug( 4015 ) << "component     MESH_RENUMBER: " <<  this->components().test( MESH_RENUMBER ) << "\n";
    Debug( 4015 ) << "component MESH_UPDATE_EDGES: " <<  this->components().test( MESH_UPDATE_EDGES ) << "\n";
    Debug( 4015 ) << "component MESH_UPDATE_FACES: " <<  this->components().test( MESH_UPDATE_FACES ) << "\n";
    Debug( 4015 ) << "component    MESH_PARTITION: " <<  this->components().test( MESH_PARTITION ) << "\n";

    if ( this->numElements() == 0 )
    {
        Debug( 4015 ) << "No elements in Mesh?\n";
        return;
    }

    boost::timer ti;

    Debug( 4015 ) << "is already updated? : " << this->isUpdatedForUse() << "\n";
    if ( !this->isUpdatedForUse() )
    {
        if ( this->components().test( MESH_RENUMBER ) )
        {

            this->renumber();
            Debug( 4015 ) << "[Mesh::updateForUse] renumber : " << ti.elapsed() << "\n";
        }

        //
        // compute the Adjacency graph
        //
        ti.restart();
        Debug( 4015 ) << "Compute adjacency graph\n";
        element_iterator iv,  en;
        std::map<size_type,boost::tuple<size_type, uint16_type, size_type> > f2e;

        _M_e2f.resize( boost::extents[this->numElements()][this->numLocalFaces()] );
        std::map<std::set<int>, size_type > _faces;
        typename std::map<std::set<int>, size_type >::iterator _faceit;
        int next_face = 0;

        boost::tie( iv, en ) = this->elementsRange();

        for ( ; iv != en; ++iv )
        {
            for ( size_type j = 0; j < this->numLocalFaces(); j++ )
            {
                std::set<int> s;

                for ( int f = 0; f < face_type::numVertices; ++f )
                {
                    if ( nDim == 1 )
                        s.insert( iv->point( j ).id() );

                    else
                        s.insert( iv->point( iv->fToP( j, f ) ).id() );
                }

                bool faceinserted = false;
                boost::tie( _faceit, faceinserted ) = _faces.insert( std::make_pair( s, next_face ) );

                if ( faceinserted )
                    ++next_face;

#if !defined ( NDEBUG )
                Debug( 4015 ) << "------------------------------------------------------------\n";
                Debug( 4015 ) << "Element id: " << iv->id() << " local face id: " << j << " process id:" << iv->processId() << "\n";
#endif
                //e = _be.addIfNotThere( baremaker( j ) );


                if ( faceinserted )
                {
#if !defined ( NDEBUG )
                    Debug( 4015 ) << " new face " << _faceit->second << " is now in store with elt " << f2e[_faceit->second].template get<0>()  << " and local face id " <<  f2e[_faceit->second].template get<1>()  << "\n";
#endif

                    f2e[_faceit->second].template get<0>() = iv->id();
                    f2e[_faceit->second].template get<1>() = j;
                    _M_e2f[iv->id()][j]=boost::make_tuple( _faceit->second, invalid_size_type_value );
                }

                else // already stored
                {
#if !defined ( NDEBUG )
                    Debug( 4015 ) << "old face " << _faceit->second << " was already in store with elt " << f2e[_faceit->second].template get<0>() << " and local face id " <<  f2e[_faceit->second].template get<1>() << "\n";
#endif

                    f2e[_faceit->second].template get<2>() = iv->id();
                    _M_e2f[iv->id()][j]=boost::make_tuple( _faceit->second, f2e[_faceit->second].template get<0>() );
                    _M_e2f[f2e[_faceit->second].template get<0>()] [f2e[_faceit->second].template get<1>()] =
                        boost::make_tuple( _faceit->second, iv->id() );
                }


            } // local face
        } // element loop

        Debug( 4015 ) << "Compute adjacency graph done in " << ti.elapsed() << "\n";
#if 0

        // partition mesh
        if ( this->components().test( MESH_PARTITION ) || ( this->worldComm().localSize() > 1 ) )
        {
            boost::timer ti1;
            this->partition();
            Debug( 4015 ) << "[Mesh::updateForUse] partition time : " << ti1.elapsed() << "\n";
        }

#endif

        if ( this->components().test( MESH_UPDATE_FACES ) )
        {
            ti.restart();
            // update connectivities of entities of co dimension 1
            this->updateEntitiesCoDimensionOne();
            // update permutation of entities of co-dimension 1
            this->updateEntitiesCoDimensionOnePermutation();
            // update in ghost cells of entities of co-dimension 1
            //if (this->worldComm().localSize()>1)
            //    this->updateEntitiesCoDimensionOneGhostCell();

            Debug( 4015 ) << "[Mesh::updateForUse] update entities of codimension 1 : " << ti.elapsed() << "\n";

        }

        if ( this->components().test( MESH_UPDATE_EDGES ) )
        {
            ti.restart();
            // update connectivities of entities of co dimension 2
            // (edges in 3D)
            this->updateEntitiesCoDimensionTwo();
            Debug( 4015 ) << "[Mesh::updateForUse] update edges : " << ti.elapsed() << "\n";
        }

        if ( this->components().test( MESH_UPDATE_FACES ) ||
             this->components().test( MESH_UPDATE_EDGES )
            )
        {
            updateOnBoundary( mpl::int_<nDim>() );
        }
        this->setUpdatedForUse( true );
    }

    {
        element_iterator iv,  en;
        boost::tie( iv, en ) = this->elementsRange();

        for ( ; iv != en; ++iv )
        {
            this->elements().modify( iv, typename super_elements::ElementConnectPointToElement() );
        }

        boost::tie( iv, en ) = this->elementsRange();
        auto pc = _M_gm->preCompute( _M_gm, _M_gm->referenceConvex().vertices() );
        auto pcf =  _M_gm->preComputeOnFaces( _M_gm, _M_gm->referenceConvex().barycenterFaces() );
        M_meas = 0;
        M_measbdy = 0;


        for ( ; iv != en; ++iv )
        {
            this->elements().modify( iv,
                                     lambda::bind( &element_type::setMeshAndGm,
                                                   lambda::_1,
                                                   this, _M_gm, _M_gm1 ) );

            this->elements().modify( iv,
                                     lambda::bind( &element_type::updateWithPc,
                                                   lambda::_1, pc, boost::ref( pcf ) ) );

            M_meas += iv->measure();
            auto _faces = iv->faces();

            if ( nDim == 1 )
                M_measbdy = 0;
            else
                for ( ; _faces.first != _faces.second; ++_faces.first )
                    if ( ( *_faces.first ) && ( *_faces.first )->isOnBoundary() )
                        M_measbdy += ( *_faces.first )->measure();
        }

        // now that all elements have been updated, build inter element
        // data such as the measure of point element neighbors
        boost::tie( iv, en ) = this->elementsRange();

        for ( ; iv != en; ++iv )
        {
            value_type meas = 0;
            BOOST_FOREACH( auto _elt, iv->pointElementNeighborIds() )
            {
                if ( this->hasElement( _elt ) )
                    meas += this->element( _elt ).measure();
            }
            this->elements().modify( iv,
                                     lambda::bind( &element_type::setMeasurePointElementNeighbors,
                                                   lambda::_1, meas ) );
        }

        typedef typename super::face_const_iterator face_const_iterator;
        face_iterator itf = this->beginFace();
        face_iterator ite = this->endFace();

        for ( ; itf != ite; ++ itf )
        {
            this->faces().modify( itf,
                                  lambda::bind( &face_type::setMesh,
                                                lambda::_1,
                                                this ) );
        }
    }
    //std::cout<<"this->worldComm().localSize()=     "<< this->worldComm().localSize() << std::endl;
#if defined(FEELPP_ENABLE_MPI_MODE)

    if ( this->components().test( MESH_UPDATE_FACES ) && this->worldComm().localSize()>1 )
    {
        this->updateEntitiesCoDimensionOneGhostCell();
    }
#endif

    propagateMarkers(mpl::int_<nDim>() );

    // check mesh connectivity
    this->check();
    //std::cout<<"pass hier\n";

    _M_gm->initCache( this );
    _M_gm1->initCache( this );

    M_tool_localization->setMesh( this->shared_from_this(),false );



    Debug( 4015 ) << "[Mesh::updateForUse] total time : " << ti.elapsed() << "\n";
}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::propagateMarkers( mpl::int_<3> )
{
    // first propagate top-down  marker from edges if points have not been marked
    std::for_each( this->beginEdge(), this->endEdge(),
                   [this]( edge_type const& e )
                   {
                       if ( e.marker().isOff() )
                           return;

                       for( int i = 0; i < edge_type::numPoints; ++i )
                       {
                           if ( e.point( i ).marker().isOff() )
                           {
                               // inherit marker from edge
                               this->points().modify( this->points().iterator_to( e.point(i) ),
                                                      [&e] ( point_type& p )
                                                      {
                                                          p.setMarker( e.marker().value() );
                                                      } );

                           }
                       } } );
    // then propagate top-down marker from face if edge has not been marked
    std::for_each( this->beginFace(), this->endFace(),
                   [this]( face_type const& f )
                   {
                       if ( f.marker().isOff() )
                           return;

                       // update points
                       for( int i = 0; i < face_type::numPoints; ++i )
                       {
                           if ( f.point( i ).marker().isOff() )
                           {
                               // inherit marker from edge
                               this->points().modify( this->points().iterator_to( f.point(i) ),
                                                      [&f] ( point_type& p )
                                                      {
                                                          p.setMarker( f.marker().value() );
                                                      } );

                           }
                       }
                       // update edges
                       for( int i = 0; i < face_type::numEdges; ++i )
                       {
                           if ( f.edge( i ).marker().isOff() )
                           {
                               // inherit marker from edge
                               this->edges().modify( this->edges().iterator_to( f.edge(i) ),
                                                      [&f] ( edge_type& e )
                                                      {
                                                          e.setMarker( f.marker().value() );
                                                      } );

                           }
                       }
                   } );
}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::renumber( mpl::bool_<true> )
{
    size_type next_free_node = 0;

    // map old/new ids
    std::vector<size_type> node_map( this->numPoints(), invalid_size_type_value );
    //std::map<size_type> node_map( this->numPoints(), invalid_size_type_value );
    //std::vector<point_type> pt_map( this->numPoints() );
    typedef std::map<size_type,point_type> ptmap_type;
    ptmap_type pt_map;

    typedef typename ptmap_type::iterator ptmap_iterator;
    typedef std::vector<size_type>::iterator nm_iterator;

    // first collect all point id that will be swapped in a dictionary
    for ( element_const_iterator elt = this->beginElement();
            elt != this->endElement(); ++elt )
    {
        element_type const& __element = *elt;
#if !defined( NDEBUG )
        Debug( 4015 ) << "mesh::renumber] element id " << __element.id() <<  " proc " << __element.processId() << "\n";

        for ( int i = 0; i < __element.nPoints(); ++i )
        {
            Debug( 4015 ) << "point id = " << __element.point( i ).id() << "\n";
        }

#endif

        // renumber the nodes of the element
        for ( int i = 0; i < __element.nPoints(); ++i )
        {
            point_type  __pt = __element.point( i );
            size_type __true_id = __pt.id();
            size_type __id = __true_id;

            // did we already change the id? if yes get the
            //  id it was changed to
            if ( node_map[ __true_id ] != invalid_size_type_value )
            {
                __id = node_map[__true_id];
#if !defined( NDEBUG )

                if ( __id >= next_free_node )
                    Debug( 4015 ) << "next_free_node = " << next_free_node
                                  << " point id = " << __id << "\n";

#endif
            }

            // don't renumber if already done
            if ( __id == next_free_node )
            {
                node_map[ __true_id ] = next_free_node;
                pt_map[next_free_node] = __pt;
                pt_map[next_free_node].setId( next_free_node );
                ++next_free_node;
            }

            else if ( __id > next_free_node )
            {
                // first check if next_free_node has been
                // already taken by another point, if so
                // then affect this point id to __id
                nm_iterator nm_it = std::find( node_map.begin(),
                                               node_map.end(),
                                               next_free_node );

                if ( nm_it != node_map.end() )
                    *nm_it = __id;

                else
                    node_map[ next_free_node ] = __id;

                // must use \p __true_id here as it is the
                // original point id we started with
                node_map[ __true_id ] = next_free_node;
                pt_map[next_free_node] = __pt;
                pt_map[next_free_node].setId( next_free_node );

#if !defined( NDEBUG )
                Debug( 4015 ) << "next_free_node = " << next_free_node
                              << " swapping point " << __true_id
                              << " and point " << next_free_node << "\n";
#endif

                // next free node
                ++next_free_node;
            } // if
        } // for

        //__element.setId(
    }

    Debug( 4015 ) << "[mesh::renumber] done collecting ids\n";
#if !defined(NDEBUG)
    std::vector<size_type> check_id( node_map );
    std::unique( check_id.begin(), check_id.end() );

    FEELPP_ASSERT( check_id.size() == node_map.size() )( node_map.size() )( check_id.size() ).error( "all ids must be unique" );

    // in parallel this will generate a lot of warning since all the points are
    //loaded into the mesh data structure including the ones not belonging to
    //the current processor. Do the test only in sequential
    if ( Environment::numberOfProcessors() == 1 )
    {
        FEELPP_ASSERT( std::find( node_map.begin(),node_map.end(), invalid_size_type_value ) == node_map.end() ).warn("invalid size_type value found as id " );
    }

#endif /* NDEBUG */

    // now we can replace the points id
    //for( size_type p = 0; p < pt_map.size(); ++p )
    ptmap_iterator ptmapit = pt_map.begin();
    ptmap_iterator ptmapen = pt_map.end();

    for ( ; ptmapit != ptmapen; ++ptmapit )
    {
        point_iterator ptit = this->points().find( ptmapit->first );
#if !defined( NDEBUG )
        Debug( 4015 ) << "[mesh::replace] replacing point " << ptit->id()
                      <<  " with " << ptmapit->second.id() << "\n";
#endif
        bool __rep1 = this->points().replace( ptit, ptmapit->second );
        FEELPP_ASSERT( __rep1 )( __rep1 )( ptit->id() )( ptmapit->second.id() ) .warn( "invalid point replacement" );

    }

    Debug( 4015 ) << "[mesh::renumber] done replace point ids\n";

    for ( element_iterator elt = this->beginElement();
            elt != this->endElement(); ++elt )
    {
        element_type __element = *elt;
#if !defined( NDEBUG )
        Debug( 4015 ) << "mesh::renumber] element id " << __element.id() <<  " proc " << __element.processId() << "\n";
#endif

        // renumber the nodes of the element
        for ( int i = 0; i < __element.nPoints(); ++i )
        {

            size_type __true_id =__element.point( i ).id();
            this->elements().modify( elt,
                                     typename super_elements::ElementUpdatePoint( i, this->point( node_map[__true_id] ) ) );
#if !defined( NDEBUG )
            Debug( 4015 ) << "point id = " << __true_id << " node map " <<  node_map[__true_id] << " "
                          << "new point id = " << elt->point( i ).id() << "\n";
#endif
        }
    }

    Debug( 4015 ) << "[mesh::renumber] done replace point ids in elements\n";

    for ( face_iterator elt = this->beginFace();
            elt != this->endFace(); ++elt )
    {
        face_type __face = *elt;
#if !defined( NDEBUG )
        Debug( 4015 ) << "face id: " << __face.id()
                      << " marker: " << __face.marker() << "\n";
#endif

        // renumber the nodes of the face
        for ( int i = 0; i < __face.nPoints(); ++i )
        {
            size_type __true_id =__face.point( i ).id();
            this->faces().modify( elt,
                                  lambda::bind( &face_type::setPoint,
                                                lambda::_1,
                                                lambda::constant( i ),
                                                boost::cref( this->point( node_map[__true_id] ) ) ) );
            face_type __face2 = *elt;
#if !defined( NDEBUG )
            Debug( 4015 ) << "id1= " << __face2.point( 0 ).id() << " id2= " << __face2.point( 1 ).id()<< "\n";
            Debug( 4015 ) << "point lid = " << i << " id = " << __true_id
                          << " nid = " << this->point( node_map[__true_id] ).id()
                          << " new point id = " << elt->point( i ).id() << "\n";
#endif

        }
    }
    renumber( node_map, mpl::int_<nDim>() );

    if ( Shape == SHAPE_TETRA && nOrder==1 )
    {
        localrenumber();

#if !defined(NDEBUG)
        checkLocalPermutation( mpl::bool_< ( Shape == SHAPE_TETRA ) >() );
#endif
    }

    // should we renumber also the faces and elements ?


}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::renumber( std::vector<size_type> const& node_map, mpl::int_<1> )
{
}
template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::renumber( std::vector<size_type> const& node_map, mpl::int_<2> )
{
}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::renumber( std::vector<size_type> const& node_map, mpl::int_<3> )
{

    for ( auto elt = this->beginEdge();
            elt != this->endEdge(); ++elt )
    {
        edge_type __edge = *elt;
#if !defined( NDEBUG )
        Debug( 4015 ) << "edge id: " << __edge.id()
                      << " marker: " << __edge.marker() << "\n";
#endif

        // renumber the nodes of the face
        for ( int i = 0; i < __edge.nPoints(); ++i )
        {
            size_type __true_id =__edge.point( i ).id();
            this->edges().modify( elt,
                                  lambda::bind( &edge_type::setPoint,
                                                lambda::_1,
                                                lambda::constant( i ),
                                                boost::cref( this->point( node_map[__true_id] ) ) ) );
            edge_type __edge2 = *elt;
#if !defined( NDEBUG )
            Debug( 4015 ) << "renumber edge: id1= " << __edge2.point( 0 ).id() << " id2= " << __edge2.point( 1 ).id()<< "\n";
            Debug( 4015 ) << "renumber edge: point lid = " << i << " id = " << __true_id
                          << " nid = " << this->point( node_map[__true_id] ).id()
                          << " new point id = " << elt->point( i ).id() << "\n";
#endif

        }
    }

}
template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::localrenumber()
{
    for ( element_const_iterator elt = this->beginElement();
            elt != this->endElement(); ++elt )
    {
        element_type __element = *elt;

        int16_type n_swap = 0;

        for ( uint16_type i=Shape::numVertices-1; i > 1; --i )
        {
            uint16_type min_index = 0;

            for ( uint16_type j=1; j < __element.nPoints()+i-( Shape::numVertices-1 ); ++j )
            {
                if ( __element.point( j ).id() < __element.point( min_index ).id() )
                    min_index = j;
            }

            if ( i != min_index )
            {
                __element.swapPoints( i, min_index );
                n_swap += int16_type( 1 );
            }
        }

        // Adds consistent permutation anti-clockwise for the two last points
        if ( n_swap == 1 )
            __element.swapPoints( 0,1 );

        // Modify the points ids
        if ( n_swap > 0 )
        {
            for ( uint16_type i=0; i < __element.nPoints(); ++i )
            {
                this->elements().modify( elt, typename super_elements::ElementUpdatePoint( i , this->point( __element.point( i ).id() ) ) );
            }
        }
    }
} /** void LocalRenumber **/

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::updateEntitiesCoDimensionOne()
{
    //updateEntitiesCoDimensionOne( mpl::bool_<nDim==nRealDim>() );
    updateEntitiesCoDimensionOne( mpl::bool_<true>() );
}
template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::updateEntitiesCoDimensionOne( mpl::bool_<false> )
{}
template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::updateEntitiesCoDimensionOne( mpl::bool_<true> )
{
    boost::timer ti;
    face_type face;
    //face.setWorldComm(this->worldComm());
    face.setProcessIdInPartition( this->worldComm().localRank() );

    std::map<std::set<int>, size_type > _faces;
    typename std::map<std::set<int>, size_type >::iterator _faceit;
    int next_face = 0;
    element_type ele;

    // First We check if we have already Faces stored
    if ( ! this->faces().empty() )
    {
        // dump all faces in the container, to maintain the correct numbering
        // if everything is correct the numbering in the bareface structure
        // will reflect the actual face numbering However, if I want to create
        // the internal faces I need to make sure that I am processing only the
        // boundary ones. So I resize the container!
        //if ( cf ) this->faces().resize( _numBFaces );


        face_iterator __it = this->beginFace();
        face_iterator __en = this->endFace();

        for ( ; __it!=__en; )
        {
            std::set<int> s;

            for ( int f = 0; f < face_type::numVertices; ++f )
            {
                s.insert( __it->point( f ).id() );
            }

            bool faceinserted = false;
            boost::tie( _faceit, faceinserted ) = _faces.insert( std::make_pair( s, next_face ) );

            if ( faceinserted )
                ++next_face;

#if !defined( NDEBUG )

            if ( faceinserted  )
                Debug( 4015 ) << "added face with id " << __it->id () << "\n";

            else
                Debug( 4015 ) << "not added face with id " << __it->id ()
                              << " was already face with id = " << _faceit->second << "\n";

            FEELPP_ASSERT( faceinserted )
            ( _faceit->second )
            ( __it->id() ).warn( "duplicated face" );
#endif

            if ( faceinserted == false )
            {
                // here we get the next face or \c end()
                size_type theid = __it->id();

                // find the other face
                face_iterator __other = this->faces().find( face_type( _faceit->second ) );
                FEELPP_ASSERT( __other->id() != theid )
                ( __other->id() )
                ( theid ).error( "faces should have different ids " );

                if ( __it->marker() != __other->marker() )
                {
                    this->faces().modify( __other, [__it]( face_type& f ) { f.setMarker2( __it->marker().value() ); } );
                }

                __it = this->eraseFace( __it );



            }

            else
            {
                // ensure that item handler ids are in sync with
                // faces ids
                face_type __f = *__it;
                __f.setId( _faceit->second );

                // set the id for this partition
                __f.setIdInPartition( this->worldComm().localRank(),__f.id() );

#if !defined( NDEBUG )
                Debug( 4015 ) << "set face id " << __f.id()
                              << " iterator id = " << __it->id()
                              << " check id = " << _faceit->second << "\n";
#endif
                this->faces().replace( __it, __f );
                ++__it;
            }
        }
    }

    Debug( 4015 ) << "[Mesh::updateFaces] adding faces : " << ti.elapsed() << "\n";
    ti.restart();

    Debug( 4015 ) << "[Mesh::updateFaces] numLocalFaces : " << this->numLocalFaces() << "\n";
    Debug( 4015 ) << "[Mesh::updateFaces] face_type::numVertices : " << face_type::numVertices << "\n";
    element_iterator iv,  en;
    boost::tie( iv, en ) = this->elementsRange();

    for ( ; iv != en; ++iv )
    {
        element_type const& __element = *iv;
        size_type __element_id = __element.id();

        //MakeBareEntity<element_type,nDim> baremaker( __element );
        for ( size_type j = 0; j < this->numLocalFaces(); j++ )
        {
            std::set<int> s;

            for ( int f = 0; f < face_type::numVertices; ++f )
            {
                uint16_type pt_localid = ( nDim==1 )?j:iv->fToP( j, f );
                s.insert( iv->point( pt_localid ).id() );
                Debug( 4015 ) << "add point local id " << f << " to face " << j  << " " << iv->fToP( j, f )
                              << " global id " << iv->point( pt_localid ).id() << "\n";
            }

            bool faceinserted = false;
            boost::tie( _faceit, faceinserted ) = _faces.insert( std::make_pair( s, next_face ) );

            if ( faceinserted )
                ++next_face;

#if !defined( NDEBUG )
            Debug( 4015 ) << "------------------------------------------------------------\n";
            Debug( 4015 ) << "Element id: " << iv->id() << " local face id: " << j << "\n";
#endif

            if ( faceinserted )
            {
#if !defined( NDEBUG )
                Debug( 4015 ) << "creating the face:" << _faceit->second << "\n";
#endif
                // set face id
                face.setId( _faceit->second );
                face.disconnect();
                face.addElement( __element_id );

                // set the process id from element
                face.setProcessId( __element.processId() );

                // set the id for this partition
                face.setIdInPartition( this->worldComm().localRank(),face.id() );

                // set the vertices of the face
                for ( size_type k = 0; k < face_type::numPoints; ++k )
                    face.setPoint( k, iv->point( ele.fToP( j, k ) ) );

                // set the connection with the element
                face.setConnection0( boost::make_tuple( boost::addressof( __element ), __element_id, j, __element.processId() ) );
                face.setOnBoundary( true );

                // adding the face
                bool inserted = false;
                face_iterator __fit;
                boost::tie( __fit, inserted )  = this->addFace( face );
                FEELPP_ASSERT( inserted && __fit != this->endFace() )
                ( _faceit->second )
                ( iv->id() )
                ( __fit->id() )
                ( face.id() ).error( "invalid face iterator" );
                this->elements().modify( iv, detail::UpdateFace<face_type>( boost::cref( *__fit ) ) );

#if !defined( NDEBUG )
                Debug( 4015 ) << "Adding [new] face info : \n";
                Debug( 4015 ) << "element id: " << __element_id << "\n";
                Debug( 4015 ) << "process id: " << __fit->processId() << "\n";
                Debug( 4015 ) << "id: " << __fit->id() << "\n";
                Debug( 4015 ) << "bdy: " << __fit->isOnBoundary() << "\n";
                Debug( 4015 ) << "marker: " << __fit->marker() << "\n";
                Debug( 4015 ) << "ad_first: " << __fit->ad_first() << "\n";
                Debug( 4015 ) << "pos_first: " << __fit->pos_first() << "\n";
                Debug( 4015 ) << "proc_first: " << __fit->proc_first() << "\n";
                Debug( 4015 ) << "ad_second: " << __fit->ad_second() << "\n";
                Debug( 4015 ) << "pos_second: " << __fit->pos_second() << "\n";
                Debug( 4015 ) << "proc_second: " << __fit->proc_second() << "\n";
#endif
            }

            else
            {
#if !defined( NDEBUG )
                Debug( 4015 ) << "found the face:" << _faceit->second << " in element " << __element_id << " and local face: " << j << "\n";
#endif

                // look in the face table for the face
                face_iterator __fit = this->faces().find( face_type( _faceit->second ) );
                FEELPP_ASSERT( __fit != this->endFace() )( _faceit->second ).error( "face is not in face container" );


                face_type face = *__fit;
                face.addElement( __element_id );

                // the three conditions below typically arise after reading a serialized mesh
                if ( __fit->isConnectedTo0() && __fit->connection0().template get<0>() == 0 && ( __element.id() == __fit->ad_first() ) )
                {
                    Debug( 4015 ) << "fixing connection 0 in face\n";
                    // reconnect the elements to the face
                    auto connect0 = __fit->connection0();
                    connect0.template get<0>() = boost::addressof( __element );
                    face.setConnection0( connect0 );

                    // need to reconstruct the neighbors
                    this->elements().modify( iv, detail::UpdateFace<face_type>( boost::cref( *__fit ) ) );

                }
                if ( __fit->isConnectedTo1() && __fit->connection1().template get<0>() == 0 && ( __element.id() == __fit->ad_second() ) )
                {
                    Debug( 4015 ) << "fixing connection 1 in face\n";

                    // reconnect the elements to the face
                    auto connect1 = __fit->connection1();
                    connect1.template get<0>() = boost::addressof( __element );
                    face.setConnection1( connect1 );

                    // need to reconstruct the neighbors
                    element_iterator elt1 = this->elementIterator( __fit->ad_first(), __fit->proc_first() );
                    this->elements().modify( elt1, detail::UpdateFace<face_type>( boost::cref( *__fit ) ) );
                    this->elements().modify( iv, detail::UpdateFace<face_type>( boost::cref( *__fit ) ) );

                }
                if ( __fit->isConnectedTo0() && __fit->isConnectedTo1() )
                {
                    FEELPP_ASSERT( face.isConnectedTo0() && face.isConnectedTo1() )
                        ( face.isConnectedTo0() )( face.isConnectedTo1() ).error ("inconsistent data structure" );
                    if ( face.processId()!=this->worldComm().localRank() )
                    {
                        if ( ( face.element0().processId()==this->worldComm().localRank() ) || ( face.element1().processId()==this->worldComm().localRank() ) )
                            face.setProcessId( this->worldComm().localRank() );
                    }
                }


                // the face could have been entered apriori given by
                // the mesh generator, so just set the connection0
                // properly .
                if ( !__fit->isConnectedTo0() )
                {
#if !defined( NDEBUG )
                    Debug( 4015 ) << "[updateFaces][boundary] element: " << __element_id
                                  << " face: " << j << " id: " << _faceit->second << "\n";
#endif

                    // set the connection with the element
                    face.setConnection0( boost::make_tuple( boost::addressof( __element ), __element_id, j, __element.processId() ) );
                    // set the process id from element
                    face.setProcessId( __element.processId() );

                    //this->faces().modify( __fit,
                    //detail::UpdateFaceConnection0<typename face_type::element_connectivity_type>( boost::make_tuple( boost::addressof( __element ), __element_id, j, __element.processId() ) ) );

                    this->faces().replace( __fit, face );

                    this->elements().modify( iv, detail::UpdateFace<face_type>( boost::cref( *__fit ) ) );
#if !defined( NDEBUG )
                    Debug( 4015 ) << "adding [!isConnectedTo0] face info : \n";
                    Debug( 4015 ) << "id: " << __fit->id() << "\n";
                    Debug( 4015 ) << "process id: " << __fit->processId() << "\n";
                    Debug( 4015 ) << "bdy: " << __fit->isOnBoundary() << "\n";
                    Debug( 4015 ) << "marker: " << __fit->marker() << "\n";
                    Debug( 4015 ) << "ad_first: " << __fit->ad_first() << "\n";
                    Debug( 4015 ) << "pos_first: " << __fit->pos_first() << "\n";
                    Debug( 4015 ) << "proc_first: " << __fit->proc_first() << "\n";
                    Debug( 4015 ) << "ad_second: " << __fit->ad_second() << "\n";
                    Debug( 4015 ) << "pos_second: " << __fit->pos_second() << "\n";
                    Debug( 4015 ) << "proc_second: " << __fit->proc_second() << "\n";
                    Debug( 4015 ) << "element process id: " << iv->processId() << "\n";
#endif
                }

                // we found an internal face
                else if ( !__fit->isConnectedTo1() )
                {


#if 0
                    this->faces().modify( __fit,
                                          detail::UpdateFaceConnection1<typename face_type::element_connectivity_type>( boost::make_tuple( boost::addressof( __element ), __element_id, j, __element.processId() ) ) );


                    FEELPP_ASSERT( __fit->isConnectedTo0() && __fit->isConnectedTo1() )
                    ( __fit->isConnectedTo0() )( __fit->isConnectedTo1() ).error( "invalid face connection" );
#endif
                    detail::UpdateFaceConnection1<typename face_type::element_connectivity_type> update1( boost::make_tuple( boost::addressof( __element ), __element_id, j, __element.processId() ) );
                    update1( face );
                    face.setOnBoundary( false );

                    // force processId equal to M_comm.rank() if face on interprocessfaces
                    if ( face.processId()!=this->worldComm().localRank() )
                    {
                        if ( ( face.element0().processId()==this->worldComm().localRank() ) || ( face.element1().processId()==this->worldComm().localRank() ) )
                            face.setProcessId( this->worldComm().localRank() );
                    }

                    this->faces().replace( __fit, face );
#if 0
                    // update neighbors for each element and replace in element container
                    element_iterator elt1 = this->elementIterator( __fit->ad_first(), __fit->proc_first() );
                    this->elements().modify( elt1, update_element_neighbor_type( __fit->pos_first(),
                                             __fit->ad_second() ) );
                    this->elements().modify( iv, update_element_neighbor_type( __fit->pos_second(),
                                             __fit->ad_first() ) );
#endif
                    // be careful: this step is crucial to set proper neighbor
                    // connectivity
                    element_iterator elt1 = this->elementIterator( __fit->ad_first(), __fit->proc_first() );
                    this->elements().modify( elt1, detail::UpdateFace<face_type>( boost::cref( *__fit ) ) );
                    this->elements().modify( iv, detail::UpdateFace<face_type>( boost::cref( *__fit ) ) );

#if !defined( NDEBUG )
                    Debug( 4015 ) << "adding face info : \n";
                    Debug( 4015 ) << "id: " << __fit->id() << "\n";
                    Debug( 4015 ) << "process id: " << __fit->processId() << "\n";
                    Debug( 4015 ) << "bdy: " << __fit->isOnBoundary() << "\n";
                    Debug( 4015 ) << "marker: " << __fit->marker() << "\n";
                    Debug( 4015 ) << "ad_first: " << __fit->ad_first() << "\n";
                    Debug( 4015 ) << "pos_first: " << __fit->pos_first() << "\n";
                    Debug( 4015 ) << "proc_first: " << __fit->proc_first() << "\n";
                    Debug( 4015 ) << "ad_second: " << __fit->ad_second() << "\n";
                    Debug( 4015 ) << "pos_second: " << __fit->pos_second() << "\n";
                    Debug( 4015 ) << "proc_second: " << __fit->proc_second() << "\n";
                    Debug( 4015 ) << "element1 process id: " << elt1->processId() << "\n";
                    Debug( 4015 ) << "element2 process id: " << iv->processId() << "\n";
#endif

                }

                FEELPP_ASSERT( (__fit->processId() == __fit->proc_first()) ||
                               (__fit->processId() == __fit->proc_second()) )
                    ( __fit->processId() )( __fit->proc_first() )( __fit->proc_second() ).error( "invalid process id" );
            }

            FEELPP_ASSERT( iv->facePtr( j ) )( j )( iv->id() ).error( "invalid element face error" );
        } // face loop
    } // element loop

#if 0
    face_iterator f_it = this->beginFace();
    face_iterator f_en = this->endFace();

    for ( ; f_it!=f_en; ++f_it )
    {
        // cleanup the face data structure :

        if ( !f_it->isConnectedTo0() )
        {
            // remove all faces that are not connected to any elements
            this->faces().erase( f_it );
        }

    }

#endif


    Debug( 4015 ) << "[Mesh::updateFaces] element/face connectivity : " << ti.elapsed() << "\n";
    ti.restart();
}
template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::updateOnBoundary( mpl::int_<1> )
{
    element_iterator iv,en;
    boost::tie( iv, en ) = this->elementsRange();
    for ( ; iv != en; ++iv )
    {
        bool isOnBoundary = false;

        for ( size_type j = 0; j < this->numLocalFaces(); j++ )
        {
            isOnBoundary |= iv->face( j ).isOnBoundary();
        }

        // an element on the boundary means that is shares a face
        // with the boundary
        this->elements().modify( iv, detail::OnBoundary( isOnBoundary ) );
    }
}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::updateOnBoundary( mpl::int_<2> )
{
    element_iterator iv,en;
    boost::tie( iv, en ) = this->elementsRange();
    for ( ; iv != en; ++iv )
    {
        bool isOnBoundary = false;

        for ( size_type j = 0; j < this->numLocalFaces(); j++ )
        {
            isOnBoundary |= iv->face( j ).isOnBoundary();
        }

        // an element on the boundary means that is shares a face
        // with the boundary
        this->elements().modify( iv, detail::OnBoundary( isOnBoundary ) );
    }
}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::updateOnBoundary( mpl::int_<3> )
{
    element_iterator iv,en;
    boost::tie( iv, en ) = this->elementsRange();
    for ( ; iv != en; ++iv )
    {
        bool isOnBoundary = false;

        for ( size_type j = 0; j < this->numLocalFaces(); j++ )
        {
            isOnBoundary |= iv->face( j ).isOnBoundary();
        }
        for ( size_type j = 0; j < this->numLocalEdges(); j++ )
        {
            isOnBoundary |= iv->edge( j ).isOnBoundary();
        }

        // an element on the boundary means that is shares a face
        // with the boundary
        this->elements().modify( iv, detail::OnBoundary( isOnBoundary ) );
    }
}
#if defined(FEELPP_ENABLE_MPI_MODE)
template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::updateEntitiesCoDimensionOneGhostCell()
{
    //std::cout << "[Mesh::updateEntitiesCoDimensionOneGhostCell] start on god rank "<< this->worldComm().godRank() << std::endl;
    std::vector<int> nbMsgToSend( this->worldComm().localSize() );
    std::fill( nbMsgToSend.begin(),nbMsgToSend.end(),0 );

    std::vector< std::map<int,int> > mapMsg( this->worldComm().localSize() );

    auto iv = this->beginGhostElement();
    auto en = this->endGhostElement();

    for ( ; iv != en; ++iv )
    {
        element_type const& __element = *iv;
        int IdProcessOfGhost = __element.processId();
        int idInPartition = __element.idInPartition( IdProcessOfGhost );
        // send
        this->worldComm().localComm().send( IdProcessOfGhost , nbMsgToSend[IdProcessOfGhost], idInPartition );
#if 0
        std::cout<< "I am the proc" << this->worldComm().localRank()<<" , I send to proc " << IdProcessOfGhost
                 <<" with tag "<< nbMsgToSend[IdProcessOfGhost]
                 << " idSend " << idInPartition
                 << " it_ghost->G() " << __element.G()
                 << std::endl;
#endif
        // save tag of request
        mapMsg[IdProcessOfGhost].insert( std::make_pair( nbMsgToSend[IdProcessOfGhost],__element.id() ) );
        // update nb send
        ++nbMsgToSend[IdProcessOfGhost];
    }

    //------------------------------------------------------------------------------------------------//
    // counter of msg received for each process
    std::vector<int> nbMsgToRecv;
    mpi::all_to_all( this->worldComm().localComm(),
                     nbMsgToSend,
                     nbMsgToRecv );

    //------------------------------------------------------------------------------------------------//
    // recv id asked and re-send set of face id
    for ( int proc=0; proc<this->worldComm().localSize(); ++proc )
    {
        for ( int cpt=0; cpt<nbMsgToRecv[proc]; ++cpt )
        {
            int idRecv;
            //recv
            this->worldComm().localComm().recv( proc, cpt, idRecv );
#if 0
            std::cout<< "I am the proc" << this->worldComm().localRank()<<" I receive to proc " << proc
                     <<" with tag "<< cpt
                     << " idRecv " << idRecv
                     << " it_ghost->G() " << this->element( idRecv ).G()
                     << std::endl;
#endif
            // get faces id
            std::vector<float/*int*/> idFaces( ( 1+nDim )*this->numLocalFaces() );

            for ( size_type j = 0; j < this->numLocalFaces(); j++ )
            {
                idFaces[( 1+nDim )*j]=this->element( idRecv ).face( j ).id();
                auto const& theface = this->element( idRecv ).face( j );
                auto const& theGj = theface.G();
#if 1
                //compute face barycenter
                typename Localization::matrix_node_type v( theGj.size1(), 1 );
                ublas::scalar_vector<T> avg( theGj.size2(), T( 1 ) );
                T n_val = int( theGj.size2() );

                for ( size_type i = 0; i < theGj.size1(); ++i )
                    v( i, 0 ) = ublas::inner_prod( ublas::row( theGj, i ), avg )/n_val;

                auto baryFace = ublas::column( v,0 );
#else // doesn't work
                /*auto GLASbaryFace =*/ //Feel::glas::average(jkj);
                auto baryFace = ublas::column( glas::average( theface.G() ),0 );
#endif

                // save facebarycenter by components
                for ( uint16_type comp = 0; comp<nDim; ++comp )
                {
                    //idFaces[(1+nDim)*j+comp+1]=this->element(idRecv).faceBarycenter(j)[comp];
                    idFaces[( 1+nDim )*j+comp+1]=baryFace[comp];
                }
            }

            // response
            this->worldComm().localComm().send( proc, cpt, idFaces );
        }
    }

    //------------------------------------------------------------------------------------------------//
    // get response to initial request and update Feel::Mesh::Faces data
    for ( int proc=0; proc<this->worldComm().localSize(); ++proc )
    {
        for ( int cpt=0; cpt<nbMsgToSend[proc]; ++cpt )
        {
            std::vector<float/*int*/> idFacesRecv( ( 1+nDim )*this->numLocalFaces() );
            //recv
            this->worldComm().localComm().recv( proc, cpt, idFacesRecv );
#if 0
            std::cout<< "I am the proc " << this->worldComm().localRank()<<" I receive to proc " << proc
                     <<" with tag "<< cpt
                     << " idFacesRecv " << idFacesRecv[0] << " " << idFacesRecv[1] << " "<< idFacesRecv[2]
                     << std::endl;
#endif
            //update data
            auto const& theelt = this->element( mapMsg[proc][cpt],proc );

            for ( size_type j = 0; j < this->numLocalFaces(); j++ )
            {
                //objective : find  face_it (hence jBis in theelt ) (permutations would be necessary)
                uint16_type jBis = invalid_uint16_type_value;
                bool hasFind=false;

                for ( uint16_type/*size_type*/ j2 = 0; j2 < this->numLocalFaces(); j2++ )
                {

                    auto const& thefacej2 = theelt.face( j2 );
                    auto const& theGj2 = thefacej2.G();
#if 1
                    //compute face barycenter
                    typename Localization::matrix_node_type v( theGj2.size1(), 1 );
                    ublas::scalar_vector<T> avg( theGj2.size2(), T( 1 ) );
                    T n_val = int( theGj2.size2() );

                    for ( size_type i = 0; i < theGj2.size1(); ++i )
                        v( i, 0 ) = ublas::inner_prod( ublas::row( theGj2, i ), avg )/n_val;

                    auto baryFace = ublas::column( v,0 );
#else // doesn't compile (I don't now why)
                    auto const& thefacej2 =theelt.face( j2 );
                    auto baryFace = ublas::column( glas::average( thefacej2.G() ),0 );
#endif

                    if ( nRealDim==1 )
                    {
                        if ( std::abs( /*theelt.faceBarycenter(j2)*/baryFace[0] - idFacesRecv[( 1+nDim )*j+1] ) < 1e-8 )
                        {
                            hasFind=true;
                            jBis=j2;
                        }
                    }

                    else if ( nRealDim==2 )
                    {
#if 0
                        std::cout << "rank " << this->worldComm().localRank() << " bary j2 "<< j2 << " " << theelt.barycenter() << " " << theelt.faceBarycenter( j2 ) << " "
                                  << ublas::column( glas::average( hola.G() ),0 ) << " "
                                  << hola/*theelt.face(j2)*/.barycenter()[0] << " " << idFacesRecv[( 1+nDim )*j+1] << " "
                                  << theelt.face( j2 ).barycenter()[1] << " " << idFacesRecv[( 1+nDim )*j+2] << std::endl;
#endif

                        if ( ( std::abs( /*theelt.faceBarycenter(j2)*/baryFace[0] - idFacesRecv[( 1+nDim )*j+1] ) < 1e-5 ) &&
                                ( std::abs( /*theelt.faceBarycenter(j2)*/baryFace[1] - idFacesRecv[( 1+nDim )*j+2] ) < 1e-5 ) )
                        {
                            hasFind=true;
                            jBis=j2;
                        }
                    }

                    else if ( nRealDim==3 )
                    {
                        if ( ( std::abs( /*theelt.faceBarycenter(j2)*/baryFace[0] - idFacesRecv[( 1+nDim )*j+1] ) < 1e-5 ) &&
                                ( std::abs( /*theelt.faceBarycenter(j2)*/baryFace[1] - idFacesRecv[( 1+nDim )*j+2] ) < 1e-5 ) &&
                                ( std::abs( /*theelt.faceBarycenter(j2)*/baryFace[2] - idFacesRecv[( 1+nDim )*j+3] ) < 1e-5 ) )
                        {
                            hasFind=true;
                            jBis=j2;
                        }
                    }
                }

                if ( !hasFind ) std::cout << "[mesh::updateEntitiesCoDimensionOneGhostCell] : PROBLEM NOT FIND" << std::endl;

                // get the good face
                //auto face_it = faceIterator(theelt.face(j).id());
                auto face_it = this->faceIterator( theelt.face( jBis ).id() );
                //update the face
                this->faces().modify( face_it, detail::update_id_in_partition_type( proc, idFacesRecv[( 1+nDim )*j] ) );

            } // for ( size_type j = 0; j < this->numLocalFaces(); j++ )
        } // for ( int cpt=0;cpt<nbMsgToSend[proc];++cpt)
    } // for (int proc=0; proc<M_comm.size();++proc)

    //------------------------------------------------------------------------------------------------//

    //std::cout << "[Mesh::updateEntitiesCoDimensionOneGhostCell] finish" << std::endl;

} // updateEntitiesCoDimensionOneGhostCell
#endif // defined(FEELPP_ENABLE_MPI_MODE)


template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::check() const
{
    if ( nDim != nRealDim )
        return;
#if !defined( NDEBUG )
    Debug( 4015 ) << "[Mesh::check] numLocalFaces = " << this->numLocalFaces() << "\n";
    element_iterator iv = this->beginElementWithProcessId( this->worldComm().localRank() );
    element_iterator en = this->endElementWithProcessId( this->worldComm().localRank() );

    //boost::tie( iv, en ) = this->elementsRange();
    for ( ; iv != en; ++iv )
    {
        element_type const& __element = *iv;

        for ( size_type j = 0; j < this->numLocalFaces(); j++ )
        {
            FEELPP_ASSERT( iv->facePtr( j ) )( j )( iv->id() ).error( "invalid element face check" );
            Debug( 4015 ) << "------------------------------------------------------------\n";
            Debug( 4015 ) << "Element : " << iv->id() << " face lid: " << j << " face gid:  " << iv->face( j ).id() << "\n";

        }

        size_type counter = 0;

        for ( uint16_type ms=0; ms < __element.nNeighbors(); ms++ )
        {
            if ( __element.neighbor( ms ).first != invalid_size_type_value )
                ++counter;

        }

        Debug( 4015 ) << "[Mesh::check] element " << __element.id() << " number of neighbors: " << counter << "\n";
        FEELPP_ASSERT( counter >= 1 )( __element.id() )( __element.nNeighbors() )( counter ).warn( "invalid neighboring data" );
#if 0

        for ( size_type j = 0; j < ( size_type )element_type::numEdges; ++j )
        {
            FEELPP_ASSERT( iv->edgePtr( j ) )( j )( iv->id() ).error( "invalid element edge check" );
            Debug( 4015 ) << "------------------------------------------------------------\n";
            Debug( 4015 ) << "Element : " << iv->id() << " edge lid: " << j << " edge gid:  " << iv->edge( j ).id() << "\n";

        }

#endif

    }

    // faces check
    typedef typename super::location_face_const_iterator location_face_const_iterator;
    location_face_const_iterator itf = this->beginFaceOnBoundary();
    location_face_const_iterator ite = this->endFaceOnBoundary();

    for ( ; itf != ite; ++ itf )
    {
        face_type const& __face = *itf;
        FEELPP_ASSERT( __face.isConnectedTo0() )
        ( __face.id() )( __face.G() )( __face.ad_first() )( __face.pos_first() )( __face.proc_first() ).warn( "invalid face" );

        if ( __face.isConnectedTo0() )
        {
            FEELPP_ASSERT( __face.element( 0 ).facePtr( __face.pos_first() ) )
            ( __face.ad_first() )( __face.pos_first() )( __face.proc_first() )( __face.element( 0 ).id() ).warn( "invalid face in element" );
        }

        FEELPP_ASSERT( !__face.isConnectedTo1() )
        ( __face.ad_first() )( __face.pos_first() )( __face.proc_first() ).warn( "invalid boundary face" );

    }

#endif
}
template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::findNeighboringProcessors()
{
    // Don't need to do anything if there is
    // only one processor.
    if ( this->worldComm().localSize() == 1 )
        return;

#ifdef FEELPP_HAS_MPI

    _M_neighboring_processors.clear();

    // Get the bounding sphere for the local processor
    Sphere bounding_sphere = processorBoundingSphere ( *this, this->worldComm().localRank() );

    // Just to be sure, increase its radius by 10%.  Sure would suck to
    // miss a neighboring processor!
    bounding_sphere.setRadius( bounding_sphere.radius()*1.1 );

    // Collect the bounding spheres from all processors, test for intersection
    {
        std::vector<float>
        send ( 4,                         0 ),
             recv ( 4*this->worldComm().localSize(), 0 );

        send[0] = bounding_sphere.center()( 0 );
        send[1] = bounding_sphere.center()( 1 );
        send[2] = bounding_sphere.center()( 2 );
        send[3] = bounding_sphere.radius();

        MPI_Allgather ( &send[0], send.size(), MPI_FLOAT,
                        &recv[0], send.size(), MPI_FLOAT,
                        this->worldComm().localComm() );


        for ( unsigned int proc=0; proc<this->worldComm().localSize(); proc++ )
        {
            const Point center ( recv[4*proc+0],
                                 recv[4*proc+1],
                                 recv[4*proc+2] );

            const Real radius = recv[4*proc+3];

            const Sphere proc_sphere ( center, radius );

            if ( bounding_sphere.intersects( proc_sphere ) )
                _M_neighboring_processors.push_back( proc );
        }

        // Print out the _neighboring_processors list
        Debug( 4015 ) << "Processor " << this->worldComm().localRank() << " intersects:\n";

        for ( unsigned int p=0; p< _M_neighboring_processors.size(); p++ )
            Debug( 4015 ) << " - proc " << _M_neighboring_processors[p] << "\n";
    }

#endif
}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::checkLocalPermutation( mpl::bool_<true> ) const
{
    bool mesh_well_oriented = true;
    std::vector<size_type> list_of_bad_elts;

    for ( element_const_iterator elt = this->beginElement();
            elt != this->endElement(); ++elt )
    {
        element_type const& __element = *elt;

        if ( !__element.isAnticlockwiseOriented() )
        {
            mesh_well_oriented = false;
            list_of_bad_elts.push_back( elt->id() );
        }
    }

    if ( mesh_well_oriented )
        Debug( 4015 ) << "Local numbering in the elements is OK . \n";

    else
    {
        Debug( 4015 ) << "Local numbering in the elements is not anticlockwise oriented. \n";
        std::for_each( list_of_bad_elts.begin(),
                       list_of_bad_elts.end(),
                       std::cout << lambda::constant( "bad element " ) << lambda::_1 << lambda::constant( "\n" ) );

    }
}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::checkAndFixPermutation(  )
{
    element_type __element = *this->beginElement();

    if ( ! ( __element.isATriangleShape() ||
             __element.isATetrahedraShape() ) )
        return;

    static const uint16_type otn_triangle[ 10 ] =
    {
        1, 0, 2, 3, 5, 4
    };
    static const uint16_type otn_tetra[ 10 ] =
    {
        1, 0, 2, 3, 4, 6, 5, 8, 7, 9
    };

    for ( element_iterator elt = this->beginElement();
            elt != this->endElement(); ++elt )
    {
        element_type __element = *elt;
        bool is_anticlockwise = __element.isAnticlockwiseOriented();
        // --verbose
#if 0
        FEELPP_ASSERT( is_anticlockwise == true )
        ( is_anticlockwise )
        ( __element.id() )
        ( __element.G() ).warn( "invalid element permutation, will fix it" );
#endif // 0

        // fix permutation
        if ( !is_anticlockwise )
        {
            if ( __element.isATriangleShape() )
                __element.exchangePoints( otn_triangle );

            else if ( __element.isATetrahedraShape() )
                __element.exchangePoints( otn_tetra );

            else
            {
                FEELPP_ASSERT( 0 )
                ( __element.id() )
                ( __element.G() ).error( "invalid element type" );
                throw std::logic_error( "invalid element type" );
            }

            is_anticlockwise = __element.isAnticlockwiseOriented();
            FEELPP_ASSERT( is_anticlockwise == true )
            ( is_anticlockwise )
            ( __element.id() )
            ( __element.G() ).error( "invalid element permutation" );

            this->elements().replace( elt, __element );
        }

    }
}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::send(int p, int tag)
{
    encode();
    Debug() << "sending markername\n";
    //this->comm().send( p, tag, M_markername.size() );
    Debug() << "sending markername size: "<< M_markername.size() << "\n";
    BOOST_FOREACH(auto m, M_markername )
    {
        Debug() << "sending key: "<< m.first << "\n";
        //this->comm().send( p, tag, m.first );
        Debug() << "sending value\n";
        //this->comm().send( p, tag, m.second );
    }
}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::recv(int p, int tag)
{
    Debug() << "receiving markername\n";
    //this->comm().recv( p, tag, M_markername );
    int s;
    //this->comm().recv( p, tag, s );
    Debug() << "receiving markername size: "<< s << "\n";
    for( int i = 0; i < s; ++i )
    {
        std::string k;
        Debug() << "receiving key\n";
        //this->comm().recv( p, tag, k );
        Debug() << "receiving key:"<< k << "\n";
        std::vector<int> v;
        Debug() << "receiving value\n";
        //this->comm().recv( p, tag, v );
        Debug() << "receiving value: "<< v[0] << ","<< v[1] <<"\n";
        //M_markername[k]=v;
    }


    //decode();

}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::encode()
{
    //std::cout<<"encode=   " << this->worldComm().localSize() << std::endl;

    M_enc_pts.clear();
    for( auto pt_it = this->beginPoint(), pt_en = this->endPoint(); pt_it != pt_en; ++pt_it )
    {
        std::vector<double> pts(3,0);
        pts[0] = pt_it->node()[0];
        if ( mesh_type::nRealDim >= 2 )
            pts[1] = pt_it->node()[1];
        if ( mesh_type::nRealDim >= 3 )
            pts[2] = pt_it->node()[2];
        M_enc_pts[pt_it->id()+1] = boost::make_tuple( pt_it->isOnBoundary(), pt_it->tags(), pts );

    }
    M_enc_elts.clear();


    auto allmarkedfaces = boundaryfaces( *this );
    auto face_it = allmarkedfaces.template get<1>();
    auto face_end = allmarkedfaces.template get<2>();

    int elem_number=1;

    GmshOrdering<element_type> ordering;

    GmshOrdering<face_type> ordering_face;
    // save the faces

    for ( ; face_it != face_end; ++face_it )
    {
        std::vector<int> faces;
        faces.push_back( ordering_face.type() );
        faces.push_back( 4 + face_it->numberOfPartitions() );
        faces.push_back( face_it->marker().value() );
        faces.push_back( face_it->marker2().value() );
        faces.push_back( face_it->numberOfPartitions() );
        faces.push_back( face_it->processId() );
        for ( size_type i=0 ; i<face_it->numberOfNeighborPartitions(); ++i )
            faces.push_back( -( face_it->neighborPartitionIds()[i] ) );
        for ( uint16_type p=0; p<face_type::numPoints; ++p )
            faces.push_back( face_it->point( ordering_face.fromGmshId( p ) ).id()+1 );

        M_enc_faces[face_it->id()] = faces;
    } // faces


    auto eltOnProccess = elements( *this );
    auto elt_it = eltOnProccess.template get<1>();
    auto elt_en = eltOnProccess.template get<2>();

    for ( ; elt_it != elt_en; ++elt_it )
    {
        std::vector<int> elts;
        elts.push_back( ordering.type() );
        elts.push_back( 4 + elt_it->numberOfPartitions() );
        elts.push_back( elt_it->marker().value() );
        elts.push_back( elt_it->marker2().value() );
        elts.push_back( elt_it->numberOfPartitions() );
        elts.push_back( elt_it->processId() );
        for ( size_type i=0 ; i<elt_it->numberOfNeighborPartitions(); ++i )
            elts.push_back( -( elt_it->neighborPartitionIds()[i] ) );
        for ( uint16_type p=0; p<element_type::numPoints; ++p )
            elts.push_back( elt_it->point( ordering.fromGmshId( p ) ).id()+1 );

        M_enc_elts[elt_it->id()] = elts;
    } // elements
    //std::cout<<"encode=   " << this->worldComm().localSize() << std::endl;
}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::decode()
{
#if 0
    std::vector<int> mapWorld(this->worldComm().size());
    for(int cpu=0; cpu < this->worldComm().size(); ++cpu)
        mapWorld[cpu] = cpu;
    WorldComm worldcomm(mapWorld);




    // std::cout<<"decode=   " << this->worldComm().size() << std::endl;
    //std::cout<<"decode=   " << worldcomm.subWorldComm().localRank() << std::endl;
    this->setWorldComm(worldcomm.subWorldComm());
#else
    Log() <<"decode=   " << this->worldComm().size() << "\n" ;
    Log() <<"decode=   " << this->worldComm().subWorldComm().localRank() << "\n";
#endif
    static const uint16_type npoints_per_face = ( face_type::numVertices*face_type::nbPtsPerVertex+
            face_type::numEdges*face_type::nbPtsPerEdge+
            face_type::numFaces*face_type::nbPtsPerFace );

    static const uint16_type npoints_per_element = element_type::numPoints;
    for( auto pt_it = M_enc_pts.begin(), pt_en = M_enc_pts.end();
         pt_it != pt_en; ++pt_it )
    {
        node_type __n( nRealDim );

        for ( uint16_type j = 0; j < nRealDim; ++j )
            __n[j] = pt_it->second.template get<2>()[j];

        point_type __pt( pt_it->first-1,__n,  pt_it->second.template get<0>() );
        __pt.setOnBoundary( pt_it->second.template get<0>() );
        __pt.setTags( pt_it->second.template get<1>() );

        this->addPoint( __pt );
    }
    for( auto face_it = M_enc_faces.begin(), face_en = M_enc_faces.end();
         face_it != face_en; ++face_it )
    {
        face_type pf;
        GmshOrdering<face_type> ordering;
        FEELPP_ASSERT( ordering.type() == face_it->second[0] )
            ( ordering.type() ) ( face_it->second[0] ).error( "invalid mesh face type" );
        std::vector<int> tags( face_it->second[1] );
        for(int i = 0; i < tags.size(); ++i ) tags[i] = face_it->second[2+i];
        pf.setTags(  tags  );
        pf.setId( this->numFaces() );
        pf.setProcessIdInPartition( this->worldComm().localRank() );
        pf.setProcessId( this->worldComm().localRank() );
        pf.setIdInPartition( this->worldComm().localRank(),pf.id() );

        const int shift = face_it->second[1]+1;
        for ( uint16_type jj = 0; jj < npoints_per_face; ++jj )
        {
            pf.setPoint( ordering.fromGmshId( jj ), this->point( face_it->second[shift+jj]-1 ) );
        }
        this->addFace( pf );
    }

    for( auto elt_it = M_enc_elts.begin(), elt_en = M_enc_elts.end();
         elt_it != elt_en; ++elt_it )
    {
        element_type pv;
        GmshOrdering<element_type> ordering;
        FEELPP_ASSERT( ordering.type() == elt_it->second[0] )
            ( ordering.type() ) ( elt_it->second[0] ).error( "invalid mesh element type" );
        std::vector<int> tags( elt_it->second[1] );
        for(int i = 0; i < tags.size(); ++i ) tags[i] = elt_it->second[2+i];
        pv.setTags(  tags  );
        pv.setProcessIdInPartition( this->worldComm().localRank() );
        pv.setProcessId( this->worldComm().localRank() );
        pv.setIdInPartition( this->worldComm().localRank(),pv.id() );

        const int shift = elt_it->second[1]+1;
        for ( uint16_type jj = 0; jj < npoints_per_element; ++jj )
        {
            pv.setPoint( ordering.fromGmshId( jj ), this->point( elt_it->second[shift+jj]-1 ) );
        }
        this->addElement( pv );
#if 0
        __idGmshToFeel=pv.id();
        auto theelt = mesh->elementIterator( pv.id(), pv.partitionId() );
        mesh->elements().modify( theelt, detail::update_id_in_partition_type( this->worldComm().localRank(), pv.id() ) );
#endif
    }
    Log() << "distance  elts: "<< std::distance( this->beginElement(), this->endElement() ) << "\n";
    Log() << "distance faces: "<< std::distance( this->beginFace(), this->endFace() ) << "\n";
    Log() << "distance marker faces: "<< std::distance( this->beginFaceWithMarker(), this->endFaceWithMarker() ) << "\n";
    Log() << "distance marker2 faces: "<< std::distance( this->beginFaceWithMarker2(), this->endFaceWithMarker2() ) << "\n";

    //this->components().set ( MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
    //this->updateForUse();
    //std::cout<<"decode=   " << this->worldComm().localSize() << std::endl;
}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::Inverse::distribute( bool extrapolation )
{
    typename self_type::element_iterator el_it;
    typename self_type::element_iterator el_en;
    boost::tie( boost::tuples::ignore, el_it, el_en ) = Feel::elements( *M_mesh );

    typedef typename self_type::element_type element_type;
    typedef typename gm_type::template Context<vm::JACOBIAN|vm::KB|vm::POINT, element_type> gmc_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
    BoundingBox<> bb;

    typename gm_type::reference_convex_type refelem;
    typename gm_type::precompute_ptrtype __geopc( new typename gm_type::precompute_type( M_mesh->gm(),
            refelem.points() ) );
    std::vector<bool> npt( this->nPoints() );

    M_dist.resize( this->nPoints() );
    M_ref_coords.resize( this->nPoints() );
    M_cvx_pts.resize( this->nPoints() );
    M_pts_cvx.clear();
    M_pts_cvx.resize( M_mesh->numElements() );

    KDTree::points_type boxpts;
    gmc_ptrtype __c( new gmc_type( M_mesh->gm(),
                                   *el_it,
                                   __geopc ) );
    Debug( 4015 ) << "[Mesh::Inverse] distribute mesh points ion kdtree\n";

    for ( ; el_it != el_en; ++el_it )
    {
        // get geometric transformation
        __c->update( *el_it );
        gic_type gic( M_mesh->gm(), *el_it );

        // create bounding box
        //bb.make( el_it->points() );
        bb.make( el_it->G() );

        for ( size_type k=0; k < bb.min.size(); ++k )
        {
            bb.min[k] -= 1e-10;
            bb.max[k] += 1e-10;
        }

#if !defined( NDEBUG )
        Debug( 4015 ) << "G = " << el_it->G() << " min = " << bb.min << ", max = " << bb.max << "\n";
#endif /* NDEBUG */


        // check if the points
        this->pointsInBox( boxpts, bb.min, bb.max );
#if !defined( NDEBUG )
        Debug( 4015 ) << "boxpts size = " << boxpts.size() << "\n";
#endif /*  */

        for ( size_type i = 0; i < boxpts.size(); ++i )
        {
            size_type index = boost::get<1>( boxpts[i] );

            if ( ( !npt[index] ) || M_dist[index] < 0 )
            {
                // check if we are in
                gic.setXReal( boost::get<0>( boxpts[i] ) );
                bool isin;
                value_type dmin;
                boost::tie( isin, dmin ) = refelem.isIn( gic.xRef() );
                bool tobeadded = extrapolation || isin;

#if !defined( NDEBUG )
                Debug( 4015 ) << "i = " << i << " index = " << index << " isin = " << ( isin >= -1e-10 )  << " xref = " << gic.xRef() << " xreal = " << boost::get<0>( boxpts[i] ) << " tobeadded= " << tobeadded << " dist=" << dmin<< "\n";
#endif

                if ( tobeadded && npt[index] )
                {
                    if ( dmin > M_dist[index] )
                        M_pts_cvx[M_cvx_pts[index]].erase( index );

                    else
                        tobeadded = false;
                }

                if  ( tobeadded )
                {
                    M_ref_coords[index] = gic.xRef();
                    M_dist[index] = dmin;
                    M_cvx_pts[index]=el_it->id();
                    M_pts_cvx[el_it->id()][index]=boost::get<2>( boxpts[i] );
                    npt[index]=true;
                }
            }
        }

    }

    Debug( 4015 ) << "[Mesh::Inverse] distribute mesh points in kdtree done\n";
}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::Localization::init()
{
    if ( !M_mesh ) return;

#if !defined( NDEBUG )
    FEELPP_ASSERT( IsInit == false )
    ( IsInit ).warn( "You have already initialized the tool of localization" );
#endif


    //clear data
    M_geoGlob_Elts.clear();
    M_kd_tree->clear();

    typename self_type::element_iterator el_it;
    typename self_type::element_iterator el_en;
    boost::tie( boost::tuples::ignore, el_it, el_en ) = Feel::elements( *M_mesh );

    for ( ; el_it != el_en; ++el_it )
    {
        for ( int i=0; i<el_it->nPoints(); ++i )
        {
            if ( boost::get<1>( M_geoGlob_Elts[el_it->point( i ).id()] ).size()==0 )
            {
                boost::get<0>( M_geoGlob_Elts[el_it->point( i ).id()] ) = el_it->point( i ).node();
                M_kd_tree->addPoint( el_it->point( i ).node(),el_it->point( i ).id() );
            }

            boost::get<1>( M_geoGlob_Elts[el_it->point( i ).id()] ).push_back( el_it->id() );
        }
    }

    IsInit=true;
    IsInitBoundaryFaces=false;

}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::Localization::initBoundaryFaces()
{
    if ( !M_mesh ) return;

#if !defined( NDEBUG )
    FEELPP_ASSERT( IsInitBoundaryFaces == false )
    ( IsInitBoundaryFaces ).warn( "You have already initialized the tool of localization" );
#endif


    //clear data
    M_geoGlob_Elts.clear();
    M_kd_tree->clear()
;
    typename self_type::location_face_iterator face_it;
    typename self_type::location_face_iterator face_en;
    boost::tie( boost::tuples::ignore, face_it, face_en ) = Feel::boundaryfaces( M_mesh );

    for ( ; face_it != face_en; ++face_it )
    {
        for ( int i=0; i<face_it->nPoints(); ++i )
        {
            if ( face_it->isConnectedTo0() )
                {
                    if ( boost::get<1>( M_geoGlob_Elts[face_it->point( i ).id()] ).size()==0 )
                        {
                            boost::get<0>( M_geoGlob_Elts[face_it->point( i ).id()] ) = face_it->point( i ).node();
                            M_kd_tree->addPoint( face_it->point( i ).node(),face_it->point( i ).id() );
                        }
                    boost::get<1>( M_geoGlob_Elts[face_it->point( i ).id()] ).push_back( face_it->element( 0 ).id() );
                }
        }
    }

    IsInitBoundaryFaces=true;
    IsInit=false;

}



template<typename Shape, typename T, int Tag>
boost::tuple<bool,typename Mesh<Shape, T, Tag>::node_type,double>
Mesh<Shape, T, Tag>::Localization::isIn( size_type _id, const node_type & _pt ) const
{
    bool isin=false;
    double dmin;
    node_type x_ref;

    //get element with the id
    auto const& elt = M_mesh->element( _id );

    if ( elt.isOnBoundary() )
        {
            // get inverse geometric transformation
            gmc_inverse_type gic( M_mesh->gm(), elt, this->mesh()->worldComm().subWorldCommSeq() );
            //apply the inverse geometric transformation for the point p
            gic.setXReal( _pt);
            x_ref=gic.xRef();
            // the point is in the reference element ?
            boost::tie( isin, dmin ) = M_refelem.isIn( gic.xRef() );
        }
    else
        {
            // get inverse geometric transformation
            gmc1_inverse_type gic( M_mesh->gm1(), elt, mpl::int_<1>(), this->mesh()->worldComm().subWorldCommSeq() );
            //apply the inverse geometric transformation for the point p
            gic.setXReal( _pt);
            x_ref=gic.xRef();
            // the point is in the reference element ?
            boost::tie( isin, dmin ) = M_refelem1.isIn( gic.xRef() );
        }

    return boost::make_tuple(isin,x_ref,dmin);
}

template<typename Shape, typename T, int Tag>
boost::tuple<uint16_type,std::vector<bool> >
Mesh<Shape, T, Tag>::Localization::isIn( std::vector<size_type> _ids, const node_type & _pt )
{
    typedef typename self_type::gm_type::reference_convex_type ref_convex_type;
    typedef typename self_type::gm1_type::reference_convex_type ref_convex1_type;

    uint16_type nbId = _ids.size();
    std::vector<bool> isin( _ids.size(),false );
    bool isin2=false;
    double dmin;
    node_type __x_ref;

    uint16_type nbIsIn=0;

    for ( uint16_type i = 0; i< nbId ; ++i )
    {
        //get element with the id
        auto const& elt = M_mesh->element( _ids[i] );

        if ( elt.isOnBoundary() )
            {
                // get inverse geometric transformation
                gmc_inverse_type gic( M_mesh->gm(), elt );
                //apply the inverse geometric transformation for the point p
                gic.setXReal( _pt);
                __x_ref=gic.xRef();
                // the point is in the reference element ?
                boost::tie( isin2, dmin ) = M_refelem.isIn( gic.xRef() );
                isin[i] = isin2;
            }
        else
            {
                // get inverse geometric transformation
                gmc1_inverse_type gic( M_mesh->gm1(), elt, mpl::int_<1>() );
                //apply the inverse geometric transformation for the point p
                gic.setXReal( _pt);
                __x_ref=gic.xRef();
                // the point is in the reference element ?
                boost::tie( isin2, dmin ) = M_refelem1.isIn( gic.xRef() );
                isin[i] = isin2;
            }
        if (isin[i]) ++nbIsIn;
    }

    return boost::make_tuple( nbIsIn,isin );
}


template<typename Shape, typename T, int Tag>
boost::tuple<bool, size_type, typename Mesh<Shape, T, Tag>::node_type>
Mesh<Shape, T, Tag>::Localization::searchElement( const node_type & p )
{

#if !defined( NDEBUG )
    FEELPP_ASSERT( IsInit == true )
    ( IsInit ).warn( "You don't have initialized the tool of localization" );
#endif
    bool isin=false;double dmin=0;
    node_type x_ref;
    size_type idEltFound = this->mesh()->beginElementWithId(this->mesh()->worldComm().localRank())->id();

    std::list< std::pair<size_type, uint> > ListTri;
    this->searchInKdTree(p,ListTri);

    //research the element which contains the point p
    auto itLT=ListTri.begin();
    auto itLT_end=ListTri.end();

#if !defined( NDEBUG )
    //if(std::distance(itLT,itLT_end)==0) std::cout<<"\nListTri vide\n";
    FEELPP_ASSERT( std::distance( itLT,itLT_end )>0 ).error( " problem in list localization : is empty" );
#endif

    while ( itLT != itLT_end && !isin  )
    {
        //get element with the id
        //elt = M_mesh->element( itLT->first );

        // search point in this elt
        boost::tie(isin,x_ref,dmin) = this->isIn(itLT->first,p);
        // if not inside, continue the research with an other element
        if (!isin) ++itLT;
        else idEltFound=itLT->first;
    }


    if (!isin)
        {
            if( this->doExtrapolation() )
                {
                    //std::cout << "WARNING EXTRAPOLATION for the point" << p << std::endl;
                    //std::cout << "W";
                    auto const& eltUsedForExtrapolation = this->mesh()->element(ListTri.begin()->first);
                    gmc_inverse_type gic( this->mesh()->gm(), eltUsedForExtrapolation, this->mesh()->worldComm().subWorldCommSeq() );
                    //apply the inverse geometric transformation for the point p
                    gic.setXReal( p);
                    boost::tie(isin,idEltFound,x_ref) = boost::make_tuple(true,eltUsedForExtrapolation.id(),gic.xRef());
                }
            else
                {
                    idEltFound = this->mesh()->beginElementWithId(this->mesh()->worldComm().localRank())->id();
                    isin = false;
                    //x_ref=?
                }
        }

    return boost::make_tuple( isin, idEltFound, x_ref);

}

template<typename Shape, typename T, int Tag>
boost::tuple<std::vector<bool>, size_type>
Mesh<Shape, T, Tag>::Localization::run_analysis( const matrix_node_type & m,
                                            const size_type & eltHypothetical )
{

#if !defined( NDEBUG )
    FEELPP_ASSERT( IsInit == true )
    ( IsInit ).warn( "You don't have initialized the tool of localization" );
#endif

    bool find_x;
    size_type cv_id=eltHypothetical;
    node_type x_ref;double dmin=0;
    std::vector<bool> hasFindPts(m.size2(),false);

    M_resultAnalysis.clear();

    bool doExtrapolationAtStart = this->doExtrapolation();
    auto nPtMaxNearNeighborAtStart = this->kdtree()->nPtMaxNearNeighbor();


    // first step : no extrapolation
    if ( doExtrapolationAtStart ) this->setExtrapolation( false );

    // first currentEltHypothetical
    auto currentEltHypothetical = eltHypothetical;//this->mesh()->beginElement()->id();//eltHypothetical;
    for ( size_type i=0; i< m.size2(); ++i )
        {
            bool testHypothetical_find = false;

            if ( eltHypothetical!=invalid_size_type_value )
                {
                    boost::tie( testHypothetical_find,x_ref,dmin ) = this->isIn( currentEltHypothetical,ublas::column( m, i ) );
                }
            if ( testHypothetical_find )
                {
                    cv_id = currentEltHypothetical;
                    M_resultAnalysis[cv_id].push_back( boost::make_tuple(i,x_ref) );
                    hasFindPts[i]=true;
                }
            else // search kdtree
                {
                    // kdtree call
                    boost::tie( find_x, cv_id, x_ref ) = this->searchElement(ublas::column( m, i ));
                    // traitement
                    if (find_x) // if find : OK
                        {
                            M_resultAnalysis[cv_id].push_back( boost::make_tuple(i,x_ref) );
                            currentEltHypothetical = cv_id;
                            hasFindPts[i]=true;
                        }
                    else// if (false) // try an other method (no efficient but maybe a solution)
                        {
                            // search in all element
                            this->kdtree()->nbNearNeighbor(2*nPtMaxNearNeighborAtStart /*15*/ /*this->mesh()->numElements()*/);

                            boost::tie( find_x, cv_id, x_ref ) = this->searchElement(ublas::column( m, i ));
                            //revert parameter
                            this->kdtree()->nbNearNeighbor(nPtMaxNearNeighborAtStart);
                            // if find : OK (but strange!)
                            if (find_x)
                                {
                                    M_resultAnalysis[cv_id].push_back( boost::make_tuple(i,x_ref) );
                                    currentEltHypothetical = cv_id;
                                    hasFindPts[i]=true;
                                }
                            else if (doExtrapolationAtStart)// && this->mesh()->worldComm().localSize()==1)
                                {
                                    this->setExtrapolation(true);

                                    boost::tie( find_x, cv_id, x_ref ) = this->searchElement(ublas::column( m, i ));
                                    // normaly is find
                                    if (find_x)
                                        {
                                            M_resultAnalysis[cv_id].push_back( boost::make_tuple(i,x_ref) );
                                            currentEltHypothetical = cv_id;
                                            hasFindPts[i]=true;
                                        }

                                    this->setExtrapolation(false);
                                }
                            else
                                {
                                    //std::cout<<"\n Il y a un GROS Probleme de Localization\n";
                                }
                        }
                } // search kdtree

        } // for (size_type i=0;i< m.size2();++i)

    //revert parameter
    this->setExtrapolation(doExtrapolationAtStart);

    return boost::make_tuple(hasFindPts,cv_id);

} //run_analysis



template<typename Shape, typename T, int Tag>
boost::tuple<bool, std::list<boost::tuple<size_type, typename Mesh<Shape, T, Tag>::node_type> > >
Mesh<Shape, T, Tag>::Localization::searchElements( const node_type & p )
{

#if !defined( NDEBUG )
    FEELPP_ASSERT( IsInit == true )
    ( IsInit ).warn( "You don't have initialized the tool of localization" );
#endif

    //this->kdtree()->nbNearNeighbor(this->mesh()->numElements());

    std::list< std::pair<size_type, uint> > ListTri;
    searchInKdTree( p,ListTri );

    typename self_type::element_type elt;
    typename self_type::gm_type::reference_convex_type refelem;
    typename self_type::gm1_type::reference_convex_type refelem1;

    bool isin=false;
    double dmin;
    node_type x_ref;

    //research the element which contains the point p
    auto itLT=ListTri.begin();
    auto itLT_end=ListTri.end();

#if !defined( NDEBUG )
    //if(std::distance(itLT,itLT_end)==0) std::cout<<"\nListTri vide\n";
    FEELPP_ASSERT( std::distance( itLT,itLT_end )>0 ).error( " problem in list localization : is empty" );
#endif

    std::list<boost::tuple<size_type,node_type> > newlistelts;
    newlistelts.clear();
    bool find = false;
    bool finishSearch=false;
    while ( itLT != itLT_end && !finishSearch /*&& !isin*/  )
    {
#if 0
        //get element with the id
        elt= M_mesh->element( itLT->first );

        if ( elt.isOnBoundary() )
        {
            // get inverse geometric transformation
            typename self_type::Inverse::gic_type gic( M_mesh->gm(), elt );

            //apply the inverse geometric transformation for the point p
            gic.setXReal( p );
            __x_ref=gic.xRef();

            // the point is in the reference element ?
            boost::tie( isin, dmin ) = refelem.isIn( gic.xRef() );
        }

        else
        {
            // get inverse geometric transformation
            typename self_type::Inverse::gic1_type gic( M_mesh->gm1(), elt,mpl::int_<1>() );

            //apply the inverse geometric transformation for the point p
            gic.setXReal( p );
            __x_ref=gic.xRef();

            // the point is in the reference element ?
            boost::tie( isin, dmin ) = refelem1.isIn( gic.xRef() );
            //std::cout << "gic.xRef()" << gic.xRef() << std::endl;
        }
#else
        boost::tie(isin,x_ref, dmin) = this->isIn(itLT->first,p);
#endif
        if ( isin )
        {
            newlistelts.push_back( boost::make_tuple( itLT->first,x_ref ) );
            find = true;
            // if not on boundary -> finish for this point
            if (dmin>1e-7) finishSearch=true;
        }

        //if (find) std::cout << elt.G() << std::endl;

        //if not inside, continue the research with an other element
        //if (!isin) ++itLT;
        ++itLT;
    }

    if ( !find ) std::cout << "\n WARNING EXTRAPOLATION IN SEARCHELEMENTS!!!"<<std::endl;

    if ( find )
        return boost::make_tuple( true,newlistelts );

    else if ( !find && !M_doExtrapolation )
        return boost::make_tuple( false,newlistelts );

    else
    {
        //std::cout << "\n WARNING EXTRAPOLATION \n";
        itLT=ListTri.begin();
        elt= M_mesh->element( itLT->first );
        typename self_type::Inverse::gic_type gic( M_mesh->gm(), elt );
        //apply the inverse geometric transformation for the point p
        //gic.setXReal(boost::get<0>(*ptsNN.begin()));
        gic.setXReal( p );
        x_ref=gic.xRef();
        //return boost::make_tuple( true, itLT->first, __x_ref);
        newlistelts.push_back( boost::make_tuple( itLT->first,x_ref ) );
        find = true;
        return boost::make_tuple( true,newlistelts );
    }
} // searchElements


template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::Localization::searchInKdTree( const node_type & p,
        std::list< std::pair<size_type, uint> > & ListTri )
{
    //search for nearest points
    M_kd_tree->search( p );

    //get the results of research
    typename KDTree::points_search_type ptsNN = M_kd_tree->pointsNearNeighbor();

    typename KDTree::points_search_const_iterator itNN = ptsNN.begin();
    typename KDTree::points_search_const_iterator itNN_end = ptsNN.end();

#if !defined( NDEBUG )
    FEELPP_ASSERT( std::distance( itNN,itNN_end )>0 ).error( "none Near Neighbor Points are find" );
#endif

    //iterator on a l(ist index element
    typename std::list<size_type>::iterator itL;
    typename std::list<size_type>::iterator itL_end;

    //ListTri will contain the indices of elements (size_type)
    //and the number of occurence(uint)
    //std::list< std::pair<size_type, uint> > ListTri;
    std::list< std::pair<size_type, uint> >::iterator itLT;
    std::list< std::pair<size_type, uint> >::iterator itLT_end;

    //create of ListTri : sort largest to smallest occurrences
    //In case of equality : if the point is closer than another then it will be before
    //                      if it is in the same point then  the lowest index will be before
    for ( ; itNN != itNN_end; ++itNN )
    {
        itL= boost::get<1>( M_geoGlob_Elts[boost::get<3>( *itNN )] ).begin();
        itL_end= boost::get<1>( M_geoGlob_Elts[boost::get<3>( *itNN )] ).end();

        for ( ; itL != itL_end; ++itL )
        {
            itLT=ListTri.begin();
            itLT_end=ListTri.end();
            bool find=false;

            while ( itLT != itLT_end && !find )
            {
                if ( itLT->first == *itL ) find=true;

                else ++itLT;
            }

            if ( find )
            {
                uint nb=itLT->second+1;
                size_type numEl=itLT->first;
                ListTri.remove( *itLT );
                itLT=ListTri.begin();
                itLT_end=ListTri.end();
                bool find=false;

                while ( itLT != itLT_end && !find )
                {
                    if ( itLT->second < nb ) find=true;

                    else ++itLT;
                }

                ListTri.insert( itLT,std::make_pair( numEl,nb ) );
            }

            else ListTri.push_back( std::make_pair( *itL,1 ) );
        }
    }

}


template<typename Shape, typename T, int Tag>
boost::tuple<bool, size_type, typename Mesh<Shape, T, Tag>::node_type>
Mesh<Shape, T, Tag>::Localization::searchElement( const node_type & p,
        const matrix_node_type & setPoints,
        mpl::int_<1> /**/ )
{
    typename self_type::element_type elt;
    typename self_type::gm_type::reference_convex_type refelem;
    typename self_type::gm1_type::reference_convex_type refelem1;

    bool isin=false,isin2=false;
    double dmin;
    node_type x_ref;
    size_type idEltFound = this->mesh()->beginElementWithId(this->mesh()->worldComm().localRank())->id();

    std::list< std::pair<size_type, uint> > ListTri;
    searchInKdTree( p,ListTri );

    auto itLT=ListTri.begin();
    auto itLT_end=ListTri.end();

#if !defined( NDEBUG )
    //if(std::distance(itLT,itLT_end)==0) std::cout<<"\nListTri vide\n";
    FEELPP_ASSERT( std::distance( itLT,itLT_end )>0 ).error( " problem in list localization : is empty" );
#endif


    //research the element which contains the point p
    while ( itLT != itLT_end && !isin  )
    {
        //get element with the id
        elt= M_mesh->element( itLT->first );

        auto eltG = elt.G();
        std::vector<bool> find( setPoints.size2() );
        std::fill( find.begin(),find.end(),false );

        for ( size_type i=0; i< setPoints.size2(); ++i )
            {
                auto thePt = ublas::column( setPoints,i );
                find[i]=false;

                for ( size_type j=0; j<eltG.size2(); ++j )
                    {
                        auto ptjeltG = ublas::column( eltG,j );

                        if ( ptjeltG.size()==1 )
                            {
                                if ( std::abs( thePt( 0 )-ptjeltG( 0 ) )<1e-5 )
                                    find[i]=true;
                            }

                        else if ( ptjeltG.size()==2 )
                            {
                                if ( std::abs( thePt( 0 )-ptjeltG( 0 ) )<1e-5 &&
                                     std::abs( thePt( 1 )-ptjeltG( 1 ) )<1e-5 )
                                    find[i]=true;
                            }

                        else if ( ptjeltG.size()==3 )
                            {
                                if ( std::abs( thePt( 0 )-ptjeltG( 0 ) )<1e-5 &&
                                     std::abs( thePt( 1 )-ptjeltG( 1 ) )<1e-5 &&
                                     std::abs( thePt( 2 )-ptjeltG( 2 ) )<1e-5 )
                                    find[i]=true;
                            }
                    }
            }
        // check if all points are found or not
        bool isOK=true;

        for ( size_type i=0; i< setPoints.size2(); ++i )
            {
                isOK &= find[i];
            }

        if ( !isOK ) isin=false;
        else isin=true;

        if ( isin ) // just a check
            {
                boost::tie(isin2,x_ref,dmin) = this->isIn(itLT->first,p);
                if ( isin!=isin2) std::cout << "Bug Mesh::Localization::searchElement<true>" << std::endl;
            }

        //if not inside, continue the research with an other element
        if ( !isin ) ++itLT;
        else idEltFound=itLT->first;

    } //while ( itLT != itLT_end && !isin  )

    if (!isin)
        {
            if( this->doExtrapolation() )
                {
                    std::cout << "WARNING EXTRAPOLATION for the point" << p << std::endl;
                    //std::cout << "W";
                    auto const& eltUsedForExtrapolation = this->mesh()->element(ListTri.begin()->first);
                    gmc_inverse_type gic( this->mesh()->gm(), eltUsedForExtrapolation, this->mesh()->worldComm().subWorldCommSeq() );
                    //apply the inverse geometric transformation for the point p
                    gic.setXReal( p);
                    boost::tie(isin,idEltFound,x_ref) = boost::make_tuple(true,eltUsedForExtrapolation.id(),gic.xRef());
                }
            else
                {
                    idEltFound = this->mesh()->beginElementWithId(this->mesh()->worldComm().localRank())->id();
                    isin = false;
                    //x_ref=?
                }
        }

    return boost::make_tuple( isin, idEltFound, x_ref);

} //searchElement


template<typename Shape, typename T, int Tag>
boost::tuple<std::vector<bool>, size_type>
Mesh<Shape, T, Tag>::Localization::run_analysis( const matrix_node_type & m,
        const size_type & eltHypothetical,
        const matrix_node_type & setPoints,
        mpl::int_<1> /**/ )
{
#if !defined( NDEBUG )
    FEELPP_ASSERT( IsInit == true )
    ( IsInit ).warn( "You don't have initialized the tool of localization" );
#endif

    bool find_x=false;
    size_type cv_id=eltHypothetical;
    node_type x_ref;
    std::vector<bool> hasFindPts(setPoints.size2(),false);

    M_resultAnalysis.clear();

    for ( size_type i=0; i< m.size2(); ++i )
    {
        boost::tie( find_x, cv_id, x_ref ) = this->searchElement( ublas::column( m, i ),setPoints,mpl::int_<1>() );

        if ( find_x )
        {
            M_resultAnalysis[cv_id].push_back( boost::make_tuple( i,x_ref ) );
            hasFindPts[i]=true;
        }
        //else std::cout<<"\nNew Probleme Localization\n" << std::endl;

    }

    return boost::make_tuple(hasFindPts,cv_id);

} // run_analysis


template<typename Shape, typename T, int Tag>
typename Mesh<Shape, T, Tag>::node_type
Mesh<Shape, T, Tag>::Localization::barycenter() const
{
    return this->barycenter(mpl::int_<nRealDim>());
}

template<typename Shape, typename T, int Tag>
typename Mesh<Shape, T, Tag>::node_type
Mesh<Shape, T, Tag>::Localization::barycenter(mpl::int_<1> /**/) const
{
    node_type res(1);
    res(0)=0;
    for (size_type i = 0 ; i<this->kdtree()->nPoints() ; ++i)
        {
            auto const& pt = this->kdtree()->points()[i].template get<0>();
            res(0)+=pt(0);
        }
    res(0)/=this->kdtree()->nPoints();
    return res;
}
template<typename Shape, typename T, int Tag>
typename Mesh<Shape, T, Tag>::node_type
Mesh<Shape, T, Tag>::Localization::barycenter(mpl::int_<2> /**/) const
{
    node_type res(2);
    res(0)=0;res(1)=0;
    for (size_type i = 0 ; i<this->kdtree()->nPoints() ; ++i)
        {
            auto const& pt = this->kdtree()->points()[i].template get<0>();
            res(0)+=pt(0);res(1)+=pt(1);
        }
    res(0)/=this->kdtree()->nPoints();res(1)/=this->kdtree()->nPoints();
    return res;
}
template<typename Shape, typename T, int Tag>
typename Mesh<Shape, T, Tag>::node_type
Mesh<Shape, T, Tag>::Localization::barycenter(mpl::int_<3> /**/) const
{
    node_type res(3);
    res(0)=0;res(1)=0;res(2)=0;
    for (size_type i = 0 ; i<this->kdtree()->nPoints() ; ++i)
        {
            auto const& pt = this->kdtree()->points()[i].template get<0>();
            res(0)+=pt(0);res(1)+=pt(1);res(2)+=pt(2);
        }
    res(0)/=this->kdtree()->nPoints();res(1)/=this->kdtree()->nPoints();res(2)/=this->kdtree()->nPoints();
    return res;
}





#if 0
#if defined( FEELPP_INSTANTIATION_MODE )

# define DIMS1 BOOST_PP_TUPLE_TO_LIST(1,(1))
# define RDIMS1 BOOST_PP_TUPLE_TO_LIST(3,(1,2,3))
# define ORDERS1 BOOST_PP_TUPLE_TO_LIST(5,(1,2,3,4,5))

# define DIMS2 BOOST_PP_TUPLE_TO_LIST(1,(2))
# define RDIMS2 BOOST_PP_TUPLE_TO_LIST(2,(2,3))
# define ORDERS2 BOOST_PP_TUPLE_TO_LIST(5,(1,2,3,4,5))

# define DIMS3 BOOST_PP_TUPLE_TO_LIST(1,(3))
# define RDIMS3 BOOST_PP_TUPLE_TO_LIST(1,(3))
# define ORDERS3 BOOST_PP_TUPLE_TO_LIST(4,(1,2,3,4))

# define FACTORY_SIMPLEX(LDIM,LORDER,RDIM) template class Mesh<Simplex<LDIM, LORDER, RDIM> >;
# define FACTORY_HYPERCUBE(LDIM,LORDER,RDIM) template class Mesh<Hypercube<LDIM, LORDER, RDIM> >;

# define FACTORY_SIMPLEX_OP(_, GDO) FACTORY_SIMPLEX GDO
# define FACTORY_HYPERCUBE_OP(_, GDO) FACTORY_HYPERCUBE GDO


# define FACTORY_SIMPLEX_E(LDIM,LORDER,RDIM) extern template class Mesh<Simplex<LDIM, LORDER, RDIM> >;
# define FACTORY_HYPERCUBE_E(LDIM,LORDER,RDIM) extern template class Mesh<Hypercube<LDIM, LORDER, RDIM> >;

# define FACTORY_SIMPLEX_OP_E(_, GDO) FACTORY_HYPERCUBE GDO
# define FACTORY_HYPERCUBE_OP_E(_, GDO) FACTORY_HYPERCUBE_E GDO

#if !defined( FEELPP_MESH_IMPL_NOEXTERN )

BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_SIMPLEX_OP_E, 3, ( DIMS1, BOOST_PP_LIST_FIRST_N( FEELPP_MESH_MAX_ORDER, ORDERS1 ), RDIMS1 ) )
BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_HYPERCUBE_OP_E, 3, ( DIMS1, BOOST_PP_LIST_FIRST_N( FEELPP_MESH_MAX_ORDER, ORDERS1 ), RDIMS1 ) )

BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_SIMPLEX_OP_E, 3, ( DIMS2, BOOST_PP_LIST_FIRST_N( FEELPP_MESH_MAX_ORDER, ORDERS2 ), RDIMS2 ) )
BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_HYPERCUBE_OP_E, 3, ( DIMS2, BOOST_PP_LIST_FIRST_N( FEELPP_MESH_MAX_ORDER, ORDERS2 ), RDIMS2 ) )

BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_SIMPLEX_OP_E, 3, ( DIMS3, BOOST_PP_LIST_FIRST_N( FEELPP_MESH_MAX_ORDER, ORDERS3 ), RDIMS3 ) )
BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_HYPERCUBE_OP_E, 3, ( DIMS3, BOOST_PP_LIST_FIRST_N( FEELPP_MESH_MAX_ORDER, ORDERS3 ), RDIMS3 ) )

#endif // FEELPP_MESH_IMPL_NOEXTERN
#endif // FEELPP_INSTANTIATION_MODE
#endif
} // namespace Feel





#endif // __MESHIMPL_HPP
