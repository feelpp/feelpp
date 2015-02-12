/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
    M_numGlobalElements( 0 ),
    M_gm( new gm_type ),
    M_gm1( new gm1_type ),
    M_meas( 0 ),
    M_measbdy( 0 ),
    M_substructuring( false ),
    //M_part(),
    M_tool_localization( new Localization() )
{
    VLOG(2) << "[Mesh] constructor called\n";
}
template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::partition ( const uint16_type n_parts )
{
    //M_part->partition( *this, n_parts );
}

template<typename Shape, typename T, int Tag>
flag_type
Mesh<Shape, T, Tag>::markerId ( boost::any const& __marker )
{
    flag_type theflag = -1;
    if ( boost::any_cast<flag_type>( &__marker ) )
    {
        theflag = boost::any_cast<flag_type>( __marker);
    }
    else if ( boost::any_cast<int>( &__marker ) )
    {
        theflag = boost::any_cast<int>( __marker);
    }
    else if ( boost::any_cast<size_type>( &__marker ) )
    {
        theflag = boost::any_cast<size_type>( __marker);
    }
    else if ( boost::any_cast<uint16_type>( &__marker ) )
    {
        theflag = boost::any_cast<uint16_type>( __marker);
    }
    else if ( boost::any_cast<std::string>( &__marker ) )
    {
        theflag = this->markerName( boost::any_cast<std::string>( __marker) );
    }
    else
        CHECK( theflag != -1 ) << "invalid flag type\n";
    return theflag;
}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::updateForUse()
{
    VLOG(2) << "component     MESH_RENUMBER: " <<  this->components().test( MESH_RENUMBER ) << "\n";
    VLOG(2) << "component MESH_UPDATE_EDGES: " <<  this->components().test( MESH_UPDATE_EDGES ) << "\n";
    VLOG(2) << "component MESH_UPDATE_FACES: " <<  this->components().test( MESH_UPDATE_FACES ) << "\n";
    VLOG(2) << "component    MESH_PARTITION: " <<  this->components().test( MESH_PARTITION ) << "\n";

    if ( this->numElements() == 0 )
    {
        VLOG(2) << "No elements in Mesh? (with process rank " << this->worldComm().rank() <<")\n";
        if ( this->worldComm().localSize()==1 ) return;
    }

    boost::timer ti;

    VLOG(2) << "is already updated? : " << this->isUpdatedForUse() << "\n";
    if ( !this->isUpdatedForUse() )
    {
        if ( this->components().test( MESH_RENUMBER ) )
        {

            //this->renumber();
            VLOG(2) << "[Mesh::updateForUse] renumber : " << ti.elapsed() << "\n";
        }

        //
        // compute the Adjacency graph
        //
        ti.restart();
        VLOG(2) << "Compute adjacency graph\n";
        element_iterator iv,  en;
        boost::unordered_map<size_type,boost::tuple<size_type, uint16_type, size_type> > f2e;

        boost::unordered_map<std::set<int>, size_type > _faces;
        typename boost::unordered_map<std::set<int>, size_type >::iterator _faceit;
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

                DVLOG(2) << "------------------------------------------------------------\n";
                DVLOG(2) << "Element id: " << iv->id() << " local face id: " << j << " process id:" << iv->processId() << "\n";

                if ( faceinserted )
                {
                    DVLOG(2) << " new face " << _faceit->second << " is now in store with elt "
                             << f2e[_faceit->second].template get<0>()  << " and local face id "
                             <<  f2e[_faceit->second].template get<1>()  << "\n";


                    f2e[_faceit->second].template get<0>() = iv->id();
                    f2e[_faceit->second].template get<1>() = j;
                    M_e2f[std::make_pair(iv->id(),j)]=boost::make_tuple( _faceit->second, invalid_size_type_value );
                }

                else // already stored
                {
                    DVLOG(2) << "old face " << _faceit->second << " was already in store with elt " << f2e[_faceit->second].template get<0>() << " and local face id " <<  f2e[_faceit->second].template get<1>() << "\n";

                    f2e[_faceit->second].template get<2>() = iv->id();
                    M_e2f[std::make_pair(iv->id(),j)]=boost::make_tuple( _faceit->second, f2e[_faceit->second].template get<0>() );
                    M_e2f[std::make_pair(f2e[_faceit->second].template get<0>(), f2e[_faceit->second].template get<1>())] =
                        boost::make_tuple( _faceit->second, iv->id() );
                }


            } // local face
        } // element loop
        f2e.clear();
        VLOG(1) << "Compute adjacency graph done in " << ti.elapsed() << "\n";
#if 0

        // partition mesh
        if ( this->components().test( MESH_PARTITION ) || ( this->worldComm().localSize() > 1 ) )
        {
            boost::timer ti1;
            this->partition();
            VLOG(2) << "[Mesh::updateForUse] partition time : " << ti1.elapsed() << "\n";
        }

#endif

        if ( this->components().test( MESH_UPDATE_FACES ) )
        {
            ti.restart();
            // update connectivities of entities of co dimension 1
            this->updateEntitiesCoDimensionOne();
            // update permutation of entities of co-dimension 1
            this->updateEntitiesCoDimensionOnePermutation();

            VLOG(1) << "[Mesh::updateForUse] update entities of codimension 1 : " << ti.elapsed() << "\n";

        }

        if ( this->components().test( MESH_UPDATE_EDGES ) )
        {
            ti.restart();
            // update connectivities of entities of co dimension 2
            // (edges in 3D)
            this->updateEntitiesCoDimensionTwo();
            VLOG(1) << "[Mesh::updateForUse] update edges : " << ti.elapsed() << "\n";
        }

        if ( this->worldComm().localSize()>1 )
        {
            if ( this->components().test( MESH_UPDATE_FACES ) )
            {
                // warning : this function change the isOnBoundary for ghost faces
                // so need call before updateOnBoundary()
                if ( false )
                    this->updateEntitiesCoDimensionOneGhostCellByUsingBlockingComm();
                else
                    this->updateEntitiesCoDimensionOneGhostCellByUsingNonBlockingComm();

                auto ipfRange = this->interProcessFaces();
                for ( auto itf = ipfRange.first, enf = ipfRange.second ; itf!=enf ; ++itf )
                    this->addFaceNeighborSubdomain( itf->partition2() );
            }

            auto iv = this->beginGhostElement();
            auto const en = this->endGhostElement();
            for ( ; iv != en; ++iv )
                this->addNeighborSubdomain( iv->processId() );
        }

        if ( this->components().test( MESH_UPDATE_FACES ) ||
             this->components().test( MESH_UPDATE_EDGES )
             )
        {
            ti.restart();
            updateOnBoundary();
            VLOG(1) << "[Mesh::updateForUse] update on boundary : " << ti.elapsed() << "\n";
        }



        if ( this->components().test( MESH_ADD_ELEMENTS_INFO ) )
        {
            ti.restart();
            boost::tie( iv, en ) = this->elementsRange();
            for ( ; iv != en; ++iv )
            {
                this->elements().modify( iv,
                                         []( element_type& e )
                                         {
                                             for ( int i = 0; i < e.numPoints; ++i )
                                                 e.point( i ).addElement( e.id() );
                                         });
            }
            VLOG(1) << "[Mesh::updateForUse] update add element info : " << ti.elapsed() << "\n"; 
        }

        this->updateNumGlobalElements();

        if ( this->components().test( MESH_PROPAGATE_MARKERS ) )
        {
            ti.restart();
            propagateMarkers(mpl::int_<nDim>() );
            VLOG(1) << "[Mesh::updateForUse] update propagate markers : " << ti.elapsed() << "\n"; 
        }

        {
            ti.restart();
            for ( auto itf = this->beginFace(), ite = this->endFace(); itf != ite; ++ itf )
            {
                this->faces().modify( itf,[this]( face_type& f ) { f.setMesh( this ); } );
            }
            VLOG(1) << "[Mesh::updateForUse] update face mesh : " << ti.elapsed() << "\n"; 
        }


#if 1
        {
            ti.restart();
            element_iterator iv,  en;
            boost::tie( iv, en ) = this->elementsRange();
            auto pc = M_gm->preCompute( M_gm, M_gm->referenceConvex().vertices() );
            auto pcf =  M_gm->preComputeOnFaces( M_gm, M_gm->referenceConvex().barycenterFaces() );

            for ( ; iv != en; ++iv )
            {
                this->elements().modify( iv,
                                         [=,&pc,&pcf]( element_type& e )
                                         {
                                             e.setMeshAndGm( this, M_gm, M_gm1 );
                                             e.updateWithPc(pc, boost::ref( pcf) );
                                         } );
#if 0
                lambda::bind( &element_type::setMeshAndGm,
                              lambda::_1,
                              this, M_gm, M_gm1 );

                this->elements().modify( iv,
                                         lambda::bind( &element_type::updateWithPc,
                                                       lambda::_1, pc, boost::ref( pcf ) ) );
#endif
            }
        }

        // update hAverage, hMin, hMax, measure of the mesh and measure of the boundary mesh
        this->updateMeasures();

        // attach mesh in the localisation tool
        M_tool_localization->setMesh( this->shared_from_this(),false );

    } // isUpdatedForUse
#endif

    {
        ti.restart();
        // check mesh connectivity
        this->check();
        VLOG(1) << "[Mesh::updateForUse] check : " << ti.elapsed() << "\n";
    }

    if ( M_is_gm_cached == false )
    {
        ti.restart();
        M_gm->initCache( this );
        M_gm1->initCache( this );
        M_is_gm_cached = true;
        VLOG(1) << "[Mesh::updateForUse] update geomap : " << ti.elapsed() << "\n"; 
    }

    this->setUpdatedForUse( true );

    if (Environment::isMasterRank() && FLAGS_v >= 1)
    {
        std::cout << "[Mesh::updateForUse] total time : " << ti.elapsed() << "\n";
        auto mem  = Environment::logMemoryUsage("memory usage after update for use");
        std::cout << "[Mesh::updateForUse] resident memory : " << mem.memory_usage/1.e9 << "\n";
    }
}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::updateMeasures()
{
    boost::timer ti;

    M_local_meas = 0;
    M_local_measbdy = 0;
    element_iterator iv,  en;
    boost::tie( iv, en ) = this->elementsRange();
    for ( ; iv != en; ++iv )
    {
        // only compute meas for active element (no ghost)
        if ( !iv->isGhostCell() )
            M_local_meas += iv->measure();

        auto _faces = iv->faces();

        if ( nDim == 1 )
            M_local_measbdy = 0;
        else
            for ( ; _faces.first != _faces.second; ++_faces.first )
                if ( ( *_faces.first ) && ( *_faces.first )->isOnBoundary() )
                    M_local_measbdy += ( *_faces.first )->measure();
    }
#if BOOST_VERSION >= 105500
    std::vector<value_type> gmeas{ M_local_meas, M_local_measbdy };
    mpi::all_reduce(this->worldComm(), mpi::inplace(gmeas.data()), 2, std::plus<value_type>());
#else
    std::vector<value_type> lmeas{ M_local_meas, M_local_measbdy };
    std::vector<value_type> gmeas( 2, 0.0 );
    mpi::all_reduce(this->worldComm(), lmeas.data(), 2, gmeas.data(), std::plus<value_type>());
#endif
    M_meas = gmeas[0];
    M_measbdy = gmeas[1];

    // now that all elements have been updated, build inter element
    // data such as the measure of point element neighbors
    if ( this->components().test( MESH_ADD_ELEMENTS_INFO ) )
    {
        boost::tie( iv, en ) = this->elementsRange();
        for ( ; iv != en; ++iv )
        {
            value_type meas = 0;
            for( auto const& _elt: iv->pointElementNeighborIds() )
            {
                // warning : only compute meas for active element (no ghost)
                if ( this->hasElement( _elt ) )
                    meas += this->element( _elt ).measure();
            }
            this->elements().modify( iv, [meas]( element_type& e ){ e.setMeasurePointElementNeighbors( meas ); } );
        }
    }
    VLOG(1) << "[Mesh::updateMeasures] update measures : " << ti.elapsed() << "\n"; 

    // compute h information: average, min and max
    {
        ti.restart();
        M_h_avg = 0;
        M_h_min = std::numeric_limits<value_type>::max();
        M_h_max = 0;
        for ( const element_type& elt : allelements( this->shared_from_this() ) )
        {
            M_h_avg += elt.h();
            M_h_min = std::min(M_h_min, elt.h());
            M_h_max = std::max(M_h_max, elt.h());
        }
        M_h_avg /= this->numGlobalElements();
        value_type reduction[3] = { M_h_avg, M_h_min, M_h_max };
        MPI_Op op;
        MPI_Op_create((MPI_User_function*)(Functor::AvgMinMax<value_type, WorldComm::communicator_type>), 1, &op);
        MPI_Allreduce(MPI_IN_PLACE, reduction, 3, mpi::get_mpi_datatype<value_type>(), op, this->worldComm());
        MPI_Op_free(&op);

        M_h_avg = reduction[0];
        M_h_min = reduction[1];
        M_h_max = reduction[2];

        LOG(INFO) << "h average : " << this->hAverage() << "\n";
        LOG(INFO) << "    h min : " << this->hMin() << "\n";
        LOG(INFO) << "    h max : " << this->hMax() << "\n";
        VLOG(1) << "[Mesh::updateForUse] update measures : " << ti.elapsed() << "\n"; 
    }

}


template<typename Shape, typename T, int Tag>
typename Mesh<Shape, T, Tag>::self_type&
Mesh<Shape, T, Tag>::operator+=( self_type const& m )
{
    std::map<std::vector<double>, size_type> mapDel;
    for( auto it = m.beginPoint(), en = m.endPoint();
         it != en;
         ++it )
    {
        bool found = false;
        if( it->isOnBoundary() )
        {
            for( auto pit = this->beginPoint(), pen = this->endPoint();
                 pit != pen;
                 ++pit )
            {
                if ( pit->isOnBoundary() )
                {
                    if ( ublas::norm_2( it->node() - pit->node() ) < 1e-10 )
                    {
                        found = true;
                        std::vector<double> ublasCopy(nDim);
                        std::copy(it->node().begin(), it->node().end(), ublasCopy.begin());
                        // std::cout << "Attention ! Le point " << it->id() << " doit etre remplace par le point " << pit->id() << " (" << it->node() << " ~= " << pit->node() << ")" << std::endl;
                        mapDel.insert(std::pair<std::vector<double>, size_type>(ublasCopy, pit->id()));
                        break;
                    }
                }
            }
        }
        if(!found)
            this->addPoint( *it );
    }
    for( auto it = m.beginFace(), en = m.endFace();it != en;++it )
    {
        // eventually need work
        if ( !it->isOnBoundary() )
        {
            face_type f = *it;
            f.disconnect();
            this->addFace( f );
        }
    }
    for( auto it = m.beginElement(), en = m.endElement();it != en;++it )
    {
        element_type e = *it;
        for( uint16_type p = 0; p < e.numPoints; ++p ) {
            if ( e.point(p).isOnBoundary() ) {
                std::vector<double> ublasCopy(nDim);
                std::copy(e.point(p).node().begin(), e.point(p).node().end(), ublasCopy.begin());
                auto itMap = mapDel.find(ublasCopy);
                if (itMap != mapDel.end() )
                {
                //    std::cout << "Attention ! Le point " << it->point(p).id() << " va etre remplace par le point " << this->point(itMap->second).id() << std::endl;
                    e.setPoint(p, this->point(itMap->second));
                }
            }
        }
        this->addElement( e );
    }
    this->setUpdatedForUse( false );
    this->updateForUse();
    return *this;
}
template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::propagateMarkers( mpl::int_<2> )
{
    // propagate top-down marker from face if edge has not been marked
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
                   } );
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

        DVLOG(2) << "mesh::renumber] element id " << __element.id() <<  " proc " << __element.processId() << "\n";
#if !defined(NDEBUG)
        for ( int i = 0; i < __element.nPoints(); ++i )
        {
            DVLOG(2) << "point id = " << __element.point( i ).id() << "\n";
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

                DVLOG_IF(2,(__id >= next_free_node )) << "next_free_node = " << next_free_node
                                                      << " point id = " << __id << "\n";
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

                DVLOG(2) << "next_free_node = " << next_free_node
                         << " swapping point " << __true_id
                         << " and point " << next_free_node << "\n";

                // next free node
                ++next_free_node;
            } // if
        } // for

        //__element.setId(
    }

    VLOG(2) << "[mesh::renumber] done collecting ids\n";
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
        VLOG(2) << "[mesh::replace] replacing point " << ptit->id()
                      <<  " with " << ptmapit->second.id() << "\n";
#endif
        bool __rep1 = this->points().replace( ptit, ptmapit->second );
        FEELPP_ASSERT( __rep1 )( __rep1 )( ptit->id() )( ptmapit->second.id() ) .warn( "invalid point replacement" );

    }

    VLOG(2) << "[mesh::renumber] done replace point ids\n";

    for ( element_iterator elt = this->beginElement();
            elt != this->endElement(); ++elt )
    {
        element_type __element = *elt;
        DVLOG(2) << "mesh::renumber] element id " << __element.id() <<  " proc " << __element.processId() << "\n";

        // renumber the nodes of the element
        for ( int i = 0; i < __element.nPoints(); ++i )
        {

            size_type __true_id =__element.point( i ).id();
            auto& new_pt = this->point( node_map[__true_id] );
            this->elements().modify( elt,[i,&new_pt]( element_type& e ) { e.setPoint( i, new_pt ); } );

            DVLOG(2) << "point id = " << __true_id << " node map " <<  node_map[__true_id] << " "
                     << "new point id = " << elt->point( i ).id() << "\n";
        }
    }

    VLOG(2) << "[mesh::renumber] done replace point ids in elements\n";

    for ( face_iterator elt = this->beginFace();
            elt != this->endFace(); ++elt )
    {
        face_type __face = *elt;
        DVLOG(2) << "face id: " << __face.id()
                 << " marker: " << __face.marker() << "\n";

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

            DVLOG(2) << "id1= " << __face2.point( 0 ).id() << " id2= " << __face2.point( 1 ).id()<< "\n";
            DVLOG(2) << "point lid = " << i << " id = " << __true_id
                     << " nid = " << this->point( node_map[__true_id] ).id()
                     << " new point id = " << elt->point( i ).id() << "\n";
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

        DVLOG(2) << "edge id: " << __edge.id()
                 << " marker: " << __edge.marker() << "\n";

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

            DVLOG(2) << "renumber edge: id1= " << __edge2.point( 0 ).id() << " id2= " << __edge2.point( 1 ).id()<< "\n";
            DVLOG(2) << "renumber edge: point lid = " << i << " id = " << __true_id
                     << " nid = " << this->point( node_map[__true_id] ).id()
                     << " new point id = " << elt->point( i ).id() << "\n";

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
                auto& new_pt = this->point( __element.point( i ).id() );
                this->elements().modify( elt, [i,&new_pt]( element_type& e ) { e.setPoint( i, new_pt ); } );
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
        LOG(INFO) << "We have " << std::distance( __it, __en ) << " faces in the database";
        auto itb = this->beginFaceOnBoundary();
        auto enb = this->beginFaceOnBoundary();
        LOG(INFO) << "We have " << std::distance( itb, enb ) << " faces on the boundary in the database";
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

            DVLOG_IF(2,faceinserted) << "added face with id " << __it->id () << "\n";
            DVLOG_IF(2,!faceinserted) << "not added face with id " << __it->id ()
                                      << " was already face with id = " << _faceit->second << "\n";


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
                //__f.setIdInPartition( this->worldComm().localRank(),__f.id() );

                DVLOG(2) << "set face id " << __f.id()
                         << " iterator id = " << __it->id()
                         << " check id = " << _faceit->second << "\n";

                this->faces().replace( __it, __f );
                ++__it;
            }
        }
    }

    VLOG(2) << "[Mesh::updateFaces] adding faces : " << ti.elapsed() << "\n";
    ti.restart();

    VLOG(2) << "[Mesh::updateFaces] numLocalFaces : " << this->numLocalFaces() << "\n";
    VLOG(2) << "[Mesh::updateFaces] face_type::numVertices : " << face_type::numVertices << "\n";
    element_iterator iv,  en;
    boost::tie( iv, en ) = this->elementsRange();

    for ( ; iv != en; ++iv )
    {
        element_type const& __element = *iv;
        size_type __element_id = __element.id();

        //MakeBareEntity<element_type,nDim> baremaker( __element );
        for ( size_type j = 0; j < this->numLocalFaces(); j++ )
        {
            DVLOG(2) << "------------------------------------------------------------\n";
            DVLOG(2) << "Element id: " << iv->id() << " local face id: " << j << "\n";

            std::set<int> s;

            for ( int f = 0; f < face_type::numVertices; ++f )
            {
                uint16_type pt_localid = ( nDim==1 )?j:iv->fToP( j, f );
                s.insert( iv->point( pt_localid ).id() );
                VLOG(3) << "add point local id " << f << " to face " << j  << " " << iv->fToP( j, f )
                        << " global id " << iv->point( pt_localid ).id() << "\n";
            }

            bool faceinserted = false;
            boost::tie( _faceit, faceinserted ) = _faces.insert( std::make_pair( s, next_face ) );

            if ( faceinserted )
                ++next_face;

            DVLOG(2) << "------------------------------------------------------------\n";
            DVLOG(2) << "Element id: " << iv->id() << " local face id: " << j << "\n";

            if ( faceinserted )
            {
                DVLOG(2) << "creating the face:" << _faceit->second << "\n";

                // set face id
                face.setId( _faceit->second );
                face.disconnect();
                if ( this->components().test( MESH_ADD_ELEMENTS_INFO ) )
                    face.addElement( __element_id );

                // set the process id from element
                face.setProcessId( __element.processId() );

                // set the id for this partition
                //face.setIdInPartition( this->worldComm().localRank(),face.id() );

                // set the vertices of the face
                for ( size_type k = 0; k < face_type::numPoints; ++k )
                    face.setPoint( k, iv->point( ele.fToP( j, k ) ) );

                // set the connection with the element
                face.setConnection0( boost::make_tuple( boost::addressof( __element ), __element_id, j, __element.processId() ) );
                face.setOnBoundary( true, face_type::nDim );

                // adding the face
                bool inserted = false;
                face_iterator __fit;
                boost::tie( __fit, inserted )  = this->addFace( face );
                FEELPP_ASSERT( inserted && __fit != this->endFace() )
                ( _faceit->second )
                ( iv->id() )
                ( __fit->id() )
                ( face.id() ).error( "invalid face iterator" );
                this->elements().modify( iv, Feel::detail::UpdateFace<face_type>( boost::cref( *__fit ) ) );


                DVLOG(2) << "Adding [new] face info : \n";
                DVLOG(2) << "element id: " << __element_id << "\n";
                DVLOG(2) << "process id: " << __fit->processId() << "\n";
                DVLOG(2) << "id: " << __fit->id() << "\n";
                DVLOG(2) << "bdy: " << __fit->isOnBoundary() << "\n";
                DVLOG(2) << "marker: " << __fit->marker() << "\n";
                DVLOG(2) << "ad_first: " << __fit->ad_first() << "\n";
                DVLOG(2) << "pos_first: " << __fit->pos_first() << "\n";
                DVLOG(2) << "proc_first: " << __fit->proc_first() << "\n";
                DVLOG(2) << "ad_second: " << __fit->ad_second() << "\n";
                DVLOG(2) << "pos_second: " << __fit->pos_second() << "\n";
                DVLOG(2) << "proc_second: " << __fit->proc_second() << "\n";

            }

            else
            {
                DVLOG(2) << "found the face:" << _faceit->second << " in element " << __element_id << " and local face: " << j << "\n";

                // look in the face table for the face
                face_iterator __fit = this->faces().find( face_type( _faceit->second ) );
                FEELPP_ASSERT( __fit != this->endFace() )( _faceit->second ).error( "face is not in face container" );


                face_type face = *__fit;
                DVLOG(2) << "the face id :"  << __fit->id() << "\n";

                if ( this->components().test( MESH_ADD_ELEMENTS_INFO ) )
                    face.addElement( __element_id );

                // the three conditions below typically arise after reading a serialized mesh
                if ( __fit->isConnectedTo0() && __fit->connection0().template get<0>() == 0 && ( __element.id() == __fit->ad_first() ) )
                {
                    DVLOG(2) << "fixing connection 0 in face\n";
                    // reconnect the elements to the face
                    auto connect0 = __fit->connection0();
                    connect0.template get<0>() = boost::addressof( __element );
                    face.setConnection0( connect0 );

                    // need to reconstruct the neighbors
                    this->elements().modify( iv, Feel::detail::UpdateFace<face_type>( boost::cref( *__fit ) ) );

                }
                if ( __fit->isConnectedTo1() && __fit->connection1().template get<0>() == 0 && ( __element.id() == __fit->ad_second() ) )
                {
                    DVLOG(2) << "fixing connection 1 in face\n";

                    // reconnect the elements to the face
                    auto connect1 = __fit->connection1();
                    connect1.template get<0>() = boost::addressof( __element );
                    face.setConnection1( connect1 );

                    // need to reconstruct the neighbors
                    element_iterator elt1 = this->elementIterator( __fit->ad_first(), __fit->proc_first() );
                    this->elements().modify( elt1, Feel::detail::UpdateFace<face_type>( boost::cref( *__fit ) ) );
                    this->elements().modify( iv, Feel::detail::UpdateFace<face_type>( boost::cref( *__fit ) ) );

                }
                if ( __fit->isConnectedTo0() && __fit->isConnectedTo1() )
                {
                    DVLOG(2) << "internal face, fixing process id if necessary\n";
                    if ( !iv->facePtr( j ) )
                        this->elements().modify( iv, Feel::detail::UpdateFace<face_type>( boost::cref( *__fit ) ) );
                    FEELPP_ASSERT( iv->facePtr( j ) )( j )( iv->id() ).warn( "invalid element face error" );
                    FEELPP_ASSERT( face.isConnectedTo0() && face.isConnectedTo1() )
                        ( face.isConnectedTo0() )( face.isConnectedTo1() ).error ("inconsistent data structure" );
                    if ( face.processId()!=this->worldComm().localRank() )
                    {
                        if ( ( face.element0().processId()==this->worldComm().localRank() ) ||
                             ( face.element1().processId()==this->worldComm().localRank() ) )
                            face.setProcessId( this->worldComm().localRank() );
                    }
                }


                // the face could have been entered apriori given by
                // the mesh generator, so just set the connection0
                // properly .
                if ( !__fit->isConnectedTo0() )
                {
                    DVLOG(2) << "[updateFaces][boundary] element: " << __element_id
                             << " face: " << j << " id: " << _faceit->second << "\n";

                    // set the connection with the element
                    face.setConnection0( boost::make_tuple( boost::addressof( __element ), __element_id, j, __element.processId() ) );
                    // set the process id from element
                    face.setProcessId( __element.processId() );
                    face.setOnBoundary( true, face_type::nDim );
                    //this->faces().modify( __fit,
                    //Feel::detail::UpdateFaceConnection0<typename face_type::element_connectivity_type>( boost::make_tuple( boost::addressof( __element ), __element_id, j, __element.processId() ) ) );

                    this->faces().replace( __fit, face );

                    this->elements().modify( iv, Feel::detail::UpdateFace<face_type>( boost::cref( *__fit ) ) );

                    DVLOG(2) << "adding [!isConnectedTo0] face info : \n";
                    DVLOG(2) << "id: " << __fit->id() << "\n";
                    DVLOG(2) << "process id: " << __fit->processId() << "\n";
                    DVLOG(2) << "bdy: " << __fit->isOnBoundary() << "\n";
                    DVLOG(2) << "marker: " << __fit->marker() << "\n";
                    DVLOG(2) << "ad_first: " << __fit->ad_first() << "\n";
                    DVLOG(2) << "pos_first: " << __fit->pos_first() << "\n";
                    DVLOG(2) << "proc_first: " << __fit->proc_first() << "\n";
                    DVLOG(2) << "ad_second: " << __fit->ad_second() << "\n";
                    DVLOG(2) << "pos_second: " << __fit->pos_second() << "\n";
                    DVLOG(2) << "proc_second: " << __fit->proc_second() << "\n";
                    DVLOG(2) << "element process id: " << iv->processId() << "\n";
                }

                // we found an internal face
                else if ( !__fit->isConnectedTo1() )
                {


#if 0
                    this->faces().modify( __fit,
                                          Feel::detail::UpdateFaceConnection1<typename face_type::element_connectivity_type>( boost::make_tuple( boost::addressof( __element ), __element_id, j, __element.processId() ) ) );


                    FEELPP_ASSERT( __fit->isConnectedTo0() && __fit->isConnectedTo1() )
                    ( __fit->isConnectedTo0() )( __fit->isConnectedTo1() ).error( "invalid face connection" );
#endif
                    Feel::detail::UpdateFaceConnection1<typename face_type::element_connectivity_type> update1( boost::make_tuple( boost::addressof( __element ), __element_id, j, __element.processId() ) );
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
                    this->elements().modify( elt1, Feel::detail::UpdateFace<face_type>( boost::cref( *__fit ) ) );
                    this->elements().modify( iv, Feel::detail::UpdateFace<face_type>( boost::cref( *__fit ) ) );

                    // fix duplication of point in connection1 with 3d mesh at order 3 and 4
                    this->fixPointDuplicationInHOMesh( iv,__fit, mpl::bool_< nDim == 3 && nOrder >= 3 >() );

                    DVLOG(2) << "adding face info : \n";
                    DVLOG(2) << "id: " << __fit->id() << "\n";
                    DVLOG(2) << "process id: " << __fit->processId() << "\n";
                    DVLOG(2) << "bdy: " << __fit->isOnBoundary() << "\n";
                    DVLOG(2) << "marker: " << __fit->marker() << "\n";
                    DVLOG(2) << "ad_first: " << __fit->ad_first() << "\n";
                    DVLOG(2) << "pos_first: " << __fit->pos_first() << "\n";
                    DVLOG(2) << "proc_first: " << __fit->proc_first() << "\n";
                    DVLOG(2) << "ad_second: " << __fit->ad_second() << "\n";
                    DVLOG(2) << "pos_second: " << __fit->pos_second() << "\n";
                    DVLOG(2) << "proc_second: " << __fit->proc_second() << "\n";
                    DVLOG(2) << "element1 process id: " << elt1->processId() << "\n";
                    DVLOG(2) << "element2 process id: " << iv->processId() << "\n";
                }

                CHECK( (__fit->processId() == __fit->proc_first()) ||
                       (__fit->processId() == __fit->proc_second()) )
                    << "invalid process id " << __fit->processId() << " with element proc first = " <<  __fit->proc_first()
                    << " and element proc second " << __fit->proc_second();
            }

            LOG_IF( WARNING, !iv->facePtr( j ) ) << "invalid element " << iv->id() << " with local face id " << j;
        } // face loop
    } // element loop


    face_iterator f_it = this->beginFace();
    face_iterator f_en = this->endFace();

    for ( ; f_it!=f_en;  )
    {
        // cleanup the face data structure :

        if ( !f_it->isConnectedTo0() )
        {
            DLOG(INFO) << "removing face id : " << f_it->id()
                       << " marker : " << f_it->marker();
            // remove all faces that are not connected to any elements
            f_it = this->faces().erase( f_it );
            //++f_it;
        }
        else
        {
            ++f_it;
        }

    }

    LOG(INFO) << "We have now " << nelements(boundaryfaces(this)) << " faces on the boundary in the database";
    VLOG(1) << "[Mesh::updateEntitiesCoDimensionOne] element/face connectivity : " << ti.elapsed() << "\n";
    ti.restart();
}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::fixPointDuplicationInHOMesh( element_iterator iv,face_iterator __fit, mpl::true_ )
{
    CHECK( nOrder == 3 || nOrder == 4 ) << "this fix with nOrder " << nOrder << " not implement, only order 3 and 4\n";
    if ( nOrder == 3 )
    {
        uint16_type startPt = face_type::numVertices*face_type::nbPtsPerVertex + face_type::numEdges*face_type::nbPtsPerEdge;
        CHECK( startPt == 9 && face_type::numPoints==10 ) << "invalid case\n";

        // easy case : only one point in face
        uint16_type idx = element_type::fToP( __fit->pos_second(), startPt );
        auto const& ptOnFace = __fit->point( startPt );
        size_type ptIdToCheck = iv->point( idx ).id();
        if ( ptIdToCheck != ptOnFace.id() )
        {
            // erase duplicate point
            this->points().erase( this->pointIterator( ptIdToCheck ) );
            // fix point coordinate
            this->elements().modify( iv, Feel::detail::UpdateEltPoint<point_type>( idx, ptOnFace ) );
        }
    }
    else if ( nOrder == 4 )
    {
        const uint16_type startPt = face_type::numVertices*face_type::nbPtsPerVertex + face_type::numEdges*face_type::nbPtsPerEdge;
        CHECK( startPt == 12 && face_type::numPoints==15 ) << "invalid case\n";

        // get a correspondance between the face and elt point numbering
        std::set<uint16_type> idxDone;
        std::vector<uint16_type> idxDistribution(3,invalid_uint16_type_value);
        for ( uint16_type k = startPt; k < face_type::numPoints; ++k )
        {
            auto const& ptOnFaceElt = iv->point( element_type::fToP( __fit->pos_second(), k ) );

            double dist=0, distMin= 1000*iv->hFace( __fit->pos_second() );// INT_MAX;
            uint16_type idxFaceNear = invalid_uint16_type_value;

            for ( uint16_type idxV = startPt ; idxV < face_type::numPoints ; ++idxV )
            {
                auto const& ptOnFace = __fit->point( idxV );
                dist=0;
                for ( uint16_type d=0 ; d< nRealDim ; ++d )
                    dist += std::pow( ptOnFace(d)- ptOnFaceElt(d),2);
                if ( dist<distMin ) { distMin=dist; idxFaceNear=idxV; }
            }

            CHECK( idxDone.find( idxFaceNear ) == idxDone.end() ) << " idxVSearch is already done\n";
            idxDone.insert( idxFaceNear );
            idxDistribution[k-12] = idxFaceNear;
        }

        // fix point if necessary
        for ( uint16_type i = 0;i<3;++i)
        {
            // point index in element
            uint16_type idx = element_type::fToP( __fit->pos_second(), startPt + i );
            size_type ptIdToCheck = iv->point( idx ).id();
            auto const& ptOnFace = __fit->point( idxDistribution[i] );

            if ( ptIdToCheck != ptOnFace.id() )
            {
                // erase duplicate point
                this->points().erase( this->pointIterator( ptIdToCheck ) );
                // fix point coordinate
                this->elements().modify( iv, Feel::detail::UpdateEltPoint<point_type>( idx, ptOnFace ) );
            }
        }
    }

}
template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::fixPointDuplicationInHOMesh( element_iterator iv, face_iterator __fit, mpl::false_ )
{}


template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::modifyEdgesOnBoundary( face_iterator& it , mpl::bool_<true> )
{
    // loop over face edges
    for ( int f = 0; f < face_type::numEdges; ++f )
    {
        if ( it->edge( f ).isOnBoundary() == false )
        {
            auto eit = this->edgeIterator( it->edge(f).id() );
            this->edges().modify( eit,
                                  []( edge_type& e )
                                  {e.setOnBoundary(true, 1 );} );
        }

    }
    DVLOG(3) << "We have " << nelements(boundaryedges(this)) <<  " boundary edges";

}
template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::modifyEdgesOnBoundary( face_iterator& f, mpl::bool_<false> )
{
}

template<typename Shape, typename T, int Tag>
bool
Mesh<Shape, T, Tag>::modifyElementOnBoundaryFromEdge( element_iterator& e, mpl::bool_<false> )
{
    return false;
}
template<typename Shape, typename T, int Tag>
bool
Mesh<Shape, T, Tag>::modifyElementOnBoundaryFromEdge( element_iterator& iv, mpl::bool_<true> )
{
    // in 3D check if the edges of the element touch the boundary
    bool isOnBoundary = false;
    for ( size_type j = 0; j < iv->nEdges(); j++ )
    {
        isOnBoundary |= iv->edge( j ).isOnBoundary();
    }
    if ( isOnBoundary )
    {
        // the element touches the boundary with just an edge
        this->elements().modify( iv, [isOnBoundary]( element_type& e ) { e.setOnBoundary( isOnBoundary, 1 ); } );
        return isOnBoundary;
    }
    return isOnBoundary;
}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::updateOnBoundary()
{
    // first go through all the faces and set the points of the boundary
    // faces to be on the boundary
    LOG(INFO) << "update boundary points...";
    LOG(INFO) << "Initially we have " << nelements(boundarypoints(this)) <<  " boundary points";
    for( auto it = this->beginFace(), en = this->endFace(); it != en; ++ it )
    {
        if ( it->isOnBoundary() == true )
        {
            modifyEdgesOnBoundary( it, mpl::bool_<(nDim >= 3)>() );
            // loop over face points
            for ( int f = 0; f < face_type::numPoints; ++f )
            {
                if ( it->point( f ).isOnBoundary() == false )
                {
                    auto pit = this->pointIterator( it->point(f).id() );
                    this->points().modify( pit,
                                           []( point_type& p )
                                           {p.setOnBoundary(true, 0 );} );
                }

            }
        }
    }
    LOG(INFO) << "We have now " << nelements(boundarypoints(this)) <<  " boundary points";
    LOG(INFO) << "update boundary elements...";
    // loop through faces to set the elements having a face on the boundary
    for ( auto iv = this->beginElement(), en = this->endElement();
          iv != en; ++iv )
    {
        bool isOnBoundary = false;

        // first check if a face is on the boundary
        for ( size_type j = 0; j < iv->nTopologicalFaces(); j++ )
        {
            isOnBoundary |= iv->face( j ).isOnBoundary();
        }

        if ( isOnBoundary )
        {
            VLOG(3) << "checking " << iv->nTopologicalFaces() << " faces, isOnBoundary: " << isOnBoundary << " face_type::nDim: " << face_type::nDim;
            this->elements().modify( iv, [isOnBoundary]( element_type& e ) { e.setOnBoundary( isOnBoundary, face_type::nDim ); } );
            // go to the next element, no need to look further
            continue;
        }
        bool e_modified = modifyElementOnBoundaryFromEdge( iv, mpl::bool_<(nDim>=3)>() );
        // go to next element if element is on boundary
        if ( e_modified )
        {
            continue;
        }

        // finally check if a point of the element touches the boundary
        for ( size_type j = 0; j < iv->nPoints(); j++ )
        {
            isOnBoundary |= iv->point( j ).isOnBoundary();
        }

        if ( isOnBoundary )
        {
            VLOG(3) << "checking " << iv->nPoints() << " points, isOnBoundary: " << isOnBoundary;
            this->elements().modify( iv, [isOnBoundary]( element_type& e ) { e.setOnBoundary( isOnBoundary, 0 ); } );
        }
    } // loop over the elements
    LOG(INFO) << "[updateOnBoundary] We have " << nelements(boundaryelements(this))
              << " elements sharing a point, a edge or a face with the boundary in the database";
    BOOST_FOREACH( auto e, this->boundaryElements( 0, 2, 0 ) )
    {
        VLOG(3) << "boundary element : " << e.id()
                << " entity on boundary max dim  " << e.boundaryEntityDimension()
                << " process id : " << e.processId();
    }

}
template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::removeFacesFromBoundary( std::initializer_list<uint16_type> markers )
{
    std::for_each( markers.begin(), markers.end(),
                   [=]( uint16_type marker )
                   {
                       auto range=markedfaces( this->shared_from_this(), boost::any(marker) );
                       LOG(INFO) << "removing " << nelements(range) << " faces marked "  << marker << " from boundary faces\n";

                       for( auto it = range.template get<1>(), en = range.template get<2>(); it != en; ++it )
                       {
                           if ( it->isOnBoundary() )
                           {
                               DVLOG(3) << "removing face "  << it->id() << "\n";
                               auto it2 = this->faces().template project<0>( it );
                               this->faces().modify( it2, []( face_type & f ) { f.setOnBoundary( false ); } );


                               //face_type f = *it2;
                               //f.setOnBoundary( false );
                               //this->faces().replace( it2, f );
                               CHECK( it2->isOnBoundary() ==false ) << " face should not be on the boundary anymore\n";
                               CHECK( it->isOnBoundary() ==false ) << " face should not be on the boundary anymore\n";
                           }

                       }
                   } );
}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::updateEntitiesCoDimensionOneGhostCellByUsingBlockingComm()
{
    typedef std::vector< boost::tuple<size_type, std::vector<double> > > resultghost_type;

    VLOG(2) << "[Mesh::updateEntitiesCoDimensionOneGhostCell] start on god rank "<< this->worldComm().godRank() << "\n";

    std::vector<int> nbMsgToSend( this->worldComm().localSize(), 0 );
    std::vector<int> nbMsgToRecv( this->worldComm().localSize(), 0 );
    std::vector< std::map<int,int> > mapMsg( this->worldComm().localSize() );

    auto iv = this->beginGhostElement();
    auto en = this->endGhostElement();
    for ( ; iv != en; ++iv )
    {
        element_type const& __element = *iv;
        const int IdProcessOfGhost = __element.processId();
        const size_type idInPartition = __element.idInOthersPartitions( IdProcessOfGhost );

        this->elements().modify( this->elementIterator(iv->id(),IdProcessOfGhost) , typename super_elements::ElementGhostConnectPointToElement() );
        if ( nDim==3 )
            this->elements().modify( this->elementIterator(iv->id(),IdProcessOfGhost) , typename super_elements::ElementGhostConnectEdgeToElement() );

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

    auto itEltActif = this->beginElementWithProcessId( this->worldComm().localRank() );
    auto const enEltActif = this->endElementWithProcessId( this->worldComm().localRank() );
    for ( ; itEltActif!=enEltActif ; ++itEltActif )
    {
        if (itEltActif->numberOfNeighborPartitions() == 0 ) continue;
        auto itneighbor = itEltActif->neighborPartitionIds().begin();
        auto const enneighbor = itEltActif->neighborPartitionIds().end();
        for ( ; itneighbor!=enneighbor ; ++itneighbor )
            nbMsgToRecv[*itneighbor]++;
    }

    //------------------------------------------------------------------------------------------------//

#if !defined( NDEBUG )
    // check nbMsgToRecv computation
    std::vector<int> nbMsgToRecv2;
    mpi::all_to_all( this->worldComm().localComm(),
                     nbMsgToSend,
                     nbMsgToRecv2 );
    for ( int proc=0; proc<this->worldComm().localSize(); ++proc )
    {
        CHECK( nbMsgToRecv[proc]==nbMsgToRecv2[proc] ) << "paritioning data incorect "
                                                       << "myrank " << this->worldComm().localRank() << " proc " << proc
                                                       << " nbMsgToRecv[proc] " << nbMsgToRecv[proc]
                                                       << " nbMsgToRecv2[proc] " << nbMsgToRecv2[proc] << "\n";
    }
#endif

    //------------------------------------------------------------------------------------------------//
    // recv id asked and re-send set of face id
    for ( int proc=0; proc<this->worldComm().localSize(); ++proc )
    {
        for ( int cpt=0; cpt<nbMsgToRecv[proc]; ++cpt )
        {
            //recv
            size_type idRecv;
            this->worldComm().localComm().recv( proc, cpt, idRecv );
#if 0
            std::cout<< "I am the proc" << this->worldComm().localRank()<<" I receive to proc " << proc
                     <<" with tag "<< cpt << " idRecv " << idRecv
                     << " it_ghost->G() " << this->element( idRecv ).G()
                     << std::endl;
#endif
            auto const& theelt = this->element( idRecv );


            // get faces id and bary
            resultghost_type idFacesWithBary(this->numLocalFaces(),boost::make_tuple(0, std::vector<double>(nRealDim) ));
            for ( size_type j = 0; j < this->numLocalFaces(); j++ )
            {
                auto const& theface = theelt.face( j );
                idFacesWithBary[j].template get<0>() = theface.id();
                //auto const& theGj = theface.G();
                auto const& theGj = theface.vertices();
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
                for ( uint16_type comp = 0; comp<nRealDim; ++comp )
                {
                    idFacesWithBary[j].template get<1>()[comp]=baryFace[comp];
                }
            }

            // get points id and nodes
            resultghost_type idPointsWithNode(element_type::numLocalVertices,boost::make_tuple(0, std::vector<double>(nRealDim) ));
            for ( size_type j = 0; j < element_type::numLocalVertices; j++ )
            {
                auto const& thepoint = theelt.point( j );
                idPointsWithNode[j].template get<0>() = thepoint.id();
                for ( uint16_type comp = 0; comp<nRealDim; ++comp )
                {
                    idPointsWithNode[j].template get<1>()[comp]=thepoint(comp);
                }
            }

            std::vector<resultghost_type> theresponse(2);
            theresponse[0] = idPointsWithNode;
            theresponse[1] = idFacesWithBary;
            //auto theresponse = boost::make_tuple( idPointsWithNode, idFacesWithBary );
            // send response
            //this->worldComm().localComm().send( proc, cpt, idFacesWithBary );
            this->worldComm().localComm().send( proc, cpt, theresponse );
        }
    }

    //------------------------------------------------------------------------------------------------//
    // get response to initial request and update Feel::Mesh::Faces data
    for ( int proc=0; proc<this->worldComm().localSize(); ++proc )
    {
        for ( int cpt=0; cpt<nbMsgToSend[proc]; ++cpt )
        {
            //recv
#if 0
            std::vector< boost::tuple<size_type, std::vector<double> > > idFacesWithBaryRecv(this->numLocalFaces());
            this->worldComm().localComm().recv( proc, cpt, idFacesWithBaryRecv );
#elif 0
            typedef std::vector< boost::tuple<size_type, std::vector<double> > > resultghost_face_type;
            typedef std::vector< boost::tuple<size_type, std::vector<double> > > resultghost_point_type;
            boost::tuple<resultghost_point_type,resultghost_face_type> requestRecv;
            this->worldComm().localComm().recv( proc, cpt, requestRecv );
            auto const& idPointsWithNodeRecv = requestRecv.template get<0>();
            auto const& idFacesWithBaryRecv = requestRecv.template get<1>();
#else
            std::vector<resultghost_type> requestRecv;
            this->worldComm().localComm().recv( proc, cpt, requestRecv );
            auto const& idPointsWithNodeRecv = requestRecv[0];
            auto const& idFacesWithBaryRecv = requestRecv[1];
#endif

#if 0
            std::cout<< "I am the proc " << this->worldComm().localRank()<<" I receive to proc " << proc
                     <<" with tag "<< cpt << std::endl;
#endif
            auto const& theelt = this->element( mapMsg[proc][cpt],proc );

            //update faces data
            for ( size_type j = 0; j < this->numLocalFaces(); j++ )
            {
                auto const& idFaceRecv = idFacesWithBaryRecv[j].template get<0>();
                auto const& baryFaceRecv = idFacesWithBaryRecv[j].template get<1>();

                //objective : find  face_it (hence jBis in theelt ) (permutations would be necessary)
                uint16_type jBis = invalid_uint16_type_value;
                bool hasFind=false;
                for ( uint16_type j2 = 0; j2 < this->numLocalFaces() && !hasFind; j2++ )
                {

                    auto const& thefacej2 = theelt.face( j2 );
                    //auto const& theGj2 = thefacej2.G();
                    auto const& theGj2 = thefacej2.vertices();
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
                    // compare barycenters
                    bool find2=true;
                    for (uint16_type d=0;d<nRealDim;++d)
                        {
                            find2 = find2 && ( std::abs( baryFace[d]-baryFaceRecv[d] )<1e-9 );
                        }
                    if (find2) { hasFind=true;jBis=j2; }

                } //for ( uint16_type j2 = 0; j2 < this->numLocalFaces() && !hasFind; j2++ )

                CHECK ( hasFind ) << "[mesh::updateEntitiesCoDimensionOneGhostCell] : invalid partitioning data, ghost face cells are not available\n";

                // get the good face
                auto face_it = this->faceIterator( theelt.face( jBis ).id() );
                //update the face
                this->faces().modify( face_it, Feel::detail::updateIdInOthersPartitions( proc, idFaceRecv ) );

            } // for ( size_type j = 0; j < this->numLocalFaces(); j++ )


            for ( size_type j = 0; j < element_type::numLocalVertices; j++ )
            {
                auto const& idPointRecv = idPointsWithNodeRecv[j].template get<0>();
                auto const& nodePointRecv = idPointsWithNodeRecv[j].template get<1>();

                uint16_type jBis = invalid_uint16_type_value;
                bool hasFind=false;
                for ( uint16_type j2 = 0; j2 < element_type::numLocalVertices && !hasFind; j2++ )
                {
                    auto const& thepointj2 = theelt.point( j2 );
                    // compare barycenters
                    bool find2=true;
                    for (uint16_type d=0;d<nRealDim;++d)
                        {
                            find2 = find2 && ( std::abs( thepointj2(d)-nodePointRecv[d] )<1e-9 );
                        }
                    if (find2) { hasFind=true;jBis=j2; }
                }

                CHECK ( hasFind ) << "[mesh::updateEntitiesCoDimensionOneGhostCell] : invalid partitioning data, ghost point cells are not available\n";
                // get the good face
                auto point_it = this->pointIterator( theelt.point( jBis ).id() );
                //update the face
                this->points().modify( point_it, Feel::detail::updateIdInOthersPartitions( proc, idPointRecv ) );
            } // for ( size_type j = 0; j < element_type::numLocalVertices; j++ )

            /*for ( size_type j = 0; j < element_type::numLocalEdges; j++ )
            {
            }*/

        } // for ( int cpt=0;cpt<nbMsgToSend[proc];++cpt)
    } // for (int proc=0; proc<M_comm.size();++proc)

    //------------------------------------------------------------------------------------------------//

    //std::cout << "[Mesh::updateEntitiesCoDimensionOneGhostCell] finish" << std::endl;

} // updateEntitiesCoDimensionOneGhostCell

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::updateEntitiesCoDimensionOneGhostCellByUsingNonBlockingComm()
{
    typedef std::vector< boost::tuple<size_type, std::vector<double> > > resultghost_point_type;
    typedef std::vector< boost::tuple<size_type, bool, std::vector<double> > > resultghost_face_type;

    DVLOG(1) << "updateEntitiesCoDimensionOneGhostCellByUsingNonBlockingComm : start on rank " << this->worldComm().localRank() << "\n";

    const rank_type nProc = this->worldComm().localSize();

    //------------------------------------------------------------------------------------------------//
    // compute size of container to send and update Point and Edge info for parallelism
    std::map< rank_type, int > nDataInVecToSend;
    auto iv = this->beginGhostElement();
    auto const en = this->endGhostElement();
    for ( ; iv != en; ++iv )
    {
        const rank_type IdProcessOfGhost = iv->processId();
        // update info for parallelism
        this->elements().modify( this->elementIterator(iv->id(),IdProcessOfGhost) , typename super_elements::ElementGhostConnectPointToElement() );
        if ( nDim==3 )
            this->elements().modify( this->elementIterator(iv->id(),IdProcessOfGhost) , typename super_elements::ElementGhostConnectEdgeToElement() );

        // be sure that the counter start to 0
        if ( nDataInVecToSend.find(IdProcessOfGhost) == nDataInVecToSend.end() )
            nDataInVecToSend[IdProcessOfGhost]=0;
        // update counter
        nDataInVecToSend[IdProcessOfGhost]++;
    }
    //------------------------------------------------------------------------------------------------//
    // init and resize the container to send
    std::map< rank_type, std::vector<size_type> > dataToSend;
    auto itNDataInVecToSend = nDataInVecToSend.begin();
    auto const enNDataInVecToSend = nDataInVecToSend.end();
    for ( ; itNDataInVecToSend!=enNDataInVecToSend ; ++itNDataInVecToSend )
    {
        const rank_type idProc = itNDataInVecToSend->first;
        const int nData = itNDataInVecToSend->second;
        dataToSend[idProc].resize( nData );
    }
    //------------------------------------------------------------------------------------------------//
    // prepare container to send
    std::map< rank_type, std::map<int,size_type> > memoryMsgToSend;
    std::map< rank_type, int > nDataInVecToSendBis;
    iv = this->beginGhostElement();
    for ( ; iv != en; ++iv )
    {
        element_type const& __element = *iv;
        const rank_type idProc = __element.processId();
        const size_type idEltInOtherPartition = __element.idInOthersPartitions( idProc );
        // be sure that the counter start to 0
        if ( nDataInVecToSendBis.find(idProc) == nDataInVecToSendBis.end() )
            nDataInVecToSendBis[idProc]=0;
        // save request
        memoryMsgToSend[idProc][nDataInVecToSendBis[idProc]] = __element.id();
        // update container
        dataToSend[idProc][nDataInVecToSendBis[idProc]] = idEltInOtherPartition;
        // update counter
        nDataInVecToSendBis[idProc]++;
    }
    //------------------------------------------------------------------------------------------------//
    // compute nbMsgToRecv
    std::set<rank_type> procToRecv;
    auto itEltActif = this->beginElementWithProcessId( this->worldComm().localRank() );
    auto const enEltActif = this->endElementWithProcessId( this->worldComm().localRank() );
    for ( ; itEltActif!=enEltActif ; ++itEltActif )
    {
        if (itEltActif->numberOfNeighborPartitions() == 0 ) continue;
        auto itneighbor = itEltActif->neighborPartitionIds().begin();
        auto const enneighbor = itEltActif->neighborPartitionIds().end();
        for ( ; itneighbor!=enneighbor ; ++itneighbor )
            procToRecv.insert(*itneighbor);
    }
    //------------------------------------------------------------------------------------------------//
    // counter of request
    int nbRequest=0;
    for ( rank_type proc=0; proc<nProc; ++proc )
    {
        if ( dataToSend.find(proc) != dataToSend.end() )
            ++nbRequest;
        if ( procToRecv.find(proc) != procToRecv.end() )
            ++nbRequest;
    }
    if ( nbRequest == 0 ) return;

    mpi::request * reqs = new mpi::request[nbRequest];
    int cptRequest=0;
    //------------------------------------------------------------------------------------------------//
    // first send
    auto itDataToSend = dataToSend.begin();
    auto const enDataToSend = dataToSend.end();
    for ( ; itDataToSend!=enDataToSend ; ++itDataToSend )
    {
        reqs[cptRequest] = this->worldComm().localComm().isend( itDataToSend->first , 0, itDataToSend->second );
        ++cptRequest;
    }
    //------------------------------------------------------------------------------------------------//
    // first recv
    std::map<rank_type,std::vector<size_type> > dataToRecv;
    auto itProcToRecv = procToRecv.begin();
    auto const enProcToRecv = procToRecv.end();
    for ( ; itProcToRecv != enProcToRecv ; ++itProcToRecv )
    {
        const rank_type idProc = *itProcToRecv;
        reqs[cptRequest] = this->worldComm().localComm().irecv( idProc , 0, dataToRecv[idProc] );
        ++cptRequest;
    }
    //------------------------------------------------------------------------------------------------//
    // wait all requests
    mpi::wait_all(reqs, reqs + nbRequest);
    //------------------------------------------------------------------------------------------------//
    // build the container to ReSend
    std::map<rank_type, std::vector< boost::tuple<resultghost_point_type,resultghost_face_type>  > > dataToReSend;
    auto itDataRecv = dataToRecv.begin();
    auto const enDataRecv = dataToRecv.end();
    for ( ; itDataRecv!=enDataRecv ; ++itDataRecv )
    {
        const rank_type idProc = itDataRecv->first;
        const int nDataRecv = itDataRecv->second.size();
        dataToReSend[idProc].resize( nDataRecv );
        for ( int k=0; k<nDataRecv; ++k )
        {
            auto const& theelt = this->element( itDataRecv->second[k] );
            //---------------------------//
            // get faces id and bary
            resultghost_face_type idFacesWithBary(this->numLocalFaces(),boost::make_tuple(0, false, std::vector<double>(nRealDim) ));
            for ( size_type j = 0; j < this->numLocalFaces(); j++ )
            {
                auto const& theface = theelt.face( j );
                idFacesWithBary[j].template get<0>() = theface.id();
                idFacesWithBary[j].template get<1>() = theface.isOnBoundary();

                auto const& theGj = theface.vertices();//theface.G();
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
                for ( uint16_type comp = 0; comp<nRealDim; ++comp )
                {
                    idFacesWithBary[j].template get<2>()[comp]=baryFace[comp];
                }
            }
            //---------------------------//
            // get points id and nodes
            resultghost_point_type idPointsWithNode(element_type::numLocalVertices,boost::make_tuple(0, std::vector<double>(nRealDim) ));
            for ( size_type j = 0; j < element_type::numLocalVertices; j++ )
            {
                auto const& thepoint = theelt.point( j );
                idPointsWithNode[j].template get<0>() = thepoint.id();
                for ( uint16_type comp = 0; comp<nRealDim; ++comp )
                {
                    idPointsWithNode[j].template get<1>()[comp]=thepoint(comp);
                }
            }
            //---------------------------//
            // update container to ReSend
            dataToReSend[idProc][k] = boost::make_tuple(idPointsWithNode,idFacesWithBary);

        } // for ( int k=0; k<nDataRecv; ++k )
    }
    //------------------------------------------------------------------------------------------------//
    // send respond to the request
    cptRequest=0;
    auto itDataToReSend = dataToReSend.begin();
    auto const enDataToReSend = dataToReSend.end();
    for ( ; itDataToReSend!=enDataToReSend ; ++itDataToReSend )
    {
        reqs[cptRequest] = this->worldComm().localComm().isend( itDataToReSend->first , 0, itDataToReSend->second );
        ++cptRequest;
    }
    //------------------------------------------------------------------------------------------------//
    // recv the initial request
    std::map<rank_type, std::vector< boost::tuple<resultghost_point_type,resultghost_face_type>  > > finalDataToRecv;

    itDataToSend = dataToSend.begin();
    for ( ; itDataToSend!=enDataToSend ; ++itDataToSend )
    {
        const rank_type idProc = itDataToSend->first;
        reqs[cptRequest] = this->worldComm().localComm().irecv( idProc, 0, finalDataToRecv[idProc] );
        ++cptRequest;
    }
    //------------------------------------------------------------------------------------------------//
    // wait all requests
    mpi::wait_all(reqs, reqs + nbRequest);
    // delete reqs because finish comm
    delete [] reqs;
    //------------------------------------------------------------------------------------------------//
    // update mesh : id in other partitions for the ghost cells
    auto itFinalDataToRecv = finalDataToRecv.begin();
    auto const enFinalDataToRecv = finalDataToRecv.end();
    for ( ; itFinalDataToRecv!=enFinalDataToRecv ; ++itFinalDataToRecv)
    {
        const rank_type idProc = itFinalDataToRecv->first;
        const int nDataRecv = itFinalDataToRecv->second.size();
        for ( int k=0; k<nDataRecv; ++k )
        {
            auto const& idPointsWithNodeRecv = itFinalDataToRecv->second[k].template get<0>();
            auto const& idFacesWithBaryRecv = itFinalDataToRecv->second[k].template get<1>();
            auto const& theelt = this->element( memoryMsgToSend[idProc][k], idProc );

            //update faces data
            for ( size_type j = 0; j < this->numLocalFaces(); j++ )
            {
                auto const& idFaceRecv = idFacesWithBaryRecv[j].template get<0>();
                const bool faceOnBoundaryRecv = idFacesWithBaryRecv[j].template get<1>();
                auto const& baryFaceRecv = idFacesWithBaryRecv[j].template get<2>();

                //objective : find  face_it (hence jBis in theelt ) (permutations would be necessary)
                uint16_type jBis = invalid_uint16_type_value;
                bool hasFind=false;
                for ( uint16_type j2 = 0; j2 < this->numLocalFaces() && !hasFind; j2++ )
                {

                    auto const& thefacej2 = theelt.face( j2 );
                    //auto const& theGj2 = thefacej2.G();
                    auto const& theGj2 = thefacej2.vertices();
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
                    // compare barycenters
                    bool find2=true;
                    for (uint16_type d=0;d<nRealDim;++d)
                        {
                            find2 = find2 && ( std::abs( baryFace[d]-baryFaceRecv[d] )<1e-9 );
                        }
                    if (find2) { hasFind=true;jBis=j2; }

                } //for ( uint16_type j2 = 0; j2 < this->numLocalFaces() && !hasFind; j2++ )

                CHECK ( hasFind ) << "[mesh::updateEntitiesCoDimensionOneGhostCell] : invalid partitioning data, ghost face cells are not available\n";

                // get the good face
                auto face_it = this->faceIterator( theelt.face( jBis ).id() );
                // update id face in other partition
                this->faces().modify( face_it, Feel::detail::updateIdInOthersPartitions( idProc, idFaceRecv ) );
                // maybe the face is not really on boundary
                if ( face_it->isOnBoundary() && !faceOnBoundaryRecv )
                    this->faces().modify( face_it, Feel::detail::UpdateFaceOnBoundary( false ) );

            } // for ( size_type j = 0; j < this->numLocalFaces(); j++ )


            for ( size_type j = 0; j < element_type::numLocalVertices; j++ )
            {
                auto const& idPointRecv = idPointsWithNodeRecv[j].template get<0>();
                auto const& nodePointRecv = idPointsWithNodeRecv[j].template get<1>();

                uint16_type jBis = invalid_uint16_type_value;
                bool hasFind=false;
                for ( uint16_type j2 = 0; j2 < element_type::numLocalVertices && !hasFind; j2++ )
                {
                    auto const& thepointj2 = theelt.point( j2 );
                    // compare barycenters
                    bool find2=true;
                    for (uint16_type d=0;d<nRealDim;++d)
                        {
                            find2 = find2 && ( std::abs( thepointj2(d)-nodePointRecv[d] )<1e-9 );
                        }
                    if (find2) { hasFind=true;jBis=j2; }
                }

                CHECK ( hasFind ) << "[mesh::updateEntitiesCoDimensionOneGhostCell] : invalid partitioning data, ghost point cells are not available\n";
                // get the good face
                auto point_it = this->pointIterator( theelt.point( jBis ).id() );
                //update the face
                this->points().modify( point_it, Feel::detail::updateIdInOthersPartitions( idProc, idPointRecv ) );

            } // for ( size_type j = 0; j < element_type::numLocalVertices; j++ )
        } // for ( int k=0; k<nDataRecv; ++k )
    } // for ( ; itFinalDataToRecv!=enFinalDataToRecv ; ++itFinalDataToRecv)
    //------------------------------------------------------------------------------------------------//

    DVLOG(1) << "updateEntitiesCoDimensionOneGhostCellByUsingNonBlockingComm : finish on rank " << this->worldComm().localRank() << "\n";

}



template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::check() const
{
    if ( nDim != nRealDim )
        return;
#if !defined( NDEBUG )
    VLOG(2) << "[Mesh::check] numLocalFaces = " << this->numLocalFaces() << "\n";
    element_iterator iv = this->beginElementWithProcessId( this->worldComm().localRank() );
    element_iterator en = this->endElementWithProcessId( this->worldComm().localRank() );
    size_type nEltInMesh = std::distance( iv,en );

    //boost::tie( iv, en ) = this->elementsRange();
    for ( ; iv != en; ++iv )
    {
        element_type const& __element = *iv;

        for ( size_type j = 0; j < this->numLocalFaces(); j++ )
        {
            FEELPP_ASSERT( iv->facePtr( j ) )( j )( iv->id() ).error( "invalid element face check" );
            VLOG(2) << "------------------------------------------------------------\n";
            VLOG(2) << "Element : " << iv->id() << " face lid: " << j << " face gid:  " << iv->face( j ).id() << "\n";

        }

        size_type counter = 0;

        for ( uint16_type ms=0; ms < __element.nNeighbors(); ms++ )
        {
            if ( __element.neighbor( ms ).first != invalid_size_type_value )
                ++counter;

        }

        VLOG(2) << "[Mesh::check] element " << __element.id() << " number of neighbors: " << counter << "\n";
        FEELPP_ASSERT( counter >= 1 || nEltInMesh==1 )( __element.id() )( __element.nNeighbors() )( counter ).warn( "invalid neighboring data" );
#if 0

        for ( size_type j = 0; j < ( size_type )element_type::numEdges; ++j )
        {
            FEELPP_ASSERT( iv->edgePtr( j ) )( j )( iv->id() ).error( "invalid element edge check" );
            VLOG(2) << "------------------------------------------------------------\n";
            VLOG(2) << "Element : " << iv->id() << " edge lid: " << j << " edge gid:  " << iv->edge( j ).id() << "\n";

        }

#endif

    }

    // faces check
    typedef typename super::location_face_const_iterator location_face_const_iterator;
    location_face_const_iterator itf = this->beginFaceOnBoundary();
    location_face_const_iterator ite = this->endFaceOnBoundary();
    std::map<int,int> nf;
    for ( ; itf != ite; ++ itf )
    {
        face_type const& __face = *itf;
        DLOG_IF( WARNING, !__face.isConnectedTo0() && !__face.isConnectedTo1() ) << "face not connected to an element face:" << __face << "\n";
        if ( !__face.isConnectedTo0() && !__face.isConnectedTo1() )
        {
            auto it = nf.find( (int)__face.marker().value() );
            if (  it == nf.end() )
                nf[(int)__face.marker().value()] = 1;
            else
                it->second += 1;
        }
        DLOG_IF( WARNING, !__face.isConnectedTo0() && __face.isConnectedTo1() ) << "face  connected to element 1 but not element 0. face:" << __face << "\n";
        DLOG_IF( WARNING, __face.isConnectedTo0() && __face.isConnectedTo1() )<< "invalid boundary face (connected to 2 elements):" << __face << "\n";

        DLOG_IF( WARNING, __face.isConnectedTo0() &&  !__face.element( 0 ).facePtr( __face.pos_first() ) ) << "invalid face in element, face: " << __face << " in element " << __face.element(0) << "\n";

    }
    auto itt = nf.begin();
    auto ent = nf.end();
    for(; itt != ent; ++itt )
    {
        LOG(INFO ) << "face with marker " << itt->first << " not attached " << itt->second << "\n";
    }
#endif
}

#if 0
template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::findNeighboringProcessors()
{
    // Don't need to do anything if there is
    // only one processor.
    if ( this->worldComm().localSize() == 1 )
        return;

#ifdef FEELPP_HAS_MPI

    M_neighboring_processors.clear();

    // Get the bounding sphere for the local processor
    Sphere bounding_sphere = processorBoundingSphere ( *this, this->worldComm().localRank() );

    // Just to be sure, increase its radius by 10%.  Sure would suck to
    // miss a neighboring processor!
    bounding_sphere.setRadius( bounding_sphere.radius()*1.1 );

    // Collect the bounding spheres from all processors, test for intersection
    {
        std::vector<double>
        send ( 4,                         0 ),
             recv ( 4*this->worldComm().localSize(), 0 );

        send[0] = bounding_sphere.center()( 0 );
        send[1] = bounding_sphere.center()( 1 );
        send[2] = bounding_sphere.center()( 2 );
        send[3] = bounding_sphere.radius();

        MPI_Allgather ( &send[0], send.size(), MPI_DOUBLE,
                        &recv[0], send.size(), MPI_DOUBLE,
                        this->worldComm().localComm() );


        for ( unsigned int proc=0; proc<this->worldComm().localSize(); proc++ )
        {
            const Point center ( recv[4*proc+0],
                                 recv[4*proc+1],
                                 recv[4*proc+2] );

            const Real radius = recv[4*proc+3];

            const Sphere proc_sphere ( center, radius );

            if ( bounding_sphere.intersects( proc_sphere ) )
                M_neighboring_processors.push_back( proc );
        }

        // Print out the _neighboring_processors list
        VLOG(2) << "Processor " << this->worldComm().localRank() << " intersects:\n";

        for ( unsigned int p=0; p< M_neighboring_processors.size(); p++ )
            VLOG(2) << " - proc " << M_neighboring_processors[p] << "\n";
    }

#endif
}
#endif

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
        VLOG(2) << "Local numbering in the elements is OK . \n";

    else
    {
        std::for_each( list_of_bad_elts.begin(),
                       list_of_bad_elts.end(),
                       []( size_type const& e )
                       {
                           LOG_FIRST_N(WARNING,10) << "element is not anticlockwise oriented(wrong local numbering): " << e << "\n";
                       });


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
    VLOG(1) << "sending markername\n";
    //this->comm().send( p, tag, M_markername.size() );
    VLOG(1) << "sending markername size: "<< M_markername.size() << "\n";
    BOOST_FOREACH(auto m, M_markername )
    {
        VLOG(1) << "sending key: "<< m.first << "\n";
        //this->comm().send( p, tag, m.first );
        VLOG(1) << "sending value\n";
        //this->comm().send( p, tag, m.second );
    }
}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::recv(int p, int tag)
{
    VLOG(1) << "receiving markername\n";
    //this->comm().recv( p, tag, M_markername );
    int s = 0;
    //this->comm().recv( p, tag, s );
    VLOG(1) << "receiving markername size: "<< s << "\n";
    for( int i = 0; i < s; ++i )
    {
        std::string k;
        VLOG(1) << "receiving key\n";
        //this->comm().recv( p, tag, k );
        VLOG(1) << "receiving key:"<< k << "\n";
        std::vector<int> v;
        VLOG(1) << "receiving value\n";
        //this->comm().recv( p, tag, v );
        VLOG(1) << "receiving value: "<< v[0] << ","<< v[1] <<"\n";
        //M_markername[k]=v;
    }


    //decode();

}

template<typename Shape, typename T, int Tag>
typename Mesh<Shape, T, Tag>::element_iterator
Mesh<Shape, T, Tag>::eraseElement( element_iterator position, bool modify )
{
    for(int i = 0; i < element_type::numTopologicalFaces; ++i )
    {
        auto fit = this->faces().iterator_to( position->face(i) );
        auto modified = this->faces().modify( fit, [position]( face_type& f ) { f.disconnect(*position); } );
        if ( fit->isOnBoundary() && !fit->isConnectedTo0() )
        {
            std::cout << "erase boundary face...\n";
            this->eraseFace( fit );
        }
        if ( !fit->isOnBoundary() && fit->isConnectedTo0() && !fit->isConnectedTo1() )
        {
            std::cout << "found boundary face...\n";
            this->faces().modify( fit, []( face_type& f ){ f.setOnBoundary( true ); } );
        }
    }
    this->elements().modify( position,
                             [] ( element_type& e )
                             {
                                 for ( int i = 0; i < e.numPoints; ++i )
                                     e.point( i ).elements().erase( e.id() );
                             } );
    auto eit = this->elements().erase( position );
    //this->setUpdatedForUse( false );
    //this->updateForUse();
    return eit;
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
    LOG(INFO) <<"decode=   " << this->worldComm().size() << "\n" ;
    LOG(INFO) <<"decode=   " << this->worldComm().subWorldComm().localRank() << "\n";
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
        //pf.setIdInPartition( this->worldComm().localRank(),pf.id() );

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
        //pv.setIdInPartition( this->worldComm().localRank(),pv.id() );

        const int shift = elt_it->second[1]+1;
        for ( uint16_type jj = 0; jj < npoints_per_element; ++jj )
        {
            pv.setPoint( ordering.fromGmshId( jj ), this->point( elt_it->second[shift+jj]-1 ) );
        }
        this->addElement( pv );
#if 0
        __idGmshToFeel=pv.id();
        auto theelt = mesh->elementIterator( pv.id(), pv.partitionId() );
        mesh->elements().modify( theelt, Feel::detail::updateIdInOthersPartitions( this->worldComm().localRank(), pv.id() ) );
#endif
    }
    LOG(INFO) << "distance  elts: "<< std::distance( this->beginElement(), this->endElement() ) << "\n";
    LOG(INFO) << "distance faces: "<< std::distance( this->beginFace(), this->endFace() ) << "\n";
    LOG(INFO) << "distance marker faces: "<< std::distance( this->beginFaceWithMarker(), this->endFaceWithMarker() ) << "\n";
    LOG(INFO) << "distance marker2 faces: "<< std::distance( this->beginFaceWithMarker2(), this->endFaceWithMarker2() ) << "\n";

    //this->components().set ( MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
    //this->updateForUse();
    //std::cout<<"decode=   " << this->worldComm().localSize() << std::endl;
    this->setSubStructuring(true);
}

template<typename Shape, typename T, int Tag1, int Tag2=Tag1, int TheTag=Tag1>
boost::shared_ptr<Mesh<Shape, T, TheTag> >
merge( boost::shared_ptr<Mesh<Shape, T, Tag1> > m1,
       boost::shared_ptr<Mesh<Shape, T, Tag2> > m2 )
{
    boost::shared_ptr<Mesh<Shape, T, TheTag> > m( new Mesh<Shape, T, TheTag> );

    // add the points

    size_type shift_p = std::distance( m1->beginPoint(),  m1->endPoint() );
    for( auto it = m1->beginPoint(),  en = m1->endPoint(); it != en; ++it  )
    {
        m->addPoint( *it );
    }
    for( auto it = m2->beginPoint(),  en = m2->endPoint(); it != en; ++it )
    {
        auto p = *it;
        // shift id
        p.setId( shift_p+p.id() );
        m->addPoint( p );
    }
    typedef typename Mesh<Shape, T, TheTag>::face_iterator face_iterator;
    typedef typename Mesh<Shape, T, TheTag>::element_iterator element_iterator;
    typedef typename Mesh<Shape, T, TheTag>::face_type face_type;
    typedef typename Mesh<Shape, T, TheTag>::element_type element_type;

    auto addface = [&]( boost::shared_ptr<Mesh<Shape, T, TheTag> > m, face_iterator it, size_type shift )
        {
            face_type pf = *it;
            pf.disconnect();
            pf.setOnBoundary( true );
            const uint16_type npoints_per_face = ( face_type::numVertices*face_type::nbPtsPerVertex+
                                                   face_type::numEdges*face_type::nbPtsPerEdge+
                                                   face_type::numFaces*face_type::nbPtsPerFace );
            for ( uint16_type jj = 0; jj < npoints_per_face; ++jj )
            {
                pf.setPoint( jj, m->point( shift+it->point(jj).id() ) );
            }
            m->addFace( pf );
        };
    // add the faces
    size_type shift_f = std::distance( m1->beginFace(),  m1->endFace() );
    for( auto it = m1->beginFace(),  en = m1->endFace(); it != en; ++it )
    {
        addface( m, it, 0 );
    }
    for( auto it = m2->beginFace(),  en = m2->endFace(); it != en; ++it )
    {
        // don't forget to shift the point id in the face
        addface( m, it, shift_p );
    }
    // add the elements
    auto addelement = [&]( boost::shared_ptr<Mesh<Shape, T, TheTag> > m, element_iterator it, size_type shift_f, size_type shift_p )
        {
            element_type pf = *it;

            static const uint16_type npoints_per_element = element_type::numPoints;

#if 0
            for ( uint16_type jj = 0; jj < pf.numLocalFaces; ++jj )
            {
                pf.setFace( jj, m->face( shift_f+it->face(jj).id() ) );
            }
#endif
            for ( uint16_type jj = 0; jj < npoints_per_element; ++jj )
            {
                pf.setPoint( jj, m->point( shift_p+it->point(jj).id() ) );
            }
            m->addElement( pf );
        };
    for( auto it = m1->beginElement(),  en = m1->endElement(); it != en; ++it )
    {
        addelement( m, it, 0, 0 );
    }
    for( auto it = m2->beginElement(),  en = m2->endElement(); it != en; ++it )
    {
        // don't forget to shift the point id in the face
        addelement( m, it, shift_f, shift_p );
    }
    m->components().set ( MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
    m->updateForUse();
    return m;
}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::Inverse::distribute( bool extrapolation )
{
    typename self_type::element_iterator el_it;
    typename self_type::element_iterator el_en;
    boost::tie( boost::tuples::ignore, el_it, el_en ) = Feel::elements( *M_mesh );
    const size_type nActifElt = std::distance( el_it, el_en );
    if ( nActifElt==0 ) return;

    typedef typename self_type::element_type element_type;
    typedef typename gm_type::template Context<vm::JACOBIAN|vm::KB|vm::POINT, element_type> gmc_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
    BoundingBox<> bb;

    typename gm_type::reference_convex_type refelem;
    typename gm_type::precompute_ptrtype __geopc( new typename gm_type::precompute_type( M_mesh->gm(),
            refelem.points() ) );
    boost::unordered_map<size_type,bool> npt;

    M_dist.clear();
    M_ref_coords.clear();
    M_cvx_pts.clear();
    M_pts_cvx.clear();
    M_pts_cvx.clear();

    KDTree::points_type boxpts;
    gmc_ptrtype __c( new gmc_type( M_mesh->gm(),
                                   *el_it,
                                   __geopc ) );
    VLOG(2) << "[Mesh::Inverse] distribute mesh points ion kdtree\n";

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

        DVLOG(2) << "G = " << el_it->G() << " min = " << bb.min << ", max = " << bb.max << "\n";

        // check if the points
        this->pointsInBox( boxpts, bb.min, bb.max );

        DVLOG(2) << "boxpts size = " << boxpts.size() << "\n";


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

                DVLOG(2) << "i = " << i << " index = " << index << " isin = " << ( isin >= -1e-10 )
                         << " xref = " << gic.xRef() << " xreal = " << boost::get<0>( boxpts[i] )
                         << " tobeadded= " << tobeadded << " dist=" << dmin<< "\n";


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

    VLOG(2) << "[Mesh::Inverse] distribute mesh points in kdtree done\n";
}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::Localization::init()
{
    auto mesh = M_mesh.lock();
    if ( !mesh ) return;

    DLOG_IF( WARNING, this->isInit() == false ) << "You have already initialized the tool of localization\n";


    //clear data
    M_geoGlob_Elts.clear();
    M_kd_tree->clear();

    typename self_type::element_iterator el_it;
    typename self_type::element_iterator el_en;
    boost::tie( boost::tuples::ignore, el_it, el_en ) = Feel::elements( *mesh );
    if ( el_it != el_en )
    {
        M_gic.reset( new gmc_inverse_type( mesh->gm(), *el_it, mesh->worldComm().subWorldCommSeq() ) );
        M_gic1.reset( new gmc1_inverse_type( mesh->gm1(), *el_it, mpl::int_<1>(), mesh->worldComm().subWorldCommSeq() ) );
    }

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

    this->computeBarycenter();

    M_isInit=true;
    M_isInitBoundaryFaces=false;

}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::Localization::initBoundaryFaces()
{
    auto mesh = M_mesh.lock();
    if ( !mesh ) return;

    DLOG_IF( WARNING, this->isInitBoundaryFaces() == false ) << "You have already initialized the tool of localization\n";

    //clear data
    M_geoGlob_Elts.clear();
    M_kd_tree->clear()
;
    typename self_type::location_face_iterator face_it;
    typename self_type::location_face_iterator face_en;
    boost::tie( boost::tuples::ignore, face_it, face_en ) = Feel::boundaryfaces( mesh );
    bool hasInitGic=false;
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

                    if ( !hasInitGic )
                    {
                        M_gic.reset( new gmc_inverse_type( mesh->gm(), face_it->element( 0 ), mesh->worldComm().subWorldCommSeq() ) );
                        M_gic1.reset( new gmc1_inverse_type( mesh->gm1(), face_it->element( 0 ), mpl::int_<1>(), mesh->worldComm().subWorldCommSeq() ) );
                        hasInitGic=true;
                    }
                }
        }
    }

    this->computeBarycenter();

    M_isInitBoundaryFaces=true;
    M_isInit=false;

}



template<typename Shape, typename T, int Tag>
boost::tuple<bool,typename Mesh<Shape, T, Tag>::node_type,double>
Mesh<Shape, T, Tag>::Localization::isIn( size_type _id, const node_type & _pt ) const
{
    bool isin=false;
    double dmin;
    node_type x_ref;

    auto mesh = M_mesh.lock();
    //get element with the id
    auto const& elt = mesh->element( _id );

    if ( elt.isOnBoundary() )
        {
#if 0
            // get inverse geometric transformation
            gmc_inverse_type gic( mesh->gm(), elt, mesh->worldComm().subWorldCommSeq() );
            //apply the inverse geometric transformation for the point p
            gic.setXReal( _pt);
            x_ref=gic.xRef();
            // the point is in the reference element ?
            boost::tie( isin, dmin ) = M_refelem.isIn( gic.xRef() );
#else
            M_gic->update( elt );
            M_gic->setXReal( _pt);
            x_ref=M_gic->xRef();
            if ( nDim == nRealDim )
                isin = M_gic->isIn();
            else // in this case, result given with gic->isIn() seems not work (see geomap.hpp)
                boost::tie( isin, dmin ) = M_refelem.isIn( M_gic->xRef() );
#endif
        }
    else
        {
#if 0
            // get inverse geometric transformation
            gmc1_inverse_type gic( mesh->gm1(), elt, mpl::int_<1>(), mesh->worldComm().subWorldCommSeq() );
            //apply the inverse geometric transformation for the point p
            gic.setXReal( _pt);
            x_ref=gic.xRef();
            // the point is in the reference element ?
            boost::tie( isin, dmin ) = M_refelem1.isIn( gic.xRef() );
#else
            M_gic1->update( elt, mpl::int_<1>() );
            M_gic1->setXReal( _pt);
            x_ref=M_gic1->xRef();
            if ( nDim == nRealDim )
                isin = M_gic1->isIn();
            else // in this case, result given with gic->isIn() seems not work (see geomap.hpp)
                boost::tie( isin, dmin ) = M_refelem1.isIn( M_gic1->xRef() );
#endif
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

    auto mesh = M_mesh.lock();
    uint16_type nbIsIn=0;

    for ( uint16_type i = 0; i< nbId ; ++i )
    {
        //get element with the id
        auto const& elt = mesh->element( _ids[i] );

        if ( elt.isOnBoundary() )
            {
#if 0
                // get inverse geometric transformation
                gmc_inverse_type gic( mesh->gm(), elt );
                //apply the inverse geometric transformation for the point p
                gic.setXReal( _pt);
                __x_ref=gic.xRef();
                // the point is in the reference element ?
                boost::tie( isin2, dmin ) = M_refelem.isIn( gic.xRef() );
#else
                M_gic->update( elt );
                M_gic->setXReal( _pt);
                __x_ref=M_gic->xRef();
                if ( nDim == nRealDim )
                    isin2 = M_gic->isIn();
                else // in this case, result given with gic->isIn() seems not work (see geomap.hpp)
                    boost::tie( isin2, dmin ) = M_refelem.isIn( M_gic->xRef() );
#endif
                isin[i] = isin2;
            }
        else
            {
#if 0
                // get inverse geometric transformation
                gmc1_inverse_type gic( mesh->gm1(), elt, mpl::int_<1>() );
                //apply the inverse geometric transformation for the point p
                gic.setXReal( _pt);
                __x_ref=gic.xRef();
                // the point is in the reference element ?
                boost::tie( isin2, dmin ) = M_refelem1.isIn( gic.xRef() );
#else
                M_gic1->update( elt, mpl::int_<1>() );
                M_gic1->setXReal( _pt);
                __x_ref=M_gic1->xRef();
                if ( nDim == nRealDim )
                    isin2 = M_gic1->isIn();
                else // in this case, result given with gic->isIn() seems not work (see geomap.hpp)
                    boost::tie( isin2, dmin ) = M_refelem1.isIn( M_gic1->xRef() );
#endif
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

    DCHECK( this->isInit() || this->isInitBoundaryFaces() ) << "You don't have initialized the tool of localization\n";

    auto mesh = M_mesh.lock();
    bool isin=false;double dmin=0;
    node_type x_ref;
    size_type idEltFound = 0;//mesh->beginElementWithId(mesh->worldComm().localRank())->id();

    std::list< std::pair<size_type, uint> > ListTri;
    this->searchInKdTree(p,ListTri);

    //research the element which contains the point p
    auto itLT=ListTri.begin();
    auto const itLT_end=ListTri.end();

    DCHECK( std::distance( itLT,itLT_end )>0 ) << "problem in localization : listTri is empty\n";

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


    if( !isin && this->doExtrapolation() )
        {
            // first elt
            idEltFound = ListTri.begin()->first;
            auto const& eltUsedForExtrapolation = mesh->element(idEltFound);
            DVLOG(1) << "localisation tool use extrapolation for the point" << p
                     << " with elt.id() " << eltUsedForExtrapolation.id()
                     << " and elt.G() " << eltUsedForExtrapolation.G()
                     << "\n";
            gmc_inverse_type gic( mesh->gm(), eltUsedForExtrapolation, mesh->worldComm().subWorldCommSeq() );
            //apply the inverse geometric transformation for the point p
            gic.setXReal( p);
            boost::tie(isin,idEltFound,x_ref) = boost::make_tuple(true,eltUsedForExtrapolation.id(),gic.xRef());
        }

    return boost::make_tuple( isin, idEltFound, x_ref);

}

template<typename Shape, typename T, int Tag>
boost::tuple<std::vector<bool>, size_type>
Mesh<Shape, T, Tag>::Localization::run_analysis( const matrix_node_type & m,
                                            const size_type & eltHypothetical )
{
    // if no init then init with all geo point
    if ( !( this->isInit() || this->isInitBoundaryFaces() ) )
        this->init();

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
                    DVLOG(1) << "localisation tool (hypthetical elt) : has found pt " << ublas::column( m, i )
                             << "in elt id " << cv_id
                             << "with G() " << M_mesh.lock()->element(cv_id).G()
                             << "\n";
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
                            DVLOG(1)  << "localisation tool (first pass) : has found pt " << ublas::column( m, i )
                                      << " in elt id " << cv_id
                                      << " with G() " << M_mesh.lock()->element(cv_id).G()
                                      << " and ref point " << x_ref
                                      << "\n";
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
                                    DVLOG(1)  << "localisation tool (second pass) : has found pt " << ublas::column( m, i )
                                              << " in elt id " << cv_id
                                              << " with G() " << M_mesh.lock()->element(cv_id).G()
                                              << " and ref point " << x_ref
                                              << "\n";
                                }
                            else if (doExtrapolationAtStart)
                                {
                                    this->setExtrapolation(true);

                                    boost::tie( find_x, cv_id, x_ref ) = this->searchElement(ublas::column( m, i ));

                                    CHECK( find_x ) << "localisation tool : invalid extrapolation \n";

                                    M_resultAnalysis[cv_id].push_back( boost::make_tuple(i,x_ref) );
                                    currentEltHypothetical = cv_id;
                                    hasFindPts[i]=true;

                                    this->setExtrapolation(false);
                                }
                        }
                } // search kdtree

            DLOG_IF(WARNING, !hasFindPts[i]) << "the localisation tool fails to find the point " << ublas::column( m, i ) << "\n";

        } // for (size_type i=0;i< m.size2();++i)

    //revert parameter
    this->setExtrapolation(doExtrapolationAtStart);

    return boost::make_tuple(hasFindPts,cv_id);

} //run_analysis



template<typename Shape, typename T, int Tag>
boost::tuple<bool, std::list<boost::tuple<size_type, typename Mesh<Shape, T, Tag>::node_type> > >
Mesh<Shape, T, Tag>::Localization::searchElements( const node_type & p )
{

    DCHECK( this->isInit() || this->isInitBoundaryFaces() ) << "You don't have initialized the tool of localization\n";

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
        auto mesh = M_mesh->lock();
        //std::cout << "\n WARNING EXTRAPOLATION \n";
        itLT=ListTri.begin();
        elt= mesh->element( itLT->first );
        typename self_type::Inverse::gic_type gic( mesh->gm(), elt );
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
boost::tuple<bool,typename Mesh<Shape, T, Tag>::node_type,double>
Mesh<Shape, T, Tag>::Localization::isIn( size_type _id,
                                         const node_type & _pt,
                                         const matrix_node_type & setPoints,
                                         mpl::int_<1> /**/ ) const
{
    bool isin=true; // warning : start with true
    double dmin=0.;
    node_type x_ref;
    auto mesh = M_mesh.lock();
    //get element with the id
    auto const& elt= mesh->element( _id );
    auto const& eltG = elt.G();

    // check conformity between setPoints (given) and eltG (localize)
    std::vector<bool> find( setPoints.size2(),false );
    for ( size_type i=0; i< setPoints.size2() && isin; ++i )
    {
        auto const& thePt = ublas::column( setPoints,i );

        for ( size_type j=0; j<eltG.size2(); ++j )
        {
            auto ptjeltG = ublas::column( eltG,j );

            if ( ptjeltG.size()==1 )
            {
                if ( std::abs( thePt( 0 )-ptjeltG( 0 ) )<1e-8 )
                    find[i]=true;
            }
            else if ( ptjeltG.size()==2 )
            {
                if ( std::abs( thePt( 0 )-ptjeltG( 0 ) )<1e-8 &&
                     std::abs( thePt( 1 )-ptjeltG( 1 ) )<1e-8 )
                    find[i]=true;
            }
            else if ( ptjeltG.size()==3 )
            {
                if ( std::abs( thePt( 0 )-ptjeltG( 0 ) )<1e-8 &&
                     std::abs( thePt( 1 )-ptjeltG( 1 ) )<1e-8 &&
                     std::abs( thePt( 2 )-ptjeltG( 2 ) )<1e-8 )
                    find[i]=true;
            }
        }
        // up checking
        isin &= find[i];
    }

    // if find -> get ref point and check
    if ( isin )
    {
        bool isin2=false;
        boost::tie(isin2,x_ref,dmin) = this->isIn(_id,_pt);
        LOG_IF(ERROR, !isin2) << "Mesh::Localization::isIn<Conformal> : check fail -> maybe x_ref is not correct";
    }

    return boost::make_tuple(isin,x_ref,dmin);
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
    double dmin=0.;
    node_type x_ref;
    auto mesh=this->mesh().lock();
    size_type idEltFound = mesh->beginElementWithId(mesh->worldComm().localRank())->id();

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
        boost::tie(isin,x_ref,dmin) = this->isIn(itLT->first,p,setPoints,mpl::int_<1>());

        //if not inside, continue the research with an other element
        if ( !isin ) ++itLT;
        else idEltFound=itLT->first;
    } //while ( itLT != itLT_end && !isin  )

    if (!isin)
        {
            if( this->doExtrapolation() )
                {
                    LOG(WARNING) << "WARNING EXTRAPOLATION for the point" << p;
                    //std::cout << "W";
                    auto const& eltUsedForExtrapolation = mesh->element(ListTri.begin()->first);
                    gmc_inverse_type gic( mesh->gm(), eltUsedForExtrapolation, mesh->worldComm().subWorldCommSeq() );
                    //apply the inverse geometric transformation for the point p
                    gic.setXReal( p);
                    boost::tie(isin,idEltFound,x_ref) = boost::make_tuple(true,eltUsedForExtrapolation.id(),gic.xRef());
                }
            else
                {
                    idEltFound = mesh->beginElementWithId(mesh->worldComm().localRank())->id();
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

    DCHECK( this->isInit() || this->isInitBoundaryFaces() ) << "You don't have initialized the tool of localization\n";

    bool find_x=false;
    size_type cv_id=eltHypothetical;
    node_type x_ref;
    double dmin;
    std::vector<bool> hasFindPts(m.size2(),false);

    M_resultAnalysis.clear();
    auto currentEltHypothetical = eltHypothetical;
    for ( size_type i=0; i< m.size2(); ++i )
    {

        bool testHypothetical_find = false;
        if ( eltHypothetical!=invalid_size_type_value )
        {
            boost::tie( testHypothetical_find,x_ref,dmin ) = this->isIn( currentEltHypothetical,ublas::column( m, i ),setPoints,mpl::int_<1>() );
        }
        if ( testHypothetical_find )
        {
            cv_id = currentEltHypothetical;
            M_resultAnalysis[cv_id].push_back( boost::make_tuple(i,x_ref) );
            hasFindPts[i]=true;
        }
        else
        {
            boost::tie( find_x, cv_id, x_ref ) = this->searchElement( ublas::column( m, i ),setPoints,mpl::int_<1>() );

            if ( find_x )
            {
                M_resultAnalysis[cv_id].push_back( boost::make_tuple( i,x_ref ) );
                hasFindPts[i]=true;
                currentEltHypothetical = cv_id;
            }
        }
        //else std::cout<<"\nNew Probleme Localization\n" << std::endl;
    }

    return boost::make_tuple(hasFindPts,cv_id);

} // run_analysis


template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::Localization::computeBarycenter()
{
    M_barycenter = this->computeBarycenter(mpl::int_<nRealDim>());
}

template<typename Shape, typename T, int Tag>
typename Mesh<Shape, T, Tag>::node_type
Mesh<Shape, T, Tag>::Localization::computeBarycenter(mpl::int_<1> /**/) const
{
    node_type res(1);
    res(0)=0;
    if ( this->kdtree()->nPoints()==0 ) return res;
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
Mesh<Shape, T, Tag>::Localization::computeBarycenter(mpl::int_<2> /**/) const
{
    node_type res(2);
    res(0)=0;res(1)=0;
    if ( this->kdtree()->nPoints()==0 )
        return res;
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
Mesh<Shape, T, Tag>::Localization::computeBarycenter(mpl::int_<3> /**/) const
{
    node_type res(3);
    res(0)=0;res(1)=0;res(2)=0;
    if ( this->kdtree()->nPoints()==0 )
        return res;
    for (size_type i = 0 ; i<this->kdtree()->nPoints() ; ++i)
        {
            auto const& pt = this->kdtree()->points()[i].template get<0>();
            res(0)+=pt(0);res(1)+=pt(1);res(2)+=pt(2);
        }
    res(0)/=this->kdtree()->nPoints();res(1)/=this->kdtree()->nPoints();res(2)/=this->kdtree()->nPoints();
    return res;
}


template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::Localization::computeBarycentersWorld()
{
    LOG(INFO) << "computeBarycentersWorld : run mpi::all_gather\n";
    auto mesh = M_mesh.lock();
    M_barycentersWorld = std::vector<boost::tuple<bool,node_type> >( mesh->worldComm().localSize() );
    mpi::all_gather( mesh->worldComm().localComm(),
                     boost::make_tuple( this->kdtree()->nPoints() > 0 , this->barycenter() ),
                     *M_barycentersWorld );
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
