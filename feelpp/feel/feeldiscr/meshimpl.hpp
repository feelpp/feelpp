/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-08-07

  Copyright (C) 2007-2010 Université Joseph Fourier (Grenoble I)
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
   \file meshimpl.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-08-07
 */
#ifndef FEELPP_MESHIMPL_HPP
#define FEELPP_MESHIMPL_HPP 1

#include <boost/preprocessor/comparison/greater_equal.hpp>
#include <boost/preprocessor/list/first_n.hpp>
#include <unordered_map>
#include <unordered_set>

#include <feel/feelconfig.h>
#include <feel/feeltiming/tic.hpp>

#include <feel/feelalg/glas.hpp>
#include <feel/feelalg/solvernonlinearpetsc.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/gmshenums.hpp>
#include <feel/feelmesh/geoentity.hpp>
#include <feel/feelmesh/hypercube.hpp>
#include <feel/feelmesh/simplex.hpp>
#include <feel/feeldiscr/measures.hpp>

namespace Feel
{
#if defined( FEELPP_INSTANTIATION_MODE )

#define DIMS1 BOOST_PP_TUPLE_TO_LIST( 1, ( 1 ) )
#define RDIMS1 BOOST_PP_TUPLE_TO_LIST( 3, ( 1, 2, 3 ) )
#define ORDERS1 BOOST_PP_TUPLE_TO_LIST( 5, ( 1, 2, 3, 4, 5 ) )

#define DIMS2 BOOST_PP_TUPLE_TO_LIST( 1, ( 2 ) )
#define RDIMS2 BOOST_PP_TUPLE_TO_LIST( 2, ( 2, 3 ) )
#define ORDERS2 BOOST_PP_TUPLE_TO_LIST( 5, ( 1, 2, 3, 4, 5 ) )

#define DIMS3 BOOST_PP_TUPLE_TO_LIST( 1, ( 3 ) )
#define RDIMS3 BOOST_PP_TUPLE_TO_LIST( 1, ( 3 ) )
#define ORDERS3 BOOST_PP_TUPLE_TO_LIST( 4, ( 1, 2, 3, 4 ) )

#define FACTORY_SIMPLEX_CODE( LDIM, LORDER, RDIM ) template class Mesh<Simplex<LDIM, LORDER, RDIM>>;
#define FACTORY_HYPERCUBE_CODE( LDIM, LORDER, RDIM ) template class Mesh<Hypercube<LDIM, LORDER, RDIM>>;

#define FACTORY_SIMPLEX_OP( _, GDO ) FACTORY_SIMPLEX_CODE GDO
#define FACTORY_HYPERCUBE_OP( _, GDO ) FACTORY_HYPERCUBE_CODE GDO

#define FACTORY_SIMPLEX_E( LDIM, LORDER, RDIM ) extern template class Mesh<Simplex<LDIM, LORDER, RDIM>>;
#define FACTORY_HYPERCUBE_E( LDIM, LORDER, RDIM ) extern template class Mesh<Hypercube<LDIM, LORDER, RDIM>>;

#define FACTORY_SIMPLEX_OP_E( _, GDO ) FACTORY_SIMPLEX_E GDO
#define FACTORY_HYPERCUBE_OP_E( _, GDO ) FACTORY_HYPERCUBE_E GDO

#if !defined( FEELPP_MESH_IMPL_NOEXTERN_1D )

BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_SIMPLEX_OP_E, 3, ( DIMS1, BOOST_PP_LIST_FIRST_N( FEELPP_MESH_MAX_ORDER, ORDERS1 ), RDIMS1 ) )
BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_HYPERCUBE_OP_E, 3, ( DIMS1, BOOST_PP_LIST_FIRST_N( FEELPP_MESH_MAX_ORDER, ORDERS1 ), RDIMS1 ) )

#endif

#if !defined( FEELPP_MESH_IMPL_NOEXTERN_2D )

BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_SIMPLEX_OP_E, 3, ( DIMS2, BOOST_PP_LIST_FIRST_N( FEELPP_MESH_MAX_ORDER, ORDERS2 ), RDIMS2 ) )
BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_HYPERCUBE_OP_E, 3, ( DIMS2, BOOST_PP_LIST_FIRST_N( FEELPP_MESH_MAX_ORDER, ORDERS2 ), RDIMS2 ) )

#endif

#if !defined( FEELPP_MESH_IMPL_NOEXTERN_3D )

BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_SIMPLEX_OP_E, 3, ( DIMS3, BOOST_PP_LIST_FIRST_N( FEELPP_MESH_MAX_ORDER, ORDERS3 ), RDIMS3 ) )
BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_HYPERCUBE_OP_E, 3, ( DIMS3, BOOST_PP_LIST_FIRST_N( FEELPP_MESH_MAX_ORDER, ORDERS3 ), RDIMS3 ) )

#endif // FEELPP_MESH_IMPL_NOEXTERN_3D
#endif // FEELPP_INSTANTIATION_MODE

namespace meshdetail
{
template <typename MeshType>
typename MeshType::gm1_ptrtype
initGm1( typename MeshType::gm_ptrtype gm, mpl::true_ /**/ )
{
    return gm;
}
template <typename MeshType>
typename MeshType::gm1_ptrtype
initGm1( typename MeshType::gm_ptrtype gm, mpl::false_ /**/ )
{
    return typename MeshType::gm1_ptrtype( new typename MeshType::gm1_type );
}
} // namespace meshdetail

template <typename MeshT>
void print( MeshT* m )
{
    std::cout << "Mesh: npoints: " << m->numPoints() << std::endl;
    std::for_each( m->beginPoint(), m->endPoint(),
                   []( auto const& p ) {
                       std::cout << p.first << " " << p.second.id() << " ";
                       std::for_each( p.second.begin(), p.second.end(),
                                      []( auto const& c ) { std::cout << c << " "; } );
                       std::cout << std::endl;
                   } );
    std::cout << "Mesh nElements: " << m->numElements() << std::endl;
    std::for_each( m->beginElement(), m->endElement(),
                   []( auto const& ep ) {
                       auto const& e = ep.second;
                       std::cout << e.id() << " ";
                       std::for_each( e.beginPoint(), e.endPoint(),
                                      []( auto const& c ) { std::cout << c->id() << " "; } );
                       //{ std::cout << c->id() << " " << c->hasMarker() << " "; } );
                       std::cout << std::endl;
                   } );
}

// Constructor.
template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
 Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::Mesh( std::string const& name, worldcomm_ptr_t const& worldComm, std::string const& props )
    : super( name, worldComm ),
      M_numGlobalElements( 0 ),
      M_gm( new gm_type ),
      M_gm1( Feel::meshdetail::initGm1<type>( M_gm, mpl::bool_<nOrder == 1>() ) ),
      M_meas( 0 ),
      M_measbdy( 0 ),
      M_substructuring( false ),
      M_structure_property( props ),
      //M_part(),
      M_tool_localization( std::make_shared<Localization<self_type>>() )
{
    VLOG( 2 ) << "[Mesh] constructor called\n";
    CHECK( this->hasWorldComm() ) << "Invalid mesh worldComm";
}
template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::partition( const uint16_type n_parts )
{
    //M_part->partition( *this, n_parts );
}

template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::updateForUse()
{
    VLOG( 2 ) << "component                  MESH_RENUMBER: " << this->components().test( MESH_RENUMBER ) << "\n";
    VLOG( 2 ) << "component MESH_UPDATE_ELEMENTS_ADJACENCY: " << this->components().test( MESH_UPDATE_ELEMENTS_ADJACENCY ) << "\n";
    VLOG( 2 ) << "component              MESH_UPDATE_EDGES: " << this->components().test( MESH_UPDATE_EDGES ) << "\n";
    VLOG( 2 ) << "component              MESH_UPDATE_FACES: " << this->components().test( MESH_UPDATE_FACES ) << "\n";
    VLOG( 2 ) << "component                 MESH_PARTITION: " << this->components().test( MESH_PARTITION ) << "\n";
    VLOG( 2 ) << "component        MESH_NO_UPDATE_MEASURES: " << this->components().test( MESH_NO_UPDATE_MEASURES ) << "\n";
    VLOG( 2 ) << "component         MESH_GEOMAP_NOT_CACHED: " << this->components().test( MESH_GEOMAP_NOT_CACHED ) << "\n";

    tic();

    VLOG( 2 ) << "is already updated? : " << this->isUpdatedForUse() << "\n";
    if ( !this->isUpdatedForUse() )
    {
        if ( this->components().test( MESH_RENUMBER ) )
        {
            //this->renumber();
            //VLOG( 2 ) << "[Mesh::updateForUse] renumber : " << ti.elapsed() << "\n";
        }

        element_iterator iv, en;

        if ( this->components().test( MESH_UPDATE_FACES ) || this->components().test( MESH_UPDATE_FACES_MINIMAL ) )
        {
            tic();
            // update connectivities of entities of co dimension 1
            this->updateEntitiesCoDimensionOne();
            // update permutation of entities of co-dimension 1
            this->updateEntitiesCoDimensionOnePermutation();
            if ( this->components().test( MESH_CHECK ) )
                this->check();
            toc( "Mesh::updateForUse update entities of codimension 1", FLAGS_v > 0 );
            VLOG( 1 ) << "[Mesh::updateForUse] update entities of codimension 1";
        }
        else if ( this->components().test( MESH_UPDATE_ELEMENTS_ADJACENCY ) )
        {
            tic();
            this->updateAdjacencyElements();
            toc( "Mesh::updateForUse update adjacency elements", FLAGS_v > 0 );
        }

        if ( this->components().test( MESH_UPDATE_EDGES ) )
        {
            tic();
            // update connectivities of entities of co dimension 2 (edges in 3D)
            this->updateEntitiesCoDimensionTwo();
            toc( "Mesh::updateForUse update entities of codimension 2", FLAGS_v > 0 );
            VLOG( 1 ) << "[Mesh::updateForUse] update entities of codimension 2";
        }

        if ( MeshBase<IndexT>::worldComm().localSize() > 1 )
        {
            tic();
            // update mesh entities with parallel data
            this->updateParallelData();
            toc( "Mesh::updateForUse update parallel data", FLAGS_v > 0 );
        }

        if ( ( this->components().test( MESH_UPDATE_FACES ) || this->components().test( MESH_UPDATE_FACES_MINIMAL ) ) &&
             ( nDim < 3 || this->components().test( MESH_UPDATE_EDGES ) ) )
        {
            tic();
            updateOnBoundary();
            toc( "Mesh::updateForUse update on boundary", FLAGS_v > 0 );
            VLOG( 1 ) << "[Mesh::updateForUse] update on boundary ";
        }

        if ( this->components().test( MESH_ADD_ELEMENTS_INFO ) )
        {
            tic();
            boost::tie( iv, en ) = this->elementsRange();
            for ( ; iv != en; ++iv )
            {
                auto& eltModified = iv->second;
                size_type eltId = eltModified.id();
                for ( uint16_type i = 0; i < eltModified.numPoints; ++i )
                    eltModified.point( i ).addElement( eltId, i );
            }
            toc( "Mesh::updateForUse update add element info", FLAGS_v > 0 );
            VLOG( 1 ) << "[Mesh::updateForUse] update add element info";
        }

        if ( true )
        {
            VLOG( 2 ) << "Adding element information to marked points";
            tic();
            boost::tie( iv, en ) = this->elementsRange();
            for ( ; iv != en; ++iv )
            {
                auto& eltModified = iv->second;
                VLOG( 3 ) << "looking into element " << eltModified.id();

                DCHECK( eltModified.isValid() ) << "Element " << eltModified << " with invalid point";
                // first look if the element has a point with a marker
                // if not we skip this element
                if ( eltModified.hasPointWithMarker() )
                {
                    int index_pt = 0;
                    std::for_each( std::begin( eltModified.points() ), std::end( eltModified.points() ),
                                   [&eltModified, &index_pt]( auto* p ) {
                                       if ( p->hasMarker() )
                                       {
                                           p->addElement( eltModified.id(), index_pt );
                                       }
                                       ++index_pt;
                                   } );
                }
            }
#if 0
            //auto& mpts = this->pointsByMarker();
            //auto& mpts = this->pointsRange();
            for( auto p_it = this->beginPoint(), p_en=this->endPoint(); p_it != p_en; ++p_it )
            {
                auto & pt = p_it->second;
                if ( pt.hasMarker() )
                {
                    if (!pt.elements().size())
                    {
                        pt.setProcessId( invalid_rank_type_value );
                    }
                }
            }
#endif
            toc( "Mesh::updateForUse register elements associated to marked points", FLAGS_v > 0 );
            VLOG( 1 ) << "[Mesh::updateForUse] update add element info for marked points";
        }

        // update number of vertices
        if ( this->numVertices() == 0 )
        {
            if ( nOrder == 1 )
                this->setNumVertices( this->numPoints() );
            else
            {
                tic();
                boost::tie( iv, en ) = this->elementsRange();
                std::unordered_set<size_type> vertexIds;
                vertexIds.reserve( this->numPoints() );
                for ( ; iv != en; ++iv )
                {
                    for ( int p = 0; p < element_type::numVertices; ++p )
                        vertexIds.insert( iv->second.point( p ).id() );
                }
                this->setNumVertices( vertexIds.size() );
                vertexIds.clear();
                std::unordered_set<size_type>().swap( vertexIds );
                toc( "Mesh::updateForUse update number of vertices", FLAGS_v > 0 );
            }
        }

        this->updateMeshFragmentation();

        this->updateNumGlobalElements<mesh_type>();

        if ( this->components().test( MESH_PROPAGATE_MARKERS ) )
        {
            tic();
            propagateMarkers( mpl::int_<nDim>() );
            toc( "Mesh::updateForUse update propagate markers", FLAGS_v > 0 );
            VLOG( 1 ) << "[Mesh::updateForUse] update propagate markers";
        }

        tic();
        this->updateCommonDataInEntities();
        toc( "Mesh::updateForUse update setMesh in elements and faces", FLAGS_v > 0 );
        VLOG( 1 ) << "[Mesh::updateForUse] update setMesh in elements and faces";

        if ( !this->components().test( MESH_NO_UPDATE_MEASURES ) )
        {
            tic();
            // update meas, measface, hAverage, hMin, hMax, measure of the mesh and measure of the boundary mesh
            this->updateMeasures();
            toc( "Mesh::updateForUse update mesh measures", FLAGS_v > 0 );
            VLOG( 1 ) << "[Mesh::updateForUse] update mesh measures";
        }

        // attach mesh in the localisation tool
        M_tool_localization->setMesh( this->shared_from_this(), false );

        if ( this->components().test( MESH_CHECK ) )
        {
            tic();
            // check mesh connectivity
            this->check();
            toc( "Mesh::updateForUse check", FLAGS_v > 0 );
            VLOG( 1 ) << "[Mesh::updateForUse] check";
        }

        this->applyUpdateInformationObject();

    } // isUpdatedForUse

    if ( !this->components().test( MESH_GEOMAP_NOT_CACHED ) && !M_is_gm_cached )
    {
        tic();
        M_gm->initCache( this );
        if ( nOrder > 1 )
            M_gm1->initCache( this );
        M_is_gm_cached = true;
        toc( "Mesh::updateForUse update geomap cache", FLAGS_v > 0 );
        VLOG( 1 ) << "[Mesh::updateForUse] update geomap cache";
    }

    this->setUpdatedForUse( true );

    toc( "[Mesh::updateForUse] total time", FLAGS_v > 0 );
    if ( Environment::isMasterRank() && FLAGS_v >= 1 )
    {
        auto mem = Environment::logMemoryUsage( "memory usage after update for use" );
        std::cout << "[Mesh::updateForUse] resident memory : " << mem.memory_usage / 1.e9 << "GBytes\n";
    }
}

template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::updateForUseAfterMovingNodes( bool upMeasures )
{
    // reset geomap cache
    if ( this->gm()->isCached() )
    {
        this->gm()->initCache( this/*imesh.get()*/ );
        if ( mesh_type::nOrder > 1 )
            this->gm1()->initCache( this/*imesh.get()*/ );
    }

    // update measures
    if ( upMeasures )
        this->updateMeasures();

    // reset localisation tool
    this->tool_localization()->reset();

    for ( MeshSupportBase* ms : this->meshSupportsAttached() )
        ms->resetLocalizationTool();


#if !defined( __INTEL_COMPILER )
    // notify observers that the mesh has changed
    LOG(INFO) << "Notify observers that the mesh has changed\n";
    this->meshChanged( MESH_CHANGES_POINTS_COORDINATES );
#endif
}

template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::updateMeasures()
{
    if ( this->components().test( MESH_NO_UPDATE_MEASURES ) )
        return;

    tic();

    M_local_meas = 0;
    M_local_measbdy = 0;

    typedef typename element_type::quad_meas_type quad_meas_type;
    typedef typename element_type::quad_meas1_type quad_meas1_type;
    quad_meas_type thequad = element_type::imMeasure();
    quad_meas1_type thequad1 = element_type::imMeasureOrder1();
    typename gm_type::precompute_ptrtype pc;
    typename gm_type::faces_precompute_type pcf;
    typename gm1_type::precompute_ptrtype pc1;
    typename gm1_type::faces_precompute_type pcf1;
    // assume that a high order mesh is straightened (probably we need to store this info in the mesh)
    bool meshIsStraightened = ( nOrder > 1 );
    bool updateMeasureWithPc = true;
    if ( updateMeasureWithPc )
    {
        pc = M_gm->preCompute( M_gm, thequad.points() );
        if ( nDim == nRealDim ) // else not works
            pcf = M_gm->preComputeOnFaces( M_gm, thequad.allfpoints() );
        if ( meshIsStraightened )
        {
            pc1 = M_gm1->preCompute( M_gm1, thequad1.points() );
            if ( nDim == nRealDim )
                pcf1 = M_gm1->preComputeOnFaces( M_gm1, thequad1.allfpoints() );
        }
    }

    element_iterator iv, en;
    boost::tie( iv, en ) = this->elementsRange();

    if ( iv != en )
    {
        auto const& eltInit = iv->second;
        std::shared_ptr<typename gm_type::template Context<element_type>> ctx;
        std::shared_ptr<typename gm_type::template Context<element_type,1>> ctxf;
        std::shared_ptr<typename gm1_type::template Context<element_type>> ctx1;
        std::shared_ptr<typename gm1_type::template Context<element_type,1>> ctxf1;
        if ( pc )
            ctx = M_gm->template context<vm::JACOBIAN>( eltInit, pc );
        if ( !pcf.empty() )
            ctxf = M_gm->template context</*vm::POINT|*/ vm::NORMAL | vm::KB | vm::JACOBIAN>( eltInit, pcf, 0 );
        if ( pc1 )
            ctx1 = M_gm1->template context<vm::JACOBIAN>( eltInit, pc1 );
        if ( !pcf1.empty() )
            ctxf1 = M_gm1->template context</*vm::POINT|*/ vm::NORMAL | vm::KB | vm::JACOBIAN>( eltInit, pcf1, 0 );

        for ( ; iv != en; ++iv )
        {
            auto& elt = iv->second;

            if ( meshIsStraightened && !elt.isOnBoundary() )
            {
                ctx1->template update<vm::JACOBIAN>( elt );
                elt.updateWithCtx( thequad1, ctx1, ctxf1 );
            }
            else
            {
                ctx->template update<vm::JACOBIAN>( elt );
                elt.updateWithCtx( thequad, ctx, ctxf );
            }
            // only compute meas for active element (no ghost)
            if ( !elt.isGhostCell() )
                M_local_meas += elt.measure();

            if constexpr ( nDim == 1 || nDim != nRealDim )
                M_local_measbdy = 0;
            else
                for ( int f = 0; f < elt.numTopologicalFaces; ++f )
                    M_local_measbdy += elt.faceMeasure( f );
        } // for element loop
    }     // if ()
#if BOOST_VERSION >= 105500
    std::vector<value_type> gmeas{M_local_meas, M_local_measbdy};
    mpi::all_reduce( MeshBase<IndexT>::worldComm(), mpi::inplace( gmeas.data() ), 2, std::plus<value_type>() );
#else
    std::vector<value_type> lmeas{M_local_meas, M_local_measbdy};
    std::vector<value_type> gmeas( 2, 0.0 );
    mpi::all_reduce( MeshBase::worldComm(), lmeas.data(), 2, gmeas.data(), std::plus<value_type>() );
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
            auto& eltModified = iv->second;
            value_type meas = 0;
            for ( auto const& _elt : eltModified.pointElementNeighborIds() )
            {
                // warning : only compute meas for active element (no ghost)
                if ( this->hasElement( _elt ) )
                    meas += this->element( _elt ).measure();
            }
            eltModified.setMeasurePointElementNeighbors( meas );
        }
    }
    toc( "[Mesh::updateMeasures] update entity measures", FLAGS_v > 0 );

    // compute h information: average, min and max
    {
        tic();
        auto [ havg, hmin, hmax ] = hMeasures( *this );
        M_h_avg = havg;
        M_h_min = hmin;
        M_h_max = hmax;

        LOG( INFO ) << "h average : " << this->hAverage() << "\n";
        LOG( INFO ) << "    h min : " << this->hMin() << "\n";
        LOG( INFO ) << "    h max : " << this->hMax() << "\n";
        toc( "[Mesh::updateMeasures] update average, min and max", FLAGS_v > 0 );
    }
}

template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
typename  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::self_type&
 Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::operator+=( self_type const& m )
{
    std::map<std::vector<double>, size_type> mapDel;
    for ( auto it = m.beginPoint(), en = m.endPoint();
          it != en;
          ++it )
    {
        bool found = false;
        if ( it->second.isOnBoundary() )
        {
            for ( auto pit = this->beginPoint(), pen = this->endPoint();
                  pit != pen;
                  ++pit )
            {
                auto const& thepoint = pit->second;
                if ( thepoint.isOnBoundary() )
                {
                    if ( ublas::norm_2( thepoint.node() - pit->second.node() ) < 1e-10 )
                    {
                        found = true;
                        std::vector<double> ublasCopy( nDim );
                        std::copy( thepoint.node().begin(), thepoint.node().end(), ublasCopy.begin() );
                        // std::cout << "Attention ! Le point " << it->id() << " doit etre remplace par le point " << pit->id() << " (" << it->node() << " ~= " << pit->node() << ")" << std::endl;
                        mapDel.insert( std::pair<std::vector<double>, size_type>( ublasCopy, pit->second.id() ) );
                        break;
                    }
                }
            }
        }
        if ( !found )
            this->addPoint( it->second );
    }
    for ( auto it = m.beginFace(), en = m.endFace(); it != en; ++it )
    {
        // eventually need work
        if ( !it->second.isOnBoundary() )
        {
            face_type f = it->second;
            f.disconnect();
            this->addFace( f );
        }
    }
    for ( auto it = m.beginElement(), en = m.endElement(); it != en; ++it )
    {
        element_type e = it->second;
        for ( uint16_type p = 0; p < e.numPoints; ++p )
        {
            if ( e.point( p ).isOnBoundary() )
            {
                std::vector<double> ublasCopy( nDim );
                std::copy( e.point( p ).node().begin(), e.point( p ).node().end(), ublasCopy.begin() );
                auto itMap = mapDel.find( ublasCopy );
                if ( itMap != mapDel.end() )
                {
                    //    std::cout << "Attention ! Le point " << it->point(p).id() << " va etre remplace par le point " << this->point(itMap->second).id() << std::endl;
                    e.setPoint( p, this->point( itMap->second ) );
                }
            }
        }
        this->addElement( e );
    }
    this->setUpdatedForUse( false );
    this->updateForUse();
    return *this;
}
template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void
     Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::propagateMarkers( mpl::int_<2> )
{
    // propagate top-down marker from face if edge has not been marked
    for ( auto itf = this->beginFace(), enf = this->endFace(); itf != enf; ++itf )
    {
        auto const& f = itf->second;
        if ( !f.hasMarker() )
            return;

        // update points
        for ( int i = 0; i < face_type::numPoints; ++i )
        {
            if ( !f.point( i ).hasMarker() )
            {
                // inherit marker from face
                this->pointIterator( f.point( i ).id() )->second.setMarker( f.marker().value() );
            }
        }
    }
}

template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void
     Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::propagateMarkers( mpl::int_<3> )
{
#if 0
    // first propagate top-down  marker from edges if points have not been marked
    for ( auto ite = this->beginEdge(), ene= this->endEdge(); ite !=ene; ++ite )
    {
        auto const& e = ite->second;
        if ( !e.hasMarker() )
            return;

        for( int i = 0; i < edge_type::numPoints; ++i )
        {
            if ( !e.point( i ).hasMarker() )
            {
                // inherit marker from edge
                this->pointIterator( e.point( i ).id() )->second.setMarker( e.marker().value() );
            }
        }
    }

    // then propagate top-down marker from face if edge has not been marked
    for ( auto itf = this->beginFace(), enf= this->endFace(); itf !=enf; ++itf )
    {
        auto const& f = itf->second;
        if ( !f.hasMarker() )
            return;

        // update points
        for( int i = 0; i < face_type::numPoints; ++i )
        {
            if ( !f.point( i ).hasMarker() )
            {
                // inherit marker from edge
                this->pointIterator( f.point( i ).id() )->second.setMarker( f.marker().value() );
            }
        }
        // update edges
        for( int i = 0; i < face_type::numEdges; ++i )
        {
            if ( !f.edge( i ).hasMarker() )
            {
                // inherit marker from edge
                this->edgeIterator( f.edge(i).id() )->second.setMarker( f.marker().value() );
            }
        }
    }
#else
    CHECK( 0 ) << "propagateMarkers not yet implemented";
#endif
}

template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
template <typename TheShape>
void
     Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::updateCommonDataInEntities( std::enable_if_t<TheShape::nDim == 0>* )
{
    for ( auto itp = this->beginPoint(), enp = this->endPoint(); itp != enp; ++itp )
        itp->second.setMesh( this );
}
template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
template <typename TheShape>
void
     Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::updateCommonDataInEntities( std::enable_if_t<TheShape::nDim == 1>* )
{
    //M_geondEltCommon = std::make_shared<GeoNDCommon<typename element_type::super>>( this, this->gm(), this->gm1() );
    //for ( auto iv = this->beginElement(), en = this->endElement(); iv != en; ++iv )
    //iv->second.setCommonData( M_geondEltCommon.get() );
    for ( auto itf = this->beginFace(), ite = this->endFace(); itf != ite; ++itf )
        itf->second.setMesh( this );
    for ( auto itp = this->beginPoint(), enp = this->endPoint(); itp != enp; ++itp )
        itp->second.setMesh( this );
}
template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
template <typename TheShape>
void
     Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::updateCommonDataInEntities( std::enable_if_t<TheShape::nDim == 2>* )
{
    //M_geondEltCommon = std::make_shared<GeoNDCommon<typename element_type::super>>( this, this->gm(), this->gm1() );
    M_geondFaceCommon = std::make_shared<GeoNDCommon<typename face_type::super>>( this /*,this->gm(), this->gm1()*/ );
    //for ( auto iv = this->beginElement(), en = this->endElement(); iv != en; ++iv )
    // iv->second.setCommonData( M_geondEltCommon.get() );
    for ( auto itf = this->beginFace(), ite = this->endFace(); itf != ite; ++itf )
        itf->second.setCommonData( M_geondFaceCommon.get() );
    for ( auto itp = this->beginPoint(), enp = this->endPoint(); itp != enp; ++itp )
        itp->second.setMesh( this );
}
template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
template <typename TheShape>
void
     Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::updateCommonDataInEntities( std::enable_if_t<TheShape::nDim == 3>* )
{
    //M_geondEltCommon = std::make_shared<GeoNDCommon<typename element_type::super>>( this, this->gm(), this->gm1() );
    M_geondFaceCommon = std::make_shared<GeoNDCommon<typename face_type::super>>( this /*,this->gm(), this->gm1()*/ );
    M_geondEdgeCommon = std::make_shared<GeoNDCommon<typename edge_type::super>>( this /*,this->gm(), this->gm1()*/ );
    //for ( auto iv = this->beginElement(), en = this->endElement(); iv != en; ++iv )
    //    iv->second.setCommonData( M_geondEltCommon.get() );
    for ( auto itf = this->beginFace(), ite = this->endFace(); itf != ite; ++itf )
        itf->second.setCommonData( M_geondFaceCommon.get() );
    for ( auto ite = this->beginEdge(), ene = this->endEdge(); ite != ene; ++ite )
        ite->second.setCommonData( M_geondEdgeCommon.get() );
    for ( auto itp = this->beginPoint(), enp = this->endPoint(); itp != enp; ++itp )
        itp->second.setMesh( this );
}

template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void
     Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::renumber( mpl::bool_<true> )
{
#if 0
    size_type next_free_node = 0;

    // map old/new ids
    std::vector<size_type> node_map( this->numPoints(), invalid_v<size_type> );
    //std::map<size_type> node_map( this->numPoints(), invalid_v<size_type> );
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
#if !defined( NDEBUG )
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
            if ( node_map[ __true_id ] != invalid_v<size_type> )
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
#if !defined( NDEBUG )
    std::vector<size_type> check_id( node_map );
    std::unique( check_id.begin(), check_id.end() );

    FEELPP_ASSERT( check_id.size() == node_map.size() )( node_map.size() )( check_id.size() ).error( "all ids must be unique" );

    // in parallel this will generate a lot of warning since all the points are
    //loaded into the mesh data structure including the ones not belonging to
    //the current processor. Do the test only in sequential
    if ( Environment::numberOfProcessors() == 1 )
    {
        FEELPP_ASSERT( std::find( node_map.begin(),node_map.end(), invalid_v<size_type> ) == node_map.end() ).warn("invalid size_type value found as id " );
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
        DVLOG(2) << "face id: " << __face.id() << "\n";

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

#if !defined( NDEBUG )
        checkLocalPermutation( mpl::bool_< ( Shape == SHAPE_TETRA ) >() );
#endif
    }

    // should we renumber also the faces and elements ?

#endif
}

template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::renumber( std::vector<size_type> const& node_map, mpl::int_<1> )
{
}
template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::renumber( std::vector<size_type> const& node_map, mpl::int_<2> )
{
}

template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::renumber( std::vector<size_type> const& node_map, mpl::int_<3> )
{
#if 0
    for ( auto elt = this->beginEdge();
            elt != this->endEdge(); ++elt )
    {
        edge_type __edge = elt->second;

        DVLOG(2) << "edge id: " << __edge.id() << "\n";

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
#else
    CHECK( 0 ) << "renumber not implemented yet";
#endif
}
template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::localrenumber()
{
#if 0
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
#else
    CHECK( 0 ) << "localrenumber not implemented yet";
#endif
} /** void LocalRenumber **/

template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::updateEntitiesCoDimensionOne()
{
    bool updateComponentFacesMinimal = !this->components().test( MESH_UPDATE_FACES ) && this->components().test( MESH_UPDATE_FACES_MINIMAL );
    if ( updateComponentFacesMinimal )
        this->updateEntitiesCoDimensionOneMinimal();
    else
    {
        //updateEntitiesCoDimensionOne( mpl::bool_<nDim==nRealDim>() );
        this->updateEntitiesCoDimensionOne( mpl::bool_<true>() );
    }
}
template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void
     Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::updateEntitiesCoDimensionOne( mpl::bool_<false> )
{
}

template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void
     Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::updateEntitiesCoDimensionOne( mpl::bool_<true> )
{
    const bool updateComponentAddElements = this->components().test( MESH_ADD_ELEMENTS_INFO );

    rank_type currentPid = MeshBase<IndexT>::worldComm().localRank();
    rank_type numPartition = MeshBase<IndexT>::worldComm().localSize();

    // stores local vertex ids of the faces to identify them
    std::vector<size_type> lids( face_type::numVertices );

    size_type next_face = 0;
    bool faceinserted = false;
    const uint16_type _numLocalFaces = this->numLocalFaces();
    // element_iterator iv,  en;
    // boost::tie( iv, en ) = this->elementsRange();
    auto iv = this->beginOrderedElement();
    auto en = this->endOrderedElement();

    size_type nElt = std::distance( iv, en );
    typedef std::unordered_map<std::vector /*set*/<size_type>, std::tuple<element_type*, uint16_type, face_type*>, Feel::HashTables::HasherContainers<size_type>> pointstoface_container_type;
    pointstoface_container_type _faces( nElt * _numLocalFaces );
    typename pointstoface_container_type::iterator _faceit;

    tic();
    // First We check if we have already Faces stored
    if ( !this->faces().empty() )
    {
        face_iterator __it = this->beginFace();
        face_iterator __en = this->endFace();
        LOG( INFO ) << "We have " << std::distance( __it, __en ) << " faces in the database";
        for ( ; __it != __en; )
        {
            auto& face = __it->second;
            for ( int f = 0; f < face_type::numVertices; ++f )
            {
                lids[f] = face.point( f ).id();
            }
            std::sort( lids.begin(), lids.end() );

            boost::tie( _faceit, faceinserted ) = _faces.emplace( std::piecewise_construct,
                                                                  std::forward_as_tuple( lids ),
                                                                  std::forward_as_tuple( nullptr, invalid_uint16_type_value, &face ) );

            DVLOG_IF( 2, faceinserted ) << "added face with id " << face.id() << "\n";
            DVLOG_IF( 2, !faceinserted ) << "not added face with id " << face.id()
                                         << " was already face with id = " << std::get<2>( _faceit->second )->id() << "\n";

            if ( faceinserted )
            {
                next_face = std::max( next_face, face.id() + 1 );
                ++__it;
            }
            else
            {
                // if ( __it->hasMarker() && ( !facePtr->hasMarker() || ( __it->marker() != facePtr->marker() ) ) )
                //     facePtr->setMarker( __it->marker().value() );
                __it = this->eraseFace( __it );
            }
        }
    }
    toc( "Mesh.updateEntitiesCoDimensionOne.add_faces", FLAGS_v > 1 );

    bool inserted = false;
    face_iterator __fit;

    std::vector<uint16_type> myfToP( face_type::numVertices * _numLocalFaces );
    for ( uint16_type j = 0; j < _numLocalFaces; j++ )
    {
        for ( int f = 0; f < face_type::numVertices; ++f )
            myfToP[j * face_type::numVertices + f] = ( nDim == 1 ) ? j : /*iv->*/ element_type::fToP( j, f );
    }

    // std::vector< std::set<size_type> > pointsFromFace( this->numLocalFaces() );
    // std::array<size_type,element_type::numVertices> pointIdInElt;
    Eigen::Matrix<size_type, element_type::numVertices, 1> pointIdInElt;
    // for ( uint16_type f = 0; f < element_type::numVertices; ++f )
    //     pointIdInElt[f] = iv->point( f ).id();

    // std::set<size_type> s;
    tic();

    // get an ordering by process (need in // for have a deterministic choice of face connection (0 and 1) and therfore to define the permutations )
    std::map<rank_type, std::vector<element_type*>, std::less<rank_type>> elementsByProcOrdering;
    for ( ; iv != en; ++iv )
    {
        auto& elt = unwrap_ref( *iv );
        rank_type pid = elt.processId();
        elementsByProcOrdering[pid].push_back( &elt );
    }
    // iterator over all element
    for ( auto const& elementsByProcPair : elementsByProcOrdering )
    {
        for ( element_type* eltPtr : elementsByProcPair.second )
        {
            element_type& elt = *eltPtr;
            const size_type eltId = elt.id();
            const rank_type eltPid = elt.processId();

            for ( uint16_type f = 0; f < element_type::numVertices; ++f )
                pointIdInElt[f] = elt.point( f ).id();

            for ( uint16_type j = 0; j < _numLocalFaces; j++ )
            {
#if !defined( NDEBUG )
                DVLOG( 2 ) << "------------------------------------------------------------\n";
                DVLOG( 2 ) << "Element id: " << eltId << " local face id: " << j << "\n";
#endif
                for ( uint16_type f = 0; f < face_type::numVertices; ++f )
                {
                    // uint16_type pt_localid = ( nDim==1 )?j:/*__element.*/element_type::fToP( j, f );
                    // lids[f]= __element.point( pt_localid ).id();
                    lids[f] = pointIdInElt[myfToP[j * face_type::numVertices + f]];

#if !defined( NDEBUG )
                    VLOG( 3 ) << "add point local id " << f << " to face " << j << " " << elt.fToP( j, f )
                              << " global id " << myfToP[j * face_type::numVertices + f] /*__element.point( pt_localid ).id()*/ << "\n";
#endif
                }
                std::sort( lids.begin(), lids.end() );

                //boost::tie( _faceit, faceinserted ) = _faces.insert( std::make_pair( lids, nullptr/*next_face*/ ) );
                boost::tie( _faceit, faceinserted ) = _faces.emplace( std::piecewise_construct,
                                                                      std::forward_as_tuple( lids ),
                                                                      std::forward_as_tuple( eltPtr, j, nullptr ) );

                if ( faceinserted )
                {
                    DVLOG( 2 ) << "creating a new face : " << next_face << "\n";

                    face_type newFace;
                    newFace.setId( next_face++ );
                    newFace.setProcessIdInPartition( currentPid );
                    newFace.setProcessId( eltPid );
                    newFace.setOnBoundary( true, face_type::nDim );
                    newFace.setConnection0( boost::make_tuple( &elt, j ) );
                    for ( size_type k = 0; k < face_type::numPoints; ++k )
                        newFace.setPoint( k, elt.point( /*ele.*/ element_type::fToP( j, k ) ) );
                    if ( updateComponentAddElements )
                        newFace.addElement( eltId );
                    auto res = this->addFace( newFace );
                    auto& faceInserted = res.first->second;
                    elt.setFace( j, faceInserted );
                    std::get<2>( _faceit->second ) = &faceInserted;

#if !defined( NDEBUG )
                    DVLOG( 2 ) << "Adding [new] face info : \n";
                    DVLOG( 2 ) << "element id: " << eltId << "\n";
                    DVLOG( 2 ) << "process id: " << faceInserted.processId() << "\n";
                    DVLOG( 2 ) << "id: " << faceInserted.id() << "\n";
                    DVLOG( 2 ) << "bdy: " << faceInserted.isOnBoundary() << "\n";
                    DVLOG( 2 ) << "ad_first: " << faceInserted.ad_first() << "\n";
                    DVLOG( 2 ) << "pos_first: " << faceInserted.pos_first() << "\n";
                    DVLOG( 2 ) << "proc_first: " << faceInserted.proc_first() << "\n";
#endif
                }
                else
                {
                    DVLOG( 2 ) << "found the face in element " << eltId << " and local face: " << j << "\n";

#if 0
                    face_type* facePtr =  _faceit->second;
                    face_type const& cface = *facePtr;
                    if ( updateComponentAddElements )
                        facePtr->addElement( __element_id );

                    // the three conditions below typically arise after reading a serialized mesh
                    if ( cface.isConnectedTo0() && cface.connection0().template get<0>() == 0 && ( __element_id == cface.ad_first() ) )
                    {
                        DVLOG(2) << "fixing connection 0 in face\n";
                        // reconnect the elements to the face
                        auto connect0 = facePtr->connection0();
                        connect0.template get<0>() = boost::addressof( __element );
                        facePtr->setConnection0( connect0 );
                    }
                    if ( cface.isConnectedTo1() && cface.connection1().template get<0>() == 0 && ( __element_id == cface.ad_second() ) )
                    {
                        DVLOG(2) << "fixing connection 1 in face\n";
                        // reconnect the elements to the face
                        auto connect1 = cface.connection1();
                        connect1.template get<0>() = boost::addressof( __element );
                        facePtr->setConnection1( connect1 );
                    }
                    if ( cface.isConnectedTo0() && cface.isConnectedTo1() )
                    {
                        DVLOG(2) << "internal face, fixing process id if necessary\n";
                        FEELPP_ASSERT( __element.facePtr( j ) )( j )( __element.id() ).warn( "invalid element face error" );
                        FEELPP_ASSERT( cface.isConnectedTo0() && cface.isConnectedTo1() )
                            ( cface.isConnectedTo0() )( cface.isConnectedTo1() ).error ("inconsistent data structure" );
                        if ( cface.processId()!=MeshBase<IndexT>::worldComm().localRank() )
                        {
                            if ( ( cface.element0().processId()==MeshBase<IndexT>::worldComm().localRank() ) ||
                                 ( cface.element1().processId()==MeshBase<IndexT>::worldComm().localRank() ) )
                                facePtr->setProcessId( currentPid );
                        }
                    }
#endif

                    auto& f2eVal = _faceit->second;
                    element_type* elt0 = std::get<0>( f2eVal );
                    face_type* facePtr = std::get<2>( f2eVal );
                    if ( !elt0 )
                    {
                        // the face could have been entered apriori given by
                        // the mesh generator, so just set the connection0
                        // properly .
                        CHECK( facePtr != nullptr ) << "invalid face";
                        facePtr->setConnection0( boost::make_tuple( eltPtr, j ) );
                        facePtr->setProcessId( eltPid );
                        facePtr->setOnBoundary( true, face_type::nDim );
                        if ( updateComponentAddElements )
                            facePtr->addElement( eltId );
                        std::get<0>( _faceit->second ) = &elt;
                        std::get<1>( _faceit->second ) = j;
                        elt.setFace( j, *facePtr );
#if !defined( NDEBUG )
                        DVLOG( 2 ) << "adding [!isConnectedTo0] face info : \n";
                        DVLOG( 2 ) << "id: " << facePtr->id() << "\n";
                        DVLOG( 2 ) << "process id: " << facePtr->processId() << "\n";
                        DVLOG( 2 ) << "bdy: " << facePtr->isOnBoundary() << "\n";
                        if ( facePtr->hasMarker() )
                            DVLOG( 2 ) << "marker: " << facePtr->marker() << "\n";
                        DVLOG( 2 ) << "ad_first: " << facePtr->ad_first() << "\n";
                        DVLOG( 2 ) << "pos_first: " << facePtr->pos_first() << "\n";
                        DVLOG( 2 ) << "proc_first: " << facePtr->proc_first() << "\n";
#endif
                    }
                    else
                    {
                        // found an internal faces
                        elt0->setNeighbor( std::get<1>( f2eVal ), eltId );
                        elt.setNeighbor( j, elt0->id() );

                        CHECK( facePtr && facePtr->isConnectedTo0() ) << "invalid face";
                        if ( !facePtr->isConnectedTo1() ) // add connection1 to the face
                        {
                            facePtr->setConnection1( boost::make_tuple( eltPtr, j ) );
                            facePtr->setOnBoundary( false );
                            elt.setFace( j, *facePtr );
                            if ( updateComponentAddElements )
                                facePtr->addElement( eltId );
                            // fix duplication of point in connection1 with 3d mesh at order 3 and 4
                            this->fixPointDuplicationInHOMesh( elt, *facePtr, mpl::bool_ < nDim == 3 && nOrder >= 3 > () );

                            // force processId equal to M_comm.rank() if face on interprocessfaces
                            if ( facePtr->processId() != currentPid )
                            {
                                if ( ( facePtr->element0().processId() == currentPid ) || ( facePtr->element1().processId() == currentPid ) )
                                    facePtr->setProcessId( currentPid );
                            }

#if !defined( NDEBUG )
                            DVLOG( 2 ) << "adding face info : \n";
                            DVLOG( 2 ) << "id: " << facePtr->id() << "\n";
                            DVLOG( 2 ) << "process id: " << facePtr->processId() << "\n";
                            DVLOG( 2 ) << "bdy: " << facePtr->isOnBoundary() << "\n";
                            if ( facePtr->hasMarker() )
                                DVLOG( 2 ) << "marker: " << facePtr->marker() << "\n";
                            DVLOG( 2 ) << "ad_first: " << facePtr->ad_first() << "\n";
                            DVLOG( 2 ) << "pos_first: " << facePtr->pos_first() << "\n";
                            DVLOG( 2 ) << "proc_first: " << facePtr->proc_first() << "\n";
                            DVLOG( 2 ) << "ad_second: " << facePtr->ad_second() << "\n";
                            DVLOG( 2 ) << "pos_second: " << facePtr->pos_second() << "\n";
                            DVLOG( 2 ) << "proc_second: " << facePtr->proc_second() << "\n";
#endif
                        }

#if !defined( NDEBUG )
                        CHECK( ( facePtr->processId() == facePtr->proc_first() ) ||
                               ( facePtr->processId() == facePtr->proc_second() ) )
                            << "invalid process id " << facePtr->processId() << " with element proc first = " << facePtr->proc_first()
                            << " and element proc second " << facePtr->proc_second();
#endif
                    }
                }
            } // face loop
        }     // element loop
    }
    toc( "Mesh.updateEntitiesCoDimensionOne.add_faces_from_elements", FLAGS_v > 1 );
    DVLOG( 2 ) << "[Mesh::updateFaces] finish elements loop";

    if ( !(nDim == 1 && nOrder > 1) ) // not remove faces in this case (can be internal point)
    {
        tic();
        // clean faces not connected to an element
        face_iterator f_it = this->beginFace();
        face_iterator f_en = this->endFace();
        for ( ; f_it != f_en; )
        {
            auto const& face = f_it->second;
            // cleanup the face data structure :
            if ( !face.isConnectedTo0() )
            {
                DLOG( INFO ) << "removing face id : " << face.id();
                // remove all faces that are not connected to any elements
                f_it = this->eraseFace( f_it );
                //++f_it;
            }
            else
                ++f_it;
        }
        toc( "Mesh.updateEntitiesCoDimensionOne.clean_faces", FLAGS_v > 1 );
    }

    LOG( INFO ) << "We have now " << std::distance( this->beginFace(), this->endFace() ) << " faces in the mesh";
}

template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::updateEntitiesCoDimensionOneMinimal()
{
    rank_type currentPid = MeshBase<IndexT>::worldComm().localRank();
    rank_type numPartition = MeshBase<IndexT>::worldComm().localSize();

    const uint16_type _numLocalFaces = this->numLocalFaces();
    std::vector<uint16_type> myfToP( face_type::numVertices * _numLocalFaces );
    for ( uint16_type j = 0; j < _numLocalFaces; j++ )
    {
        for ( int f = 0; f < face_type::numVertices; ++f )
            myfToP[j * face_type::numVertices + f] = ( nDim == 1 ) ? j : /*iv->*/ element_type::fToP( j, f );
    }

    const bool updateComponentAddElements = this->components().test( MESH_ADD_ELEMENTS_INFO );

    // special trick (at this time, we need also to add active faces which touch interprocess boundary)
    // else mpi dof table can be not build
    std::unordered_set<size_type> ptsTouchInterProcess;
    if ( numPartition > 1 )
    {
        auto rangeGhostElement = this->ghostElements();
        auto ghostelt_it = std::get<0>( rangeGhostElement );
        auto const ghostelt_en = std::get<1>( rangeGhostElement );
        for ( ; ghostelt_it != ghostelt_en; ++ghostelt_it )
        {
            auto const& ghostelt = boost::unwrap_ref( *ghostelt_it );
            for ( uint16_type vLocId = 0; vLocId < element_type::numVertices; ++vLocId )
            {
                auto const& thept = ghostelt.point( vLocId );
                if ( thept.processId() != invalid_rank_type_value )
                    ptsTouchInterProcess.insert( thept.id() );
            }
        }
    }

    bool faceinserted = false;

    Eigen::Matrix<size_type, element_type::numVertices, 1> pointIdInElt;
    std::vector<size_type> lids( face_type::numVertices );

    element_iterator iv, en;
    boost::tie( iv, en ) = this->elementsRange();
    size_type nElt = std::distance( iv, en );

    typedef std::unordered_map<std::vector /*set*/<size_type>, std::tuple<element_type*, uint16_type, face_type*>, Feel::HashTables::HasherContainers<size_type>> pointstoface_container_type;
    pointstoface_container_type _faces( nElt * _numLocalFaces );
    typename pointstoface_container_type::iterator _faceit;

    size_type next_face = 0;

    // First We check if we have already Faces stored
    if ( !this->faces().empty() )
    {
        face_iterator __it = this->beginFace();
        face_iterator __en = this->endFace();
        LOG( INFO ) << "We have " << std::distance( __it, __en ) << " faces in the database";
        for ( ; __it != __en; )
        {
            auto& face = __it->second;
            for ( int f = 0; f < face_type::numVertices; ++f )
            {
                lids[f] = face.point( f ).id();
            }
            std::sort( lids.begin(), lids.end() );

            boost::tie( _faceit, faceinserted ) = _faces.emplace( std::piecewise_construct,
                                                                  std::forward_as_tuple( lids ),
                                                                  std::forward_as_tuple( nullptr, invalid_uint16_type_value, &face ) );

            DVLOG_IF( 2, faceinserted ) << "added face with id " << face.id() << "\n";
            DVLOG_IF( 2, !faceinserted ) << "not added face with id " << face.id()
                                         << " was already face with id = " << std::get<2>( _faceit->second )->id() << "\n";

            if ( faceinserted )
            {
                next_face = std::max( next_face, face.id() + 1 );
                ++__it;
            }
            else
            {
                // if ( __it->hasMarker() && ( !facePtr->hasMarker() || ( __it->marker() != facePtr->marker() ) ) )
                //     facePtr->setMarker( __it->marker().value() );
                __it = this->eraseFace( __it );
            }
        }
    }

    // get an ordering by process (need in // for have a deterministic choice of face connection (0 and 1) and therfore to define the permutations )
    std::map<rank_type, std::vector<element_type*>, std::less<rank_type>> elementsByProcOrdering;
    for ( ; iv != en; ++iv )
    {
        auto& elt = iv->second;
        rank_type pid = elt.processId();
        elementsByProcOrdering[pid].push_back( &elt );
    }
    // iterator over all element
    for ( auto const& elementsByProcPair : elementsByProcOrdering )
    {
        for ( element_type* eltPtr : elementsByProcPair.second )
        {
            element_type& elt = *eltPtr;

            const size_type eltId = elt.id();
            const rank_type eltPid = elt.processId();

            for ( uint16_type f = 0; f < element_type::numVertices; ++f )
                pointIdInElt[f] = elt.point( f ).id();

            for ( uint16_type j = 0; j < _numLocalFaces; j++ )
            {
                for ( uint16_type f = 0; f < face_type::numVertices; ++f )
                    lids[f] = pointIdInElt[myfToP[j * face_type::numVertices + f]];
                std::sort( lids.begin(), lids.end() );

                boost::tie( _faceit, faceinserted ) = _faces.emplace( std::piecewise_construct,
                                                                      std::forward_as_tuple( lids ),
                                                                      std::forward_as_tuple( eltPtr, j, nullptr ) );
                if ( faceinserted )
                {
                    //DVLOG(2) << "Connection0 face id: " << next_face << " to element id: " << eltId << " local face id: " << j << "\n";
                }
                else // already stored
                {
                    //DVLOG(2) << "Connection1 face id: " << faceIdRegistered << " to element id: " << eltId << " local face id: " << j << "\n";
                    auto& f2eVal = _faceit->second;
                    element_type* elt0 = std::get<0>( f2eVal );
                    face_type* facePtr = std::get<2>( f2eVal );
                    if ( !elt0 )
                    {
                        // found an face already stored by the mesh importer
                        CHECK( facePtr != nullptr ) << "invalid face";
                        facePtr->setConnection0( boost::make_tuple( eltPtr, j ) );
                        facePtr->setProcessId( eltPid );
                        facePtr->setOnBoundary( true, face_type::nDim );
                        std::get<0>( _faceit->second ) = &elt;
                        std::get<1>( _faceit->second ) = j;
                        elt.setFace( j, *facePtr );
                    }
                    else
                    {
                        // found an internal faces
                        elt0->setNeighbor( std::get<1>( f2eVal ), eltId );
                        elt.setNeighbor( j, elt0->id() );

                        if ( facePtr ) // face already built
                        {
                            facePtr->setOnBoundary( false );
                            // force processId equal to M_comm.rank() if face on interprocessfaces
                            if ( ( elt0->processId() == currentPid ) || ( elt.processId() == currentPid ) )
                                facePtr->setProcessId( currentPid );
                            facePtr->setConnection1( boost::make_tuple( eltPtr, j ) );
                            elt.setFace( j, *facePtr );
                            if ( updateComponentAddElements )
                                facePtr->addElement( eltId );
                            // fix duplication of point in connection1 with 3d mesh at order 3 and 4
                            this->fixPointDuplicationInHOMesh( elt, *facePtr, mpl::bool_ < nDim == 3 && nOrder >= 3 > () );
                        }
                        else //maybe a new face to build
                        {
                            // add internal faces in ghost element and connected to interprocess faces
                            if ( numPartition > 1 )
                            {
                                bool isFaceOnGhostElt = ( elt0->processId() != currentPid ) || ( eltPid != currentPid );
                                bool isFaceConnectedToInterProcess = isFaceOnGhostElt;
                                if ( !isFaceConnectedToInterProcess )
                                {
                                    for ( uint16_type vLocId = 0; vLocId < face_type::numVertices; ++vLocId )
                                    {
                                        if ( ptsTouchInterProcess.find( lids[vLocId] ) != ptsTouchInterProcess.end() )
                                        {
                                            isFaceConnectedToInterProcess = true;
                                            break;
                                        }
                                    }
                                }
                                if ( isFaceConnectedToInterProcess )
                                {
                                    face_type newFace;
                                    newFace.setId( next_face++ );
                                    newFace.setProcessIdInPartition( currentPid );
                                    // force processId equal to M_comm.rank() if face on interprocessfaces
                                    if ( ( elt0->processId() == currentPid ) || ( elt.processId() == currentPid ) )
                                        newFace.setProcessId( currentPid );
                                    newFace.setOnBoundary( false );
                                    uint16_type posElt0 = std::get<1>( f2eVal );
                                    newFace.setConnection0( boost::make_tuple( elt0, posElt0 ) );
                                    newFace.setConnection1( boost::make_tuple( eltPtr, j ) );
                                    for ( size_type k = 0; k < face_type::numPoints; ++k )
                                        newFace.setPoint( k, elt0->point( /*ele.*/ element_type::fToP( posElt0, k ) ) );

                                    if ( updateComponentAddElements )
                                    {
                                        newFace.addElement( elt0->id() );
                                        newFace.addElement( eltId );
                                    }

                                    auto res = this->addFace( newFace );
                                    auto& faceInserted = res.first->second;
                                    elt0->setFace( posElt0, faceInserted );
                                    elt.setFace( j, faceInserted );
                                    std::get<2>( _faceit->second ) = &faceInserted;
                                    // fix duplication of point in connection1 with 3d mesh at order 3 and 4
                                    this->fixPointDuplicationInHOMesh( elt, faceInserted, mpl::bool_ < nDim == 3 && nOrder >= 3 > () );
                                }
                            } // numPartition > 1
                        }
                    }
                }
            } // local face
        }     // element loop
    }

    // add boundary faces if not present
    boost::tie( iv, en ) = this->elementsRange();
    for ( ; iv != en; ++iv )
    {
        element_type& elt = iv->second;
        const size_type eltId = elt.id();
        const rank_type eltPid = elt.processId();
        for ( uint16_type j = 0; j < _numLocalFaces; j++ )
        {
            if ( elt.facePtr( j ) )
                continue;
            if ( elt.hasNeighbor( j ) )
                continue;
            face_type newFace;
            newFace.setId( next_face++ );
            newFace.setProcessIdInPartition( currentPid );
            newFace.setProcessId( eltPid );
            newFace.setOnBoundary( true, face_type::nDim );
            newFace.setConnection0( boost::make_tuple( &elt, j ) );
            for ( size_type k = 0; k < face_type::numPoints; ++k )
                newFace.setPoint( k, elt.point( /*ele.*/ element_type::fToP( j, k ) ) );
            if ( updateComponentAddElements )
                newFace.addElement( eltId );
            auto res = this->addFace( newFace );
            auto& faceInserted = res.first->second;
            elt.setFace( j, faceInserted );
        }
    }

    // clean faces not connected to an element
    face_iterator f_it = this->beginFace();
    face_iterator f_en = this->endFace();
    for ( ; f_it != f_en; )
    {
        auto const& face = f_it->second;
        // cleanup the face data structure :
        if ( !face.isConnectedTo0() )
        {
            DLOG( INFO ) << "removing face id : " << face.id();
            // remove all faces that are not connected to any elements
            f_it = this->eraseFace( f_it );
            //++f_it;
        }
        else
            ++f_it;
    }
}

#if 0
// this version is a lit bit less efficiently
template<typename Shape, typename T, int Tag, typename IndexT>
void
 Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::updateAdjacencyElements()
{
    VLOG(2) << "Compute adjacency graph\n";
    //boost::unordered_map<std::set<int>, size_type > _faces;
    //typename boost::unordered_map<std::set<int>, size_type >::iterator _faceit;
    boost::unordered_map<std::set<size_type>, size_type > _faces;
    typename boost::unordered_map<std::set<size_type>, size_type >::iterator _faceit;
    //std::map<std::set<size_type>, size_type > _faces;
    //typename std::map<std::set<size_type>, size_type >::iterator _faceit;
    std::vector< std::vector<std::tuple<size_type,rank_type,uint16_type> > > f2e;
    size_type next_face = 0;
    bool faceinserted = false;

    element_iterator iv,  en;
    boost::tie( iv, en ) = this->elementsRange();

    const uint16_type _numLocalFaces = this->numLocalFaces();
    for ( ; iv != en; ++iv )
    {
        const size_type eltId = iv->id();
        const rank_type eltPid = iv->processId();
        for ( uint16_type j = 0; j < _numLocalFaces; j++ )
        {
            std::set<size_type> s;
            //std::set<int> s;
            for ( uint16_type f = 0; f < face_type::numVertices; ++f )
            {
                if ( nDim == 1 )
                    s.insert( iv->point( j ).id() );
                else
                    s.insert( iv->point( iv->fToP( j, f ) ).id() );
            }

            boost::tie( _faceit, faceinserted ) = _faces.insert( std::make_pair( s, next_face ) );

            auto const& faceIdRegistered = _faceit->second;
            DVLOG(2) << "Connect face id: " << faceIdRegistered << " to element id: " << eltId << " local face id: " << j << " process id:" << eltPid << "\n";

            if ( faceinserted )
            {
                f2e.push_back( std::vector<std::tuple<size_type,rank_type,uint16_type> >(1,std::make_tuple(eltId,eltPid,j) ) );
                ++next_face;
            }
            else // already stored
            {
                f2e[faceIdRegistered].reserve(2);
                f2e[faceIdRegistered].push_back( std::make_tuple(eltId,eltPid,j) );

                // erase face desc in map (maybe improve other access)
                //_faces.erase( _faceit );
            }

        } // local face
    } // element loop

    for ( auto const& faceDatas : f2e )
    {
        int nConnectedElts = faceDatas.size();
        if ( nConnectedElts <= 1 ) continue;
        for ( int k1=0;k1<nConnectedElts;++k1 )
        {
            auto const& data1 = faceDatas[k1];
            element_iterator eltIt1 = this->elementIterator( std::get<0>( data1 ), std::get<1>( data1 ) );
            uint16_type j = std::get<2>( data1 );
            for ( int k2=0;k2<nConnectedElts;++k2 )
            {
                if ( k1==k2 ) continue;
                auto const& data2 = faceDatas[k2];
                size_type eltId2 = std::get<0>( data2 );
                rank_type eltPid2 = std::get<1>( data2 );
                this->elements().modify( eltIt1, [&j,&eltId2,&eltPid2]( element_type& elt ) { elt.setNeighbor(j,eltId2 ); });
            }
        }
    }
}
#else
template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::updateAdjacencyElements()
{
    VLOG( 2 ) << "Compute adjacency graph\n";
    //typedef std::unordered_map<std::vector/*set*/<size_type>, size_type, Feel::HashTables::HasherContainers<size_type> > pointstoface_container_type;
    //typedef std::unordered_map<std::vector/*set*/<size_type>, std::tuple<size_type,uint16_type>, Feel::HashTables::HasherContainers<size_type> > pointstoface_container_type;
    typedef std::unordered_map<std::vector /*set*/<size_type>, std::tuple<element_type*, uint16_type>, Feel::HashTables::HasherContainers<size_type>> pointstoface_container_type;

    std::vector<uint16_type> myfToP( face_type::numVertices * this->numLocalFaces() );
    for ( uint16_type j = 0; j < this->numLocalFaces(); j++ )
    {
        for ( int f = 0; f < face_type::numVertices; ++f )
            myfToP[j * face_type::numVertices + f] = ( nDim == 1 ) ? j : /*iv->*/ element_type::fToP( j, f );
    }

    bool faceinserted = false;

    Eigen::Matrix<size_type, element_type::numVertices, 1> pointIdInElt;
    std::vector<size_type> lids( face_type::numVertices );

    const uint16_type _numLocalFaces = this->numLocalFaces();

    element_iterator iv, en;
    boost::tie( iv, en ) = this->elementsRange();
    size_type nElt = std::distance( iv, en );

    pointstoface_container_type _faces;
    _faces.reserve( nElt * _numLocalFaces );
    typename pointstoface_container_type::iterator _faceit;

    for ( ; iv != en; ++iv )
    {
        element_type& elt = iv->second;
        const size_type eltId = elt.id();

        for ( uint16_type f = 0; f < element_type::numVertices; ++f )
            pointIdInElt[f] = elt.point( f ).id();

        for ( uint16_type j = 0; j < _numLocalFaces; j++ )
        {
            for ( uint16_type f = 0; f < face_type::numVertices; ++f )
                lids[f] = pointIdInElt[myfToP[j * face_type::numVertices + f]];
            std::sort( lids.begin(), lids.end() );

            auto [_faceit, faceinserted] = _faces.try_emplace( lids, &elt, j );

            if ( faceinserted )
            {
                //DVLOG(2) << "Connection0 face id: " << next_face << " to element id: " << eltId << " local face id: " << j << "\n";
            }
            else // already stored
            {
                //DVLOG(2) << "Connection1 face id: " << faceIdRegistered << " to element id: " << eltId << " local face id: " << j << "\n";
                auto const& f2eVal = _faceit->second;
                element_type* elt2 = std::get<0>( f2eVal );
                elt2->setNeighbor( std::get<1>( f2eVal ), eltId );
                elt.setNeighbor( j, elt2->id() );
            }
        } // local face
    }     // element loop
}

#endif

template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::fixPointDuplicationInHOMesh( element_type& elt, face_type const& face, mpl::true_ )
{
#if 0
    CHECK( nOrder == 3 || nOrder == 4 ) << "this fix with nOrder " << nOrder << " not implement, only order 3 and 4\n";
    if ( nOrder == 3 )
    {
        uint16_type startPt = face_type::numVertices*face_type::nbPtsPerVertex + face_type::numEdges*face_type::nbPtsPerEdge;
        CHECK( startPt == 9 && face_type::numPoints==10 ) << "invalid case\n";

        // easy case : only one point in face
        uint16_type idx = element_type::fToP( face.pos_second(), startPt );
        auto const& ptOnFace = face.point( startPt );
        size_type ptIdToCheck = elt.point( idx ).id();
        if ( ptIdToCheck != ptOnFace.id() )
        {
            // erase duplicate point
            this->erasePoint( this->pointIterator( ptIdToCheck ) );
            // fix point coordinate
            elt.setPoint( idx, ptOnFace );
        }
    }
    else if ( nOrder == 4 )
    {
        const uint16_type startPt = face_type::numVertices*face_type::nbPtsPerVertex + face_type::numEdges*face_type::nbPtsPerEdge;
        CHECK( startPt == 12 && face_type::numPoints==15 ) << "invalid case\n";

        // get a correspondence between the face and elt point numbering
        std::set<uint16_type> idxDone;
        std::vector<uint16_type> idxDistribution(3,invalid_uint16_type_value);
        for ( uint16_type k = startPt; k < face_type::numPoints; ++k )
        {
            auto const& ptOnFaceElt = elt.point( element_type::fToP( face.pos_second(), k ) );

            double dist=0, distMin= 1000*elt.hFace( face.pos_second() );// INT_MAX;
            uint16_type idxFaceNear = invalid_uint16_type_value;

            for ( uint16_type idxV = startPt ; idxV < face_type::numPoints ; ++idxV )
            {
                auto const& ptOnFace = face.point( idxV );
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
            uint16_type idx = element_type::fToP( face.pos_second(), startPt + i );
            size_type ptIdToCheck = elt.point( idx ).id();
            auto const& ptOnFace = face.point( idxDistribution[i] );

            if ( ptIdToCheck != ptOnFace.id() )
            {
                // erase duplicate point
                this->erasePoint( this->pointIterator( ptIdToCheck ) );
                // fix point coordinate
                elt.setPoint( idx, ptOnFace );
            }
        }
    }
#else
    CHECK( 0 ) << "fixPointDuplicationInHOMesh not implemented yet ";
#endif
}
template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::fixPointDuplicationInHOMesh( element_type& elt, face_type const& face, mpl::false_ )
{
}

template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
template <typename TheShape>
void  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::modifyEdgesOnBoundary( face_type& face, std::enable_if_t<TheShape::nDim == 3>* )
{
    // loop over face edges
    for ( int f = 0; f < face_type::numEdges; ++f )
    {
        auto& edge = face.edge( f );
        if ( !edge.isOnBoundary() )
        {
            edge.setOnBoundary( true, 1 );
        }
    }
    DVLOG( 3 ) << "We have " << nelements( boundaryedges( this ) ) << " boundary edges";
}
template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
template <typename TheShape>
void  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::modifyEdgesOnBoundary( face_type& f, std::enable_if_t<TheShape::nDim != 3>* )
{
}

template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
template <typename TheShape>
bool  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::modifyElementOnBoundaryFromEdge( element_type& elt, std::enable_if_t<TheShape::nDim != 3>* )
{
    return false;
}
template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
template <typename TheShape>
bool  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::modifyElementOnBoundaryFromEdge( element_type& elt, std::enable_if_t<TheShape::nDim == 3>* )
{
    // in 3D check if the edges of the element touch the boundary
    bool isOnBoundary = false;
    for ( size_type j = 0; j < elt.nEdges(); j++ )
    {
        auto edgePtr = elt.edgePtr( j );
        if ( edgePtr )
            isOnBoundary |= edgePtr->isOnBoundary();
    }
    if ( isOnBoundary )
    {
        // the element touches the boundary with just an edge
        elt.setOnBoundary( isOnBoundary, 1 );
        return isOnBoundary;
    }
    return isOnBoundary;
}

template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::updateOnBoundary()
{
    // first go through all the faces and set the points of the boundary
    // faces to be on the boundary
    LOG( INFO ) << "update boundary points...";
    DVLOG( 2 ) << "Initially we have " << nelements( boundarypoints( *this ) ) << " boundary points";
    for ( auto it = this->beginFace(), en = this->endFace(); it != en; ++it )
    {
        auto& face = it->second;
        if ( face.isOnBoundary() )
        {
            modifyEdgesOnBoundary( face );
            // loop over face points
            for ( int f = 0; f < face_type::numPoints; ++f )
            {
                auto& point = face.point( f );
                if ( !point.isOnBoundary() )
                {
                    point.setOnBoundary( true, 0 );
                }
            }
        }
    }
    DVLOG( 2 ) << "We have now " << nelements( boundarypoints( *this ) ) << " boundary points";
    LOG( INFO ) << "update boundary elements...";
    // loop through faces to set the elements having a face on the boundary
    for ( auto iv = this->beginElement(), en = this->endElement();
          iv != en; ++iv )
    {
        bool isOnBoundary = false;
        auto& elt = iv->second;
        // first check if a face is on the boundary
        for ( size_type j = 0; j < elt.nTopologicalFaces(); j++ )
        {
            auto facePtr = elt.facePtr( j );
            if ( facePtr )
                isOnBoundary |= facePtr->isOnBoundary();
        }

        if ( isOnBoundary )
        {
            VLOG( 3 ) << "checking " << elt.nTopologicalFaces() << " faces, isOnBoundary: " << isOnBoundary << " face_type::nDim: " << face_type::nDim;
            elt.setOnBoundary( isOnBoundary, face_type::nDim );
            ;
            // go to the next element, no need to look further
            continue;
        }
        bool e_modified = modifyElementOnBoundaryFromEdge( elt );
        // go to next element if element is on boundary
        if ( e_modified )
        {
            continue;
        }

        // finally check if a point of the element touches the boundary
        for ( size_type j = 0; j < elt.nPoints(); j++ )
        {
            isOnBoundary |= elt.point( j ).isOnBoundary();
        }

        if ( isOnBoundary )
        {
            VLOG( 3 ) << "checking " << elt.nPoints() << " points, isOnBoundary: " << isOnBoundary;
            elt.setOnBoundary( isOnBoundary, 0 );
        }
    } // loop over the elements
#if !defined( NDEBUG )
    DVLOG( 2 ) << "[updateOnBoundary] We have " << nelements( boundaryelements( this->shared_from_this() ) )
               << " elements sharing a point, a edge or a face with the boundary in the database";
    // for( auto const& e : this->boundaryElements( 0, 2, 0 ) )
    auto rangebe = this->boundaryElements( 0, 2, 0 );
    auto itbe = std::get<0>( rangebe );
    auto enbe = std::get<1>( rangebe );
    for ( ; itbe != enbe; ++itbe )
    {
        auto const& e = boost::unwrap_ref( *itbe );
        VLOG( 3 ) << "boundary element : " << e.id()
                  << " entity on boundary max dim  " << e.boundaryEntityDimension()
                  << " process id : " << e.processId();
    }
#endif
}
template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::removeFacesFromBoundary( std::initializer_list<uint16_type> markers )
{
    std::for_each( markers.begin(), markers.end(),
                   [this]( uint16_type marker ) {
                       auto range = markedfaces( this->shared_from_this(), boost::any( marker ) );
                       LOG( INFO ) << "removing " << nelements( range ) << " faces marked " << marker << " from boundary faces\n";

                       for ( auto it = range.begin(), en = range.end(); it != en; ++it )
                       {
                           auto const& face = boost::unwrap_ref( *it );
                           if ( face.isOnBoundary() )
                           {
                               DVLOG( 3 ) << "removing face " << face.id() << "\n";

                               auto& faceModified = this->faceIterator( face )->second;
                               faceModified.setOnBoundary( false );

                               //face_type f = *it2;
                               //f.setOnBoundary( false );
                               //this->faces().replace( it2, f );
                               CHECK( face.isOnBoundary() == false ) << " face should not be on the boundary anymore\n";
                           }
                       }
                   } );
}



template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void
Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::updateParallelData()
{
    // some contexts and data used in update
    using container_elements_type = std::vector<std::tuple<size_type,size_type>>; // (id of ghost elt, id of active element)
    using container_faces_type = std::map<size_type,std::tuple<std::tuple<size_type,uint16_type>,
                                                               std::map<uint16_type, typename face_type::marker_type>,
                                                               bool >>;// (id of face -> ( (id of element,face id in element), markers, isOnboundary )
    using container_edges_type = std::map<size_type,std::tuple<std::tuple<size_type,uint16_type>,
                                                               std::map<uint16_type, typename edge_type::marker_type> >>;// (id of edge -> ((id of element,edge id in element), markers )
    using container_points_type = std::map<size_type,std::tuple<std::tuple<size_type,uint16_type>,
                                                                std::map<uint16_type, typename point_type::marker_type> >>;// (id of point -> ( (id of element,point id in element0), markers )

    std::map<rank_type, std::tuple<container_elements_type,
                                   container_faces_type,
                                   container_edges_type,
                                   container_points_type>> dataToSend, dataToRecv;
    static constexpr int tuple_id_info_elements = 0;
    static constexpr int tuple_id_info_faces = 1;
    static constexpr int tuple_id_info_edges = 2;
    static constexpr int tuple_id_info_points = 3;

    bool meshHasUpdateFaces = this->components().test( MESH_UPDATE_FACES ) || this->components().test( MESH_UPDATE_FACES_MINIMAL );
    bool meshHasUpdateEdges = nDim == 3 && this->components().test( MESH_UPDATE_EDGES );


    // update interprocess entities and neighbor subdomains
    // TODO optimisation if we have interprocessfaces

    std::unordered_map<index_type,std::tuple<bool,std::set<rank_type>>> pointsInterprocessDetection; // ( pt id -> ( isOnActiveElt, isOnGhostEltRanks ) )
    std::unordered_map<index_type,std::tuple<bool,std::set<rank_type>>> edgesInterprocessDetection; // ( edge id -> ( isOnActiveElt, isOnGhostEltRanks ) )
    pointsInterprocessDetection.reserve( std::distance( this->beginOrderedPoint(),
                                                        this->endOrderedPoint() ) );
    if constexpr ( nDim == 3 )
        edgesInterprocessDetection.reserve( std::distance( this->beginOrderedEdge(),
                                                           this->endOrderedEdge() ) );
    auto itPointIpDetect = pointsInterprocessDetection.begin();
    auto itEdgeIpDetect = edgesInterprocessDetection.begin();
    auto itoe = this->beginOrderedElement();
    auto enoe = this->endOrderedElement();
    for ( ; itoe!=enoe;++itoe )
    {
        auto const& elt = unwrap_ref( *itoe );

        // nothing to do if no neighbor process
        if ( elt.idInOthersPartitions().empty() )
            continue;
        // update neighbor subdomains
        for ( auto const& [pid,eid] : elt.idInOthersPartitions() )
            this->addNeighborSubdomain( pid );

        for ( uint16_type n=0; n < elt.nPoints(); n++ )
        {
            auto const& point = elt.point( n );
            std::tie( itPointIpDetect,std::ignore ) = pointsInterprocessDetection.try_emplace( point.id(), false, std::set<rank_type>{} );
            if ( elt.isGhostCell() )
                std::get<1>( itPointIpDetect->second ).insert( elt.processId() );
            else
                std::get<0>( itPointIpDetect->second ) = true;
        }

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
    }

    //this->updateInterprocessPoints( std::move( pointsInterprocessDetection ) );
    this->updateInterprocessPoints( pointsInterprocessDetection );
    if constexpr ( nDim == 3 )
    {
        if ( meshHasUpdateEdges )
            this->updateInterprocessEdges( edgesInterprocessDetection );
    }
    //------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------------------//

    auto initialNeighborSubdomains = this->neighborSubdomains();
    rank_type currentProcessId = MeshBase<IndexT>::worldComm().localRank();

    auto iv = this->beginOrderedElement();
    auto const en = this->endOrderedElement();
    for ( ; iv != en; ++iv )
    {
        auto& ghostelt = boost::unwrap_ref( *iv );
        if ( !ghostelt.isGhostCell() )
            continue;
        const rank_type ghosteltPid = ghostelt.processId();

        // elements
        const size_type idEltInOtherPartition = ghostelt.idInOthersPartitions( ghosteltPid );
        std::get<tuple_id_info_elements>( dataToSend[ghosteltPid] ).push_back( std::make_tuple( ghostelt.id(), idEltInOtherPartition ) );

        // faces
        auto & dataToSendFaces = std::get<tuple_id_info_faces>( dataToSend[ghosteltPid] );
        for ( uint16_type j = 0; j < ghostelt.nTopologicalFaces(); j++ )
        {
            if ( !ghostelt.facePtr( j ) )
                continue;
            auto & face = ghostelt.face( j );
            // do nothing if already updated
            auto itFindFace = dataToSendFaces.find(face.id());
            if ( itFindFace != dataToSendFaces.end() )
                continue;
            // set process id at interprocess (use min rank)
            if ( face.isInterProcessDomain( currentProcessId ) )
                face.setProcessId( std::min( face.element0().processId(), face.element1().processId() ) );
            // update face data
            std::tie( itFindFace, std::ignore ) = dataToSendFaces.emplace( face.id(), typename container_faces_type::mapped_type{} );
            std::get<0>( itFindFace->second ) = std::make_tuple(idEltInOtherPartition, j );
            std::get<1>( itFindFace->second ) = face.markers();
            std::get<2>( itFindFace->second ) = face.isOnBoundary();
        }

        // edges
        if constexpr ( nDim == 3 )
        {
            auto & dataToSendEdges = std::get<tuple_id_info_edges>( dataToSend[ghosteltPid] );
            for ( size_type j = 0; j < ghostelt.nEdges(); j++ )
            {
                if ( !ghostelt.edgePtr( j ) )
                    continue;
                auto & edge = ghostelt.edge( j );
                // do nothing if already updated
                auto itFindEdge = dataToSendEdges.find( edge.id() );
                if ( itFindEdge != dataToSendEdges.end() )
                    continue;
                // set process id at interprocess (use min rank)
                auto [isInterprocess,itInterprocess] = this->findInterprocessEdges( edge.id() );
                if ( isInterprocess )
                    edge.setProcessId( std::min( currentProcessId, *itInterprocess->second.begin() ) );
                // update edge data
                std::tie( itFindEdge, std::ignore ) = dataToSendEdges.emplace( edge.id(), typename container_edges_type::mapped_type{} );
                std::get<0>( itFindEdge->second ) = std::make_tuple(idEltInOtherPartition, j );
                std::get<1>( itFindEdge->second ) = edge.markers();
            }
        }

        // points
        if ( nDim > 1 || !meshHasUpdateFaces )
        {
            auto & dataToSendPoints = std::get<tuple_id_info_points>( dataToSend[ghosteltPid] );
            for ( size_type j = 0; j < ghostelt.nPoints(); j++ )
            {
                auto & point =  ghostelt.point( j );
                // do nothing if already updated
                auto itFindPoint = dataToSendPoints.find(point.id());
                if ( itFindPoint != dataToSendPoints.end() )
                    continue;
                // set process id at interprocess (use min rank)
                auto [isInterprocess,itInterprocess] = this->findInterprocessPoints( point.id() );
                if ( isInterprocess )
                    point.setProcessId( std::min( currentProcessId, *itInterprocess->second.begin() ) );
                // update point data
                std::tie( itFindPoint, std::ignore ) = dataToSendPoints.emplace( point.id(), typename container_points_type::mapped_type{} );
                std::get<0>( itFindPoint->second ) = std::make_tuple(idEltInOtherPartition, j );
                std::get<1>( itFindPoint->second ) = point.markers();
            }
        }
    }

    //------------------------------------------------------------------------------------------------//
    //int neighborSubdomains = this->neighborSubdomains().size();
    int nbMaxRequest = (nDim+1)*2*initialNeighborSubdomains.size();//neighborSubdomains;

    std::vector<mpi::request> reqs( nbMaxRequest );
    int countRequest = 0;

    // get size of data to transfer
    std::map<rank_type,std::array<std::size_t,nDim==1?3:nDim+1>> sizeRecv, sizeSended;
    uint16_type sizeDataPointIndex = nDim==1?2:nDim;
    for ( rank_type neighborRank : initialNeighborSubdomains )
    {
        sizeSended[neighborRank][0] = std::get<tuple_id_info_elements>( dataToSend[neighborRank] ).size();
        sizeSended[neighborRank][1] = std::get<tuple_id_info_faces>( dataToSend[neighborRank] ).size();
        if constexpr ( nDim == 3 )
            sizeSended[neighborRank][2] = std::get<tuple_id_info_edges>( dataToSend[neighborRank] ).size();
        //if constexpr ( nDim > 2  )
        sizeSended[neighborRank][/*nDim*/sizeDataPointIndex] = std::get<tuple_id_info_points>( dataToSend[neighborRank] ).size();
        reqs[countRequest++] = MeshBase<IndexT>::worldComm().localComm().isend( neighborRank, 0, sizeSended[neighborRank] );
        reqs[countRequest++] = MeshBase<IndexT>::worldComm().localComm().irecv( neighborRank, 0, sizeRecv[neighborRank] );
    }
    // wait all requests
    mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
    countRequest = 0;

    for ( rank_type neighborRank : initialNeighborSubdomains )
    {
        std::size_t nRecvDataElement = sizeRecv[neighborRank][0];
        std::get<0>( dataToRecv[neighborRank] ).resize( nRecvDataElement );
        if ( nRecvDataElement > 0 )
            reqs[countRequest++] = MeshBase<IndexT>::worldComm().localComm().irecv( neighborRank, tuple_id_info_elements, std::get<tuple_id_info_elements>( dataToRecv[neighborRank] ).data(), nRecvDataElement );
        std::size_t nSendDataElement = sizeSended[neighborRank][0];
        if ( nSendDataElement > 0 )
            reqs[countRequest++] = MeshBase<IndexT>::worldComm().localComm().isend( neighborRank, tuple_id_info_elements, std::get<tuple_id_info_elements>( dataToSend[neighborRank] ).data(), nSendDataElement );

        std::size_t nRecvDataFaces = sizeRecv[neighborRank][1];
        if ( nRecvDataFaces > 0 )
            reqs[countRequest++] = MeshBase<IndexT>::worldComm().localComm().irecv( neighborRank, tuple_id_info_faces, std::get<tuple_id_info_faces>( dataToRecv[neighborRank] ) );
        std::size_t nSendDataFaces = sizeSended[neighborRank][1];
        if ( nSendDataFaces > 0 )
            reqs[countRequest++] = MeshBase<IndexT>::worldComm().localComm().isend( neighborRank, tuple_id_info_faces, std::get<tuple_id_info_faces>( dataToSend[neighborRank] ) );

        if constexpr ( nDim == 3 )
        {
            std::size_t nRecvDataEdges = sizeRecv[neighborRank][2];
            if ( nRecvDataEdges > 0 )
                reqs[countRequest++] = MeshBase<IndexT>::worldComm().localComm().irecv( neighborRank, tuple_id_info_edges, std::get<tuple_id_info_edges>( dataToRecv[neighborRank] ) );
            std::size_t nSendDataEdges = sizeSended[neighborRank][2];
            if ( nSendDataEdges > 0 )
                reqs[countRequest++] = MeshBase<IndexT>::worldComm().localComm().isend( neighborRank, tuple_id_info_edges, std::get<tuple_id_info_edges>( dataToSend[neighborRank] ) );
        }

        //if constexpr ( nDim > 2 )
        if ( nDim > 1 || !meshHasUpdateFaces )
        {
            std::size_t nRecvDataPoints = sizeRecv[neighborRank][sizeDataPointIndex/*nDim*/];
            if ( nRecvDataPoints > 0 )
                reqs[countRequest++] = MeshBase<IndexT>::worldComm().localComm().irecv( neighborRank, tuple_id_info_points, std::get<tuple_id_info_points>( dataToRecv[neighborRank] ) );
            std::size_t nSendDataPoints = sizeSended[neighborRank][sizeDataPointIndex/*nDim*/];
            if ( nSendDataPoints > 0 )
                reqs[countRequest++] = MeshBase<IndexT>::worldComm().localComm().isend( neighborRank, tuple_id_info_points, std::get<tuple_id_info_points>( dataToSend[neighborRank] ) );
        }

    }
    // wait all requests
    mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
    countRequest = 0;

    //------------------------------------------------------------------------------------------------//
    M_neighbor_processors.clear();
    // update mesh data (step 1)
    std::map<size_type,element_type*> eltUpdated;
    std::map<size_type,face_type*> facesUpdated;
    std::map<size_type,edge_type*> edgesUpdated;
    std::map<size_type,point_type*> pointsUpdated;
    for ( auto const& [rankRecv,dataEntities] : dataToRecv )
    {
        this->addNeighborSubdomain( rankRecv );
        // elements
        for ( auto const& [ eltIdRecv,eltIdCurrent ] : std::get<tuple_id_info_elements>( dataEntities ) )
        {
            CHECK( this->hasElement( eltIdCurrent ) ) << "invalid elt id" << eltIdCurrent;
            auto& theelt = this->elementIterator( eltIdCurrent )->second;
            theelt.setIdInOtherPartitions( rankRecv, eltIdRecv );
            eltUpdated.emplace( theelt.id(), &theelt );
        }

        // faces
        for ( auto const& [faceIdRecv,dataFace] : std::get<tuple_id_info_faces>( dataEntities ) )
        {
            auto const& [eltIdCurrent,faceIdInElt] = std::get<0>( dataFace );
            auto const& mapMarkers = std::get<1>( dataFace );
            bool faceIsOnBoundary = std::get<2>( dataFace );

            CHECK( this->hasElement( eltIdCurrent ) ) << "invalid elt id" << eltIdCurrent;
            auto& theelt = this->elementIterator( eltIdCurrent )->second;
            if ( !theelt.facePtr( faceIdInElt ) )
                continue;
            auto & theface = theelt.face( faceIdInElt );
            theface.setIdInOtherPartitions( rankRecv, faceIdRecv );
            for ( auto [markerType,markerValues] : mapMarkers )
                theface.addMarker( markerType, markerValues );
            if ( !faceIsOnBoundary )
                theface.setOnBoundary( false );
            facesUpdated.emplace( theface.id(), &theface );
        }
        if constexpr ( nDim == 3 )
        {
            // edges
            for ( auto const& [edgeIdRecv,dataEdge] : std::get<tuple_id_info_edges>( dataEntities ) )
            {
                auto const& [eltIdCurrent,edgeIdInElt] = std::get<0>( dataEdge );
                auto const& mapMarkers = std::get<1>( dataEdge );

                CHECK( this->hasElement( eltIdCurrent ) ) << "invalid elt id" << eltIdCurrent;
                auto& theelt = this->elementIterator( eltIdCurrent )->second;
                if ( !theelt.edgePtr( edgeIdInElt ) )
                    continue;
                auto & edge = theelt.edge( edgeIdInElt );
                edge.setIdInOtherPartitions( rankRecv, edgeIdRecv );
                for ( auto [markerType,markerValues] : mapMarkers )
                    edge.addMarker( markerType, markerValues );
                edgesUpdated.emplace( edge.id(), &edge );
            }
        }

        if ( nDim > 1 || !meshHasUpdateFaces )
        {
            // points
            for ( auto const& [pointIdRecv,dataPoint] : std::get<tuple_id_info_points>( dataEntities ) )
            {
                auto const& [eltIdCurrent,pointIdInElt] = std::get<0>( dataPoint );
                auto const& mapMarkers = std::get<1>( dataPoint );

                CHECK( this->hasElement( eltIdCurrent ) ) << "invalid elt id" << eltIdCurrent;
                auto& theelt = this->elementIterator( eltIdCurrent )->second;
                if ( !theelt.pointPtr( pointIdInElt ) )
                    continue;
                auto & point = theelt.point( pointIdInElt );
                point.setIdInOtherPartitions( rankRecv, pointIdRecv );
                for ( auto [markerType,markerValues] : mapMarkers )
                    point.addMarker( markerType, markerValues );
                pointsUpdated.emplace( point.id(), &point );
            }
        }
    }

    // Step 2 : active elements and subentity included are properly built, we can send info to ghost elements

    // prepare data to send (second step)
    using container_step2_elements_type = std::map<size_type, std::tuple<std::map<rank_type,size_type>,
                                                                         std::vector<uint16_type> > >; // (id of elt, (id in other part mapping), (face permutation in elt (only 3D)))
    using container_step2_faces_type = std::map<size_type, std::tuple<std::map<rank_type,size_type>,
                                                                      std::map<uint16_type, typename face_type::marker_type>,
                                                                      rank_type,bool> >; // (id of face, (id in other part mapping),markers,processId,isOnBoundary)
    using container_step2_edges_type = std::map<size_type, std::tuple<std::map<rank_type,size_type>,
                                                                      std::map<uint16_type, typename edge_type::marker_type>,
                                                                      rank_type> >; // (id of edge, (id in other part mapping),markers,processId)
    using container_step2_points_type = std::map<size_type, std::tuple<std::map<rank_type,size_type>,
                                                                       std::map<uint16_type, typename point_type::marker_type>,
                                                                       rank_type> >; // (id of point, (id in other part mapping),markers,processId)
    std::map<rank_type, std::tuple<container_step2_elements_type,
                                   container_step2_faces_type,
                                   container_step2_edges_type,
                                   container_step2_points_type>> dataToSendStep2, dataToRecvStep2;

    auto prepareDataToSendStep2 =
        [&dataToSendStep2]( auto const& mapEntities, auto tupleId ){
            static constexpr int tuple_id = std::decay_t<decltype(tupleId)>::value;
            for ( auto & [entityId,entityPtr] : mapEntities )
            {
                for ( auto const& [pidGhost,eidGhost] : entityPtr->idInOthersPartitions() )
                {
                    auto & dataToSendStep2Elt = std::get<tuple_id>( dataToSendStep2[pidGhost] )[eidGhost];

                    // idInOthersPartitions
                    std::map<rank_type,size_type> tmpIiop;
                    for ( auto [r,e] : entityPtr->idInOthersPartitions() )
                        if ( r != pidGhost )
                            tmpIiop.emplace( r, e );
                    // NOTE: not need to add current id with element entities (already registered in ghost element)
                    if constexpr( tuple_id != tuple_id_info_elements )
                        tmpIiop.emplace( entityPtr->pidInPartition(), entityPtr->id() );
                    std::get<0>( dataToSendStep2Elt ) = std::move( tmpIiop );

                    if constexpr( tuple_id == tuple_id_info_elements )
                    {
                        if constexpr ( nDim == 3 )
                        {
                            std::vector<uint16_type> facePermutations( entityPtr->nTopologicalFaces(), invalid_v<uint16_type> );
                            for ( size_type j = 0; j < entityPtr->nTopologicalFaces(); j++ )
                            {
                                if ( !entityPtr->facePtr( j ) )
                                    continue;
                                auto const& face = entityPtr->face( j );
                                facePermutations[j] = entityPtr->facePermutation( j ).value();
                            }
                            std::get<1>( dataToSendStep2Elt ) = std::move( facePermutations );
                        }
                    }
                    else if constexpr( tuple_id == tuple_id_info_faces )
                    {
                        std::get<1>( dataToSendStep2Elt ) = entityPtr->markers();
                        std::get<3>( dataToSendStep2Elt ) = entityPtr->isOnBoundary();
                        std::get<2>( dataToSendStep2Elt ) = entityPtr->processId();
                    }
                    else if constexpr( tuple_id == tuple_id_info_edges || tuple_id == tuple_id_info_points )
                    {
                        std::get<1>( dataToSendStep2Elt ) = entityPtr->markers();
                        std::get<2>( dataToSendStep2Elt ) = entityPtr->processId();
                    }
                }
            }
        };
    prepareDataToSendStep2( eltUpdated, std::integral_constant<int,tuple_id_info_elements>{} );
    prepareDataToSendStep2( facesUpdated,std::integral_constant<int,tuple_id_info_faces>{} );
    if constexpr ( nDim == 3 )
        prepareDataToSendStep2( edgesUpdated,std::integral_constant<int,tuple_id_info_edges>{} );
    //if constexpr ( nDim > 2 )
    if ( nDim > 1 || !meshHasUpdateFaces )
        prepareDataToSendStep2( pointsUpdated,std::integral_constant<int,tuple_id_info_points>{} );
    for ( rank_type neighborRank : initialNeighborSubdomains )
    {
        reqs[countRequest++] = MeshBase<IndexT>::worldComm().localComm().isend( neighborRank, 0, dataToSendStep2[neighborRank] );
        reqs[countRequest++] = MeshBase<IndexT>::worldComm().localComm().irecv( neighborRank, 0, dataToRecvStep2[neighborRank] );
    }
    // wait all requests
    mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
    countRequest = 0;

    //------------------------------------------------------------------------------------------------//
    // update ghost elements and subentities (step 2)
    for ( auto const& [rankRecv,dataEntities] : dataToRecvStep2 )
    {
        // elements
        for ( auto const& [eltId,dataElt] : std::get<tuple_id_info_elements>( dataEntities ) )
        {
            CHECK( this->hasElement( eltId ) ) << "invalid elt id" << eltId;
            auto& elt = this->elementIterator( eltId )->second;
            for ( auto [r,e] : std::get<0>( dataElt ) )
                elt.setIdInOtherPartitions( r, e );
            if constexpr ( nDim == 3 )
            {
                for ( uint16_type j = 0; j < elt.nTopologicalFaces(); j++ )
                {
                    if ( !elt.facePtr( j ) )
                        continue;
                    auto & face = elt.face( j );
                    elt.setFacePermutation( j, typename element_type::face_permutation_type( std::get<1>( dataElt ).at(j) ) );
                }
            }
        }
        // faces
        for ( auto const& [faceId,dataFace] : std::get<1>( dataEntities ) )
        {
            CHECK( this->hasFace( faceId ) ) << "invalid face id" << faceId;
            auto& face = this->faceIterator( faceId )->second;
            for ( auto [r,e] : std::get<0>( dataFace ) )
                face.setIdInOtherPartitions( r, e );
            for ( auto const& [markerType,markerValues] : std::get<1>( dataFace ) )
                face.addMarker( markerType, markerValues );
            // // update process id
            face.setProcessId( std::get<2>( dataFace ) );

            //face.setOnBoundary( std::get<3>( dataFace ) );
            if ( !std::get<3>( dataFace ) )
                face.setOnBoundary( false );
            if constexpr ( nDim == 1 )
            {
                for ( auto const& [pidOther,eidOther] : face.idInOthersPartitions() )
                    this->addNeighborSubdomain( pidOther );
            }
        }
        // edges
        if constexpr ( nDim == 3 )
        {
            for ( auto const& [edgeId,dataEdge] : std::get<tuple_id_info_edges>( dataEntities ) )
            {
                CHECK( this->hasEdge( edgeId ) ) << "invalid edge id" << edgeId;
                auto& edge = this->edgeIterator( edgeId )->second;
                for ( auto [r,e] : std::get<0>( dataEdge ) )
                    edge.setIdInOtherPartitions( r, e );
                for ( auto const& [markerType,markerValues] : std::get<1>( dataEdge ) )
                    edge.addMarker( markerType, markerValues );
                // update process id
                edge.setProcessId( std::get<2>( dataEdge ) );
            }
        }
        // points
        if ( nDim > 1 || !meshHasUpdateFaces )
        {
            for ( auto const& [pointId,dataPoint] : std::get<tuple_id_info_points>( dataEntities ) )
            {
                CHECK( this->hasPoint( pointId ) ) << "invalid point id" << pointId;
                auto& point = this->pointIterator( pointId )->second;
                for ( auto [r,e] : std::get<0>( dataPoint ) )
                    point.setIdInOtherPartitions( r, e );
                for ( auto const& [markerType,markerValues] : std::get<1>( dataPoint ) )
                    point.addMarker( markerType, markerValues );
                // update process id
                point.setProcessId( std::get<2>( dataPoint ) );
                // TODO check interprocess point : process id  should be already good or put invalid_value at the begin for ghost non interprocess
                for ( auto const& [pidOther,eidOther] : point.idInOthersPartitions() )
                    this->addNeighborSubdomain( pidOther );
            }
        }
    }

}





template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::check() const
{
    if ( nDim != nRealDim )
        return;
#if !defined( NDEBUG )
    VLOG( 2 ) << "[Mesh::check] numLocalFaces = " << this->numLocalFaces() << "\n";
    //element_iterator iv = this->beginElementWithProcessId( MeshBase<IndexT>::worldComm().localRank() );
    //element_iterator en = this->endElementWithProcessId( MeshBase<IndexT>::worldComm().localRank() );
    auto iv = this->beginElement();
    auto en = this->endElement();
    size_type nEltInMesh = std::distance( iv, en );
    VLOG( 2 ) << "[Mesh::check] nEltInMesh = " << nEltInMesh;

    //boost::tie( iv, en ) = this->elementsRange();
    for ( ; iv != en; ++iv )
    {
        element_type const& elt = iv->second;

        if ( this->components().test( MESH_UPDATE_FACES ) )
        {
            for ( size_type j = 0; j < this->numLocalFaces(); j++ )
            {
                //FEELPP_ASSERT( elt.facePtr( j ) )( j )( elt.id() ).error( "invalid element face check" );
                VLOG( 2 ) << "------------------------------------------------------------\n";
                if ( elt.hasFace( j ) )
                    VLOG( 2 ) << "Element : " << elt.id() << " face lid: " << j << " face gid:  "
                              << elt.face( j ).id();
                else
                    VLOG( 2 ) << "Element : " << elt.id() << " face lid: " << j;
            }
        }

        if ( this->components().test( MESH_UPDATE_FACES ) || this->components().test( MESH_UPDATE_FACES_MINIMAL ) ||
             this->components().test( MESH_UPDATE_ELEMENTS_ADJACENCY ) )
        {
            size_type counter = 0;
            for ( uint16_type ms = 0; ms < elt.nNeighbors(); ms++ )
            {
                if ( elt.neighbor( ms ) != invalid_v<size_type> )
                    ++counter;
            }
            VLOG( 2 ) << "[Mesh::check] element " << elt.id() << " number of neighbors: " << counter << "\n";
            if ( elt.nNeighbors() > 0 )
            {
                //FEELPP_ASSERT( counter >= 1 || nEltInMesh==1 )( elt.id() )( elt.nNeighbors() )( counter ).warn( "invalid neighboring data" );
            }
        }

#if 0

        for ( size_type j = 0; j < ( size_type )element_type::numEdges; ++j )
        {
            FEELPP_ASSERT( elt.edgePtr( j ) )( j )( elt.id() ).error( "invalid element edge check" );
            VLOG(2) << "------------------------------------------------------------\n";
            VLOG(2) << "Element : " << elt.id() << " edge lid: " << j << " edge gid:  " << elt.edge( j ).id() << "\n";

        }

#endif
    }

    // faces check
    auto rangeBoundaryFaces = this->facesOnBoundary();
    auto itf = std::get<0>( rangeBoundaryFaces );
    auto ite = std::get<1>( rangeBoundaryFaces );
    std::map<int, int> nf;
    for ( ; itf != ite; ++itf )
    {
        face_type const& __face = boost::unwrap_ref( *itf );
        DLOG_IF( WARNING, !__face.isConnectedTo0() && !__face.isConnectedTo1() ) << "face not connected to an element face:" << __face << "\n";
        if ( !__face.isConnectedTo0() && !__face.isConnectedTo1() )
        {
            auto it = nf.find( (int)__face.id() );
            if ( it == nf.end() )
                nf[(int)__face.id()] = 1;
            else
                it->second += 1;
        }
        DLOG_IF( WARNING, !__face.isConnectedTo0() && __face.isConnectedTo1() ) << "face  connected to element 1 but not element 0. face:" << __face << "\n";
        DLOG_IF( WARNING, __face.isConnectedTo0() && __face.isConnectedTo1() ) << "invalid boundary face (connected to 2 elements):" << __face << "\n";

        DLOG_IF( WARNING, __face.isConnectedTo0() && !__face.element( 0 ).facePtr( __face.pos_first() ) ) << "invalid face in element, face: " << __face << " in element " << __face.element( 0 ) << "\n";
    }
    auto itt = nf.begin();
    auto ent = nf.end();
    for ( ; itt != ent; ++itt )
    {
        LOG( INFO ) << "face with id " << itt->first << " not attached " << itt->second << "\n";
    }
#endif
}

#if 0
template<typename Shape, typename T, int Tag, typename IndexT>
void
 Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::findNeighboringProcessors()
{
    // Don't need to do anything if there is
    // only one processor.
    if ( MeshBase<IndexT>::worldComm().localSize() == 1 )
        return;

#ifdef FEELPP_HAS_MPI

    M_neighboring_processors.clear();

    // Get the bounding sphere for the local processor
    Sphere bounding_sphere = processorBoundingSphere ( *this, MeshBase<IndexT>::worldComm().localRank() );

    // Just to be sure, increase its radius by 10%.  Sure would suck to
    // miss a neighboring processor!
    bounding_sphere.setRadius( bounding_sphere.radius()*1.1 );

    // Collect the bounding spheres from all processors, test for intersection
    {
        std::vector<double>
        send ( 4,                         0 ),
             recv ( 4*MeshBase<IndexT>::worldComm().localSize(), 0 );

        send[0] = bounding_sphere.center()( 0 );
        send[1] = bounding_sphere.center()( 1 );
        send[2] = bounding_sphere.center()( 2 );
        send[3] = bounding_sphere.radius();

        MPI_Allgather ( &send[0], send.size(), MPI_DOUBLE,
                        &recv[0], send.size(), MPI_DOUBLE,
                        MeshBase<IndexT>::worldComm().localComm() );


        for ( unsigned int proc=0; proc<MeshBase<IndexT>::worldComm().localSize(); proc++ )
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
        VLOG(2) << "Processor " << MeshBase<IndexT>::worldComm().localRank() << " intersects:\n";

        for ( unsigned int p=0; p< M_neighboring_processors.size(); p++ )
            VLOG(2) << " - proc " << M_neighboring_processors[p] << "\n";
    }

#endif
}
#endif

template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void
     Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::checkLocalPermutation( mpl::bool_<true> ) const
{
    bool mesh_well_oriented = true;
    std::vector<size_type> list_of_bad_elts;

    for ( element_const_iterator iv = this->beginElement();
          iv != this->endElement(); ++iv )
    {
        element_type const& elt = iv->second;

        if ( !elt.isAnticlockwiseOriented() )
        {
            mesh_well_oriented = false;
            list_of_bad_elts.push_back( elt.id() );
        }
    }

    if ( mesh_well_oriented )
        VLOG( 2 ) << "Local numbering in the elements is OK . \n";

    else
    {
        std::for_each( list_of_bad_elts.begin(),
                       list_of_bad_elts.end(),
                       []( size_type const& e ) {
                           LOG_FIRST_N( WARNING, 10 ) << "element is not anticlockwise oriented(wrong local numbering): " << e << "\n";
                       } );
    }
}

#if 0
template<typename Shape, typename T, int Tag, typename IndexT>
void
 Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::checkAndFixPermutation(  )
{
    element_type __element = *this->beginElement()->second;

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
        element_type __element = elt->second;
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
#endif

template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::send( int p, int tag )
{
    encode();
    VLOG( 1 ) << "sending markername\n";
    //this->comm().send( p, tag, M_markername.size() );
    VLOG( 1 ) << "sending markername size: " << this->M_markername.size() << "\n";
    for ( auto const& m : this->M_markername )
    {
        VLOG( 1 ) << "sending key: " << m.first << "\n";
        //this->comm().send( p, tag, m.first );
        VLOG( 1 ) << "sending value\n";
        //this->comm().send( p, tag, m.second );
    }
}

template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::recv( int p, int tag )
{
    VLOG( 1 ) << "receiving markername\n";
    //this->comm().recv( p, tag, this->M_markername );
    int s = 0;
    //this->comm().recv( p, tag, s );
    VLOG( 1 ) << "receiving markername size: " << s << "\n";
    for ( int i = 0; i < s; ++i )
    {
        std::string k;
        VLOG( 1 ) << "receiving key\n";
        //this->comm().recv( p, tag, k );
        VLOG( 1 ) << "receiving key:" << k << "\n";
        std::vector<int> v;
        VLOG( 1 ) << "receiving value\n";
        //this->comm().recv( p, tag, v );
        VLOG( 1 ) << "receiving value: " << v[0] << "," << v[1] << "\n";
        // this->M_markername[k]=v;
    }

    //decode();
}

template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
typename  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::element_iterator
 Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::eraseElement( element_iterator position, bool modify )
{
#if 0
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
#else
    element_iterator eit;
    return eit;

#endif
}
template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::encode()
{
    //std::cout<<"encode=   " << MeshBase<IndexT>::worldComm().localSize() << std::endl;

    M_enc_pts.clear();
    for ( auto pt_it = this->beginPoint(), pt_en = this->endPoint(); pt_it != pt_en; ++pt_it )
    {
        auto const& pt = pt_it->second;
        std::vector<double> pts( 3, 0 );
        pts[0] = pt.node()[0];
        if ( mesh_type::nRealDim >= 2 )
            pts[1] = pt.node()[1];
        if ( mesh_type::nRealDim >= 3 )
            pts[2] = pt.node()[2];
        M_enc_pts[pt.id() + 1] = boost::make_tuple( pt.isOnBoundary(), pt.tags(), pts );
    }
    M_enc_elts.clear();

    auto allmarkedfaces = boundaryfaces( this->shared_from_this() );
    auto face_it = allmarkedfaces.begin();
    auto face_end = allmarkedfaces.end();

    int elem_number = 1;

    GmshOrdering<element_type> ordering;

    GmshOrdering<face_type> ordering_face;
    // save the faces

    for ( ; face_it != face_end; ++face_it )
    {
        auto const& curface = boost::unwrap_ref( *face_it );
        std::vector<int> faces;
        faces.push_back( ordering_face.type() );
        faces.push_back( 4 + curface.numberOfPartitions() );
        faces.push_back( curface.hasMarker() ? curface.marker().value() : 0 );
        faces.push_back( curface.hasMarker2() ? curface.marker2().value() : 0 );
        faces.push_back( curface.numberOfPartitions() );
        faces.push_back( curface.processId() );
        for ( auto const& [pid,eid] : curface.neighborProcessIds() )
            faces.push_back(-pid);
        for ( uint16_type p = 0; p < face_type::numPoints; ++p )
            faces.push_back( curface.point( ordering_face.fromGmshId( p ) ).id() + 1 );

        M_enc_faces[curface.id()] = faces;
    } // faces

    auto rangeElements = this->elementsWithProcessId();
    auto elt_it = std::get<0>( rangeElements );
    auto elt_en = std::get<1>( rangeElements );
    for ( ; elt_it != elt_en; ++elt_it )
    {
        auto const& curelt = boost::unwrap_ref( *elt_it );
        std::vector<int> elts;
        elts.push_back( ordering.type() );
        elts.push_back( 4 + curelt.numberOfPartitions() );
        elts.push_back( curelt.hasMarker() ? curelt.marker().value() : 0 );
        elts.push_back( curelt.hasMarker2() ? curelt.marker2().value() : 0 );
        elts.push_back( curelt.numberOfPartitions() );
        elts.push_back( curelt.processId() );
        for ( auto const& [pid,eid] : curelt.neighborProcessIds() )
            elts.push_back(-pid);
        for ( uint16_type p = 0; p < element_type::numPoints; ++p )
            elts.push_back( curelt.point( ordering.fromGmshId( p ) ).id() + 1 );

        M_enc_elts[curelt.id()] = elts;
    } // elements
    //std::cout<<"encode=   " << MeshBase<IndexT>::worldComm().localSize() << std::endl;
}

template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void  Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::decode()
{
#if 0
    std::vector<int> mapWorld(MeshBase<IndexT>::worldComm().size());
    for(int cpu=0; cpu < MeshBase<IndexT>::worldComm().size(); ++cpu)
        mapWorld[cpu] = cpu;
    WorldComm worldcomm(mapWorld);




    // std::cout<<"decode=   " << MeshBase<IndexT>::worldComm().size() << std::endl;
    //std::cout<<"decode=   " << worldcomm.subWorldComm().localRank() << std::endl;
    this->setWorldComm(worldcomm.subWorldComm());
#else
    LOG( INFO ) << "decode=   " << MeshBase<IndexT>::worldComm().size() << "\n";
    LOG( INFO ) << "decode=   " << MeshBase<IndexT>::worldComm().subWorldComm().localRank() << "\n";
#endif
    static const uint16_type npoints_per_face = ( face_type::numVertices * face_type::nbPtsPerVertex +
                                                  face_type::numEdges * face_type::nbPtsPerEdge +
                                                  face_type::numFaces * face_type::nbPtsPerFace );

    static const uint16_type npoints_per_element = element_type::numPoints;
    for ( auto pt_it = M_enc_pts.begin(), pt_en = M_enc_pts.end();
          pt_it != pt_en; ++pt_it )
    {
        node_type __n( nRealDim );

        for ( uint16_type j = 0; j < nRealDim; ++j )
            __n[j] = pt_it->second.template get<2>()[j];

        point_type __pt( pt_it->first - 1, __n, pt_it->second.template get<0>() );
        __pt.setOnBoundary( pt_it->second.template get<0>() );
        __pt.setTags( pt_it->second.template get<1>() );

        this->addPoint( __pt );
    }
    for ( auto face_it = M_enc_faces.begin(), face_en = M_enc_faces.end();
          face_it != face_en; ++face_it )
    {
        face_type pf;
        GmshOrdering<face_type> ordering;
        FEELPP_ASSERT( ordering.type() == face_it->second[0] )
        ( ordering.type() )( face_it->second[0] ).error( "invalid mesh face type" );
        std::vector<int> tags( face_it->second[1] );
        for ( int i = 0; i < tags.size(); ++i )
            tags[i] = face_it->second[2 + i];
        pf.setTags( tags );
        pf.setId( this->numFaces() );
        pf.setProcessIdInPartition( MeshBase<IndexT>::worldComm().localRank() );
        pf.setProcessId( MeshBase<IndexT>::worldComm().localRank() );
        //pf.setIdInPartition( MeshBase<IndexT>::worldComm().localRank(),pf.id() );

        const int shift = face_it->second[1] + 1;
        for ( uint16_type jj = 0; jj < npoints_per_face; ++jj )
        {
            pf.setPoint( ordering.fromGmshId( jj ), this->point( face_it->second[shift + jj] - 1 ) );
        }
        this->addFace( std::move( pf ) );
    }

    for ( auto elt_it = M_enc_elts.begin(), elt_en = M_enc_elts.end();
          elt_it != elt_en; ++elt_it )
    {
        element_type pv;
        GmshOrdering<element_type> ordering;
        FEELPP_ASSERT( ordering.type() == elt_it->second[0] )
        ( ordering.type() )( elt_it->second[0] ).error( "invalid mesh element type" );
        std::vector<int> tags( elt_it->second[1] );
        for ( int i = 0; i < tags.size(); ++i )
            tags[i] = elt_it->second[2 + i];
        pv.setTags( tags );
        pv.setProcessIdInPartition( MeshBase<IndexT>::worldComm().localRank() );
        pv.setProcessId( MeshBase<IndexT>::worldComm().localRank() );
        //pv.setIdInPartition( MeshBase<IndexT>::worldComm().localRank(),pv.id() );

        const int shift = elt_it->second[1] + 1;
        for ( uint16_type jj = 0; jj < npoints_per_element; ++jj )
        {
            pv.setPoint( ordering.fromGmshId( jj ), this->point( elt_it->second[shift + jj] - 1 ) );
        }
        this->addElement( pv );
#if 0
        __idGmshToFeel=pv.id();
        auto theelt = mesh->elementIterator( pv.id(), pv.partitionId() );
        mesh->elements().modify( theelt, Feel::detail::updateIdInOthersPartitions( MeshBase<IndexT>::worldComm().localRank(), pv.id() ) );
#endif
    }
    LOG( INFO ) << "distance  elts: " << std::distance( this->beginElement(), this->endElement() ) << "\n";
    LOG( INFO ) << "distance faces: " << std::distance( this->beginFace(), this->endFace() ) << "\n";
    // LOG(INFO) << "distance marker faces: "<< std::distance( this->beginFaceWithMarker(), this->endFaceWithMarker() ) << "\n";
    // LOG(INFO) << "distance marker2 faces: "<< std::distance( this->beginFaceWithMarker2(), this->endFaceWithMarker2() ) << "\n";

    //this->components().set ( MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
    //this->updateForUse();
    //std::cout<<"decode=   " << MeshBase<IndexT>::worldComm().localSize() << std::endl;
    this->setSubStructuring( true );
}

template <typename Shape, typename T, int Tag1, int Tag2 = Tag1, int TheTag = Tag1>
std::shared_ptr<Mesh<Shape, T, TheTag>>
merge( std::shared_ptr<Mesh<Shape, T, Tag1>> m1,
       std::shared_ptr<Mesh<Shape, T, Tag2>> m2 )
{
    std::shared_ptr<Mesh<Shape, T, TheTag>> m( new Mesh<Shape, T, TheTag> );

    // add the points

    size_type shift_p = std::distance( m1->beginPoint(), m1->endPoint() );
    for ( auto it = m1->beginPoint(), en = m1->endPoint(); it != en; ++it )
    {
        m->addPoint( *it );
    }
    for ( auto it = m2->beginPoint(), en = m2->endPoint(); it != en; ++it )
    {
        auto p = *it;
        // shift id
        p.setId( shift_p + p.id() );
        m->addPoint( p );
    }
    typedef typename Mesh<Shape, T, TheTag>::face_iterator face_iterator;
    typedef typename Mesh<Shape, T, TheTag>::element_iterator element_iterator;
    typedef typename Mesh<Shape, T, TheTag>::face_type face_type;
    typedef typename Mesh<Shape, T, TheTag>::element_type element_type;

    auto addface = [&]( std::shared_ptr<Mesh<Shape, T, TheTag>> m, face_iterator it, size_type shift ) {
        face_type pf = *it;
        pf.disconnect();
        pf.setOnBoundary( true );
        const uint16_type npoints_per_face = ( face_type::numVertices * face_type::nbPtsPerVertex +
                                               face_type::numEdges * face_type::nbPtsPerEdge +
                                               face_type::numFaces * face_type::nbPtsPerFace );
        for ( uint16_type jj = 0; jj < npoints_per_face; ++jj )
        {
            pf.setPoint( jj, m->point( shift + it->point( jj ).id() ) );
        }
        m->addFace( pf );
    };
    // add the faces
    size_type shift_f = std::distance( m1->beginFace(), m1->endFace() );
    for ( auto it = m1->beginFace(), en = m1->endFace(); it != en; ++it )
    {
        addface( m, it, 0 );
    }
    for ( auto it = m2->beginFace(), en = m2->endFace(); it != en; ++it )
    {
        // don't forget to shift the point id in the face
        addface( m, it, shift_p );
    }
    // add the elements
    auto addelement = [&]( std::shared_ptr<Mesh<Shape, T, TheTag>> m, element_iterator it, size_type shift_f, size_type shift_p ) {
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
            pf.setPoint( jj, m->point( shift_p + it->point( jj ).id() ) );
        }
        m->addElement( pf );
    };
    for ( auto it = m1->beginElement(), en = m1->endElement(); it != en; ++it )
    {
        addelement( m, it, 0, 0 );
    }
    for ( auto it = m2->beginElement(), en = m2->endElement(); it != en; ++it )
    {
        // don't forget to shift the point id in the face
        addelement( m, it, shift_f, shift_p );
    }
    m->components().set( MESH_UPDATE_EDGES | MESH_UPDATE_FACES | MESH_CHECK );
    m->updateForUse();
    return m;
}

namespace details
{
template<typename IndexT = uint32_type>
struct RemoveMarkerNameWithoutEntityForAllReduce : std::function<std::vector<std::vector<IndexT>>(
                                                                 std::vector<std::vector<IndexT>> const&,
                                                                 std::vector<std::vector<IndexT>> const& )>
{
    typedef std::vector<std::vector<IndexT>> cont_type;
    cont_type operator()( cont_type const& x, cont_type const& y ) const
    {
        cont_type res( x.size() );
        for ( int i = 0; i < x.size(); ++i )
        {
            int nMarkerByEntity = x[i].size();
            CHECK( y[i].size() == nMarkerByEntity ) << "size should be equal";
            res[i].resize( nMarkerByEntity );
            for ( int j = 0; j < nMarkerByEntity; ++j )
            {
                res[i][j] = x[i][j] + y[i][j];
            }
        }
        return res;
    }
};

} // namespace details


template <typename Shape, typename T, int Tag, typename IndexT, bool EnableSharedFromThis>
void
Mesh<Shape, T, Tag, IndexT, EnableSharedFromThis>::removeMarkerNameWithoutEntity()
{
    std::vector<std::vector<std::string>> markersEntity( nDim + 1 );
    std::vector<std::vector<size_type>> nMarkedEntity( nDim + 1 );
    std::vector<std::string> markersToRemove;
    std::vector<std::map<size_type,size_type>> mapMarkerFlagToContainer( nDim + 1 ); // (flag->index in container)
    for ( auto const& m : this->markerNames() )
    {
        std::string const& markerName = m.first;
        int markerDim = m.second[1];
        index_type markerFlag = m.second[0];

        mapMarkerFlagToContainer[markerDim][markerFlag] = nMarkedEntity[markerDim].size();
        markersEntity[markerDim].push_back( markerName );
        nMarkedEntity[markerDim].push_back( 0 );
    }

    auto updateMarkedEntityCounter =
        [&nMarkedEntity,&mapMarkerFlagToContainer]( auto itBegin, auto itEnd )
            {
                std::for_each(itBegin,itEnd, [&nMarkedEntity,&mapMarkerFlagToContainer]( auto const& pairIdElt )
                              {
                                  auto const& elt = pairIdElt.second;
                                  auto & nMarkedEntityCurrent = nMarkedEntity[elt.nDim];
                                  auto & mapMarkerFlagToContainerCurrent = mapMarkerFlagToContainer[elt.nDim];
                                  for ( auto const& [mtag,mflags] : elt.markers() )
                                  {
                                      for ( auto mflag : mflags )
                                      {
                                          auto itFind = mapMarkerFlagToContainerCurrent.find( mflag );
                                          if ( itFind != mapMarkerFlagToContainerCurrent.end() )
                                              ++(nMarkedEntityCurrent[itFind->second]);
                                      }
                                  }
                              });
            };

    if constexpr ( nDim >= 1 )
        updateMarkedEntityCounter( this->beginElement(), this->endElement() );
    if constexpr ( nDim == 3 )
        updateMarkedEntityCounter( this->beginEdge(), this->endEdge() );
    if constexpr ( nDim >= 2 )
        updateMarkedEntityCounter( this->beginFace(), this->endFace() );
    updateMarkedEntityCounter( this->beginPoint(), this->endPoint() );

    mpi::all_reduce( this->worldComm(), mpi::inplace( nMarkedEntity ), Feel::details::RemoveMarkerNameWithoutEntityForAllReduce<index_type>() );

    for ( int d = 0; d <= nDim; ++d )
        for ( int k = 0; k < nMarkedEntity[d].size(); ++k )
            if ( nMarkedEntity[d][k] == 0 )
                markersToRemove.push_back( markersEntity[d][k] );

    for ( std::string const& markerName : markersToRemove )
        this->M_markername.erase( markerName );
}

} // namespace Feel

#include <feel/feeldiscr/meshio.hpp>

namespace Feel
{

} // namespace Feel

#endif // __MESHIMPL_HPP
