/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-08-10

  Copyright (C) 2011-2020 Feel++ Consortium
  Copyright (C) 2005,2006 EPFL
  

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
   \file geoentity.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-08-10
 */
#ifndef __GeoEntity_H
#define __GeoEntity_H 1

#include <feel/feelmesh/simplex.hpp>
#include <feel/feelmesh/hypercube.hpp>
#include <feel/feelmesh/refentity.hpp>

namespace Feel
{

/**
   \class GeoEntity
   \brief base class for all geometric entities

   @author Christophe Prud'homme
   @see
*/
template<typename Entity, typename T = double, typename IndexT = uint32_type>
class GeoEntity
    :
        boost::equality_comparable<GeoEntity<Entity,T,IndexT> >,
    boost::less_than_comparable<GeoEntity<Entity,T,IndexT> >,
    boost::less_than_comparable<GeoEntity<Entity,T,IndexT>, IndexT>,
    public Entity
{
    static inline const uint16_type nBitShiftedGeoEntityContext = 0;
    static inline const uint16_type nBitShiftedReferenceGeometry = 2;
    static inline const uint16_type nBitShiftedReferenceShapes = 8;
public:


    /** @name Typedefs
     */
    //@{

    typedef Entity super;
    typedef GeoEntity<Entity,T,IndexT> GeoShape;
    typedef GeoEntity<Entity,T,IndexT> self_type;
    typedef T value_type;
    using index_type = IndexT;
    using size_type = index_type;
    typedef typename super::topological_face_type face_type;
    typedef face_type GeoBShape;
    typedef typename Entity::edge_permutation_type edge_permutation_type;
    typedef typename Entity::face_permutation_type face_permutation_type;

    static inline const size_type Shape = super::Shape;
    static inline const size_type Geometry = super::Geometry;

    static inline const uint16_type nDim = super::nDim;
    static inline const uint16_type nOrder = super::nOrder;
    static inline const uint16_type nRealDim = super::nRealDim;


    static inline const uint16_type numVertices = super::numVertices;
    static inline const uint16_type numFaces = super::numFaces;
    static inline const uint16_type numGeometricFaces = super::numGeometricFaces;
    static inline const uint16_type numTopologicalFaces = super::numTopologicalFaces;
    static inline const uint16_type numEdges = super::numEdges;
    static inline const uint16_type numNormals = super::numNormals;

    static inline const uint16_type numPoints = super::numPoints;
    static inline const uint16_type nbPtsPerVertex = super::nbPtsPerVertex;
    static inline const uint16_type nbPtsPerEdge = super::nbPtsPerEdge;
    static inline const uint16_type nbPtsPerFace = super::nbPtsPerFace;
    static inline const uint16_type nbPtsPerVolume = super::nbPtsPerVolume;

    typedef Entity convex_type;

    static inline const bool is_simplex = super::is_simplex;
    static inline const bool is_hypercube = super::is_hypercube;

    /**
     * helper class to construct the associated reference convex.
     */
    template<typename TT = double>
    struct reference_convex
    {
        typedef Reference<Entity, nDim, nOrder, nRealDim, TT> type;
    };
    template<typename TT = double>
    using reference_convex_type =  Reference<Entity, nDim, nOrder, nRealDim, TT>;

    using marker_type = Marker<flag_type/*uint16_type*/>;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    GeoEntity()
        :
        super(),
        M_id( 0 ),
        M_entity( (MESH_ENTITY_INTERNAL<<nBitShiftedGeoEntityContext) | (Geometry<<nBitShiftedReferenceGeometry) | (Shape<<nBitShiftedReferenceShapes) ),
        M_boundaryEntityDimension( invalid_uint16_type_value ),
        M_pid( invalid_rank_type_value ),
        M_pidInPartition( invalid_rank_type_value ),
        M_neighor_pids(),
        M_idInOtherPartitions(),
        M_elist(),
        M_elistGhost()
    {}

    explicit GeoEntity( size_type i,
                        size_type geometry = Geometry,
                        size_type shape = Shape,
                        size_type context = MESH_ENTITY_INTERNAL )
        :
        super(),
        M_id( i ),
        M_entity( (context<<nBitShiftedGeoEntityContext) | (geometry<<nBitShiftedReferenceGeometry) | (shape<<nBitShiftedReferenceShapes) ),
        M_boundaryEntityDimension( invalid_uint16_type_value ),
        M_pid( invalid_rank_type_value ),
        M_pidInPartition( invalid_rank_type_value ),
        M_neighor_pids(),
        M_idInOtherPartitions(),
        M_elist(),
        M_elistGhost()
    {}

    GeoEntity( GeoEntity const& __me ) = default;
    GeoEntity( GeoEntity && __me )
        :
        super( std::move( __me ) ),
        M_id( std::move( __me.M_id ) ),
        M_entity( std::move( __me.M_entity ) ),
        M_boundaryEntityDimension( std::move( __me.M_boundaryEntityDimension ) ),
        M_pid( std::move( __me.M_pid ) ),
        M_pidInPartition( std::move( __me.M_pidInPartition ) ),
        M_neighor_pids( std::move( __me.M_neighor_pids ) ),
        M_idInOtherPartitions( std::move( __me.M_idInOtherPartitions ) ),
        M_elist( std::move( __me.M_elist ) ),
        M_elistGhost( std::move( __me.M_elistGhost ) ),
        M_markers( std::move( __me.M_markers ) )
        {
            //std::cout << "GeoEntity moved ctor\n";
        }
    GeoEntity& operator=( GeoEntity const& __me ) = default;
    GeoEntity& operator=( GeoEntity && __me )
        {
            super::operator=( std::move( __me ) );
            M_id= std::move( __me.M_id );
            M_entity= std::move( __me.M_entity );
            M_boundaryEntityDimension= std::move( __me.M_boundaryEntityDimension );
            M_pid= std::move( __me.M_pid );
            M_pidInPartition= std::move( __me.M_pidInPartition );
            M_neighor_pids= std::move( __me.M_neighor_pids );
            M_idInOtherPartitions= std::move( __me.M_idInOtherPartitions );
            M_elist= std::move( __me.M_elist );
            M_elistGhost= std::move( __me.M_elistGhost );
            M_markers = std::move( __me.M_markers );
            //std::cout << "GeoEntity moved assign\n";
            return *this;
        }

    ~GeoEntity() override
    {}

    //@}

    /** @name Operator overloads
     */
    //@{
    bool operator==( GeoEntity const& e ) const
    {
        return M_id == e.id();
    }
    bool operator<( GeoEntity const& e ) const
    {
        return M_id < e.id();
    }

    bool operator<( size_type __i ) const
    {
        return M_id < __i;
    }

    //@}

    /** @name Accessors
     */
    //@{

    size_type id() const noexcept
    {
        return M_id;
    }


    /**
     * the dimension of the reference shape
     *
     *
     * @return the dimension of the reference shape
     */
    constexpr uint16_type refDim() const
    {
        return super::nDim;
    }

    /**
     * number of points on the reference shape
     *
     * @return the number of points on the reference shape
     */
    constexpr uint16_type nPoints() const
    {
        return super::numPoints;
    }

    /**
     * number of vertices on the reference shape
     *
     * @return the number of vertices on the reference shape
     */
    constexpr uint16_type nVertices() const
    {
        return super::numVertices;
    }

    /**
     * number of edges on the reference shape
     *
     * @return the number of edges on the reference shape
     */
    constexpr uint16_type nEdges() const
    {
        return super::numEdges;
    }

    /**
     * number of faces on the reference shape
     *
     * @return the number of edges on the reference shape
     */
    constexpr uint16_type nFaces() const
    {
        return super::numFaces;
    }

    /**
     * number of topological faces on the reference shape
     *
     * @return the number of topological faces on the reference shape
     */
    constexpr uint16_type nTopologicalFaces() const
    {
        return super::numTopologicalFaces;
    }

    /**
     * number of faces on the reference shape
     *
     * @return the number of edges on the reference shape
     */
    constexpr uint16_type nGeometricFaces() const
    {
        return super::numGeometricFaces;
    }

    /**
     * number of normals on the reference shape
     *
     * @return the number of normals on the reference shape
     */
    constexpr uint16_type nNormals() const
    {
        return super::numNormals;
    }


    /**
     *
     *
     *
     * @return true if the entoty has the shape \c __shape, false otherwise
     */
    bool hasShape( size_type __shape ) const
    {
        return M_entity.test( __shape << nBitShiftedReferenceShapes );
    }

    /**
     * @return true of the entity is a volume
     */
    bool isAVolume() const
    {
        return M_entity.test( GEOMETRY_VOLUME << nBitShiftedReferenceGeometry );
    }

    /**
     * @return true of the entity is a surface
     */
    bool isASurface() const
    {
        return M_entity.test( GEOMETRY_SURFACE << nBitShiftedReferenceGeometry );
    }

    /**
     * @return true of the entity is a line
     */
    bool isALine() const
    {
        return M_entity.test( GEOMETRY_LINE << nBitShiftedReferenceGeometry );
    }

    /**
     * @return true of the entity is a point
     */
    bool isAPoint() const
    {
        return M_entity.test( GEOMETRY_POINT << nBitShiftedReferenceGeometry );
    }

    /**
     * @return true of the entity is a shape point
     */
    bool isAPointShape() const
    {
        return M_entity.test( SHAPE_POINT << nBitShiftedReferenceShapes );
    }

    /**
     * @return true of the entity is a shape line
     */
    bool isALineShape() const
    {
        return M_entity.test( SHAPE_LINE << nBitShiftedReferenceShapes );
    }

    /**
     * @return true of the entity is a triangle shape
     */
    bool isATriangleShape() const
    {
        return M_entity.test( SHAPE_TRIANGLE << nBitShiftedReferenceShapes );
    }

    /**
     * @return true of the entity is a quadrangle
     */
    bool isAQuadrangleShape() const
    {
        return M_entity.test( SHAPE_QUAD << nBitShiftedReferenceShapes );
    }

    /**
     * @return true of the entity is a tetrahedra shape
     */
    bool isATetrahedraShape() const
    {
        return M_entity.test( SHAPE_TETRA << nBitShiftedReferenceShapes );
    }

    /**
     * @return true of the entity is a hexahedra
     */
    bool isAHexahedraShape() const
    {
        return M_entity.test( SHAPE_HEXA << nBitShiftedReferenceShapes );
    }

    /**
     * @return true if the shape is linear, false otherwise
     */
    bool isLinear() const
    {
        return M_entity.test( SHAPE_LINEAR << nBitShiftedReferenceShapes );
    }

    /**
     * @return true if the shape is bilinear, false otherwise
     */
    bool isBilinear() const
    {
        return M_entity.test( SHAPE_BILINEAR << nBitShiftedReferenceShapes );
    }

    /**
     * @return true if the shape is quadratic, false otherwise
     */
    bool isQuadratic() const
    {
        return M_entity.test( SHAPE_QUADRATIC << nBitShiftedReferenceShapes );
    }

    /**
     * @return true if the entity is internal, false otherwise
     */
    bool isInternal() const
    {
        return M_entity.test( MESH_ENTITY_INTERNAL << nBitShiftedGeoEntityContext );
    }


    /**
     * Tells if  item is on the boundary
     * @return true if on boundary, false otherwise
     */
    bool isOnBoundary() const noexcept
    {
        return M_entity.test( MESH_ENTITY_BOUNDARY << nBitShiftedGeoEntityContext );
    }

    /**
     * maximum dimension of the entity of the element touching the boundary
     */
    uint16_type boundaryEntityDimension() const noexcept
    {
        return M_boundaryEntityDimension;
    }
    /**
     * \return \c true if ghost cell, \c false otherwise
     */
    bool isGhostCell() const noexcept
    {
        //return (this->worldComm().localRank()!=M_pid);
        //mpi::communicator world;
        //return (world.rank()!=M_pid);
        return ( M_pidInPartition!=M_pid );
    }

    /**
     * \return the processor id of the entity
     */
    rank_type processId() const noexcept
    {
        return M_pid;
    }

    /**
     * set the processor id of the entity
     & \param pid processor id
     */
    void setProcessId( rank_type pid ) noexcept
    {
        M_pid = pid ;
    }

    /**
     * \return the processor id of the entity
     */
    rank_type pidInPartition() const noexcept
    {
        return M_pidInPartition;
    }
    /**
     * set the processor id of the entity
     & \param pid processor id
     */
    void setProcessIdInPartition( rank_type pid ) noexcept
    {
        M_pidInPartition = pid ;
    }

    /**
     * \return the partition id
     */
    rank_type partitionId() const noexcept
    {
        return M_pid;
    }

    /**
     * \return the number of partition the element is linked to including the
     * partition to which it belongs
     */
    rank_type numberOfPartitions() const noexcept
    {
        return static_cast<rank_type>(M_neighor_pids.size()+1);
    }

    /**
     * \return the number of partition the element is linked to
     */
    rank_type numberOfNeighborPartitions() const
    {
        return static_cast<rank_type>(M_neighor_pids.size());
    }

    /**
     * \return true if the element is linked to other partitions through one of
     * more of its faces
     */
    bool isLinkedToOtherPartitions() const
    {
        return M_neighor_pids.size() > 0;
    }

    /**
     * \return the number of partition the element is linked to
     */
    std::vector<rank_type> const& neighborPartitionIds() const
    {
        return M_neighor_pids;
    }
    /**
     * \return the number of partition the element is linked to
     */
    std::vector<rank_type> & neighborPartitionIds()
    {
        return M_neighor_pids;
    }
    /**
     * clear the neighbor partition ids container
     */
    void clearNeighborPartitionIds()
    {
        M_neighor_pids.clear();
    }

    /**
     * set id in a partition pid of the entity
     */
    FEELPP_DEPRECATED void setIdInOthersPartitions( rank_type pid, size_type id )
    {
        M_idInOtherPartitions.insert( std::make_pair( pid, id ) );
    }
    void setIdInOtherPartitions( rank_type pid, size_type id )
    {
        M_idInOtherPartitions.insert( std::make_pair( pid, id ) );
    }

    /**
     * set (partition,id) in other partitions of the entity
     */
    void setIdInOtherPartitions( std::map<rank_type,size_type> const& iop )
        {
            M_idInOtherPartitions = iop;
        }

    
    /**
     * set (partition,id) in other partitions of the entity
     */
    void setIdInOtherPartitions( std::map<rank_type,size_type>&& iop )
        {
            M_idInOtherPartitions = iop;
        }

    /**
     * \return the id of the entity in a partition pid
     */
    size_type idInOthersPartitions( rank_type pid ) const
    {
        DCHECK( M_idInOtherPartitions.find( pid )!=M_idInOtherPartitions.end() ) 
            << " local id " << this->id() << " is unknown for this partition " << pid << "\n";
        return M_idInOtherPartitions.find( pid )->second;
    }

    /**
     * \return idInOthersPartitions map
     */
    std::map<rank_type, size_type> const& idInOthersPartitions() const
    {
        return M_idInOtherPartitions;
    }
    /**
     * clear id in others partitions container
     */
    void clearIdInOthersPartitions()
    {
        M_idInOtherPartitions.clear();
    }

    /**
     * \return \c true if active, \c false otherwise
     *
     * \note for now it is a dummy function that returns always true,
     * will change when work on AMR starts
     */
    bool active() const
    {
        return true;
    }

    /**
     * \return the measure of the entity
     */
    virtual value_type measure() const { return value_type{0}; }

    //@}

    /** @name  Mutators
     */
    //@{
    void setId( size_type id )
    {
        M_id = id;
    }

    /**
     * set the boundary flag
     * @param b true if the item is on the boundary, false otherwise
     */
    void setOnBoundary( bool b, uint16_type ent_d = invalid_uint16_type_value )
    {
        if ( b )
        {
            M_entity.set( MESH_ENTITY_BOUNDARY << nBitShiftedGeoEntityContext );
            M_entity.clear( MESH_ENTITY_INTERNAL << nBitShiftedGeoEntityContext );
        }
        else
        {
            M_entity.clear( MESH_ENTITY_BOUNDARY << nBitShiftedGeoEntityContext );
            M_entity.set( MESH_ENTITY_INTERNAL << nBitShiftedGeoEntityContext );
        }
        M_boundaryEntityDimension = ent_d;
    }

    /**
     * \return the number of partition the element is linked to including the
     * partition to which it belongs
     */
    FEELPP_DEPRECATED void setNumberOfPartitions( uint16_type np )
    {
        CHECK( 0 ) << "Invalid call to setNumberOfPartitions()";
    }

    /**
     * set the number of partition the element is linked to
     */
    void setNumberOfNeighborPartitions( uint16_type nep )
    {
        M_neighor_pids.resize( nep );
    }

    /**
     * \return the number of partition the element is linked to
     */
    void setNeighborPartitionIds( std::vector<rank_type> const& npids )
    {
        M_neighor_pids = npids;
    }

    void addNeighborPartitionId( rank_type p )
    {
        if ( std::find( M_neighor_pids.begin(), M_neighor_pids.end(), p) == M_neighor_pids.end() )
        {
            M_neighor_pids.push_back(p);
        }
    }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * eToP(i,j) = localId of jth point on ith local edge
     */
    static uint16_type eToP( uint16_type const __localEdge, uint16_type const __point )
    {
        return super::e2p( __localEdge, __point );
    }

    /**
     * fToP(i,j) = localId of jth point on ith local edge
     */
    static uint16_type fToP( uint16_type const __localFace, uint16_type const __point )
    {
        return super::f2p( __localFace, __point );
    }

    /**
     * fToE(i,j) = localId of jth edge on ith local face
     */
    static uint16_type fToE( uint16_type const __localFace, uint16_type const __edge )
    {
        return super::f2e( __localFace, __edge );
    }

    /**
     * add a new element to which the point belongs
     */
    self_type& addElement( size_type e, int id_in_element = 0 )
    {
        M_elist.insert( std::make_pair(e,id_in_element) );
        return *this;
    }

    /**
     * \return the number of elements whom the point belongs to
     */
    size_type numberOfElements() const
    {
        return M_elist.size();
    }

    /**
     * \return the set of ids of elements whom the point belongs to
     */
    std::set<std::pair<size_type,uint16_type>> const& elements() const
    {
        return M_elist;
    }
    std::set<std::pair<size_type,uint16_type>> & elements()
    {
        return M_elist;
    }


    /**
     * add a new ghost element to which the point belongs
     */
    self_type& addElementGhost( rank_type proc, size_type e  )
    {
        M_elistGhost[proc].insert(e);
        return *this;
    }

    /**
     * \return the number of ghost elements whom the point belongs to
     */
    rank_type numberOfProcGhost() const
    {
        return M_elistGhost.size();
    }

    /**
     * \return the set of ids of ghost elements whom the point belongs to
     */
    std::map<rank_type,std::set<size_type> > const& elementsGhost() const
    {
        return M_elistGhost;
    }

    /**
     * \return true if geoentity is a subentity of a ghost elt belongs to partition p
     */
    bool isSubEntityOfGhostElement( rank_type p ) const
    {
        return M_elistGhost.find( p ) != M_elistGhost.end();
    }

    void clearElementsGhost()
    {
        M_elistGhost.clear();
    }

    //@}

    /**
     * set the tags associated to the points
     * - tags[0] physical region
     * - tags[1] elementary region
     * - tags[2] particular region
     */
    virtual void setTags( std::vector<int> const& tags )
    {
        M_markers[1].assign( tags[0] );
        if ( tags.size() > 1 )
            M_markers[2].assign( tags[1] );

        if ( tags.size() > 2 )
        {
            this->setProcessId( tags[3] );

            if ( tags[2] > 1 )
            {
                // ghosts
                std::vector<rank_type> p( tags[2] - 1 );

                for ( size_type i = 0; i < p.size(); ++i )
                {
                    p[i] = tags[4 + i];
                }

                this->setNeighborPartitionIds( p );
            }
        }
    }

    //! return all markers
    std::map<uint16_type, marker_type> const&
    markers() const
    {
        return M_markers;
    }
    //! set all markers
    void setMarkers( std::map<uint16_type, marker_type> const& markers )
    {
        M_markers = markers;
    }

    //! return true if has marker type id k
    bool hasMarkerType( uint16_type k ) const
    {
        auto itFindMarker = M_markers.find( k );
        if ( itFindMarker == M_markers.end() )
            return false;
        if ( itFindMarker->second.isOff() )
            return false;
        return true;
    }
    //! return true if has marker type id k
    bool hasMarker( uint16_type k ) const { return this->hasMarkerType( k ); }

    //! return marker with type id k
    marker_type const& marker( uint16_type k ) const
    {
        DCHECK( this->hasMarkerType( k ) ) << "no marker type " << k;
        return M_markers.find( k )->second;
    }
    //! set marker type id k with flag v
    void setMarker( uint16_type k, flag_type v )
    {
        M_markers[k].assign( v );
    }
    //! add flag v to the marker type id k
    void addMarker( uint16_type k, flag_type v )
    {
        M_markers[k].insert( v );
    }
    //! set marker type id k with flag v
    template <typename TT,std::enable_if_t< is_iterable_of_v<TT,flag_type>, bool> = true >
    void setMarker( uint16_type k,TT const& vs )
    {
        M_markers[k].assign( vs );
    }
    //! add flags vs to the marker type id k
    template <typename TT,std::enable_if_t< is_iterable_of_v<TT,flag_type>, bool> = true >
    void addMarker( uint16_type k, TT const& vs )
    {
        for ( auto const& v : vs )
            this->addMarker( k, v );
    }

    // methods for marker type id  1
    bool hasMarker() const
    {
        return this->hasMarkerType( 1 );
    }
    marker_type const& marker() const
    {
        return this->marker( 1 );
    }
    void setMarker( flag_type v )
    {
        this->setMarker( 1, v );
    }
    void addMarker( flag_type v )
    {
        this->addMarker( 1, v );
    }
    template <typename TT,std::enable_if_t< is_iterable_of_v<TT,flag_type>, bool> = true >
    void setMarker( TT const& vs )
    {
        this->setMarker( 1, vs );
    }
    template <typename TT,std::enable_if_t< is_iterable_of_v<TT,flag_type>, bool> = true >
    void addMarker( TT const& vs )
    {
        this->addMarker( 1, vs );
    }

    // methods for marker type id  2
    bool hasMarker2() const
    {
        return this->hasMarkerType( 2 );
    }
    marker_type const& marker2() const
    {
        return this->marker( 2 );
    }
    void setMarker2( flag_type v )
    {
        this->setMarker( 2, v );
    }

    // methods for marker type id  3
    bool hasMarker3() const
    {
        return this->hasMarkerType( 3 );
    }
    marker_type const& marker3() const
    {
        return this->marker( 3 );
    }
    void setMarker3( flag_type v )
    {
        this->setMarker( 3, v );
    }


protected:

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            DVLOG(2) << "Serializing GeoEntity...\n";
            DVLOG(2) << "  - id...\n";
            ar & M_id;
            DVLOG(2) << "  - id:" << M_id << "\n";
            DVLOG(2) << "  - entity...\n";
            ar & M_entity;
            DVLOG(2) << "  - entity:" << M_entity.context() << "\n";
            DVLOG(2) << "  - geometry...\n";
            ar & M_pid;
            DVLOG(2) << "  - pid:" << M_pid << "\n";
            ar & M_pidInPartition;
            ar & M_neighor_pids;
            ar & M_idInOtherPartitions;
            DVLOG( 2 ) << "  - markers...\n";
            ar& M_markers;
        }

private:


    size_type M_id;

    //! 2 bits GeoEntityContext, 6 bits ReferenceGeometry, 13 bits ReferenceShapes
    meta::Context<uint32_type> M_entity;

    //! maximum dimension of the entity touching the boundary within the element
    uint16_type M_boundaryEntityDimension;

    rank_type M_pid;
    rank_type M_pidInPartition;
    std::vector<rank_type> M_neighor_pids;
    std::map<rank_type, size_type> M_idInOtherPartitions;

    //! element list to which the point belongs
    std::set<std::pair<size_type,uint16_type>>  M_elist;
    //! ghost elements which share the entity
    std::map<rank_type,std::set<size_type > > M_elistGhost;

    //! mapping from marker index to marker flag
    std::map<uint16_type, marker_type> M_markers;
};

typedef GeoEntity<Simplex<0, 1> > GeoPoint;

// simplices
typedef GeoEntity<Simplex<1, 1> > LinearLine;
typedef GeoEntity<Simplex<2, 1> > LinearTriangle;
typedef GeoEntity<Simplex<3, 1> > LinearTetra;
typedef GeoEntity<Simplex<1, 2> > QuadraticLine;
typedef GeoEntity<Simplex<2, 2> > QuadraticTriangle;
typedef GeoEntity<Simplex<3, 2> > QuadraticTetra;

// tensor products
typedef GeoEntity<Hypercube<2, 1> > LinearQuad;
typedef GeoEntity<Hypercube<3, 1> > LinearHexa;
typedef GeoEntity<Hypercube<2, 2> > QuadraticQuad;
typedef GeoEntity<Hypercube<3, 2> > QuadraticHexa;

} // Feel

#endif /* __GeoEntity_H */
