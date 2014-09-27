/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-08-10

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
template<typename Entity>
class GeoEntity
    :
    boost::equality_comparable<GeoEntity<Entity> >,
    boost::less_than_comparable<GeoEntity<Entity> >,
    boost::less_than_comparable<GeoEntity<Entity>, size_type>,
public Entity

{
public:


    /** @name Typedefs
     */
    //@{

    typedef Entity super;
    typedef GeoEntity<Entity> GeoShape;
    typedef GeoEntity<Entity> self_type;
    typedef typename super::topological_face_type face_type;
    typedef face_type GeoBShape;
    typedef typename Entity::edge_permutation_type edge_permutation_type;
    typedef typename Entity::face_permutation_type face_permutation_type;

    static const size_type Shape = super::Shape;
    static const size_type Geometry = super::Geometry;

    static const uint16_type nDim = super::nDim;
    static const uint16_type nOrder = super::nOrder;
    static const uint16_type nRealDim = super::nRealDim;


    static const uint16_type numVertices = super::numVertices;
    static const uint16_type numFaces = super::numFaces;
    static const uint16_type numGeometricFaces = super::numGeometricFaces;
    static const uint16_type numTopologicalFaces = super::numTopologicalFaces;
    static const uint16_type numEdges = super::numEdges;
    static const uint16_type numNormals = super::numNormals;

    static const uint16_type numPoints = super::numPoints;
    static const uint16_type nbPtsPerVertex = super::nbPtsPerVertex;
    static const uint16_type nbPtsPerEdge = super::nbPtsPerEdge;
    static const uint16_type nbPtsPerFace = super::nbPtsPerFace;
    static const uint16_type nbPtsPerVolume = super::nbPtsPerVolume;

    typedef Entity convex_type;

    static const bool is_simplex = super::is_simplex;
    static const bool is_hypercube = super::is_hypercube;

    /**
     * helper class to construct the associated reference convex.
     */
    template<typename T = double>
    struct reference_convex
    {
        typedef Reference<Entity, nDim, nOrder, nRealDim, T> type;
    };
    //@}

    /** @name Constructors, destructor
     */
    //@{

    GeoEntity()
        :
        super(),
        M_id( 0 ),
        M_entity( MESH_ENTITY_INTERNAL ),
        M_geometry( Geometry ),
        M_shape( Shape ),
        M_boundaryEntityDimension( invalid_uint16_type_value ),
        M_npids( 1 ),
        M_pid( invalid_rank_type_value ),
        M_pidInPartition( invalid_rank_type_value ),
        M_neighor_pids(),
        M_idInOthersPartitions(),
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
        M_entity( context ),
        M_geometry( geometry ),
        M_shape( shape ),
        M_boundaryEntityDimension( invalid_uint16_type_value ),
        M_npids( 1 ),
        M_pid( invalid_rank_type_value ),
        M_pidInPartition( invalid_rank_type_value ),
        M_neighor_pids(),
        M_idInOthersPartitions(),
        M_elist(),
        M_elistGhost()
    {}

    GeoEntity( GeoEntity const& __me )
        :
        super(),
        M_id( __me.M_id ),
        M_entity( __me.M_entity ),
        M_geometry( __me.M_geometry ),
        M_shape( __me.M_shape ),
        M_boundaryEntityDimension( __me.M_boundaryEntityDimension ),
        M_npids( __me.M_npids ),
        M_pid( __me.M_pid ),
        M_pidInPartition( __me.M_pidInPartition ),
        M_neighor_pids( __me.M_neighor_pids ),
        M_idInOthersPartitions( __me.M_idInOthersPartitions ),
        M_elist( __me.M_elist ),
        M_elistGhost( __me.M_elistGhost )
    {}

    GeoEntity& operator=( GeoEntity const& __me )
    {
        if ( this != &__me )
        {
            M_id = __me.M_id;
            M_entity = __me.M_entity;
            M_geometry = __me.M_geometry;
            M_shape = __me.M_shape;
            M_boundaryEntityDimension = __me.M_boundaryEntityDimension;
            M_npids = __me.M_npids;
            M_pid = __me.M_pid;
            M_pidInPartition = __me.M_pidInPartition;
            M_neighor_pids = __me.M_neighor_pids;
            M_idInOthersPartitions = __me.M_idInOthersPartitions;
            M_elist = __me.M_elist;
            M_elistGhost = __me.M_elistGhost;
        }

        return *this;
    }

    virtual ~GeoEntity()
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

    size_type id() const
    {
        return M_id;
    }


    /**
     * the dimension of the reference shape
     *
     *
     * @return the dimension of the reference shape
     */
    uint16_type refDim() const
    {
        return super::nDim;
    }

    /**
     * number of points on the reference shape
     *
     * @return the number of points on the reference shape
     */
    uint16_type nPoints() const
    {
        return super::numPoints;
    }

    /**
     * number of vertices on the reference shape
     *
     * @return the number of vertices on the reference shape
     */
    uint16_type nVertices() const
    {
        return super::numVertices;
    }

    /**
     * number of edges on the reference shape
     *
     * @return the number of edges on the reference shape
     */
    uint16_type nEdges() const
    {
        return super::numEdges;
    }

    /**
     * number of faces on the reference shape
     *
     * @return the number of edges on the reference shape
     */
    uint16_type nFaces() const
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
    uint16_type nGeometricFaces() const
    {
        return super::numGeometricFaces;
    }

    /**
     * number of normals on the reference shape
     *
     * @return the number of normals on the reference shape
     */
    uint16_type nNormals() const
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
        return M_shape.test( __shape );
    }

    /**
     * @return true of the entity is a volume
     */
    bool isAVolume() const
    {
        return M_geometry.test( GEOMETRY_VOLUME );
    }

    /**
     * @return true of the entity is a surface
     */
    bool isASurface() const
    {
        return M_geometry.test( GEOMETRY_SURFACE );
    }

    /**
     * @return true of the entity is a line
     */
    bool isALine() const
    {
        return M_geometry.test( GEOMETRY_LINE );
    }

    /**
     * @return true of the entity is a point
     */
    bool isAPoint() const
    {
        return M_geometry.test( GEOMETRY_POINT );
    }

    /**
     * @return true of the entity is a shape point
     */
    bool isAPointShape() const
    {
        return M_shape.test( SHAPE_POINT );
    }

    /**
     * @return true of the entity is a shape line
     */
    bool isALineShape() const
    {
        return M_shape.test( SHAPE_LINE );
    }

    /**
     * @return true of the entity is a triangle shape
     */
    bool isATriangleShape() const
    {
        return M_shape.test( SHAPE_TRIANGLE );
    }

    /**
     * @return true of the entity is a quadrangle
     */
    bool isAQuadrangleShape() const
    {
        return M_shape.test( SHAPE_QUAD );
    }

    /**
     * @return true of the entity is a tetrahedra shape
     */
    bool isATetrahedraShape() const
    {
        return M_shape.test( SHAPE_TETRA );
    }

    /**
     * @return true of the entity is a hexahedra
     */
    bool isAHexahedraShape() const
    {
        return M_shape.test( SHAPE_HEXA );
    }

    /**
     * @return true if the shape is linear, false otherwise
     */
    bool isLinear() const
    {
        return M_shape.test( SHAPE_LINEAR );
    }

    /**
     * @return true if the shape is bilinear, false otherwise
     */
    bool isBilinear() const
    {
        return M_shape.test( SHAPE_BILINEAR );
    }

    /**
     * @return true if the shape is quadratic, false otherwise
     */
    bool isQuadratic() const
    {
        return M_shape.test( SHAPE_QUADRATIC );
    }

    /**
     * @return true if the entity is internal, false otherwise
     */
    bool isInternal() const
    {
        return M_entity.test( MESH_ENTITY_INTERNAL );
    }


    /**
     * Tells if  item is on the boundary
     * @return true if on boundary, false otherwise
     */
    bool isOnBoundary() const
    {
        return M_entity.test( MESH_ENTITY_BOUNDARY );
    }

    /**
     * maximum dimension of the entity of the element touching the boundary
     */
    uint16_type boundaryEntityDimension() const
    {
        return M_boundaryEntityDimension;
    }
    /**
     * \return \c true if ghost cell, \c false otherwise
     */
    bool isGhostCell() const
    {
        //return (this->worldComm().localRank()!=M_pid);
        //mpi::communicator world;
        //return (world.rank()!=M_pid);
        return ( M_pidInPartition!=M_pid );
    }

    /**
     * \return the processor id of the entity
     */
    rank_type processId() const
    {
        return M_pid;
    }

    /**
     * set the processor id of the entity
     & \param pid processor id
     */
    void setProcessId( rank_type pid )
    {
        M_pid = pid ;
    }

    /**
     * \return the processor id of the entity
     */
    rank_type pidInPartition() const
    {
        return M_pidInPartition;
    }
    /**
     * set the processor id of the entity
     & \param pid processor id
     */
    void setProcessIdInPartition( rank_type pid )
    {
        M_pidInPartition = pid ;
    }

    /**
     * \return the partition id
     */
    rank_type partitionId() const
    {
        return M_pid;
    }

    /**
     * \return the number of partition the element is linked to including the
     * partition to which it belongs
     */
    rank_type numberOfPartitions() const
    {
        return M_npids;
    }

    /**
     * \return the number of partition the element is linked to
     */
    size_type numberOfNeighborPartitions() const
    {
        return M_neighor_pids.size();
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
    void setIdInOthersPartitions( rank_type pid, size_type id )
    {
        M_idInOthersPartitions.insert( std::make_pair( pid, id ) );
    }

    /**
     * \return the id of the entity in a partition pid
     */
    size_type idInOthersPartitions( rank_type pid ) const
    {
        DCHECK( M_idInOthersPartitions.find( pid )!=M_idInOthersPartitions.end() ) << " id is unknow for this pid " << pid << "\n";
        return M_idInOthersPartitions.find( pid )->second;
    }

    /**
     * \return idInOthersPartitions map
     */
    std::map<rank_type, size_type> const& idInOthersPartitions() const
    {
        return M_idInOthersPartitions;
    }
    /**
     * clear id in others partitions container
     */
    void clearIdInOthersPartitions()
    {
        M_idInOthersPartitions.clear();
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
    virtual double measure() const = 0;

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
            M_entity.set( MESH_ENTITY_BOUNDARY );
            M_entity.clear( MESH_ENTITY_INTERNAL );
        }
        else
        {
            M_entity.clear( MESH_ENTITY_BOUNDARY );
            M_entity.set( MESH_ENTITY_INTERNAL );
        }
        M_boundaryEntityDimension = ent_d;
    }

    /**
     * \return the number of partition the element is linked to including the
     * partition to which it belongs
     */
    void setNumberOfPartitions( uint16_type np )
    {
        M_npids = np;
    }

    /**
     * set the number of partition the element is linked to
     */
    void setNumberOfNeighborPartitions( uint16_type nep )
    {
        FEELPP_ASSERT( M_npids -1 == M_neighor_pids.size() )( M_npids )( M_neighor_pids ).error( "invalid partitioning data" );
        M_neighor_pids.size();
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
            ++M_npids;
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
    self_type& addElement( size_type e )
    {
        M_elist.insert( e );
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
    std::set<size_type> const& elements() const
    {
        return M_elist;
    }
    std::set<size_type>& elements()
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
            ar & M_geometry;
            DVLOG(2) << "  - geometry:" << M_geometry.context() << "\n";
            DVLOG(2) << "  - shape...\n";
            ar & M_shape;
            DVLOG(2) << "  - shape:" << M_shape.context() << "\n";
            DVLOG(2) << "  - npids...\n";
            ar & M_npids;
            DVLOG(2) << "  - npids:" << M_npids << "\n";
            DVLOG(2) << "  - pid...\n";
            ar & M_pid;
            DVLOG(2) << "  - pid:" << M_pid << "\n";
            ar & M_pidInPartition;
            ar & M_neighor_pids;
            ar & M_idInOthersPartitions;
        }

private:


    size_type M_id;

    Context M_entity;
    Context M_geometry;
    Context M_shape;

    //! maximum dimension of the entity touching the boundary within the element
    uint16_type M_boundaryEntityDimension;

    rank_type M_npids;
    rank_type M_pid;
    rank_type M_pidInPartition;
    std::vector<rank_type> M_neighor_pids;
    std::map<uint16_type, size_type> M_idInOthersPartitions;

    //! element list to which the point belongs
    std::set<size_type> M_elist;
    //! ghost elements which share the entity
    std::map<rank_type,std::set<size_type > > M_elistGhost;

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
