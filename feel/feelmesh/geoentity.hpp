/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
    static const bool is_simplex_product = super::is_simplex_product;

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
        _M_id( 0 ),
        _M_entity( MESH_ENTITY_INTERNAL ),
        _M_geometry( Geometry ),
        _M_shape( Shape ),
        _M_pid( 0 )
        {}

    explicit GeoEntity( size_type i,
                        size_type geometry = Geometry,
                        size_type shape = Shape,
                        size_type context = MESH_ENTITY_INTERNAL )
        :
        super(),
        _M_id( i ),
        _M_entity( context ),
        _M_geometry( geometry ),
        _M_shape( shape ),
        _M_pid( 0 )
        {}

    GeoEntity( GeoEntity const& __me )
        :
        super(),
        _M_id( __me._M_id ),
        _M_entity( __me._M_entity ),
        _M_geometry( __me._M_geometry ),
        _M_shape( __me._M_shape ),
        _M_pid( __me._M_pid )
        {}

    GeoEntity& operator=( GeoEntity const& __me )
        {
            if ( this != &__me )
            {
                _M_id = __me._M_id;
                _M_entity = __me._M_entity;
                _M_geometry = __me._M_geometry;
                _M_shape = __me._M_shape;
                _M_pid = __me._M_pid;
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
            return _M_id == e.id();
        };
    bool operator<( GeoEntity const& e ) const
        {
            return _M_id < e.id();
        };

    bool operator<( size_type __i ) const
        {
            return _M_id < __i;
        };

    //@}

    /** @name Accessors
     */
    //@{

    size_type id() const
        {
            return _M_id;
        }


    /**
     * the dimension of the reference shape
     *
     *
     * @return the dimension of the reference shape
     */
    uint16_type refDim() const { return super::nDim; }

    /**
     * number of points on the reference shape
     *
     * @return the number of points on the reference shape
     */
    uint16_type nPoints() const { return super::numPoints; }

    /**
     * number of vertices on the reference shape
     *
     * @return the number of vertices on the reference shape
     */
    uint16_type nVertices() const { return super::numVertices; }

    /**
     * number of edges on the reference shape
     *
     * @return the number of edges on the reference shape
     */
    uint16_type nEdges() const { return super::numEdges; }

    /**
     * number of faces on the reference shape
     *
     * @return the number of edges on the reference shape
     */
    uint16_type nFaces() const { return super::numFaces; }

    /**
     * number of faces on the reference shape
     *
     * @return the number of edges on the reference shape
     */
    uint16_type nGeometricFaces() const { return super::numGeometricFaces; }

    /**
     * number of normals on the reference shape
     *
     * @return the number of normals on the reference shape
     */
    uint16_type nNormals() const { return super::numNormals; }


    /**
     *
     *
     *
     * @return true if the entoty has the shape \c __shape, false otherwise
     */
    bool hasShape( size_type __shape ) const { return _M_shape.test( __shape ); }

    /**
     * @return true of the entity is a volume
     */
    bool isAVolume() const { return _M_geometry.test( GEOMETRY_VOLUME ); }

    /**
     * @return true of the entity is a surface
     */
    bool isASurface() const { return _M_geometry.test( GEOMETRY_SURFACE ); }

    /**
     * @return true of the entity is a line
     */
    bool isALine() const { return _M_geometry.test( GEOMETRY_LINE ); }

    /**
     * @return true of the entity is a point
     */
    bool isAPoint() const { return _M_geometry.test( GEOMETRY_POINT ); }

    /**
     * @return true of the entity is a shape point
     */
    bool isAPointShape() const { return _M_shape.test( SHAPE_POINT ); }

    /**
     * @return true of the entity is a shape line
     */
    bool isALineShape() const { return _M_shape.test( SHAPE_LINE ); }

    /**
     * @return true of the entity is a triangle shape
     */
    bool isATriangleShape() const { return _M_shape.test( SHAPE_TRIANGLE ); }

    /**
     * @return true of the entity is a quadrangle
     */
    bool isAQuadrangleShape() const { return _M_shape.test( SHAPE_QUAD ); }

    /**
     * @return true of the entity is a tetrahedra shape
     */
    bool isATetrahedraShape() const { return _M_shape.test( SHAPE_TETRA ); }

    /**
     * @return true of the entity is a hexahedra
     */
    bool isAHexahedraShape() const { return _M_shape.test( SHAPE_HEXA ); }

    /**
     * @return true if the shape is linear, false otherwise
     */
    bool isLinear() const { return _M_shape.test( SHAPE_LINEAR ); }

    /**
     * @return true if the shape is bilinear, false otherwise
     */
    bool isBilinear() const { return _M_shape.test( SHAPE_BILINEAR ); }

    /**
     * @return true if the shape is quadratic, false otherwise
     */
    bool isQuadratic() const { return _M_shape.test( SHAPE_QUADRATIC ); }

    /**
     * @return true if the entity is internal, false otherwise
     */
    bool isInternal() const { return _M_entity.test( MESH_ENTITY_INTERNAL ); }


    /**
     * Tells if  item is on the boundary
     * @return true if on boundary, false otherwise
     */
    bool isOnBoundary() const
        {
            return _M_entity.test( MESH_ENTITY_BOUNDARY );
        };

    /**
     * \return the processor id of the entity
     */
    uint16_type processId() const { return _M_pid; }

    /**
     * set the processor id of the entity
     & \param pid processor id
     */
    void setProcessId( uint16_type pid )  { _M_pid = pid ; }

    /**
     * \return \c true if active, \c false otherwise
     *
     * \note for now it is a dummy function that returns always true,
     * will change when work on AMR starts
     */
    bool active() const { return true; }

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
            _M_id = id;
        }

    /**
     * set the boundary flag
     * @param b true if the item is on the boundary, false otherwise
     */
    void setOnBoundary( bool b )
        {
            if ( b )
            {
                _M_entity.set( MESH_ENTITY_BOUNDARY );
                _M_entity.clear( MESH_ENTITY_INTERNAL );
            }
            else
            {
                _M_entity.clear( MESH_ENTITY_BOUNDARY );
                _M_entity.set( MESH_ENTITY_INTERNAL );
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

    //@}



protected:

private:

    size_type _M_id;

    Context _M_entity;
    Context _M_geometry;
    Context _M_shape;

    uint16_type _M_pid;
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
