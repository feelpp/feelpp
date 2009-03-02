/*
  This file is part of the Life library

  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano
  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
  \file geond.hpp
*/
#ifndef _GEOND_HH_
#define _GEOND_HH_

#include <boost/numeric/ublas/storage.hpp>

#include <life/lifecore/application.hpp>

#include <life/lifemesh/geoentity.hpp>
#include <life/lifemesh/geo0d.hpp>
#include <life/lifepoly/geomap.hpp>
#include <life/lifemesh/marker.hpp>
#include <life/lifemesh/meshutil.hpp>
#include <life/lifemesh/meshbase.hpp>

namespace Life
{
class MeshBase;

template<int Dim, int Order,  template<uint16_type,uint16_type,uint16_type> class Entity, typename T> struct GT_Lagrange;

/// \cond detail
namespace detail
{
/**
 * \class ReversePoint
 * An utility to  invert point numbering on a GeoShape
 */
template <typename GeoShape>
class ReversePoint
{
public:
    uint16_type operate( uint16_type const & point )
        {
            return ( point < GeoShape::numVertices ) ?
                GeoShape::numVertices - point :
                GeoShape::numPoints - point + GeoShape::numVertices;
        }

};
}
/// \endcond


/**
 * @class GeoND
 * @brief Base class for Multi-dimensional basis Geometrical Entities.
 *
 */
template <uint16_type Dim,
          typename GEOSHAPE,
          typename T = double,
          typename POINTTYPE = Geo0D<Dim, T> >
class GeoND
    :
    public GeoEntity<GEOSHAPE>
{
    typedef GeoEntity<GEOSHAPE> super;
public:

    typedef T value_type;
    typedef GeoND<Dim, GEOSHAPE, value_type, POINTTYPE> self_type;
    typedef self_type element_type;

    typedef GEOSHAPE GeoShape;
    typedef POINTTYPE PointType;

    typedef PointType point_type;
    typedef typename super::face_type face_type;

    static const size_type Shape = super::Shape;
    static const uint16_type numPoints = super::numPoints;
    static const uint16_type numVertices = super::numVertices;
    static const uint16_type numLocalPoints = super::numPoints;
    static const uint16_type numLocalVertices = super::numVertices;
    static const int numFaces = super::numFaces;
    static const int numEdges = super::numEdges;
    static const int numTopologicalFaces = super::numTopologicalFaces;
    static const uint16_type numNeighbors = super::numTopologicalFaces;


    typedef typename ublas::bounded_array<point_type*, numPoints>::iterator point_iterator;
    typedef typename ublas::bounded_array<point_type*, numPoints>::const_iterator point_const_iterator;

    typedef typename matrix_node<value_type>::type matrix_node_type;
    typedef typename node<value_type>::type node_type;

    static const uint16_type nDim = super::nDim;
    static const uint16_type nOrder = super::nOrder;
    static const uint16_type nRealDim = super::nRealDim;

    typedef typename mpl::if_<mpl::bool_<GeoShape::is_simplex_product>,
                              mpl::identity<GT_Lagrange<nDim, nOrder, SimplexProduct, T> >,
                              mpl::identity<GT_Lagrange<nDim, nOrder, Simplex, T> > >::type::type gm_type;
    typedef boost::shared_ptr<gm_type> gm_ptrtype;

    /**
     * default constructor
     */
    GeoND()
        :
        super( 0 ),
        _M_points( numPoints ),
        _M_face_points( numTopologicalFaces ),
        _M_G( nRealDim, numPoints ),
        _M_has_points( false ),
        _M_neighbors( numNeighbors, std::make_pair( invalid_size_type_value, 0 ) ),
        _M_marker1(),
        _M_marker2(),
        _M_marker3()
    {
    }

    /**
     * constructor from an id
     *
     * @param id identifier for the element to store
     *
     */
    explicit GeoND( size_type id )
        :
        super( id ),
        _M_points( numPoints ),
        _M_face_points( numTopologicalFaces ),
        _M_G( nRealDim, numPoints ),
        _M_has_points( false ),
        _M_neighbors( numNeighbors, std::make_pair( invalid_size_type_value, 0 ) ),
        _M_marker1(),
        _M_marker2(),
        _M_marker3()
    {
    }

    /**
     * copy constructor
     */
    GeoND( GeoND const& G )
        :
        super( G ),
        _M_points( G._M_points ),
        _M_face_points( G._M_face_points ),
        _M_G( G._M_G ),
        _M_has_points( G._M_has_points ),
        _M_neighbors( G._M_neighbors ),
        _M_marker1( G._M_marker1 ),
        _M_marker2( G._M_marker2 ),
        _M_marker3( G._M_marker3 )
    {
#if 0
        for ( uint16_type i = 0; i < numLocalPoints; ++i )
            _M_points[ i ] = G._M_points[ i ];

        _M_G = G._M_G;
#endif
    }

    /**
     * destructor, make it virtual for derived classes
     */
    virtual ~GeoND()
    {
    }

    /**
     * set the mesh to which this geometric entity belongs to
     */
    void setMeshAndGm( MeshBase const* m, gm_ptrtype const& gm ) const
    {
        M_mesh = m;
        M_gm = gm;
    }

    void setMesh( MeshBase const* m ) const
    {
        M_mesh = m;
    }

    //! return the geometric mapping if a mesh was set
    gm_ptrtype gm() const { return M_gm; }

    /**
     * \return the mesh to which this geometric entity belongs to
     */
    MeshBase const* mesh() const { return M_mesh; }

    /**
     * \return true if points have been inserted in elements, false
     * otherwise
     */
    bool hasPoints() const { return _M_has_points; }

    /**
     * assignment operator
     *
     * @param G the object to assign
     *
     * @return the object that was assigned
     */
    GeoND & operator=( GeoND const & G )
    {
        if ( this != &G )
        {
            super::operator=( G );
            for ( uint16_type i = 0; i < numLocalPoints; ++i )
                _M_points[ i ] = G._M_points[ i ];

            _M_face_points = G._M_face_points;
            _M_G = G._M_G;

            _M_has_points = G._M_has_points;

            _M_neighbors = G._M_neighbors;

            _M_marker1 = G._M_marker1;
            _M_marker2 = G._M_marker2;
            _M_marker3 = G._M_marker3;
        }
        return *this;
    }

    /**
     * \return the number of points in convex
     */
    uint16_type nPoints() const { return numPoints; }

    /**
     * the number of neighbors is equal to the number of
     * faces. Neighbors are stored as pointers and some of them can be
     * null pointers if the corresponding face is on the boundary.
     *
     *\return the number of neighbors
     */
    uint16_type nNeighbors() const { return numNeighbors; }

    /**
     * Neighbors are stored as pointers and some of them can be null
     * pointers if the corresponding face is on the boundary.
     *
     * \return the pair neighbor \p n index and process \p id it belongs to
     */
    std::pair<size_type,size_type> const& neighbor( uint16_type n ) const { return _M_neighbors[n]; }

    /**
     * set the \p n -th neighbor with \p neigh
     */
    void setNeighbor( uint16_type n, size_type neigh_id, size_type proc_id = Application::processId() )
    {
        _M_neighbors[n] = std::make_pair( neigh_id, proc_id );
    }

    /**
     * It returns the reference to an point object (possibly derived from
     * Geo0D)
     */
    PointType & point( uint16_type i )
    {
        return *( static_cast<POINTTYPE *>( _M_points[ i ] ) );
    }

    /**
     * It returns the reference to an point object (possibly derived from
     * Geo0D)
     */
    PointType const & point ( uint16_type i ) const
    {
        return *( static_cast<POINTTYPE *>( _M_points[ i ] ) );
    }
    /**
     * It returns the pointer to an point object (possibly derived from
     * Geo0D)
     */
    PointType* pointPtr( uint16_type i )
    {
        return _M_points[ i ];
    }

    /**
     * It returns the pointer to an point object (possibly derived from
     * Geo0D)
     */
    PointType const* pointPtr ( uint16_type i ) const
    {
        return _M_points[ i ];
    }

    /**
     */
    PointType const & facePoint ( uint16_type __f, uint16_type const __i ) const { return _M_face_points[__f][__i]; }


    /**
     *  The ith point (starting from the end)
     *
     *  It returns the reference to an point object (possibly derived from
     *  Geo0D). It starts from the last point, yet it follows the rule: vertices
     *  first. It may be used to access the points of a Geometry Element in a
     *  reverse way (i.e. with the opposite GeoElement permutation)
     */
    PointType & reversepoint( uint16_type const i )
    {
        return *( static_cast<POINTTYPE *>( _M_points[ detail::ReversePoint<GEOSHAPE>::operate( i ) ] ) );
    }


    /**
     *  The ith point (starting from the end)
     *
     *  It returns the reference to an point object (possibly derived from
     *  Geo0D). It starts from the last point, yet it follows the rule: vertices
     *  first. It may be used to access the points of a Geometry Element in a
     *  reverse way (i.e. with the opposite GeoElement permutation)
     */
    PointType const & reversepoint ( uint16_type const i ) const
    {
        return *( static_cast<POINTTYPE *>( _M_points[ detail::ReversePoint<GEOSHAPE>::operate( i ) ] ) );
    }


    /**
     * Inserts a point.  Uses point references
     * put point
     */
    void setPoint( uint16_type const i, point_type const & p );

#if 0
    /**
     * Inserts a point. Uses pointers
     * put point
     */
    void setPoint( uint16_type const i, point_type const * p );

    /**
     * Inserts a point Uses point references (bounds check)
     * with forced bound check
     */
    bool setPointBD( uint16_type const i, point_type const & p );


    /**
     * Inserts a point. Uses pointers (bounds check)
     * with forced bound check
     */
    bool setPointBD( uint16_type const i, point_type const * p );
#endif
    /**
     * show information about the geoND
     *
     * @param verbose true if verbose mode, false otherwise
     * @param c output stream
     *
     * @return the output stream
     */
    std::ostream & showMe( bool verbose = false, std::ostream & c = std::cout ) const;

    /**
     * Swap Points
     *
     *  This is a member function to be used ONLY by routines for checking or
     *  amending meshes. You must give the local id
     *
     * @param pt1 1st point to swap with 2nd point
     * @param pt2 2nd point to swap with 1st point
     */
    void swapPoints( const uint16_type & pt1, const uint16_type & pt2 );

    /**
     *  Exchange Points
     *
     *  Exchanges points according to a list of old2new local id numbering !
     *  old2new[i] is the new local id of a point whose old local id was ! i+1
     *  (remeber the numbering from 1 of the id's!. This is a member function
     *  to be used ONLY by routines for checking or amending meshes. You must
     *  give uint16_type (which start
     */
    void exchangePoints( const uint16_type otn[ numPoints ] );

    /**
     * matrix of geometric nodes
     * retrieve the matrix of geometric nodes (Dim x NumPoints) the
     * matrix is column oriented, the column i contains the coordinate
     * of the i-th geometric node of the element
     *
     * \return the matrix of geometric nodes
     */
    matrix_node_type const& G() const { return _M_G; }

    /**
     * matrix of geometric nodes
     * retrieve the matrix of geometric nodes (Dim x NumPoints) the
     * matrix is column oriented, the column i contains the coordinate
     * of the i-th geometric node of the element
     *
     * \return the matrix of geometric nodes
     */
    matrix_node_type & G() { return _M_G; }

    point_iterator beginPoint() { return _M_points.begin(); }
    point_const_iterator beginPoint() const { return _M_points.begin(); }
    point_iterator endPoint() { return _M_points.end(); }
    point_const_iterator endPoint() const { return _M_points.end(); }


    struct tt
    {
        static uint16_type fToP( uint16_type const _localFace, uint16_type const _point )
        {
            return super::eToP( _localFace, _point );
        }
    };

    /**
     * Get the local id of the point in the element
     *
     * @param _localFace local id of a face in the element
     * @param _point local id of a point in the face
     *
     * @return the local id of the point in the element
     */
    static uint16_type fToP( uint16_type const _localFace, uint16_type const _point )
    {
#if 1
        typedef typename mpl::if_<mpl::greater<mpl::int_<super::nDim>, mpl::int_<2> >,
            mpl::identity<super>,
            mpl::identity<tt> >::type the_type;
        return the_type::type::fToP( _localFace, _point );
#else
        return super::fToP( _localFace, _point );
#endif

    }


    /**
     * get the number of opposite points per face
     *
     * @return the number of opposite points per face
     */
    uint16_type nOppositePointsPerFace() const
    {
        return super::nbOppositePointsPerFace;
    }

    /**
     * faceToOppositePoint(i,j) = localId of jth opposite point to ith local face
     *
     * @return the localId of _point-th opposite point to _localFace-th local face
     */
    uint16_type faceToOppositePoint( uint16_type const _localFace, uint16_type const _point ) const
    {
        return super::faceToOppositePoint( _localFace, _point );
    }

    /**
     * Determines if the local numbering of a 2D element
     * is oriented anticlockwise
     */
    bool isAnticlockwiseOriented() const
    {

        // Calculate vectors originating from vertex zero
        ublas::matrix<T> orientation_matrix (nRealDim,nRealDim);

        for (int i = 0; i < nRealDim ; ++i)
            {
                ublas::row( orientation_matrix, i ) = ( ublas::column( this->G(), i+1 ) -
                                                        ublas::column( this->G(),   0 ) );

            }
        LU< ublas::matrix<T> > lu(orientation_matrix);
        T sgn=lu.det();

        return (sgn > 0) ? 1 : 0;
    }

    void applyDisplacement( int i, ublas::vector<double> const& u )
    {
        ublas::column( _M_G, i ) += u;
        (*_M_points[ i ]) += u;
    }
    void applyDisplacementG( int i, ublas::vector<double> const& u )
    {
        ublas::column( _M_G, i ) += u;
    }
    Marker1 const& marker() const { return _M_marker1; }
    Marker1& marker() { return _M_marker1; }
    void setMarker( flag_type v ) { return _M_marker1.assign( v ); }

    Marker2 const& marker2() const { return _M_marker2; }
    Marker2& marker2() { return _M_marker2; }
    void setMarker2( flag_type v ) { return _M_marker2.assign( v ); }

    Marker3 const& marker3() const { return _M_marker3; }
    Marker3& marker3() { return _M_marker3; }
    void setMarker3( flag_type v ) { return _M_marker3.assign( v ); }

private:

    /** geometric nodes of the element */
    ublas::bounded_array<point_type*, numPoints> _M_points;

    /** geometric nodes of the faces of the element */
    std::vector<ublas::bounded_array<point_type*, numPoints> > _M_face_points;

    /**< matrix of the geometric nodes */
    matrix_node_type _M_G;

    bool _M_has_points;

    /**
     * store neighbor element id
     */
    std::vector<std::pair<size_type,size_type> > _M_neighbors;


    Marker1 _M_marker1;
    Marker2 _M_marker2;
    Marker3 _M_marker3;

    // mesh to which the geond element belongs to
    mutable MeshBase const* M_mesh;
    mutable gm_ptrtype M_gm;
};

template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
const uint16_type GeoND<Dim,GEOSHAPE, T, POINTTYPE>::numLocalPoints;

template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
const uint16_type GeoND<Dim,GEOSHAPE, T, POINTTYPE>::numLocalVertices;

template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
inline
void
GeoND<Dim,GEOSHAPE, T, POINTTYPE>::setPoint( uint16_type const i, point_type const & p )
{
    _M_points[ i ] = const_cast<point_type *>( &p );
    //Debug() << "[setPoint] üpdate point index " << i << " with "<< _M_points[i]->id() << "\n";
    LIFE_ASSERT( const_cast<point_type *>( &p ) != 0 ).error( "invalid Geo0D<>" );
    ublas::column( _M_G, i ) = _M_points[i]->node();
    _M_has_points = true;
    //Debug() << "[setPoint] üpdate point index " << i << " with "<< _M_points[i]->id() << "\n";
}

#if 0
template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
bool GeoND<Dim,GEOSHAPE, T, POINTTYPE>::setPointBD( uint16_type const i, point_type const & p )
{
    // if not assert we need anyway to avoid under/overflows
    if ( i > numLocalVertices )
        return false;

    _M_points[ i ] = const_cast<point_type *>( &p );
    ublas::column( _M_G, i ) = _M_points[i]->node();
    _M_has_points = true;
    return true;
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
inline
void GeoND<Dim,GEOSHAPE, T, POINTTYPE>::setPoint( uint16_type const i, point_type const * p )
{
    _M_points[ i ] = const_cast<point_type *>( p );
    ublas::column( _M_G, i ) = _M_points[i]->node();
    _M_has_points = true;
}


template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
bool GeoND<Dim,GEOSHAPE, T, POINTTYPE>::setPointBD( uint16_type const i, point_type const * p )
{
    // if not assert we need anyway to avoid under/overflows
    if ( i > numLocalVertices )
        return false;

    _M_points[ i ] = const_cast<point_type *>( p );
    ublas::column( _M_G, i ) = _M_points[i]->node();
    _M_has_points = true;
    return true;
}
#endif

template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
std::ostream &
GeoND<Dim,GEOSHAPE, T, POINTTYPE>::showMe( bool verbose, std::ostream & out ) const
{
    out << "----- BEGIN OF GeoND data ---" << std::endl << std::endl;
    out << " GeoND object of shape " << Shape << std::endl;
    out << " Number of Vertices = " << numVertices << std::endl;
    out << "   Number of Points = " << numPoints << std::endl;
    out << "                 id = " << this->id() << std::endl;
    out << "                  G = " << _M_G << "\n";
    for ( int i = 0; i < numVertices; i++ )
    {
        out << "POINT id = " << i << std::endl;
        out << point( i ).showMe( verbose, out );
    }
    out << "----- END OF GeoND data ---" << std::endl << std::endl;
    return out;
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
void GeoND<Dim,GEOSHAPE, T, POINTTYPE>::swapPoints( const uint16_type & pt1, const uint16_type & pt2 )
{
    point_type * tmp( _M_points[ pt1 ] );
    _M_points[ pt1 ] = _M_points[ pt2 ];
    _M_points[ pt2 ] = tmp;

    // swap also the entries in G
    ublas::column( _M_G, pt1 ).swap( ublas::column( _M_G, pt2 ) );
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
void GeoND<Dim,GEOSHAPE, T, POINTTYPE>::exchangePoints( const uint16_type otn[ numPoints ] )
{
    point_type * tmp[ numPoints ];
    for ( unsigned int i = 0;i < numPoints;++i )
    {
        tmp[ i ] = _M_points[ i ];
    }
    for ( unsigned int i = 0;i < numPoints;++i )
    {
        _M_points[ i ] = tmp[ otn[ i ] ];
        ublas::column( _M_G, i ) = _M_points[i]->node();
    }
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
inline
DebugStream&
operator<<( DebugStream& __os, GeoND<Dim,GEOSHAPE, T, POINTTYPE> const& __n )
{
    if ( __os.doPrint() )
    {
        std::ostringstream __str;

        __str << __n.showMe( true, __str );

        __os << __str.str() << "\n";
    }
    return __os;
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename POINTTYPE>
inline
NdebugStream&
operator<<( NdebugStream& __os, GeoND<Dim,GEOSHAPE, T, POINTTYPE> const& __n )
{
    return __os;
}


} // Life
#endif
