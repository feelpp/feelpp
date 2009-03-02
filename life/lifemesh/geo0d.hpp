/*
 This file is part of the Life library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

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
 * \file geo0D.hpp
 */
#ifndef _GEO0D_HH_
#define _GEO0D_HH_

#include <boost/operators.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifealg/glas.hpp>
#include <life/lifemesh/geoentity.hpp>
#include <life/lifemesh/meshbase.hpp>
#include <life/lifemesh/marker.hpp>

namespace Life
{
class MeshBase;

/**
 *  \defgroup GeoXD Basis Geometrical Entities Geo0D and GeoND.
 *
 *  They are intermediate classes used to build the actual Geometry classes
 *
 *  \warning Geo1D/2D/3D are template classes because some of the info is not
 *  known a priori and I want all vector dimensions determined at compile time
 *  to enhance memory access time.
 */

/*@{*/

/**
 * \class Geo0D
 *
 *  Zero dimensional entity.
 *
 * \ingroup GeoXD
 */
template<uint16_type Dim, typename T = double>
class Geo0D
    :
        public boost::equality_comparable<Geo0D<Dim,T> >,
        public boost::less_than_comparable<Geo0D<Dim,T> >,
        public boost::less_than_comparable<Geo0D<Dim,T>, size_type>,
        public GeoPoint
{
    typedef GeoPoint super;
public:

    typedef Geo0D<Dim,T> self_type;
    static const uint16_type nDim = Dim;
    typedef T value_type;
    typedef typename matrix_node<value_type>::type matrix_node_type;
    typedef typename node<T,  nDim>::type node_type;
    /**
     * default constructor
     *
     */
    Geo0D();

    /**
     * constructor where I give the id and declare if Geo0D object is on a
     * boundary
     *
     * @param id identifier of the Geo0D
     * @param boundary true if on the boundary, false otherwise
     * @param is_vertex true if the point is a vertex
     *
     * @return
     */
    explicit Geo0D( size_type id, bool boundary = false, bool is_vertex = false );

    /**
     * constructor where I give the id, the point coordinate and I declare
     * if the Geo0D object is on a boundary
     *
     * @param id identifier of the node
     * @param x x-coordinate of the node
     * @param y y-coordinate of the node
     * @param z z-coordinate of the node
     * @param boundary true if on the boundary, false otherwise
     * @param is_vertex true if the point is a vertex
     *
     */
    explicit Geo0D( size_type id, value_type x, value_type y, value_type z, bool boundary = false, bool is_vertex = false );

    /**
     * give the point coordinate
     *
     * @param x x-coordinate of the node
     * @param y y-coordinate of the node
     * @param z z-coordinate of the node
     */
    Geo0D( value_type x, value_type y, value_type z )
        :
        super( 0, MESH_ENTITY_INTERNAL ),
        _M_coor(Dim),
        _M_is_vertex( false )
    {
        if ( Dim < 2 )
            _M_coor[ 0 ] = x;
        if (  Dim < 3 )
            _M_coor[ 1 ] = y;
        if ( Dim == 3 )
            _M_coor[ 2 ] = z;
    }


    /**
     * constructor where I give the id, the point coordinate and I declare
     * if the Geo0D object is on a boundary
     *
     * @param id identifier
     * @param __x node coordinate
     * @param boundary true if on the boundary, false otherwise
     * @param is_vertex true if the point is a vertex
     *
     */
    Geo0D( size_type id, node_type const& __x, bool boundary = false, bool is_vertex = false );

    /**
     * the point coordinate
     *
     * @param __x node coordinate
     *
     */
    Geo0D( node_type const& __x )
        :
        super( 0, MESH_ENTITY_INTERNAL ),
        _M_coor(__x),
        _M_is_vertex( false )
    {
    }

    /**
     * the point coordinate expression
     *
     * @param __expr expression for the coordinates
     *
     */
    template<typename AE>
    Geo0D( ublas::vector_expression<AE> const& __expr )
        :
        super( 0, MESH_ENTITY_INTERNAL ),
        _M_coor(__expr),
        _M_is_vertex( false )
    {
    }

    /**
     * copy constructor
     * @param G the Geo0D to copy
     */
    Geo0D( Geo0D const & G );

    /**
     * assignement operator
     *
     * @param G the geo0D to copy
     *
     * @return the newly assigned Geo0D
     */
    Geo0D & operator=( Geo0D const & G );

    Geo0D & operator+=( node_type const & G )
    {
        _M_coor += G;
        return *this;
    }


    value_type& operator()( int i )
    {
        return _M_coor[i];
    }

    value_type  operator()( int i ) const
    {
        return _M_coor[i];
    }

    /**
     * \return \p true if point is a vertex, \p false otherwise
     * \attention DO NOT USE YET, returns always false
     */
    bool isVertex() const { return _M_is_vertex; }

    /**
     * set the point as a vertex or not using \p v
     */
    void  setAsVertex( bool v )
    {
        _M_is_vertex = v;
    }

    /**
     * set the mesh to which this geometric entity belongs to
     */
    void setMesh( MeshBase const* m )
    {
        M_mesh = m;
    }

    /**
     * \return the mesh to which this geometric entity belongs to
     */
    MeshBase const* mesh() const { return M_mesh; }

    /**
     * @return the \c begin() iterator of the coordinate container
     */
    typename node_type::iterator begin() { return _M_coor.begin(); }

    /**
     * @return the \c begin() iterator of the coordinate container
     */
    typename node_type::const_iterator begin() const { return _M_coor.begin(); }

    /**
     * @return the \c end() iterator of the coordinate container
     */
    typename node_type::iterator end() { return _M_coor.end(); }

    /**
     * @return the \c end() iterator of the coordinate container
     */
    typename node_type::const_iterator end() const { return _M_coor.end(); }

    /**
     * @return the node data structure
     */
    node_type const& node() const { return _M_coor; }

    /**
     * @return the node data structure
     */
    matrix_node_type G() const { matrix_node_type __G( Dim, 1 ); ublas::column( __G, 0 ) = _M_coor; return __G; }

    /**
     * set the node coordinates
     *
     * @param __n the node coordinates
     */
    void setNode( node_type const& __n )
    {
        _M_coor = __n;
        //Debug() << "point " << id() << " coords: " << _M_coor << "\n";
    }

    /**
     * \return true if points are equal, false otherwise
     */
    bool operator==( Geo0D const& geo0d ) const
    {
        return this->id() == geo0d.id();//this->isEqual( geo0d, mpl::int_<Dim>() );
    }

    bool operator<( Geo0D const& e ) const
    {
        return this->id() < e.id();
    };

    bool operator<( size_type __i ) const
    {
        return this->id() < __i;
    };

    /**
     * show the information about the Geo0D
     *
     * @param verbose true if verbose mode, false otherwise
     * @param c the output stream to use
     *
     * @return the output stream
     */
    std::ostream & showMe( bool verbose = false, std::ostream & c = std::cout ) const;

    /**
     * set the point coordinates. This will typically be called when
     * creating faces (points) in 1D.
     */
    void setPoint( uint16_type const /*i*/, self_type const & p )
    {
        _M_coor = p._M_coor;
    }

    /**
     * translate the point by \p trans
     */
    self_type& translate( node_type const& trans )
    {
        _M_coor += trans;
        return *this;
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
    bool isEqual( Geo0D const& geo0d, mpl::int_<1> ) const
    {
        return std::abs( _M_coor[0] - geo0d._M_coor[0] ) < 1e-10;
    }
    bool isEqual( Geo0D const& geo0d, mpl::int_<2> ) const
    {
        return ( std::abs( _M_coor[0] - geo0d._M_coor[0] ) < 1e-10 &&
                 std::abs( _M_coor[1] - geo0d._M_coor[1] ) < 1e-10 );
    }
    bool isEqual( Geo0D const& geo0d, mpl::int_<3> ) const
    {
        return ( std::abs( _M_coor[0] - geo0d._M_coor[0] ) < 1e-10 &&
                 std::abs( _M_coor[1] - geo0d._M_coor[1] ) < 1e-10 &&
                 std::abs( _M_coor[2] - geo0d._M_coor[2] ) < 1e-10 );

    }
private:
    node_type _M_coor;

    bool _M_is_vertex;

    // mesh to which the geond element belongs to
    MeshBase const* M_mesh;


    Marker1 _M_marker1;
    Marker2 _M_marker2;
    Marker3 _M_marker3;

};

// Alias for Geo0D<3>
typedef Geo0D<3> Point;


/*--------------------------------------------------------------
  Geo0D
  ---------------------------------------------------------------*/
template<uint16_type Dim, typename T>
Geo0D<Dim, T>::Geo0D()
    :
    super( 0, MESH_ENTITY_INTERNAL ),
    _M_coor(Dim),
    _M_is_vertex( false ),
    _M_marker1(),
    _M_marker2(),
    _M_marker3()
{
    _M_coor.clear();
}

template<uint16_type Dim, typename T>
Geo0D<Dim, T>::Geo0D( size_type id, bool boundary, bool is_vertex )
    :
    super( id, MESH_ENTITY_INTERNAL ),
    _M_coor(Dim),
    _M_is_vertex( is_vertex ),
    _M_marker1(),
    _M_marker2(),
    _M_marker3()
{
    _M_coor.clear();
    this->setOnBoundary( boundary );
}

template<uint16_type Dim, typename T>
Geo0D<Dim, T>::Geo0D( size_type id, value_type x, value_type y, value_type z, bool boundary, bool is_vertex )
    :
    super( id, MESH_ENTITY_INTERNAL ),
    _M_coor(Dim),
    _M_is_vertex( is_vertex ),
    _M_marker1(),
    _M_marker2(),
    _M_marker3()
{
    if ( Dim < 2 )
        _M_coor[ 0 ] = x;
    if (  Dim < 3 )
        _M_coor[ 1 ] = y;
    if ( Dim == 3 )
        _M_coor[ 2 ] = z;

    this->setOnBoundary( boundary );
}

template<uint16_type Dim, typename T>
Geo0D<Dim, T>::Geo0D( size_type id, node_type const& __p, bool boundary, bool is_vertex )
    :
    super( id, MESH_ENTITY_INTERNAL ),
    _M_coor(__p ),
    _M_is_vertex( is_vertex ),
    _M_marker1(),
    _M_marker2(),
    _M_marker3()
{
    LIFE_ASSERT( __p.size() == Dim )( __p )( Dim ).error( "invalid node" );

    this->setOnBoundary( boundary );
}

template<uint16_type Dim, typename T>
Geo0D<Dim, T>::Geo0D( Geo0D const & G )
    :
    super( G ),
    _M_coor( G._M_coor ),
    _M_is_vertex( G._M_is_vertex ),
    _M_marker1( G._M_marker1 ),
    _M_marker2( G._M_marker2 ),
    _M_marker3( G._M_marker3 )
{
}

template<uint16_type Dim, typename T>
Geo0D<Dim, T> &
Geo0D<Dim, T>::operator=( Geo0D<Dim, T> const & G )
{
    if (  this == &G )
        return *this;
    super::operator=( G );
    _M_coor = G._M_coor;
    _M_is_vertex = G._M_is_vertex;
    _M_marker1 = G._M_marker1;
    _M_marker2 = G._M_marker2;
    _M_marker3 = G._M_marker3;
    return *this;
}

template<uint16_type Dim, typename T>
std::ostream &
Geo0D<Dim, T>::showMe( bool /*verbose*/, std::ostream & out ) const
{
    out.setf( std::ios::scientific, std::ios::floatfield );
    out << "----- BEGIN of Geo0D ---\n";
    out << "id = " << id() << " node:" << node() << "\n";
    out << "is a vertex = " << isVertex() << "\n";
    out << "----- END OF Geo0D ---\n";
    return out;
}

template<uint16_type Dim, typename T>
inline
DebugStream&
operator<<( DebugStream& __os, Geo0D<Dim, T> const& __n )
{
    if ( __os.doPrint() )
    {
        std::ostringstream __str;

        __str << __n.showMe( true, __str );

        __os << __str.str() << "\n";
    }
    return __os;
}
template<uint16_type Dim, typename T>
inline
NdebugStream&
operator<<( NdebugStream& __os, Geo0D<Dim, T> const& __n ) { return __os; }


template<typename T>
inline
T
distance( Geo0D<1,T> const& p1, Geo0D<1,T> const& p2 )
{
    return ublas::norm_2( p1.node()-p2.node() );
}

template<typename T>
inline
T
distance( Geo0D<2,T> const& p1, Geo0D<2,T> const& p2 )
{
    return ublas::norm_2( p1.node()-p2.node() );
}

template<typename T>
inline
T
distance( Geo0D<3,T> const& p1, Geo0D<3,T> const& p2 )
{
    return ublas::norm_2( p1.node()-p2.node() );
}
template<typename T>
inline
Geo0D<1,T>
middle( Geo0D<1,T> const& p1, Geo0D<1,T> const& p2 )
{
    return Geo0D<1,T>( ( p1.node()+p2.node() )/2 );
}
template<typename T>
inline
Geo0D<2,T>
middle( Geo0D<2,T> const& p1, Geo0D<2,T> const& p2 )
{
    return Geo0D<2,T>( ( p1.node()+p2.node() )/2 );
}
template<typename T>
inline
Geo0D<3,T>
middle( Geo0D<3,T> const& p1, Geo0D<3,T> const& p2 )
{
    return Geo0D<3,T>( ( p1.node()+p2.node() )/2 );
}

} // Life

#endif

