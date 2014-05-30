/*
 This file is part of the Feel library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

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
 * \file geo0D.hpp
 */
#ifndef _GEO0D_HH_
#define _GEO0D_HH_

#include <boost/operators.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelmesh/geoentity.hpp>
#include <feel/feelmesh/meshbase.hpp>
#include <feel/feelmesh/marker.hpp>

namespace Feel
{
class MeshBase;

/**
 *  \defgroup GeoXD Basis Geometrical Entities Geo0D and GeoND.
 *  \ingroup Obsolet_Groups
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
public GeoEntity<Simplex<0, 1, Dim> >,
public node<T, Dim>::type
{
    typedef GeoEntity<Simplex<0, 1, Dim> > super;
    typedef typename node<T, Dim>::type super2;
public:

    typedef Geo0D<Dim,T> self_type;
    static const uint16_type nDim = Dim;
    typedef T value_type;
    typedef typename matrix_node<value_type>::type matrix_node_type;
    typedef super2 node_type;
    typedef typename node<T, 2>::type parametric_node_type;
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
        super2( Dim ),
        M_is_vertex( false )
    {
        this->operator[]( 0 ) = x;

        if (  Dim >= 2 )
            this->operator[]( 1 ) = y;

        if ( Dim == 3 )
            this->operator[]( 2 ) = z;
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
        super2( __x ),
        M_is_vertex( false )
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
        super2( __expr ),
        M_is_vertex( false )
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

    template<typename AE>
    Geo0D & operator=( ublas::vector_expression<AE> const& expr )
    {
        super2::operator=( expr );
        return *this;
    }
    template<typename AE>
    Geo0D & operator+=( ublas::vector_expression<AE> const& expr )
    {
        super2::operator+=( expr );
        return *this;
    }

    value_type& operator()( int i )
    {
        return this->operator[]( i );
    }

    value_type  operator()( int i ) const
    {
        return this->operator[]( i );
    }

    /**
     * the master id is different from the id in the context of periodic
     * boundary conditions on the slave side
     * @return the master id
     */
    size_type masterId() const
    {
        return M_master_id;
    }

    self_type const* masterVertex() const
    {
        return M_master_vertex;
    }

    /**
     * @return true if the entity is periodic, false otherwise
     */
    bool isPeriodic() const { return M_master_id != this->M_id; }

    /**
     * \return \p true if point is a vertex, \p false otherwise
     * \attention DO NOT USE YET, returns always false
     */
    bool isVertex() const
    {
        return M_is_vertex;
    }

    /**
     * set the point as a vertex or not using \p v
     */
    void  setAsVertex( bool v )
    {
        M_is_vertex = v;
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
    MeshBase const* mesh() const
    {
        return M_mesh;
    }

    /**
     * @return the node data structure
     */
    Geo0D const& node() const
    {
        return *this;
    }

    /**
     * @return the node data structure
     */
    matrix_node_type G() const
    {
        matrix_node_type __G( Dim, 1 );
        ublas::column( __G, 0 ) = *this;
        return __G;
    }

    /**
     * @return the node data structure
     */
    matrix_node_type vertices() const
    {
        return this->G();
    }

    /**
     * \return the measure of a point
     */
    value_type measure() const
    {
        return 0;
    }

    /**
     * \return \p true if the point has parametric coordinates
     */
    bool isParametric() const
    {
        return M_is_parametric;
    }

    /**
     * \return the geometric dimension of the entity the point belongs to
     */
    int gDim() const
    {
        return M_gdim;
    }

    /**
     * \return the geometric dimension of the entity the point belongs to
     */
    int gTag() const
    {
        return M_gtag;
    }

    /**
     * \return the parametric coordinates
     */
    parametric_node_type const& parametricCoordinates() const
    {
        return M_uv;
    }

    /**
     * \return the first parametric coordinates
     */
    value_type u() const
    {
        return M_uv[0];
    }

    /**
     * \return the first parametric coordinates
     */
    value_type v() const
    {
        return M_uv[1];
    }

    /**
     * set the node coordinates
     *
     * @param __n the node coordinates
     */
    void setNode( node_type const& __n )
    {
        *this = __n;
    }

    /**
     * \return true if points are equal, false otherwise
     */
    bool operator==( Geo0D const& geo0d ) const
    {
        return this->id() == geo0d.id();
    }

    bool operator<( Geo0D const& e ) const
    {
        return this->id() < e.id();
    }

    bool operator<( size_type __i ) const
    {
        return this->id() < __i;
    }

    /**
     * show the information about the Geo0D
     *
     * @param verbose true if verbose mode, false otherwise
     * @param c the output stream to use
     *
     * @return the output stream
     */
    std::ostream & showMe( bool verbose = false, std::ostream & c = std::cout ) const;


    void setMasterId( size_type id )
    {
        M_master_id = id;
    }
    void setMasterVertex( self_type const* m )
    {
        M_master_vertex = m;
    }

    /**
     * set the point coordinates. This will typically be called when
     * creating faces (points) in 1D.
     */
    void setPoint( uint16_type const /*i*/, self_type const & p )
    {
        *this = p;
    }
    void setPointCoordG( int i, ublas::vector<double> const& p )
    {
        *this = p;
    }

    /**
     * translate the point by \p trans
     */
    self_type& translate( node_type const& trans )
    {
        *this += trans;
        return *this;
    }

    Marker1 const& marker() const
    {
        return M_marker1;
    }
    Marker1& marker()
    {
        return M_marker1;
    }
    void setMarker( flag_type v )
    {
        return M_marker1.assign( v );
    }

    Marker2 const& marker2() const
    {
        return M_marker2;
    }
    Marker2& marker2()
    {
        return M_marker2;
    }
    void setMarker2( flag_type v )
    {
        return M_marker2.assign( v );
    }

    Marker3 const& marker3() const
    {
        return M_marker3;
    }
    Marker3& marker3()
    {
        return M_marker3;
    }
    void setMarker3( flag_type v )
    {
        return M_marker3.assign( v );
    }

    /**
     * set the tags associated to the points
     * - tags[0] physical region
     * - tags[1] elementary region
     * - tags[2] particular region
     */
    void setTags( std::vector<int> const& tags )
    {
        this->setMarker( tags[0] );

        if ( tags.size() > 1 )
            this->setMarker2( tags[1] );

        if ( tags.size() > 2 )
            this->setProcessId( tags[2] );
    }

    std::vector<int> tags() const
        {
            std::vector<int> thetags(3);
            thetags[0] = M_marker1.value();
            thetags[1] = M_marker2.value();
            thetags[2] = this->processId();
            return thetags;
        }
    /**
     * set the geometric dimension of the entity the points belongs to
     */
    void setGDim( int gdim )
    {
        M_gdim = gdim;
        M_is_parametric = true;
    }

    /**
     * set the geometric tag of the entity the points belongs to
     */
    void setGTag( int gtag )
    {
        M_gtag = gtag;
        M_is_parametric = true;
    }

    /**
     * set the parametric coordinates of the node (if it is on an point, edge or
     * surface geometric entity)
     */
    void setParametricCoordinates( parametric_node_type const& x )
    {
        M_uv = x;
        M_is_parametric = true;
    }

    /**
     * set the parametric coordinates of the node (if it is on an point, edge or
     * surface geometric entity)
     */
    void setParametricCoordinates( value_type u, value_type v )
    {
        M_uv[0] = u;
        M_uv[1] = v;
        M_is_parametric = true;
    }

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            ar & boost::serialization::base_object<super>( *this );
            ar & boost::serialization::base_object<super2>( *this );
            //ar & M_is_vertex;
            //ar & M_is_parametric;
            ar & M_marker1;
            ar & M_marker2;
            ar & M_marker3;
            /*
            ar & M_gdim;
            ar & M_gtag;
            ar & M_uv;
            */
        }

private:

    size_type M_master_id;
    self_type const* M_master_vertex;

    bool M_is_vertex;
    bool M_is_parametric;

    // mesh to which the geond element belongs to
    MeshBase const* M_mesh;


    Marker1 M_marker1;
    Marker2 M_marker2;
    Marker3 M_marker3;

    int M_gdim;
    int M_gtag;
    parametric_node_type M_uv;
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
    super2( Dim ),
    M_master_id( 0 ),
    M_is_vertex( false ),
    M_is_parametric( false ),
    M_marker1(),
    M_marker2(),
    M_marker3(),
    M_gdim( 0 ),
    M_gtag( 0 ),
    M_uv( 2 )
{
    this->clear();
}

template<uint16_type Dim, typename T>
Geo0D<Dim, T>::Geo0D( size_type id, bool boundary, bool is_vertex )
    :
    super( id, MESH_ENTITY_INTERNAL ),
    super2( Dim ),
    M_master_id( id ),
    M_is_vertex( is_vertex ),
    M_marker1(),
    M_marker2(),
    M_marker3(),
    M_gdim( 0 ),
    M_gtag( 0 ),
    M_uv( 2 )
{
    this->clear();
    this->setOnBoundary( boundary );
}

template<uint16_type Dim, typename T>
Geo0D<Dim, T>::Geo0D( size_type id, value_type x, value_type y, value_type z, bool boundary, bool is_vertex )
    :
    super( id, MESH_ENTITY_INTERNAL ),
    super2( Dim ),
    M_master_id( id ),
    M_is_vertex( is_vertex ),
    M_is_parametric( false ),
    M_marker1(),
    M_marker2(),
    M_marker3(),
    M_gdim( 0 ),
    M_gtag( 0 ),
    M_uv( 2 )
{
    this->operator[]( 0 ) = x;

    if (  Dim >= 2 )
        this->operator[]( 1 ) = y;

    if ( Dim == 3 )
        this->operator[]( 2 ) = z;

    this->setOnBoundary( boundary );
}

template<uint16_type Dim, typename T>
Geo0D<Dim, T>::Geo0D( size_type id, node_type const& __p, bool boundary, bool is_vertex )
    :
    super( id, MESH_ENTITY_INTERNAL ),
    super2( __p ),
    M_master_id( id ),
    M_is_vertex( is_vertex ),
    M_is_parametric( false ),
    M_marker1(),
    M_marker2(),
    M_marker3(),
    M_gdim( 0 ),
    M_gtag( 0 ),
    M_uv( 2 )
{
    FEELPP_ASSERT( __p.size() == Dim )( __p )( Dim ).error( "invalid node" );

    this->setOnBoundary( boundary );
}

template<uint16_type Dim, typename T>
Geo0D<Dim, T>::Geo0D( Geo0D const & G )
    :
    super( G ),
    super2( G ),
    M_master_id( G.id() ),
    M_is_vertex( G.M_is_vertex ),
    M_is_parametric( G.M_is_parametric ),
    M_marker1( G.M_marker1 ),
    M_marker2( G.M_marker2 ),
    M_marker3( G.M_marker3 ),
    M_gdim( G.M_gdim ),
    M_gtag( G.M_gtag ),
    M_uv( G.M_uv )
{
}

template<uint16_type Dim, typename T>
Geo0D<Dim, T> &
Geo0D<Dim, T>::operator=( Geo0D<Dim, T> const & G )
{
    if (  this == &G )
        return *this;

    super::operator=( G );
    super2::operator=( G );
    M_master_id = G.masterId();
    M_is_vertex = G.M_is_vertex;
    M_is_parametric = G.M_is_parametric;
    M_marker1 = G.M_marker1;
    M_marker2 = G.M_marker2;
    M_marker3 = G.M_marker3;
    M_gdim = G.M_gdim;
    M_gtag = G.M_gtag;
    M_uv = G.M_uv;
    return *this;
}

template<uint16_type Dim, typename T>
std::ostream &
Geo0D<Dim, T>::showMe( bool /*verbose*/, std::ostream & out ) const
{
    out.setf( std::ios::scientific, std::ios::floatfield );
    out << "----- BEGIN of Geo0D ---\n";
    out << "id = " << this->id() << " node:" << this->node() << "\n";
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
operator<<( NdebugStream& __os, Geo0D<Dim, T> const& __n )
{
    return __os;
}


template<typename T>
inline
T
distance( Geo0D<1,T> const& p1, Geo0D<1,T> const& p2 )
{
    return ublas::norm_2( p1-p2 );
}

template<typename T>
inline
T
distance( Geo0D<2,T> const& p1, Geo0D<2,T> const& p2 )
{
    return ublas::norm_2( p1-p2 );
}

template<typename T>
inline
T
distance( Geo0D<3,T> const& p1, Geo0D<3,T> const& p2 )
{
    return ublas::norm_2( p1-p2 );
}
template<typename T>
inline
Geo0D<1,T>
middle( Geo0D<1,T> const& p1, Geo0D<1,T> const& p2 )
{
    return ( p1+p2 )/2;
}

template<typename T>
inline
Geo0D<2,T>
middle( Geo0D<2,T> const& p1, Geo0D<2,T> const& p2 )
{
    return ( p1+p2 )/2;
}

template<typename T>
inline
Geo0D<3,T>
middle( Geo0D<3,T> const& p1, Geo0D<3,T> const& p2 )
{
    return ( p1+p2 )/2;
}

template<typename E1,typename E2>
inline
ublas::vector<double>
cross( ublas::vector_expression<E1>  _p1,
       ublas::vector_expression<E2>  _p2 )
{
    ublas::vector<double> v( 3 );
    ublas::vector<double> p1( _p1 );
    ublas::vector<double> p2( _p2 );
    v( 0 ) = p1( 1 )*p2( 2 )-p1( 2 )*p2( 1 );
    v( 1 ) = p1( 2 )*p2( 0 )-p1( 0 )*p2( 2 );
    v( 2 ) = p1( 0 )*p2( 1 )-p1( 1 )*p2( 0 );
    return v;
}

template<typename T>
inline
ublas::vector<double>
cross( Geo0D<3,T> p1,
       Geo0D<3,T> p2 )
{
    ublas::vector<double> v( 3 );
    v( 0 ) = p1( 1 )*p2( 2 ) - p1( 2 )*p2( 1 );
    v( 1 ) = p1( 2 )*p2( 0 ) - p1( 0 )*p2( 2 );
    v( 2 ) = p1( 0 )*p2( 1 ) - p1( 1 )*p2( 0 );
    return v;
}


} // Feel

#endif
