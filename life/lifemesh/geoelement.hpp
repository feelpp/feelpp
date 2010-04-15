/*
 This file is part of the Life library
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
/*!
  \file geoElement.hpp
  \brief Geometric elements
  Introduces all the geometric elements
*/

#ifndef _GEOELEMENT_HH_
#define _GEOELEMENT_HH_

#include <boost/optional.hpp>
#include <boost/tuple/tuple.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifealg/matrix.hpp>
#include <life/lifemesh/marker.hpp>
#include <life/lifemesh/geond.hpp>
#include <life/lifemesh/bareitems.hpp>
#include <life/lifealg/lu.hpp>

namespace Life
{

class SubFaceOfNone
{
public:
    template<typename ET>
    struct Element
    {
        typedef ET type;
    };
    typedef boost::tuple<size_type, size_type, uint16_type> element_connectivity_type;

    boost::none_t element( uint16_type /* e */ ) const { return boost::none_t(); }

    SubFaceOfNone(){}
    SubFaceOfNone( SubFaceOfNone const& ){}

    virtual ~SubFaceOfNone() {}

    template<typename SFO>
    SubFaceOfNone( SFO const& /*sf*/ )
    {
    }
    bool
    isInterProcessDomain( size_type /*p*/ ) const
    {
        return false;
    }
};

template<typename ElementType>
class SubFaceOf
{
public:

    template<typename ET>
    struct Element
    {
        typedef ElementType type;
    };
    typedef ElementType entity_type;
    typedef boost::tuple<ElementType const*, size_type, uint16_type, size_type> element_connectivity_type;

    SubFaceOf()
        :
        _M_element0( 0, invalid_size_type_value, invalid_uint16_type_value, invalid_size_type_value ),
        _M_element1( 0, invalid_size_type_value, invalid_uint16_type_value, invalid_size_type_value )
    {}

    SubFaceOf( element_connectivity_type const& connect0 )
        :
        _M_element0( connect0 ),
        _M_element1( 0, invalid_size_type_value, invalid_uint16_type_value, invalid_uint16_type_value )
    {}
    SubFaceOf( element_connectivity_type const& connect0,
               element_connectivity_type const& connect1 )
        :
        _M_element0( connect0 ),
        _M_element1( connect1 )
    {}

    SubFaceOf( SubFaceOf const& sf )
        :
        _M_element0( sf._M_element0 ),
        _M_element1( sf._M_element1 )
    {
    }
    SubFaceOf( SubFaceOfNone const& /*sf*/ )
        :
        _M_element0( 0, invalid_size_type_value, invalid_uint16_type_value, invalid_size_type_value ),
        _M_element1( 0, invalid_size_type_value, invalid_uint16_type_value, invalid_size_type_value )
    {
    }
    virtual ~SubFaceOf() {}

    SubFaceOf& operator=( SubFaceOf const& sf )
    {
        if ( this != &sf )
        {
            _M_element0 = sf._M_element0;
            _M_element1 = sf._M_element1;

        }
        return *this;
    }
    SubFaceOf& operator=( SubFaceOfNone const& /*sf*/ )
    {
        return *this;
    }
    entity_type const& element( uint16_type e ) const
    {
        if ( e == 0 )
            return *boost::get<0>( _M_element0 );
        else
            return *boost::get<0>( _M_element1 );
    }

    entity_type const& element0() const { return *boost::get<0>( _M_element0 ); }
    entity_type const& element1() const { return *boost::get<0>( _M_element1 ); }

    size_type ad_first() const { return boost::get<1>( _M_element0 ); }
    uint16_type pos_first() const { return boost::get<2>( _M_element0 ); }
    size_type proc_first() const { return boost::get<3>( _M_element0 ); }

    size_type ad_second() const { return boost::get<1>( _M_element1 ); }
    uint16_type pos_second() const { return boost::get<2>( _M_element1 ); }
    size_type proc_second() const { return boost::get<3>( _M_element1 ); }

    element_connectivity_type const& connection0() const { return _M_element0; }
    element_connectivity_type const& connection1() const { return _M_element1; }

    void setConnection( uint16_type f, element_connectivity_type const& connect )
    {
        if ( f == 0 )
            _M_element0 = connect;
        else
            _M_element1 = connect;

    }

    void setConnection0( element_connectivity_type const& connect )
    {
        _M_element0 = connect;
    }
    void setConnection1( element_connectivity_type const& connect ) { _M_element1 = connect; }

    bool isConnectedTo0() const { return ( boost::get<1>( _M_element0 ) != invalid_size_type_value &&
                                           boost::get<2>( _M_element0 ) != invalid_uint16_type_value &&
                                           boost::get<3>( _M_element0 ) != invalid_size_type_value ); }
    bool isConnectedTo1() const { return ( boost::get<1>( _M_element1 ) != invalid_size_type_value &&
                                           boost::get<2>( _M_element1 ) != invalid_uint16_type_value &&
                                           boost::get<3>( _M_element1 ) != invalid_size_type_value ); }


    bool
    isInterProcessDomain( size_type p ) const
    {
        return ( ( boost::get<3>( _M_element1 ) != invalid_size_type_value ) &&
                 ( boost::get<3>( _M_element0 ) == p ) &&
                 ( boost::get<3>( _M_element0 ) != boost::get<3>( _M_element1 )) );
    }
    bool
    isIntraProcessDomain( size_type p ) const
    {
        return ( ( boost::get<3>( _M_element0 ) == p ) &&
                 ( boost::get<3>( _M_element1 ) == p ) );
    }

    void disconnect()
    {
        _M_element0 = boost::make_tuple( (ElementType const*)0,
                                         invalid_size_type_value,
                                         invalid_uint16_type_value,
                                         invalid_size_type_value );
        _M_element1 = boost::make_tuple( (ElementType const*)0,
                                         invalid_size_type_value,
                                         invalid_uint16_type_value,
                                         invalid_size_type_value );
    }

private:

    element_connectivity_type _M_element0;
    element_connectivity_type _M_element1;

};

//     *********** Geometrical Elements *****************
//! \defgroup GeoEle Geometry Element classes
/*@{*/

/**
 * Class for Points and Vertices
 *
 * \bug: in the 1D mesh case the points are subfaces of the segments
 *       but this is not handled yet, fixes in regionmesh1D needed
 */
template <uint16_type Dim,
          typename SubFace = SubFaceOfNone,
          typename T = double>
class GeoElement0D
    :
    public Geo0D<Dim,T>,
    public SubFace
{
public:

    static const uint16_type nDim = 0;
    static const uint16_type nRealDim = Dim;

    typedef Geo0D<Dim,T> geo0d_type;
    typedef typename geo0d_type::node_type node_type;

    typedef geo0d_type super;
    typedef SubFace super2;

    typedef GeoElement0D<Dim,SubFace,T> self_type;
    typedef typename SubFace::template Element<self_type>::type element_type;
    typedef self_type point_type;


    GeoElement0D()
        :
        super(),
        super2(),
        _M_marker1(),
        _M_marker2(),
        _M_marker3(),
        _M_pt()
    {}


    //! Declares item id and if it is on boundary
    GeoElement0D( size_type id, bool boundary = false )
        :
        super( id, boundary ),
        super2(),
        _M_marker1(),
        _M_marker2(),
        _M_marker3(),
        _M_pt()
    {}

    GeoElement0D( size_type id, node_type const& n,  bool boundary = false )
        :
        super( id, n, boundary ),
        super2(),
        _M_marker1(),
        _M_marker2(),
        _M_marker3(),
        _M_pt()
    {}

    //! Declares item id and if it is on boundary, and provides coordinate
    //! data.
    GeoElement0D( size_type id, Real x, Real y, Real z, bool boundary = false )
        :
        super( id, x, y, z, boundary ),
        super2(),
        _M_marker1(),
        _M_marker2(),
        _M_marker3(),
        _M_pt()
    {}

    template<typename SF>
    GeoElement0D( GeoElement0D<Dim,SF,T> const & g )
        :
        super( g ),
        super2( g ),
        _M_marker1( g.marker() ),
        _M_marker2( g.marker2() ),
        _M_marker3( g.marker3() ),
        _M_pt( g._M_pt )
    {
    }

    ~GeoElement0D()
    {}

    template<typename SF>
    GeoElement0D & operator = ( GeoElement0D<Dim,SF,T> const & g )
    {
        super::operator=( g );
        super2::operator=( g );
        _M_marker1 = g.marker();
        _M_marker2 = g.marker2();
        _M_marker3 = g.marker3();
        _M_pt = g._M_pt;
        return *this;
    }

    //void setMesh( MeshBase const* m ) { super::setMesh( m ); }
    MeshBase const* mesh() const { return super::mesh(); }

    /**
     * \return process id
     */
    uint16_type processId() const { return super::processId(); }

    /**
     * \return \p true if interprocess domain face, \p false otherwise
     */
    bool isInterProcessDomain() const
    {
        return super2::isInterProcessDomain( super::processId()  );
    }


    /**
     * \return \c true if on the boundary, \c false otherwise
     */
    bool isOnBoundary() const { return super::isOnBoundary(); }

    /**
     * \return the point associated to the face
     */
    geo0d_type const& point( uint16_type /*i*/ ) const { return _M_pt; }

    /**
     * set the geometrical point associated to the face
     */
    //void setPoint( uint16_type /*i*/, GeoElement0D<Dim,SubFaceOfNone,T> const& e ) { _M_pt = e; }
    void setPoint( uint16_type /*i*/, geo0d_type const& e ) { _M_pt = e; }

    /**
     * \return marker1() index
     */
    Marker1 const& marker() const { return _M_marker1; }
    /**
     * \return marker1() index
     */
    Marker1& marker() { return _M_marker1; }
    /**
     * \return marker2() index
     */
    Marker2 const& marker2() const { return _M_marker2; }
    /**
     * \return marker2() index
     */
    Marker2& marker2() { return _M_marker2; }
    /**
     * \return marker3() index
     */
    Marker3 const& marker3() const { return _M_marker3; }
    /**
     * \return marker3() index
     */
    Marker3& marker3() { return _M_marker3; }

private:
    Marker1 _M_marker1;
    Marker2 _M_marker2;
    Marker3 _M_marker3;

    geo0d_type _M_pt;
};




/**
 * \class GeoElement1D
 * \brief class for 1D elements
 *
 * In the 2D case, we store the size_types of the adjacent 2D elements
 * and their relative position.
 */
template<uint16_type Dim,
         typename GEOSHAPE,
         typename SubFace = SubFaceOfNone,
         typename T = double>
class GeoElement1D
    :
    public GeoND<Dim, GEOSHAPE, T, GeoElement0D<Dim, SubFaceOfNone, T> >,
    public SubFace
{
public:

    //enum { nDim = Dim };



    typedef GeoND<Dim, GEOSHAPE, T, GeoElement0D<Dim, SubFaceOfNone, T> > super;
    typedef SubFace super2;

    static const uint16_type nDim = super::nDim;
    static const uint16_type nOrder = super::nOrder;
    static const uint16_type nRealDim = super::nRealDim;

    static const bool condition = (Dim==nRealDim);
    BOOST_MPL_ASSERT_MSG( (condition), INVALID_ELEMENT_REAL_DIMENSION, ( mpl::int_<Dim>, mpl::int_<nRealDim>, GEOSHAPE ) );

    typedef GEOSHAPE GeoShape;
    typedef GeoElement1D<Dim, GEOSHAPE, SubFace, T > self_type;
    //typedef typename SubFace::template Element<self_type>::type element_type;
    typedef self_type element_type;
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<nDim>,mpl::int_<1> >,
                              mpl::identity<GeoElement0D<Dim, SubFaceOf<self_type>, T> >,
                              mpl::identity<GeoElement0D<Dim, SubFaceOfNone, T> > >::type::type point_type;
    typedef point_type GeoBElement;

    static const uint16_type numLocalVertices = super::numLocalVertices;
    static const uint16_type numLocalEdges = super::numEdges;
    static const uint16_type numLocalFaces = super::numLocalVertices;


    typedef typename super::node_type node_type;
    typedef typename super::vertex_permutation_type vertex_permutation_type;
    typedef typename super::edge_permutation_type edge_permutation_type;
    typedef typename super::face_permutation_type face_permutation_type;

    /**
     * default constructor, make it explicit to avoid implict
     * inversion to \c size_type
     */
    explicit GeoElement1D( size_type id = 0 )
        :
        super( id ),
        super2(),
        _M_h(),
        _M_h_face( 2 ),
        _M_vertices( numLocalVertices, 0 )
    {
    }
    /**
     * copy consttructor
     */
    GeoElement1D( GeoElement1D const& g )
        :
        super( g ),
        super2( g ),
        _M_h( g._M_h ),
        _M_h_face( g._M_h_face ),
        _M_vertices( g._M_vertices )
    {}

    /**
     * destructor
     */
    ~GeoElement1D()
    {}

    /**
     * copy operator
     */
    GeoElement1D& operator=( GeoElement1D const& g )
    {
        if ( this != &g )
        {
            super::operator=( g );
            super2::operator=( g );
            _M_h = g._M_h;
            _M_h_face = g._M_h_face;
            _M_vertices = g._M_vertices;
        }
        return *this;
    }

    /**
     * get the max length of the edges of the element
     *
     *
     * @return the max length of the edges of the element
     */
    double h() const
    {
        if ( !_M_h && this->hasPoints() )
        {
            double __h = 0;
            for ( uint16_type __e = 0;__e < numLocalEdges;++__e )
            {
                node_type const& __x1 = this->point( this->eToP( __e, 0 ) ).node();
                node_type const& __x2 = this->point( this->eToP( __e, 1 ) ).node();
                double __l = ublas::norm_2( __x1-__x2 );
                __h = ( __h > __l )?__h:__l;
            }
            _M_h = __h;

        }
        return _M_h.get();

    }

    //void setMesh( MeshBase const* m ) { super::setMesh( m ); }
    MeshBase const* mesh() const { return super::mesh(); }

    /**
     * get the max length of the edge in the local face \c f
     *
     * @param f local id of the face
     *
     * @return the max length of the edges of the local face
     */
    double hFace( uint16_type f ) const
    {
        if ( !_M_h_face[f] )
        {
            _M_h_face[f] = 1;
        }
        return _M_h_face[f].get();
    }

    /**
     * \return \c true if on the boundary, \c false otherwise
     */
    size_type id() const { return super::id(); }

    /**
     * \return \p true if interprocess domain face, \p false otherwise
     */
    bool isInterProcessDomain() const
    {
        return super2::isInterProcessDomain( super::processId() );
    }

    /**
     * \return \c true if on the boundary, \c false otherwise
     */
    bool isOnBoundary() const { return super::isOnBoundary(); }

    /**
     * \return process id
     */
    uint16_type processId() const { return super::processId(); }

    void setMap( uint8_type k_1, uint8_type k_2 )
    {
        _M_map[k_1] = k_2;
    }

    uint8_type map( uint8_type k_1 ) const
    {
        return _M_map[ k_1 ];
    }

    Marker1 const& marker() const { return super::marker(); }
    Marker1& marker() { return super::marker(); }
    Marker2 const& marker2() const { return super::marker2(); }
    Marker3 const& marker3() const { return super::marker3(); }

    /**
     * Inserts a point as face of the edge geometric element
     */
    void setFace( uint16_type const i, point_type const & p )
    {
        LIFE_ASSERT( i < numLocalVertices )( i ).error( "invalid local point index" );
        _M_vertices[i] = const_cast<point_type*>( boost::addressof( p ) );
    }

    edge_permutation_type permutation( uint16_type /*i*/ ) const
    {
        return edge_permutation_type();
    }

    point_type const& face( uint16_type i ) const
    {
        return *_M_vertices[i];
    }
    point_type const* facePtr( uint16_type i ) const
    {
        LIFE_ASSERT( i < numLocalVertices )( this->id() )( i ).error( "invalid local vertex index" );
        return _M_vertices[i];
    }
    point_type* facePtr( uint16_type i )
    {
        LIFE_ASSERT( i < numLocalVertices )( this->id() )( i ).error( "invalid local vertex index" );
        return _M_vertices[i];
    }

    typedef typename ublas::bounded_array<point_type*, numLocalVertices>::iterator face_iterator;
    typedef typename ublas::bounded_array<point_type*, numLocalVertices>::const_iterator face_const_iterator;

    /**
     * \return the iterator pair (begin,end) of faces
     */
    std::pair<face_iterator,face_iterator>
    faces()
    {
        return std::make_pair( _M_vertices.begin(), _M_vertices.end() );
    }

    /**
     * \return the iterator pair (begin,end) of faces
     */
    std::pair<face_const_iterator,face_const_iterator>
    faces() const
    {
        return std::make_pair( _M_vertices.begin(), _M_vertices.end() );
    }
private:

    std::vector<uint8_type> _M_map;
    mutable boost::optional<double> _M_h;
    mutable ublas::vector<boost::optional<double>,ublas::bounded_array<boost::optional<double>, 2> > _M_h_face;
    ublas::bounded_array<point_type*, numLocalVertices> _M_vertices;

};


/**
 * \class GeoElement2D
 * \brief  Class for 2D elements.
 *
 * In the 3D case, we store the size_types of the adjacent 3D elements
 * and their relative position.
 */
template<uint16_type Dim,
         typename GEOSHAPE,
         typename SubFace = SubFaceOfNone,
         typename T = double>
class GeoElement2D
    :
    public GeoND<Dim, GEOSHAPE, T, GeoElement0D<Dim, SubFaceOfNone, T> >,
        public SubFace
{
public:


    typedef GeoND<Dim, GEOSHAPE, T, GeoElement0D<Dim, SubFaceOfNone, T> > super;
    typedef SubFace super2;

    static const uint16_type nDim = super::nDim;
    static const uint16_type nOrder = super::nOrder;
    static const uint16_type nRealDim = super::nRealDim;

    static const bool condition = (Dim==nRealDim);
    BOOST_MPL_ASSERT_MSG( (condition), INVALID_ELEMENT_REAL_DIMENSION, ( mpl::int_<Dim>, mpl::int_<nRealDim>, GEOSHAPE ) );

    //! Number of element edges
    static const uint16_type numLocalEdges = super::numEdges;
    static const uint16_type numLocalFaces = super::numFaces;

    typedef GEOSHAPE GeoShape;
    typedef typename super::face_type entity_face_type;
    typedef GeoElement2D<Dim, GEOSHAPE,SubFace, T> self_type;
    //typedef typename SubFace::template Element<self_type>::type element_type;
    typedef self_type element_type;
    typedef GeoElement1D<Dim, entity_face_type, SubFaceOf<self_type>, T > edge_type;
    typedef GeoElement0D<Dim, SubFaceOfNone, T> point_type;
#if 0
    BOOST_MPL_ASSERT_MSG( (boost::is_same<point_type,typename edge_type::point_type>::value),
                          INCOMPATIBLE_POINT_TYPE,
                          (point_type, typename edge_type::point_type, edge_type, element_type, self_type ) );
    BOOST_STATIC_ASSERT( (boost::is_same<point_type,typename edge_type::point_type>::value) );
#endif // 0
    typedef typename super::node_type node_type;

    typedef typename super::vertex_permutation_type vertex_permutation_type;
    typedef typename super::edge_permutation_type edge_permutation_type;
    typedef typename super::edge_permutation_type permutation_type;
    typedef typename super::face_permutation_type face_permutation_type;
    typedef typename super2::element_connectivity_type element_connectivity_type;

    /**
     * default constructor, make it explicit to avoid implict
     * inversion to \c size_type
     */
    explicit GeoElement2D( size_type id = 0 )
        :
        super( id ),
        super2(),
        _M_h(),
        _M_h_face(numLocalEdges),
        _M_edges( numLocalEdges ),
        _M_edge_permutation( numLocalEdges )
    {
        std::fill( _M_edges.begin(), _M_edges.end(), (edge_type*)0 );
        std::fill( _M_edge_permutation.begin(),
                   _M_edge_permutation.end(),
                   edge_permutation_type(edge_permutation_type::IDENTITY) );
    }

    /**
     * copy consttructor
     */
    GeoElement2D( GeoElement2D const& g )
        :
        super( g ),
        super2( g ),
        _M_h( g._M_h ),
        _M_h_face( g._M_h_face ),
        _M_edges( g._M_edges ),
        _M_edge_permutation( g._M_edge_permutation )
    {}

    /**
     * destructor
     */
    ~GeoElement2D()
    {}

    /**
     * copy operator
     */
    GeoElement2D& operator=( GeoElement2D const& g )
    {
        if ( this != &g )
        {
            super::operator=( g );
            super2::operator=( g );
            _M_h = g._M_h;
            _M_h_face = g._M_h_face;
            _M_edges = g._M_edges;
            _M_edge_permutation = g._M_edge_permutation;
        }
        return *this;
    }

    //void setMesh( MeshBase const* m ) { super::setMesh( m ); }
    MeshBase const* mesh() const { return super::mesh(); }
    /**
     * \return \c true if on the boundary, \c false otherwise
     */
    size_type id() const { return super::id(); }

    Marker1 const& marker() const { return super::marker(); }
    Marker1& marker() { return super::marker(); }
    Marker2 const& marker2() const { return super::marker2(); }
    Marker3 const& marker3() const { return super::marker3(); }

    /**
     * \return \p true if interprocess domain face, \p false otherwise
     */
    bool isInterProcessDomain() const
    {
        return super2::isInterProcessDomain( super::processId() );
    }

    /**
     * \return \c true if on the boundary, \c false otherwise
     */
    bool isOnBoundary() const { return super::isOnBoundary(); }

    /**
     * \return process id
     */
    uint16_type processId() const { return super::processId(); }

    /**
     * get the max length of the edges of the element
     *
     *
     * @return the max length of the edges of the element
     */
    double h() const
    {
        if ( !_M_h && this->hasPoints() )
        {
            double __h = 0;
            for ( uint16_type __e = 0;__e < numLocalEdges;++__e )
            {
                node_type const& __x1 = this->point( this->fToP( __e, 0 ) ).node();
                node_type const& __x2 = this->point( this->fToP( __e, 1 ) ).node();
                double __l = ublas::norm_2( __x1-__x2 );
                __h = ( __h > __l )?__h:__l;
            }
            _M_h = __h;

        }
        return _M_h.get();

    }

    /**
     * get the max length of the edge in the local face \c f
     *
     * @param f local id of the face
     *
     * @return the max length of the edges of the local face
     */
    double hFace( uint16_type f ) const
    {
        if ( !_M_h_face[f] && this->hasPoints() )
        {
            node_type const& __x1 = this->point( this->fToP( f, 0 ) ).node();
            node_type const& __x2 = this->point( this->fToP( f, 1 ) ).node();
            _M_h_face[f] = ublas::norm_2( __x1-__x2 );
        }
        return _M_h_face[f].get();
    }

    /**
     * \sa face()
     */
    edge_type const& edge( uint16_type i ) const
    {
        LIFE_ASSERT( i < numLocalEdges )( i )( numLocalEdges ).error( "invalid local edge index" );
        LIFE_ASSERT( _M_edges[i] )( i ).error( "invalid edge (null pointer)" );
        return boost::cref( *_M_edges[i] );
    }

    /**
     * \sa edge()
     */
    edge_type const& face( uint16_type i ) const
    {
        LIFE_ASSERT( i < numLocalEdges )( i )( numLocalEdges ).error( "invalid local edge index" );
        LIFE_ASSERT( _M_edges[i] )( i ).error( "invalid edge (null pointer)" );
        return boost::cref( *_M_edges[i] );
    }

    edge_type const* facePtr( uint16_type i ) const
    {
        LIFE_ASSERT( i < numLocalEdges )( this->id() )( i ).error( "invalid local edge index" );
        return _M_edges[i];
    }

    /**
     * Inserts an edge.
     * \sa setEdge()
     */
    void setFace( uint16_type const i, edge_type const & p )
    {
        LIFE_ASSERT( i < numLocalEdges )( i ).error( "invalid local edge index" );
        _M_edges[i] = const_cast<edge_type*>( boost::addressof( p ) );
    }

    /**
     * \sa facePermutation(), permutation()
     */
    edge_permutation_type edgePermutation( uint16_type i ) const
    {
        LIFE_ASSERT( i < numLocalEdges )( i )( numLocalEdges ).error( "invalid local edge index" );
        return _M_edge_permutation[i];
    }
    /**
     * \sa edgePermutation(), permutation()
     */
    edge_permutation_type facePermutation( uint16_type i ) const
    {
        LIFE_ASSERT( i < numLocalEdges )( i )( numLocalEdges ).error( "invalid local edge index" );
        return _M_edge_permutation[i];
    }

    /**
     * \sa edgePermutation(), facePermutation()
     */
    edge_permutation_type permutation( uint16_type i ) const
    {
        LIFE_ASSERT( i < numLocalEdges )( i )( numLocalEdges ).error( "invalid local edge index" );
        return _M_edge_permutation[i];
    }

    /**
     * Inserts a point.  Uses point references put point
     * \sa setFace()
     */
    void setEdge( uint16_type i, edge_type const & p )
    {
        LIFE_ASSERT( i < numLocalEdges )( i ).error( "invalid local edge index" );
        _M_edges[i] = const_cast<edge_type*>( boost::addressof( p ) );
    }

    void setEdgePermutation( uint16_type i, edge_permutation_type o )
    {
        LIFE_ASSERT( i < numLocalEdges )( i ).error( "invalid local edge index" );
        _M_edge_permutation[i] = o;
    }

    typedef typename ublas::bounded_array<edge_type*, numLocalEdges>::iterator face_iterator;
    typedef typename ublas::bounded_array<edge_type*, numLocalEdges>::const_iterator face_const_iterator;

    /**
     * \return the iterator pair (begin,end) of faces
     */
    std::pair<face_iterator,face_iterator>
    faces()
    {
        return std::make_pair( _M_edges.begin(), _M_edges.end() );
    }

    /**
     * \return the iterator pair (begin,end) of faces
     */
    std::pair<face_const_iterator,face_const_iterator>
    faces() const
    {
        return std::make_pair( _M_edges.begin(), _M_edges.end() );
    }

private:

    mutable boost::optional<double> _M_h;
    mutable ublas::vector<boost::optional<double>,ublas::bounded_array<boost::optional<double>, numLocalEdges> > _M_h_face;

    ublas::bounded_array<edge_type*, numLocalEdges> _M_edges;
    ublas::bounded_array<edge_permutation_type, numLocalEdges> _M_edge_permutation;
};

/*-------------------------------------------------------------------------
  GeoElement2D
  --------------------------------------------------------------------------*/
template <uint16_type Dim, typename GEOSHAPE, typename SFO, typename T>
const uint16_type GeoElement2D<Dim, GEOSHAPE, SFO, T>::numLocalEdges;


/**
 * \class GeoElement3D
 * \brief Class for 3D elements
 *
 */
template<uint16_type Dim,
         typename GEOSHAPE,
         typename T = double>
class GeoElement3D
    :
    public GeoND<Dim, GEOSHAPE, T, GeoElement0D<Dim, SubFaceOfNone, T> >,
    public SubFaceOfNone
{
public:

    enum { nDim = Dim };

    typedef GeoND<Dim, GEOSHAPE, T, GeoElement0D<Dim, SubFaceOfNone, T> > super;
    typedef SubFaceOfNone super2;

    typedef GEOSHAPE GeoShape;

    typedef typename super::face_type entity_face_type;

    typedef GeoElement3D<Dim, GEOSHAPE,T> self_type;
    typedef self_type element_type;
    typedef GeoElement2D<Dim, entity_face_type, SubFaceOf<self_type>, T > face_type;
    typedef GeoElement1D<Dim, typename entity_face_type::topological_face_type, SubFaceOfNone, T> edge_type;
    typedef GeoElement0D<Dim, SubFaceOfNone, T> point_type;

    typedef typename super::node_type node_type;

    typedef typename super::vertex_permutation_type vertex_permutation_type;
    typedef typename super::edge_permutation_type edge_permutation_type;
    typedef typename super::face_permutation_type face_permutation_type;
    typedef typename super::face_permutation_type permutation_type;

    //! Number of local Vertices
    static const uint16_type numLocalVertices = super::numVertices;
    //! Number of local Faces
    static const uint16_type numLocalFaces = super::numFaces;
    //! Number of local Edges (using Euler Formula)
    static const uint16_type numLocalEdges = super::numEdges;



    /**
     *
     *
     * @param id identifier of the element
     */
    explicit GeoElement3D( size_type id = 0 )
        :
        super( id ),
        super2(),
        _M_h(),
        _M_h_face( numLocalEdges ),
        _M_edges( numLocalEdges ),
        _M_faces( numLocalFaces ),
        _M_edge_permutation( numLocalEdges ),
        _M_face_permutation( numLocalFaces )
    {
        std::fill( _M_edge_permutation.begin(), _M_edge_permutation.end(), edge_permutation_type( edge_permutation_type::IDENTITY ) );
        std::fill( _M_face_permutation.begin(), _M_face_permutation.end(), face_permutation_type( face_permutation_type::IDENTITY ) );
    }

    /**
     * copy consttructor
     */
    GeoElement3D( GeoElement3D const& g )
        :
        super( g ),
        super2( g ),
        _M_h( g._M_h ),
        _M_h_face( g._M_h_face ),
        _M_edges( numLocalEdges ),
        _M_faces( numLocalFaces ),
        _M_edge_permutation( g._M_edge_permutation ),
        _M_face_permutation( g._M_face_permutation )
    {}

    /**
     * destructor
     */
    ~GeoElement3D()
    {}

    /**
     * copy operator
     */
    GeoElement3D& operator=( GeoElement3D const& g )
    {
        if ( this != &g )
            {
                super::operator=( g );
                _M_h = g._M_h;
                _M_h_face = g._M_h_face;
                _M_edges = g._M_edges;
                _M_faces = g._M_faces;
                _M_edge_permutation = g._M_edge_permutation;
                _M_face_permutation = g._M_face_permutation;
            }
        return *this;
    }

    //void setMesh( MeshBase const* m ) { super::setMesh( m ); }
    MeshBase const* mesh() const { return super::mesh(); }
    /**
     * \return \c true if on the boundary, \c false otherwise
     */
    size_type id() const { return super::id(); }

    Marker1 const& marker() const { return super::marker(); }
    Marker1& marker() { return super::marker(); }
    Marker2 const& marker2() const { return super::marker2(); }
    Marker3 const& marker3() const { return super::marker3(); }

    /**
     * \return \p true if interprocess domain face, \p false otherwise
     */
    bool isInterProcessDomain() const
    {
        return super2::isInterProcessDomain( super::processId() );
    }

    /**
     * \return \c true if on the boundary, \c false otherwise
     */
    bool isOnBoundary() const { return super::isOnBoundary(); }

    /**
     * \return process id
     */
    uint16_type processId() const { return super::processId(); }

    /**
     * get the max length of the edges of the element
     *
     *
     * @return the max length of the edges of the element
     */
    double h() const
    {
        if ( !_M_h && this->hasPoints() )
            {
                double __h = 0;
                for ( uint16_type __e = 0;__e < numLocalEdges;++__e )
                    {
                        node_type const& __x1 = this->point( this->eToP( __e, 0 ) ).node();
                        node_type const& __x2 = this->point( this->eToP( __e, 1 ) ).node();
                        double __l = ublas::norm_2( __x1-__x2 );
                        __h = ( __h > __l )?__h:__l;
                    }
                _M_h = __h;

            }
        return _M_h.get();

    }

    /**
     * get the max length of the edge in the local face \c f
     *
     * @param f local id of the face
     *
     * @return the max length of the edges of the local face
     */
    double hFace( uint16_type f ) const
    {
        if ( !_M_h_face[f] && this->hasPoints() )
            {
#if 0
                double __h = 0;
                for ( uint16_type __e = 0;__e < FaceShape::numEdges;++__e )
                    {
                        node_type const& __x1 = this->point( super::eToP( super::fToE( f, __e ).first, 0 ) ).node();
                        node_type const& __x2 = this->point( super::eToP( super::fToE( f, __e ).first, 1 ) ).node();
                        double __l = ublas::norm_2( __x1-__x2 );
                        __h = ( __h > __l )?__h:__l;
                    }
                _M_h_face[f] = __h;
#else
                _M_h_face[f] = 1;
#endif
            }
        return _M_h_face[f].get();
    }

    size_type ad_first() const { return invalid_size_type_value; }
    uint16_type pos_first() const { return invalid_uint16_type_value; }
    size_type ad_second() const { return invalid_size_type_value; }
    uint16_type pos_second() const { return invalid_uint16_type_value; }

    edge_type const& edge( uint16_type i ) const
    {
        LIFE_ASSERT( i < numLocalEdges )( i ).error( "invalid local edge index" );
        LIFE_ASSERT( _M_edges[i] )( i ).error( "invalid edge (null pointer)" );
        return *_M_edges[i];
    }

    edge_type const* edgePtr( uint16_type i ) const
    {
        LIFE_ASSERT( i < numLocalEdges )( i ).error( "invalid local edge index" );
        LIFE_ASSERT( _M_edges[i] )( i ).error( "invalid edge (null pointer)" );
        return _M_edges[i];
    }

    edge_permutation_type edgePermutation( uint16_type i ) const
    {
        LIFE_ASSERT( i < numLocalEdges )( i ).error( "invalid local edge index" );

        LIFE_ASSERT( _M_edges[i] )( i ).warn( "invalid edge (null pointer)" );

        return _M_edge_permutation[i];
    }

    /**
     * Inserts an edge
     */
    void setEdge( uint16_type const i, edge_type const & p )
    {
        LIFE_ASSERT( boost::addressof( p ) )( i ).error( "invalid edge (null pointer)" );
        _M_edges[i] = const_cast<edge_type*>( boost::addressof( p ) );
        LIFE_ASSERT( _M_edges[i] )( i ).error( "invalid edge (null pointer)" );

    }

    void setEdgePermutation( uint16_type i, edge_permutation_type o )
    {
        LIFE_ASSERT( i < numLocalEdges )( i ).error( "invalid local edge index" );
        _M_edge_permutation[i] = o;
    }

    face_type const& face( uint16_type i ) const
    {
        LIFE_ASSERT( i < numLocalFaces )( this->id() )( i ).error( "invalid local edge index" );
        LIFE_ASSERT( _M_faces[i] )( this->id() )( i ).error( "invalid edge (null pointer)" );
        return *_M_faces[i];
    }

    face_type const* facePtr( uint16_type i ) const
    {
        LIFE_ASSERT( i < numLocalFaces )( this->id() )( i ).error( "invalid local edge index" );
        //LIFE_ASSERT( _M_faces[i] )( i ).error( "invalid edge (null pointer)" );
        return _M_faces[i];
    }

    face_permutation_type facePermutation( uint16_type i ) const
    {
        LIFE_ASSERT( i < numLocalFaces )( this->id() )( i ).error( "invalid local face index" );
        LIFE_ASSERT( _M_faces[i] )( this->id() )( i ).error( "invalid face (null pointer)" );
        return _M_face_permutation[i];
    }
    face_permutation_type permutation( uint16_type i ) const
    {
        LIFE_ASSERT( i < numLocalFaces )( this->id() )( i ).error( "invalid local face index" );
        LIFE_ASSERT( _M_faces[i] )( this->id() )( i ).error( "invalid face (null pointer)" );
        return _M_face_permutation[i];
    }

    /**
     * Inserts a face.
     */
    void setFace( uint16_type const i, face_type const & p ) { _M_faces[i] = const_cast<face_type*>( boost::addressof( p ) ); }

    void setFacePermutation( uint16_type i, face_permutation_type o )
    {
        LIFE_ASSERT( i < numLocalFaces )( this->id() )( i ).error( "invalid local face index" );
        _M_face_permutation[i] = o;
    }

    typedef typename ublas::bounded_array<face_type*, numLocalFaces>::iterator face_iterator;
    typedef typename ublas::bounded_array<face_type*, numLocalFaces>::const_iterator face_const_iterator;

    /**
     * \return the iterator pair (begin,end) of faces
     */
    std::pair<face_iterator,face_iterator>
    faces()
    {
        return std::make_pair( _M_faces.begin(), _M_faces.end() );
    }

    /**
     * \return the iterator pair (begin,end) of faces
     */
    std::pair<face_const_iterator,face_const_iterator>
    faces() const
    {
        return std::make_pair( _M_faces.begin(), _M_faces.end() );
    }
private:

    mutable boost::optional<double> _M_h;
    mutable ublas::vector<boost::optional<double>,ublas::bounded_array<boost::optional<double>, numLocalEdges> > _M_h_face;

    ublas::bounded_array<edge_type*, numLocalEdges> _M_edges;
    ublas::bounded_array<face_type*, numLocalFaces> _M_faces;

    ublas::bounded_array<edge_permutation_type, numLocalEdges> _M_edge_permutation;
    ublas::bounded_array<face_permutation_type, numLocalFaces> _M_face_permutation;
};
/*@}*/




/*-------------------------------------------------------------------------
  GeoElement3D
  --------------------------------------------------------------------------*/
template <uint16_type Dim, typename GEOSHAPE, typename T>
const uint16_type GeoElement3D<Dim, GEOSHAPE, T>::numLocalVertices;
template <uint16_type Dim, typename GEOSHAPE, typename T>
const uint16_type GeoElement3D<Dim, GEOSHAPE, T>::numLocalFaces;
template <uint16_type Dim, typename GEOSHAPE, typename T>
const uint16_type GeoElement3D<Dim, GEOSHAPE, T>::numLocalEdges;

} // Life
#endif

