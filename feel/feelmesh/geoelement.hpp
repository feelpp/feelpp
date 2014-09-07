/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano
   Copyright (C) 2006-2010 Université Joseph Fourier (Grenoble I)

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

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/matrix.hpp>
#include <feel/feelmesh/marker.hpp>
#include <feel/feelmesh/geond.hpp>

#include <feel/feelalg/lu.hpp>

namespace Feel
{

class SubFaceOfNone
{
public:
    static const uint16_type nDim = 1;
    template<typename ET>
    struct Element
    {
        typedef ET type;
    };
    typedef boost::tuple<size_type, size_type, uint16_type> element_connectivity_type;

    boost::none_t element( uint16_type /* e */ ) const
    {
        return boost::none_t();
    }

    SubFaceOfNone() {}
    SubFaceOfNone( SubFaceOfNone const& ) {}

    virtual ~SubFaceOfNone() {}

    template<typename SFO>
    SubFaceOfNone( SFO const& /*sf*/ )
    {
    }
    bool
    isGhostFace( rank_type /*p*/ ) const
    {
        return false;
    }
    bool
    isInterProcessDomain( rank_type /*p*/ ) const
    {
        return false;
    }
    rank_type partition1( rank_type /*p*/ ) const { return invalid_rank_type_value; }
    rank_type partition2( rank_type /*p*/ ) const { return invalid_rank_type_value; }
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
        }
};

template<typename ElementType>
class SubFaceOf
{
public:
    static const uint16_type nDim = ElementType::nDim;
    template<typename ET>
    struct Element
    {
        typedef ElementType type;
    };
    typedef ElementType entity_type;
    typedef boost::tuple<ElementType const*, size_type, uint16_type, rank_type> element_connectivity_type;

    SubFaceOf()
        :
        M_element0( 0, invalid_size_type_value, invalid_uint16_type_value, invalid_rank_type_value ),
        M_element1( 0, invalid_size_type_value, invalid_uint16_type_value, invalid_rank_type_value )
    {}

    SubFaceOf( element_connectivity_type const& connect0 )
        :
        M_element0( connect0 ),
        M_element1( 0, invalid_size_type_value, invalid_uint16_type_value, invalid_rank_type_value )
    {}
    SubFaceOf( element_connectivity_type const& connect0,
               element_connectivity_type const& connect1 )
        :
        M_element0( connect0 ),
        M_element1( connect1 )
    {}

    SubFaceOf( SubFaceOf const& sf )
        :
        M_element0( sf.M_element0 ),
        M_element1( sf.M_element1 )
    {
    }
    SubFaceOf( SubFaceOfNone const& /*sf*/ )
        :
        M_element0( 0, invalid_size_type_value, invalid_uint16_type_value, invalid_rank_type_value ),
        M_element1( 0, invalid_size_type_value, invalid_uint16_type_value, invalid_rank_type_value )
    {
    }
    virtual ~SubFaceOf() {}

    SubFaceOf& operator=( SubFaceOf const& sf )
    {
        if ( this != &sf )
        {
            M_element0 = sf.M_element0;
            M_element1 = sf.M_element1;

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
            return *boost::get<0>( M_element0 );

        else
            return *boost::get<0>( M_element1 );
    }

    entity_type const& element0() const
    {
        return *boost::get<0>( M_element0 );
    }
    entity_type const& element1() const
    {
        return *boost::get<0>( M_element1 );
    }

    size_type ad_first() const
    {
        return boost::get<1>( M_element0 );
    }
    uint16_type pos_first() const
    {
        return boost::get<2>( M_element0 );
    }
    rank_type proc_first() const
    {
        return boost::get<3>( M_element0 );
    }
    rank_type partition1( rank_type p ) const
    {
        return p;
    }

    size_type ad_second() const
    {
        return boost::get<1>( M_element1 );
    }
    uint16_type pos_second() const
    {
        return boost::get<2>( M_element1 );
    }
    rank_type proc_second() const
    {
        return boost::get<3>( M_element1 );
    }
    rank_type partition2( size_type p ) const
    {
        return ( p == this->proc_first() )? this->proc_second() : this->proc_first();
    }

    element_connectivity_type const& connection0() const
    {
        return M_element0;
    }
    element_connectivity_type const& connection1() const
    {
        return M_element1;
    }

    void setConnection( uint16_type f, element_connectivity_type const& connect )
    {
        if ( f == 0 )
            M_element0 = connect;

        else
            M_element1 = connect;

    }

    void setConnection0( element_connectivity_type const& connect )
    {
        M_element0 = connect;
    }
    void setConnection1( element_connectivity_type const& connect )
    {
        M_element1 = connect;
    }

    bool isConnected() const { return isConnectedTo0() && isConnectedTo1(); }

    bool isConnectedTo0() const
    {
        return ( boost::get<1>( M_element0 ) != invalid_size_type_value &&
                 boost::get<2>( M_element0 ) != invalid_uint16_type_value &&
                 boost::get<3>( M_element0 ) != invalid_rank_type_value );
    }
    bool isConnectedTo1() const
    {
        return ( boost::get<1>( M_element1 ) != invalid_size_type_value &&
                 boost::get<2>( M_element1 ) != invalid_uint16_type_value &&
                 boost::get<3>( M_element1 ) != invalid_rank_type_value );
    }

    bool
    isGhostFace( rank_type p ) const
    {
        return ( ( boost::get<3>( M_element1 ) != invalid_rank_type_value ) &&
                 ( ( ( boost::get<3>( M_element0 ) == p ) && ( boost::get<3>( M_element1 ) < p ) ) ||
                   ( ( boost::get<3>( M_element0 ) < p ) && ( boost::get<3>( M_element1 ) == p ) ) ) );
    }

    bool
    isInterProcessDomain( rank_type p ) const
    {
        return ( ( boost::get<3>( M_element1 ) != invalid_rank_type_value ) &&
                 ( ( boost::get<3>( M_element0 ) == p ) || ( boost::get<3>( M_element1 ) == p ) ) &&
                 ( boost::get<3>( M_element0 ) != boost::get<3>( M_element1 ) ) );
    }
    bool
    isIntraProcessDomain( rank_type p ) const
    {
        return ( ( boost::get<3>( M_element0 ) == p ) &&
                 ( boost::get<3>( M_element1 ) == p ) );
    }

    void disconnect0()
    {
        M_element0 = boost::make_tuple( ( ElementType const* )0,
                                        invalid_size_type_value,
                                        invalid_uint16_type_value,
                                        invalid_rank_type_value );
    }

    void disconnect1()
    {
        M_element1 = boost::make_tuple( ( ElementType const* )0,
                                        invalid_size_type_value,
                                        invalid_uint16_type_value,
                                        invalid_rank_type_value );
    }

    void disconnect()
    {
        disconnect0();
        disconnect1();
    }

    void disconnect( ElementType const& elem )
    {
        if(boost::get<0>( M_element0 ) == boost::addressof(elem))
        {
            DVLOG(2) << "connecting 1 to 0 and disconnecting 1..\n";
            M_element0 = M_element1;
            disconnect1();
        }
        else
        {
            DVLOG(2) << "disconnecting 1..\n";
            disconnect1();
        }
    }

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
#if 1
            //ar & M_element0.template get<0>();
            ar & M_element0.template get<1>();
            ar & M_element0.template get<2>();
            ar & M_element0.template get<3>();
            M_element0.template get<0>() = 0;


            //ar & M_element1.template get<0>();
            ar & M_element1.template get<1>();
            ar & M_element1.template get<2>();
            ar & M_element1.template get<3>();
            M_element1.template get<0>() = 0;
#endif
        }
private:

    element_connectivity_type M_element0;
    element_connectivity_type M_element1;

};
#if 0
template<typename ElementType>
class SubFaceOfMany
{
public:
    static const uint16_type nDim = ElementType::nDim;
    static const uint16_type nRealDim = ElementType::nRealDim;
    template<typename ET>
    struct Element
    {
        typedef ElementType type;
    };
    typedef ElementType entity_type;
    typedef boost::tuple<ElementType const*, size_type, uint16_type, size_type> element_connectivity_type;

    SubFaceOfMany()
        :
        M_elements()
        {}

    SubFaceOfMany( SubFaceOfMany const& sf )
        :
        M_elements( sf.M_elements )
        {
        }
    SubFaceOfMany( SubFaceOfNone const& /*sf*/ )
        :
        M_elements()
        {
        }
    virtual ~SubFaceOfMany() {}

    SubFaceOfMany& operator=( SubFaceOfMany const& sf )
        {
            if ( this != &sf )
            {
                M_elements = sf.M_elements;
            }

            return *this;
        }
    SubFaceOfMany& operator=( SubFaceOfNone const& /*sf*/ )
    {
        return *this;
    }

    entity_type const& element0() const
    {
        return *boost::get<0>( *M_elements.begin() );
    }
    entity_type const& element1() const
    {
        return *boost::get<0>( *boost::next(M_elements.begin()) );
    }

    void setConnection( element_connectivity_type const& connect )
    {
        this->insert( connect );
    }

    size_type ad_first() const
    {
        return boost::get<1>( *M_elements.begin() );
    }
    uint16_type pos_first() const
    {
        return boost::get<2>( *M_elements.begin() );
    }
    size_type proc_first() const
    {
        return boost::get<3>( *M_elements.begin() );
    }

    size_type ad_second() const
    {
        return boost::get<1>( *boost::next(M_elements.begin()) );
    }
    uint16_type pos_second() const
    {
        return boost::get<2>( *boost::next(M_elements.begin()) );
    }
    size_type proc_second() const
    {
        return boost::get<3>( *boost::next(M_elements.begin()) );
    }


    void setConnection0( element_connectivity_type const& connect )
    {
        M_elements.insert( connect );
    }
    void setConnection1( element_connectivity_type const& connect )
    {
        M_elements.insert( connect );
    }

    element_connectivity_type const& connection0() const
    {
        return *M_elements.begin();
    }
    element_connectivity_type const& connection1() const
    {
        return *boost::next(M_elements.begin());
    }

    bool isConnectedTo0() const
    {
        return ( boost::get<1>( *M_elements.begin() ) != invalid_size_type_value &&
                 boost::get<2>( *M_elements.begin() ) != invalid_uint16_type_value &&
                 boost::get<3>( *M_elements.begin() ) != invalid_size_type_value );
    }
    bool isConnectedTo1() const
    {
        return ( boost::get<1>( *boost::next(M_elements.begin()) ) != invalid_size_type_value &&
                 boost::get<2>( *boost::next(M_elements.begin()) ) != invalid_uint16_type_value &&
                 boost::get<3>( *boost::next(M_elements.begin()) ) != invalid_size_type_value );
    }
    bool
    isInterProcessDomain( size_type p ) const
    {
        return false;
    }
    bool
    isIntraProcessDomain( size_type p ) const
    {
        return true;
    }

    entity_type const& element( uint16_type e ) const
    {
        return *boost::get<0>( *M_elements.begin() );
    }
    void disconnect()
    {
        M_elements.clear();
    }

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
        }
private:

    std::set<element_connectivity_type> M_elements;
};
#else
//template<typename ElementType> class SubFaceOfMany: public SubFaceOf<ElementType> {};
#define SubFaceOfMany SubFaceOf

#endif

//     *********** Geometrical Elements *****************
/** \defgroup GeoEle Geometry Element classes
    \ingroup Obsolet_Groups */
/*@{*/

/**
 * Class for Points and Vertices
 */
template <uint16_type Dim,
         typename SubFace = SubFaceOfNone,
         typename T = double>
class GeoElement0D
    :
    //public GeoND<Dim, GEOSHAPE, T, GeoElement0D<Dim, SubFaceOfNone, T> >,
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
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<SubFace::nDim>, mpl::int_<0> >, mpl::identity<self_type>, mpl::identity<typename SubFace::template Element<self_type>::type> >::type::type element_type;
    typedef self_type point_type;

    typedef typename super::matrix_node_type matrix_node_type;

    static const uint16_type numLocalVertices = super::numVertices;

    GeoElement0D()
        :
        super(),
        super2()
        //M_facept()
    {}


    //! Declares item id and if it is on boundary
    GeoElement0D( size_type id, bool boundary = false )
        :
        super( id, boundary ),
        super2()
        //M_facept()
    {}

    GeoElement0D( size_type id, node_type const& n,  bool boundary = false )
        :
        super( id, n, boundary ),
        super2()
        //M_facept()
    {}

    //! Declares item id and if it is on boundary, and provides coordinate
    //! data.
    GeoElement0D( size_type id, Real x, Real y, Real z, bool boundary = false )
        :
        super( id, x, y, z, boundary ),
        super2()
        //M_facept()
    {}

    template<typename SF>
    GeoElement0D( GeoElement0D<Dim,SF,T> const & g )
        :
        super( g ),
        super2( g ),
        M_facept( g.M_facept )
    {
    }

    ~GeoElement0D()
    {}

    template<typename SF>
    GeoElement0D & operator = ( GeoElement0D<Dim,SF,T> const & g )
    {
        super::operator=( g );
        super2::operator=( g );
        M_facept = g.M_facept;
        return *this;
    }

    //void setMesh( MeshBase const* m ) { super::setMesh( m ); }
    MeshBase const* mesh() const
    {
        return super::mesh();
    }

    /**
     * \return id
     */
    size_type id() const
    {
        return super::id();
    }

    /**
     * \return process id
     */
    rank_type processId() const
    {
        return super::processId();
    }

    /**
     * \return process id
     */
    rank_type partition1() const
    {
        return super2::partition1( super::processId() );
    }

    /**
     * \return process id
     */
    rank_type partition2() const
    {
        return super2::partition2( super::processId() );
    }

    bool isGhostFace() const
    {
        return super2::isGhostFace( super::processId()  );
    }

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
    bool isOnBoundary() const
    {
        return super::isOnBoundary();
    }

    /**
     * \return maximum \c dimension of the sub-entity touching the boundary of the element
     */
    uint16_type boundaryEntityDimension() const
    {
        return 0;
    }

    /**
     * \return \c true if ghost cell, \c false otherwise
     */
    bool isGhostCell() const
    {
        return super::isGhostCell();
    }

    /**
     * \return the point associated to the face
     */
    geo0d_type const& point( uint16_type /*i*/ ) const
    {
        return *( static_cast<geo0d_type *>( M_facept ) );
        //return M_facept;
    }

    /**
     * set the geometrical point associated to the face
     */
    void setPoint( uint16_type /*i*/, geo0d_type const& e )
    {
        //M_facept = e;
        M_facept= const_cast<geo0d_type *>( &e );
    }

    matrix_node_type /*const&*/ G() const
    {
        return this->G( mpl::bool_<boost::is_same<SubFace,SubFaceOfNone>::value>() );
    }

    matrix_node_type /*const&*/ G( mpl::bool_<true> /**/ ) const
    {
        return super::G();
    }

    matrix_node_type /*const&*/ G( mpl::bool_<false> /**/ ) const
    {
        return M_facept->G();
    }

    matrix_node_type /*const&*/ vertices() const
    {
        return this->vertices( mpl::bool_<boost::is_same<SubFace,SubFaceOfNone>::value>() );
    }

    matrix_node_type /*const&*/ vertices( mpl::bool_<true> /**/ ) const
    {
        return super::vertices();
    }

    matrix_node_type /*const&*/ vertices( mpl::bool_<false> /**/ ) const
    {
        return M_facept->vertices();
    }

    /**
     * \return marker1() index
     */
    Marker1 const& marker() const
    {
        return super::marker();
    }
    /**
     * \return marker1() index
     */
    Marker1& marker()
    {
        return super::marker();
    }
    /**
     * \return marker2() index
     */
    Marker2 const& marker2() const
    {
        return super::marker2();
    }
    /**
     * \return marker2() index
     */
    Marker2& marker2()
    {
        return super::marker2();
    }
    /**
     * \return marker3() index
     */
    Marker3 const& marker3() const
    {
        return super::marker3();
    }
    /**
     * \return marker3() index
     */
    Marker3& marker3()
    {
        return super::marker3();
    }

//private:
    geo0d_type* M_facept;

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            ar & boost::serialization::base_object<super>( *this );
            ar & boost::serialization::base_object<super2>( *this );
        }


};

template<uint16_type Dim,
         typename SubFace,
         typename T>
const uint16_type GeoElement0D<Dim,SubFace,T>::numLocalVertices;



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

    typedef GeoND<Dim, GEOSHAPE, T, GeoElement0D<Dim, SubFaceOfNone, T> > super;
    typedef SubFace super2;

    static const uint16_type nDim = super::nDim;
    static const uint16_type nOrder = super::nOrder;
    static const uint16_type nRealDim = super::nRealDim;

    static const bool condition = ( Dim==nRealDim );
    BOOST_MPL_ASSERT_MSG( ( condition ), INVALID_ELEMENT_REAL_DIMENSION, ( mpl::int_<Dim>, mpl::int_<nRealDim>, GEOSHAPE ) );

    typedef GEOSHAPE GeoShape;
    typedef GeoElement1D<Dim, GEOSHAPE, SubFace, T > self_type;
    //typedef typename SubFace::template Element<self_type>::type element_type;
    typedef self_type element_type;
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<nRealDim>,mpl::int_<1> >,
                              mpl::identity<GeoElement0D<Dim, SubFaceOf<self_type>, T> >,
                              mpl::identity<GeoElement0D<Dim, SubFaceOfMany<self_type>, T> > >::type::type point_type;
    typedef point_type GeoBElement;

    static const uint16_type numLocalVertices = super::numLocalVertices;
    static const uint16_type numLocalEdges = super::numEdges;
    static const uint16_type numLocalFaces = super::numLocalVertices;


    typedef typename super::node_type node_type;
    typedef typename super::vertex_permutation_type vertex_permutation_type;
    typedef typename super::edge_permutation_type edge_permutation_type;
    typedef typename super::face_permutation_type face_permutation_type;
    typedef typename super::permutation_type permutation_type;

    /**
     * default constructor, make it explicit to avoid implict
     * inversion to \c size_type
     */
    explicit GeoElement1D( size_type id = 0 )
        :
        super( id ),
        super2(),
        M_vertices( numLocalVertices, 0 ),
        M_vertex_permutation( numLocalVertices )
    {
    }
    /**
     * copy consttructor
     */
    GeoElement1D( GeoElement1D const& g )
        :
        super( g ),
        super2( g ),
        M_vertices( g.M_vertices ),
        M_vertex_permutation( g.M_vertex_permutation )
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
            M_vertices = g.M_vertices;
            M_vertex_permutation = g.M_vertex_permutation;
        }

        return *this;
    }


    //void setMesh( MeshBase const* m ) { super::setMesh( m ); }
    MeshBase const* mesh() const
    {
        return super::mesh();
    }


    /**
     * \return \c true if on the boundary, \c false otherwise
     */
    size_type id() const
    {
        return super::id();
    }

    bool isGhostFace() const
    {
        return super2::isGhostFace( super::processId()  );
    }


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
    bool isOnBoundary() const
    {
        return super::isOnBoundary();
    }

    /**
     * \return maximum \c dimension of the sub-entity touching the boundary of the element
     */
    uint16_type boundaryEntityDimension() const
    {
        return 0;
    }

    /**
     * \return \c true if ghost cell, \c false otherwise
     */
    bool isGhostCell() const
    {
        return super::isGhostCell();
    }

    /**
     * \return process id
     */
    rank_type processId() const
    {
        return super::processId();
    }

    /**
     * \return process id
     */
    rank_type partition1() const
    {
        return super2::partition1( super::processId() );
    }

    /**
     * \return process id
     */
    rank_type partition2() const
    {
        return super2::partition2( super::processId() );
    }

    void setMap( uint8_type k_1, uint8_type k_2 )
    {
        M_map[k_1] = k_2;
    }

    uint8_type map( uint8_type k_1 ) const
    {
        return M_map[ k_1 ];
    }

    Marker1 const& marker() const
    {
        return super::marker();
    }
    Marker1& marker()
    {
        return super::marker();
    }
    Marker2 const& marker2() const
    {
        return super::marker2();
    }
    Marker3 const& marker3() const
    {
        return super::marker3();
    }

    /**
     * Inserts a point as face of the edge geometric element
     */
    void setFace( uint16_type const i, point_type const & p )
    {
        FEELPP_ASSERT( i < numLocalVertices )( i ).error( "invalid local point index" );
        M_vertices[i] = const_cast<point_type*>( boost::addressof( p ) );
    }

    edge_permutation_type permutation( uint16_type /*i*/ ) const
    {
        return edge_permutation_type();
    }

    point_type const& edge( uint16_type i ) const
    {
        return *M_vertices[i];
    }
    point_type const& face( uint16_type i ) const
    {
        return *M_vertices[i];
    }
    point_type const* facePtr( uint16_type i ) const
    {
        FEELPP_ASSERT( i < numLocalVertices )( this->id() )( i ).error( "invalid local vertex index" );
        return M_vertices[i];
    }
    point_type* facePtr( uint16_type i )
    {
        FEELPP_ASSERT( i < numLocalVertices )( this->id() )( i ).error( "invalid local vertex index" );
        return M_vertices[i];
    }

    typedef typename ublas::bounded_array<point_type*, numLocalVertices>::iterator face_iterator;
    typedef typename ublas::bounded_array<point_type*, numLocalVertices>::const_iterator face_const_iterator;

    /**
     * \return the iterator pair (begin,end) of faces
     */
    std::pair<face_iterator,face_iterator>
    faces()
    {
        return std::make_pair( M_vertices.begin(), M_vertices.end() );
    }

    /**
     * \return the iterator pair (begin,end) of faces
     */
    std::pair<face_const_iterator,face_const_iterator>
    faces() const
    {
        return std::make_pair( M_vertices.begin(), M_vertices.end() );
    }

    /**
     * \sa edgePermutation(), permutation()
     */
    vertex_permutation_type facePermutation( uint16_type i ) const
    {
        FEELPP_ASSERT( i < numLocalVertices )( i )( numLocalVertices ).error( "invalid local vertex index" );
        return M_vertex_permutation[i];
    }

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            DVLOG(2) << "Serializing Geoelement1D id: " << this->id() << "...\n";
            ar & boost::serialization::base_object<super>( *this );
            ar & boost::serialization::base_object<super2>( *this );
        }



private:

    std::vector<uint8_type> M_map;
    ublas::bounded_array<point_type*, numLocalVertices> M_vertices;
    ublas::bounded_array<vertex_permutation_type, numLocalVertices> M_vertex_permutation;

};
template<uint16_type Dim,
         typename GEOSHAPE,
         typename SubFace,
         typename T>
const uint16_type GeoElement1D<Dim,GEOSHAPE,SubFace,T>::numLocalVertices;

template<uint16_type Dim,
         typename GEOSHAPE,
         typename SubFace,
         typename T>
const uint16_type GeoElement1D<Dim,GEOSHAPE,SubFace,T>::nDim;

template<uint16_type Dim,
         typename GEOSHAPE,
         typename SubFace,
         typename T>
const uint16_type GeoElement1D<Dim,GEOSHAPE,SubFace,T>::nRealDim;

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

    static const bool condition = ( Dim==nRealDim );
    BOOST_MPL_ASSERT_MSG( ( condition ), INVALID_ELEMENT_REAL_DIMENSION, ( mpl::int_<Dim>, mpl::int_<nRealDim>, GEOSHAPE ) );

    //! Number of element edges
    static const uint16_type numLocalEdges = super::numEdges;
    static const uint16_type numLocalFaces = super::numFaces;

    typedef GEOSHAPE GeoShape;
    typedef typename super::face_type entity_face_type;
    typedef GeoElement2D<Dim, GEOSHAPE,SubFace, T> self_type;
    //typedef typename SubFace::template Element<self_type>::type element_type;
    typedef self_type element_type;
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<nRealDim>,mpl::int_<2> >,
                              mpl::identity<GeoElement1D<Dim, entity_face_type, SubFaceOf<self_type>, T> >,
                              mpl::identity<GeoElement1D<Dim, entity_face_type, SubFaceOfMany<self_type>, T> > >::type::type edge_type;
    //typedef GeoElement1D<Dim, entity_face_type, SubFaceOf<self_type>, T > edge_type;
    typedef GeoElement0D<Dim, SubFaceOfNone, T> point_type;
#if 0
    BOOST_MPL_ASSERT_MSG( ( boost::is_same<point_type,typename edge_type::point_type>::value ),
                          INCOMPATIBLE_POINT_TYPE,
                          ( point_type, typename edge_type::point_type, edge_type, element_type, self_type ) );
    BOOST_STATIC_ASSERT( ( boost::is_same<point_type,typename edge_type::point_type>::value ) );
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
        M_edges( numLocalEdges, nullptr ),
        M_edge_permutation( numLocalEdges, edge_permutation_type( edge_permutation_type::IDENTITY ) )
    {
    }

    /**
     * copy consttructor
     */
    GeoElement2D( GeoElement2D const& g )
        :
        super( g ),
        super2( g ),
        M_edges( g.M_edges ),
        M_edge_permutation( g.M_edge_permutation )
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
            M_edges = g.M_edges;
            M_edge_permutation = g.M_edge_permutation;
        }

        return *this;
    }

    //void setMesh( MeshBase const* m ) { super::setMesh( m ); }
    MeshBase const* mesh() const
    {
        return super::mesh();
    }
    /**
     * \return \c true if on the boundary, \c false otherwise
     */
    size_type id() const
    {
        return super::id();
    }

    Marker1 const& marker() const
    {
        return super::marker();
    }
    Marker1& marker()
    {
        return super::marker();
    }
    Marker2 const& marker2() const
    {
        return super::marker2();
    }
    Marker3 const& marker3() const
    {
        return super::marker3();
    }

    bool isGhostFace() const
    {
        return super2::isGhostFace( super::processId()  );
    }


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
    bool isOnBoundary() const
    {
        return super::isOnBoundary();
    }

    /**
     * \return maximum \c dimension of the sub-entity touching the boundary of the element
     */
    uint16_type boundaryEntityDimension() const
    {
        return super::boundaryEntityDimension();
    }

    /**
     * \return \c true if ghost cell, \c false otherwise
     */
    bool isGhostCell() const
    {
        return super::isGhostCell();
    }

    /**
     * \return process id
     */
    rank_type processId() const
    {
        return super::processId();
    }

    /**
     * \return process id
     */
    rank_type partition1() const
    {
        return super2::partition1( super::processId() );
    }

    /**
     * \return process id
     */
    rank_type partition2() const
    {
        return super2::partition2( super::processId() );
    }




    /**
     * \sa face()
     */
    edge_type const& edge( uint16_type i ) const
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;
        DCHECK( M_edges[i] != nullptr ) << "invalid edge (null pointer) for edge local id " << i << " in element " << this->id();
        return boost::cref( *M_edges[i] );
    }

    /**
     * \sa face()
     */
    edge_type& edge( uint16_type i )
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;
        DCHECK( M_edges[i] != nullptr ) << "invalid edge (null pointer) for edge local id " << i << " in element " << this->id();
        return boost::ref( *M_edges[i] );
    }

    edge_type & face( uint16_type i )
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;
        DCHECK( M_edges[i] != nullptr ) << "invalid edge (null pointer) for edge local id " << i << " in element " << this->id();

        return boost::ref( *M_edges[i] );
    }

    /**
     * \sa edge()
     */
    edge_type const& face( uint16_type i ) const
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;
        DCHECK( M_edges[i] != nullptr ) << "invalid edge (null pointer) for edge local id " << i << " in element " << this->id();

        return boost::cref( *M_edges[i] );
    }

    edge_type const* facePtr( uint16_type i ) const
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;
        DCHECK( M_edges[i] != nullptr ) << "invalid edge (null pointer) for edge local id " << i << " in element " << this->id();

        return M_edges[i];
    }

    /**
     * Inserts an edge.
     * \sa setEdge()
     */
    void setFace( uint16_type const i, edge_type const & p )
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;

        M_edges[i] = const_cast<edge_type*>( boost::addressof( p ) );
    }

    /**
     * \sa facePermutation(), permutation()
     */
    edge_permutation_type edgePermutation( uint16_type i ) const
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;

        return M_edge_permutation[i];
    }
    /**
     * \sa edgePermutation(), permutation()
     */
    edge_permutation_type facePermutation( uint16_type i ) const
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;

        return M_edge_permutation[i];
    }

    /**
     * \sa edgePermutation(), facePermutation()
     */
    edge_permutation_type permutation( uint16_type i ) const
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges << " in element id " << this->id();
        return M_edge_permutation[i];
    }

    /**
     * Inserts a point.  Uses point references put point
     * \sa setFace()
     */
    void setEdge( uint16_type i, edge_type const & p )
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;
        M_edges[i] = const_cast<edge_type*>( boost::addressof( p ) );
    }

    void setEdgePermutation( uint16_type i, edge_permutation_type o )
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;

        M_edge_permutation[i] = o;
    }

    typedef typename std::vector<edge_type*>::iterator face_iterator;
    typedef typename std::vector<edge_type*>::const_iterator face_const_iterator;

    /**
     * \return the iterator pair (begin,end) of faces
     */
    std::pair<face_iterator,face_iterator>
    faces()
    {
        return std::make_pair( M_edges.begin(), M_edges.end() );
    }

    /**
     * \return the iterator pair (begin,end) of faces
     */
    std::pair<face_const_iterator,face_const_iterator>
    faces() const
    {
        return std::make_pair( M_edges.begin(), M_edges.end() );
    }

    void disconnectSubEntities()
    {
        for(unsigned int i = 0; i<numLocalEdges;++i)
        {
            M_edges[i]->disconnect(*this);
        }
    }

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            DVLOG(2) << "Serializing Geoelement2D id: " << this->id() << "...\n";
            ar & boost::serialization::base_object<super>( *this );
            ar & boost::serialization::base_object<super2>( *this );
            ar & M_edges;
        }


private:

    std::vector<edge_type*> M_edges;
    std::vector<edge_permutation_type> M_edge_permutation;
};

/*-------------------------------------------------------------------------
  GeoElement2D
  --------------------------------------------------------------------------*/
template <uint16_type Dim, typename GEOSHAPE, typename SFO, typename T>
const uint16_type GeoElement2D<Dim, GEOSHAPE, SFO, T>::numLocalEdges;
template <uint16_type Dim, typename GEOSHAPE, typename SFO, typename T>
const uint16_type GeoElement2D<Dim, GEOSHAPE, SFO, T>::nDim;
template <uint16_type Dim, typename GEOSHAPE, typename SFO, typename T>
const uint16_type GeoElement2D<Dim, GEOSHAPE, SFO, T>::nRealDim;


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

    static const uint16_type nDim = Dim;

    typedef GeoND<Dim, GEOSHAPE, T, GeoElement0D<Dim, SubFaceOfNone, T> > super;
    typedef SubFaceOfNone super2;

    typedef GEOSHAPE GeoShape;

    typedef typename super::face_type entity_face_type;

    typedef GeoElement3D<Dim, GEOSHAPE,T> self_type;
    typedef self_type element_type;
    typedef GeoElement2D<Dim, entity_face_type, SubFaceOf<self_type>, T > face_type;
    typedef GeoElement1D<Dim, typename entity_face_type::topological_face_type, SubFaceOfMany<face_type>, T> edge_type;
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
        M_edges( numLocalEdges, nullptr ),
        M_faces( numLocalFaces ),
        M_edge_permutation( numLocalEdges, edge_permutation_type( edge_permutation_type::IDENTITY ) ),
        M_face_permutation( numLocalFaces )
    {
        std::fill( M_faces.begin(), M_faces.end(), ( face_type* )0 );

        std::fill( M_face_permutation.begin(), M_face_permutation.end(), face_permutation_type( face_permutation_type::IDENTITY ) );
    }

    /**
     * copy consttructor
     */
    GeoElement3D( GeoElement3D const& g )
        :
        super( g ),
        super2( g ),
        M_edges( g.M_edges ),
        M_faces( g.M_faces ),
        M_edge_permutation( g.M_edge_permutation ),
        M_face_permutation( g.M_face_permutation )
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
            M_edges = g.M_edges;
            M_faces = g.M_faces;
            M_edge_permutation = g.M_edge_permutation;
            M_face_permutation = g.M_face_permutation;
        }

        return *this;
    }

    //void setMesh( MeshBase const* m ) { super::setMesh( m ); }
    MeshBase const* mesh() const
    {
        return super::mesh();
    }
    /**
     * \return \c true if on the boundary, \c false otherwise
     */
    size_type id() const
    {
        return super::id();
    }

    Marker1 const& marker() const
    {
        return super::marker();
    }
    Marker1& marker()
    {
        return super::marker();
    }
    Marker2 const& marker2() const
    {
        return super::marker2();
    }
    Marker3 const& marker3() const
    {
        return super::marker3();
    }

    bool isGhostFace() const
    {
        return super2::isGhostFace( super::processId()  );
    }


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
    bool isOnBoundary() const
    {
        return super::isOnBoundary();
    }

    /**
     * \return maximum \c dimension of the sub-entity touching the boundary of the element
     */
    uint16_type boundaryEntityDimension() const
    {
        return super::boundaryEntityDimension();
    }

    /**
     * \return \c true if ghost cell, \c false otherwise
     */
    bool isGhostCell() const
    {
        return super::isGhostCell();
    }

    /**
     * \return process id
     */
    rank_type processId() const
    {
        return super::processId();
    }

    /**
     * \return process id
     */
    rank_type partition1() const
    {
        return super2::partition1( super::processId() );
    }

    /**
     * \return process id
     */
    rank_type partition2() const
    {
        return super2::partition2( super::processId() );
    }


    size_type ad_first() const
    {
        return invalid_size_type_value;
    }
    uint16_type pos_first() const
    {
        return invalid_uint16_type_value;
    }
    size_type ad_second() const
    {
        return invalid_size_type_value;
    }
    uint16_type pos_second() const
    {
        return invalid_uint16_type_value;
    }

    edge_type const& edge( uint16_type i ) const
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;
        DCHECK( M_edges[i] != nullptr ) << "invalid edge (null pointer) for edge local id " << i << " in element " << this->id();

        return *M_edges[i];
    }

    edge_type& edge( uint16_type i )
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;
        DCHECK( M_edges[i] != nullptr ) << "invalid edge (null pointer) for edge local id " << i << " in element " << this->id();

        return *M_edges[i];
    }

    edge_type const* edgePtr( uint16_type i ) const
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;
        DCHECK( M_edges[i] != nullptr ) << "invalid edge (null pointer) for edge local id " << i << " in element " << this->id();

        return M_edges[i];
    }

    edge_permutation_type edgePermutation( uint16_type i ) const
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;
        DCHECK( M_edges[i] != nullptr ) << "invalid edge (null pointer) for edge local id " << i << " in element " << this->id();


        return M_edge_permutation[i];
    }

    /**
     * Inserts an edge
     */
    void setEdge( uint16_type const i, edge_type const & p )
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;
        DCHECK( boost::addressof( p ) ) << "invalid edge (null pointer) for edge local id " << i << " in element " << this->id();
        M_edges[i] = const_cast<edge_type*>( boost::addressof( p ) );
    }

    void setEdgePermutation( uint16_type i, edge_permutation_type o )
    {
        DCHECK( i < numLocalEdges ) << "invalid local edge index " << i << " should be less than " << numLocalEdges ;

        M_edge_permutation[i] = o;
    }

    face_type const& face( uint16_type i ) const
    {
        FEELPP_ASSERT( i < numLocalFaces )( this->id() )( i ).error( "invalid local edge index" );
        FEELPP_ASSERT( M_faces[i] )( this->id() )( i ).error( "invalid edge (null pointer)" );
        return *M_faces[i];
    }

    face_type& face( uint16_type i )
    {
        FEELPP_ASSERT( i < numLocalFaces )( this->id() )( i ).error( "invalid local edge index" );
        FEELPP_ASSERT( M_faces[i] )( this->id() )( i ).error( "invalid edge (null pointer)" );
        return *M_faces[i];
    }

    face_type const* facePtr( uint16_type i ) const
    {
        FEELPP_ASSERT( i < numLocalFaces )( this->id() )( i ).error( "invalid local edge index" );
        //FEELPP_ASSERT( M_faces[i] )( i ).error( "invalid edge (null pointer)" );
        return M_faces[i];
    }

    face_permutation_type facePermutation( uint16_type i ) const
    {
        FEELPP_ASSERT( i < numLocalFaces )( this->id() )( i ).error( "invalid local face index" );
        FEELPP_ASSERT( M_faces[i] )( this->id() )( i ).error( "invalid face (null pointer)" );
        return M_face_permutation[i];
    }
    face_permutation_type permutation( uint16_type i ) const
    {
        FEELPP_ASSERT( i < numLocalFaces )( this->id() )( i ).error( "invalid local face index" );
        FEELPP_ASSERT( M_faces[i] )( this->id() )( i ).error( "invalid face (null pointer)" );
        return M_face_permutation[i];
    }

    /**
     * Inserts a face.
     */
    void setFace( uint16_type const i, face_type const & p )
    {
        M_faces[i] = const_cast<face_type*>( boost::addressof( p ) );
    }

    void setFacePermutation( uint16_type i, face_permutation_type o )
    {
        FEELPP_ASSERT( i < numLocalFaces )( this->id() )( i ).error( "invalid local face index" );
        M_face_permutation[i] = o;
    }

    typedef typename ublas::bounded_array<face_type*, numLocalFaces>::iterator face_iterator;
    typedef typename ublas::bounded_array<face_type*, numLocalFaces>::const_iterator face_const_iterator;

    /**
     * \return the iterator pair (begin,end) of faces
     */
    std::pair<face_iterator,face_iterator>
    faces()
    {
        return std::make_pair( M_faces.begin(), M_faces.end() );
    }

    /**
     * \return the iterator pair (begin,end) of faces
     */
    std::pair<face_const_iterator,face_const_iterator>
    faces() const
    {
        return std::make_pair( M_faces.begin(), M_faces.end() );
    }
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            ar & boost::serialization::base_object<super>( *this );
        }


private:

    ublas::bounded_array<edge_type*, numLocalEdges> M_edges;
    ublas::bounded_array<face_type*, numLocalFaces> M_faces;

    ublas::bounded_array<edge_permutation_type, numLocalEdges> M_edge_permutation;
    ublas::bounded_array<face_permutation_type, numLocalFaces> M_face_permutation;
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
template <uint16_type Dim, typename GEOSHAPE, typename T>
const uint16_type GeoElement3D<Dim, GEOSHAPE, T>::nDim;

} // Feel
#endif
