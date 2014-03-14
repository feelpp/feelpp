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
/*! file GeoEntity.h */
#ifndef _GEOENTITY_HH_
#define _GEOENTITY_HH_

#include <vector>

#include <boost/version.hpp>
#if (BOOST_VERSION >= 103400)
#include <boost/none.hpp>
#else
#include <boost/none_t.hpp>
#endif /* BOOST_VERSION >= 103400 */

#include <boost/mpl/vector.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/operators.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/context.hpp>

namespace Feel
{

enum FaceLocation { INTERNAL = false, ON_BOUNDARY = true };

/**
 * \enum GeoEntityContext
 */
enum GeoEntityContext
{
    MESH_ENTITY_INTERNAL          = ( 1<<0 ), /**< internal entity */
    MESH_ENTITY_BOUNDARY          = ( 1<<1 )  /**< boundary entity */
};

/**
 * \enum ReferenceGeometry
 *
 */
enum ReferenceGeometry
{
    GEOMETRY_POINT    = ( 1<<0 ), /**< point entity */
    GEOMETRY_LINE     = ( 1<<1 ), /**< line entity */
    GEOMETRY_SURFACE  = ( 1<<2 ), /**< surface entity */
    GEOMETRY_VOLUME   = ( 1<<3 ), /**< volume entity */
    GEOMETRY_4        = ( 1<<4 ), /**< hypercube entity */
    GEOMETRY_5        = ( 1<<5 )  /**< hypercube entity */
};

/**
 * \enum ReferenceShapes
 *
 */
enum ReferenceShapes
{
    SHAPE_LINEAR   = ( 1<<0 ),
    SHAPE_BILINEAR = ( 1<<1 ),
    SHAPE_QUADRATIC= ( 1<<2 ),
    SHAPE_NONE     = ( 1<<3 ),
    SHAPE_POINT    = ( 1<<4 ),
    SHAPE_LINE     = ( 1<<5 ),
    SHAPE_TRIANGLE = ( 1<<6 ),
    SHAPE_QUAD     = ( 1<<7 ),
    SHAPE_HEXA     = ( 1<<8 ),
    SHAPE_PRISM    = ( 1<<9 ),
    SHAPE_TETRA    = ( 1<<10 ),
    SHAPE_SP4      = ( 1<<11 ),
    SHAPE_SP5      = ( 1<<12 )
};

/**
 * @defgroup GeoEntites Basis Reference Shapes
 \ingroup Obsolet_Groups */
//@{

/// \cond detail
/**
 * Represents the range of entities of topological dimension d for a given Shape
 */
template<typename E>
class EntityRange
{
    template<uint16_type td>
    struct num
    {
        typedef typename mpl::if_<mpl::equal_to<mpl::int_<td>, mpl::int_<0> >,
                mpl::int_<E::numVertices>,
                typename mpl::if_<mpl::equal_to<mpl::int_<td>, mpl::int_<1> >,
                mpl::int_<E::numEdges>,
                typename mpl::if_<mpl::equal_to<mpl::int_<td>, mpl::int_<2> >,
                mpl::int_<E::numGeometricFaces>,
                mpl::int_<E::numVolumes>
                >::type // int_<2>
                >::type // int_<1>
                >::type type;
        static const uint16_type value = type::value;
    };
    uint16_type d;
public:
    EntityRange( uint16_type td = 0 )
        :
        d( td )
    {
        check_invariant();
    }
    EntityRange( EntityRange const& r ) : d( r.d ) {}
    ~EntityRange() {}

    uint16_type topologicalDimension() const
    {
        return d;
    }

    void setTopologicalDimension( uint16_type td )
    {
        d = td;
        check_invariant();
    }

    uint16_type begin() const
    {
        return 0;
    }
    uint16_type end() const
    {
        if ( d == 0 )
            return num<0>::value;

        if ( d == 1 )
            return num<1>::value;

        if ( d == 2 )
            return num<2>::value;

        if ( d == 3 )
            return num<3>::value;

        return num<0>::value;
    }
private:
    void check_invariant()
    {
        //FEELPP_ASSERT(  d <= E::topological_dimension )( d )( E::topological_dimension ).error( "invalid topological dimension" );
    }
};


namespace details
{
template<int N, int P>
struct pow
{
    static const size_type value = N*pow<N, P-1>::value;
};
template<int N>
struct pow<N, 0>
{
    static const size_type value = 1;
};
}

/**
 * \typedef no_permutation_type
 *
 * defines an enum for the possible edge permutation with respect to
 * its parent elements
 */
struct no_permutation: public boost::detail::identifier<uint16_type, no_permutation>
{
    static const uint16_type NO_PERMUTATION       = 0; //!< no permutation
    static const uint16_type IDENTITY             = 1; //!< no permutation
    static const uint16_type N_PERMUTATIONS       = 2; //!< number of permutations
    typedef boost::detail::identifier<uint16_type, no_permutation> super;
    typedef super::value_type value_type;
    no_permutation()                           : super( IDENTITY ) {}
    explicit no_permutation( value_type v )    : super( v ) {}
    no_permutation & operator=( value_type v )
    {
        this->assign( v );
        return *this;
    }
    no_permutation& operator++()
    {
        this->assign( this->value()+1 );
        return *this;
    }
};


/**
 * \typedef edge_permutation_type
 *
 * defines an enum for the possible edge permutation with respect to
 * its parent elements
 */
struct line_permutations: public boost::detail::identifier<uint16_type, line_permutations>
{
    static const uint16_type NO_PERMUTATION       = 0; //!< the edge has no permutation defined
    static const uint16_type IDENTITY             = 1; //!< the edge is oriented according to the standard one
    static const uint16_type REVERSE_PERMUTATION  = 2; //!< the edge is orientated reversely wrt parent element
    static const uint16_type N_PERMUTATIONS       = 3; //!< number of permutations
    typedef boost::detail::identifier<uint16_type, line_permutations> super;
    typedef super::value_type value_type;
    line_permutations()                           : super( IDENTITY ) {}
    explicit line_permutations( value_type v )    : super( v ) {}
    operator value_type() const { return this->value(); }
    line_permutations& operator=( value_type v )
    {
        this->assign( v );
        return *this;
    }
    line_permutations& operator++()
    {
        this->assign( this->value()+1 );
        return *this;
    }
};

enum line_permutations_dummy {};

struct triangular_faces_type: public boost::detail::identifier<uint16_type, triangular_faces_type>
{
    static const uint16_type NO_PERMUTATION       = 0; //!< the face has no permutation associated
    static const uint16_type IDENTITY             = 1; //!< the face has the identity permutation associated (standard one)
    static const uint16_type REVERSE_HEIGHT       = 2; //!< reflection according to the height direction
    static const uint16_type REVERSE_BASE         = 3; //!< reflection according to the base direction
    static const uint16_type REVERSE_HYPOTENUSE   = 4; //!< reflection according to the hypotenuse direction
    static const uint16_type ROTATION_ANTICLOCK   = 5; //!< rotates the points once in the anticlockwise sense
    static const uint16_type ROTATION_CLOCKWISE   = 6; //!< rotates the points once in the clockwise sense
    static const uint16_type N_PERMUTATIONS       = 7; //!< number of permutations

    static const uint16_type PRINCIPAL_DIAGONAL       = 4; //!< reflection according to the principal diagonal
    static const uint16_type SECOND_DIAGONAL          = 7; //!< reflection according to the second diagonal
    static const uint16_type ROTATION_TWICE_CLOCKWISE = 8; //!< rotates the points twice in the (anti)clockwise sense

    typedef boost::detail::identifier<uint16_type, triangular_faces_type> super;
    typedef super::value_type value_type;
    triangular_faces_type()                           : super( IDENTITY ) {}
    explicit triangular_faces_type( value_type v )    : super( v ) {}
    ///triangular_faces_type( triangular_faces_type const& tft ) : super( tft ) { value( tft.value() ); }
    triangular_faces_type & operator=( value_type v )
    {
        this->assign( v );
        return *this;
    }
    triangular_faces_type& operator++()
    {
        this->assign( this->value()+1 );
        return *this;
    }
};

struct quadrangular_faces: public boost::detail::identifier<uint16_type, quadrangular_faces>
{
    static const uint16_type NO_PERMUTATION           = 0; //!< the face has no permutation associated
    static const uint16_type IDENTITY                 = 1; //!< the face has the identity permutation associated (standard one)
    static const uint16_type REVERSE_HEIGHT           = 2; //!< reflection according to the height direction
    static const uint16_type REVERSE_BASE             = 3; //!< reflection according to the base direction
    static const uint16_type PRINCIPAL_DIAGONAL       = 4; //!< reflection according to the principal diagonal
    static const uint16_type ROTATION_ANTICLOCK       = 5; //!< rotates the points once in the anticlockwise sense
    static const uint16_type ROTATION_CLOCKWISE       = 6; //!< rotates the points once in the clockwise sense
    static const uint16_type SECOND_DIAGONAL          = 7; //!< reflection according to the second diagonal
    static const uint16_type ROTATION_TWICE_CLOCKWISE = 8; //!< rotates the points twice in the (anti)clockwise sense

    static const uint16_type N_PERMUTATIONS           = 9; //!< number of permutations

    static const uint16_type REVERSE_HYPOTENUSE   = 4; //!< reflection according to the hypotenuse direction

    typedef boost::detail::identifier<uint16_type, quadrangular_faces> super;
    typedef super::value_type value_type;
    quadrangular_faces()                           : super( IDENTITY ) {}
    explicit quadrangular_faces( value_type v )    : super( v ) {}
    //quadrangular_faces( quadrangular_faces const& tft ) : super( tft ) { value( tft.value() ); }
    quadrangular_faces & operator=( value_type v )
    {
        this->assign( v );
        return *this;
    }
    quadrangular_faces& operator++()
    {
        this->assign( this->value()+1 );
        return *this;
    }
};
/// \endcond

} // Feel
#endif
