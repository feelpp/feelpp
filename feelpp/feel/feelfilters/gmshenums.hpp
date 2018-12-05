/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-07-15

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
   \file gmshenums.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-07-15
 */
#ifndef __GmshEnums_H
#define __GmshEnums_H 1

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wredeclared-class-member"
#endif
#include <boost/bimap.hpp>
#if defined(__clang__)
#pragma clang diagnostic pop
#endif
#include <boost/assign/list_of.hpp>
namespace Feel
{
enum GMSH_PARTITIONER
{
    GMSH_PARTITIONER_CHACO = 1,
    GMSH_PARTITIONER_METIS = 2

};
#if 1// defined( FEELPP_HAS_GMSH_H )
extern const GMSH_PARTITIONER GMSH_PARTITIONER_DEFAULT;


#endif

extern const char* FEELPP_GMSH_FORMAT_VERSION;

enum GMSH_ORDER
{
    GMSH_ORDER_ONE = 1,
    GMSH_ORDER_TWO = 2,
    GMSH_ORDER_THREE = 3,
    GMSH_ORDER_FOUR = 4,
    GMSH_ORDER_FIVE = 5
};

enum GMSH_FORMAT
{
    GMSH_FORMAT_ASCII = 0,
    GMSH_FORMAT_BINARY = 1
};


/**
 * \enum GMSH_ENTITY
 *
 * enumerate the various elements available in gmsh
 */
enum GMSH_ENTITY
{
    GMSH_LINE = 1, //!< Line (2 nodes).
    GMSH_TRIANGLE = 2, //!< Triangle (3 nodes).
    GMSH_QUADRANGLE = 3,  //!< Quadrangle (4 nodes).
    GMSH_TETRAHEDRON = 4, //!< Tetrahedron (4 nodes).
    GMSH_HEXAHEDRON = 5, //!< Hexahedron (8 nodes).
    GMSH_PRISM = 6,  //!< Prism (6 nodes).
    GMSH_PYRAMID = 7, //!< Pyramid (5 nodes).
    GMSH_LINE_2 = 8,  //!< Second order line (3 nodes: 2 associated
    //with the vertices and 1 with the edge).
    GMSH_TRIANGLE_2 = 9, //!< Second order triangle (6 nodes: 3
    //associated with the vertices and 3 with the
    //edges).
    GMSH_QUADRANGLE_2 = 10, //!<Second order quadrangle (9 nodes: 4
    //associated with the vertices, 4 with the
    //edges and 1 with the face).
    GMSH_TETRAHEDRON_2 = 11, //!< Second order tetrahedron (10 nodes:
    //4 associated with the vertices and 6
    //with the edges).
    GMSH_HEXAHEDRON_2 = 12,  //!< Second order hexahedron (27 nodes: 8
    //associated with the vertices, 12 with
    //the edges, 6 with the faces and 1 with
    //the volume).
    GMSH_PRISM_2 = 13, //!<Second order prism (18 nodes: 6 associated
    //with the vertices, 9 with the edges and 3 with
    //the quadrangular faces).
    GMSH_PYRAMID_2 = 14, //!<Second order pyramid (14 nodes: 5
    //associated with the vertices, 8 with the
    //edges and 1 with the quadrangular face).
    GMSH_POINT = 15, //!< Point (1 node).

    GMSH_TRIANGLE_INCOMPLETE_3=20, //!< triangle of order 3
    GMSH_TRIANGLE_3=21, //!< triangle of order 3
    GMSH_TRIANGLE_INCOMPLETE_4=22, //!< triangle of order 4
    GMSH_TRIANGLE_4=23, //!< triangle of order 4
    GMSH_TRIANGLE_INCOMPLETE_5=24, //!< triangle of order 5
    GMSH_TRIANGLE_5=25, //!< triangle of order 5

    GMSH_LINE_3=26, //!< line of order 3
    GMSH_LINE_4=27, //!< line of order 4
    GMSH_LINE_5=28, //!< line of order 5

    GMSH_TETRAHEDRON_3=29, //!< tetra of order 3
    GMSH_TETRAHEDRON_4=30, //!< tetra of order 4
    GMSH_TETRAHEDRON_5=31  //!< tetra of order 5
};

template<typename ConvexType>
class GmshOrdering
{
public:
    typedef boost::bimap<int,int> id_type;
    GmshOrdering();
    ~GmshOrdering() {}

    /**
     * \return the type of ids
     */
    int type() const
    {
        return M_type;
    }

    /**
     * \return the number of ids
     */
    int size() const
    {
        return M_id.size();
    }

    /**
     * \return the Feel id from the Gmsh id \p p
     */
    int fromGmshId( int p ) const
    {
        return M_id.right.at( p );
    }

    /**
     * \return the Gmsh id from the Feel id \p p
     */
    int toGmshId( int p ) const
    {
        return M_id.left.at( p );
    }

    /**
     * \return the Gmsh id from the Feel id \p p
     */
    int id( int p ) const
    {
        return M_id.left.at( p );
    }

private:

    int M_type;
    id_type M_id;
};

/// \cond detail
namespace detail
{
const int line_type[6] = { 0, GMSH_LINE, GMSH_LINE_2, GMSH_LINE_3, GMSH_LINE_4, GMSH_LINE_5 };
const int triangle_type[6] = { 0, GMSH_TRIANGLE, GMSH_TRIANGLE_2, GMSH_TRIANGLE_3, GMSH_TRIANGLE_4, GMSH_TRIANGLE_5 };
const int tetrahedron_type[6] = { 0, GMSH_TETRAHEDRON, GMSH_TETRAHEDRON_2, GMSH_TETRAHEDRON_3, GMSH_TETRAHEDRON_4, GMSH_TETRAHEDRON_5 };



const int quad_type[6] = { 0, GMSH_QUADRANGLE, GMSH_QUADRANGLE_2, 0, 0, 0 };
const int hexa_type[6] = { 0, GMSH_HEXAHEDRON, GMSH_HEXAHEDRON_2, 0, 0, 0 };

//int gmshquadtype[5] = { 0, 4, 10, 0, 0 };
}
/// \endcond detail

template<typename ConvexType>
GmshOrdering<ConvexType>::GmshOrdering()
{
    using namespace boost::assign;
    typedef typename id_type::relation relation;

    if ( ConvexType::nDim == 0 )
    {
        M_type = GMSH_POINT;

        for ( int i = 0; i < ConvexType::numPoints; ++i )
            M_id.insert( id_type::value_type( i, i ) );
    }

    else if ( ConvexType::is_simplex )
    {
        if ( ConvexType::nDim == 1 )
        {
            M_type = Feel::detail::line_type[ConvexType::nOrder];

            for ( int i = 0; i < ConvexType::numPoints; ++i )
                M_id.insert( id_type::value_type( i, i ) );
        }

        if ( ConvexType::nDim == 2 )
        {
            M_type = Feel::detail::triangle_type[ConvexType::nOrder];

            if ( ConvexType::nOrder == 1 )
                M_id = list_of<relation>( 0,0 )( 1,1 )( 2,2 );

            //M_id+=0,1,2;
            if ( ConvexType::nOrder == 2 )
                M_id = list_of<relation>( 0,0 )( 1,1 )( 2,2 )( 3,4 )( 4,5 )( 5,3 );

            //M_id += 0,1,2,5,3,4;
            if ( ConvexType::nOrder == 3 )
                M_id = list_of<relation>( 0,0 )( 1,1 )( 2,2 )( 3,5 )( 4,6 )( 5,7 )( 6,8 )( 7,3 )( 8,4 )( 9,9 );

            //M_id += 0,1,2,7,8,3,4,5,6,9;
            if ( ConvexType::nOrder == 4 )
                M_id = list_of<relation>( 0,0 )( 1,1 )( 2,2 )( 3,6 )( 4,7 )( 5,8 )( 6,9 )( 7,10 )( 8,11 )( 9,3 )( 10,4 )( 11,5 )( 12,12 )( 13,13 )( 14,14 );

            //M_id += 0,1,2,9,10,11,3,4,5,6,7,8,12,13,14;
            if ( ConvexType::nOrder == 5 )
                M_id = list_of<relation>( 0,0 )( 1,1 )( 2,2 )( 3,7 )( 4,8 )( 5,9 )( 6,10 )( 7,11 )( 8,12 )( 9,13 )( 10,14 )( 11,3 )( 12,4 )( 13,5 )( 14,6 )( 15,15 )( 16,18 )( 17,16 )( 18,20 )( 19,19 )( 20,17 );

            //M_id += 0,1,2,11,12,13,14,3,4,5,6,7,8,9,10,15,16,17,19,20,18;
        }

        if ( ConvexType::nDim == 3 )
        {
            M_type = Feel::detail::tetrahedron_type[ConvexType::nOrder];

            if ( ConvexType::nOrder == 1 )
                M_id = list_of<relation>
                       ( 0,0 )( 1,1 )( 2,2 )( 3,3 ) // vertices
                       ;

            //M_id+=0,1,2,3;
            if ( ConvexType::nOrder == 2 )
                M_id = list_of<relation>
                       ( 0,0 )( 1,1 )( 2,2 )( 3,3 ) // vertices
                       ( 4,5 ) // edge 0
                       ( 5,6 ) // edge 1
                       ( 6,4 ) // edge 2
                       ( 7,7 ) // edge 3
                       ( 8,9 ) // edge 4
                       ( 9,8 ) // edge 5
                       ;

            //M_id+=0,1,2,3,6,4,5,7,9,8;
            if ( ConvexType::nOrder == 3 )
                M_id = list_of<relation>
                       ( 0,0 )( 1,1 )( 2,2 )( 3,3 ) // vertices
                       ( 4,6 )( 5,7 ) // edge 0
                       ( 6,8 )( 7,9 ) // edge 1
                       ( 8,4 )( 9,5 ) // edge 2
                       ( 10,11 )( 11,10 ) // edge 3
                       ( 12,15 )( 13,14 ) // edge 4
                       ( 14,13 )( 15,12 ) // edge 5
                       ( 16,19 )      // face 0
                       ( 17,18 )      // face 1
                       ( 18,17 )      // face 2
                       ( 19,16 )      // face 3
                       ;

            if ( ConvexType::nOrder == 4 )
                M_id = list_of<relation>
                       ( 0,0 )( 1,1 )( 2,2 )( 3,3 ) // vertices
                       ( 4,7 )( 5,8 )( 6,9 ) // edge 0
                       ( 7,10 )( 8,11 )( 9,12 ) // edge 1
                       ( 10,4 )( 11,5 )( 12,6 ) // edge 2
                       ( 13,15 )( 14,14 )( 15,13 ) // edge 3
                       ( 16,21 )( 17,20 )( 18,19 ) // edge 4
                       ( 19,18 )( 20,17 )( 21,16 ) // edge 5
                       ( 22,32 )( 23,33 )( 24,31 ) // face 0
                       ( 25,28 )( 26,30 )( 27,29 ) // face 1
                       ( 28,25 )( 29,26 )( 30,27 ) // face 2
                       ( 31,22 )( 32,24 )( 33,23 ) // face 3
                       ( 34,34 )             // interior point
                       ;

            if ( ConvexType::nOrder > 4 )
                for ( int i = 0; i < ConvexType::numPoints; ++i )
                    M_id.insert( id_type::value_type( i, i ) );

        }
    }

    else
    {

        if ( ConvexType::nDim == 1 )
        {
            M_type = Feel::detail::line_type[ConvexType::nOrder];

            for ( int i = 0; i < ConvexType::numPoints; ++i )
                M_id.insert( id_type::value_type( i, i ) );
        }

        if ( ConvexType::nDim == 2 )
        {
            M_type = Feel::detail::quad_type[ConvexType::nOrder];

            if ( ConvexType::nOrder == 1 )
                M_id = list_of<relation>( 0,0 )( 1,1 )( 2,2 )( 3,3 );

            //M_id+=0,1,2,3;
            if ( ConvexType::nOrder == 2 )
                M_id = list_of<relation>( 0,0 )( 1,1 )( 2,2 )( 3,3 ) // vertices
                       ( 4,4 ) // edge 0
                       ( 5,5 ) // edge 1
                       ( 6,6 ) // edge 2
                       ( 7,7 ) // edge 3
                       ( 8,8 ) // face 0
                       ;
        }

        if ( ConvexType::nDim == 3 )
        {
            M_type = Feel::detail::hexa_type[ConvexType::nOrder];

            if ( ConvexType::nOrder == 1 )
                M_id = list_of<relation>( 0,0 )( 1,1 )( 2,2 )( 3,3 )( 4,4 )( 5,5 )( 6,6 )( 7,7 );

            if ( ConvexType::nOrder == 2 )
                M_id = list_of<relation>
                       ( 0,0 )( 1,1 )( 2,2 )( 3,3 )( 4,4 )( 5,5 )( 6,6 )( 7,7 ) // vertices
                       ( 8,8 ) // edge 0
                       ( 9,11 ) // edge 1
                       ( 10,13 ) // edge 2
                       ( 11,9 ) // edge 3
                       ( 12,12 ) // edge 4
                       ( 13,16 ) // edge 5
                       ( 14,10 ) // edge 6
                       ( 15,14 ) // edge 7
                       ( 16,18 ) // edge 8
                       ( 17,15 ) // edge 9
                       ( 18,19 ) // edge 10
                       ( 19,17 ) // edge 11
                       ( 20,20 ) // face 0
                       ( 21,21 ) // face 1
                       ( 22,23 ) // face 2
                       ( 23,24 ) // face 3
                       ( 24,22 ) // face 4
                       ( 25,25 ) // face 5
                       ( 26,26 ) // volume 0
                       ;

            //M_id+=0,1,2,3,4,5,6,7;
        }

    }
}

} // Feel

#if !defined( FEELPP_HAS_GMSH_H )
// Element types in .msh file format (numbers should not be changed)
#define MSH_LIN_2    1
#define MSH_TRI_3    2
#define MSH_QUA_4    3
#define MSH_TET_4    4
#define MSH_HEX_8    5
#define MSH_PRI_6    6
#define MSH_PYR_5    7
#define MSH_LIN_3    8
#define MSH_TRI_6    9
#define MSH_QUA_9    10
#define MSH_TET_10   11
#define MSH_HEX_27   12
#define MSH_PRI_18   13
#define MSH_PYR_14   14
#define MSH_PNT      15
#define MSH_QUA_8    16
#define MSH_HEX_20   17
#define MSH_PRI_15   18
#define MSH_PYR_13   19
#define MSH_TRI_9    20
#define MSH_TRI_10   21
#define MSH_TRI_12   22
#define MSH_TRI_15   23
#define MSH_TRI_15I  24
#define MSH_TRI_21   25
#define MSH_LIN_4    26
#define MSH_LIN_5    27
#define MSH_LIN_6    28
#define MSH_TET_20   29
#define MSH_TET_35   30
#define MSH_TET_56   31
#define MSH_TET_22   32
#define MSH_TET_28   33
#define MSH_POLYG_   34
#define MSH_POLYH_   35
#define MSH_QUA_16   36
#define MSH_QUA_25   37
#define MSH_QUA_36   38
#define MSH_QUA_12   39
#define MSH_QUA_16I  40
#define MSH_QUA_20   41
#define MSH_TRI_28   42
#define MSH_TRI_36   43
#define MSH_TRI_45   44
#define MSH_TRI_55   45
#define MSH_TRI_66   46
#define MSH_QUA_49   47
#define MSH_QUA_64   48
#define MSH_QUA_81   49
#define MSH_QUA_100  50
#define MSH_QUA_121  51
#define MSH_TRI_18   52
#define MSH_TRI_21I  53
#define MSH_TRI_24   54
#define MSH_TRI_27   55
#define MSH_TRI_30   56
#define MSH_QUA_24   57
#define MSH_QUA_28   58
#define MSH_QUA_32   59
#define MSH_QUA_36I  60
#define MSH_QUA_40   61
#define MSH_LIN_7    62
#define MSH_LIN_8    63
#define MSH_LIN_9    64
#define MSH_LIN_10   65
#define MSH_LIN_11   66
#define MSH_LIN_B    67
#define MSH_TRI_B    68
#define MSH_POLYG_B  69
#define MSH_LIN_C    70
// TETS COMPLETE (6->10)
#define MSH_TET_84   71
#define MSH_TET_120  72
#define MSH_TET_165  73
#define MSH_TET_220  74
#define MSH_TET_286  75
// TETS INCOMPLETE (6->10)
#define MSH_TET_34   79
#define MSH_TET_40   80
#define MSH_TET_46   81
#define MSH_TET_52   82
#define MSH_TET_58   83
//
#define MSH_LIN_1    84
#define MSH_TRI_1    85
#define MSH_QUA_1    86
#define MSH_TET_1    87
#define MSH_HEX_1    88
#define MSH_PRI_1    89
#define MSH_PRI_40   90
#define MSH_PRI_75   91
// HEXES COMPLETE (3->9)
#define MSH_HEX_64   92
#define MSH_HEX_125  93
#define MSH_HEX_216  94
#define MSH_HEX_343  95
#define MSH_HEX_512  96
#define MSH_HEX_729  97
#define MSH_HEX_1000 98
// HEXES INCOMPLETE (3->9)
#define MSH_HEX_32   99
#define MSH_HEX_44  100
#define MSH_HEX_56  101
#define MSH_HEX_68  102
#define MSH_HEX_80  103
#define MSH_HEX_92  104
#define MSH_HEX_104 105
// PRISMS COMPLETE (5->9)
#define MSH_PRI_126 106
#define MSH_PRI_196 107
#define MSH_PRI_288 108
#define MSH_PRI_405 109
#define MSH_PRI_550 110
// PRISMS INCOMPLETE (3->9)
#define MSH_PRI_24  111
#define MSH_PRI_33  112
#define MSH_PRI_42  113
#define MSH_PRI_51  114
#define MSH_PRI_60  115
#define MSH_PRI_69  116
#define MSH_PRI_78  117
// PYRAMIDS COMPLETE (3->9)
#define MSH_PYR_30  118
#define MSH_PYR_55  119
#define MSH_PYR_91  120
#define MSH_PYR_140 121
#define MSH_PYR_204 122
#define MSH_PYR_285 123
#define MSH_PYR_385 124
// PYRAMIDS INCOMPLETE (3->9)
#define MSH_PYR_21  125
#define MSH_PYR_29  126
#define MSH_PYR_37  127
#define MSH_PYR_45  128
#define MSH_PYR_53  129
#define MSH_PYR_61  130
#define MSH_PYR_69  131
// Additional types
#define MSH_PYR_1   132
#define MSH_PNT_SUB 133
#define MSH_LIN_SUB 134
#define MSH_TRI_SUB 135
#define MSH_TET_SUB 136
#define MSH_TET_16  137
#define MSH_TRI_MINI 138
#define MSH_TET_MINI 139
#define MSH_TRIH_4   140
#define MSH_NUM_TYPE 140
// Element types
#define TYPE_PNT     1
#define TYPE_LIN     2
#define TYPE_TRI     3
#define TYPE_QUA     4
#define TYPE_TET     5
#define TYPE_PYR     6
#define TYPE_PRI     7
#define TYPE_HEX     8
#define TYPE_POLYG   9
#define TYPE_POLYH   10
#define TYPE_XFEM    11
#define TYPE_MINI    12
#define TYPE_TRIH    13

#endif
#endif /* __GmshEnums_H */
