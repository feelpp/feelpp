/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2010-07-15
 */
#ifndef __GmshEnums_H
#define __GmshEnums_H 1

#include <boost/bimap.hpp>
#include <boost/assign/list_of.hpp>
namespace Life
{
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
    GMSH_TETRAHEDRON_5=31, //!< tetra of order 5
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
    int type() const { return M_type; }

    /**
     * \return the number of ids
     */
    int size() const { return M_id.size(); }

    /**
     * \return the Life id from the Gmsh id \p p
     */
    int fromGmshId( int p ) const { return M_id.right.at(p); }

    /**
     * \return the Gmsh id from the Life id \p p
     */
    int toGmshId( int p ) const { return M_id.left.at(p); }

    /**
     * \return the Gmsh id from the Life id \p p
     */
    int id( int p ) const { return M_id.left.at(p); }

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
        for( int i = 0; i < ConvexType::numPoints; ++i )
            M_id.insert( id_type::value_type( i, i ) );
    }
    else if ( ConvexType::is_simplex )
    {
        if ( ConvexType::nDim == 1 )
        {
            M_type = detail::line_type[ConvexType::nOrder];
            for( int i = 0; i < ConvexType::numPoints; ++i )
                M_id.insert( id_type::value_type( i, i ) );
        }
        if ( ConvexType::nDim == 2 )
        {
            M_type = detail::triangle_type[ConvexType::nOrder];
            if ( ConvexType::nOrder == 1 )
                M_id = list_of<relation>(0,0)(1,1)(2,2);
            //M_id+=0,1,2;
            if ( ConvexType::nOrder == 2 )
                M_id = list_of<relation>(0,0)(1,1)(2,2)(3,4)(4,5)(5,3);
            //M_id += 0,1,2,5,3,4;
            if ( ConvexType::nOrder == 3 )
                M_id = list_of<relation>(0,0)(1,1)(2,2)(3,5)(4,6)(5,7)(6,8)(7,3)(8,4)(9,9);
            //M_id += 0,1,2,7,8,3,4,5,6,9;
            if ( ConvexType::nOrder == 4 )
                M_id = list_of<relation>(0,0)(1,1)(2,2)(3,6)(4,7)(5,8)(6,9)(7,10)(8,11)(9,3)(10,4)(11,5)(12,12)(13,13)(14,14);
            //M_id += 0,1,2,9,10,11,3,4,5,6,7,8,12,13,14;
            if ( ConvexType::nOrder == 5 )
                M_id = list_of<relation>(0,0)(1,1)(2,2)(3,11)(4,12)(5,13)(6,14)(7,3)(8,4)(9,5)(10,6)(11,7)(12,8)(13,9)(14,10)(15,15)(16,16)(17,17)(18,19)(19,20)(20,18);
            //M_id += 0,1,2,11,12,13,14,3,4,5,6,7,8,9,10,15,16,17,19,20,18;
        }
        if ( ConvexType::nDim == 3 )
        {
            M_type = detail::tetrahedron_type[ConvexType::nOrder];
            if ( ConvexType::nOrder == 1 )
                M_id = list_of<relation>
                    (0,0)(1,1)(2,2)(3,3) // vertices
                    ;
            //M_id+=0,1,2,3;
            if ( ConvexType::nOrder == 2 )
                M_id = list_of<relation>
                    (0,0)(1,1)(2,2)(3,3)  // vertices
                    (4,5) // edge 0
                    (5,6) // edge 1
                    (6,4) // edge 2
                    (7,7) // edge 3
                    (8,8) // edge 4
                    (9,9) // edge 5
                    ;
            //M_id+=0,1,2,3,6,4,5,7,9,8;
            if ( ConvexType::nOrder == 3 )
                M_id = list_of<relation>
                    (0,0)(1,1)(2,2)(3,3) // vertices
                    (4,6)(5,7)     // edge 0
                    (6,8)(7,9)    // edge 1
                    (8,4)(9,5)     // edge 2
                    (10,10)(11,11) // edge 3
                    (12,12)(13,13) // edge 4
                    (14,14)(15,15) // edge 5
                    (16,18)        // face 0
                    (17,17)        // face 1
                    (18,19)        // face 2
                    (19,16)        // face 3
                    ;
            if ( ConvexType::nOrder == 4 )
                M_id = list_of<relation>
                    (0,0)(1,1)(2,2)(3,3) // vertices
                    (4,7)(5,8)(6,9)      // edge 0
                    (7,10)(8,11)(9,12)   // edge 1
                    (10,4)(11,5)(12,6)   // edge 2
                    (13,13)(14,14)(15,15)// edge 3
                    (16,16)(17,17)(18,18)// edge 4
                    (19,19)(20,20)(21,21)// edge 5
                    (22,30)(23,28)(24,29)// face 0
                    (25,27)(26,25)(27,26)// face 1
                    (28,33)(29,32)(30,31)// face 2
                    (31,22)(32,23)(33,24)// face 4
                    (34,34)              // interior point
                    ;
            if ( ConvexType::nOrder > 4 )
                for( int i = 0; i < ConvexType::numPoints; ++i )
                    M_id.insert( id_type::value_type( i, i ) );

        }
    }
    else
    {

        if ( ConvexType::nDim == 1 )
        {
            M_type = detail::line_type[ConvexType::nOrder];
            for( int i = 0; i < ConvexType::numPoints; ++i )
                M_id.insert( id_type::value_type( i, i ) );
        }
        if ( ConvexType::nDim == 2 )
        {
            M_type = detail::quad_type[ConvexType::nOrder];
            if ( ConvexType::nOrder == 1 )
                M_id = list_of<relation>(0,0)(1,1)(2,2)(3,3);
            //M_id+=0,1,2,3;
            if ( ConvexType::nOrder == 2 )
                M_id = list_of<relation>(0,0)(1,1)(2,2)(3,3)(4,7)(5,4)(6,5)(7,6);
            //M_id += 0,1,2,3,7,4,5,6;
        }
        if ( ConvexType::nDim == 3 )
        {
            M_type = detail::hexa_type[ConvexType::nOrder];
            if ( ConvexType::nOrder == 1 )
                M_id = list_of<relation>(0,0)(1,1)(2,2)(3,3)(4,4)(5,5)(6,6)(7,7);
            //M_id+=0,1,2,3,4,5,6,7;
        }

    }
    std::cout << "There are " << M_id.size() << "relations" << std::endl;

    for( auto iter = M_id.begin(), iend = M_id.end(); iter != iend; ++iter )
    {
        // iter->left  : data : int
        // iter->right : data : std::string

        std::cout << iter->left << " <--> " << iter->right << std::endl;
    }
}

} // Life

#endif /* __GmshEnums_H */
