/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@imag.fr>
       Date: 2011-03-03

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file geotool.hpp
   \author Vincent Chabannes <vincent.chabannes@imag.fr>
   \date 2011-03-03
 */


#ifndef FEELPP_GEOTOOL_HPP
#define FEELPP_GEOTOOL_HPP 1

#include <iostream>
#include <string>
#include <sstream>
#include <list>
#include <map>

#include <boost/preprocessor/tuple/elem.hpp>

#include <feel/feelalg/glas.hpp>
//#include <boost/parameter/keyword.hpp>
//#include <boost/parameter/preprocessor.hpp>
#include <feel/feelcore/parameter.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/straightenmesh.hpp>
#include <feel/feelfilters/importergmsh.hpp>
#include <feel/feelfilters/detail/mesh.hpp>





/*_________________________________________________*/
/*_________________________________________________*/
/*_________________________________________________*/

# define GEOTOOL_SHAPE                                                  \
    ( 19, ( ( Line          , 1, 0, 0, "line"         , 2, LINE        ), \
            ( Triangle      , 2, 1, 0, "triangle"     , 3, TRIANGLE    ), \
            ( Rectangle     , 2, 1, 0, "rectangle"    , 2, RECTANGLE   ), \
            ( Quadrangle    , 2, 1, 0, "quadrangle"   , 4, QUADRANGLE  ), \
            ( Pentagon      , 2, 1, 0, "pentagon"     , 5, PENTAGON    ), \
            ( Hexagon       , 2, 1, 0, "hexagon"      , 6, HEXAGON     ), \
            ( Circle        , 2, 1, 0, "circle"       , 2, CIRCLE      ), \
            ( Ellipse       , 2, 1, 0, "ellipse"      , 3, ELLIPSE     ), \
            ( Pie           , 2, 1, 0, "pie"          , 3, PIE         ), \
            ( Special_1a    , 2, 2, 0, "special_1a"   , 1, SPECIAL_1A  ), \
            ( Special_1b    , 2, 1, 0, "special_1b"   , 1, SPECIAL_1B  ), \
            ( Peanut        , 2, 1, 0, "peanut"       , 4, PEANUT      ), \
            ( Tetrahedron   , 3, 4, 1, "tetrahedron"  , 4, TETRAHEDRON ), \
            ( Hexahedron    , 3, 6, 1, "hexahedron"   , 8, HEXAHEDRON  ), \
            ( Cube          , 3, 6, 1, "cube"         , 2, CUBE        ), \
            ( Cylindre      , 3, 6, 1, "cylindre"     , 4, CYLINDRE    ), \
            ( Sphere        , 3, 8, 1, "sphere"       , 2, SPHERE      ), \
            ( Tube          , 3,20, 4, "tube"         , 5, TUBE        ), \
            ( Special3D_1   , 3, 9, 1, "special3D_1"  , 1, SPECIAL3D_1 ) \
            )                                                           \
      )                                                                 \
    /**/

/*_________________________________________________*/

# define GEOTOOL_MARKER_POINT_LINE              \
    ( 2, ( ( 1, 1, ( 1 ) ),                     \
           ( 2, 1, ( 2 ) ) )                    \
      )                                         \
    /**/

# define GEOTOOL_MARKER_LINE_LINE               \
    ( 1, ( ( 1, 1, ( 1 ) ) )                    \
      )                                         \
    /**/

/*_________________________________________________*/

# define GEOTOOL_MARKER_POINT_TRIANGLE          \
    ( 3, ( ( 1, 1, ( 1 ) ),                     \
           ( 2, 1, ( 2 ) ),                     \
           ( 3, 1, ( 3 ) ) )                    \
      )                                         \
    /**/

# define GEOTOOL_MARKER_LINE_TRIANGLE           \
    ( 3, ( ( 1, 1, ( 1 ) ),                     \
           ( 2, 1, ( 2 ) ),                     \
           ( 3, 1, ( 3 ) ) )                    \
      )                                         \
    /**/
# define GEOTOOL_MARKER_SURFACE_TRIANGLE        \
    ( 1, ( ( 1, 1, ( 1 ) ) )                    \
      )                                         \
    /**/

/*_________________________________________________*/

# define GEOTOOL_MARKER_POINT_RECTANGLE         \
    ( 4, ( ( 1, 1, ( 1 ) ),                     \
           ( 2, 1, ( 2 ) ),                     \
           ( 3, 1, ( 3 ) ),                     \
           ( 4, 1, ( 4 ) ) )                    \
      )                                         \
    /**/
# define GEOTOOL_MARKER_LINE_RECTANGLE          \
    ( 4, ( ( 1, 1, ( 1 ) ),                     \
           ( 2, 1, ( 2 ) ),                     \
           ( 3, 1, ( 3 ) ),                     \
           ( 4, 1, ( 4 ) ) )                    \
      )                                         \
    /**/
# define GEOTOOL_MARKER_SURFACE_RECTANGLE       \
    ( 1, ( ( 1, 1, ( 1 ) ) )                    \
      )                                         \
    /**/

/*_________________________________________________*/

# define GEOTOOL_MARKER_POINT_QUADRANGLE        \
    ( 4, ( ( 1, 1, ( 1 ) ),                     \
           ( 2, 1, ( 2 ) ),                     \
           ( 3, 1, ( 3 ) ),                     \
           ( 4, 1, ( 4 ) ) )                    \
      )                                         \
    /**/
# define GEOTOOL_MARKER_LINE_QUADRANGLE         \
    ( 4, ( ( 1, 1, ( 1 ) ),                     \
           ( 2, 1, ( 2 ) ),                     \
           ( 3, 1, ( 3 ) ),                     \
           ( 4, 1, ( 4 ) ) )                    \
      )                                         \
    /**/
# define GEOTOOL_MARKER_SURFACE_QUADRANGLE      \
    ( 1, ( ( 1, 1, ( 1 ) ) )                    \
      )                                         \
    /**/

/*_________________________________________________*/

# define GEOTOOL_MARKER_POINT_PENTAGON          \
    ( 5, ( ( 1, 1, ( 1 ) ),                     \
           ( 2, 1, ( 2 ) ),                     \
           ( 3, 1, ( 3 ) ),                     \
           ( 4, 1, ( 4 ) ),                     \
           ( 5, 1, ( 5 ) ) )                    \
      )                                         \
    /**/
# define GEOTOOL_MARKER_LINE_PENTAGON           \
    ( 5, ( ( 1, 1, ( 1 ) ),                     \
           ( 2, 1, ( 2 ) ),                     \
           ( 3, 1, ( 3 ) ),                     \
           ( 4, 1, ( 4 ) ),                     \
           ( 5, 1, ( 5 ) ) )                    \
      )                                         \
    /**/
# define GEOTOOL_MARKER_SURFACE_PENTAGON        \
    ( 1, ( ( 1, 1, ( 1 ) ) )                    \
      )                                         \
    /**/

/*_________________________________________________*/


# define GEOTOOL_MARKER_POINT_HEXAGON           \
    ( 6, ( ( 1, 1, ( 1 ) ),                     \
           ( 2, 1, ( 2 ) ),                     \
           ( 3, 1, ( 3 ) ),                     \
           ( 4, 1, ( 4 ) ),                     \
           ( 5, 1, ( 5 ) ),                     \
           ( 6, 1, ( 6 ) ) )                    \
      )                                         \
    /**/
# define GEOTOOL_MARKER_LINE_HEXAGON            \
    ( 6, ( ( 1, 1, ( 1 ) ),                     \
           ( 2, 1, ( 2 ) ),                     \
           ( 3, 1, ( 3 ) ),                     \
           ( 4, 1, ( 4 ) ),                     \
           ( 5, 1, ( 5 ) ),                     \
           ( 6, 1, ( 6 ) ) )                    \
      )                                         \
    /**/
# define GEOTOOL_MARKER_SURFACE_HEXAGON         \
    ( 1, ( ( 1, 1, ( 1 ) ) )                    \
      )                                         \
    /**/

/*_________________________________________________*/

# define GEOTOOL_MARKER_POINT_CIRCLE            \
    ( 2, ( ( 1, 2, ( 1,3 ) ),                   \
           ( 2, 1, ( 2 ) ) )                    \
      )                                         \
/**/
# define GEOTOOL_MARKER_LINE_CIRCLE             \
    ( 1, ( ( 1, 2, ( 1,2 ) ) )                  \
      )                                         \
    /**/

# define GEOTOOL_MARKER_SURFACE_CIRCLE          \
    ( 1, ( ( 1, 1, ( 1 ) ) )                    \
      )                                         \
    /**/

/*_________________________________________________*/

# define GEOTOOL_MARKER_POINT_ELLIPSE           \
    ( 2, ( ( 1, 4, ( 2,3,4,5 ) ),               \
           ( 2, 1, ( 1 ) ) )                    \
      )                                         \
/**/
# define GEOTOOL_MARKER_LINE_ELLIPSE            \
    ( 1, ( ( 1, 4, ( 1,2,3,4 ) ) )              \
      )                                         \
    /**/

# define GEOTOOL_MARKER_SURFACE_ELLIPSE         \
    ( 1, ( ( 1, 1, ( 1 ) ) )                    \
      )                                         \
    /**/

/*_________________________________________________*/
//false : must be fix!
# define GEOTOOL_MARKER_POINT_PIE               \
    ( 2, ( ( 1, 2, ( 1,3 ) ),                   \
           ( 2, 1, ( 2 ) ) )                    \
      )                                         \
    /**/
# define GEOTOOL_MARKER_LINE_PIE                \
    ( 2, ( ( 1, 4, ( 1,2,3,4 ) ),               \
           ( 2, 1, (    5    ) ) )              \
      )                                         \
    /**/

# define GEOTOOL_MARKER_SURFACE_PIE             \
    ( 1, ( ( 1, 1, ( 1 ) ) )                    \
      )                                         \
    /**/

/*_________________________________________________*/

# define GEOTOOL_MARKER_POINT_SPECIAL_1A          \
    ( 1, ( ( 1, 24, ( 1,2,3,4,5,6,                \
                      7,8,9,10,11,12,             \
                      13,14,15,16,17,18,          \
                      19,20,21,22,23,24 ) ) ) )   \
    /**/
# define GEOTOOL_MARKER_LINE_SPECIAL_1A            \
    ( 4, ( ( 1, 2, ( 1,5 ) ),                      \
           ( 2, 2, ( 2,6 ) ),                      \
           ( 3, 2, ( 3,7 ) ),                      \
           ( 4, 2, ( 4,8 ) ) )                     \
      )                                            \
    /**/
# define GEOTOOL_MARKER_SURFACE_SPECIAL_1A      \
    ( 1, ( ( 1, 2, ( 1,2 ) ) )                  \
      )                                         \
    /**/

/*_________________________________________________*/

# define GEOTOOL_MARKER_POINT_SPECIAL_1B          \
    ( 1, ( ( 1, 12, ( 1,2,3,4,5,6,                \
                      7,8,9,10,11,12 ) ) ) )      \
    /**/
# define GEOTOOL_MARKER_LINE_SPECIAL_1B            \
    ( 3, ( ( 1, 2, ( 1,2 ) ),                      \
           ( 2, 1, ( 3   ) ),                      \
           ( 3, 1, ( 4   ) )                       \
           )                                       \
      )                                            \
    /**/
# define GEOTOOL_MARKER_SURFACE_SPECIAL_1B      \
    ( 1, ( ( 1, 1, ( 1 ) ) )                    \
      )                                         \
    /**/
/*_________________________________________________*/

# define GEOTOOL_MARKER_POINT_PEANUT              \
    ( 1, ( ( 1, 8, ( 1,2,3,4,5,6,7,8 ) ) ) )      \
    /**/
# define GEOTOOL_MARKER_LINE_PEANUT                \
    ( 1, ( ( 1, 1, ( 1 ) ) ) )                     \
    /**/
# define GEOTOOL_MARKER_SURFACE_PEANUT          \
    ( 1, ( ( 1, 1, ( 1 ) ) ) )                  \
    /**/
/*_________________________________________________*/

# define GEOTOOL_MARKER_POINT_CYLINDRE            \
    ( 1, ( ( 1, 10, ( 1,2,3,4,5,6,                \
                      7,8,9,10 ) ) ) )            \
    /**/
# define GEOTOOL_MARKER_LINE_CYLINDRE           \
    ( 12, ( (  1, 1, (  1 ) ),                  \
            (  2, 1, (  2 ) ),                  \
            (  3, 1, (  3 ) ),                  \
            (  4, 1, (  4 ) ),                  \
            (  5, 1, (  5 ) ),                  \
            (  6, 1, (  6 ) ),                  \
            (  7, 1, (  7 ) ),                  \
            (  8, 1, (  8 ) ),                  \
            (  9, 1, (  9 ) ),                  \
            ( 10, 1, ( 10 ) ),                  \
            ( 11, 1, ( 11 ) ),                  \
            ( 12, 1, ( 12 ) )                   \
            )                                   \
      )                                         \
    /**/
# define GEOTOOL_MARKER_SURFACE_CYLINDRE        \
    ( 3, ( ( 1, 1, ( 1 ) ),                     \
           ( 2, 1, ( 2 ) ),                     \
           ( 3, 4, ( 3,4,5,6 ) )                \
           )                                    \
      )                                         \
    /**/
# define GEOTOOL_MARKER_VOLUME_CYLINDRE         \
    ( 1, ( ( 1, 1, ( 1 ) ) )                    \
      )                                         \
    /**/
/*_________________________________________________*/

# define GEOTOOL_MARKER_POINT_TUBE                \
    ( 1, ( ( 1, 18, ( 1,2,3,4,5,6,                \
                      7,8,9,10,11,12,             \
                      13,14,15,16,17,18 ) ) ) )   \
    /**/
# define GEOTOOL_MARKER_LINE_TUBE               \
    ( 12, ( (  1, 1, (  1 ) ),                  \
            (  2, 1, (  2 ) ),                  \
            (  3, 1, (  3 ) ),                  \
            (  4, 1, (  4 ) ),                  \
            (  5, 1, (  5 ) ),                  \
            (  6, 1, (  6 ) ),                  \
            (  7, 1, (  7 ) ),                  \
            (  8, 1, (  8 ) ),                  \
            (  9, 1, (  9 ) ),                  \
            ( 10, 1, ( 10 ) ),                  \
            ( 11, 1, ( 11 ) ),                  \
            ( 12, 1, ( 12 ) )                   \
            )                                   \
      )                                         \
    /**/
# define GEOTOOL_MARKER_SURFACE_TUBE            \
    ( 5, ( ( 1, 4, ( 1,2,3,4 ) ),                \
           ( 2, 4, ( 5,6,7,8 ) ),                \
           ( 3, 4, ( 9,10,11,12 ) ),\
           ( 4, 4, ( 13,14,15,16 ) ),           \
           ( 5, 4, ( 17,18,19,20 ) )            \
           )                                    \
      )                                         \
    /**/
# define GEOTOOL_MARKER_VOLUME_TUBE         \
    ( 1, ( ( 1, 4, ( 1,2,3,4 ) ) )              \
      )                                         \
    /**/
/*_________________________________________________*/

# define GEOTOOL_MARKER_POINT_SPHERE              \
    ( 1, ( ( 1, 7, ( 1,2,3,4,5,6,7 ) ) ) )        \
    /**/
# define GEOTOOL_MARKER_LINE_SPHERE             \
    ( 12, ( (  1, 1, (  1 ) ),                  \
            (  2, 1, (  2 ) ),                  \
            (  3, 1, (  3 ) ),                  \
            (  4, 1, (  4 ) ),                  \
            (  5, 1, (  5 ) ),                  \
            (  6, 1, (  6 ) ),                  \
            (  7, 1, (  7 ) ),                  \
            (  8, 1, (  8 ) ),                  \
            (  9, 1, (  9 ) ),                  \
            ( 10, 1, ( 10 ) ),                  \
            ( 11, 1, ( 11 ) ),                  \
            ( 12, 1, ( 12 ) )                   \
            )                                   \
      )                                         \
    /**/
# define GEOTOOL_MARKER_SURFACE_SPHERE          \
    ( 1, ( ( 1, 8, ( 1,2,3,4,5,6,7,8 ) )        \
           )                                    \
      )                                         \
    /**/
# define GEOTOOL_MARKER_VOLUME_SPHERE           \
    ( 1, ( ( 1, 1, ( 1 ) ) )                    \
      )                                         \
    /**/

/*_________________________________________________*/

# define GEOTOOL_MARKER_POINT_HEXAHEDRON          \
    ( 1, ( ( 1, 8, ( 1,2,3,4,5,6,7,8 ) ) ) )      \
    /**/
# define GEOTOOL_MARKER_LINE_HEXAHEDRON         \
    ( 12, ( (  1, 1, (  1 ) ),                  \
            (  2, 1, (  2 ) ),                  \
            (  3, 1, (  3 ) ),                  \
            (  4, 1, (  4 ) ),                  \
            (  5, 1, (  5 ) ),                  \
            (  6, 1, (  6 ) ),                  \
            (  7, 1, (  7 ) ),                  \
            (  8, 1, (  8 ) ),                  \
            (  9, 1, (  9 ) ),                  \
            ( 10, 1, ( 10 ) ),                  \
            ( 11, 1, ( 11 ) ),                  \
            ( 12, 1, ( 12 ) )                   \
            )                                   \
      )                                         \
    /**/
# define GEOTOOL_MARKER_SURFACE_HEXAHEDRON      \
    ( 6, ( ( 1, 1, ( 1 ) ),                     \
           ( 2, 1, ( 2 ) ),                     \
           ( 3, 1, ( 3 ) ),                     \
           ( 4, 1, ( 4 ) ),                     \
           ( 5, 1, ( 5 ) ),                     \
           ( 6, 1, ( 6 ) )                      \
           )                                    \
      )                                         \
    /**/
# define GEOTOOL_MARKER_VOLUME_HEXAHEDRON       \
    ( 1, ( ( 1, 1, ( 1 ) ) )                    \
      )                                         \
    /**/
/*_________________________________________________*/

# define GEOTOOL_MARKER_POINT_TETRAHEDRON          \
    ( 1, ( ( 1, 4, ( 1,2,3,4 ) ) ) )               \
    /**/
# define GEOTOOL_MARKER_LINE_TETRAHEDRON       \
    ( 6, ( (  1, 1, (  1 ) ),                  \
            (  2, 1, (  2 ) ),                  \
            (  3, 1, (  3 ) ),                  \
            (  4, 1, (  4 ) ),                  \
            (  5, 1, (  5 ) ),                  \
            (  6, 1, (  6 ) )                   \
            )                                   \
      )                                         \
    /**/
# define GEOTOOL_MARKER_SURFACE_TETRAHEDRON     \
    ( 4, ( ( 1, 1, ( 1 ) ),                     \
           ( 2, 1, ( 2 ) ),                     \
           ( 3, 1, ( 3 ) ),                     \
           ( 4, 1, ( 4 ) )                      \
           )                                    \
      )                                         \
    /**/
# define GEOTOOL_MARKER_VOLUME_TETRAHEDRON      \
    ( 1, ( ( 1, 1, ( 1 ) ) )                    \
      )                                         \
    /**/
/*_________________________________________________*/

# define GEOTOOL_MARKER_POINT_CUBE            \
    ( 1, ( ( 1, 8, ( 1,2,3,4,5,6,7,8 ) ) ) )      \
    /**/
# define GEOTOOL_MARKER_LINE_CUBE           \
    ( 12, ( (  1, 1, (  1 ) ),                  \
            (  2, 1, (  2 ) ),                  \
            (  3, 1, (  3 ) ),                  \
            (  4, 1, (  4 ) ),                  \
            (  5, 1, (  5 ) ),                  \
            (  6, 1, (  6 ) ),                  \
            (  7, 1, (  7 ) ),                  \
            (  8, 1, (  8 ) ),                  \
            (  9, 1, (  9 ) ),                  \
            ( 10, 1, ( 10 ) ),                  \
            ( 11, 1, ( 11 ) ),                  \
            ( 12, 1, ( 12 ) )                   \
            )                                   \
      )                                         \
    /**/
# define GEOTOOL_MARKER_SURFACE_CUBE        \
    ( 6, ( ( 1, 1, ( 1 ) ),                     \
           ( 2, 1, ( 2 ) ),                     \
           ( 3, 1, ( 3 ) ),                     \
           ( 4, 1, ( 4 ) ),                     \
           ( 5, 1, ( 5 ) ),                     \
           ( 6, 1, ( 6 ) )                      \
           )                                    \
      )                                         \
    /**/
# define GEOTOOL_MARKER_VOLUME_CUBE         \
    ( 1, ( ( 1, 1, ( 1 ) ) )                    \
      )                                         \
    /**/

/*_________________________________________________*/
//special3D_1

# define GEOTOOL_MARKER_POINT_SPECIAL3D_1                 \
    ( 1, ( ( 1, 12, ( 1,2,3,4,5,6,7,8,9,10,11,12 ) ) ) )  \
    /**/
# define GEOTOOL_MARKER_LINE_SPECIAL3D_1                                \
    ( 1, ( ( 1, 17, ( 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17 ) ) ) ) \
    /**/
# define GEOTOOL_MARKER_SURFACE_SPECIAL3D_1        \
    ( 2, ( ( 1, 5, ( 1,2,3,4,5 ) ),                \
           ( 3, 4, ( 6,7,8,9 ) )                   \
           )                                       \
      )                                            \
    /**/
# define GEOTOOL_MARKER_VOLUME_SPECIAL3D_1         \
    ( 1, ( ( 1, 1, ( 1 ) ) )                       \
      )                                            \
    /**/



/*_________________________________________________*/

#if 0
# define GEOTOOL_MARKER_SURFACE_DEFAULT         \
    ( 1, ( ( 1, 1, ( 1 ) ) )                    \
      )                                         \
    /**/




# define GEOTOOL_MARKER_VOLUME_DEFAULT          \
    ( 1, ( ( 1, 1, ( 1 ) ) )                    \
      )                                         \
    /**/

#endif

// Accessors

# define GEOTOOL_SHAPE_NAME_CLASS(i) BOOST_PP_TUPLE_ELEM(7, 0, BOOST_PP_ARRAY_ELEM(i,GEOTOOL_SHAPE))
# define GEOTOOL_SHAPE_DIM(i) BOOST_PP_TUPLE_ELEM(7, 1, BOOST_PP_ARRAY_ELEM(i,GEOTOOL_SHAPE))
# define GEOTOOL_SHAPE_NBSURFACE(i) BOOST_PP_TUPLE_ELEM(7, 2, BOOST_PP_ARRAY_ELEM(i,GEOTOOL_SHAPE))
# define GEOTOOL_SHAPE_NBVOLUME(i) BOOST_PP_TUPLE_ELEM(7, 3, BOOST_PP_ARRAY_ELEM(i,GEOTOOL_SHAPE))
# define GEOTOOL_SHAPE_NAME_STR(i) BOOST_PP_TUPLE_ELEM(7, 4, BOOST_PP_ARRAY_ELEM(i,GEOTOOL_SHAPE))
# define GEOTOOL_SHAPE_NBPARAM(i) BOOST_PP_TUPLE_ELEM(7, 5, BOOST_PP_ARRAY_ELEM(i,GEOTOOL_SHAPE))
# define GEOTOOL_SHAPE_NAME_MACRO(i) BOOST_PP_TUPLE_ELEM(7, 6, BOOST_PP_ARRAY_ELEM(i,GEOTOOL_SHAPE))

# define GEOTOOL_MARKER_POINT_INDICE(O) BOOST_PP_TUPLE_ELEM(3, 0, O)
# define GEOTOOL_MARKER_POINT_NBMARK(F,i) BOOST_PP_TUPLE_ELEM(3, 1, BOOST_PP_ARRAY_ELEM(i,F))
# define GEOTOOL_MARKER_POINT_ARRAYMARK(O) BOOST_PP_TUPLE_ELEM(3, 2, O)
# define GEOTOOL_MARKER_POINT_MARKVALUE(F,i,j)                           \
    BOOST_PP_TUPLE_ELEM( GEOTOOL_MARKER_POINT_NBMARK(F,i),j,GEOTOOL_MARKER_POINT_ARRAYMARK(BOOST_PP_ARRAY_ELEM(i, F)))

# define GEOTOOL_MARKER_LINE_INDICE(O) BOOST_PP_TUPLE_ELEM(3, 0, O)
# define GEOTOOL_MARKER_LINE_NBMARK(F,i) BOOST_PP_TUPLE_ELEM(3, 1, BOOST_PP_ARRAY_ELEM(i,F))
# define GEOTOOL_MARKER_LINE_ARRAYMARK(O) BOOST_PP_TUPLE_ELEM(3, 2, O)
# define GEOTOOL_MARKER_LINE_MARKVALUE(F,i,j)                           \
    BOOST_PP_TUPLE_ELEM( GEOTOOL_MARKER_LINE_NBMARK(F,i),j,GEOTOOL_MARKER_LINE_ARRAYMARK(BOOST_PP_ARRAY_ELEM(i, F)))

# define GEOTOOL_MARKER_SURFACE_INDICE(O) BOOST_PP_TUPLE_ELEM(3, 0, O)
# define GEOTOOL_MARKER_SURFACE_NBMARK(F,i) BOOST_PP_TUPLE_ELEM(3, 1, BOOST_PP_ARRAY_ELEM(i,F))
# define GEOTOOL_MARKER_SURFACE_ARRAYMARK(O) BOOST_PP_TUPLE_ELEM(3, 2, O)
# define GEOTOOL_MARKER_SURFACE_MARKVALUE(F,i,j)                        \
    BOOST_PP_TUPLE_ELEM( GEOTOOL_MARKER_SURFACE_NBMARK(F,i),j,GEOTOOL_MARKER_SURFACE_ARRAYMARK(BOOST_PP_ARRAY_ELEM(i, F)))

# define GEOTOOL_MARKER_VOLUME_INDICE(O) BOOST_PP_TUPLE_ELEM(3, 0, O)
# define GEOTOOL_MARKER_VOLUME_NBMARK(F,i) BOOST_PP_TUPLE_ELEM(3, 1, BOOST_PP_ARRAY_ELEM(i,F))
# define GEOTOOL_MARKER_VOLUME_ARRAYMARK(O) BOOST_PP_TUPLE_ELEM(3, 2, O)
# define GEOTOOL_MARKER_VOLUME_MARKVALUE(F,i,j)                         \
    BOOST_PP_TUPLE_ELEM( GEOTOOL_MARKER_VOLUME_NBMARK(F,i),j,GEOTOOL_MARKER_VOLUME_ARRAYMARK(BOOST_PP_ARRAY_ELEM(i, F)))

/*_________________________________________________*/
/*_________________________________________________*/
/*_________________________________________________*/


namespace Feel
{

namespace GeoTool
{

typedef node<double>::type node_type;

class GeoGMSHTool;
typedef boost::shared_ptr< GeoGMSHTool> GeoGMSHTool_ptrtype;

typedef std::map<uint16_type,uint16_type> map_data_type;
typedef std::vector<map_data_type> vec_map_data_type;
typedef boost::shared_ptr<vec_map_data_type> vec_map_data_ptrtype;

//if bool=true => surface stoker dans un tableau gmsh
typedef std::vector<std::map<uint16_type,bool> > vec_map_data_surf1_type;
typedef boost::shared_ptr<vec_map_data_surf1_type> vec_map_data_surf1_ptrtype;
//=> la string est le nom de ce tableau
typedef std::vector<std::map<uint16_type,std::string> > vec_map_data_surf2_type;
typedef boost::shared_ptr<vec_map_data_surf2_type> vec_map_data_surf2_ptrtype;
// list of pt define in more in the surface
typedef std::vector<std::map<uint16_type,std::list<uint16_type> > > vec_map_data_ptsinsurf_type;
typedef boost::shared_ptr<vec_map_data_ptsinsurf_type> vec_map_data_ptsinsurf_ptrtype;

typedef std::map<int,std::list<int> > map_surfaceLoop_type;
//typedef boost::shared_ptr<map_surfaceLoop_type> map_surfaceLoop_ptrtype;


typedef boost::tuple< GeoGMSHTool_ptrtype,
        vec_map_data_ptrtype,
        std::string,
        std::string,
        vec_map_data_surf1_ptrtype,
        vec_map_data_surf2_ptrtype,
        vec_map_data_surf1_ptrtype,
        vec_map_data_ptsinsurf_ptrtype,
        map_surfaceLoop_type > data_geo_type;
typedef boost::shared_ptr<data_geo_type> data_geo_ptrtype;


void run( data_geo_ptrtype __dg );



#define GEOTOOL_INSTANTIATES_FOR_COMP(r, state)                         \
        BOOST_PP_NOT_EQUAL( BOOST_PP_TUPLE_ELEM(2, 0, state),           \
                            BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(2, 1, state)) \
                            )                                           \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#define GEOTOOL_INSTANTIATES_FOR_INCR(r, state)             \
        (                                                   \
         BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(2, 0, state)),	\
         BOOST_PP_TUPLE_ELEM(2, 1, state) )                 \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#define GEOTOOL_INSTANTIATES_FOR(r,state)                               \
        void BOOST_PP_CAT(run,GEOTOOL_SHAPE_NAME_CLASS(BOOST_PP_TUPLE_ELEM(2,0,state))) (data_geo_ptrtype dg); \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
BOOST_PP_FOR( ( 0, BOOST_PP_SUB( BOOST_PP_ARRAY_SIZE( GEOTOOL_SHAPE ),1 ) ),
              GEOTOOL_INSTANTIATES_FOR_COMP,
              GEOTOOL_INSTANTIATES_FOR_INCR,
              GEOTOOL_INSTANTIATES_FOR )




} // namespace GeoTool

namespace GeoTool
{

/*_________________________________________________*
 *_________________________________________________*
 *_________________________________________________*
 * GeoGMSHTool :                                   *
 *_________________________________________________*
 *_________________________________________________*
 *_________________________________________________*/


class Node
{
public :

    Node()
        :
        M_node( new node_type() )
    {}

    Node( double __x ) :
        M_node( new node_type( 1 ) )
    {
        ( *M_node )( 0 )=__x;
    }

    Node( double __x, double __y ) :
        M_node( new node_type( 2 ) )
    {
        ( *M_node )( 0 )=__x;
        ( *M_node )( 1 )=__y;
    }

    Node( double __x, double __y, double __z ) :
        M_node( new node_type( 3 ) )
    {
        ( *M_node )( 0 )=__x;
        ( *M_node )( 1 )=__y;
        ( *M_node )( 2 )=__z;
    }

    Node( Node const & m )
        :
        M_node( m.M_node )
    {}

    Node operator=( Node const & m )
    {
        M_node.reset( new node_type( *( m.M_node ) ) );
        return *this;
    }

    double operator()( uint16_type n ) const
    {
        return this->getNode()( n );
    }

    double & operator()( uint16_type n )
    {
        return ( *M_node )( n );
    }

    node_type
    getNode() const
    {
        return *M_node;
    }

    node_type &
    getNode()
    {
        return *M_node;
    }

    boost::shared_ptr<node_type> M_node;
};

/*_________________________________________________*/

class Loop
{
public :

    Loop( Loop const & L ) : M_loop( L.M_loop ) {}

    Loop()
    {
        M_loop.clear();
    }

    void  operator=( Loop m )
    {
        this->M_loop=m.M_loop;
    }
    Loop  operator>>( int __n )
    {
        M_loop.push_back( __n );
        return *this;
    }

    uint16_type size()
    {
        return M_loop.size();
    }

    std::list<int>::const_iterator begin() const
    {
        return M_loop.begin();
    }
    std::list<int>::const_iterator end() const
    {
        return M_loop.end();
    }

    std::list<int> M_loop;
};


/*_________________________________________________*
 *_________________________________________________*
 *_________________________________________________*
 * GeoGMSHTool :                                   *
 *_________________________________________________*
 *_________________________________________________*
 *_________________________________________________*/


class GeoGMSHTool
{
public:

    typedef node<double>::type node_type;

    /*            // list de < nameMesh, meshSize >
    typedef boost::tuple<std::string,double> names_base_type;
    typedef std::list< names_base_type > names_type;
    typedef std::map< std::string, names_type > map_shape_names_type;
    typedef names_type::const_iterator names_const_iterator_type;
    typedef map_shape_names_type::const_iterator map_shape_names_const_iterator_type;
    */
    typedef boost::tuple<std::string,std::string,uint16_type> marker_base_type;
    typedef std::map<std::string,std::list<marker_base_type > > marker_markerName_type;
    typedef std::map< std::string, marker_markerName_type > marker_type_type;
    typedef std::map< std::string, marker_type_type > marker_name_type;
    typedef std::map< std::string, marker_type_type > marker_shape_type;

    typedef marker_markerName_type::const_iterator marker_markerName_const_iterator_type;
    typedef marker_type_type::const_iterator marker_type_const_iterator_type;
    typedef marker_name_type::const_iterator marker_name_const_iterator_type;
    typedef marker_shape_type::const_iterator marker_shape_const_iterator_type;

    typedef std::vector<node_type> parameter_rectangle_type;
    typedef std::map<std::string, parameter_rectangle_type > parameter_name_type;
    typedef std::map<std::string, parameter_name_type > parameter_shape_type;
    typedef parameter_name_type::const_iterator parameter_name_const_iterator_type;
    typedef parameter_shape_type::const_iterator parameter_shape_const_iterator_type;


    // gestion des lignes : shape,name,value,meshSize
    typedef boost::tuple<std::string,std::string,uint16_type,double > ligne_type;
    typedef std::list< ligne_type > ligne_type_type;
    typedef std::list< ligne_type_type > ligne_name_type;
    typedef ligne_type_type::const_iterator ligne_type_const_iterator_type;
    typedef ligne_name_type::const_iterator ligne_name_const_iterator_type;

    // gestion des surfaces : shape,name,(numGlobSurface,valueOfLineloop),meshSize
    typedef boost::tuple<std::string,std::string,std::pair<int,int>,double > surface_type;
    typedef std::list< surface_type > surface_type_type;
    typedef std::list< surface_type_type > surface_name_type;
    typedef surface_type_type::const_iterator surface_type_const_iterator_type;
    typedef surface_name_type::const_iterator surface_name_const_iterator_type;

    // gestion des volumes : shape,name,(numGlobVolume,value),meshSize
    typedef boost::tuple<std::string,std::string, std::pair<int,int>,double > volume_type;
    typedef std::list< volume_type > volume_type_type;
    typedef std::list< volume_type_type > volume_name_type;
    typedef volume_type_type::const_iterator volume_type_const_iterator_type;
    typedef volume_name_type::const_iterator volume_name_const_iterator_type;

    // gestion des surfaceLoop : shape,name, numLoopLoc->list<value>
    typedef boost::tuple<std::string,std::string, std::map< int, std::list<int> > > surfaceloop_type;
    typedef std::list< surfaceloop_type > surfaceloop_type_type;
    typedef std::list< surfaceloop_type_type > surfaceloop_name_type;
    typedef surfaceloop_type_type::const_iterator surfaceloop_type_const_iterator_type;
    typedef surfaceloop_name_type::const_iterator surfaceloop_name_const_iterator_type;


    GeoGMSHTool( uint16_type __dim, std::string __shape="NO_SHAPE", std::string __name="NO_NAME", double __meshSize=0.1 );

    GeoGMSHTool( uint16_type __dim,  std::string const & geoUserStr, double __meshSize=0.1, std::string __shape="NO_SHAPE", std::string __name="NO_NAME" );

    GeoGMSHTool( GeoGMSHTool const & m );

    void zeroCpt();

    void operator=( GeoGMSHTool const & m );

    GeoGMSHTool operator+( const GeoGMSHTool & m );
    GeoGMSHTool operator-( const GeoGMSHTool & m );

    GeoGMSHTool opFusion( const GeoGMSHTool & m,int __typeop );

    void init( int orderGeo,
               std::string gmshFormatVersion,
               double hmin=0,double hmax=1e22,
               int refine=0,
               bool optimize3dNetgen=true,
               GMSH_PARTITIONER partitioner=GMSH_PARTITIONER_CHACO,
               int partitions=1,
               bool partition_file=false );

    /*
     *
     */
    void initData( std::string __shape,
                   std::string __name,
                   double __meshSize,
                   std::vector<GeoTool::Node> & __param,
                   uint16_type dim,
                   uint16_type __nbligne,
                   uint16_type __nbsurface,
                   uint16_type __nbvolume );

    /*
     *Utile pour la fct geoStr()
     *Pas de maj pour cptSurface et cptVolume car traitement different
     */
    void updateData( GeoGMSHTool const & m );

    /*
     * Update the output stringstream wich generate the gmsh code
     */
    void updateOstr( std::string __str )
    {
        *M_ostr << __str;
    }

    /*
     * Generate the gmsh code
     */
    void geoStr();

    /*
     * Clean
     */
    void cleanOstr()
    {
        M_ostr.reset( new std::ostringstream() );
    }





    BOOST_PARAMETER_MEMBER_FUNCTION(
        ( typename Feel::detail::mesh<Args>::ptrtype ), // return type
        createMesh, // function name
        tag,
        ( required
          ( mesh, * )
          ( name, ( std::string ) )
        ) //required
        ( optional
          ( format,         *, option(_name="gmsh.format").template as<int>() )
          ( straighten,     *( boost::is_integral<mpl::_> ), 1 )
          ( refine,          *( boost::is_integral<mpl::_> ), 0 )
          ( partitions,   *( boost::is_integral<mpl::_> ), Environment::worldComm().size() )
          ( partition_file,   *( boost::is_integral<mpl::_> ), 0 )
          ( partitioner,   *( boost::is_integral<mpl::_> ), GMSH_PARTITIONER_CHACO )
          ( worldcomm,      *, Environment::worldComm() )
          ( hmin,     ( double ), 0 )
          ( hmax,     ( double ), 1e22 )
          ( optimize3d_netgen, *( boost::is_integral<mpl::_> ), true )
        ) //optional
    )
    {
        typedef typename Feel::detail::mesh<Args>::type _mesh_type;
        typedef typename Feel::detail::mesh<Args>::ptrtype _mesh_ptrtype;

        _mesh_ptrtype _mesh( mesh );
        _mesh->setWorldComm( worldcomm );

        if ( worldcomm.isActive() )
        {

            this->cleanOstr();
            this->zeroCpt();
            Gmsh gmsh( _mesh_type::nDim, _mesh_type::nOrder, worldcomm );
            gmsh.setRecombine( _mesh_type::shape_type::is_hypercube );
            gmsh.setRefinementLevels( refine );
            gmsh.setFileFormat( (GMSH_FORMAT)format );
            gmsh.setNumberOfPartitions( partitions );
            gmsh.setPartitioner( partitioner );
            gmsh.setMshFileByPartition( partition_file );
            this->init( _mesh_type::nOrder,gmsh.version(),
                        hmin,hmax,refine,
                        optimize3d_netgen,
                        partitioner,partitions,partition_file );

            std::string geostring;

            if ( M_geoIsDefineByUser )
            {
                geostring= M_ostrDefineByUser->str();
            }

            else
            {
                this->geoStr();
                geostring = M_ostr->str();
            }


            std::string fname;
            bool gen;
            boost::tie( fname, gen ) = gmsh.generate( name,
                                                      geostring,
                                                      false,false,false );

            ImporterGmsh<_mesh_type> import( fname, FEELPP_GMSH_FORMAT_VERSION, worldcomm );
            _mesh->accept( import );
            _mesh->components().set ( MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
            _mesh->updateForUse();

            if ( straighten && _mesh_type::nOrder > 1 )
                return straightenMesh( _mesh=_mesh,
                                       _worldcomm=worldcomm.subWorldComm() );

        } // if (worldcomm.isActive())

        return _mesh;
    }


#if 0 // old function
    template<typename mesh_type>
    boost::shared_ptr<mesh_type>
    createMesh( std::string name, int straighten = 1, WorldComm const& worldcomm=Environment::worldComm() )
    {
        boost::shared_ptr<mesh_type> mesh( new mesh_type );
        mesh->setWorldComm( worldcomm );

        if ( worldcomm.isActive() )
        {
            this->cleanOstr();
            this->zeroCpt();

            Gmsh gmsh( mesh_type::nDim,mesh_type::nOrder, worldcomm );
            gmsh.setOrder( mesh_type::nOrder );
            gmsh.setRecombine( mesh_type::shape_type::is_hypercube );

            this->init( mesh_type::nOrder,gmsh.version() );

            std::string geostring;

            if ( M_geoIsDefineByUser )
            {
                geostring= M_ostrDefineByUser->str();
            }

            else
            {
                this->geoStr();
                geostring = M_ostr->str();
            }

            std::string fname = gmsh.generate( name,
                                               geostring,false,false,false );

            ImporterGmsh<mesh_type> import( fname, FEELPP_GMSH_FORMAT_VERSION, worldcomm );
            mesh->accept( import );
            mesh->components().set ( MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
            mesh->updateForUse();

            if ( straighten && mesh_type::nOrder > 1 )
                return straightenMesh( mesh );
        } // if (worldcomm.isActive())

        return mesh;
    }
#endif
    /*_________________________________________________*
     *_________________________________________________*
     * Accessor                                        *
     *_________________________________________________*
     *_________________________________________________*/

    uint16_type dim() const
    {
        return  M_dim;
    }
    uint16_type cptPt() const
    {

        return M_cptPt;
    }
    uint16_type cptLine() const
    {
        return M_cptLine;
    }
    uint16_type cptLineLoop() const
    {
        return M_cptLineLoop;
    }
    uint16_type cptSurface() const
    {
        return M_cptSurface;
    }
    uint16_type cptTableau() const
    {
        return M_cptTableau;   //voir les extrudes par exemple
    }
    uint16_type cptSurfaceLoop() const
    {
        return M_cptSurfaceLoop;
    }
    uint16_type cptVolume() const
    {
        return M_cptVolume;
    }

    /*_________________________________________________*
     * Parameter
     *_________________________________________________*/

    parameter_shape_const_iterator_type
    paramShapeBegin() const
    {
        return M_paramShape->begin();
    }

    parameter_shape_const_iterator_type paramShapeEnd() const
    {
        return M_paramShape->end();
    }

    parameter_name_const_iterator_type
    paramNameBegin( std::string __shape ) const
    {
        return M_paramShape->find( __shape )->second.begin();
    }

    parameter_name_const_iterator_type
    paramNameEnd( std::string __shape ) const
    {
        return M_paramShape->find( __shape )->second.end();
    }

    parameter_rectangle_type
    getParameter( std::string __shape, std::string __name ) const
    {
        return M_paramShape->find( __shape )->second.find( __name )->second;
    }

    /*_________________________________________________*
     * Marker
     *_________________________________________________*/
    /*
    marker_shape_const_iterator_type
    markShapeBegin() const
    {
        return M_markShape->begin();
    }

    marker_shape_const_iterator_type
    markShapeEnd() const
    {
        return M_markShape->end();
    }*/


    marker_type_const_iterator_type
    markerTypeBegin( /*std::string __shape*/ ) const
    {
        //return M_markShape->find(__shape)->second.begin();
        return M_markShape->begin();
    }

    marker_type_const_iterator_type
    markerTypeEnd( /*std::string __shape*/ ) const
    {
        //return M_markShape->find(__shape)->second.end();
        return M_markShape->end();
    }
    /*
    marker_type_type
    markerType(std::string __shape) const
    {
        //return M_markShape->find(__shape)->second;
        }*/

    marker_markerName_const_iterator_type
    markerMarkerNameBegin( /*std::string __shape,*/ std::string __type ) const
    {
        //return M_markShape->find(__shape)->second.find(__type)->second.begin();
        return M_markShape->find( __type )->second.begin();
    }

    marker_markerName_const_iterator_type
    markerMarkerNameEnd( /*std::string __shape,*/ std::string __type ) const
    {
        //return M_markShape->find(__shape)->second.find(__type)->second.end();
        return M_markShape->find( __type )->second.end();
    }

    marker_markerName_type
    markerMarkerName( /*std::string __shape,*/ std::string __type ) const
    {
        //return M_markShape->find(__shape)->second.find(__type)->second;
        return M_markShape->find( __type )->second;
    }

    std::list<marker_base_type>::const_iterator
    markerListIndiceBegin( /*std::string __shape,*/ std::string __type ,std::string __markerName ) const
    {
        return M_markShape->find( __type )->second.find( __markerName )->second.begin();
    }

    std::list<marker_base_type>::const_iterator
    markerListIndiceEnd( /*std::string __shape,*/ std::string __type ,std::string __markerName ) const
    {
        //return M_markShape->find(__shape)->second.find(__type)->second.find(__markerName)->second.end();
        return M_markShape->find( __type )->second.find( __markerName )->second.end();
    }


    std::list<marker_base_type>
    getMarkerName( /*std::string __shape,*/ std::string __type ,std::string __markerName ) const
    {
        //return M_markShape->find(__shape)->second.find(__type)->second.find(__markerName)->second;
        return M_markShape->find( __type )->second.find( __markerName )->second;
    }

    /*_________________________________________________*
     *_________________________________________________*
     * Members                                         *
     *_________________________________________________*
     *_________________________________________________*/

    uint16_type M_dim;
    // memory
    uint16_type M_cptPt;
    uint16_type M_cptLine;
    uint16_type M_cptLineLoop;
    uint16_type M_cptSurface;
    uint16_type M_cptTableau;
    uint16_type M_cptSurfaceLoop;
    uint16_type M_cptVolume;

    // gestion des surface : shape,name,value
    // value is the marker associated to the planeSurface (init to 0 and to use when call geoStr())
    //std::list< std::list< boost::tuple<std::string,std::string, uint16_type > > > M_surfaceList;
    boost::shared_ptr<ligne_name_type> M_ligneList;
    boost::shared_ptr<surface_name_type> M_surfaceList;
    boost::shared_ptr<volume_name_type> M_volumeList;
    boost::shared_ptr<surfaceloop_name_type> M_surfaceLoopList;

    boost::shared_ptr<std::ostringstream> M_ostrExtrude;
    boost::shared_ptr<std::ostringstream> M_ostrSurfaceLoop;


    // data containers
    //boost::shared_ptr<map_shape_names_type> M_map_Shape;
    boost::shared_ptr<parameter_shape_type> M_paramShape;
    //boost::shared_ptr<marker_shape_type> M_markShape;
    boost::shared_ptr<marker_type_type> M_markShape;

    // output string
    boost::shared_ptr<std::ostringstream> M_ostr;

    boost::shared_ptr<std::ostringstream> M_ostrDefineByUser;
    bool M_geoIsDefineByUser;
};

/*_________________________________________________*
 *_________________________________________________*
 *_________________________________________________*
 * Function on the namespace                       *
 *_________________________________________________*
 *_________________________________________________*
 *_________________________________________________*/


void run( data_geo_ptrtype __dg );

template <uint16_type Numero>
node_type
param( data_geo_ptrtype __dg );


void
writePoint( uint16_type __numLoc, data_geo_ptrtype __dg ,double __x1,double __x2=0, double __x3=0 );

void
writeLine( uint16_type __numLoc, data_geo_ptrtype __dg ,uint16_type __n1, uint16_type __n2 );

void
writeCircle( uint16_type __numLoc, data_geo_ptrtype __dg ,uint16_type __n1, uint16_type __n2, uint16_type __n3 );

void
writeEllipse( uint16_type __numLoc, data_geo_ptrtype __dg ,uint16_type __n1, uint16_type __n2, uint16_type __n3, uint16_type __n4 );

void
writeSpline( uint16_type __numLoc, data_geo_ptrtype __dg ,Loop __loop );

void
writeBSpline( uint16_type __numLoc, data_geo_ptrtype __dg ,Loop __loop );

void
writeLineLoop( uint16_type __numLoc, data_geo_ptrtype __dg , Loop /*const*/ __loop );

void
writePlaneSurface( uint16_type __numLoc, data_geo_ptrtype __dg , uint16_type __ind );

void
writeRuledSurface( uint16_type __numLoc, data_geo_ptrtype __dg , uint16_type __ind );

void
writeExtrudeSurface( uint16_type __numLoc,data_geo_ptrtype __dg , uint16_type __ind,Loop /*const*/ __loop );

void
writePtInSurface( data_geo_ptrtype __dg , uint16_type __indPt,uint16_type __indSurf );

void
writeSurfaceLoop( uint16_type __numLoc, data_geo_ptrtype __dg , Loop /*const*/ __loop );

void
writeVolume( uint16_type __numLoc, data_geo_ptrtype __dg , uint16_type __ind );

boost::tuple<Node,Node,Node>
computeBasisOrthogonal( node_type dir,node_type centre );


/*_________________________________________________*
 *_________________________________________________*
 *_________________________________________________*
 * PREPROCESSOR METHODS                            *
 *_________________________________________________*
 *_________________________________________________*
 *_________________________________________________*/




#define GEOTOOL_FOR_COMP2(r, state)                                     \
        BOOST_PP_NOT_EQUAL( BOOST_PP_TUPLE_ELEM(4, 0, state),           \
                            BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(4, 1, state)) \
                            )                                           \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#define GEOTOOL_FOR_INCR2(r, state)                         \
        (                                                   \
         BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(4, 0, state)),    \
         BOOST_PP_TUPLE_ELEM(4, 1, state),                  \
         BOOST_PP_TUPLE_ELEM(4, 2, state),                  \
         BOOST_PP_TUPLE_ELEM(4, 3, state) )                 \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#define GEOTOOL_FOR_MARKER_POINT_MACRO2(r, state)                       \
        __listMarker.push_back(boost::make_tuple(this->shape(),this->name(), \
                                                 GEOTOOL_MARKER_POINT_MARKVALUE( BOOST_PP_CAT(GEOTOOL_MARKER_POINT_, \
                                                                                             GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(4,3,state) )), \
                                                                                BOOST_PP_TUPLE_ELEM(4, 2, state), \
                                                                                BOOST_PP_TUPLE_ELEM(4, 0, state) ) \
                                                 )                      \
                               );                                       \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#define GEOTOOL_FOR_MARKER_LINE_MACRO2(r, state)                        \
        __listMarker.push_back(boost::make_tuple(this->shape(),this->name(), \
                                                 GEOTOOL_MARKER_LINE_MARKVALUE( BOOST_PP_CAT(GEOTOOL_MARKER_LINE_, \
                                                                                             GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(4,3,state) )), \
                                                                                BOOST_PP_TUPLE_ELEM(4, 2, state), \
                                                                                BOOST_PP_TUPLE_ELEM(4, 0, state) ) \
                                                 )                      \
                               );                                       \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#if 1
#define GEOTOOL_FOR_MARKER_SURFACE_MACRO2(r, state)                     \
        __listMarker.push_back(boost::make_tuple(this->shape(),this->name(), \
                                                 GEOTOOL_MARKER_SURFACE_MARKVALUE( BOOST_PP_CAT(GEOTOOL_MARKER_SURFACE_, \
                                                                                                GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(4,3,state) )), \
                                                                                   BOOST_PP_TUPLE_ELEM(4, 2, state), \
                                                                                   BOOST_PP_TUPLE_ELEM(4, 0, state) ) \
                                                 )                      \
                               );                                       \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#else
#define GEOTOOL_FOR_MARKER_SURFACE_MACRO2(r, state)                     \
        __listMarker.push_back(boost::make_tuple(this->shape(),this->name(), \
                                                 GEOTOOL_MARKER_SURFACE_MARKVALUE( BOOST_PP_CAT(GEOTOOL_MARKER_SURFACE_, \
                                                                                                BOOST_PP_IF(BOOST_PP_GREATER(GEOTOOL_SHAPE_NBSURFACE(BOOST_PP_TUPLE_ELEM(4,3,state)) ,0), \
                                                                                                GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(4,3,state) ), \
                                                                                                DEFAULT)), \
                                                                                   BOOST_PP_TUPLE_ELEM(4, 2, state), \
                                                                                   BOOST_PP_TUPLE_ELEM(4, 0, state) ) \
                                                 )                      \
                               );                                       \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/

#endif


#if 0
#define GEOTOOL_FOR_MARKER_VOLUME_MACRO2(r, state)                      \
        __listMarker.push_back(boost::make_tuple(this->shape(),this->name(), \
                                                 GEOTOOL_MARKER_VOLUME_MARKVALUE( BOOST_PP_CAT(GEOTOOL_MARKER_VOLUME_, \
                                                                                               BOOST_PP_IF(BOOST_PP_GREATER(GEOTOOL_SHAPE_NBVOLUME(BOOST_PP_TUPLE_ELEM(4,3,state)) ,0), \
                                                                                                           GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(4,3,state)), \
                                                                                                           DEFAULT)), \
                                                                                  BOOST_PP_TUPLE_ELEM(4, 2, state), \
                                                                                  BOOST_PP_TUPLE_ELEM(4, 0, state) ) \
                                                 )                      \
                               );                                       \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#else
#define GEOTOOL_FOR_MARKER_VOLUME_MACRO2(r, state)                      \
        __listMarker.push_back(boost::make_tuple(this->shape(),this->name(), \
                                                 GEOTOOL_MARKER_VOLUME_MARKVALUE( BOOST_PP_CAT(GEOTOOL_MARKER_VOLUME_, \
                                                                                               GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(4,3,state) )), \
                                                                                  BOOST_PP_TUPLE_ELEM(4, 2, state), \
                                                                                  BOOST_PP_TUPLE_ELEM(4, 0, state) ) \
                                                 )                      \
                               );                                       \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#endif


#define GEOTOOL_FOR_MARKER_POINT_MACRO(r, state)                         \
        if (BOOST_PP_CAT(marker,                                        \
                         BOOST_PP_ADD(BOOST_PP_TUPLE_ELEM(3, 0, state),	\
                                      1 ) ) )                           \
            {                                                           \
                BOOST_PP_FOR( (0,                                       \
                               BOOST_PP_SUB(GEOTOOL_MARKER_POINT_NBMARK(BOOST_PP_CAT(GEOTOOL_MARKER_POINT_, \
                                                                                    GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(3,2,state))), \
                                                                       BOOST_PP_TUPLE_ELEM(3, 0, state) ),1), \
                               BOOST_PP_TUPLE_ELEM(3, 0, state),		\
                               BOOST_PP_TUPLE_ELEM(3, 2, state)			\
                               ),                                       \
                              GEOTOOL_FOR_COMP2, GEOTOOL_FOR_INCR2, GEOTOOL_FOR_MARKER_POINT_MACRO2) \
                    }                                                   \

/**/
/*_________________________________________________*/
/*                                                 */
/**/
#define GEOTOOL_FOR_MARKER_LINE_MACRO(r, state)                         \
        if (BOOST_PP_CAT(marker,                                        \
                         BOOST_PP_ADD(BOOST_PP_TUPLE_ELEM(3, 0, state),	\
                                      1 ) ) )                           \
            {                                                           \
                BOOST_PP_FOR( (0,                                       \
                               BOOST_PP_SUB(GEOTOOL_MARKER_LINE_NBMARK(BOOST_PP_CAT(GEOTOOL_MARKER_LINE_, \
                                                                                    GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(3,2,state))), \
                                                                       BOOST_PP_TUPLE_ELEM(3, 0, state) ),1), \
                               BOOST_PP_TUPLE_ELEM(3, 0, state),		\
                               BOOST_PP_TUPLE_ELEM(3, 2, state)			\
                               ),                                       \
                              GEOTOOL_FOR_COMP2, GEOTOOL_FOR_INCR2, GEOTOOL_FOR_MARKER_LINE_MACRO2) \
                    }                                                   \

/**/
/*_________________________________________________*/
/*                                                 */
/**/
#if 1
#define GEOTOOL_FOR_MARKER_SURFACE_MACRO(r, state)                      \
        if (BOOST_PP_CAT(marker,                                        \
                         BOOST_PP_ADD(BOOST_PP_TUPLE_ELEM(3, 0, state),	\
                                      1 ) ) )                           \
            {                                                           \
                BOOST_PP_FOR( (0,                                       \
                               BOOST_PP_SUB(GEOTOOL_MARKER_SURFACE_NBMARK(BOOST_PP_CAT(GEOTOOL_MARKER_SURFACE_, \
                                                                                       GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(3,2,state)) ), \
                                                                          BOOST_PP_TUPLE_ELEM(3, 0, state) ),1), \
                               BOOST_PP_TUPLE_ELEM(3, 0, state),		\
                               BOOST_PP_TUPLE_ELEM(3, 2, state)			\
                               ),                                       \
                              GEOTOOL_FOR_COMP2, GEOTOOL_FOR_INCR2, GEOTOOL_FOR_MARKER_SURFACE_MACRO2) \
                    }                                                   \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#else
#define GEOTOOL_FOR_MARKER_SURFACE_MACRO(r, state)                      \
        if (BOOST_PP_CAT(marker,                                        \
                         BOOST_PP_ADD(BOOST_PP_TUPLE_ELEM(3, 0, state),	\
                                      1 ) ) )                           \
            {                                                           \
                BOOST_PP_FOR( (0,                                       \
                               BOOST_PP_SUB(GEOTOOL_MARKER_SURFACE_NBMARK(BOOST_PP_CAT(GEOTOOL_MARKER_SURFACE_, \
                                                                                       BOOST_PP_IF(BOOST_PP_GREATER(GEOTOOL_SHAPE_NBSURFACE(BOOST_PP_TUPLE_ELEM(3,2,state)),0), \
                                                                                                   GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(3,2,state)), \
                                                                                                   DEFAULT )), \
                                                                          BOOST_PP_TUPLE_ELEM(3, 0, state) ),1) , \
                               BOOST_PP_TUPLE_ELEM(3, 0, state),		\
                               BOOST_PP_TUPLE_ELEM(3, 2, state)			\
                               ),                                       \
                              GEOTOOL_FOR_COMP2, GEOTOOL_FOR_INCR2, GEOTOOL_FOR_MARKER_SURFACE_MACRO2) \
                   }                                                  \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/

#endif

#if 0
#define GEOTOOL_FOR_MARKER_VOLUME_MACRO(r, state)                       \
        if (BOOST_PP_CAT(marker,                                        \
                         BOOST_PP_ADD(BOOST_PP_TUPLE_ELEM(3, 0, state),	\
                                      1 ) ) )                           \
            {                                                           \
                BOOST_PP_FOR( (0,                                       \
                               BOOST_PP_SUB(GEOTOOL_MARKER_VOLUME_NBMARK(BOOST_PP_CAT(GEOTOOL_MARKER_VOLUME_, \
                                                                                      BOOST_PP_IF(BOOST_PP_GREATER(GEOTOOL_SHAPE_NBVOLUME(BOOST_PP_TUPLE_ELEM(3,2,state)) ,0), \
                                                                                                  GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(3,2,state)),\
                                                                                                  DEFAULT )), \
                                                                          BOOST_PP_TUPLE_ELEM(3, 0, state) ),1), \
                               BOOST_PP_TUPLE_ELEM(3, 0, state),		\
                               BOOST_PP_TUPLE_ELEM(3, 2, state)			\
                               ),                                       \
                              GEOTOOL_FOR_COMP2, GEOTOOL_FOR_INCR2, GEOTOOL_FOR_MARKER_VOLUME_MACRO2) \
                    }                                                   \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#else
#define GEOTOOL_FOR_MARKER_VOLUME_MACRO(r, state)                       \
        if (BOOST_PP_CAT(marker,                                        \
                         BOOST_PP_ADD(BOOST_PP_TUPLE_ELEM(3, 0, state),	\
                                      1 ) ) )                           \
            {                                                           \
                BOOST_PP_FOR( (0,                                       \
                               BOOST_PP_SUB(GEOTOOL_MARKER_VOLUME_NBMARK(BOOST_PP_CAT(GEOTOOL_MARKER_VOLUME_, \
                                                                                      GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(3,2,state)) ), \
                                                                         BOOST_PP_TUPLE_ELEM(3, 0, state) ),1), \
                               BOOST_PP_TUPLE_ELEM(3, 0, state),		\
                               BOOST_PP_TUPLE_ELEM(3, 2, state)			\
                               ),                                       \
                              GEOTOOL_FOR_COMP2, GEOTOOL_FOR_INCR2, GEOTOOL_FOR_MARKER_VOLUME_MACRO2) \
                    }                                                   \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/

#endif


#define GEOTOOL_FOR_COMP1(r, state)                                     \
        BOOST_PP_NOT_EQUAL( BOOST_PP_TUPLE_ELEM(3, 0, state),           \
                            BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(3, 1, state)) \
                            )                                           \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#define GEOTOOL_FOR_INCR1(r, state)                         \
        (                                                   \
         BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(3, 0, state)),	\
         BOOST_PP_TUPLE_ELEM(3, 1, state),                  \
         BOOST_PP_TUPLE_ELEM(3, 2, state) )                 \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#define GEOTOOL_FOR_COMP(r, state)                                      \
        BOOST_PP_LESS( BOOST_PP_TUPLE_ELEM(2, 0, state),           \
                            BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(2, 1, state)) \
                            )                                           \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#define GEOTOOL_FOR_INCR(r, state)                          \
        (                                                   \
         BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(2, 0, state)),	\
         BOOST_PP_TUPLE_ELEM(2, 1, state) )                 \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#define GEOTOOL_SHAPE_PARAM(r, state)                                   \
        M_param[BOOST_PP_TUPLE_ELEM(2,0,state)] = BOOST_PP_CAT( __param, \
                                                                 BOOST_PP_TUPLE_ELEM(2,0,state) ); \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#define GEOTOOL_SHAPE_FOR_PARAM_SIGNATURE(r, state)                     \
    Node BOOST_PP_CAT( __param, BOOST_PP_TUPLE_ELEM(2,0,state) ) = Node(0,0,0) BOOST_PP_COMMA() \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#define GEOTOOL_SHAPE_PARAM_SIGNATURE(state)                           \
        BOOST_PP_FOR( (0, BOOST_PP_SUB(GEOTOOL_SHAPE_NBPARAM(BOOST_PP_TUPLE_ELEM(2,0,state)),1) ), \
                      GEOTOOL_FOR_COMP,                                 \
                      GEOTOOL_FOR_INCR,                                 \
                      GEOTOOL_SHAPE_FOR_PARAM_SIGNATURE)                \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#define GEOTOOL_SHAPE_CLASS(r,state)                                    \
        class GEOTOOL_SHAPE_NAME_CLASS(BOOST_PP_TUPLE_ELEM(2,0,state)) : public GeoGMSHTool \
        {                                                               \
        public :                                                        \
                                                                        \
            typedef GeoGMSHTool::node_type node_type;                   \
            typedef GeoTool::Node Node;                                 \
                                                                        \
                                                                        \
            GEOTOOL_SHAPE_NAME_CLASS(BOOST_PP_TUPLE_ELEM(2,0,state))(double __meshSize, \
                                                                     std::string __name, \
                                                                     GEOTOOL_SHAPE_PARAM_SIGNATURE(state) \
                                                                     uint16_type type = 0 ) /*Ne sert a rien, juste a cause de la virgule au dessus)*/ \
                :                                                       \
                GeoGMSHTool( GEOTOOL_SHAPE_DIM(BOOST_PP_TUPLE_ELEM(2,0,state)),shape(), __name, __meshSize), \
                M_name(__name)                                         \
                {                                                       \
                    M_param.resize( GEOTOOL_SHAPE_NBPARAM(BOOST_PP_TUPLE_ELEM(2,0,state))); \
                    BOOST_PP_FOR( (0, BOOST_PP_SUB(GEOTOOL_SHAPE_NBPARAM(BOOST_PP_TUPLE_ELEM(2,0,state)),1) ), \
                                  GEOTOOL_FOR_COMP,                     \
                                  GEOTOOL_FOR_INCR,                     \
                                  GEOTOOL_SHAPE_PARAM);                 \
                                                                        \
                    initData(shape(),                                   \
                             __name,                                    \
                             __meshSize,                                \
                             M_param,                                  \
                             GEOTOOL_SHAPE_DIM(BOOST_PP_TUPLE_ELEM(2,0,state)), \
                             1,                                         \
                             GEOTOOL_SHAPE_NBSURFACE(BOOST_PP_TUPLE_ELEM(2,0,state)), \
                             GEOTOOL_SHAPE_NBVOLUME(BOOST_PP_TUPLE_ELEM(2,0,state))); \
                }                                                       \
                                                                        \
                                                                        \
                                                                        \
                                                                        \
            GEOTOOL_SHAPE_NAME_CLASS(BOOST_PP_TUPLE_ELEM(2,0,state))(const GEOTOOL_SHAPE_NAME_CLASS(BOOST_PP_TUPLE_ELEM(2,0,state)) & m) \
                :                                                       \
                GeoGMSHTool(m),                                         \
                M_param(m.M_param)                                    \
                    {}                                                  \
                                                                        \
            BOOST_PARAMETER_MEMBER_FUNCTION(                            \
                                            (void),                     \
                                            setMarker,                  \
                                            tag,                        \
                                            (required                   \
                                             ( type, (std::string))		\
                                             ( name, (std::string)) )   \
                                            (optional                   \
                                             (markerAll, (bool), false) \
                                             (marker1, (bool), false)   \
                                             (marker2, (bool), false)   \
                                             (marker3, (bool), false)   \
                                             (marker4, (bool), false)   \
                                             (marker5, (bool), false)   \
                                             (marker6, (bool), false)   \
                                             (marker7, (bool), false)   \
                                             (marker8, (bool), false)   \
                                             (marker9, (bool), false)   \
                                             (marker10, (bool), false)  \
                                             (marker11, (bool), false)  \
                                             (marker12, (bool), false)  \
                                             ))                         \
                {                                                       \
                                                                        \
                    if (markerAll) {                                    \
                        marker1=true;                                   \
                        marker2=true;                                   \
                        marker3=true;                                   \
                        marker4=true;                                   \
                        marker5=true;                                   \
                        marker6=true;                                   \
                        marker7=true;                                   \
                        marker8=true;                                   \
                        marker9=true;                                   \
                        marker10=true;                                  \
                        marker11=true;                                  \
                        marker12=true;                                  \
                    }                                                   \
                                                                        \
                    std::list<marker_base_type > __listMarker = (*(M_markShape))[type][name]; \
                                                                        \
                                                                        \
                    if (type=="point")                                  \
                        {                                               \
                            BOOST_PP_FOR( (0, BOOST_PP_SUB(BOOST_PP_ARRAY_SIZE( BOOST_PP_CAT(GEOTOOL_MARKER_POINT_, \
                                                                                             GEOTOOL_SHAPE_NAME_MACRO( BOOST_PP_TUPLE_ELEM(2,0,state)))), \
                                                           1), BOOST_PP_TUPLE_ELEM(2,0,state)), \
                                          GEOTOOL_FOR_COMP1,            \
                                          GEOTOOL_FOR_INCR1,            \
                                          GEOTOOL_FOR_MARKER_POINT_MACRO) \
                                }                                       \
                    else if (type=="line")                              \
                        {                                               \
                            BOOST_PP_FOR( (0, BOOST_PP_SUB(BOOST_PP_ARRAY_SIZE( BOOST_PP_CAT(GEOTOOL_MARKER_LINE_, \
                                                                                             GEOTOOL_SHAPE_NAME_MACRO( BOOST_PP_TUPLE_ELEM(2,0,state)))), \
                                                           1), BOOST_PP_TUPLE_ELEM(2,0,state)), \
                                          GEOTOOL_FOR_COMP1,            \
                                          GEOTOOL_FOR_INCR1,            \
                                          GEOTOOL_FOR_MARKER_LINE_MACRO) \
                                }                                       \
                    else if (type=="surface")                           \
                        {                                               \
                            BOOST_PP_IF(BOOST_PP_NOT_EQUAL(GEOTOOL_SHAPE_NBSURFACE(BOOST_PP_TUPLE_ELEM(2,0,state)),0), \
                                        BOOST_PP_FOR( (0, BOOST_PP_SUB(BOOST_PP_ARRAY_SIZE( BOOST_PP_CAT(GEOTOOL_MARKER_SURFACE_, \
                                                                                                         GEOTOOL_SHAPE_NAME_MACRO( BOOST_PP_TUPLE_ELEM(2,0,state)))), \
                                                                       1), BOOST_PP_TUPLE_ELEM(2,0,state)), \
                                                      GEOTOOL_FOR_COMP1, \
                                                      GEOTOOL_FOR_INCR1, \
                                                      GEOTOOL_FOR_MARKER_SURFACE_MACRO), \
                                        )                               \
                                }                                       \
                    else if (type=="volume")                            \
                        {                                               \
                            BOOST_PP_IF(BOOST_PP_NOT_EQUAL(GEOTOOL_SHAPE_NBVOLUME(BOOST_PP_TUPLE_ELEM(2,0,state)),0), \
                                        BOOST_PP_FOR( (0, BOOST_PP_SUB(BOOST_PP_ARRAY_SIZE( BOOST_PP_CAT(GEOTOOL_MARKER_VOLUME_, \
                                                                                                         GEOTOOL_SHAPE_NAME_MACRO( BOOST_PP_TUPLE_ELEM(2,0,state)))), \
                                                                       1), BOOST_PP_TUPLE_ELEM(2,0,state)), \
                                                      GEOTOOL_FOR_COMP1, \
                                                      GEOTOOL_FOR_INCR1, \
                                                      GEOTOOL_FOR_MARKER_VOLUME_MACRO), \
                                        )                               \
                                }                                       \
                                                                        \
                    (*(M_markShape))[type][name] = __listMarker;       \
                }                                                       \
                                                                        \
                                                                        \
            std::string M_name;                                        \
                                                                        \
            static const std::string shape() { return GEOTOOL_SHAPE_NAME_STR(BOOST_PP_TUPLE_ELEM(2,0, state));} \
            const std::string name() const {return M_name;}			\
                                                                        \
                                                                        \
            std::vector<GeoTool::Node> M_param;                        \
                                                                        \
        };                                                              \
        /**/
/*_________________________________________________*/
/*_________________________________________________*/
/*                                                 */
/**/



//creation des classes representants les objets geotool
BOOST_PP_FOR( ( 0, BOOST_PP_SUB( BOOST_PP_ARRAY_SIZE( GEOTOOL_SHAPE ),1 ) ),
              GEOTOOL_FOR_COMP,
              GEOTOOL_FOR_INCR,
              GEOTOOL_SHAPE_CLASS )



template<typename mesh_type>
boost::shared_ptr<mesh_type>
createMeshFromGeoFile( std::string geofile,std::string name,double meshSize,int straighten = 1,
                       int partitions=1, WorldComm worldcomm=Environment::worldComm(),
                       int partition_file = 0, GMSH_PARTITIONER partitioner = GMSH_PARTITIONER_CHACO )
{

    boost::shared_ptr<mesh_type> mesh( new mesh_type );
    mesh->setWorldComm( worldcomm );

    if ( !worldcomm.isActive() ) return mesh;

    Gmsh gmsh( mesh_type::nDim,mesh_type::nOrder,worldcomm );
    gmsh.setCharacteristicLength( meshSize );
    gmsh.setNumberOfPartitions( partitions );
    gmsh.setPartitioner( partitioner );
    gmsh.setMshFileByPartition( partition_file );
    gmsh.setRecombine( mesh_type::shape_type::is_hypercube );

    std::ostringstream ostr;

    // preambule :
    ostr << "Mesh.MshFileVersion = " << gmsh.version() << ";\n"
         << "Mesh.CharacteristicLengthExtendFromBoundary=1;\n"
         << "Mesh.CharacteristicLengthFromPoints=1;\n"
         << "Mesh.ElementOrder=" << gmsh.order() << ";\n"
         << "Mesh.SecondOrderIncomplete = 0;\n"
         << "Mesh.Algorithm = 6;\n" // 2D mesh algorithm (1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad)
         << "Mesh.Algorithm3D = 4;\n" // 3D mesh algorithm (1=Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D)
         << "Mesh.OptimizeNetgen=1;\n"
         << "// partitioning data\n"
         << "Mesh.Partitioner=" << partitioner << ";\n"
         << "Mesh.NbPartitions=" << partitions << ";\n"
         << "Mesh.MshFilePartitioned=" << partition_file << ";\n";


    std::string contenu;
    std::ifstream ifstr( geofile.c_str(), std::ios::in );

    if ( ifstr )
    {
        // each line of the stream is appended in contenu
        while ( getline( ifstr, contenu ) )
            ostr << contenu<<"\n";

        ifstr.close();
    }


    std::string fname;
    bool generated;
    boost::tie( fname, generated ) = gmsh.generate( name,
                                                    ostr.str(),false,false,true );

    ImporterGmsh<mesh_type> import( fname, FEELPP_GMSH_FORMAT_VERSION, worldcomm );

    mesh->accept( import );
    mesh->components().set ( MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
    mesh->updateForUse();

    if ( straighten && mesh_type::nOrder > 1 )
        return straightenMesh( _mesh=mesh, _worldcomm=worldcomm.subWorldComm() );

    return mesh;
}



}//GeoTool

} //Feel
#endif /* FEELPP_GEOTOOL_HPP */
