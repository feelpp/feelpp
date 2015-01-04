/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2011-12-30

  Copyright (C) 2014 Université Joseph Fourier (Grenoble I)

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
   \file geotoolshape.cpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2014-12-30
 */

namespace Feel
{

namespace GeoTool
{

/*_________________________________________________*/

// Accessors GEOTOOL_MARKER_... macro

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




/*_________________________________________________*
 *_________________________________________________*
 * Function user                                   *
 *_________________________________________________*
 *_________________________________________________*/


void
runLine( data_geo_ptrtype dg )
{
    node_type PtA = param<0>( dg );
    node_type PtB = param<1>( dg );

    writePoint( 1, dg , PtA( 0 ), PtA( 1 ) );
    writePoint( 2, dg , PtB( 0 ), PtB( 1 ) );
    writeLine( 1, dg , 1 , 2 );
}


void
runTriangle( data_geo_ptrtype dg )
{
    node_type PtA = param<0>( dg );
    node_type PtB = param<1>( dg );
    node_type PtC = param<2>( dg );

    writePoint( 1, dg , PtA( 0 ), PtA( 1 ) );
    writePoint( 2, dg , PtB( 0 ), PtB( 1 ) );
    writePoint( 3, dg , PtC( 0 ), PtC( 1 ) );

    writeLine( 1, dg , 1 , 2 );
    writeLine( 2, dg , 2 , 3 );
    writeLine( 3, dg , 3 , 1 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3 );

    writePlaneSurface( 1, dg, 1 );
}



void
runRectangle( data_geo_ptrtype dg )
{

    node_type PtA = param<0>( dg );
    node_type PtB = param<1>( dg );

    writePoint( 1, dg , PtA( 0 ), PtA( 1 ) );
    writePoint( 2, dg , PtB( 0 ), PtA( 1 ) );
    writePoint( 3, dg , PtB( 0 ), PtB( 1 ) );
    writePoint( 4, dg , PtA( 0 ), PtB( 1 ) );

    writeLine( 1, dg , 1 , 2 );
    writeLine( 2, dg , 2 , 3 );
    writeLine( 3, dg , 3 , 4 );
    writeLine( 4, dg , 4 , 1 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3>>4 );

    writePlaneSurface( 1, dg, 1 );

}

void
runQuadrangle( data_geo_ptrtype dg )
{

    node_type PtA = param<0>( dg );
    node_type PtB = param<1>( dg );
    node_type PtC = param<2>( dg );
    node_type PtD = param<3>( dg );

    writePoint( 1, dg , PtA( 0 ), PtA( 1 ) );
    writePoint( 2, dg , PtB( 0 ), PtB( 1 ) );
    writePoint( 3, dg , PtC( 0 ), PtC( 1 ) );
    writePoint( 4, dg , PtD( 0 ), PtD( 1 ) );

    writeLine( 1, dg , 1 , 2 );
    writeLine( 2, dg , 2 , 3 );
    writeLine( 3, dg , 3 , 4 );
    writeLine( 4, dg , 4 , 1 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3>>4 );

    writePlaneSurface( 1, dg, 1 );

}

void
runPentagon( data_geo_ptrtype dg )
{

    node_type PtA = param<0>( dg );
    node_type PtB = param<1>( dg );
    node_type PtC = param<2>( dg );
    node_type PtD = param<3>( dg );
    node_type PtE = param<4>( dg );

    writePoint( 1, dg , PtA( 0 ), PtA( 1 ) );
    writePoint( 2, dg , PtB( 0 ), PtB( 1 ) );
    writePoint( 3, dg , PtC( 0 ), PtC( 1 ) );
    writePoint( 4, dg , PtD( 0 ), PtD( 1 ) );
    writePoint( 5, dg , PtE( 0 ), PtE( 1 ) );

    writeLine( 1, dg , 1 , 2 );
    writeLine( 2, dg , 2 , 3 );
    writeLine( 3, dg , 3 , 4 );
    writeLine( 4, dg , 4 , 5 );
    writeLine( 5, dg , 5 , 1 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3>>4>>5 );

    writePlaneSurface( 1, dg, 1 );

}


void
runHexagon( data_geo_ptrtype dg )
{

    node_type PtA = param<0>( dg );
    node_type PtB = param<1>( dg );
    node_type PtC = param<2>( dg );
    node_type PtD = param<3>( dg );
    node_type PtE = param<4>( dg );
    node_type PtF = param<5>( dg );

    writePoint( 1, dg , PtA( 0 ), PtA( 1 ) );
    writePoint( 2, dg , PtB( 0 ), PtB( 1 ) );
    writePoint( 3, dg , PtC( 0 ), PtC( 1 ) );
    writePoint( 4, dg , PtD( 0 ), PtD( 1 ) );
    writePoint( 5, dg , PtE( 0 ), PtE( 1 ) );
    writePoint( 6, dg , PtF( 0 ), PtF( 1 ) );

    writeLine( 1, dg , 1 , 2 );
    writeLine( 2, dg , 2 , 3 );
    writeLine( 3, dg , 3 , 4 );
    writeLine( 4, dg , 4 , 5 );
    writeLine( 5, dg , 5 , 6 );
    writeLine( 6, dg , 6 , 1 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3>>4>>5>>6 );

    writePlaneSurface( 1, dg, 1 );

}


void
runCircle( data_geo_ptrtype dg )
{
    node_type PtA = param<0>( dg );
    node_type PtB = param<1>( dg );

    writePoint( 1, dg , PtB( 0 ), PtB( 1 ) );
    writePoint( 2, dg , PtA( 0 ), PtA( 1 ) );
    writePoint( 3, dg , 2*PtA( 0 )-PtB( 0 ), 2*PtA( 1 )-PtB( 1 ) );

    writeCircle( 1, dg, 1, 2, 3 );
    writeCircle( 2, dg, 3, 2, 1 );

    writeLineLoop( 1, dg, Loop()>>1>>2 );

    writePlaneSurface( 1, dg, 1 );
    writePtInSurface( dg,2,1 );
}


void
runEllipse( data_geo_ptrtype dg )
{
    node_type PtC = param<0>( dg );
    node_type PtMinor = param<1>( dg );
    node_type PtMajor = param<2>( dg );

    writePoint( 1, dg , PtC( 0 ), PtC( 1 ) );
    writePoint( 2, dg , PtMinor( 0 ), PtMinor( 1 ) );
    writePoint( 3, dg , 2*PtC( 0 )-PtMinor( 0 ), 2*PtC( 1 )-PtMinor( 1 ) );
    writePoint( 4, dg , PtMajor( 0 ), PtMajor( 1 ) );
    writePoint( 5, dg , 2*PtC( 0 )-PtMajor( 0 ), 2*PtC( 1 )-PtMajor( 1 ) );

    writeEllipse( 1, dg, 2, 1, 4, 4 );
    writeEllipse( 2, dg, 4, 1, 3, 3 );
    writeEllipse( 3, dg, 3, 1, 5, 5 );
    writeEllipse( 4, dg, 5, 1, 2, 2 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3>>4 );

    writePlaneSurface( 1, dg, 1 );
    writePtInSurface( dg,2,1 );
}



void
runPie( data_geo_ptrtype dg )
{
    node_type PtA = param<0>( dg );
    node_type PtB = param<1>( dg );
    node_type PtC = param<2>( dg );

    writePoint( 1, dg , PtB( 0 ), PtB( 1 ) );
    writePoint( 2, dg , PtA( 0 ), PtA( 1 ) );
    writePoint( 3, dg , PtC( 0 ), PtC( 1 ) );

    writeCircle( 1, dg, 1, 2, 3 );
    writeLine( 2, dg , 3 , 2 );
    writeLine( 3, dg , 2 , 1 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3 );

    writePlaneSurface( 1, dg, 1 );
    writePtInSurface( dg,2,1 );

}

void
runSpecial_1a( data_geo_ptrtype dg )
{
    double yh=1.0;
    Node a1( 0.0, yh );
    Node a2( 3.0, yh );
    Node a3( 6.0, yh+0.7 );
    Node a4( 7.3, yh-0.5 );
    Node a5( 8.5, yh );
    Node a6( 11.0, yh );
    double ep=0.3;
    //_______________________________________________//
    writePoint( 1, dg , a1( 0 ), a1( 1 ) );
    writePoint( 2, dg , a2( 0 ), a2( 1 ) );
    writePoint( 3, dg , a3( 0 ), a3( 1 ) );
    writePoint( 4, dg , a4( 0 ), a4( 1 ) );
    writePoint( 5, dg , a5( 0 ), a5( 1 ) );
    writePoint( 6, dg , a6( 0 ), a6( 1 ) );

    writeSpline( 1, dg, Loop()>>1>>2>>3>>4>>5>>6 );

    writePoint( 7, dg , a1( 0 ), a1( 1 )+ep );
    writePoint( 8, dg , a2( 0 ), a2( 1 )+ep );
    writePoint( 9, dg , a3( 0 ), a3( 1 )+ep );
    writePoint( 10, dg , a4( 0 ), a4( 1 )+ep );
    writePoint( 11, dg , a5( 0 ), a5( 1 )+ep );
    writePoint( 12, dg , a6( 0 ), a6( 1 )+ep );

    writeSpline( 2, dg, Loop()>>7>>8>>9>>10>>11>>12 );

    writeLine( 3, dg, 1,7 );
    writeLine( 4, dg, 6,12 );

    writeLineLoop( 1, dg, Loop()>>1>>4>>-2>>-3 );
    writePlaneSurface( 1, dg, 1 );

    writePoint( 13, dg , a1( 0 ), -a1( 1 ) );
    writePoint( 14, dg , a2( 0 ), -a2( 1 ) );
    writePoint( 15, dg , a3( 0 ), -a3( 1 ) );
    writePoint( 16, dg , a4( 0 ), -a4( 1 ) );
    writePoint( 17, dg , a5( 0 ), -a5( 1 ) );
    writePoint( 18, dg , a6( 0 ), -a6( 1 ) );

    writeSpline( 5, dg, Loop()>>13>>14>>15>>16>>17>>18 );

    writePoint( 19, dg , a1( 0 ), -a1( 1 )-ep );
    writePoint( 20, dg , a2( 0 ), -a2( 1 )-ep );
    writePoint( 21, dg , a3( 0 ), -a3( 1 )-ep );
    writePoint( 22, dg , a4( 0 ), -a4( 1 )-ep );
    writePoint( 23, dg , a5( 0 ), -a5( 1 )-ep );
    writePoint( 24, dg , a6( 0 ), -a6( 1 )-ep );

    writeSpline( 6, dg, Loop()>>19>>20>>21>>22>>23>>24 );

    writeLine( 7, dg, 13,19 );
    writeLine( 8, dg, 18,24 );

    writeLineLoop( 2, dg, Loop()>>5>>8>>-6>>-7 );
    writePlaneSurface( 2, dg, 2 );

}

void
runSpecial_1b( data_geo_ptrtype dg )
{
    double yh=1.0;
    Node a1( 0.0, yh );
    Node a2( 3.0, yh );
    Node a3( 6.0, yh+0.7 );
    Node a4( 7.3, yh-0.5 );
    Node a5( 8.5, yh );
    Node a6( /*11.0*/16.0, yh );
    //_______________________________________________//
    writePoint( 1, dg , a1( 0 ), a1( 1 ) );
    writePoint( 2, dg , a2( 0 ), a2( 1 ) );
    writePoint( 3, dg , a3( 0 ), a3( 1 ) );
    writePoint( 4, dg , a4( 0 ), a4( 1 ) );
    writePoint( 5, dg , a5( 0 ), a5( 1 ) );
    writePoint( 6, dg , a6( 0 ), a6( 1 ) );

    writeSpline( 1, dg, Loop()>>1>>2>>3>>4>>5>>6 );

    writePoint( 7, dg , a1( 0 ), -a1( 1 ) );
    writePoint( 8, dg , a2( 0 ), -a2( 1 ) );
    writePoint( 9, dg , a3( 0 ), -a3( 1 ) );
    writePoint( 10, dg , a4( 0 ), -a4( 1 ) );
    writePoint( 11, dg , a5( 0 ), -a5( 1 ) );
    writePoint( 12, dg , a6( 0 ), -a6( 1 ) );

    writeSpline( 2, dg, Loop()>>7>>8>>9>>10>>11>>12 );

    writeLine( 3, dg, 1,7 );
    writeLine( 4, dg, 6,12 );

    writeLineLoop( 1, dg, Loop()>>1>>4>>-2>>-3 );

    writePlaneSurface( 1, dg, 1 );

}

void
runPeanut( data_geo_ptrtype dg )
{
    node_type PtCenter = param<0>( dg );
    node_type majorRadiusParam = param<1>( dg ); //2
    node_type minorRadiusParam = param<2>( dg ); //1
    node_type penautRadiusParam = param<3>( dg ); //0.1
    double majorRadius = majorRadiusParam( 0 );
    double minorRadius = minorRadiusParam( 0 );
    double penautRadius = penautRadiusParam( 0 );
    writePoint( 1, dg , PtCenter( 0 )               , PtCenter( 1 )+penautRadius, 0. );
    writePoint( 2, dg , PtCenter( 0 )               , PtCenter( 1 )-penautRadius, 0. );
    writePoint( 3, dg , PtCenter( 0 )-majorRadius   , PtCenter( 1 )             , 0. );
    writePoint( 4, dg , PtCenter( 0 )+majorRadius   , PtCenter( 1 )             , 0. );
    writePoint( 5, dg , PtCenter( 0 )-majorRadius/4., PtCenter( 1 )+minorRadius , 0. );
    writePoint( 6, dg , PtCenter( 0 )-majorRadius/4., PtCenter( 1 )-minorRadius , 0. );
    writePoint( 7, dg , PtCenter( 0 )+majorRadius/4., PtCenter( 1 )+minorRadius , 0. );
    writePoint( 8, dg , PtCenter( 0 )+majorRadius/4., PtCenter( 1 )-minorRadius , 0. );
    writeBSpline( 1, dg, Loop()>>1>>5>>3>>6>>2>>8>>4>>7>>1 );

    writeLineLoop( 1, dg, Loop()>>1 );
    writePlaneSurface( 1, dg, 1 );
}

void
runHexahedron( data_geo_ptrtype dg )
{

    node_type Pt1 = param<0>( dg );
    node_type Pt2 = param<1>( dg );
    node_type Pt3 = param<2>( dg );
    node_type Pt4 = param<3>( dg );
    node_type Pt5 = param<4>( dg );
    node_type Pt6 = param<5>( dg );
    node_type Pt7 = param<6>( dg );
    node_type Pt8 = param<7>( dg );

    writePoint( 1, dg , Pt1( 0 ), Pt1( 1 ), Pt1( 2 ) );
    writePoint( 2, dg , Pt2( 0 ), Pt2( 1 ), Pt2( 2 ) );
    writePoint( 3, dg , Pt3( 0 ), Pt3( 1 ), Pt3( 2 ) );
    writePoint( 4, dg , Pt4( 0 ), Pt4( 1 ), Pt4( 2 ) );
    writePoint( 5, dg , Pt5( 0 ), Pt5( 1 ), Pt5( 2 ) );
    writePoint( 6, dg , Pt6( 0 ), Pt6( 1 ), Pt6( 2 ) );
    writePoint( 7, dg , Pt7( 0 ), Pt7( 1 ), Pt7( 2 ) );
    writePoint( 8, dg , Pt8( 0 ), Pt8( 1 ), Pt8( 2 ) );

    writeLine( 1, dg , 1 , 2 );
    writeLine( 2, dg , 2 , 3 );
    writeLine( 3, dg , 3 , 4 );
    writeLine( 4, dg , 4 , 1 );
    writeLine( 5, dg , 5 , 6 );
    writeLine( 6, dg , 6 , 7 );
    writeLine( 7, dg , 7 , 8 );
    writeLine( 8, dg , 8 , 5 );
    writeLine( 9, dg , 1 , 5 );
    writeLine( 10, dg , 2 , 6 );
    writeLine( 11, dg , 3 , 7 );
    writeLine( 12, dg , 4 , 8 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3>>4 );
    writePlaneSurface( 1, dg, 1 );
    writeLineLoop( 2, dg, Loop()>>-5>>-8>>-7>>-6 );
    writePlaneSurface( 2, dg, 2 );
    writeLineLoop( 3, dg, Loop()>>-1>>9>>5>>-10 );
    writePlaneSurface( 3, dg, 3 );
    writeLineLoop( 4, dg, Loop()>>10>>6>>-11>>-2 );
    writePlaneSurface( 4, dg, 4 );
    writeLineLoop( 5, dg, Loop()>>11>>7>>-12>>-3 );
    writePlaneSurface( 5, dg, 5 );
    writeLineLoop( 6, dg, Loop()>>8>>-9>>-4>>12 );
    writePlaneSurface( 6, dg, 6 );

    writeSurfaceLoop( 1, dg, Loop()>>1>>2>>3>>4>>5>>6 );
    writeVolume( 1, dg, 1 );
}

void
runTetrahedron( data_geo_ptrtype dg )
{
    node_type Pt1 = param<0>( dg );
    node_type Pt2 = param<1>( dg );
    node_type Pt3 = param<2>( dg );
    node_type Pt4 = param<3>( dg );

    writePoint( 1, dg , Pt1( 0 ), Pt1( 1 ), Pt1( 2 ) );
    writePoint( 2, dg , Pt2( 0 ), Pt2( 1 ), Pt2( 2 ) );
    writePoint( 3, dg , Pt3( 0 ), Pt3( 1 ), Pt3( 2 ) );
    writePoint( 4, dg , Pt4( 0 ), Pt4( 1 ), Pt4( 2 ) );

    writeLine( 1, dg , 1 , 2 );
    writeLine( 2, dg , 2 , 3 );
    writeLine( 3, dg , 3 , 1 );
    writeLine( 4, dg , 1 , 4 );
    writeLine( 5, dg , 2 , 4 );
    writeLine( 6, dg , 3 , 4 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3 );
    writePlaneSurface( 1, dg, 1 );

    //writeLineLoop( 2, dg, Loop()>>5>>-4>>1 );
    writeLineLoop( 2, dg, Loop()>>-1>>4>>-5 );
    writePlaneSurface( 2, dg, 2 );
    //writeLineLoop( 3, dg, Loop()>>2>>6>>-5 );
    writeLineLoop( 3, dg, Loop()>>5>>-6>>-2 );
    writePlaneSurface( 3, dg, 3 );
    //writeLineLoop( 4, dg, Loop()>>3>>4>>-6 );
    writeLineLoop( 4, dg, Loop()>>6>>-4>>-3 );
    writePlaneSurface( 4, dg, 4 );

    writeSurfaceLoop( 1, dg, Loop()>>1>>2>>3>>4 );

    writeVolume( 1, dg, 1 );
}


void
runCube( data_geo_ptrtype dg )
{

    node_type Pt1 = param<0>( dg );
    node_type Pt7 = param<1>( dg );
    double lenx = Pt7( 0 )-Pt1( 0 );
    double leny = Pt7( 1 )-Pt1( 1 );
    double lenz = Pt7( 2 )-Pt1( 2 );
    node_type Pt2 = Pt1;
    Pt2( 0 )+=lenx;
    node_type Pt3 = Pt1;
    Pt3( 0 )+=lenx;
    Pt3( 1 )+=leny;
    node_type Pt4 = Pt1;
    Pt4( 1 )+=leny;
    node_type Pt5 = Pt1;
    Pt5( 2 )+=lenz;
    node_type Pt6 = Pt1;
    Pt6( 0 )+=lenx;
    Pt6( 2 )+=lenz;
    // Pt7 was give above by the user
    node_type Pt8 = Pt1;
    Pt8( 1 )+=leny;
    Pt8( 2 )+=lenz;

    writePoint( 1, dg , Pt1( 0 ), Pt1( 1 ), Pt1( 2 ) );
    writePoint( 2, dg , Pt2( 0 ), Pt2( 1 ), Pt2( 2 ) );
    writePoint( 3, dg , Pt3( 0 ), Pt3( 1 ), Pt3( 2 ) );
    writePoint( 4, dg , Pt4( 0 ), Pt4( 1 ), Pt4( 2 ) );
    writePoint( 5, dg , Pt5( 0 ), Pt5( 1 ), Pt5( 2 ) );
    writePoint( 6, dg , Pt6( 0 ), Pt6( 1 ), Pt6( 2 ) );
    writePoint( 7, dg , Pt7( 0 ), Pt7( 1 ), Pt7( 2 ) );
    writePoint( 8, dg , Pt8( 0 ), Pt8( 1 ), Pt8( 2 ) );

    writeLine( 1, dg , 1 , 2 );
    writeLine( 2, dg , 2 , 3 );
    writeLine( 3, dg , 3 , 4 );
    writeLine( 4, dg , 4 , 1 );
    writeLine( 5, dg , 5 , 6 );
    writeLine( 6, dg , 6 , 7 );
    writeLine( 7, dg , 7 , 8 );
    writeLine( 8, dg , 8 , 5 );
    writeLine( 9, dg , 1 , 5 );
    writeLine( 10, dg , 2 , 6 );
    writeLine( 11, dg , 3 , 7 );
    writeLine( 12, dg , 4 , 8 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3>>4 );
    writePlaneSurface( 1, dg, 1 );
    writeLineLoop( 2, dg, Loop()>>-5>>-8>>-7>>-6 );
    writePlaneSurface( 2, dg, 2 );
    writeLineLoop( 3, dg, Loop()>>-1>>9>>5>>-10 );
    writePlaneSurface( 3, dg, 3 );
    writeLineLoop( 4, dg, Loop()>>10>>6>>-11>>-2 );
    writePlaneSurface( 4, dg, 4 );
    writeLineLoop( 5, dg, Loop()>>11>>7>>-12>>-3 );
    writePlaneSurface( 5, dg, 5 );
    writeLineLoop( 6, dg, Loop()>>8>>-9>>-4>>12 );
    writePlaneSurface( 6, dg, 6 );

    writeSurfaceLoop( 1, dg, Loop()>>1>>2>>3>>4>>5>>6 );

    writeVolume( 1, dg, 1 );
}

void
runCylindre( data_geo_ptrtype dg )
{
#if 0
    node_type centre = param<0>( dg );
    double par = param<1>( dg )( 0 );

    Node c1( -1,-par,0 );
    Node c2( -1, 0,-par );
    Node c3( -1, par,0 );
    Node c4( -1, 0, par );
    writePoint( 1, dg , centre( 0 ), centre( 1 ), centre( 2 ) );
    writePoint( 2, dg , c1( 0 ), c1( 1 ), c1( 2 ) );
    writePoint( 3, dg , c2( 0 ), c2( 1 ), c2( 2 ) );
    writePoint( 4, dg , c3( 0 ), c3( 1 ), c3( 2 ) );
    writePoint( 5, dg , c4( 0 ), c4( 1 ), c4( 2 ) );

    writeCircle( 1, dg, 2, 1, 3 );
    writeCircle( 2, dg, 3, 1, 4 );
    writeCircle( 3, dg, 4, 1, 5 );
    writeCircle( 4, dg, 5, 1, 2 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3>>4 );
    writePlaneSurface( 1, dg, 1 );

    //volume 1, extrude surface 1, creation surface 2,3,4,5,6
    writeExtrudeSurface( 1,dg,1,Loop()>>2>>3>>4>>5>>6 );
#endif

    node_type centre = param<0>( dg );
    node_type direction = param<1>( dg );
    double rayon = param<2>( dg )( 0 );
    double longueur=param<3>( dg )( 0 );

    auto basis = computeBasisOrthogonal( direction,centre );
    Node dir( boost::get<0>( basis ) );
    Node u( boost::get<1>( basis ) );
    Node v( boost::get<2>( basis ) );

    Node geoX1( centre( 0 )+u( 0 )*rayon,
                centre( 1 )+u( 1 )*rayon,
                centre( 2 )+u( 2 )*rayon  );

    Node geoX2( centre( 0 )+v( 0 )*rayon,
                centre( 1 )+v( 1 )*rayon,
                centre( 2 )+v( 2 )*rayon  );

    Node geoX3( centre( 0 )-u( 0 )*rayon,
                centre( 1 )-u( 1 )*rayon,
                centre( 2 )-u( 2 )*rayon  );

    Node geoX4( centre( 0 )-v( 0 )*rayon,
                centre( 1 )-v( 1 )*rayon,
                centre( 2 )-v( 2 )*rayon  );

    Node centre2( centre( 0 )+longueur*dir( 0 ),
                  centre( 1 )+longueur*dir( 1 ),
                  centre( 2 )+longueur*dir( 2 ) );

    Node geoX5( geoX1( 0 )+longueur*dir( 0 ),
                geoX1( 1 )+longueur*dir( 1 ),
                geoX1( 2 )+longueur*dir( 2 ) );

    Node geoX6( geoX2( 0 )+longueur*dir( 0 ),
                geoX2( 1 )+longueur*dir( 1 ),
                geoX2( 2 )+longueur*dir( 2 ) );

    Node geoX7( geoX3( 0 )+longueur*dir( 0 ),
                geoX3( 1 )+longueur*dir( 1 ),
                geoX3( 2 )+longueur*dir( 2 ) );

    Node geoX8( geoX4( 0 )+longueur*dir( 0 ),
                geoX4( 1 )+longueur*dir( 1 ),
                geoX4( 2 )+longueur*dir( 2 ) );

    //--------------------------------------------------------------------------//

    writePoint( 1, dg, centre( 0 ), centre( 1 ) ,centre( 2 ) );
    writePoint( 2, dg, geoX1( 0 ), geoX1( 1 ) ,geoX1( 2 ) );
    writePoint( 3, dg, geoX2( 0 ), geoX2( 1 ) ,geoX2( 2 ) );
    writePoint( 4, dg, geoX3( 0 ), geoX3( 1 ) ,geoX3( 2 ) );
    writePoint( 5, dg, geoX4( 0 ), geoX4( 1 ) ,geoX4( 2 ) );

    writeCircle( 1, dg, 2,1,3 );
    writeCircle( 2, dg, 3,1,4 );
    writeCircle( 3, dg, 4,1,5 );
    writeCircle( 4, dg, 5,1,2 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3>>4 );

    writePoint( 6, dg, centre2( 0 ), centre2( 1 ) ,centre2( 2 ) );
    writePoint( 7, dg, geoX5( 0 ), geoX5( 1 ) ,geoX5( 2 ) );
    writePoint( 8, dg, geoX6( 0 ), geoX6( 1 ) ,geoX6( 2 ) );
    writePoint( 9, dg, geoX7( 0 ), geoX7( 1 ) ,geoX7( 2 ) );
    writePoint( 10, dg, geoX8( 0 ), geoX8( 1 ) ,geoX8( 2 ) );

    writeCircle( 5, dg, 7,6,8 );
    writeCircle( 6, dg,8,6,9 );
    writeCircle( 7, dg,9,6,10 );
    writeCircle( 8 ,dg,10,6,7 );

    writeLineLoop( 2, dg, Loop()>>-5>>-8>>-7>>-6 );

    writeLine( 9, dg, 4, 9 );
    writeLine( 10, dg, 5, 10 );
    writeLine( 11, dg, 2, 7 );
    writeLine( 12, dg, 8, 3 );
    writeLineLoop( 13, dg, Loop()>>-9>>-2>>-12>>6 );
    writeLineLoop( 15, dg, Loop()>>9>>7>>-10>>-3 );
    writeLineLoop( 17, dg, Loop()>>10>>8>>-11>>-4 );
    writeLineLoop( 19, dg, Loop()>>11>>5>>12>>-1 );

    writePlaneSurface( 1, dg, 1 );
    writePlaneSurface( 2, dg, 2 );
    writePtInSurface(dg,1,1 );
    writePtInSurface(dg,6,2 );

    writeRuledSurface( 3, dg, 13 );
    writeRuledSurface( 4, dg, 15 );
    writeRuledSurface( 5, dg, 17 );
    writeRuledSurface( 6, dg, 19 );

    writeSurfaceLoop( 1, dg, Loop()>>5>>4>>3>>2>>6>>1 );
    writeVolume( 1, dg, 1 );

}

void
runSphere( data_geo_ptrtype dg )
{

    node_type centre = param<0>( dg );
    double R = param<1>( dg )( 0 );

    double x_center = centre( 0 );
    double y_center = centre( 1 );
    double z_center = centre( 2 );

    writePoint( 1, dg, x_center, y_center, z_center );
    writePoint( 2, dg, x_center - R, y_center, z_center );
    writePoint( 4, dg, x_center, y_center - R, z_center );
    writePoint( 5, dg, x_center + R, y_center, z_center );
    writePoint( 8, dg, x_center, y_center, z_center - R );
    writePoint( 11, dg, x_center, y_center + R, z_center );
    writePoint( 14, dg, x_center, y_center, z_center + R );
    writeCircle( 1, dg, 2, 1, 4 );
    writeCircle( 2, dg, 4, 1, 5 );
    writeCircle( 3, dg, 2, 1, 8 );
    writeCircle( 4, dg, 4, 1, 8 );
    writeCircle( 6, dg, 2, 1, 11 );
    writeCircle( 7, dg, 8, 1, 11 );
    writeCircle( 9, dg, 2, 1, 14 );
    writeCircle( 10, dg, 11, 1, 14 );
    writeCircle( 13, dg, 14, 1, 4 );
    writeCircle( 15, dg, 8, 1, 5 );
    writeCircle( 18, dg, 11, 1, 5 );
    writeCircle( 21, dg, 14, 1, 5 );

#if 0
    writeLineLoop( 1, dg, Loop()>>1>>4>>-3 );
    writeRuledSurface( 1, dg,1 );
    writeLineLoop( 2, dg, Loop()>>3>>7>>-6 );
    writeRuledSurface( 2, dg,2 );
    writeLineLoop( 3, dg, Loop()>>6>>10>>-9 );
    writeRuledSurface( 3,dg,3 );
    writeLineLoop( 4, dg, Loop()>>9>>13>>-1 );
    writeRuledSurface( 4, dg, 4 );
    writeLineLoop( 5, dg, Loop()>>-15>>-4>>2 );
    writeRuledSurface( 5, dg, 5 );
    writeLineLoop( 6, dg, Loop()>>-18>>-7>>15 );
    writeRuledSurface( 6, dg, 6 );
    writeLineLoop( 7, dg, Loop()>>-21>>-10>>18 );
    writeRuledSurface( 7,dg,7 );
    writeLineLoop( 8, dg, Loop()>>-2>>-13>>21 );
    writeRuledSurface( 8, dg,8 );
    writeSurfaceLoop( 1, dg, Loop()>>1>>4>>3>>2>>6>>7>>8>>5 );

#else
    writeLineLoop( 1, dg, Loop()>>1>>4>>-3 );
    writeRuledSurface( 1, dg,1 );
    writeLineLoop( 2, dg, Loop()>>7>>-6>>3 );//
    writeRuledSurface( 2, dg,2 );
    writeLineLoop( 3, dg, Loop()>>6>>10>>-9 );//
    writeRuledSurface( 3,dg,3 );
    writeLineLoop( 4, dg, Loop()>>13>>-1>>9 );//
    writeRuledSurface( 4, dg, 4 );
    writeLineLoop( 5, dg, Loop()>>-4>>2>>-15 );//
    writeRuledSurface( 5, dg, 5 );
    writeLineLoop( 6, dg, Loop()>>-18>>-7>>15 );//
    writeRuledSurface( 6, dg, 6 );
    writeLineLoop( 7, dg, Loop()>>-10>>18>>-21 );//
    writeRuledSurface( 7,dg,7 );
    writeLineLoop( 8, dg, Loop()>>-2>>-13>>21 );
    writeRuledSurface( 8, dg,8 );
    writeSurfaceLoop( 1, dg, Loop()>>3>>2>>6>>7>>8>>5>>1>>4 );
#endif

    writeVolume( 1, dg, 1 );
}



void
runTube( data_geo_ptrtype dg )
{

    node_type centre = param<0>( dg );
    node_type direction = param<1>( dg );
    double rayon = param<2>( dg )( 0 );
    double longueur=param<3>( dg )( 0 );
    double epaisseur=param<4>( dg )( 0 );

    auto basis = computeBasisOrthogonal( direction,centre );
    Node dir( boost::get<0>( basis ) );
    Node u( boost::get<1>( basis ) );
    Node v( boost::get<2>( basis ) );

    Node geoX1( centre( 0 )+u( 0 )*rayon,
                centre( 1 )+u( 1 )*rayon,
                centre( 2 )+u( 2 )*rayon  );

    Node geoX2( centre( 0 )+v( 0 )*rayon,
                centre( 1 )+v( 1 )*rayon,
                centre( 2 )+v( 2 )*rayon  );

    Node geoX3( centre( 0 )-u( 0 )*rayon,
                centre( 1 )-u( 1 )*rayon,
                centre( 2 )-u( 2 )*rayon  );

    Node geoX4( centre( 0 )-v( 0 )*rayon,
                centre( 1 )-v( 1 )*rayon,
                centre( 2 )-v( 2 )*rayon  );

    Node centre2( centre( 0 )+longueur*dir( 0 ),
                  centre( 1 )+longueur*dir( 1 ),
                  centre( 2 )+longueur*dir( 2 ) );

    Node geoX5( geoX1( 0 )+longueur*dir( 0 ),
                geoX1( 1 )+longueur*dir( 1 ),
                geoX1( 2 )+longueur*dir( 2 ) );

    Node geoX6( geoX2( 0 )+longueur*dir( 0 ),
                geoX2( 1 )+longueur*dir( 1 ),
                geoX2( 2 )+longueur*dir( 2 ) );

    Node geoX7( geoX3( 0 )+longueur*dir( 0 ),
                geoX3( 1 )+longueur*dir( 1 ),
                geoX3( 2 )+longueur*dir( 2 ) );

    Node geoX8( geoX4( 0 )+longueur*dir( 0 ),
                geoX4( 1 )+longueur*dir( 1 ),
                geoX4( 2 )+longueur*dir( 2 ) );


    Node geoX1B( centre( 0 )+u( 0 )*( rayon+epaisseur ),
                 centre( 1 )+u( 1 )*( rayon+epaisseur ),
                 centre( 2 )+u( 2 )*( rayon+epaisseur )  );

    Node geoX2B( centre( 0 )+v( 0 )*( rayon+epaisseur ),
                 centre( 1 )+v( 1 )*( rayon+epaisseur ),
                 centre( 2 )+v( 2 )*( rayon+epaisseur ) );

    Node geoX3B( centre( 0 )-u( 0 )*( rayon+epaisseur ),
                 centre( 1 )-u( 1 )*( rayon+epaisseur ),
                 centre( 2 )-u( 2 )*( rayon+epaisseur ) );

    Node geoX4B( centre( 0 )-v( 0 )*( rayon+epaisseur ),
                 centre( 1 )-v( 1 )*( rayon+epaisseur ),
                 centre( 2 )-v( 2 )*( rayon+epaisseur ) );

    Node geoX5B( geoX1B( 0 )+longueur*dir( 0 ),
                 geoX1B( 1 )+longueur*dir( 1 ),
                 geoX1B( 2 )+longueur*dir( 2 ) );

    Node geoX6B( geoX2B( 0 )+longueur*dir( 0 ),
                 geoX2B( 1 )+longueur*dir( 1 ),
                 geoX2B( 2 )+longueur*dir( 2 ) );

    Node geoX7B( geoX3B( 0 )+longueur*dir( 0 ),
                 geoX3B( 1 )+longueur*dir( 1 ),
                 geoX3B( 2 )+longueur*dir( 2 ) );

    Node geoX8B( geoX4B( 0 )+longueur*dir( 0 ),
                 geoX4B( 1 )+longueur*dir( 1 ),
                 geoX4B( 2 )+longueur*dir( 2 ) );

    //--------------------------------------------------------------------------//

    writePoint( 1, dg, centre( 0 ), centre( 1 ) ,centre( 2 ) );
    writePoint( 2, dg, geoX1( 0 ), geoX1( 1 ) ,geoX1( 2 ) );
    writePoint( 3, dg, geoX2( 0 ), geoX2( 1 ) ,geoX2( 2 ) );
    writePoint( 4, dg, geoX3( 0 ), geoX3( 1 ) ,geoX3( 2 ) );
    writePoint( 5, dg, geoX4( 0 ), geoX4( 1 ) ,geoX4( 2 ) );

    writeCircle( 1, dg, 2,1,3 );
    writeCircle( 2, dg, 3,1,4 );
    writeCircle( 3, dg, 4,1,5 );
    writeCircle( 4, dg, 5,1,2 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3>>4 );

    writePoint( 6, dg, centre2( 0 ), centre2( 1 ) ,centre2( 2 ) );
    writePoint( 7, dg, geoX5( 0 ), geoX5( 1 ) ,geoX5( 2 ) );
    writePoint( 8, dg, geoX6( 0 ), geoX6( 1 ) ,geoX6( 2 ) );
    writePoint( 9, dg, geoX7( 0 ), geoX7( 1 ) ,geoX7( 2 ) );
    writePoint( 10, dg, geoX8( 0 ), geoX8( 1 ) ,geoX8( 2 ) );

    writeCircle( 5, dg, 7,6,8 );
    writeCircle( 6, dg,8,6,9 );
    writeCircle( 7, dg,9,6,10 );
    writeCircle( 8 ,dg,10,6,7 );

    writeLineLoop( 2, dg,Loop()>>5>>6>>7>>8 );

    writeLine( 9, dg, 4, 9 );
    writeLine( 10, dg, 5, 10 );
    writeLine( 11, dg, 2, 7 );
    writeLine( 12, dg, 8, 3 );

    writeLineLoop( 3, dg, Loop()>>9>>-6>>12>>2 );
    writeLineLoop( 4, dg, Loop()>>9>>7>>-10>>-3 );
    writeLineLoop( 5, dg, Loop()>>10>>8>>-11>>-4 );
    writeLineLoop( 6, dg, Loop()>>11>>5>>12>>-1 );

    //writePlaneSurface(1, dg, 1);
    //writePlaneSurface(2, dg, 2);

    writeRuledSurface( 3, dg, 3 );
    writeRuledSurface( 4, dg, 4 );
    writeRuledSurface( 5, dg, 5 );
    writeRuledSurface( 6, dg, 6 );


    writePoint( 11, dg, geoX1B( 0 ), geoX1B( 1 ) ,geoX1B( 2 ) );
    writePoint( 12, dg, geoX2B( 0 ), geoX2B( 1 ) ,geoX2B( 2 ) );
    writePoint( 13, dg, geoX3B( 0 ), geoX3B( 1 ) ,geoX3B( 2 ) );
    writePoint( 14, dg, geoX4B( 0 ), geoX4B( 1 ) ,geoX4B( 2 ) );

    writeCircle( 13, dg, 11,1,12 );
    writeCircle( 14, dg, 12,1,13 );
    writeCircle( 15, dg, 13,1,14 );
    writeCircle( 16, dg, 14,1,11 );

    writePoint( 15, dg, geoX5B( 0 ), geoX5B( 1 ) ,geoX5B( 2 ) );
    writePoint( 16, dg, geoX6B( 0 ), geoX6B( 1 ) ,geoX6B( 2 ) );
    writePoint( 17, dg, geoX7B( 0 ), geoX7B( 1 ) ,geoX7B( 2 ) );
    writePoint( 18, dg, geoX8B( 0 ), geoX8B( 1 ) ,geoX8B( 2 ) );

    writeCircle( 17, dg, 15,6,16 );
    writeCircle( 18, dg, 16,6,17 );
    writeCircle( 19, dg, 17,6,18 );
    writeCircle( 20 ,dg, 18,6,15 );

    writeLine( 21, dg, 13, 17 );
    writeLine( 22, dg, 14, 18 );
    writeLine( 23, dg, 11, 15 );
    writeLine( 24, dg, 16, 12 );

    writeLineLoop( 7, dg, Loop()>>21>>-18>>24>>14 );
    writeLineLoop( 8, dg, Loop()>>21>>19>>-22>>-15 );
    writeLineLoop( 9, dg, Loop()>>22>>20>>-23>>-16 );
    writeLineLoop( 10, dg, Loop()>>23>>17>>24>>-13 );

    writeRuledSurface( 7, dg, 7 );
    writeRuledSurface( 8, dg, 8 );
    writeRuledSurface( 9, dg, 9 );
    writeRuledSurface( 10, dg, 10 );

    writeLine( 25, dg, 2, 11 );
    writeLine( 26, dg, 3, 12 );
    writeLine( 27, dg, 4, 13 );
    writeLine( 28, dg, 5, 14 );
    writeLine( 29, dg, 7, 15 );
    writeLine( 30, dg, 8, 16 );
    writeLine( 31, dg, 9, 17 );
    writeLine( 32, dg, 10, 18 );


    writeLineLoop( 11, dg, Loop()>>19>>-32>>-7>>31 );
    writeLineLoop( 12, dg, Loop()>>20>>-29>>-8>>32 );
    writeLineLoop( 13, dg, Loop()>>6>>31>>-18>>-30 );
    writeLineLoop( 14, dg, Loop()>>5>>30>>-17>>-29 );
    writeLineLoop( 15, dg, Loop()>>15>>-28>>-3>>27 );
    writeLineLoop( 16, dg, Loop()>>27>>-14>>-26>>2 );
    writeLineLoop( 17, dg, Loop()>>26>>-13>>-25>>1 );
    writeLineLoop( 18, dg, Loop()>>25>>-16>>-28>>4 );

    // internal surface
    writeLineLoop( 19 ,dg, Loop()>>23>>-29>>-11>>25 );
    writeLineLoop( 20 ,dg, Loop()>>10>>32>>-22>>-28 );
    writeLineLoop( 21 ,dg, Loop()>>9>>31>>-21>>-27 );
    writeLineLoop( 22 ,dg, Loop()>>24>>-26>>-12>>30 );

    // internal surface
    writePlaneSurface( 19,dg,19 );
    writeRuledSurface( 20,dg,20 );
    writePlaneSurface( 21,dg,21 );
    writePlaneSurface( 22,dg,22 );

    // Inlet or outlet?
    writeRuledSurface( 11,dg,11 );
    writeRuledSurface( 12,dg,12 );
    writeRuledSurface( 13,dg,13 );
    writeRuledSurface( 14,dg,14 );

    // Inlet or outlet?
    writeRuledSurface( 15,dg,15 );
    writeRuledSurface( 16,dg,16 );
    writeRuledSurface( 17,dg,17 );
    writeRuledSurface( 18,dg,18 );



    //writeSurfaceLoop(1, dg, Loop()>>3>>4>>5>>6>>7>>8>>9>>10>>11>>12>>13>>14>>15>>16>>17>>18);
    //writeSurfaceLoop(1, dg, Loop()>>3>>4>>5>>6>>7>>8>>9>>10>>11>>12>>13>>14>>15>>16>>17>>18>>19>>20>>21>>22);
    //writeVolume(1, dg, 1);
#if 0
    writeSurfaceLoop( 1, dg, Loop()>>6>>9>>2>>13>>18>>19 );
    writeVolume( 1,dg,1 );
    writeSurfaceLoop( 2, dg, Loop()>>18>>7>>10>>3>>16>>17 );
    writeVolume( 2,dg,2 );
    writeSurfaceLoop( 37, dg, Loop()>>20>>12>>4>>15>>8>>17 );
    writeVolume( 3,dg,3 );
    writeSurfaceLoop( 4, dg, Loop()>>11>>1>>14>>5>>19>>20 );
    writeVolume( 4,dg,4 );
#else
    writeSurfaceLoop( 1, dg, Loop()>>4>>16>>8>>19>>12>>9 );
    writeVolume( 1,dg,1 );
    writeSurfaceLoop( 2, dg, Loop()>>5>>15>>1>>18>>11>>12 );
    writeVolume( 2,dg,2 );
    writeSurfaceLoop( 3, dg, Loop()>>6>>13>>2>>17>>11>>10 );
    writeVolume( 3,dg,3 );
    writeSurfaceLoop( 4, dg, Loop()>>10>>7>>14>>3>>20>>9 );
    writeVolume( 4,dg,4 );
#endif


} // runTube


void
runSpecial3D_1( data_geo_ptrtype dg )
{
    double lgstruct=0.35101;
    double xL = 0.6-lgstruct;
    double yMin = -0.12,yMax=0.12;
    writePoint(1,dg,  xL, yMin, 0.19 );
    writePoint(2,dg, 0.6, yMin, 0.19 );
    writePoint(3,dg, 0.6, yMin, 0.21 );
    writePoint(4,dg,  xL, yMin, 0.21 );
    writePoint(5,dg,  xL, yMax, 0.19 );
    writePoint(6,dg, 0.6, yMax, 0.19 );
    writePoint(7,dg, 0.6, yMax, 0.21 );
    writePoint(8,dg,  xL, yMax, 0.21 );
    // point sup
    writePoint(9,dg,  0.2     , yMin, 0.2 );// center
    writePoint(10,dg, 0.2-0.05, yMin, 0.2 ); // on circle
    writePoint(11,dg, 0.2     , yMax, 0.2 );// center
    writePoint(12,dg, 0.2-0.05, yMax, 0.2); // on circle

    writeLine(1,dg, 1,2);
    writeLine(2,dg, 2,3);
    writeLine(3,dg, 3,4);
    writeLine(4,dg, 4,1);
    writeLine(5,dg, 5,6);
    writeLine(6,dg, 6,7);
    writeLine(7,dg, 7,8);
    writeLine(8,dg, 8,5);
    writeLine(9,dg, 1,5);
    writeLine(10,dg,2,6);
    writeLine(11,dg,3,7);
    writeLine(12,dg,4,8);
    // line sup
    writeCircle(13,dg, 1,9,10 );
    writeCircle(14,dg, 10,9,4 );
    writeCircle(15,dg, 5,11,12 );
    writeCircle(16,dg, 12,11,8 );
    // line on cylinder
    writeLine(17,dg, 12, 10);

    writeLineLoop(1,dg, Loop()>>1>>2>>3>>4);
    writeLineLoop(2,dg, Loop()>>-8>>-7>>-6>>-5);
    writeLineLoop(3,dg, Loop()>>9>>5>>-10>>-1);
    writeLineLoop(4,dg, Loop()>>10>>6>>-11>>-2);
    writeLineLoop(5,dg, Loop()>>11>>7>>-12>>-3);
    //writeLineLoop(6,dg, Loop()>>9>>-8>>-12>>4);
    writePlaneSurface(1,dg,1);
    writePlaneSurface(2,dg,2);
    writePlaneSurface(3,dg,3);
    writePlaneSurface(4,dg,4);
    writePlaneSurface(5,dg,5);
    //writePlaneSurface(6,dg,6);

    writeLineLoop(6,dg,  Loop()>>17>>14>>12>>-16 );
    writeLineLoop(7,dg,  Loop()>>-15>>-9>>13>>-17);
    writeLineLoop(8,dg,  Loop()>>-13>>-4>>-14);
    writeLineLoop(9,dg, Loop()>>15>>16>>8);
    writeRuledSurface(6,dg,6);
    writeRuledSurface(7,dg,7);
    writePlaneSurface(8,dg,8);
    writePlaneSurface(9,dg,9);
    writePtInSurface(dg,9,8 );
    writePtInSurface(dg,11,9 );

#if 0
    writeSurfaceLoop( 1, dg, Loop()>>7>>10>>8>>9>>6 );
    writeVolume( 1,dg,1 );
    writeSurfaceLoop( 2, dg, Loop()>>5>>4>>3>>1>>2>>6 );
    writeVolume( 2,dg,2 );
#else
    writeSurfaceLoop( 1, dg, Loop()>>6>>9>>7>>8>>5>>4>>3>>1>>2 );
    writeVolume( 1,dg,1 );
#endif

} // runSpecial3D_1


} //namespace GeoTool

} //namespace Feel
