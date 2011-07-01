/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-06-11

  Copyright (C) 2007,2010 Université Joseph Fourier (Grenoble I)

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
   \file fin_cooler.cpp
   \author Baptiste Morin <baptistemorin@gmail.com> <baptiste.morin@e.ujf-grenoble.fr>
   \date 2011-06-28
 */
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feeldiscr/mesh.hpp>

using namespace Feel;


gmsh_ptrtype
makefin( double hsize , double deep)
{
    std::ostringstream ostr;
    if (!deep) { // 2D Mesh
        //typedef Simplex<2> convex_type;
        ostr << "Mesh.MshFileVersion = 1;\n"
             << "h=" << hsize << ";\n"
             << "Point (1) = {0, 0, 0, h};\n"
             << "Point (2) = {3.5, 0, 0, h};\n"
             << "Point (3) = {3.5, 3, 0, h};\n"
             << "Point (4) = {1, 3, 0, h};\n"
             << "Point (5) = {0, 3, 0, h};\n"
             << "Point (6) = {1, 17, 0, h};\n"
             << "Point (7) = {0, 17, 0, h};\n"
             << "Line (1) = {1, 2};\n"
             << "Line (2) = {2, 3};\n"
             << "Line (3) = {3, 4};\n"
             << "Line (4) = {4, 5};\n"
             << "Line (5) = {4, 6};\n"
             << "Line (6) = {6, 7};\n"
             << "Line (7) = {7, 5};\n"
             << "Line (8) = {5, 1};\n"
             << "Line Loop (9) = {1, 2, 3, 4, 8};\n"
             << "Line Loop (10) = {5, 6, 7, -4};\n"
             << "Plane Surface (9) = {9};\n"
             << "Plane Surface (10) = {10};\n"
            // physical entities
             << "Physical Surface (1) = {9};\n"
             << "Physical Surface (2) = {10};\n";
    }

    else { //3D Mesh
        //typedef Simplex<3> convex_type;
        ostr << "Mesh.MshFileVersion = 1;\n"
             << "h=" << hsize << ";\n"
             << "d=" << deep << ";\n"
             << "Point (1) = {0, 0, 0, h};\n"
             << "Point (2) = {3.5, 0, 0, h};\n"
             << "Point (3) = {3.5, 3, 0, h};\n"
             << "Point (4) = {1, 3, 0, h};\n"
             << "Point (5) = {0, 3, 0, h};\n"
             << "Point (6) = {1, 17, 0, h};\n"
             << "Point (7) = {0, 17, 0, h};\n"
             << "Point (8) = {0, 0, d, h};\n"
             << "Point (9) = {3.5, 0, d, h};\n"
             << "Point (10) = {3.5, 3, d, h};\n"
             << "Point (11) = {1, 3, d, h};\n"
             << "Point (12) = {0, 3, d, h};\n"
             << "Point (13) = {1, 17, d, h};\n"
             << "Point (14) = {0, 17, d, h};\n"
             << "Line (1) = {1, 2};\n"
             << "Line (2) = {2, 3};\n"
             << "Line (3) = {3, 4};\n"
             << "Line (4) = {4, 5};\n"
             << "Line (5) = {4, 6};\n"
             << "Line (6) = {6, 7};\n"
             << "Line (7) = {7, 5};\n"
             << "Line (8) = {5, 1};\n"
             << "Line (9) = {2, 9};\n"
             << "Line (10) = {1, 8};\n"
             << "Line (11) = {8, 9};\n"
             << "Line (12) = {12, 8};\n"
             << "Line (13) = {9, 10};\n"
             << "Line (14) = {10, 11};\n"
             << "Line (15) = {5, 12};\n"
             << "Line (16) = {3, 10};\n"
             << "Line (17) = {4, 11};\n"
             << "Line (18) = {11, 13};\n"
             << "Line (19) = {14, 12};\n"
             << "Line (20) = {11, 12};\n"
             << "Line (21) = {13, 14};\n"
             << "Line (22) = {6, 13};\n"
             << "Line (23) = {7, 14};\n"
             << "Line Loop (24) = {1, 9, -11, -10};\n"
             << "Line Loop (25) = {-3, 16, 14, 20, -15, -4 };\n"
             << "Line Loop (26) = {1, 2, 3, 4, 8};\n"
             << "Line Loop (27) = {5, 6, 7, -4};\n"
             << "Line Loop (28) = {11, 13, 14, 20, 12};\n"
             << "Line Loop (29) = {18, 21, 19, -20};\n"
             << "Line Loop (30) = {-6, 22, 21, -23};\n"
             << "Line Loop (31) = {5, 22, -18, -17};\n"
             << "Line Loop (32) = {15, -19, -23, 7};\n"
             << "Line Loop (33) = {-6, 22, 21, -23};\n"
             << "Line Loop (34) = {17, 20, -15, -4};\n"
             << "Line Loop (35) = {9, 13, -16, -2};\n"
             << "Line Loop (36) = {10, -12, -15, 8};\n"
             << "Line Loop (37) = {-3, 16, 14, -17};\n"
             << "Plane Surface (24) = {24};\n"
             << "Plane Surface (25) = {25};\n"
             << "Plane Surface (26) = {26};\n"
             << "Plane Surface (27) = {27};\n"
             << "Plane Surface (28) = {28};\n"
             << "Plane Surface (29) = {29};\n"
             << "Plane Surface (30) = {30};\n"
             << "Plane Surface (31) = {31};\n"
             << "Plane Surface (32) = {32};\n"
             << "Plane Surface (33) = {33};\n"
             << "Plane Surface (34) = {34};\n"
             << "Plane Surface (35) = {35};\n"
             << "Plane Surface (36) = {36};\n"
             << "Plane Surface (37) = {37};\n"
             << "Surface Loop (1) = {32, 31, 27, 29, 33, 34};\n"
             << "Surface Loop (2) = {35, 36, 34, 26, 24, 28, 37};\n"
             << "Volume (1) = {1};\n"
             << "Volume (2) = {2};\n";
    }

    std::cout<<"start "<<std::endl;
    std::ostringstream nameStr;
    nameStr.precision( 3 );
    nameStr << "fin_sink";
    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( ostr.str() );
    return gmshp;

}

