/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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


/**
 * makefin creates the .geo file and returns a gmsh pointer type.
 * Be carefull, the problem considered here is dimensional
 *
 * \arg hsize : size of the mesh
 * \arg width : the width of the fin
 * \arg deep : depth of the mesh (in meters)
 * \arg L : dimensional length of the sink (in meters)
 */
gmsh_ptrtype
makefin( double hsize, double width, double deep, double L )
{


    std::ostringstream nameStr;
    nameStr.precision( 3 );
    nameStr << "fin_sink";
    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setPrefix( nameStr.str() );
    std::ostringstream ostr;
    ostr << gmshp->preamble();
    if ( !deep ) // 2D Mesh
    {
        ostr << "Point (1) = {0, 0, 0, " << hsize << "};\n"
             << "Point (2) = {0.001, 0, 0, " << hsize << "};\n"
             << "Point (3) = {0.001, 0.001, 0, " << hsize << "};\n"
             << "Point (4) = {"<< width <<", 0.001, 0, " << hsize << "};\n"
             << "Point (5) = {0, 0.001, 0, " << hsize << "};\n"
             << "Point (6) = {"<< width <<", "<<0.0005+L<<", 0, " << hsize << "};\n"
             << "Point (7) = {0, "<<0.0005+L<<", 0, " << hsize << "};\n"
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
             << "Physical Line (\"gamma1\") = {5};\n"
             << "Physical Line (\"gamma2\") = {6};\n"
             << "Physical Line (\"gamma3\") = {4};\n"
             << "Physical Line (\"gamma4\") = {1};\n"
             << "Physical Line (\"gamma5\") = {3};\n"
             << "Physical Line (\"gamma6\") = {7};\n"
             << "Physical Line (\"gamma7\") = {2};\n"
             << "Physical Line (\"gamma8\") = {8};\n"
             << "Physical Surface (\"spreader_mesh\") = {9};\n"
             << "Physical Surface (\"fin_mesh\") = {10};\n";
    }

    else   //3D Mesh
    {
        ostr << "Point (1) = {0, 0, 0, " << hsize << "};\n"
             << "Point (2) = {0.001, 0, 0, " << hsize << "};\n"
             << "Point (3) = {0.001, 0.001, 0, " << hsize << "};\n"
             << "Point (4) = {"<< width <<", 0.001, 0, " << hsize << "};\n"
             << "Point (5) = {0, 0.001, 0, " << hsize << "};\n"
             << "Point (6) = {"<< width <<", "<<0.0005+L<<", 0, " << hsize << "};\n"
             << "Point (7) = {0, "<<0.0005+L<<", 0, " << hsize << "};\n"
             << "Point (8) = {0, 0, "<< deep <<", " << hsize << " };\n"
             << "Point (9) = {0.001, 0, "<< deep <<", " << hsize << " };\n"
             << "Point (10) = {0.001, 0.001, "<< deep <<", " << hsize << " };\n"
             << "Point (11) = {"<< width <<", 0.001, "<< deep <<", " << hsize << " };\n"
             << "Point (12) = {0, 0.001, "<< deep <<", " << hsize << " };\n"
             << "Point (13) = {"<< width <<", "<<0.0005+L<<", "<< deep <<", " << hsize << " };\n"
             << "Point (14) = {0, "<<0.0005+L<<", "<< deep <<", " << hsize << " };\n"
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
             // name the physical entities
             << "Physical Surface (\"gamma4\") = {24};\n"
             << "Physical Surface (25) = {25};\n"
             << "Physical Surface (26) = {26};\n"
             << "Physical Surface (27) = {27};\n"
             << "Physical Surface (28) = {28};\n"
             << "Physical Surface (29) = {29};\n"
             << "Physical Surface (30) = {30};\n"
             << "Physical Surface (\"gamma1\") = {31};\n"
             << "Physical Surface (32) = {32};\n"
             << "Physical Surface (33) = {33};\n"
             << "Physical Surface (34) = {34};\n"
             << "Physical Surface (35) = {35};\n"
             << "Physical Surface (36) = {36};\n"
             << "Physical Surface (37) = {37};\n"
             //volumes
             << "Surface Loop (1) = {32, 31, 27, 29, 33, 34};\n"
             << "Surface Loop (2) = {35, 36, 34, 26, 24, 28, 37};\n"
             << "Volume (3) = {1};\n"
             << "Volume (4) = {2};\n"
             << "Physical Volume (\"fin_mesh\") = {3};\n"
             << "Physical Volume (\"spreader_mesh\") = {4};\n";


    }

    gmshp->setDescription( ostr.str() );
    return gmshp;
}

