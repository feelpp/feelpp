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
   \file fin.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-06-11
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
makefin( double hsize )
{
    std::ostringstream ostr;
    ostr << "h=" << hsize << ";\n"
         << "Point (1) = {0, 0, 0, h};\n"
         << "Point (2) = {1, 0, 0, h};\n"
         << "Point (3) = {1, 1, 0, h};\n"
         << "Point (4) = {0, 1, 0, h};\n"
         << "Point (5) = {0.5, 1, 0, h};\n"
         << "Point (6) = {0, 0.5, 0, h};\n"
         << "Point (7) = {1, 0.5, 0, h};\n"
         << "Point (8) = {0.5, 0, 0, h};\n"
         << "Point (9) = {0.5, 0.5, 0, h};\n"
         << "Line (1) = {1, 8};\n"
         << "Line (2) = {8, 2};\n"
         << "Line (3) = {2, 7};\n"
         << "Line (4) = {7, 3};\n"
         << "Line (5) = {3, 5};\n"
         << "Line (6) = {5, 4};\n"
         << "Line (7) = {4, 6};\n"
         << "Line (8) = {6, 1};\n"
         << "Line (9) = {8, 9};\n"
         << "Line (10) = {9, 5};\n"
         << "Line (11) = {6, 9};\n"
         << "Line (12) = {9, 7};\n"
         << "Line Loop (14) = {7, 11, 10, 6};\n"
         << "Plane Surface (14) = {14};\n"
         << "Line Loop (16) = {8, 1, 9, -11};\n"
         << "Plane Surface (16) = {16};\n"
         << "Line Loop (18) = {2, 3, -12, -9};\n"
         << "Plane Surface (18) = {18};\n"
         << "Line Loop (20) = {12, 4, 5, -10};\n"
         << "Plane Surface (20) = {20};\n"

        // physical entities
         << "Physical Line (1) = {1, 2};\n"
         << "Physical Line (2) = {3, 4, 5, 6, 7, 8};\n"
         << "Physical Surface (1) = {16};\n"
         << "Physical Surface (2) = {18};\n"
         << "Physical Surface (3) = {14};\n"
         << "Physical Surface (4) = {20};\n";


    std::ostringstream nameStr;
    nameStr.precision( 3 );
    nameStr << "fin";

    gmsh_ptrtype gmshp( new Gmsh);
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( ostr.str() );
    return gmshp;
}

