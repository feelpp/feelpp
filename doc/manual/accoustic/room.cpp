/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-07-06

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
   \file room.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-07-06
 */
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>
#include <feel/feelfilters/gmsh.hpp>

namespace Feel
{
gmsh_ptrtype
createRoom( int Dim, double meshSize )
{
    std::ostringstream ostr;
    std::ostringstream nameStr;

    //    std::string fname;//
    switch ( Dim )
    {

        // On rajoute le cas ou la dimension = 1

    case 1:
        ostr << "hTemp=" << meshSize << ";\n"
             << "Point (1) = {0, 0, 0, hTemp};\n"
             << "Point (2) = {1, 0, 0, hTemp};\n"
             << "Line (1) = {1, 2};\n"
             << "Physical Point (1) = {1};\n"
             << "Physical Point (2) = {2};\n"
             << "Physical Line (10) = {1};\n";

        nameStr << "room." << meshSize;
        break;


    case 2:
        ostr << "hTemp=" << meshSize << ";\n"
             << "Point (1) = {0, 0, 0, hTemp};\n"
             << "Point (2) = {1, 0, 0, hTemp};\n"
             << "Point (3) = {5, 1, 0, hTemp};\n"
             << "Point (4) = {5, 3, 0, hTemp};\n"
             << "Point (5) = {3, 3, 0, hTemp};\n"
             << "Point (6) = {1, 1, 0, hTemp};\n"
             << "Point (7) = {0, 1, 0, hTemp};\n"
             << "Line (1) = {1, 2};\n"
             << "Line (2) = {2, 3};\n"
             << "Line (3) = {3, 4};\n"
             << "Line (4) = {4, 5};\n"
             << "Line (5) = {5, 6};\n"
             << "Line (6) = {6, 7};\n"
             << "Line (7) = {7, 1};\n"
             << "Line Loop (9) = {4, 5, 6, 7, 1, 2, 3};\n"
             << "Plane Surface (9) = {9};\n"
             << "Physical Line (1) = {1, 2, 3, 4, 5, 6};\n"
             << "Physical Line (2) = {7};\n"
             << "Physical Surface (13) = {9};\n";

        nameStr << "room." << meshSize;
        break;

    case 3:
        ostr << "hTemp=" << meshSize << ";\n"
             << "Point (1) = {0, 0, 0, hTemp};\n"
             << "Point (2) = {1, 0, 0, hTemp};\n"
             << "Point (3) = {5, 1, 0, hTemp};\n"
             << "Point (4) = {5, 3, 0, hTemp};\n"
             << "Point (5) = {3, 3, 0, hTemp};\n"
             << "Point (6) = {1, 1, 0, hTemp};\n"
             << "Point (7) = {0, 1, 0, hTemp};\n"
             << "Line (1) = {1, 2};\n"
             << "Line (2) = {2, 3};\n"
             << "Line (3) = {3, 4};\n"
             << "Line (4) = {4, 5};\n"
             << "Line (5) = {5, 6};\n"
             << "Line (6) = {6, 7};\n"
             << "Line (7) = {7, 1};\n"
             << "Line Loop (9) = {4, 5, 6, 7, 1, 2, 3};\n"
             << "Plane Surface (9) = {9};\n"
             << "Extrude {0,0,1} {\n"
             << "  Surface{9};\n"
             << "}\n"
             << "Surface Loop(47) = {33,9,21,25,29,46,37,41,45};\n"
             << "Volume(1) = {47};\n"
             << "Physical Surface(51) = {33};\n"
             << "Physical Surface(52) = {37,41,46,9,25,29,21,45};\n"
             << "Physical Volume(1) = {1};\n";
        nameStr << "room." << meshSize;
        break;

    default:
        std::ostringstream os;
        os << "invalid dimension: " << Dim;
        throw std::logic_error( os.str() );
    }

    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( ostr.str() );
    return gmshp;
}
} // Feel
