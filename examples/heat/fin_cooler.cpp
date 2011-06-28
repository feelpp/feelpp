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

std::pair<std::string, std::string>
makefin( double hsize )
{
    std::ostringstream ostr;
    ostr << "Mesh.MshFileVersion = 1;\n"
         << "h=" << hsize << ";\n"
		 << "Point (1) = {0, 0, 0, h};\n"
		 << "Point (2) = {1, 0, 0, h};\n"
		 << "Point (3) = {1, 1, 0, h};\n"
		 << "Point (4) = {0.333, 1, 0, h};\n"
		 << "Point (5) = {0.333, 5, 0, h};\n"
		 << "Point (6) = {0, 5, 0, h};\n"
		 << "Point (7) = {0, 1, 0, h};\n"
		 << "Line (1) = {1, 2};\n"
		 << "Line (2) = {2, 3};\n"
		 << "Line (3) = {3, 4};\n"
		 << "Line (4) = {4, 5};\n"
		 << "Line (5) = {5, 6};\n"
		 << "Line (6) = {6, 7};\n"
		 << "Line (7) = {7, 1};\n"
		 << "Line (8) = {4, 7};\n"
		 << "Line Loop (9) = {1, 2, 3, 8, 7};\n"
		 << "Plane Surface (9) = {9};\n"
		 << "Line Loop (10) = {4, 5, 6, -8};\n"
		 << "Plane Surface (10) = {10};\n"
	
		// physical entities
		 << "Physical Line (1) = {1, 2, 3, 4, 7};\n"
		 << "Physical Line (2) = {4, 5, 6, 7};\n"
		 << "Physical Surface (1) = {9};\n"
		 << "Physical Surface (2) = {10};\n";


    std::ostringstream nameStr;
    nameStr.precision( 3 );
    nameStr << "fin";

    return std::make_pair( nameStr.str(), ostr.str() );
}

