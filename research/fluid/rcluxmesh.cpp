/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-03-19

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file rcluxmesh.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-03-19
 */
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>



std::pair<std::string, std::string>
createMesh( double meshSize )
{
    std::ostringstream ostr;

    ostr << "h=" << meshSize << ";\n"
         << "hfine=" << "h/2;\n"
         << "Point(1) = {0,0,0,hfine};\n"
         << "Point(2) = {0.002,0,0,hfine};\n"
         << "Point(3) = {0.002,0.02,0,hfine};\n"
         << "Point(4) = {0.01,0.02,0,h};\n"
         << "Point(5) = {0.01,0.138,0,h};\n"
         << "Point(6) = {0.002,0.138,0,hfine};\n"
         << "Point(7) = {0.002,0.188,0,h};\n"
         << "Point(8) = {0,0.188,0,hfine};\n"
         << "Point(9) = {0,0.0138,0,hfine};\n"
         << "Point(10) = {0,0.0138,0,hfine};\n"
         << "Point(11) = {0,0.138,0,hfine};\n"
         << "Point(12) = {0,0.02,0,hfine};\n"
         << "\n"
         << "Line(1) = {1,2};\n"
         << "Line(2) = {2,3};\n"
         << "Line(3) = {3,4};\n"
         << "Line(4) = {4,5};\n"
         << "Line(5) = {5,6};\n"
         << "Line(6) = {6,7};\n"
         << "Line(7) = {7,8};\n"
         << "Line(8) = {8,11};\n"
         << "Line(9) = {11,12};\n"
         << "Line(10) = {12,1};\n"
         << "Line(11) = {3,12};\n"
         << "Line(12) = {6,11};\n"
         << "Line Loop(13) = {1,2,11,10};\n"
         << "Plane Surface(14) = {13};\n"
         << "Line Loop(15) = {3,4,5,12,9,-11};\n"
         << "Plane Surface(16) = {15};\n"
         << "Line Loop(17) = {8,-12,6,7};\n"
         << "Plane Surface(18) = {17};\n"
         << "\n"
         << "\n"
         << "Physical Line(\"inflow\") = {1};\n"
         << "Physical Line(\"wall\") = {2,3,4,5,6};\n"
         << "Physical Line(\"outflow\") = {7};\n"
         << "Physical Line(\"symmetry\") = {8,9,10};\n"
         << "\n"
         << "Physical Surface(\"dom1\") = {14};\n"
         << "Physical Surface(\"dom2\") = {16};\n"
         << "Physical Surface(\"dom3\") = {18};\n";



    std::ostringstream nameStr;
    nameStr.precision( 3 );
    nameStr << "rclux";

    return std::make_pair( nameStr.str(), ostr.str() );
}


