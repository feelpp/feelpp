/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-10-12

  Copyright (C) 2006 University Joseph Fourier (UJF)

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
   \file mesh.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-10-12
 */
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>

std::pair<std::string,std::string>
createGeometry( int Dim, double meshSize, double a, double b )
{
    std::ostringstream ostr;
    std::ostringstream nameStr;

    switch ( Dim )
    {
    case 1:
        ostr << "Mesh.MshFileVersion = 1;\n"
             << "h=" << meshSize << ";\n"
             << "a=" << a << ";\n"
             << "Point(1) = {-a,0,0,h};\n"
             << "Point(2) = {a,0,0,h};\n"
             << "Line(1) = {1,2};\n"
             << "Physical Point(21) = {1,2};\n"
             << "Physical Line(23) = {1};\n";
        nameStr << "line." << meshSize;
        break;

    case 2:
#if 0
        ostr << "Mesh.MshFileVersion = 1;\n"
             << "h=" << meshSize << ";\n"
             << "a=" << a << ";\n"
             << "b=" << b << ";\n"
             << "Point(1) = {-a,-b,0,h};\n"
             << "Point(2) = {a,-b,0,h};\n"
             << "Point(3) = {a,b,0,h};\n"
             << "Point(4) = {-a,b,0,h};\n"
             << "Line(1) = {1,2};\n"
             << "Line(2) = {2,3};\n"
             << "Line(3) = {3,4};\n"
             << "Line(4) = {4,1};\n"
             << "Line Loop(5) = {1,2,3,4};\n"
             << "Plane Surface(6) = {5};\n"
             << "Physical Line(10) = {1};\n"
             << "Physical Line(20) = {2,3,4};\n"
             << "Physical Surface(7) = {6};\n";
        nameStr << "triangle." << meshSize;
#endif

        ostr << "Mesh.MshFileVersion = 1;\n"
             << "h=" << meshSize << ";\n"
             << "b = " << a << ";\n"
             << "a = " << b << ";\n"
             << "Point(1) = {-a,-b,0,h};\n"
             << "Point(2) = {a,-b,0,h};\n"
             << "Point(3) = {a,b,0,h};\n"
             << "Point(4) = {-a,b,0,h};\n"
             << "Point(5) = {0,b,0,h};\n"
             << "Point(6) = {0,-b,0,h};\n"
             << "Point(7) = {-a,0,0,h};\n"
             << "Point(8) = {a,0,0,h};\n"
             << "Point(9) = {0,0,0,h};\n"
             << "\n" // circle
             << "Point(10) = {-a/2,-b/2,0,h};\n"
             << "Point(11) = {-a/2-a/4,-b/2,0,h};\n"
             << "Point(12) = {-a/2+a/4,-b/2,0,h};\n"
             << "Point(14) = {-a/2,-b/2+b/4,0,h};\n"
             << "Point(15) = {-a/2,-b/2-b/4,0,h};\n"
             << "\n"
             << "Line(1) = {3,5};\n"
             << "Line(2) = {5,4};\n"
             << "Line(3) = {4,7};\n"
             << "Line(4) = {7,1};\n"
             << "Line(5) = {1,6};\n"
             << "Line(6) = {6,9};\n"
             << "Line(7) = {9,7};\n"
             << "Line(8) = {5,9};\n"
             << "Line(9) = {9,8};\n"
             << "Line(10) = {8,3};\n"
             << "Line(11) = {8,2};\n"
             << "Line(12) = {2,6};\n"
             << "\n"
#if 1
             << "Circle(27) = {15,10,11};\n"
             << "Circle(28) = {11,10,14};\n"
             << "Circle(29) = {14,10,12};\n"
             << "Circle(30) = {12,10,15};\n"
             << "Line Loop(31) = {28,29,30,27};\n"
#endif
             << "\n"
             << "Line Loop(13) = {2,3,-7,-8};\n"
             << "Plane Surface(14) = {13};\n"
             << "Line Loop(15) = {1,8,9,10};\n"
             << "Plane Surface(16) = {15};\n"
             << "Line Loop(17) = {7,4,5,6};\n"
#if 1
             << "Plane Surface(18) = {17,31};\n"
#else
             << "Plane Surface(18) = {17};\n"
#endif
             << "Line Loop(19) = {9,11,12,6};\n"
             << "Plane Surface(20) = {19};\n"
             << "\n"
             << "Physical Line(21) = {4,3,2,1,10,11};\n"
             << "Physical Line(22) = {5};\n"
             << "Physical Line(23) = {12};\n"
             << "Physical Surface(23) = {18};\n"
             << "Physical Surface(24) = {20};\n"
             << "Physical Surface(25) = {14};\n"
             << "Physical Surface(26) = {16};\n"
             << "Physical Line(33) = {28,29,30,27};\n";

        nameStr << "triangle." << meshSize;

        break;

    case 3:
        ostr << "Mesh.MshFileVersion = 1;\n"
             << "h=" << meshSize << ";\n"
             << "a=" << a << ";\n"
             << "b=" << b << ";\n"
             << "Point(1) = {-a,-b,-1,h};\n"
             << "Point(2) = {-a, b,-1,h};\n"
             << "Point(3) = { a, b,-1,h};\n"
             << "Point(4) = { a,-b,-1,h};\n"
             << "Line(1) = {1,4};\n"
             << "Line(2) = {4,3};\n"
             << "Line(3) = {3,2};\n"
             << "Line(4) = {2,1};\n"
             << "Line Loop(5) = {3,4,1,2};\n"
             << "Plane Surface(6) = {5};\n"
             << "Extrude Surface {6, {0,0,2}};\n"
             << "Physical Surface(22) = {6};\n"
             << "Physical Surface(21) = {15,19,23,27,28};\n"
             << "Surface Loop(31) = {28,15,-6,19,23,27};\n"
             << "Volume(1) = {31};\n"
             << "Physical Volume(23) = {1};\n";
        nameStr << "cube." << meshSize;
        break;

    default:
        std::ostringstream os;
        os << "invalid dimension: " << Dim;
        throw std::logic_error( os.str() );
    }

    return std::make_pair( nameStr.str(), ostr.str() );
}
