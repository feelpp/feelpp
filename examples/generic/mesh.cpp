/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2006-10-12
 */
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>

std::pair<std::string,std::string>
createRing( int Dim, double meshSize )
{
    std::ostringstream ostr;
    std::ostringstream nameStr;
//    std::string fname;//
    switch( Dim ) {
    case 2:
        ostr << "h=" << meshSize << ";\n"
             << "Point(1) = {0.1,0,0,h/2};\n"
             << "Point(2) = {1,0,0,h};\n"
             << "Point(3) = {0,1,0,h};\n"
             << "Point(4) = {0,0.1,0,h/2};\n"
             << "Point(5) = {0,0,0,h/2};\n"
             << "Line(1) = {1,2};\n"
             << "Circle(2) = {2,5,3};\n"
             << "Line(3) = {3,4};\n"
             << "Circle(4) = {4,5,1};\n"
             << "Line Loop(5) = {1,2,3,4};\n"
             << "Plane Surface(6) = {5};\n"
             << "Physical Line(10) = {1};\n"
             << "Physical Line(20) = {2};\n"
             << "Physical Line(30) = {3};\n"
             << "Physical Line(40) = {4};\n"
//             << "Physical Line(20) = {1,2,4};\n"
             << "Physical Surface(7) = {6};\n";
        nameStr << "ring." << meshSize;
//        fname = __gmsh.generateSquare( "advectiondg2d", meshSize );//
        break;
// To be added for 3D something like:
/*    case 3:
        ostr << "h=" << meshSize << ";\n"
             << "Point(1) = {-1,-1,-1,h};\n"
             << "Point(2) = {-1, 1,-1,h};\n"
             << "Point(3) = { 1, 1,-1,h};\n"
             << "Point(4) = { 1,-1,-1,h};\n"
             << "Line(1) = {1,4};\n"
             << "Line(2) = {4,3};\n"
             << "Line(3) = {3,2};\n"
             << "Line(4) = {2,1};\n"
             << "Line Loop(5) = {3,4,1,2};\n"
             << "Plane Surface(6) = {5};\n"
             << "Extrude Surface {6, {0,0,2}};\n"
             << "Physical Surface(10) = {15,23,6,28};\n"
             << "Physical Surface(20) = {19,27};\n"
             << "Surface Loop(31) = {28,15,-6,19,23,27};\n"
             << "Volume(1) = {31};\n"
             << "Physical Volume(2) = {1};\n";
        nameStr << "cube." << meshSize;
        break;*/
    default:
        std::ostringstream os;
        os << "invalid dimension: " << Dim;
        throw std::logic_error( os.str() );
    }
    return std::make_pair( nameStr.str(), ostr.str() );
}
