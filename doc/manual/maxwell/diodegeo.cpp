/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-06-15

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file diode-mesh.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-06-15
 */
#include <iostream>
#include <sstream>
#include <string>

#include <feel/feelfilters/gmsh.hpp>

using namespace Feel;

gmsh_ptrtype
diodegeo( double h, int Order, std::string const& convex )
{
    std::ostringstream ostr;
    std::ostringstream nameStr;
    gmsh_ptrtype gmshp( new Gmsh( 2, Order ) );
    gmshp->setCharacteristicLength( h );
    ostr << gmshp->preamble() << "\n"
         << "Point(1) = {0, 0, 0, h};\n"
         << "Point(2) = {0, 1, 0, h};\n"
         << "Point(3) = {0, 2, 0, h};\n"
         << "Point(4) = {2, 0, 0, h};\n"
         << "Point(5) = {2, 1, 0, h};\n"
         << "Point(6) = {2, 2, 0, h};\n"
         << "Point(7) = {3, 0, 0, h};\n"
         << "Point(8) = {4, 0, 0, h};\n"
         << "Point(9) = {4, 2, 0, h};\n"
         << "Point(10) = {2.707106781, 0.707106781, 0, h};\n"
         << "\n"
         << "Line(1) = {2, 3};\n"
         << "Line(2) = {3, 6};\n"
         << "Line(3) = {9, 6};\n"
         << "Line(4) = {9, 8};\n"
         << "Line(5) = {8, 7};\n"
         << "Line(6) = {10, 9};\n"
         << "Line(7) = {6, 5};\n"
         << "Line(8) = {5, 2};\n"
         << "Circle(9) = {5, 4, 10};\n"
         << "Circle(10) = {7, 4, 10};\n"
         << "//Physical Line(8) = {1,2,3,4,5,9,10,8};\n"
         << "//Physical Line(8) = {1};\n";

    if ( convex == "simplex" )
    {
        ostr << "Line Loop(11) = {1, 2, -3, 4, 5, 10, -9, 8};\n"
             << "Plane Surface(12) = {11};\n"
             << "\n"
             << "Physical Line(\"Dirichlet\") = {1};\n"
             << "Physical Line(\"Metal\") = {2,3,4,5,8,9,10};\n"
             << "Physical Surface(\"Omega\") = {12};\n";
    }

    if ( convex == "hypercube" )
    {
        ostr << "Line Loop(1) = {1, 8, 7, 2};\n"
             << "Plane Surface(1) = {1};\n"
             << "Line Loop(2) = {7, 9, 6, 3};\n"
             << "Plane Surface(2) = {2};\n"
             << "Line Loop(3) = {4, 6, 10, 5};\n"
             << "Plane Surface(3) = {3};\n"
             << "\n"
             << "Physical Line(\"Dirichlet\") = {1};\n"
             << "Physical Line(\"Metal\") = {2,3,4,5,8,9,10};\n"
             << "Physical Surface(9) = {-1};\n"
             << "Physical Surface(10) = {2};\n"
             << "Physical Surface(11) = {-3};\n"

             << "Recombine Surface {1};\n"
             << "Transfinite Surface {1};\n"
             << "Recombine Surface {2};\n"
             << "Transfinite Surface {2};\n"
             << "Recombine Surface {3};\n"
             << "Transfinite Surface {3};\n"
             << "Transfinite Line {1, 7} = 17 Using Progression 1;\n"
             << "Transfinite Line {7, 6} = 17 Using Progression 1;\n"
             << "Transfinite Line {6, 5} = 17 Using Progression 1;\n"
             << "Transfinite Line {2, 8} = 17 Using Progression 1;\n"
             << "Transfinite Line {3, 9} = 17 Using Progression 1;\n"
             << "Transfinite Line {4, 10}= 17 Using Progression 1;\n";
        //Mesh.ElementOrder=2 ;
        //Mesh.SecondOrderIncomplete = 1 ;
    }


    nameStr << "diode-" << convex << "-" << Order;


    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( ostr.str() );
    return gmshp;
}



