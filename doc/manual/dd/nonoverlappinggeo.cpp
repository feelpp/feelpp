/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Abdoulaye Samake <Abdoulaye.Samake@imag.fr>
       Date: 2011-12-31

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
   \file fixddgeo.cpp
   \author Abdoulaye Samake <Abdoulaye.Samake@imag.fr>
   \date 2011-06-04
 */
#include <feel/feelfilters/gmsh.hpp>

namespace Feel
{
gmsh_ptrtype
nonOverlapGeometryLeft( int RDim, double hsize )
{
    std::ostringstream ostr;
    std::ostringstream nameStr;
    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setOrder( GMSH_ORDER_ONE );
    gmshp->setRecombine( false );
    gmshp->setCharacteristicLength( hsize );
    ostr << gmshp->preamble() << "\n";

    switch ( RDim )
    {
    case 2:
        ostr << "Point(1) = {0,0,0,h};\n"
             << "Point(2) = {1,0,0,h};\n"
             << "Point(3) = {0,1,0,h};\n"
             << "Point(4) = {-1,0,0,h/5};\n"
             << "Point(5) = {0,-1,0,h/5};\n"
             << "Line(1) = {2,3};\n"
             << "Circle(2) = {3,1,4};\n"
             << "Circle(3) = {4,1,5};\n"
             << "Circle(4) = {5,1,2};\n"
             << "Line Loop (5) = {1,2,3,4};\n"
             << "Plane Surface(6) = {5};\n"
             << "Physical Line(1) = {2};\n"
             << "Physical Line(2) = {3};\n"
             << "Physical Line(3) = {1};\n"
             << "Physical Line(4) = {4};\n"
             << "Physical Surface(\"Mat1\") = {6};\n";
        nameStr << "leftgeo2D";
        break;

    case 3:
        ostr << "Point(1) = {-1,-2,-1, h};\n"
             << "Point(2) = {1,-2,-1, h};\n;"
             << "Point(3) = {1,2,-1, h};\n"
             << "Point(4) = {-1,2,-1, h};\n"
             << "Point(5) = {-1,-2,1, h};\n"
             << "Point(6) = {1,-2,1, h};\n"
             << "Point(7) = {1,2,1, h};\n"
             << "Point(8) = {-1,2,1, h};\n"
             << "Point(9) = {1.0,0.0,0.0, h};\n"
             << "Point(10) = {1.0,0.5,0.0, h};\n"
             << "Point(11) = {1.0,0.0,-0.5, h};\n"
             << "Point(12) = {1.0,-0.5,0.0, h};\n"
             << "Point(13) = {1.0,0.0,0.5, h};\n"
             << "Line(1) = {1,2};\n"
             << "Line(2) = {2,3};\n"
             << "Line(3) = {3,4};\n"
             << "Line(4) = {4,1};\n"
             << "Line(5) = {5,6};\n"
             << "Line(6) = {6,7};\n"
             << "Line(7) = {7,8};\n"
             << "Line(8) = {8,5};\n"
             << "Line(9) = {1,5};\n"
             << "Line(10) = {2,6};\n"
             << "Line(11) = {3,7};\n"
             << "Line(12) = {4,8};\n"
             << "Line Loop(1) = {1,2,3,4};\n"
             << "Line Loop(2) = {5,6,7,8};\n"
             << "Line Loop(3) = {1,10,-5,-9};\n"
             << "Line Loop(5) = {11,7,-12,-3};\n"
             << "Line Loop(6) = {9,-8,-12,4};\n"
             << "Circle(13) = {10,9,11};\n"
             << "Circle(14) = {11,9,12};\n"
             << "Circle(15) = {12,9,13};\n"
             << "Circle(16) = {13,9,10};\n"
             << "Line Loop(7) = {13,14,15,16};\n"
             << "Line Loop(4) = {10,6,-11,-2,13,14,15,16};\n"
             << "Plane Surface(1) = {1};\n"
             << "Plane Surface(2) = {2};\n"
             << "Plane Surface(3) = {3};\n"
             << "Plane Surface(4) = {4};\n"
             << "Plane Surface(5) = {5};\n"
             << "Plane Surface(6) = {6};\n"
             << "Plane Surface(7) = {7};\n"
             << "Surface Loop(1) = {1,2,3,4,5,6,7};\n"
             << "Physical Surface(1) = {1}\n"
             << "Physical Surface(2) = {2};\n"
             << "Physical Surface(3) = {3};\n"
             << "Physical Surface(4) = {4};\n"
             << "Physical Surface(5) = {5};\n"
             << "Physical Surface(6) = {6};\n"
             << "Physical Surface(7) = {7};\n"
             << "Volume(1) = {1};\n"
             << "Physical Volume(\"Omega1\") = {1};\n";
        nameStr << "Parallelepiped.geo";
        break;

    default:
        std::ostringstream os;
        os << "invalid dimension: " << RDim;
        throw std::logic_error( os.str() );
    }

    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( ostr.str() );
    return gmshp;
}

gmsh_ptrtype
nonOverlapGeometryRight( int RDim, double hsize )
{
    std::ostringstream ostr;
    std::ostringstream nameStr;
    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setOrder( GMSH_ORDER_ONE );
    gmshp->setRecombine( false );
    gmshp->setCharacteristicLength( hsize );
    ostr << gmshp->preamble() << "\n";

    switch ( RDim )
    {
    case 2:
        ostr << "Point(1) = {1,0,0,h};\n"
             << "Point(2) = {2,0,0,h};\n"
             << "Point(3) = {2,1,0,h/2};\n"
             << "Point(4) = {0,1,0,h/2};\n"
             << "Line(1) = {1,2};\n"
             << "Line(2) = {2,3};\n"
             << "Line(3) = {3,4};\n"
             << "Line(4) = {4,1};\n"
             << "Line Loop (5) = {1,2,3,4};\n"
             << "Plane Surface(6) = {5};\n"
             << "Physical Line(1) = {4};\n"
             << "Physical Line(2) = {1};\n"
             << "Physical Line(3) = {2};\n"
             << "Physical Line(4) = {3};\n"
             << "Physical Surface(\"Mat1\") = {6};\n";
        nameStr << "rightgeo2D";
        break;

    case 3:
        ostr << "Point(1) = {1.,0.,0., h};\n"
             << "Point(2) = {1.,0.5,0., h};\n"
             << "Point(3) = {1.,0.,-0.5, h};\n"
             << "Point(4) = {1.,-0.5,0., h};\n"
             << "Point(5) = {1.,0.,0.5, h};\n"
             << "Point(6) = {6.,0.,0., h};\n"
             << "Point(7) = {6.,0.5,0., h};\n"
             << "Point(8) = {6.,0.,-0.5, h};\n"
             << "Point(9) = {6.,-0.5,0., h};\n"
             << "Point(10) = {6.,0.,0.5, h};\n"
             << "Circle(1) = {2,1,3};\n"
             << "Circle(2) = {3,1,4};\n"
             << "Circle(3) = {4,1,5};\n"
             << "Circle(4) = {5,1,2};\n"
             << "Line Loop(1) = {1,2,3,4};\n"
             << "Circle(5) = {7,6,8};\n"
             << "Circle(6) = {8,6,9};\n"
             << "Circle(7) = {9,6,10};\n"
             << "Circle(8) = {10,6,7};\n"
             << "Line Loop(2) = {5,6,7,8};\n"
             << "Line(9) = {4,9};\n"
             << "Line(10) = {5,10};\n"
             << "Line(11) = {2,7};\n"
             << "Line(12) = {8,3};\n"
             << "Line Loop(3) = {9,-6,12,2};\n"
             << "Line Loop(4) = {9,7,-10,-3};\n"
             << "Line Loop(5) = {10,8,-11,-4};\n"
             << "Line Loop(6) = {11,5,12,-1};\n"
             << "Plane Surface(1) = {1};\n"
             << "Plane Surface(2) = {2};\n"
             << "Point{1} In Surface{1};\n"
             << "Point{6} In Surface{2};\n"
             << "Ruled Surface(3) = {3};\n"
             << "Ruled Surface(4) = {4};\n"
             << "Ruled Surface(5) = {5};\n"
             << "Ruled Surface(6) = {6};\n"
             << "Surface Loop(1) = {5,4,3,2,6,1};\n"
             << "Volume(1) = {1};\n"
             << "Physical Surface(1) = {1};\n"
             << "Physical Surface(2) = {2};\n"
             << "Physical Surface(3) = {3};\n"
             << "Physical Surface(4) = {4};\n"
             << "Physical Surface(5) = {5};\n"
             << "Physical Surface(6) = {6};\n"
             << "Physical Volume(\"Omega2\") = {1};\n";
        nameStr << "Cylinder.geo";
        break;

    default:
        std::ostringstream os;
        os << "invalid dimension: " << RDim;
        throw std::logic_error( os.str() );
    }

    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( ostr.str() );
    return gmshp;

}

}
