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
ddmethodGeometryLeft( int RDim, double hsize )
{
    std::ostringstream ostr;
    std::ostringstream nameStr;
    gmsh_ptrtype gmshp( new Gmsh( RDim ) );
    gmshp->setOrder( GMSH_ORDER_ONE );
    gmshp->setRecombine( false );
    gmshp->setCharacteristicLength( hsize );
    // ostr << gmshp->preamble() << "\n";


    switch ( RDim )
    {
    case 2:

        ostr<< "h =" << hsize <<";\n"
            << "Point(1) = {0,0,0,h};\n"
            << "Point(2) = {1,0,0,h};\n"
            << "Point(3) = {0,1,0,h};\n"
            << "Point(4) = {-1,0,0,h/5};\n"
            << "Point(5) = {0,-1,0,h/5};\n"
            << "Circle(1) = {2,1,3};\n"
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

        ostr << "Point(1) = {a,0,0,h};\n"
             << "Point(2) = {a,1,0,h};\n"
             << "Point(3) = {a,1,1,h};\n"
             << "Point(4) = {a,0,1,h};\n"
             << "Line(1) = {1,2};\n"
             << "Line(2) = {2,3};\n"
             << "Line(3) = {3,4};\n"
             << "Line(4) = {4,1};\n"
             << "Line Loop(5) = {1,2,3,4};\n"
             << "Plane Surface(6) = {5};\n"
             << "Physical Surface(\"Mat1\") = {6};\n";
        nameStr << "interface3D";
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
ddmethodGeometryRight( int RDim, double hsize )
{
    std::ostringstream ostr;
    std::ostringstream nameStr;
    gmsh_ptrtype gmshp( new Gmsh( RDim ) );
    gmshp->setOrder( GMSH_ORDER_ONE );
    gmshp->setRecombine( false );
    gmshp->setCharacteristicLength( hsize );
    // ostr << gmshp->preamble() << "\n";


    switch ( RDim )
    {
    case 2:

        ostr << "h =" << hsize <<";\n"
             << "Point(1) = {1,0,0,h};\n"
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

        ostr << "Point(1) = {a,0,0,h};\n"
             << "Point(2) = {a,1,0,h};\n"
             << "Point(3) = {a,1,1,h};\n"
             << "Point(4) = {a,0,1,h};\n"
             << "Line(1) = {1,2};\n"
             << "Line(2) = {2,3};\n"
             << "Line(3) = {3,4};\n"
             << "Line(4) = {4,1};\n"
             << "Line Loop(5) = {1,2,3,4};\n"
             << "Plane Surface(6) = {5};\n"
             << "Physical Surface(\"Mat1\") = {6};\n";
        nameStr << "interface3D";
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
