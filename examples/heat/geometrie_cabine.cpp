/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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



    ostr << "hTemp=" << meshSize << ";\n"
         << "Point(1) = {-0.8,-0.4,0,hTemp};\n"
         << "Point(2) = {0.8,-0.4,0,hTemp};\n"
         << "Point(3) = {0.8,0.4,0,hTemp};\n"
         << "Point(4) = {-0.8,0.4,0,hTemp};\n"
         << "Line(1) = {1,2};\n"
         << "Line(2) = {2,3};\n"
         << "Line(3) = {3,4};\n"
         << "Line(4) = {4,1};\n"
         << "Line Loop(1) = {2,3,4,1};\n"
         << "Plane Surface(3) = {1};\n"
         << "Physical Line(11) = {4};\n"
         << "Physical Line(12) = {2};\n"
         << "Physical Line(13) = {1,3};\n"
         << "Physical Surface(14) = {3};\n";
    //diriclet=11
    //neumann=12
    //autre=13
    //surface=14

    nameStr << "room." << meshSize;
    //        break;

    //    default:
    //        std::ostringstream os;
    //        os << "invalid dimension: " << Dim;
    //        throw std::logic_error( os.str() );


    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( ostr.str() );
    return gmshp;
}
} // Feel
