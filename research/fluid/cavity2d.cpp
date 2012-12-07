/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-06-20

  Copyright (C) 2006 EPFL
  Copyright (C) 2007,2008 Université Joseph Fourier (Grenoble 1)

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
   \file cavity2d.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-06-20
 */
#include <examples/fluid/cavity.hpp>


std::string
createCavity( double h )
{
    std::ostringstream ostr;
    ostr << "Mesh.MshFileVersion = 2;\n"
         << "a=0;\n"
         << "b=1;\n"
         << "c=0;\n"
         << "d=1;\n"
         << "h=" << h << ";\n"
         << "Point(1) = {a,c,0.0,h};\n"
         << "Point(2) = {b,c,0.0,h};\n"
         << "Point(3) = {b,d,0.0,h};\n"
         << "Point(4) = {a,d,0.0,h};\n"
         << "Line(1) = {1,4};\n"
         << "Line(2) = {4,3};\n"
         << "Line(3) = {3,2};\n"
         << "Line(4) = {2,1};\n"
         << "Line Loop(5) = {3,4,1,2};\n"
         << "Plane Surface(6) = {5};\n"
         << "Physical Line(1) = {2};\n"
         << "Physical Line(2) = {1,3,4};\n"
         << "Physical Surface(6) = {6};\n";
    return ostr.str();
}

int
main( int argc, char** argv )
{
    /* assertions handling */
    Feel::Cavity<2> cavity( argc, argv, makeAbout(), makeOptions() );
    cavity.run();
}




