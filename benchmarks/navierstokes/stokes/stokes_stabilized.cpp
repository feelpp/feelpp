/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-06-20

  Copyright (C) 2006 EPFL

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
   \file stokes_stabilized.cpp
   \author Goncalo Pena <goncalo.pena@epfl.ch>
   \date 2007-11-20
 */
#include <feel/feelconfig.h>

#include <stokes_stabilized.hpp>

std::string
createDomain( double h )
{
    std::ostringstream ostr;
    ostr << "Mesh.MshFileVersion = 1;\n"
         << "h=" << h << ";\n"
         << "Point(1) = {-1,-1,0.0,h};\n"
         << "Point(2) = { 1,-1,0.0,h};\n"
         << "Point(3) = { 1, 1,0.0,h};\n"
         << "Point(4) = {-1, 1,0.0,h};\n"
         << "Line(1) = {1,2};\n"
         << "Line(2) = {2,3};\n"
         << "Line(3) = {3,4};\n"
         << "Line(4) = {4,1};\n"
         << "Line Loop(7) = {4,3,2,1};\n"
         << "Plane Surface(8) = {7};\n"
         << "Physical Line(1) = {1,2,3,4};\n"
         << "Physical Surface(8) = {8};\n";
    return ostr.str();
}

int
main( int argc, char** argv )
{
    MPI_Init( &argc, &argv );

    Feel::StokesStabilized<1> test( argc, argv, makeAbout(), makeOptions() );

    test.run();
}




