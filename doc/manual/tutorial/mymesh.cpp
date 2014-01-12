/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
	     Guillaume Dollé <guillaume.dolle@math.unistra.fr>

  Date 2013-02-11

  Copyright (C) 2008-2013
       Université Joseph Fourier (Grenoble I)
       Université de Strasbourg

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   \file mymesh.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>,
                 Guillaume Dollé <guillaume.dolle@math.unistra.fr>
   \date 2013-02-11
   This program show how to create a mesh and export results.
 */
#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelfilters/exporter.hpp>
using namespace Feel;

/**
 * Entry point
 */
//\code
//# marker_main #
int main( int argc, char** argv )
{
    // initialize Feel++ Environment
    Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about( _name="mymesh" ,
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ) );

    // create a mesh with GMSH using Feel++ geometry tool
    auto mesh = unitSquare();

    // export results for post processing
    auto e = exporter( _mesh=mesh );
    e->addRegions();
    e->save();

}   // main
//# endmarker_main #
//\endcode
