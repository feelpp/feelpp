/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-07

  Copyright (C) 2008-2009 Université Joseph Fourier (Grenoble I)

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
//! [all]
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/integrate.hpp>
#include <feel/feelfilters/exporter.hpp>


using namespace Feel;

int
main( int argc, char** argv )
{
    // Initialize Feel++ Environment
    Environment env( _argc=argc, _argv=argv,
                     _about=about( _name="myintegrals" ,
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ) );

    //! [mesh]
    // create the mesh (specify the dimension of geometric entity)
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    //! [mesh]

    //! [expression]
    // our function to integrate
    auto f = expr( option(_name="functions.g").as<std::string>(), Symbols{"x","y"});
    //! [expression]

    //! [integrals]
    // compute integral of f (global contribution)
    auto intf_1 = integrate( _range = elements( mesh ),
                                 _expr = f ).evaluate();

    // compute integral f on boundary
    double intf_3 = integrate( _range = boundaryfaces( mesh ),
                               _expr = f ).evaluate();

    if ( Environment::isMasterRank() == 0 )
        std::cout << "int global ; local ; boundary" << std::endl
                  << intf_1 << ";" << intf_2 << ";" << intf_3 << std::endl;
    //! [integrals]
}
//! [all]
