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
/// [all]
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/integrate.hpp>
#include <feel/feelvf/ginac.hpp>
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

    /// [mesh]
    // create the mesh (specify the dimension of geometric entity)
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    /// [mesh]

    /// [expression]
    // our function to integrate
    auto g = expr( soption(_name="functions.g") );
    /// [expression]

    /// [integrals]
    // compute integral of g (global contribution): \(\int_{\partial \Omega} g\)
    auto intf_1 = integrate( _range = elements( mesh ),
                                 _expr = g ).evaluate();

    // compute integral g on boundary: \( \int_{\partial \Omega} g \)
    auto intf_2 = integrate( _range = boundaryfaces( mesh ),
                             _expr = g ).evaluate();

    // compute integral of grad f (global contribution): \( \int_{\Omega} \nabla g \)
    auto grad_g = grad<2>(g);
    auto intgrad_f = integrate( _range = elements( mesh ),
                                _expr = grad_g ).evaluate();

    // only the process with rank 0 prints to the screen to avoid clutter
    if ( Environment::isMasterRank() )
        std::cout << "int_Omega " << g << " = " << intf_1  << std::endl
                  << "int_{boundary of Omega} " << g << " = " << intf_2 << std::endl
                  << "int_Omega grad " << g << " = "
                  << "int_Omega  " << grad_g << " = "
                  << intgrad_f  << std::endl;
    /// [integrals]
}
/// [all]
