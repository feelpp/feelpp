/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
	     Guillaume Dollé <guillaume.dolle@math.unistra.fr>
 
  Date 2013-02-25

  Copyright (C) 2013 Université de Strasbourg

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
   \file myadvection.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>,
                 Guillaume Dollé <guillaume.dolle@math.unistra.fr>
   \date 2013-02-25
   This program show how to solve the advection problem.
 */
#include <feel/feel.hpp>
using namespace Feel;

/**
 * Entry point
 */
//\code
//# marker_main #
int
main( int argc, char** argv )
{
   
    // Initialize Feel++ Environment
    Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="myadvection",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org") );
    // create mesh
    auto mesh = unitSquare();

    // function space
    auto Xh = Pch<1>( mesh );
    auto u = Xh->element( "u" );
    auto v = Xh->element( "v" );

    // diffusion coeff.
    double epsilon = 1;
    // reaction coeff.
    double mu = 1;    
    auto beta = vec( cst(1.),
                     cst(1.) );
    auto f = cst(1.);

    // left hand side
    auto a = form2( _test=Xh, _trial=Xh );
    a += integrate( _range=elements( mesh ),
                    _expr=( epsilon*gradt( u )*trans( grad( v ) )
                         + ( gradt( u )*beta )*id(v)
                         + mu*idt( u )*id( v ) ) );

    // right hand side
    auto l = form1( _test=Xh );
    l+= integrate( _range=elements( mesh ), _expr=f*id( v ) );

    // boundary condition
    a += on( _range=boundaryfaces( mesh ), _rhs=l, _element=u,
             _expr=cst(0.) );

    // solve the system
    a.solve( _rhs=l, _solution=u );

    // export results
    auto e = exporter( _mesh=mesh );
    e->add("u",u);
    e->save();
} // end main
//# endmarker_main #
//\endcode





