/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
	     Guillaume Dollé <guillaume.dolle@math.unistra.fr>
 
  Date 2013-02-19

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
   \file mystokes.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>,
                 Guillaume Dollé <guillaume.dolle@math.unistra.fr>
   \date 2013-02-19
   This program show how to solve Poiseuille flow equations (stokes).
 */

#include <feel/feel.hpp>
using namespace Feel;

/*
 * Entry point
 */
//\code
//# marker_main #
int main(int argc, char**argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="mystokes",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    // create the mesh
    auto mesh = unitSquare();

    // function space
    auto Vh = THch<2>( mesh );

    // element U=(u,p) in Vh
    auto U = Vh->element();
    auto u = U.element<0>();
    auto p = U.element<1>();

    // left hand side
    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate(_range=elements(mesh),
                  _expr=trace(gradt(u)*trans(grad(u))) );

    a+= integrate(_range=elements(mesh),
                  _expr=-div(u)*idt(p)+divt(u)*id(p));

    // right hand side
    auto l = form1( _test=Vh );

    // boundary condition
    a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u,
          _expr=Py()*(1-Py()) );

    // solve a(u,v)=l(v)
    a.solve(_rhs=l,_solution=U);

    // save results
    auto e = exporter( _mesh=mesh, _name="poiseuille" );
    e->add( "u", u );
    e->add( "p", p );
    e->save();
}
//# endmarker_main #
//\endcode

