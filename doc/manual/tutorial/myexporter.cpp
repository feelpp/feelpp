/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Guillaume Dollé <guillaume.dolle@math.unistra.fr>

   Date 2013-02-18

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
//! [all]
#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>
#include <feel/feelfilters/unitcircle.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/projectors.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/stdmathfunctors.hpp>
#include <feel/feelvf/geometricdata.hpp>

int main(int argc, char**argv )
{
  using namespace Feel;
  Environment env( _argc=argc, _argv=argv,
      _about=about(_name="myexporter",
        _author="Christophe Prud'homme",
        _email="christophe.prudhomme@feelpp.org"));

  //! [mesh]
  auto mesh = unitCircle<2>(); // circle - geometrical order: 2
  //! [mesh]
  //! [P1_mesh]
  auto meshp1 = unitCircle<1>(); // circle - geometrical order: 1
  //! [P1_mesh]

  //! [space]
  auto Xh = Pch<2>( mesh ); // \( \mathbb{p}_2 \) space
  //! [space]

  //! [function]
  auto v = project( _space=Xh, _range=elements(mesh),
      _expr=sin(pi*Px()));
  //! [function]

  //! [exporter]	
  auto exhi = exporter( _mesh=mesh, _name="exhi" );
  auto exlo = exporter( _mesh=meshp1, _name="exlo" );
  auto exhilo = exporter( _mesh=lagrangeP1(_space=Xh)->mesh(),_name="exhilo");
  //! [exporter]	

  //! [adding]
  int max = 10; double dt = 0.1;
  double time = 0;
  for (int i = 0; i<max; i++)
  {
    exhilo->step( time )->add( "vhilo", v );
    exlo->step( time )->add( "vlo", v );
    exhi->step( time )->add( "vhi", v );
    time += dt;
    //! [save]	
    exhi->save();
    exlo->save();
    exhilo->save();
    //! [save]	
  }
  //! [adding]	

}
//! [all]
