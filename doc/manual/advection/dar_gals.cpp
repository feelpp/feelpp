/* -*- coding: utf-8; mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-07-23

  Copyright (C) 2013 Université de Strasbourg

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
   \file gals.cpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-07-23
 */

#include <feel/feel.hpp>
int main( int argc, char**argv )
{
    //# marker1 #
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="dar_gals",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));
    //# endmarker1 #

    //# marker2 #
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    auto Vh = Pch<3>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();
    //# endmarker2 #

    auto vars = symbols<2>();
    auto beta_x = expr( option(_name="beta_x",_prefix="functions").as<std::string>(), vars );
    auto beta_y = expr( option(_name="beta_y",_prefix="functions").as<std::string>(), vars );
    auto beta = vec( beta_x, beta_y );
    auto epsilon = expr( option(_name="epsilon",_prefix="functions").as<std::string>(), vars );
    auto gamma = expr( option(_name="gamma",_prefix="functions").as<std::string>(), vars );

    //# marker3 #
    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),
                  _expr=id(v));


    auto a = form2( _trial=Vh, _test=Vh);
    a = integrate(_range=elements(mesh),
                  _expr=(gradt(u)*beta)*id(v)+epsilon*gradt(u)*trans(grad(v))+gamma*idt(u)*id(v) );
    a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u,
          _expr=constant(0.) );
    a.solve(_rhs=l,_solution=u);
    //# endmarker3 #

    //# marker4 #
    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    e->save();
    return 0;
    //# endmarker4 #
}
