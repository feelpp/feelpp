/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-05-28

  Copyright (C) 2014 Feel++ Consortium

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
   \file magnetostatic.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-05-28
 */
#include <feel/feel.hpp>
#include <feel/feeldiscr/ned1h.hpp>

#define FEELPP_DIM 2
int main(int argc, char**argv )
{
    /// [marker1]
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="magnetostatic",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<FEELPP_DIM>> );

    auto Nh = Ned1h<0>( mesh );
    auto a = form2( _trial=Nh, _test=Nh );
    auto c = doption("parameters.c");
    auto f = expr<FEELPP_DIM,1>(soption("functions.f"));
    auto e = expr<FEELPP_DIM,1>(soption("functions.e"));
    auto u = Nh->element();
    auto v = Nh->element();

    a = integrate(_range=elements(mesh), _expr=c*trans(idt(u))*id(v)+curlxt(u)*curlx(v));
    auto l = form1( _test=Nh );
    l = integrate( _range=elements(mesh), _expr=trans(f)*id(v));
    a.solve(_rhs=l,_solution=u);
    auto err = normL2( _range=elements(mesh), _expr=idv(u)-e );
    if ( Environment::isMasterRank() )
    {
        std::cout << "L2 error = " << err << "\n";
    }
    auto ex = exporter( _mesh=mesh );
    ex->add( "u", u );
    ex->save();

}
