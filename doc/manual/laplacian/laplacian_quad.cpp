/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-05-06

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
   \file laplacian_quad.cpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-05-06
 */
#include <feel/feel.hpp>

int main(int argc, char**argv )
{
    using namespace Feel;

	Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="laplacian_quad",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    //const std::string geofilename = option( _name="gmsh.filename" ).as<std::string>();
    //auto thegeo = geo( _filename=geofilename );
    //double alpha = thegeo->getGeoParameter("alpha");
    //double beta = thegeo->getGeoParameter("beta");
    double alpha = 0.45;
    double beta = 0.65;

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    auto Vh = Pch<1>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();
    auto vars=symbols<2>();
    auto& x = vars[0];
    auto& y = vars[1];
    ex uex = 16*x*(1-x)*y*(1-y)*sin(x + y-alpha);
    auto rhs=-laplacian(uex,vars);

    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),
                  _expr=chi(Px()+Py()-alpha > 0)*expr(rhs,vars)*id(v),
                  _quad=_Q<20>() );

    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate(_range=elements(mesh),
                  _expr=gradt(u)*trans(grad(v)) );
    a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u,
          _expr=constant(0.) );
    a.solve(_rhs=l,_solution=u);

    auto g = chi(Px()+Py()-alpha > 0)*expr(uex, vars );
    auto gradg = chi(Px()+Py()-alpha > 0)*expr<1,2,20>(grad(uex,vars), vars );
    LOG(INFO) << "grad(g)="<< grad(uex,vars) << "\n";
    LOG(INFO) << "H1 error omega : " << normH1( _range=elements(mesh),
                                                _expr=(idv(u)-g)*(idv(u)-g),
                                                _grad_expr=(gradv(u)-gradg)*trans(gradv(u)-gradg),
                                                _quad=_Q<20>() ) << "\n";

    LOG(INFO) << "H1 error omega0 : " << normH1( _range=markedelements(mesh,"Omega0"),
                                                 _expr=(idv(u)-g)*(idv(u)-g),
                                                 _grad_expr=(gradv(u)-gradg)*trans(gradv(u)-gradg),
                                                 _quad=_Q<20>() ) << "\n";
    v = project( _space=Vh, _expr=g );

    std::cout << "values with marker WEST : " << u.extractValuesWithMarker( "WEST", backend() ) << "\n";
    std::cout << "values without marker WEST : " << u.extractValuesWithoutMarker( "WEST", backend() ) << "\n";

    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    e->add( "exact", v );
    e->save();
    return 0;
}
