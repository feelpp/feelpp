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


#include <feel/feel.hpp>
using namespace Feel;

/// [marker_main]
int main(int argc, char**argv )
{
	po::options_description app_options( "MyBackend options" );
	app_options.add_options()
		( "value",
			po::value<double>() -> default_value(4.2),
			"a 'double' with default value" );
	app_options.add(Feel::feel_options());
	app_options.add(backend_options("myBackend"));
	
	Environment env( _argc=argc, _argv=argv,
			_desc = app_options, 
			_about=about(_name="mybackend",
				_author="Feel++ Consortium",
				_email="feelpp-devel@feelpp.org"));

	// create a backend
	boost::shared_ptr<Backend<double>> myBackend(backend(_name="myBackend"));

	// create the mesh
	auto mesh = loadMesh(_mesh=new Mesh<Simplex< 2 > > );

	// function space
	auto Vh = Pch<2>( mesh );

	// element in Vh
	auto u  = Vh->element();
	auto u1 = Vh->element(); //computed with built-in backend 
	auto u2 = Vh->element(); //computed with home maid backend

	// left hand side
	auto a = form2( _trial=Vh, _test=Vh );
	a = integrate(_range=elements(mesh),
			_expr=expr(soption("functions.alpha"))*trace(gradt(u)*trans(grad(u))) );

	// right hand side
	auto l = form1( _test=Vh );
	l = integrate(_range=elements(mesh),
			_expr=expr(soption("functions.f"))*id(u));

	// BC
	a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u,
			_expr=expr(soption("functions.g")));

	// solve a(u,v)=l(v)
	std::cout << "With standard backend\n";
	a.solve(_rhs=l,_solution=u1);
	std::cout << "With home maid backend\n";
	a.solveb(_rhs=l,_solution=u2, _backend=myBackend);
	//# endmarker_main #

	// save results
	auto e = exporter( _mesh=mesh );
	e->step(0) -> add( "u1", u1 );
	e->step(0) -> add( "u2", u2 );
	e->step(1) -> add( "u1", u1 );
	e->step(1) -> add( "u2", u2 );
	e->save();
}
/// [marker_main]
