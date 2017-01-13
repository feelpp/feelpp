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
//! [all]
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/geo.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/exporter.hpp>
using namespace Feel;

int main( int argc, char** argv )
{
	// initialize Feel++ Environment
	Environment env( _argc=argc, _argv=argv,
                     _about=about( _name="mymeshFromString" ,
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ) );

	//! [create]
	// create a mesh with GMSH using Feel++ geometry tool
	std::ostringstream str;
	str
		<< "h = 0.1;\n"
		<< "Point(1) = {0,0,0,h};\n"
		<< "Point(2) = {1,0,0,h};\n"
		<< "Point(3) = {0,1,0,h};\n"
		<< "Point(4) = {-1,0,0,h/5};\n"
		<< "Point(5) = {0,-1,0,h/5};\n"
		<< "Line(1) = {2,3};\n"
		<< "Circle(2) = {3,1,4};\n"
		<< "Circle(3) = {4,1,5};\n"
		<< "Circle(4) = {5,1,2};\n"
		<< "Line Loop (5) = {1,2,3,4};\n"
		<< "Plane Surface(6) = {5};\n"
		<< "Physical Line(1) = {2};\n"
		<< "Physical Line(2) = {3};\n"
		<< "Physical Line(3) = {1};\n"
		<< "Physical Line(4) = {4};\n"
		<< "Physical Surface(\"Mat1\") = {6};\n";

	auto mesh2 = createGMSHMesh(_mesh=new  Mesh<Simplex<2>>,
                                _desc = geo(_filename="aMesh",_desc=str.str()));
	//! [create]

	//! [export]
	// export results for post processing
	auto e = exporter( _mesh=mesh2 );
    e->add("pid", regionProcess(Pdh<0>(mesh2,true)) );
	e->save();
	//! [export]

}   // main
//! [all]
