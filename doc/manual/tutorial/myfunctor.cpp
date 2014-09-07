/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4

   This file is part of the Feel library

   Author(s): Vincent HUBER <vincent.huber@cemosis.fr>

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
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
using namespace Feel;
namespace Feel
{
	struct myFunctor
	{
		static const size_type context = vm::JACOBIAN|vm::POINT;
		typedef double value_type;
		typedef Feel::uint16_type uint16_type;
		static const uint16_type rank = 0;
		static const uint16_type imorder = 1;
		static const bool imIsPoly = true;
		double operator()( uint16_type, uint16_type, ublas::vector<double> const& x, ublas::vector<double> const& n ) const
		{
			return x[i];
		}
		int val = 0;
		void setVal(int i){val = i;}
	};
}

int main(int argc, char**argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="myexpression",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    //! [mesh]
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    //! [mesh]

    //! [space]
    auto spaces= Pch<1>(mesh); 
    auto spacev= Pchv<1>(mesh); 
    //! [space]
	
		//! [functors]
		myFunctor functor1, functor2;
		functor1.setVal(0);
		functor2.setVal(1);
		//! [functors]

		//! [projection]
		auto u = vf::project(spaces,elements(mesh),idf(functor1)); // Will contain x
		auto v = vf::project(spaces,elements(mesh),idf(functor2)); // will contain y
		auto V = vf::project(spacev,elements(mesh),vec(idf(functor1),idf(functor2)) );
		//! [projection]
}
//! [all]
