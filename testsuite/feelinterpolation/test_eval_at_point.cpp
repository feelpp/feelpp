/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

	 This file is part of the Feel library

	 Author(s): Vincent Chabannes <vincent.chabannes@imag.fr>
Date: 2011-06-19

Copyright (C) 2011 Universite Joseph Fourier (Grenoble I)

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
	\file test_form_interpolation.cpp
	\author Vincent Chabannes <vincent.chabannes@imag.fr>
	\date 2011-06-19
	*/

#define USE_BOOST_TEST 1

#define BOOST_TEST_MODULE eval_at_point testsuite

#include <testsuite/testsuite.hpp>

#include <feel/options.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/geotool.hpp>

using namespace Feel;
using namespace Feel::vf;

typedef Application Application_type;
typedef boost::shared_ptr<Application_type> Application_ptrtype;

namespace test_eval_at_point
{
	inline
		po::options_description
		makeOptions()
		{
			po::options_description desc_options( "test form interpolation options" );
			desc_options.add_options()
				( "x", po::value<double>()->default_value( 0.23 ), " x coord" )
				( "y", po::value<double>()->default_value( 0.52 ), " y coord" );
			;
			return desc_options.add( Feel::feel_options() );//.add( backend_options( "pressure" ) );
		}
	inline
		AboutData
		makeAbout()
		{
			AboutData about( "Test_eval_at_point" ,
					"Test_eval_at_point" ,
					"0.1",
					"test evaluation at point",
					Feel::AboutData::License_GPL,
					"Copyright (c) 2011 Universite Joseph Fourier" );

			about.addAuthor( "Vincent HUBER", "developer", "vincent.huber@cemosis.fr", "" );
			return about;
		}

	template <uint16_type OrderPoly = 1, uint16_type OrderGeo = 2>
		void run()
		{
			auto mesh = loadMesh(_mesh = new Mesh<Simplex<OrderGeo>>); 
			auto Vh = Pchv<OrderPoly>( mesh );
			auto u = Vh->element();
			u = vf::project(Vh,
					elements(mesh),
					vec( 
						Px()*Px()*Py(), 
						Px()*Py()*Py() ));
			node_type pt(2); 
			pt[0] = option(_name="x").as<double>(); 
			pt[1] = option(_name="y").as<double>();
			auto eval = u(pt);
			auto val1 = eval(0,0,0); auto sol1 = pt[0]*pt[0]*pt[1];  
			auto val2 = eval(1,0,0); auto sol2 = pt[0]*pt[1]*pt[1];
			double min = option(_name="gmsh.hsize").as<double>();

#if USE_BOOST_TEST
			BOOST_CHECK_SMALL( (val1-sol1)/sol1, min );
			BOOST_CHECK_SMALL( (val2-sol2)/sol2, min );
#else
			assert((val1-sol1)/sol1<min);
			assert((val2-sol2)/sol2<min);
#endif
		}
}

#if USE_BOOST_TEST
FEELPP_ENVIRONMENT_WITH_OPTIONS( test_eval_at_point::makeAbout(),
		test_eval_at_point::makeOptions() )

BOOST_AUTO_TEST_SUITE( form_eval_at_point )

BOOST_AUTO_TEST_CASE( form_eval_at_point )
{
	test_eval_at_point::run();
}
BOOST_AUTO_TEST_SUITE_END()

#else
int main(int argc, char**argv )
{
	Feel::Environment env( argc, argv );
	auto theApp = Application_ptrtype( new Application_type );
	if ( theApp->vm().count( "help" ) )
	{
		std::cout << theApp->optionsDescription() << "\n";
		exit( 0 );
	}
	test_eval_at_point::run<1,1>( theApp );
}

#endif
