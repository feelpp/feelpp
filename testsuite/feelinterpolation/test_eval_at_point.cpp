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
template <uint16_type OrderPoly = 1, uint16_type DimGeo = 2>
void run()
{
    BOOST_TEST_MESSAGE( "test_eval_at_point D=" << DimGeo << " P=" << OrderPoly << "..." );
    Environment::changeRepository( boost::format( "%1%/D%2%/P%3%" ) % Environment::about().appName() % DimGeo % OrderPoly );
    auto mesh = loadMesh(_mesh = new Mesh<Simplex<DimGeo>>);
    auto Vh = Pchv<OrderPoly>( mesh );
    auto u = Vh->element();
    std::string e_str;
    if ( DimGeo == 2 )
        e_str = "{x*x*y,x*y*y}:x:y";
    else if ( DimGeo == 3 )
        e_str = "{x*x*y*z,x*y*y*z,x*y*z*z}:x:y:z";

    auto e = expr<DimGeo,1>( e_str );
    u = vf::project(Vh,elements(mesh), e );

    node_type pt(DimGeo);
    pt[0] = 0.5;
    if ( DimGeo >= 2 )
        pt[1] = 0.5;
    if ( DimGeo >= 3 )
        pt[2] = 0.5;
    auto eval = u(pt)[0];

    if ( DimGeo >= 2 )
    {
        e.setParameterValues( { { "x", 0.5 },{ "y", 0.5 } } );
        if ( DimGeo >= 3 )
            e.setParameterValues( { { "x", 0.5 },{ "y", 0.5 }, { "z", 0.5 } } );
        auto sol = e.evaluate();

        double min = doption(_name="gmsh.hsize");

#if USE_BOOST_TEST
        if ( OrderPoly < 4 )
            BOOST_CHECK_SMALL( (eval-sol).cwiseQuotient(sol).norm(), min );
        else
            BOOST_CHECK_SMALL( (eval-sol).cwiseQuotient(sol).norm(), 1e-13 );
#else
        CHECK( (eval-sol).norm() < min ) << "Error in evaluation at point " << pt
                                         << " value = [ " << eval << " ] expected = [ " << sol << " ] ";
#endif
        BOOST_TEST_MESSAGE( "test_eval_at_point D=" << DimGeo << " P=" << OrderPoly << "done" );
    }
}
}

#if USE_BOOST_TEST
FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( form_eval_at_point )

BOOST_AUTO_TEST_CASE( form_eval_at_point_P1_D2 )
{
//    for(int i = 0; i < 10; ++i )
        test_eval_at_point::run<1,2>();
}
BOOST_AUTO_TEST_CASE( form_eval_at_point_P1_D3 )
{
//    for(int i = 0; i < 5; ++i )
        test_eval_at_point::run<1,3>();
}
BOOST_AUTO_TEST_CASE( form_eval_at_point_P4_D2 )
{
    test_eval_at_point::run<4,2>();
}
BOOST_AUTO_TEST_CASE( form_eval_at_point_P4_D3 )
{
    test_eval_at_point::run<4,3>();
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
