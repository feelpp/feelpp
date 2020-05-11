/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 8 may 2020

 Copyright (C) 2020 Feel++ Consortium

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

// give a name to the testsuite
#define BOOST_TEST_MODULE vf symbolic diff expr testsuite

#include <feel/feelcore/testsuite.hpp>

// #include <feel/feeldiscr/functionspace.hpp>
// #include <feel/feelfilters/loadmesh.hpp>
// #include <feel/feelfilters/exporter.hpp>
// #include <feel/feelmesh/concatenate.hpp>
// #include <feel/feelfilters/unitcube.hpp>
#include <feel/feelvf/vf.hpp>


FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( symbolic_diff )


BOOST_AUTO_TEST_CASE( test1 )
{
    using namespace Feel;

    auto e1 = expr( "u*u:u");
    e1.setParameterValues( { "u", 2 } );
    auto e2 = expr( "2*e1:e1");
    auto e3base = expr( "1/e2:e2");  // -> - e1p /e1^2
    auto e3 = expr( e3base, symbolExpr("e1",e1), symbolExpr("e2",e2) );
    BOOST_CHECK( e3.hasSymbolDependency( "u" ) );
    BOOST_CHECK( e3.hasSymbolDependency( "e1" ) );
    BOOST_CHECK( e3.hasSymbolDependency( "e2" ) );
    auto diff_e3_u = e3.diff<1>( "u" );
    BOOST_CHECK_CLOSE( diff_e3_u.evaluate()(0,0), -0.125, 1e-10 );

    auto f1 = expr<2,1>( "{cos(u),u*u}:u");
    f1.setParameterValues( { "u", 2 } );
    auto f2 = expr( "2*f1_1:f1_1");
    auto f3base = expr( "1/f2:f2");
    auto f3 = expr( f3base, symbolExpr("f1",f1,SymbolExprComponentSuffix(2,1)), symbolExpr("f2",f2) );
    auto diff_f3_u = f3.diff<1>( "u" );
    BOOST_CHECK_CLOSE( diff_f3_u.evaluate()(0,0), -0.125, 1e-10 );

    auto g1 = expr( "u:u");
    auto g2 = expr( "2*u:u");
    auto g3 = cst(1.)/(g1*g2);
    auto g4base = expr( "g3:g3");
    auto g4 = expr( g4base, symbolExpr("g3",g3) );
    g4.setParameterValues( { "u", 2 } );
    auto diff_g4_u = g4.diff<1>( "u" );
    BOOST_CHECK_CLOSE( diff_g4_u.evaluate()(0,0), -0.125, 1e-10 );

}

BOOST_AUTO_TEST_SUITE_END()
