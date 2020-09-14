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
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>
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

BOOST_AUTO_TEST_CASE( test2 )
{
    using namespace Feel;

    auto a1 = expr( "3*u^2:u" );
    auto a2 = expr( "2*u:u" );
    auto a3 = expr( "u*a1:u:a1" );

    auto my_eval_diff = [&a1,&a2,&a3]( auto && a4 )
        {
            auto a5base = expr( "5*a4:a4" );
            auto a5 = expr( a5base, symbolExpr("a1",a1), symbolExpr("a2",a3), symbolExpr("a3",a3), symbolExpr("a4",a4) );
            a5.setParameterValues( { "u", 2 } );

            auto diff_a5_u = a5.template diff<1>( "u" );
            BOOST_CHECK_CLOSE( diff_a5_u.evaluate()(0,0),1520,1e-10 );

        };

    my_eval_diff( a1*a1+a2*a2 );
    auto myvec = vec( a1, a2 );
    my_eval_diff( inner( myvec ) );
    my_eval_diff( inner( myvec, vec( a1, a2 )  ) );
    my_eval_diff( pow( inner( myvec, mpl::int_<InnerProperties::SQRT>() ), cst(2.) ) );
    my_eval_diff( pow( inner( myvec, vec( a1, a2 ), mpl::int_<InnerProperties::SQRT>() ), cst(2.) ) );
    my_eval_diff( pow(myvec(0,0),2.0) + pow(myvec(1,0),2.0) );
}

BOOST_AUTO_TEST_CASE( test3 )
{
    using namespace Feel;
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<3,1>>);
    auto Vh = Pch<2>( mesh );
    auto u = Vh->element( inner(P()) );

    auto e1 = expr( "u*u:u");
    auto e2 = expr( "2*e1:e1");
    auto e3base = expr( "3*e2:e2");

    auto e3 = expr( e3base, symbolExpr("e1",e1), symbolExpr("e2",e2), symbolExpr("u",inner(P()) ) );
    auto diff_e3_x = e3.diff<1>( "x" );
    auto diff_e3_y = e3.diff<1>( "y" );
    auto diff_e3_z = e3.diff<1>( "z" );
    auto grad_e3 = grad<3>( e3 );
    auto diff_e3_x_exact = 3*4*inner(P())*2*Px();
    auto diff_e3_y_exact = 3*4*inner(P())*2*Py();
    auto diff_e3_z_exact = 3*4*inner(P())*2*Pz();
    auto grad_e3_exact = trans(vec(diff_e3_x_exact,diff_e3_y_exact,diff_e3_z_exact));

    double error_diff_e3_x = normL2(_range=elements(mesh),_expr= diff_e3_x - diff_e3_x_exact );
    double error_diff_e3_y = normL2(_range=elements(mesh),_expr= diff_e3_y - diff_e3_y_exact );
    double error_diff_e3_z = normL2(_range=elements(mesh),_expr= diff_e3_z - diff_e3_z_exact );
    double error_grad_e3 = normL2(_range=elements(mesh),_expr= grad_e3 - grad_e3_exact );
    BOOST_CHECK_SMALL( error_diff_e3_x, 1e-10 );
    BOOST_CHECK_SMALL( error_diff_e3_y, 1e-10 );
    BOOST_CHECK_SMALL( error_diff_e3_z, 1e-10 );
    BOOST_CHECK_SMALL( error_grad_e3, 1e-10 );

    auto e3b = expr( e3base, symbolExpr("e1",e1), symbolExpr("e2",e2), symbolExpr("u",idv(u)) );
    auto diff_e3b_x = e3b.diff<1>( "x" );
    auto diff_e3b_y = e3b.diff<1>( "y" );
    auto diff_e3b_z = e3b.diff<1>( "z" );
    auto grad_e3b = grad<3>( e3b );
    double error_diff_e3b_x = normL2(_range=elements(mesh),_expr= diff_e3b_x - diff_e3_x_exact );
    double error_diff_e3b_y = normL2(_range=elements(mesh),_expr= diff_e3b_y - diff_e3_y_exact );
    double error_diff_e3b_z = normL2(_range=elements(mesh),_expr= diff_e3b_z - diff_e3_z_exact );
    double error_grad_e3b = normL2(_range=elements(mesh),_expr= grad_e3b - grad_e3_exact );
    BOOST_CHECK_SMALL( error_diff_e3b_x, 1e-10 );
    BOOST_CHECK_SMALL( error_diff_e3b_y, 1e-10 );
    BOOST_CHECK_SMALL( error_diff_e3b_z, 1e-10 );
    BOOST_CHECK_SMALL( error_grad_e3b, 1e-10 );
}


BOOST_AUTO_TEST_SUITE_END()
