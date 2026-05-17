/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 28 nov 2019

 Copyright (C) 2019 Feel++ Consortium

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
#define BOOST_TEST_MODULE expr evaluate testsuite

#include <feel/feelcore/testsuite.hpp>

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchm.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/eig.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( expr_evaluate )

BOOST_AUTO_TEST_CASE( test1 )
{
    using namespace Feel;
    double tolerance = 1e-10;

    // constant
    auto ex1 = cst(3.14);
    auto eval1 = ex1.evaluate();
    BOOST_CHECK( eval1.rows() == 1 && eval1.cols() == 1 );
    BOOST_CHECK_CLOSE( eval1(0,0), 3.14, tolerance );

    // operator *
    auto ex2 = cst(3.14)*cst(2.0);
    auto eval2 = ex2.evaluate();
    BOOST_CHECK( eval2.rows() == 1 && eval2.cols() == 1 );
    BOOST_CHECK_CLOSE( eval2(0,0), 3.14*2, tolerance );

    // idendity matrix
    auto ex3 = eye<3,3>();
    auto eval3 = ex3.evaluate();
    BOOST_CHECK( eval3.rows() == 3 && eval3.cols() == 3 );
    for (int i=0;i<eval3.rows();++i )
    {
        for (int j=0;j<eval3.cols();++j )
        {
            if ( i == j )
                BOOST_CHECK_CLOSE( eval3(i,j), 1.0, tolerance );
            else
                BOOST_CHECK_SMALL( eval3(i,j), tolerance );
        }
    }

    // symbolic expression scalar
    auto ex4 = expr("2*3");
    auto eval4 = ex4.evaluate();
    BOOST_CHECK( eval4.rows() == 1 && eval4.cols() == 1 );
    BOOST_CHECK_CLOSE( eval4(0,0), 6.0, tolerance );

    // symbolic expression vectorial
    auto ex5 = expr<3,1>("{1,2,3}");
    auto eval5 = ex5.evaluate();
    BOOST_CHECK( eval5.rows() == 3 && eval5.cols() == 1 );
    for (int i=0;i<3;++i)
        BOOST_CHECK_CLOSE( eval5(i,0), 1.0+i, tolerance );

    // expr component
    auto ex6 = ex5(1,0);
    auto eval6 = ex6.evaluate();
    BOOST_CHECK( eval6.rows() == 1 && eval6.cols() == 1 );
    BOOST_CHECK_CLOSE( eval6(0,0), 2.0, tolerance );

    // product scalar by matrix
    auto ex7 = ex4*ex3;
    auto eval7 = ex7.evaluate();
    BOOST_CHECK( eval7.rows() == 3 && eval7.cols() == 3 );
    for (int i=0;i<eval7.rows();++i )
    {
        for (int j=0;j<eval7.cols();++j )
        {
            if ( i == j )
                BOOST_CHECK_CLOSE( eval7(i,j), 6.0, tolerance );
            else
                BOOST_CHECK_SMALL( eval7(i,j), tolerance );
        }
    }

    // symbolic expression
    using mesh_t = Mesh<Simplex<2>>;
    auto mesh = loadMesh( _mesh = new mesh_t );
    auto Vh = Pch<2>( mesh );
    auto u = Vh->element();
    auto ex8base = expr("2*3+v:v");
    auto ex8 = expr( ex8base, symbolExpr("u",idv(u) ), symbolExpr("v",cst(5.) ) );
    auto eval8 = ex8.evaluate();
    BOOST_CHECK( eval8.rows() == 1 && eval8.cols() == 1 );
    BOOST_CHECK_CLOSE( eval8(0,0), 11.0, tolerance );

    // vec
    auto ex9 = vec(cst(1.),cst(2.),cst(3.));
    auto eval9 = ex9.evaluate();
    BOOST_CHECK( eval9.rows() == 3 && eval9.cols() == 1 );
    for (int i=0;i<eval9.rows();++i )
        BOOST_CHECK_CLOSE( eval9(i,0), 1.0+i, tolerance );

    auto ex9c = vec( 1., 2., 3. );
    auto eval9c = ex9c.evaluate();
    BOOST_CHECK( eval9c.rows() == 3 && eval9c.cols() == 1 );
    for ( int i = 0; i < eval9c.rows(); ++i )
        BOOST_CHECK_CLOSE( eval9c( i, 0 ), 1.0+i, tolerance );

    // matrix
    auto ex10 = mat<3,2>(cst(1.),cst(2.),
                         cst(3.),cst(4.),
                         cst(5.),cst(6.));
    auto eval10 = ex10.evaluate();
    BOOST_CHECK( eval10.rows() == 3 && eval10.cols() == 2 );
    for (int i=0;i<eval10.rows();++i )
        for (int j=0;j<eval10.cols();++j )
            BOOST_CHECK_CLOSE( eval10(i,j), 1.0+i*eval10.cols()+j, tolerance );

    auto ex10d = mat<2,2>( 1., 2.,
                           3., 4. );
    auto eval10d = ex10d.evaluate();
    BOOST_CHECK( eval10d.rows() == 2 && eval10d.cols() == 2 );
    BOOST_CHECK_CLOSE( eval10d(0,0), 1.0, tolerance );
    BOOST_CHECK_CLOSE( eval10d(0,1), 2.0, tolerance );
    BOOST_CHECK_CLOSE( eval10d(1,0), 3.0, tolerance );
    BOOST_CHECK_CLOSE( eval10d(1,1), 4.0, tolerance );

    auto ex10bLeft = mat<2,2>( cst(1.), cst(2.),
                               cst(3.), cst(4.) );
    auto ex10bRight = mat<2,2>( cst(5.), cst(6.),
                                cst(7.), cst(8.) );
    auto ex10b = ex10bLeft*ex10bRight;
    auto eval10b = ex10b.evaluate();
    BOOST_CHECK( eval10b.rows() == 2 && eval10b.cols() == 2 );
    BOOST_CHECK_CLOSE( eval10b(0,0), 19.0, tolerance );
    BOOST_CHECK_CLOSE( eval10b(0,1), 22.0, tolerance );
    BOOST_CHECK_CLOSE( eval10b(1,0), 43.0, tolerance );
    BOOST_CHECK_CLOSE( eval10b(1,1), 50.0, tolerance );

    auto ex10bHadamard = hadamard( ex10bLeft, ex10bRight );
    auto eval10bHadamard = ex10bHadamard.evaluate();
    BOOST_CHECK( eval10bHadamard.rows() == 2 && eval10bHadamard.cols() == 2 );
    BOOST_CHECK_CLOSE( eval10bHadamard(0,0), 5.0, tolerance );
    BOOST_CHECK_CLOSE( eval10bHadamard(0,1), 12.0, tolerance );
    BOOST_CHECK_CLOSE( eval10bHadamard(1,0), 21.0, tolerance );
    BOOST_CHECK_CLOSE( eval10bHadamard(1,1), 32.0, tolerance );

    auto ex10bRuntimeComponent = ex10b( 1, 0 );
    auto eval10bRuntimeComponent = ex10bRuntimeComponent.evaluate();
    BOOST_CHECK( eval10bRuntimeComponent.rows() == 1 && eval10bRuntimeComponent.cols() == 1 );
    BOOST_CHECK_CLOSE( eval10bRuntimeComponent( 0, 0 ), 43.0, tolerance );

    auto ex10bStaticComponent = component<1,0>( ex10b );
    auto eval10bStaticComponent = ex10bStaticComponent.evaluate();
    BOOST_CHECK( eval10bStaticComponent.rows() == 1 && eval10bStaticComponent.cols() == 1 );
    BOOST_CHECK_CLOSE( eval10bStaticComponent( 0, 0 ), 43.0, tolerance );

    auto ex10bChain = ex10b*ex10bLeft;
    auto eval10bChain = ex10bChain.evaluate();
    BOOST_CHECK( eval10bChain.rows() == 2 && eval10bChain.cols() == 2 );
    BOOST_CHECK_CLOSE( eval10bChain( 0, 0 ), 85.0, tolerance );
    BOOST_CHECK_CLOSE( eval10bChain( 0, 1 ), 126.0, tolerance );
    BOOST_CHECK_CLOSE( eval10bChain( 1, 0 ), 193.0, tolerance );
    BOOST_CHECK_CLOSE( eval10bChain( 1, 1 ), 286.0, tolerance );

    auto ex10bVec = ex10bLeft*vec( cst( 5.0 ), cst( 6.0 ) );
    auto eval10bVec = ex10bVec.evaluate();
    BOOST_CHECK( eval10bVec.rows() == 2 && eval10bVec.cols() == 1 );
    BOOST_CHECK_CLOSE( eval10bVec( 0, 0 ), 17.0, tolerance );
    BOOST_CHECK_CLOSE( eval10bVec( 1, 0 ), 39.0, tolerance );

    auto ex10bVecComponent = component<1,0>( ex10bVec );
    auto eval10bVecComponent = ex10bVecComponent.evaluate();
    BOOST_CHECK( eval10bVecComponent.rows() == 1 && eval10bVecComponent.cols() == 1 );
    BOOST_CHECK_CLOSE( eval10bVecComponent( 0, 0 ), 39.0, tolerance );

    auto ex10c = ex10bLeft( 0, 1 );
    auto eval10c = ex10c.evaluate();
    BOOST_CHECK( eval10c.rows() == 1 && eval10c.cols() == 1 );
    BOOST_CHECK_CLOSE( eval10c( 0, 0 ), 2.0, tolerance );

    auto ex10cStatic = component<0,1>( ex10bLeft );
    auto eval10cStatic = ex10cStatic.evaluate();
    BOOST_CHECK( eval10cStatic.rows() == 1 && eval10cStatic.cols() == 1 );
    BOOST_CHECK_CLOSE( eval10cStatic( 0, 0 ), 2.0, tolerance );

    auto ex9b = ex9( 2, 0 );
    auto eval9b = ex9b.evaluate();
    BOOST_CHECK( eval9b.rows() == 1 && eval9b.cols() == 1 );
    BOOST_CHECK_CLOSE( eval9b( 0, 0 ), 3.0, tolerance );

    auto ex9bStatic = component<2,0>( ex9 );
    auto eval9bStatic = ex9bStatic.evaluate();
    BOOST_CHECK( eval9bStatic.rows() == 1 && eval9bStatic.cols() == 1 );
    BOOST_CHECK_CLOSE( eval9bStatic( 0, 0 ), 3.0, tolerance );

    auto ex9Hadamard = hadamard( ex9, vec( cst( 4.0 ), cst( 5.0 ), cst( 6.0 ) ) );
    auto eval9Hadamard = ex9Hadamard.evaluate();
    BOOST_CHECK( eval9Hadamard.rows() == 3 && eval9Hadamard.cols() == 1 );
    BOOST_CHECK_CLOSE( eval9Hadamard( 0, 0 ), 4.0, tolerance );
    BOOST_CHECK_CLOSE( eval9Hadamard( 1, 0 ), 10.0, tolerance );
    BOOST_CHECK_CLOSE( eval9Hadamard( 2, 0 ), 18.0, tolerance );

    // transpose
    auto ex11 = trans(ex10)*ex9;
    auto eval11 = ex11.evaluate();
    BOOST_CHECK( eval11.rows() == 2 && eval11.cols() == 1 );
    BOOST_CHECK_CLOSE( eval11(0,0), 22, tolerance );
    BOOST_CHECK_CLOSE( eval11(1,0), 28, tolerance );

    auto ex11b = trans( ex9 );
    auto eval11b = ex11b.evaluate();
    BOOST_CHECK( eval11b.rows() == 1 && eval11b.cols() == 3 );
    BOOST_CHECK_CLOSE( eval11b(0,0), 1.0, tolerance );
    BOOST_CHECK_CLOSE( eval11b(0,1), 2.0, tolerance );
    BOOST_CHECK_CLOSE( eval11b(0,2), 3.0, tolerance );

    auto ex11c = trans( vec( cst(1.), cst(2.) ) )*mat<2,2>( cst(3.), cst(4.),
                                                             cst(5.), cst(6.) );
    auto eval11c = ex11c.evaluate();
    BOOST_CHECK( eval11c.rows() == 1 && eval11c.cols() == 2 );
    BOOST_CHECK_CLOSE( eval11c(0,0), 13.0, tolerance );
    BOOST_CHECK_CLOSE( eval11c(0,1), 16.0, tolerance );

    auto ex11d = ex11c*vec( cst(7.), cst(8.) );
    auto eval11d = ex11d.evaluate();
    BOOST_CHECK( eval11d.rows() == 1 && eval11d.cols() == 1 );
    BOOST_CHECK_CLOSE( eval11d(0,0), 219.0, tolerance );

    // inner
    auto ex12 = inner(ex9, /*4.0**/ones<3,1>());
    auto eval12 = ex12.evaluate();
    BOOST_CHECK( eval12.rows() == 1 && eval12.cols() == 1 );
    BOOST_CHECK_CLOSE( eval12(0,0), 6.0, tolerance );

    // trace
    auto ex13 = trace( mat<3,3>(cst(1.),cst(2.),cst(3.),
                                cst(4.),cst(5.),cst(6.),
                                cst(7.),cst(8.),cst(9.) ) );
    auto eval13 = ex13.evaluate();
    BOOST_CHECK( eval13.rows() == 1 && eval13.cols() == 1 );
    BOOST_CHECK_CLOSE( eval13(0,0), 15.0, tolerance );

    // det
    auto ex14 = det( eye<3,3>() );
    auto eval14 = ex14.evaluate();
    BOOST_CHECK( eval14.rows() == 1 && eval14.cols() == 1 );
    BOOST_CHECK_CLOSE( eval14(0,0), 1.0, tolerance );

    // inv
    auto ex15base = mat<3,3>(cst(1.),cst(0.),cst(0.),
                             cst(0.),cst(5.),cst(0.),
                             cst(0.),cst(0.),cst(9.) );
    auto ex15 = ex15base*inv(ex15base);
    auto eval15 = ex15.evaluate();
    BOOST_CHECK( eval15.rows() == 3 && eval15.cols() == 3 );
    for (int i=0;i<eval15.rows();++i )
    {
        for (int j=0;j<eval15.cols();++j )
        {
            if ( i == j )
                BOOST_CHECK_CLOSE( eval15(i,j), 1.0, tolerance );
            else
                BOOST_CHECK_SMALL( eval15(i,j), tolerance );
        }
    }

    // cross 2d
    auto ex16 = cross( vec(cst(1.),cst(2.)), vec(cst(3.),cst(4.)) );
    auto eval16 = ex16.evaluate();
    BOOST_CHECK( eval16.rows() == 1 && eval16.cols() == 1 );
    BOOST_CHECK_CLOSE( eval16(0,0), -2.0, tolerance );

    // cross 3d
    auto ex17 = cross( vec(cst(1.),cst(2.),cst(3.)), vec(cst(4.),cst(5.),cst(6.)) );
    auto eval17 = ex17.evaluate();
    BOOST_CHECK( eval17.rows() == 3 && eval17.cols() == 1 );
    BOOST_CHECK_CLOSE( eval17(0,0), 2.0*6.0-3.0*5.0, tolerance );
    BOOST_CHECK_CLOSE( eval17(1,0), 3.0*4.0-1.0*6.0, tolerance );
    BOOST_CHECK_CLOSE( eval17(2,0), 1.0*5.0-2.0*4.0, tolerance );

    // math function
    auto ex18 = sin( cst(1.5) );
    auto eval18 = ex18.evaluate();
    BOOST_CHECK( eval18.rows() == 1 && eval18.cols() == 1 );
    BOOST_CHECK_CLOSE( eval18(0,0), std::sin(1.5), tolerance );

    // pow
    auto ex19 = pow( cst(4.), cst(3.) );
    auto eval19 = ex19.evaluate();
    BOOST_CHECK( eval19.rows() == 1 && eval19.cols() == 1 );
    BOOST_CHECK_CLOSE( eval19(0,0), std::pow(4,3), tolerance );

    // eig (symmetric mat)
    auto ex20base = mat<3,3>(cst(0.),cst(1),cst(-1.),
                             cst(1.),cst(1.),cst(0.),
                             cst(-1.),cst(0.),cst(1.) );
    auto ex20 = eig( ex20base );
    auto eval20 = ex20.evaluate();
    BOOST_CHECK( eval20.rows() == 3 && eval20.cols() == 1 );
    BOOST_CHECK_CLOSE( eval20(0,0), -1.0, tolerance );
    BOOST_CHECK_CLOSE( eval20(1,0), 1.0, tolerance );
    BOOST_CHECK_CLOSE( eval20(2,0), 2.0, tolerance );

    auto VhP0m = Pchm<0>( mesh );
    auto uP0m = VhP0m->element();
    uP0m.on(_range=elements(mesh),_expr=mat<2,2>(cst(1.),cst(2.),cst(3.),cst(4.)));
    auto eval21 = idv(uP0m).evaluate();
    BOOST_CHECK( eval21.rows() == 2 && eval21.cols() == 2 );
    BOOST_CHECK_CLOSE( eval21(0,0), 1.0, tolerance );
    BOOST_CHECK_CLOSE( eval21(0,1), 2.0, tolerance );
    BOOST_CHECK_CLOSE( eval21(1,0), 3.0, tolerance );
    BOOST_CHECK_CLOSE( eval21(1,1), 4.0, tolerance );
}
BOOST_AUTO_TEST_SUITE_END()
