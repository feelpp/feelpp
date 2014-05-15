/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-01-08

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
   \file test_lambda.cpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-01-08
 */


#define USE_BOOST_TEST 1

// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE function space testsuite
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN

#include <testsuite/testsuite.hpp>

#include <boost/timer.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <feel/feelcore/environment.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/pch.hpp>

#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelvf/vf.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( lambda )

BOOST_AUTO_TEST_CASE( test_lambda_cst1 )
{
    using namespace Feel;
    auto I = _e1;
    typedef decltype(I( cst(1.) )) t1;
    typedef decltype(cst(1.)) t2;
    BOOST_MPL_ASSERT_MSG( (boost::is_same<t1,t2>::value), INVALID_TYPE,(decltype(I( cst(1.) )),decltype(cst(1.))));
}
BOOST_AUTO_TEST_CASE( test_lambda_cst2 )
{
    BOOST_TEST_MESSAGE( "test_lambda_cst2" );
    using namespace Feel;
    auto I1 = _e1;
    typedef decltype(I1( cst(1.) )) t1;
    typedef decltype(cst(1.)) t2;
    BOOST_MPL_ASSERT_MSG( (boost::is_same<t1,t2>::value), INVALID_TYPE,(decltype(I1( cst(1.) )),decltype(cst(1.))));
    t1 a(cst(1.0));
    t2 b(cst(1.0));
    BOOST_CHECK_CLOSE( a.expression().value(), 1., 1e-13 );
    BOOST_CHECK_CLOSE( b.expression().value(), 1., 1e-13 );
    BOOST_CHECK_CLOSE( a.expression().value(), b.expression().value(), 1e-13 );
    BOOST_TEST_MESSAGE( "test_lambda_cst2 done" );
}
BOOST_AUTO_TEST_CASE( test_lambda_cst3 )
{
    using namespace Feel;
    auto I = _e1*_e2;
    typedef decltype(I( cst(1.), cst(.5) )) t1;
    typedef decltype(cst(1.)*cst(.5)) t2;
    BOOST_MPL_ASSERT_MSG( (boost::is_same<t1,t2>::value), INVALID_TYPE,(decltype(I( cst(1.),cst(.5) )),decltype(cst(1.)*cst(.5))));
}

BOOST_AUTO_TEST_CASE( test_lambda_cos )
{
    BOOST_TEST_MESSAGE( "test_lambda_cos" );
    using namespace Feel;
    auto I1 = cos(_e1);
    typedef decltype(I1( cst(1.) )) t1;
    typedef decltype(cos(cst(1.))) t2;
    BOOST_MPL_ASSERT_MSG( (boost::is_same<t1,t2>::value), INVALID_TYPE,(decltype(I1( cst(1.) )),decltype(cos(cst(1.)))));
    BOOST_TEST_MESSAGE( "test_lambda_cos done" );
}
BOOST_AUTO_TEST_CASE( test_lambda_trans )
{
    BOOST_TEST_MESSAGE( "test_lambda_trans" );
    using namespace Feel;
    auto I1 = trans(_e1);
    typedef decltype(I1( cos(cst(1.)))) t1;
    typedef decltype(trans(cos(cst(1.)))) t2;
    BOOST_MPL_ASSERT_MSG( (boost::is_same<t1,t2>::value), INVALID_TYPE,(decltype(I1( cos(cst(1.)))),decltype(trans(cos(cst(1.))))));
    BOOST_TEST_MESSAGE( "test_lambda_trans done" );
}

BOOST_AUTO_TEST_CASE( test_lambda_int_cst )
{
    BOOST_TEST_MESSAGE( "test_lambda_int_cst" );
    using namespace Feel;
    auto mesh = unitSquare();
    BOOST_TEST_MESSAGE( "test_lambda_int_cst mesh generated" );
    auto I = integrate( elements(mesh), _expr=_e1, _verbose=true );
    BOOST_TEST_MESSAGE( "test_lambda_int_cst integral defined" );

    auto I1 = integrate( elements(mesh), cst(1.) );
    auto Xh = Pch<1>( mesh );
    auto u = project( _space=Xh, _range=elements(mesh), _expr=cst(1.) );

    BOOST_CHECK_CLOSE( I1.evaluate()( 0, 0 ), 1, 1e-10 );
    BOOST_CHECK_CLOSE( I( cst(1.0) ).expression().expression().expression().value(), 1, 1e-10 );
    BOOST_CHECK_CLOSE( I( cst(1.0) ).evaluate()( 0, 0 ), 1, 1e-10 );
    BOOST_CHECK_CLOSE( I( idv(u) ).evaluate()( 0, 0 ), 1, 1e-10 );

    BOOST_TEST_MESSAGE( "test_lambda_int_cst done" );
}
BOOST_AUTO_TEST_CASE( test_lambda_int_cst2 )
{
    BOOST_TEST_MESSAGE( "test_lambda_int_cst2" );
    using namespace Feel;
    auto mesh = unitSquare();
    BOOST_TEST_MESSAGE( "test_lambda_int_cst2 mesh generated" );
    auto I = integrate( elements(mesh), _expr=_e1*_e2, _verbose=true );
    auto I3 = integrate( elements(mesh), _expr=_e1*_e2*_e3, _verbose=true );
    BOOST_TEST_MESSAGE( "test_lambda_int_cst2 integral defined" );

    auto I1 = integrate( elements(mesh), cst(1.) );
    auto Xh = Pch<1>( mesh );
    auto u = project( _space=Xh, _range=elements(mesh), _expr=cst(.5) );
    auto w = project( _space=Xh, _range=elements(mesh), _expr=Px()+1 );
    auto v = project( _space=Xh, _range=elements(mesh), _expr=Px()+Py() );

    BOOST_CHECK_CLOSE( I1.evaluate()( 0, 0 ), 1, 1e-10 );
    BOOST_CHECK_CLOSE( I( cst(.5), cst(2.) ).evaluate()( 0, 0 ), 1, 1e-10 );
    BOOST_CHECK_CLOSE( I( idv(u), gradv(v)(0,0)*2.).evaluate()( 0, 0 ), 1, 1e-10 );
    BOOST_CHECK_CLOSE( I3( cst(.5), cst(2.), cst(1.0) ).evaluate()( 0, 0 ), 1, 1e-10 );
    BOOST_CHECK_CLOSE( I3( cst(2.)*idv(u), gradv(v)*trans(gradv(v)), idv(w)/(cst(2.0)*(Px()+1)) ).evaluate()( 0, 0 ), 1, 1e-10 );

    BOOST_TEST_MESSAGE( "test_lambda_int_cst2 done" );
}

BOOST_AUTO_TEST_CASE( test_lambda_int )
{
    BOOST_TEST_MESSAGE( "test_lambda_int" );
    using namespace Feel;
    auto mesh = unitSquare();
    BOOST_TEST_MESSAGE( "test_lambda_int mesh generated" );
    auto I = integrate( elements(mesh), _expr=_e1, _verbose=true );
    BOOST_TEST_MESSAGE( "test_lambda_int integral defined" );

    auto I1 = integrate( elements(mesh), Px()*Px()+Py()*Py() );
    auto Xh = Pch<2>( mesh );
    auto u = project( _space=Xh, _range=elements(mesh), _expr=Px()*Px()+Py()*Py() );

    BOOST_CHECK_CLOSE( I1.evaluate()( 0, 0 ), 2./3., 1e-10 );
    BOOST_CHECK_CLOSE( I( Px()*Px()+Py()*Py() ).evaluate()( 0, 0 ), 2./3., 1e-10 );
    BOOST_CHECK_CLOSE( I( idv(u) ).evaluate()( 0, 0 ), 2./3., 1e-10 );

    BOOST_TEST_MESSAGE( "test_lambda_int done" );
}

BOOST_AUTO_TEST_CASE( test_lambda_op )
{
    BOOST_TEST_MESSAGE( "test_lambda_op" );
    using namespace Feel;
    auto mesh = unitSquare();
    BOOST_TEST_MESSAGE( "test_lambda_op mesh generated" );
    auto I = integrate( elements(mesh), _expr=0.1*_e1, _verbose=true );
    BOOST_TEST_MESSAGE( "test_lambda_op integral defined" );

    auto I1 = integrate( elements(mesh), Px()*Px()+Py()*Py() );
    auto Xh = Pch<2>( mesh );
    auto u = project( _space=Xh, _range=elements(mesh), _expr=Px()*Px()+Py()*Py() );

    BOOST_CHECK_CLOSE( I1.evaluate()( 0, 0 ), 2./3., 1e-10 );
    BOOST_CHECK_CLOSE( I( cst(10.) ).evaluate()( 0, 0 ), 1., 1e-10 );
    BOOST_CHECK_CLOSE( I( idv(u) ).evaluate()( 0, 0 ), 0.1*2./3., 1e-10 );
    BOOST_CHECK_CLOSE( I( Px()*Px()+Py()*Py() ).evaluate()( 0, 0 ), 0.1*2./3., 1e-10 );

    BOOST_TEST_MESSAGE( "test_lambda_op done" );
}

BOOST_AUTO_TEST_CASE( test_lambda_2 )
{
    BOOST_TEST_MESSAGE( "test_lambda_2" );
    using namespace Feel;
    auto mesh = unitSquare();
    BOOST_TEST_MESSAGE( "test_lambda_2 mesh generated" );
    auto I = integrate( elements(mesh), _expr=_e1, _verbose=true );
    BOOST_TEST_MESSAGE( "test_lambda_2 integral defined" );

    auto Xh = Pch<2>( mesh );
    auto u = project( _space=Xh, _range=elements(mesh), _expr=Px() );
    auto I1 = integrate( elements(mesh), idv(u) * idv(u) );

    BOOST_CHECK_CLOSE( I1.evaluate()( 0, 0 ), I( idv(u) * idv(u) ).evaluate()( 0, 0 ) , 1e-10 );
    BOOST_CHECK_CLOSE( I( Px()*Px() ).evaluate()( 0, 0 ), 1./3., 1e-10 );
    BOOST_CHECK_CLOSE( I( idv(u)*idv(u) ).evaluate()( 0, 0 ), 1./3., 1e-10 );

    auto lambda_expr = _e1;
    double integrate_lambda = integrate( _range=elements(mesh), _expr=lambda_expr( idv(u)*idv(u) ) ).evaluate()(0,0);
    BOOST_CHECK_CLOSE( integrate_lambda , 1./3., 1e-10 );
    //use twice _e1, once with idv() and one with gradv() or trans( gradv() )
}

BOOST_AUTO_TEST_CASE( test_lambda_div )
{
    BOOST_TEST_MESSAGE( "test_lambda_div" );
    using namespace Feel;
    auto mesh = unitSquare();
    BOOST_TEST_MESSAGE( "test_lambda_div mesh generated" );
    auto I = integrate( _range=elements(mesh), _expr=cst(1.)/(1+_e1), _verbose=true );
    BOOST_TEST_MESSAGE( "test_lambda_div integral defined" );

    auto I1 = integrate( elements(mesh), cst(1.)/(1.+Px()) );
    auto Xh = Pch<1>( mesh );
    auto u = project( _space=Xh, _range=elements(mesh), _expr=Px() );

    //compilation failed
    BOOST_CHECK_CLOSE( I1.evaluate()( 0, 0 ), I( idv(u)  ).evaluate()( 0, 0 ) , 1e-10 );
    BOOST_CHECK_CLOSE( I( Px() ).evaluate()( 0, 0 ), I( idv(u) ).evaluate()( 0, 0 ) , 1e-10 );

    BOOST_TEST_MESSAGE( "test_lambda_div done" );

    double alpha=2;
    auto T = project( _space=Xh, _range=elements(mesh), _expr=cos( Px() ) );
    double T0=3;
    double sigma=0.5;
    auto expr_lambda = cst(sigma)/( cst(1)+alpha*( _e1-T0 ) );
    auto I_sigma = integrate( _range=elements(mesh), _expr=cst(sigma)/(cst(1)+alpha*(idv(T)-T0)) );
    auto I_lambda_sigma = integrate( _range=elements(mesh), _expr=expr_lambda );
    //compilation failed
    double integrate_lambda = integrate( _range=elements(mesh), _expr=expr_lambda( idv(T) ) ).evaluate()(0,0);
    BOOST_CHECK_CLOSE( I_sigma.evaluate()( 0, 0 ), I_lambda_sigma(idv(T)).evaluate()( 0, 0 ) , 1e-10 );
    BOOST_CHECK_CLOSE( integrate_lambda, I_lambda_sigma(idv(T)).evaluate()( 0, 0 ) , 1e-10 );

}


BOOST_AUTO_TEST_SUITE_END()
