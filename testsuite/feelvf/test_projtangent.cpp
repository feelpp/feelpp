/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-01-14

  Copyright (C) 2014 Feel++ Consortium

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
#if defined(USE_BOOST_TEST)
#define BOOST_TEST_MODULE test_projtangent
#include <testsuite/testsuite.hpp>
#endif

#include <feel/feelalg/backend.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/unitcircle.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelfilters/unitcube.hpp>
#include <feel/feelvf/vf.hpp>

/** use Feel namespace */
using namespace Feel;

inline
AboutData
makeAbout()
{
    AboutData about( "test_projtangent" ,
                     "test_projtangent" ,
                     "0.2",
                     "nD(n=2,3) test projtangent",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2014 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

#if defined(USE_BOOST_TEST)
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), feel_options() );
BOOST_AUTO_TEST_SUITE( projtangent_suite )

BOOST_AUTO_TEST_CASE( test_0 )
{
    BOOST_TEST_MESSAGE( "test_circle" );
    Feel::Environment::changeRepository( boost::format( "/testsuite/feelvf/%1%/test_circle/h_%2%/" )
                                         % Feel::Environment::about().appName()
                                         % option(_name="gmsh.hsize").as<double>() );
    auto mesh = unitCircle();
    double a = integrate( _range= boundaryfaces( mesh ), _expr= pT()*N() ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( a, 1e-13);
    BOOST_TEST_MESSAGE( "pT*N=" << a );

    double b = integrate( _range= boundaryfaces( mesh ), _expr= inner(vec(-Py(),Px()),vec(-Py(),Px()) ) ).evaluate()( 0,0 );
    BOOST_TEST_MESSAGE( "(-y,x)=" << b );
    BOOST_CHECK_CLOSE( b, 8*atan(1.), 1);

    double c = integrate( _range= boundaryfaces( mesh ), _expr= norm2(pT()*vec(-Py(),Px()) ) ).evaluate()( 0,0 );
    BOOST_TEST_MESSAGE( "pT(-y,x)=" << c );
    BOOST_CHECK_CLOSE( c, 8*atan(1.), 1);

    BOOST_TEST_MESSAGE( "test_circle done" );
}

BOOST_AUTO_TEST_CASE( test_1 )
{
    BOOST_TEST_MESSAGE( "test_square" );
    Feel::Environment::changeRepository( boost::format( "/testsuite/feelvf/%1%/test_square/h_%2%/" )
                                         % Feel::Environment::about().appName()
                                         % option(_name="gmsh.hsize").as<double>() );
    auto mesh = unitSquare();
    double a = integrate( _range= boundaryfaces( mesh ), _expr= norm2(pT()*N()) ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( a, 1e-13);
    BOOST_TEST_MESSAGE( "pT*N=" << a );

    double b = integrate( _range= boundaryfaces( mesh ), _expr= norm2(pT()*vec(-Ny(),Nx()) )).evaluate()( 0,0 );
    BOOST_CHECK_CLOSE( b, 4, 1e-13 );
    BOOST_TEST_MESSAGE( "pT*(-Ny,Nx)=" << b << " should be equal to perimeter of unit square : 4" );

    boost::mpi::timer ti;
    integrate( _range= boundaryfaces( mesh ), _expr= norm2(pT()*N()) ).evaluate()( 0,0 );
    BOOST_TEST_MESSAGE( "time project tangent : " << ti.elapsed() );
    ti.restart();
    integrate( _range= boundaryfaces( mesh ), _expr= norm2((eye<2,2>()-N()*trans(N()))*N()) ).evaluate()( 0,0 );
    BOOST_TEST_MESSAGE( "time project tangent by hand : " << ti.elapsed() );

    BOOST_TEST_MESSAGE( "test_square done." );

}

BOOST_AUTO_TEST_CASE( test_2 )
{
    BOOST_TEST_MESSAGE( "test_cube" );
    Feel::Environment::changeRepository( boost::format( "/testsuite/feelvf/%1%/test_cube/h_%2%/" )
                                         % Feel::Environment::about().appName()
                                         % option(_name="gmsh.hsize").as<double>() );
    auto mesh = unitCube();
    double a = integrate( _range= boundaryfaces( mesh ), _expr= pT()*N() ).evaluate()( 0,0 );
    BOOST_CHECK_SMALL( a, 1e-13);
    BOOST_TEST_MESSAGE( "pT*N=" << a );

    boost::mpi::timer ti;
    integrate( _range= boundaryfaces( mesh ), _expr= norm2(pT()*N()) ).evaluate()( 0,0 );
    BOOST_TEST_MESSAGE( "time project tangent : " << ti.elapsed() );
    ti.restart();
    integrate( _range= boundaryfaces( mesh ), _expr= norm2((eye<3,3>()-N()*trans(N()))*N()) ).evaluate()( 0,0 );
    BOOST_TEST_MESSAGE( "time project tangent by hand : " << ti.elapsed() );

    BOOST_TEST_MESSAGE( "test_cube done." );

}

BOOST_AUTO_TEST_SUITE_END()
#else
int main(int argc, char** argv )
{
    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=feel_options() );

    Feel::Environment::changeRepository( boost::format( "/testsuite/feelvf/%1%/test_cube/h_%2%/" )
                                         % Feel::Environment::about().appName()
                                         % doption(_name="gmsh.hsize") );
    auto mesh = unitCube();
    double a = integrate( _range= boundaryfaces( mesh ), _expr= norm2(pT()*N()) ).evaluate()( 0,0 );
    LOG_IF(WARNING, math::abs(a) < 1e-13 ) << "should be 0 and the value is " << a;

    boost::mpi::timer ti;
    integrate( _range= boundaryfaces( mesh ), _expr= norm2(pT()*N()) ).evaluate()( 0,0 );
    LOG(INFO) << "time project tangent : " << ti.elapsed();
    ti.restart();
    integrate( _range= boundaryfaces( mesh ), _expr= norm2((eye<3,3>()-N()*trans(N()))*N()) ).evaluate()( 0,0 );
    LOG(INFO) << "time project tangent by hand : " << ti.elapsed();
}
#endif
