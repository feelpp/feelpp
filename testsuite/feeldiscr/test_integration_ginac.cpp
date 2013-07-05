/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-07-05

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
   \file test_integration_ginac.cpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-07-05
 */
#include <sstream>
#include <boost/timer.hpp>

#define BOOST_TEST_MODULE ginac integration testsuite
#include <testsuite/testsuite.hpp>

#include <boost/mpl/list.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( ginacsuite )
//typedef boost::mpl::list<boost::mpl::int_<1>,boost::mpl::int_<2>,boost::mpl::int_<3> > dim_types;
typedef boost::mpl::list<boost::mpl::int_<2> > dim_types;
//typedef boost::mpl::list<boost::mpl::int_<2>,boost::mpl::int_<3>,boost::mpl::int_<1> > dim_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( ginacint, T, dim_types )
{
    BOOST_TEST_MESSAGE( "check 1 ginac int  = " << T::value << "\n" );
    typedef Mesh<Simplex<T::value,1> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name=( boost::format( "hypercube-%1%" )  % T::value ).str() ,
                                                _usenames=true,
                                                _addmidpoint=false,
                                                _shape="hypercube",
                                                _h=0.2 ),
                                        _update=MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES );

    auto P1h = Pch<1>( mesh );
    //auto vars = symbols<T::value>();
    symbol x("x"), y ("y");

    //auto f1g = vars[0];auto f2g = vars[1];
    ex f1g = sin(x)/y ;
    ex f2g = x*cos(y);
    std::vector<symbol> s = {x,y};
    auto f1= expr( f1g, s, "a" );
    auto f2 = expr(f2g, s, "b");
    auto xg = integrate(_range=elements(mesh), _expr=f1 ).evaluate()(0,0);
    auto yg = integrate(_range=elements(mesh), _expr=f2 ).evaluate()(0,0);

    BOOST_CHECK_CLOSE( xg, 0.5, 1e-10 );
    BOOST_CHECK_CLOSE( yg, 0.5, 1e-10 );
    BOOST_TEST_MESSAGE( "test gravity center " << xg << "," << yg );

}

BOOST_AUTO_TEST_SUITE_END()

#if 0
int BOOST_TEST_CALL_DECL
main( int argc, char* argv[] )
{
    Feel::Environment env( argc, argv );
    int ret = ::boost::unit_test::unit_test_main( &init_unit_test, argc, argv );

    return ret;
}

#endif
