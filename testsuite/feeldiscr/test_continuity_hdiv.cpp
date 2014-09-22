/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-05-20

  Copyright (C) 2014 Feel++ Consortium

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
   \file test_continuity.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-05-20
 */
#include <sstream>
#include <boost/timer.hpp>

#define BOOST_TEST_MODULE continuity testsuite
#include <testsuite/testsuite.hpp>

#include <boost/mpl/list.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/dh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/ginac.hpp>
#include <feel/feelvf/print.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( continuitysuite )
//typedef boost::mpl::list<boost::mpl::int_<1>,boost::mpl::int_<2>,boost::mpl::int_<3> > dim_types;
    typedef boost::mpl::list<boost::mpl::pair<boost::mpl::int_<2>,boost::mpl::int_<0> >,
                             // boost::mpl::pair<boost::mpl::int_<2>,boost::mpl::int_<1> >,
                             // boost::mpl::pair<boost::mpl::int_<2>,boost::mpl::int_<2> >,
                             boost::mpl::pair<boost::mpl::int_<3>,boost::mpl::int_<0> >
                             > dim_types;
//typedef boost::mpl::list<boost::mpl::int_<2>,boost::mpl::int_<3>,boost::mpl::int_<1> > dim_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( HDivRT0, T, dim_types )
{
    BOOST_TEST_MESSAGE( "check continuity for HDivRT in  " << T::first::value << "D P" << T::second::value << "\n" );
    Feel::Environment::changeRepository( boost::format( "/testsuite/feeldiscr/%1%/%2%D/P%3%/" )
                                         % Feel::Environment::about().appName()
                                         % T::first::value % T::second::value );


    typedef Mesh<Simplex<T::first::value,1> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    mesh_ptrtype mesh = loadMesh( _mesh=new mesh_type );

    auto Xh = Dh<T::second::value>( mesh );
    auto u = Xh->element();
    auto a1 = form1( _test=Xh );
    a1  = integrate( internalfaces( mesh ), (leftface(trans(id(u))*N())+rightface(-trans(id(u))*N())) );
    u.on(  _range=elements(mesh), _expr=expr<T::first::value,1>("{x*y,x+y}:x:y") );
    a1.vector().printMatlab("HDivRT.m");
    BOOST_CHECK_SMALL( a1( u ), 1e-10 );
    u.printMatlab("uRT0.m");

    BOOST_TEST_MESSAGE( "HDivRT, a1(u)=" << a1(u)  );
    BOOST_TEST_MESSAGE( "check continuity for HDivRT in  " << T::first::value << "D P" << T::second::value << " done\n" );
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
