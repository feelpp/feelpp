/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-05-20

  Copyright (C) 2014-2016 Feel++ Consortium

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

#define BOOST_TEST_MODULE continuity testsuite
#include <testsuite.hpp>

#include <boost/mpl/list.hpp>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( continuitysuite )
    typedef boost::mpl::list<boost::mpl::int_<1>,
                             boost::mpl::int_<2>,
                             boost::mpl::int_<3>
                             > dim_types;
BOOST_AUTO_TEST_CASE_TEMPLATE( LagrangeCG, T, dim_types )
{
    BOOST_TEST_MESSAGE( "check continuity for LagrangeCG in  " << T::value << "D\n" );
    typedef Mesh<Simplex<T::value,1> > mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    mesh_ptrtype mesh = loadMesh( _mesh=new mesh_type );

    auto Xh = Pch<1>( mesh, true );
    auto u = Xh->element();
    auto a1 = form1( _test=Xh );
    a1  = integrate( _range=internalfaces( mesh ), _expr=(leftface(id(u))+rightface(-id(u)) ) );
    a1.vectorPtr()->close();
    if ( mesh_type::nDim == 1 )
        u.on( _range=elements(mesh), _expr=expr("x:x") );
    else if ( mesh_type::nDim == 2 )
        u.on( _range=elements(mesh), _expr=expr("x*y:x:y") );
    else
        u.on( _range=elements(mesh), _expr=expr("x*y*z:x:y:z") );
    if ( Environment::numberOfProcessors() == 1 )
    {
        a1.vector().printMatlab("LagrangeCG.m");
        u.printMatlab("uCG.m");
    }
    BOOST_CHECK_SMALL( a1( u ), 1e-10 );

    BOOST_TEST_MESSAGE( "LagrangeCG, a1(u)=" << a1(u)  );
    BOOST_TEST_MESSAGE( "check continuity for LagrangeCG in  " << T::value << "D\n" );
}

BOOST_AUTO_TEST_SUITE_END()
