/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-10-27

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
   \file test_mortar.cpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-10-27
 */
#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE Mortar testsuite


#include <testsuite/testsuite.hpp>

#include <boost/timer.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <feel/feelcore/environment.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelvf/vf.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( mortar )

BOOST_AUTO_TEST_CASE( test_no_mortar_1 )
{
    BOOST_TEST_MESSAGE( "test_no_mortar_1" );
    using namespace Feel;
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Xh = Pch<1>(mesh);
    BOOST_CHECK_MESSAGE( Xh->is_mortar == false, "Space should not be mortar" ) ;
}
BOOST_AUTO_TEST_CASE( test_no_mortar_2 )
{
    BOOST_TEST_MESSAGE( "test_no_mortar_2" );
    using namespace Feel;
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Xh = THch<1>(mesh);
    BOOST_CHECK_MESSAGE( Xh->functionSpace<0>()->is_mortar == false, "Space 0 should not be mortar" ) ;
    BOOST_CHECK_MESSAGE( Xh->functionSpace<1>()->is_mortar == false, "Space 1 should not be mortar" ) ;

    BOOST_TEST_MESSAGE( "test_no_mortar_2 done" );
}
BOOST_AUTO_TEST_CASE( test_mortar_1 )
{
    using namespace Feel;
    BOOST_TEST_MESSAGE( "test_mortar_1" );
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Xh = Moch<1>(mesh);
    BOOST_CHECK_MESSAGE(Xh->is_mortar == true, "Space should be mortar" ) ;
    BOOST_TEST_MESSAGE( "test_mortar_1 done" );
}


BOOST_AUTO_TEST_SUITE_END()
