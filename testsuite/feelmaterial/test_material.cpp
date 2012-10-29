/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-05-08

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file test_material.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-05-08
 */

#define BOOST_TEST_MODULE material tests
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;


#include <feel/feelmaterial/materiallib.hpp>
using namespace Feel;


BOOST_AUTO_TEST_SUITE( material_testsuite )

BOOST_AUTO_TEST_CASE( test_air )
{
    MaterialLib lib;
    material_ptrtype m = MaterialLib::material( "Air" );
    BOOST_CHECK_CLOSE( m->rho(), 1.2, 1e-15 );
    BOOST_CHECK_CLOSE( m->nu(), 1.7*1e-5, 1e-15 );
}

BOOST_AUTO_TEST_SUITE_END() // material test suite


