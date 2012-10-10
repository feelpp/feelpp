/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-08-07

  Copyright (C) 2010 Universite Joseph Fourier (Grenoble I)

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
   \file test_vector.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-08-07
 */
#include <cmath>

#include <boost/timer.hpp>
#define BOOST_TEST_MODULE vector testsuite
// disable the main function creation, use our own
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


using boost::unit_test::test_suite;


#include <feel/feelcore/traits.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelalg/vectorublas.hpp>

BOOST_AUTO_TEST_SUITE( vector )

BOOST_AUTO_TEST_CASE( test1 )
{
    using namespace Feel;
    Feel::Environment env( _argc=boost::unit_test::framework::master_test_suite().argc,
                           _argv=boost::unit_test::framework::master_test_suite().argv );

    using namespace Feel;
    VectorUblas<double> v1( 100 ), v2( 100 ), v3( 100 );
    v1.setConstant( 1 );
    v2.setConstant( 2 );
    BOOST_CHECK_CLOSE( v1.sqrt().sum(), v1.size(), 1e-10 );
    BOOST_CHECK_CLOSE( v2.sqrt().sum(), sqrt( 2 )*v1.size(), 1e-10 );
    BOOST_CHECK_CLOSE( v2.pow( 2 ).sqrt().sum(), 2*v1.size(), 1e-10 );
    v3.setZero();
    v3 = element_product( v1, v2 );
    BOOST_CHECK_CLOSE( v3.sqrt().sum(), sqrt( 2*1 )*v1.size(), 1e-10 );
    BOOST_CHECK_CLOSE( v3.sum(), 2*1*v1.size(), 1e-10 );



}
BOOST_AUTO_TEST_SUITE_END()
