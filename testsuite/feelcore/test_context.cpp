/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-04-07

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file test_context.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-04-07
 */
#define BOOST_TEST_MODULE context test
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/context.hpp>
#include <feel/feelalg/enums.hpp>



BOOST_AUTO_TEST_CASE( test_context )
{
    using namespace Feel;
    Context ctx( int(OnContext::ELIMINATION) | int(OnContext::ELIMINATION_KEEP_DIAGONAL) );
    BOOST_CHECK_EQUAL( ctx.test( OnContext::PENALISATION ), false );
    BOOST_CHECK_EQUAL( ctx.test( OnContext::ELIMINATION ), true );
    BOOST_CHECK_EQUAL( ctx.test( OnContext::ELIMINATION_KEEP_DIAGONAL ), true );
    BOOST_CHECK_EQUAL( ctx.test( OnContext::ELIMINATION_SYMMETRIC ), false );


    ctx = int(OnContext::ELIMINATION) | int(OnContext::ELIMINATION_SYMMETRIC);
    BOOST_CHECK_EQUAL( ctx.test( OnContext::PENALISATION ), false );
    BOOST_CHECK_EQUAL( ctx.test( OnContext::ELIMINATION ), true );
    BOOST_CHECK_EQUAL( ctx.test( OnContext::ELIMINATION_KEEP_DIAGONAL ), false );
    BOOST_CHECK_EQUAL( ctx.test( OnContext::ELIMINATION_SYMMETRIC ), true );
    BOOST_CHECK_EQUAL( ctx.test( int(OnContext::ELIMINATION_SYMMETRIC) | int(OnContext::ELIMINATION) ), true );

#if 0
    detail::Flags<on_context_type> oct;

    //FEELPP_DECLARE_FLAGS( onct, on_context_type );


    onct |= ( OnContext::ELIMINATION|OnContext::ELIMINATION_SYMMETRIC );
    BOOST_CHECK_EQUAL( onct.testFlag( OnContext::ELIMINATION ), true );
    BOOST_CHECK_EQUAL( onct.testFlag( OnContext::ELIMINATION_KEEP_DIAGONAL ), false );
    BOOST_CHECK_EQUAL( onct.testFlag( OnContext::ELIMINATION_SYMMETRIC ), true );
#endif
}

BOOST_AUTO_TEST_CASE( test_context_matrix_properties )
{
    using namespace Feel;
    Context ctx( HERMITIAN | POSITIVE_DEFINITE );
    BOOST_CHECK_EQUAL( ctx.test( HERMITIAN | POSITIVE_DEFINITE ), true );
    BOOST_CHECK_EQUAL( ctx.test( SINGULAR ), false );

    ctx = HERMITIAN|SINGULAR;
    BOOST_CHECK_EQUAL( ctx.test( POSITIVE_DEFINITE ), false );
    BOOST_CHECK_EQUAL( ctx.test( SINGULAR ), true );
    BOOST_CHECK_EQUAL( ctx.test( HERMITIAN ), true );

}
