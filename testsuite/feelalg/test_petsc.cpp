/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2004-10-27

  Copyright (C) 2004 EPFL

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
   \file test_petsc.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2004-10-27
 */
#include <cstdlib>
#include <cassert>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/traits.hpp>


#include <feel/feelconfig.h>

#if defined(FEELPP_HAS_BOOST_TEST)
// Boost.Test
#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/bind.hpp>

using boost::unit_test_framework::test_suite;


#if defined(FEELPP_HAS_PETSC_H)
#include <feel/feelcore/application.hpp>
#include <feel/feelalg/solverlinearpetsc.hpp>

Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "test_petsc" ,
                           "test_petsc" ,
                           "0.1",
                           "test_petsc",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2007 Université Joseph Fourier Grenoble 1" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

/**
 * Test non linear solver framework
 */
void petsc_snes()
{

}


test_suite*
init_unit_test_suite( int argc, char** argv )
{
    test_suite* test= BOOST_TEST_SUITE( "PETSC Unit Test" );

    Feel::Application app( argc,argv,makeAbout() );

    // this example will pass cause we know ahead of time number of expected failures
    test->add( BOOST_TEST_CASE( &petsc_snes ), 0 );

    return test;
}
#else
test_suite*
init_unit_test_suite( int /*argc*/, char** /*argv*/ )
{
    test_suite* test= BOOST_TEST_SUITE( "PETSC Unit Test" );
    return test;
}
#endif /* FEELPP_HAS_PETSC_H */
#else
int main()
{
    return EXIT_SUCCESS;
}
#endif
