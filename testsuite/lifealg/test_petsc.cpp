/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2004-10-27
 */
#include <cstdlib>
#include <cassert>

#include <life/lifecore/life.hpp>
#include <life/lifecore/traits.hpp>


#include <lifeconfig.h>

#if defined(HAVE_BOOST_TEST)
// Boost.Test
#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/bind.hpp>

using boost::unit_test_framework::test_suite;


#if defined(HAVE_PETSC_H)
#include <life/lifecore/application.hpp>
#include <life/lifealg/solverlinearpetsc.hpp>

Life::AboutData
makeAbout()
{
    Life::AboutData about( "test_petsc" ,
                            "test_petsc" ,
                            "0.1",
                            "test_petsc",
                            Life::AboutData::License_GPL,
                            "Copyright (c) 2007 Université Joseph Fourier Grenoble 1");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
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

    Life::Application app(argc,argv,makeAbout());

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
#endif /* HAVE_PETSC_H */
#else
int main()
{
    return EXIT_SUCCESS;
}
#endif
