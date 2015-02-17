/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-02-20

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
   \file test_env.cpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-02-20
 */
#define USE_BOOST_TEST 1
// Boost.Test

// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE submesh testsuite
#include <testsuite/testsuite.hpp>

#include <feel/feelcore/environment.hpp>

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description integrationoptions( "Test Environment options" );
    integrationoptions.add_options()
        ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "h value" )
        ( "fname", Feel::po::value<std::string>()->default_value( "toto" ), "file name" )
    ;
    return integrationoptions.add( Feel::feel_options() );
}

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "test_env" ,
                           "test_env" ,
                           "0.2",
                           "Environment class tests",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2013-2015 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( submesh )

BOOST_AUTO_TEST_CASE( test_env1 )
{
    using namespace Feel;
    BOOST_TEST_MESSAGE( "test_env1" );

    double hsize = Feel::option(_name="hsize").as<double>(), hsize2;
    BOOST_CHECK_CLOSE( Feel::optionT(_name="hsize", _opt=hsize2 ), hsize, 1e-15 );
    BOOST_CHECK_CLOSE( hsize2, hsize, 1e-15 );

    std::string fname = Feel::option(_name="fname").as<std::string>(), fname2;
    double t;
    BOOST_CHECK_EQUAL( Feel::optionT(_name="fname", _opt=fname2 ), fname );
    BOOST_CHECK_EQUAL( fname2, fname );

    // if you try to do the following an exception is triggerred and a nice
    // message with a backtrace is printed
    // Feel::optionT(_name="fname", _opt=t );

    BOOST_TEST_MESSAGE( "test_env1 done" );
}

BOOST_AUTO_TEST_SUITE_END()

