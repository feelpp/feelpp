/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2022-02-20

  Copyright (C) 2022 Université de Strasbourg

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
   \file test_json.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2022-02-20
 */
#define USE_BOOST_TEST 1
// Boost.Test

#define BOOST_TEST_MODULE json testsuite
#include <feel/feelcore/testsuite.hpp>

#include <boost/algorithm/string.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/json.hpp>
inline Feel::po::options_description
makeOptions()
{
    Feel::po::options_description integrationoptions( "Test Environment options" );
    integrationoptions.add_options()( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "h value" )( "fname", Feel::po::value<std::string>()->default_value( "toto" ), "file name" );
    return integrationoptions.add( Feel::feel_options() );
}

inline Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "test_json",
                           "test_json",
                           "0.2",
                           "Environment class tests",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2013-2022 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( git )

BOOST_AUTO_TEST_CASE( test_json1 )
{
    using namespace Feel;
    nl::json j;
    j["address"]["street"] = "Fake Street";
    j["address"]["housenumber"] = "12";

    try
    {
        int housenumber = j["address"]["housenumber"];
    }
    catch ( nl::json::exception& e )
    {
        std::string estr{ e.what() };
        bool found_extra_info = boost::algorithm::find_first(estr, "/address/housenumber" );
        BOOST_CHECK( !JSON_DIAGNOSTICS || found_extra_info );
        handleExceptions();
    }
}

BOOST_AUTO_TEST_SUITE_END()
