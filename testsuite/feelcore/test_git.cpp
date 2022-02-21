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
   \file test_git.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2022-02-20
 */
#define USE_BOOST_TEST 1
// Boost.Test

#define BOOST_TEST_MODULE git testsuite
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/git.hpp>

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
    Feel::AboutData about( "test_git",
                           "test_git",
                           "0.2",
                           "Environment class tests",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2013-2016 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( git )

BOOST_AUTO_TEST_CASE( test_git1 )
{
    using namespace Feel;

    if ( GitMetadata::populated() )
    {
        if ( GitMetadata::anyUncommittedChanges() )
        {
            std::cerr << "WARN: there were uncommitted changes at build-time." << std::endl;
        }
        std::cout << "commit " << GitMetadata::commitSHA1() << " (HEAD)\n"
                  << "describe " << GitMetadata::describe() << "\n"
                  << "Author: " << GitMetadata::authorName() << " <" << GitMetadata::authorEmail() << ">\n"
                  << "Date: " << GitMetadata::commitDate() << "\n\n"
                  << GitMetadata::commitSubject() << "\n"
                  << GitMetadata::commitBody() << std::endl;
    }
    else
    {
        std::cerr << "WARN: failed to get the current git state. Is this a git repo?" << std::endl;
    }
}

BOOST_AUTO_TEST_SUITE_END()
