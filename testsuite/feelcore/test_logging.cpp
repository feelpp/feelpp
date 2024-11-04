/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2024-10-30

  Copyright (C) 2024 Université de Strasbourg

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
   \file test_logging.cpp
   \author Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
   \date 2024-10-30
 */
#define USE_BOOST_TEST 1
// Boost.Test

#define BOOST_TEST_MODULE logging testsuite


#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/fmt.hpp>
#include <feel/feelcore/testsuite.hpp>


inline Feel::po::options_description
makeOptions()
{
    return Feel::feel_options();
}

inline Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "test_logging",
                           "test_logging",
                           "0.2",
                           "Environment class tests",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2013-2016 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@cemosis.fr", "" );
    return about;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( logging )

BOOST_AUTO_TEST_CASE( test_logging1 )
{
    using namespace Feel;
    LOG(INFO) << "This is an info message";
    LOG(WARNING) << "This is a warning message";
    LOG(ERROR) << "This is an error message";
    // demonstrate log every n times
    for (int i = 0; i < 100; ++i) {
        LOG_EVERY_N(INFO, 10) << fmt::format("This is an info message every 10 times, this is the {}th time", i);
    }
    // conditional logging
    for (int i = 0; i < 100; ++i) {
        LOG_IF(INFO, i % 10 == 0) << fmt::format("This is an info message every 10 times, this is the {}th time", i);
    }
    // combine conditional and every n times logging
    for (int i = 0; i < 100; ++i) {
        LOG_IF_EVERY_N(INFO, i % 15 == 0, 10) << fmt::format("This is an info message every 10 times, this is the {}th time",i);
    }
}

BOOST_AUTO_TEST_SUITE_END()
