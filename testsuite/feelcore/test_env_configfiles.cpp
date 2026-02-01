/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/*
  This file is part of the Feel++ library

  Author(s): Feel++ Contributors
       Date: 2025-10-28

  Copyright (C) 2025 Feel++ Consortium

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

#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE environment_config_files testsuite

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/testsuite.hpp>

#include <fstream>
#include <random>

using namespace Feel;

inline po::options_description makeOptions()
{
    po::options_description opts( "Environment config files test options" );
    opts.add_options()
        ( "single", po::value<std::string>()->default_value( "default-single" ), "single file option" )
        ( "override", po::value<std::string>()->default_value( "default-override" ), "overridden option" );
    return opts.add( feel_options() );
}

inline AboutData makeAbout()
{
    AboutData about( "test_env_configfiles",
                     "test_env_configfiles",
                     "0.1",
                     "Environment config files tests",
                     AboutData::License_GPL,
                     "Copyright (C) 2025 Feel++ Consortium" );
    about.addAuthor( "Feel++ Contributors", "developer", "feelpp-devel@feelpp.org", "" );
    return about;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( environment_config_files )

BOOST_AUTO_TEST_CASE( multiple_config_files_respect_priority )
{
    // Generate a unique directory name using random number
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1000, 9999);
    
    std::string unique_suffix = "feelpp-config-" + std::to_string(dis(gen)) + "-" + std::to_string(dis(gen));
    fs::path const tempDir = fs::temp_directory_path() / unique_suffix;
    fs::create_directories( tempDir );

    fs::path const firstCfg = tempDir / "first.cfg";
    fs::path const secondCfg = tempDir / "second.cfg";

    {
        std::ofstream ofs( firstCfg.string() );
        ofs << "single=from-first\n";
        ofs << "override=from-first\n";
    }

    {
        std::ofstream ofs( secondCfg.string() );
        ofs << "override=from-second\n";
    }

    Environment::setConfigFiles( { firstCfg.string(), secondCfg.string() } );

    BOOST_CHECK_EQUAL( option(_name="single").as<std::string>(), "from-first" );
    BOOST_CHECK_EQUAL( option(_name="override").as<std::string>(), "from-second" );

    fs::remove_all( tempDir );
}

BOOST_AUTO_TEST_SUITE_END()

