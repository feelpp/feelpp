/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 19 Feb 2016

 Copyright (C) 2016 Feel++ Consortium

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
#define BOOST_TEST_MODULE test_feelio
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/logger.hpp>

/** use Feel namespace */
using namespace Feel;

inline
po::options_description makeOptions()
{
    po::options_description options( "Test Feelio  Options" );
    return options;
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_feelio" ,
                     "test_feelio" ,
                     "0.2",
                     "nD(n=2,3) test feelio",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2016 Feel++ Consortium" );

    about.addAuthor( "C Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}



FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
BOOST_AUTO_TEST_SUITE( feelio_suite )

BOOST_AUTO_TEST_CASE( test_cout_cerr_clog )
{
    BOOST_MESSAGE( "test_cout" );
    BOOST_CHECK( cout.worldCommPtr() );
    BOOST_CHECK( cerr.worldCommPtr() );
    BOOST_CHECK( clog.worldCommPtr() );
    cout << "cout:: Hello World from process " << Environment::rank() << std::endl;
    cerr << "cerr:: Hello World from process " << Environment::rank() << std::endl;
    clog << "clog:: Hello World from process " << Environment::rank() << std::endl;
    BOOST_MESSAGE( "test_cout done." );
}

BOOST_AUTO_TEST_CASE( test_stringstream )
{
    BOOST_MESSAGE( "test_stringstream" );
    std::ostringstream os;
    MasterStream fs( os ) ;
    fs << "Hello World from process " << Environment::rank();
    
    // MasterStream now writes on all ranks (behavior changed with new logging)
    std::string expected = "Hello World from process " + std::to_string(Environment::rank());
    BOOST_TEST_MESSAGE( "str:: --" << fs.str() << "--\n");
    BOOST_CHECK_EQUAL( fs.str(), expected );

    BOOST_MESSAGE( "test_fstream done." );
}

BOOST_AUTO_TEST_CASE( test_log_console )
{
    BOOST_MESSAGE( "test_log_console" );
    
    // Test console logger - respects log.mpi setting
    Logger::console()->info("log::console:: Hello from rank {}", Environment::rank());
    Logger::console()->warn("log::console:: Warning from rank {}", Environment::rank());
    Logger::console()->error("log::console:: Error from rank {}", Environment::rank());
    Logger::console()->flush();
    
    // Test file logger (default logger)
    Logger::file()->info("log::file:: This goes to file from rank {}", Environment::rank());
    Logger::file()->debug("log::file:: Debug message from rank {}", Environment::rank());
    Logger::file()->flush();
    
    // Test null logger
    Logger::null()->info("log::null:: This is discarded from rank {}", Environment::rank());
    
    // Verify loggers exist
    BOOST_CHECK( Logger::console() != nullptr );
    BOOST_CHECK( Logger::file() != nullptr );
    BOOST_CHECK( Logger::null() != nullptr );
    
    BOOST_MESSAGE( "test_log_console done." );
}


BOOST_AUTO_TEST_SUITE_END()
