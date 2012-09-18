/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-04-27

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

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
   \file test_backend.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-04-27
 */
#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE backend testsuite

#if defined(USE_BOOST_TEST)
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;
#include <boost/test/floating_point_comparison.hpp>
#endif

#include <boost/timer.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>

#define FEELAPP( argc, argv, about, options )                           \
    Feel::Application app( argc, argv, about, options );                \
    if ( app.vm().count( "help" ) )                                     \
    {                                                                   \
        std::cout << app.optionsDescription() << "\n";                  \
    }



namespace Feel
{
inline
po::options_description
makeOptions()
{
    po::options_description simgetoptions( "test_simget options" );
    simgetoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.5 ), "mesh size" )
    ;
    return simgetoptions.add( Feel::feel_options() );
}


inline
AboutData
makeAbout()
{
    AboutData about( "test_simget" ,
                     "test_simget" ,
                     "0.1",
                     "SimGet tests",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2010 Université Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

/**
 *
 */
class sim : public Simget
{
public:
    sim( po::variables_map const& vm, AboutData const& about )
        :
        Simget( vm, about ),
        meshSize( this->vm()["hsize"].as<double>() )
        {
        }
    void run()

        {
            BOOST_CHECK( mpi::environment::initialized() );
            BOOST_CHECK_EQUAL( detail::BackendManager::instance().size(), 0 );
            auto b = backend(_name="test1");
            BOOST_CHECK_EQUAL( detail::BackendManager::instance().size(), 1 );
            BOOST_CHECK_EQUAL( detail::BackendManager::instance().begin()->first.first, BACKEND_PETSC );
            BOOST_CHECK_EQUAL( detail::BackendManager::instance().begin()->first.second, "test1" );
            BOOST_FOREACH( auto b, detail::BackendManager::instance() )
            {
                BOOST_TEST_MESSAGE( "[test_backend] backend name =" << b.first.second << " kind =" << b.first.first );;
            }
            auto v = b->newVector( 10, 10 );
        }
    void run( const double*, long unsigned int, double*, long unsigned int ) {}
private:
    double meshSize;
};
} // Feel

using namespace Feel;

//BOOST_GLOBAL_FIXTURE( Environment )

BOOST_AUTO_TEST_SUITE( backendsuite )

BOOST_AUTO_TEST_CASE( test_backend1 )
{
    Environment env;

    BOOST_TEST_MESSAGE( "test_backend1" );
    BOOST_CHECK( Environment::initialized() );
    //BOOST_CHECK( mpi::environment::initialized() );
    BOOST_TEST_MESSAGE( "initializing the Application" );

    FEELAPP( 0,//boost::unit_test::framework::master_test_suite().argc,
             boost::unit_test::framework::master_test_suite().argv,
             makeAbout(), makeOptions() );
    BOOST_CHECK( mpi::environment::initialized() );
    BOOST_TEST_MESSAGE( "adding simget" );
    app.add( new sim( app.vm(), app.about() ) );
    BOOST_TEST_MESSAGE( "run simget" );
    app.run();


    BOOST_TEST_MESSAGE( "test_backend1 done" );

}

BOOST_AUTO_TEST_SUITE_END()

