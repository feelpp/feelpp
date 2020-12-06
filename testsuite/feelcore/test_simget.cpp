/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-09-12

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
   \file test_simget.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-09-12
 */

// give a name to the testsuite
#define BOOST_TEST_MODULE simget testsuite
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/options.hpp>


using namespace Feel;

namespace test_simget
{
inline
po::options_description
makeOptions()
{
    po::options_description simgetoptions( "test_simget options" );
    simgetoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.5 ), "mesh size" )
    ;
    return simgetoptions;
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
    sim()
        :
        Simget(),
        meshSize( doption(_name="hsize") )
    {
    }
    void run() override {}
    void run( const double*, long unsigned int, double*, long unsigned int ) override {}
private:
    double meshSize;
};
} // test_simget

BOOST_AUTO_TEST_SUITE( simget )

FEELPP_ENVIRONMENT_WITH_OPTIONS( test_simget::makeAbout(), test_simget::makeOptions() )

BOOST_AUTO_TEST_CASE( test_sim1 )
{
    BOOST_TEST_MESSAGE( "test_sim1" );
    BOOST_CHECK( Feel::Environment::initialized() );
    Feel::Application app;
    app.add( new test_simget::sim );
    app.run();
    BOOST_TEST_MESSAGE( "test_sim1 done" );
}

BOOST_AUTO_TEST_SUITE_END()

