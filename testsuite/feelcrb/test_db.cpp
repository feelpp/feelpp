/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Alexandre Ancel <alexandre.ancel@cemosis.fr>
       Date: 2016-04-06

  Copyright (C) 2016 Universite de Strasbourg

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
   \file test_db.cpp
   \author Alexandre Ancel <alexandre.ancel@cemosis.fr>
   \date 2016-04-06
 */
#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE eim testsuite

#include <vector>

#include <testsuite/testsuite.hpp>

#include <boost/timer.hpp>
#include <boost/smart_ptr/enable_shared_from_this.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/unitsquare.hpp>

#include <feel/feelcrb/opusapp.hpp>
#include <applications/crb/heat1d/heat1d.hpp>

namespace Feel
{
inline
po::options_description
makeOptions()
{
    po::options_description dboptions( "test_db options" );
    /*
    dboptions.add_options()
    ( "cvg-study" , po::value<bool>()->default_value( false ), "run a convergence study if true" )
    ;
    */
    return dboptions;
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_db" ,
                     "test_db" ,
                     "0.1",
                     "Database tests",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2016 Universite de Strasbourg" );

    about.addAuthor( "Alexandre Ancel", "developer", "alexandre.ancel@cemosis.fr", "" );
    return about;

}
};

using namespace Feel;

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( db )

BOOST_AUTO_TEST_CASE( test_db1 )
{
    BOOST_TEST_MESSAGE( "test_db starts..." );

    using namespace Feel;

    char ** argv;
    int argc = 3;
    argv = new char*[argc + 1];
    argv[0] = strdup("test_db");
    argv[1] = strdup("--config-file");
    argv[2] = strdup("/ssd/ancel/feelpp/clang/build/testsuite/feelcrb/heat1d.cfg");
    argv[3] = NULL;

    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=opusapp_options("heat1d")
                           .add(crbOptions())
                           .add(crbSEROptions())
                           .add(makeHeat1DOptions())
                           .add(eimOptions())
                           .add(podOptions())
                           .add(backend_options("backend-primal"))
                           .add(backend_options("backend-dual"))
                           .add(backend_options("backend-l2"))
                           .add(bdf_options("Heat1D")),
                           _about=makeHeat1DAbout( "heat1d" ) );

    std::vector<int> wnSize;
    std::vector<int> wnduSize;

    // Setup initial options
    Environment::setOptionValue("crb.results-repo-name", std::string("test_db"));
    Environment::setOptionValue("crb.rebuild-database", true);

    // Execute the Heat1d app with both boost and hdf5 databases
    Environment::setOptionValue("crb.db.format", std::string("boost"));
    Feel::OpusApp<Feel::Heat1D > * app = new Feel::OpusApp<Feel::Heat1D >();
    app->run();

    auto crb = app->getCRB();
    auto WN = crb->wn();
    auto WNdu = crb->wndu();
    wnSize.push_back(WN.size());
    wnduSize.push_back(WNdu.size());

    delete app;

    Environment::setOptionValue("crb.db.format", std::string("hdf5"));
    app = new Feel::OpusApp<Feel::Heat1D >();
    app->run();

    crb = app->getCRB();
    WN = crb->wn();
    WNdu = crb->wndu();
    wnSize.push_back(WN.size());
    wnduSize.push_back(WNdu.size());

    delete app;

    // Execute the Heat1d app with both boost and hdf5 databases
    // Except we are reloading the databases here
    Environment::setOptionValue("crb.rebuild-database", false);

    Environment::setOptionValue("crb.db.format", std::string("boost"));
    app = new Feel::OpusApp<Feel::Heat1D >();
    app->run();

    crb = app->getCRB();
    WN = crb->wn();
    WNdu = crb->wndu();
    wnSize.push_back(WN.size());
    wnduSize.push_back(WNdu.size());

    delete app;

    Environment::setOptionValue("crb.db.format", std::string("hdf5"));
    app = new Feel::OpusApp<Feel::Heat1D >();
    app->run();

    crb = app->getCRB();
    WN = crb->wn();
    WNdu = crb->wndu();
    wnSize.push_back(WN.size());
    wnduSize.push_back(WNdu.size());

    delete app;

    /* Check that we have the correct basis sizes */
    for(int i = 1; i < wnSize.size(); i++)
    {
        BOOST_CHECK(wnSize[0] == wnSize[i]);
    }
 
    for(int i = 1; i < wnSize.size(); i++)
    {
        BOOST_CHECK(wnSize[0] == wnSize[i]);
    }   

    BOOST_TEST_MESSAGE( "test_db done" );
}

BOOST_AUTO_TEST_SUITE_END()
