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

#include <vector>

#include <testsuite/testsuite.hpp>

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

int main(int argc, char **argv)
{
    char ** argv1;
    int argc1 = 3;
    argv1 = new char*[argc1 + 1];
    argv1[0] = strdup("test_db");
    argv1[1] = strdup("--config-file");
    argv1[2] = strdup("/ssd/ancel/feelpp/clang/build/testsuite/feelcrb/heat1d.cfg");
    argv1[3] = NULL;

    Feel::Environment env( _argc=argc1, _argv=argv1,
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

    int sampleSize = 4;
    std::vector<double> wnSample;

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

    if(WN.size() > 0 && WN[0].size() > sampleSize)
    {
        for(int i = 0; i < sampleSize; i++)
        { wnSample.push_back(WN[0][i]); }
    }

    delete app;

    Environment::setOptionValue("crb.db.format", std::string("hdf5"));
    app = new Feel::OpusApp<Feel::Heat1D >();
    app->run();

    crb = app->getCRB();
    WN = crb->wn();
    WNdu = crb->wndu();
    wnSize.push_back(WN.size());
    wnduSize.push_back(WNdu.size());

    if(WN.size() > 0 && WN[0].size() > sampleSize)
    {
        for(int i = 0; i < sampleSize; i++)
        { wnSample.push_back(WN[0][i]); }
    }

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

    if(WN.size() > 0 && WN[0].size() > sampleSize)
    {
        for(int i = 0; i < sampleSize; i++)
        { wnSample.push_back(WN[0][i]); }
    }

    delete app;

    Environment::setOptionValue("crb.db.format", std::string("hdf5"));
    app = new Feel::OpusApp<Feel::Heat1D >();
    app->run();

    crb = app->getCRB();
    WN = crb->wn();
    WNdu = crb->wndu();
    wnSize.push_back(WN.size());
    wnduSize.push_back(WNdu.size());

    if(WN.size() > 0 && WN[0].size() > sampleSize)
    {
        for(int i = 0; i < sampleSize; i++)
        { wnSample.push_back(WN[0][i]); }
    }

    /*
    std::cout << wnSample.size() << std::endl;

    for(int i = 0; i < wnSample.size(); i++)
    {
        if(i != 0 && (i % sampleSize) == 0)
        {
            std::cout << std::endl;
        }
        std::cout << wnSample.at(i) << " ";
    }
    std::cout << std::endl;
    */

    delete app;

    /* Check that we have the correct basis sizes */
    for(int i = 1; i < wnSize.size(); i++)
    {
        if(wnSize[0] == wnSize[i])
        {
            return 1;
        }
    }
 
    for(int i = 1; i < wnSize.size(); i++)
    {
        if(wnSize[0] != wnSize[i])
        {
            return 1;
        }
    }   

    return 0;
}
