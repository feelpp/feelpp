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
   \file test_taskflow.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2022-02-20
 */
#define USE_BOOST_TEST 1
// Boost.Test

#define BOOST_TEST_MODULE taskflow testsuite
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <taskflow/taskflow.hpp>
#include <taskflow/algorithm/reduce.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/chrono.h>

inline Feel::po::options_description
makeOptions()
{
    Feel::po::options_description options( "Test Environment options" );
    options.add_options()
        ( "nt", Feel::po::value<int>()->default_value( 2 ), "number of threads" )
        ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "h value" )
        ( "fname", Feel::po::value<std::string>()->default_value( "toto" ), "file name" );
    return options.add( Feel::feel_options() );
}

inline Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "test_taskflow",
                           "test_taskflow",
                           "0.2",
                           "Environment class tests",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2013-2022 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( git )

BOOST_AUTO_TEST_CASE( test_taskflow1 )
{
    using namespace Feel;

    using mesh_t = Mesh<Simplex<3, 1, 3>>;
    using mesh_ptr_t = std::shared_ptr<mesh_t>;
    auto beg = std::chrono::high_resolution_clock::now();
    auto mesh = loadMesh( _mesh = new mesh_t );
    auto end = std::chrono::high_resolution_clock::now();
    BOOST_TEST_MESSAGE( fmt::format( "load mesh in {}", std::chrono::duration_cast<std::chrono::seconds>( end - beg ) ) );

    tf::Executor executor( ioption( "nt" ) );
    tf::Taskflow taskflow;
    // auto [r,meas] = std::tuple{boundaryfaces(mesh),6};
    // auto [r, meas] = std::tuple{markedfaces(mesh,"Walls"), 6};
    // auto [r, meas] = std::tuple{markedfaces(mesh,"1"), 1};
    auto [r, meas] = std::tuple{ elements( mesh ), 1 };
    double area_loc = 0;
    beg = std::chrono::high_resolution_clock::now();
    taskflow.transform_reduce(
                r.get<1>(), r.get<2>(), area_loc,
                []( double a, double b )
                {
                    return a + b;
                },
                []( auto const& wf )
                {
                    auto const& f = boost::unwrap_ref( wf );
                    return f.measure();
                } )
        .name( "compute_measure" );
    executor.run( taskflow ).get();
    end = std::chrono::high_resolution_clock::now();
    double area= 0;
    mpi::all_reduce( Environment::worldComm().comm(), area_loc, area, std::plus<double>() );
    BOOST_TEST_MESSAGE( fmt::format("area domain tf({} threads): {} in {}", ioption( "nt" ), area, std::chrono::duration_cast<std::chrono::milliseconds>( end - beg ) ) );
    beg = std::chrono::high_resolution_clock::now();
    double area_seq_loc = 0, area_seq = 0;
    for ( auto const& wf : r )
    {
        auto const& f = boost::unwrap_ref( wf );
        area_seq_loc += f.measure();
    }
    end = std::chrono::high_resolution_clock::now();
    mpi::all_reduce(Environment::worldComm().comm(), area_seq_loc, area_seq, std::plus<double>());
    BOOST_TEST_MESSAGE( fmt::format("area domain seq: {} in {}", area_seq, std::chrono::duration_cast<std::chrono::milliseconds>( end - beg ) ) );

    BOOST_CHECK_CLOSE( area, area_seq, 1e-10);
    BOOST_CHECK_CLOSE( area, meas, 1e-10 );
}

BOOST_AUTO_TEST_SUITE_END()
