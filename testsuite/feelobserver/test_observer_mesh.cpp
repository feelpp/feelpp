// This file is part of the Feel library
//
// Author(s): Feel++ Contortium
//      Date: 2017-07-10
//
// @copyright (C) 2017 University of Strasbourg
// @copyright (C) 2012-2017 Feel++ Consortium
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

#define BOOST_TEST_MODULE test_observers
#include <feel/feelcore/testsuite.hpp>

#include <feel/feel.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( observers )

BOOST_AUTO_TEST_CASE( journal_mesh )
{
    tic();
    const auto& sigptr = Environment::journalSignal();
    int na = sigptr->num_slots();

    auto mesh1 = loadMesh( _mesh=new Mesh<Simplex<1,1>>, _name="Mesh1" );
    auto mesh2 = loadMesh( _mesh=new Mesh<Simplex<2,1>>, _name="Mesh2" );
    auto mesh3 = loadMesh( _mesh=new Mesh<Simplex<3,1>>, _name="Mesh3" );
    // Mesh 4 and 5 have no name. It will be saved as "undefined", mesh5
    // overide info!
    auto mesh4 = loadMesh( _mesh=new Mesh<Simplex<2,1>> );
    auto mesh5 = loadMesh( _mesh=new Mesh<Simplex<2,1>> );

    auto mesh6 = loadMesh( _mesh=new Mesh<Simplex<3,1>> );
    mesh6->journalWatcherName( "Mesh6", false);
    toc("loadMesh");

    // By default, mesh is watchable. You can connect the mesh
    // to be observed and its data will be retrieved by the observer
    // manager (Environment class by default).

    // only mesh1 and mesh3 will be observed.
    mesh1->journalConnect();
    mesh2->journalConnect();
    //mesh3 automatic connection.
    mesh4->journalConnect();
    mesh5->journalConnect();
    mesh6->journalConnect();

    // Mesh 2 won't be observed anymore.
    mesh2->journalDisconnect();

    int nb = sigptr->num_slots();
    CHECK( nb-na == 5 );
    // This create a checkpoint and save the result in a json file.
    Environment::journalFinalize();
}

BOOST_AUTO_TEST_SUITE_END()

// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
