/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 16 september 2019

 Copyright (C) 2019 Feel++ Consortium

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

// give a name to the testsuite
#define BOOST_TEST_MODULE importer mesh testsuite

#include <feel/feelcore/testsuite.hpp>

#include <feel/feelfilters/loadmesh.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( test_importer_mesh )

typedef boost::mpl::list<boost::mpl::int_<1>,boost::mpl::int_<2>,boost::mpl::int_<3> > dim_types;
using namespace Feel;

BOOST_AUTO_TEST_CASE_TEMPLATE( test_scale, T, dim_types )
{
    typedef Mesh<Simplex<T::value> > mesh_type;
    auto mesh = loadMesh( _mesh=new mesh_type);
    double scalingUsed = 10;
    auto meshWithScale = loadMesh( _mesh=new mesh_type,
                                   _scale=scalingUsed );

    double volume = mesh->measure();
    double volumeWithScale = meshWithScale->measure();
    BOOST_CHECK_CLOSE( volume*std::pow(scalingUsed, mesh_type::nDim), volumeWithScale, 1e-10 );

    // save mesh with hdf5 format
    double scalingUsed2 = 1./3.;
    std::string filenameHDF5 = (fs::current_path()/(boost::format("test_scale_mesh_%1%d.json")%mesh_type::nDim).str()).string();
    mesh->saveHDF5( filenameHDF5, scalingUsed2 );
    // load mesh with hdf5 format
    double scalingUsed3 = 5;
    auto meshHDF5 = loadMesh( _mesh=new mesh_type,_filename=filenameHDF5,_scale=scalingUsed3);
    double volumeHDF5 = meshHDF5->measure();
    BOOST_CHECK_CLOSE( volume*std::pow(scalingUsed2*scalingUsed3, mesh_type::nDim), volumeHDF5, 1e-10 );
}

BOOST_AUTO_TEST_SUITE_END()
