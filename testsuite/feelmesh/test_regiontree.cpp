/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-02-01

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2010-2014 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file test_regiontree.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-02-01
 */
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <iostream>

#include <feel/feelcore/environment.hpp>
#include <feel/feelalg/boundingbox.hpp>
#include <feel/feelmesh/regiontree.hpp>

#include <feel/feeldiscr/mesh.hpp>

#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>



using namespace Feel;


int main(int argc, char** argv)
{
    Feel::Environment env(argc, argv);
    using mesh_type = Mesh<Simplex<3, 1>>;
    constexpr uint16_type nDim = mesh_type::nDim;

    std::string shape = "hypercube";
    auto aMesh = createGMSHMesh(_mesh = new mesh_type,
                                _desc = domain(_name = (boost::format("%1%-%2%") % shape % nDim).str(),
                                               _usenames = true,
                                               _shape = shape,
                                               _dim = nDim,
                                               _h = 1));

    RegionTree regionTree;
    using node_type = node<double>::type;
    node_type min, max;
    scalar_type EPS = 1E-13;

    regionTree.clear();

    BoundingBox<> boundingBox(true);

    for (const auto& element : elements(aMesh))
    {
        boundingBox.make(boost::unwrap_ref(element).G());

        for (auto& min_val : boundingBox.min) 
        {
            min_val -= EPS;
        }

        for (auto& max_val : boundingBox.max) 
        {
            max_val += EPS;
        }

        regionTree.addBox(boundingBox.min, boundingBox.max, boost::unwrap_ref(element).id());
    }

    regionTree.dump();

    node_type point(3);
    point[0] = 0;
    point[1] = 0;
    point[2] = 0;

    // size_type cv;
    // FIXME: implementation of function below must have been lost
    // during some refactoring
    // find_a_point(regionTree, point, point, cv);
}