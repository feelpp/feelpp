/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-02-01

  Copyright (C) 2005,2006 EPFL

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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-02-01
 */
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <iostream>
#include <boost/lambda/lambda.hpp>

#include <life/lifealg/boundingbox.hpp>
#include <life/lifemesh/regiontree.hpp>

#include <life/lifediscr/mesh.hpp>

#include <life/lifefilters/importergmsh.hpp>
#include <life/lifefilters/gmsh.hpp>



using namespace Life;
namespace lambda = boost::lambda;


int main( int /*argc*/,  char** /*argv*/ )
{
    typedef Mesh<Simplex<3,1> > mesh_type;
    mesh_type aMesh;

    Gmsh __gmsh;
    std::string fname = __gmsh.generateCube( "cube", 1 );

    ImporterGmsh<mesh_type> gmsh_import( fname );
    aMesh.accept( gmsh_import );

    RegionTree __rt;
    typedef node<double>::type node_type;
    node_type min, max;
    scalar_type EPS=1E-13;

    __rt.clear();

    BoundingBox<> bb( true );
    for ( size_type __i = 0; __i < aMesh.numElements(); ++__i )
    {
        bb.make( aMesh.element( __i ).G() );
        for (unsigned k=0; k < min.size(); ++k)
        {
            bb.min[k]-=EPS;
            bb.max[k]+=EPS;
        }
        __rt.addBox(bb.min, bb.max, __i );
    }

    __rt.dump();

    node_type __n( 3 );

    __n[0]=0;
    __n[1]=0;
    __n[2]=0;

    //size_type __cv;
    // FIXME: implementation of function below must have been lost
    //during some refactoring
    // find_a_point( __rt, __n, __n, __cv );
}
