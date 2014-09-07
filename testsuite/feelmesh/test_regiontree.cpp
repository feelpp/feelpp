/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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
#include <boost/lambda/lambda.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelalg/boundingbox.hpp>
#include <feel/feelmesh/regiontree.hpp>

#include <feel/feeldiscr/mesh.hpp>

#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>



using namespace Feel;
namespace lambda = boost::lambda;


int main( int argc,  char** argv )
{
    Feel::Environment env( argc, argv );
    typedef Mesh<Simplex<3,1> > mesh_type;

    static const uint16_type nDim = mesh_type::nDim;

    auto shape = "hypercube";
    auto aMesh = createGMSHMesh( _mesh=new mesh_type,
                                 _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % nDim ).str() ,
                                         _usenames=true,
                                         _shape=shape,
                                         _dim=nDim,
                                         _h=1 ) );

    RegionTree __rt;
    typedef node<double>::type node_type;
    node_type min, max;
    scalar_type EPS=1E-13;

    __rt.clear();

    BoundingBox<> bb( true );

    for ( auto const& elt : elements(aMesh) )
    {
        bb.make( elt.G() );

        for ( unsigned k=0; k < min.size(); ++k )
        {
            bb.min[k]-=EPS;
            bb.max[k]+=EPS;
        }

        __rt.addBox( bb.min, bb.max, elt.id() );
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
