/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

   This file is part of the Feel++ library

   Author(s): Christophe Prud'homme
   Date     : Fri Mar 21 18:10:49 2014

   Copyright (C) 2014-2016 Feel++ Consortium

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
#include <feel/options.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/product.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/loadgmshmesh.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>


int main(int argc, char**argv )
{
    using namespace Feel;
    // initialize feel++
    Environment env( _argc=argc, _argv=argv,
                     _desc_lib=feel_options().add( backend_options( "u" ) ).add ( backend_options( "gradu" ) ),
                     _about=about(_name="productspace_dofs",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    double meshSize = doption("gmsh.hsize");
    Feel::cout << "hsize: " << meshSize << std::endl;

    std::string geofile = soption("gmsh.filename");
    Feel::cout << "geofile: " << geofile << std::endl;

    auto mesh = loadMesh(_mesh = new Mesh<Simplex<2>>,
                         _filename = geofile, // geofile = "/home/LNCMI-G/trophime/feelpp_build/B_Map/clang-3.7/testsuite/feelvf/cube.geo"
                         _savehdf5=false,
                         _h = meshSize,
                         _force_rebuild = true,
                         _update=MESH_UPDATE_EDGES|MESH_UPDATE_FACES );
    typedef FunctionSpace<Mesh<Simplex<2> >, bases<Lagrange<1, Scalar>, Lagrange<0, Scalar> > > space_type;

    auto Vh = space_type::New( mesh );
    auto U = Vh->element();
    auto u = U.element<0>() ;
    auto l = U.element<1>() ;

    Feel::cout << "Vh Dofs: " << Vh->nDof() << "[" << u.functionSpace()->nDof() << "," << l.functionSpace()->nDof() <<"]" << std::endl;


    Feel::cout << "hsize: " << meshSize/2.0 << std::endl;
    if ( geofile.empty() || geofile == "untitled.geo" )
    {
        std::string filenameExpand = Environment::expand("hypercube.geo");
        fs::path mesh_name=fs::path(Environment::findFile(filenameExpand));
        geofile =  mesh_name.string(); //"hypercube.geo";
    }

    gmsh_ptrtype desc_geo=  geo(_filename = geofile, _h = meshSize/2.0 );
    mesh = createGMSHMesh( _mesh=new Mesh<Simplex<2>>,
                           _desc = desc_geo,
                           _force_rebuild = true,
                           _update=MESH_UPDATE_FACES | MESH_UPDATE_EDGES);

    Vh = space_type::New( mesh );
    U = Vh->element();
    u = U.element<0>() ;
    l = U.element<1>() ;

    Feel::cout << "Vh Dofs: " << Vh->nDof() << "[" << u.functionSpace()->nDof() << "," << l.functionSpace()->nDof() <<"]" << std::endl;
}
