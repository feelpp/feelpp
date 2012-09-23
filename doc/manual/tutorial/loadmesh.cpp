/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2011-06-21

   Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
 \file testload.cpp
 \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 \date 16-06-2011
 */
#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelvf/vf.hpp>

int main( int argc, char** argv )
{
    // Declare the supported options.
    namespace po = boost::program_options;
    po::options_description desc( "Allowed options" );


    using namespace Feel;

    Environment env( _argc=argc, _argv=argv,
                     _desc=desc,
                     _about=about(_name="loadmesh",
                                  _author="Christophe Prud'homme",
                                  _email="christophe.prudhomme@feelpp.org") );

    typedef Mesh<Simplex<3> > mesh_type;
    fs::path mesh_name=Environment::vm()["filename"].as<std::string>();
    auto mesh =  mesh_type::New();
    if ( mesh_name.extension() == ".geo" )
        mesh = createGMSHMesh( _mesh=mesh,
                               _desc=geo( _filename=mesh_name.string(),_depends=Environment::vm()["depends"].as<std::string>() ),
                               _physical_are_elementary_regions=true,
                               _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES );
    else if ( mesh_name.extension() == ".msh" )
        mesh = loadGMSHMesh( _mesh=mesh,
                             _filename=mesh_name.string(),
                             _rebuild_partitions=true,
                             _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES );
    else
    {
        std::cerr << "invalid file name : " << mesh_name << "\n";
        return 0;
    }

    std::cout << "mesh " << mesh_name << " loaded\n" << std::endl;

    std::cout << "volume =" << std::endl
              << integrate( elements( mesh ), cst( 1. ) ).evaluate() << "\n";
    std::cout << "surface =" << std::endl
              << integrate( boundaryfaces( mesh ), cst( 1. ) ).evaluate() << "\n";


}
