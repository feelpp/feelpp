/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Author(s): Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
   Date: 2011-16-12

   Copyright (C) 2008-2010 Universite Joseph Fourier (Grenoble I)
   Copyright (C) CNRS

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
   \file mesh_initializer.cpp
   \author Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
   \date 2011-16-12
*/

#ifndef __MESH_INIT_HPP
#define __MESH_INIT_HPP 1

#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/mesh.hpp>

#include <feel/feelfilters/geo.hpp>
#include <feel/feelfilters/loadmesh.hpp>

namespace Feel
{
    template<int Dim, int G_order>
    class mesh_initializer
    {
    public :
        //! geometry entities type composing the mesh, here Simplex in Dimension Dim (Order G_order)
        typedef Simplex<Dim,G_order> convex_type;
        //! mesh type
        typedef Mesh<convex_type> mesh_type;
        //! mesh shared_ptr<> type
        typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

        mesh_ptrtype initializeMesh(po::variables_map const& _vm);
    };
}

template<int Dim, int G_order>
typename Feel::mesh_initializer<Dim, G_order>::mesh_ptrtype
Feel::mesh_initializer<Dim, G_order>::initializeMesh(Feel::po::variables_map const& _vm)
{
    int proc_rank = Environment::worldComm().globalRank();

    mesh_ptrtype mesh = mesh_type::New();

    // Collect general options: hsize defined in gmsh.hsize...
    std::list<std::string> general_options = boost::assign::list_of("hsize")("geofile")("geo_depends")("geofile-path")("meshadapt_method")("meshadapt_type");

    for (auto& option: general_options)
        {
            if( !_vm.count( option) )
                throw std::logic_error( "Option" + option + "have to be defined to initialize the mesh");
        }

    double meshSize = doption("gmsh.hsize");
    std::string geofile = soption("geofile");
    std::string geo_depends = soption("geo_depends");
    std::string geofile_path = soption("geofile-path");
    std::string meshadapt_method = soption("meshadapt_method");
    std::string meshadapt_type = soption("meshadapt_type");
    bool repart = boption("repart");
    bool meshadapt = boption("meshadapt");

    // Feel::cout << "MeshInitializer:\n";
    // Feel::cout << "geofile =" << geofile << "\n";
    // Feel::cout << "geo_depends =" << geo_depends << "\n";
    // Feel::cout << "geofile_path =" << geofile_path << "\n";
    // Feel::cout << "meshadapt =" << meshadapt << "\n";
    // Feel::cout << "meshadapt_method =" << meshadapt_method << "\n";
    // Feel::cout << "meshadapt_type =" << meshadapt_type << "\n" << std::endl;

    std::string access_geofile = geofile;
    if ( !geofile_path.empty() )
        access_geofile = ( boost::format( "%1%/%2%" ) % Environment::expand( geofile_path ) % geofile ).str();
    // Feel::cout << "access_geofile =" << access_geofile << std::endl;
    // Feel::cout << "gmsh.verbosity =" << ioption("gmsh.verbosity") << std::endl;
    LOG(INFO) << "Loading " << access_geofile << "\n" << std::flush;
    if( fs::exists( access_geofile ) )
        {
            mesh = loadMesh( _mesh=new mesh_type,
                             _filename=access_geofile,
                             _h=meshSize,
                             _rebuild_partitions=repart,
                             _savehdf5=boption("mesh.save.enable"),
                             _depends=geo_depends,
                             _verbose=ioption("gmsh.verbosity"),
                             _update=MESH_UPDATE_EDGES|MESH_UPDATE_FACES );
        }
    else
        throw std::logic_error( "initializeMesh: " + access_geofile + " no such file");

    // Mesh adaptation needs geofile and his dependances in current directory
    if( meshadapt && proc_rank == 0)
        {
            // Copy geofile into current directory
            fs::path mesh_name=fs::path(access_geofile);
            std::string geofileWE = mesh_name.stem().string()+".geo";
            std::string geofile_completePath = (boost::format( "%1%/%2%" ) % Environment::expand( geofile_path ) % geofileWE ).str();
            std::string geofile_newPath = (boost::format( "./%1%" ) % geofileWE ).str();

            fs::path file_path( geofile_completePath ); // current geofile path
            fs::path new_path( geofile_newPath ); // to copy geofile in current directory
            try
                {
                    boost::system::error_code ec;
                    if ( !( fs::exists( file_path ) && fs::is_regular_file( file_path ) ) )
                        std::cout << "[mesh_initializer] File : " << file_path << " doesn't exist or is not a regular file" << std::endl;
                    else if ( !fs::exists( new_path )  )
                        fs::copy_file( file_path, fs::path( geofileWE ), fs::copy_option::none );
                }
            catch ( const fs::filesystem_error& e )
                {
                    std::cerr << "mesh_initializer Error: " << e.what() << std::endl;
                }

            // Copy geofile dependances
            std::vector<std::string> depends_on_files;
            algorithm::split( depends_on_files, geo_depends, algorithm::is_any_of( ":,; " ), algorithm::token_compress_on );

            for (auto& _filename: depends_on_files)
                {
                    std::string depends_completePath = (boost::format( "%1%/%2%" ) % Environment::expand( geofile_path ) % _filename ).str();
                    std::string depends_newPath = (boost::format( "./%1%" ) % _filename ).str();
                    fs::path depends_file_path( depends_completePath );
                    fs::path depends_new_path( depends_newPath );

                    try
                        {
                            boost::system::error_code ec;
                            if ( !( fs::exists( depends_file_path ) && fs::is_regular_file( depends_file_path ) ) )
                                std::cout << "[mesh_initializer] File : " << depends_file_path << " doesn't exist or is not a regular file" << std::endl;
                            else if ( !fs::exists( depends_new_path )  )
                                fs::copy_file( depends_file_path, fs::path( _filename ), fs::copy_option::none );
                        }

                    catch ( const fs::filesystem_error& e )
                        {
                            std::cerr << "Error: " << e.what() << std::endl;
                        }
                }
            // Wait for copy of files by proc of rank 0
            if ( mpi::environment::initialized() )
                Environment::worldComm().barrier();

        }


    LOG(INFO) << "Mesh built \n";

    // Check for bad partitionning
    if( nelements( elements(mesh) ) <= 0 )
        {
            std::cout << "Mesh proc[" << proc_rank << "] bad partition - no elements" << std::endl;
        }
    LOG(INFO) << "Mesh check for bad partitioning \n";

    return mesh;
}

#endif // __MESH_INIT_HPP
