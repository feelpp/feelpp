/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2016-01-29

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
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelmesh/partitionmesh.hpp>
#include <feel/feelfilters/partitionio.hpp>
//#include <feel/feelfilters/savegmshmesh.hpp>

using namespace Feel;

template <typename ShapeType>
void run( std::vector<int> const& nParts)
{
    typedef Mesh<ShapeType> mesh_type;

    // only master rank because only metis is supported at this time
    if ( Environment::isMasterRank() )
    {
        std::cout << "run partioner with shape " << ShapeType::name() << "\n";
        fs::path inputPathMesh = fs::system_complete( soption("ifile") );

        tic();
        auto mesh = loadMesh(_mesh=new mesh_type(Environment::worldCommSeq()), _savehdf5=0,
                             _filename=inputPathMesh.string(),
                             _update=size_type(MESH_UPDATE_ELEMENTS_ADJACENCY|MESH_NO_UPDATE_MEASURES));
                             //_update=size_type(MESH_UPDATE_FACES_MINIMAL|MESH_NO_UPDATE_MEASURES));
                             //_update=size_type(MESH_UPDATE_FACES|MESH_UPDATE_EDGES));
        toc("loading mesh done",FLAGS_v>0);

        std::cout << "      number of elements in memory : " << mesh->numGlobalElements() << std::endl;
        std::cout << "      number of faces in memory : " << mesh->numGlobalFaces() << std::endl;
        if ( mesh_type::nDim == 3 )
            std::cout << "      number of edges in memory : " << mesh->numGlobalEdges() << std::endl;
        std::cout << "      number of points  in memory : " << mesh->numGlobalPoints() << std::endl;
        for( auto marker: mesh->markerNames() )
        {
            auto name = marker.first;
            auto data = marker.second;
            if ( data[1] == mesh->dimension()-1 )
            {
                size_type nelts = nelements( markedfaces(mesh, name ), true );
                std::cout << "      number of marked faces " << name << " with tag " << data[0] << " : " << nelts << std::endl;
            }
        }

        //saveGMSHMesh(_mesh=mesh,_filename="tototi.msh");

        fs::path inputDir = inputPathMesh.parent_path();
        std::string inputFilenameWithoutExt = inputPathMesh.stem().string();

        for ( int nPartition : nParts )
        {

            std::string outputFilenameWithoutExt = (boost::format("%1%_p%2%")%inputFilenameWithoutExt %nPartition ).str();
            std::string outputFilenameWithExt = outputFilenameWithoutExt + ".json";
            std::string outputPathMesh = ( fs::current_path() / fs::path(outputFilenameWithExt) ).string();
            if ( Environment::vm().count("ofile") )
            {
                outputPathMesh = fs::system_complete( soption("ofile") ).string();
                if ( fs::path(outputPathMesh).extension() != ".json" )
                    outputPathMesh = outputPathMesh + ".json";
            }
            else if ( Environment::vm().count("odir") )
            {
                fs::path odir = fs::system_complete( soption("odir") ).string();
                outputPathMesh = (odir / fs::path(outputFilenameWithExt) ).string();
            }

            fs::path outputDir = fs::path(outputPathMesh).parent_path();
            if ( !fs::exists( outputDir ) )
                fs::create_directories( outputDir );


            std::cout << "start mesh partitioning and save on disk : " << outputPathMesh << "\n";

            // build a MeshPartitionSet based on a mesh partition that will feed a
            // partition io data structure to generate a parallel hdf5 file from which
            // the parallel mesh can be loaded
            tic();
            using io_t = PartitionIO<mesh_t<decltype(mesh)>>;
            io_t io( outputPathMesh );
            io.write( partitionMesh( mesh, nPartition ) );
            toc("paritioning and save on disk done",FLAGS_v>0);
        }
    }

}



int main( int argc, char** argv )
{
	po::options_description meshpartoptions( "Mesh Partitioner options" );
	meshpartoptions.add_options()
        ( "dim", po::value<int>()->default_value( 3 ), "mesh dimension" )
        ( "shape", po::value<std::string>()->default_value( "simplex" ), "mesh dimension" )
        ( "part", po::value<std::vector<int> >()->multitoken(), "number of partition" )
        ( "ifile", po::value<std::string>(), "input mesh filename" )
        ( "odir", po::value<std::string>(), "output directory [optional]" )
        ( "ofile", po::value<std::string>(), "output mesh filename [optional]" )
		;

    Environment env( _argc=argc, _argv=argv,
                     _desc=meshpartoptions,
                     _about=about( _name="mesh_partitioner" ,
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ) );

    int dim = ioption(_name="dim");
    std::string shape = soption(_name="shape");

    std::vector<int> nParts;
    if ( Environment::vm().count("part"))
        nParts = Environment::vm()["part"].as<std::vector<int> >();

    if ( nParts.empty() )
    {
        if ( Environment::isMasterRank() )
            std::cout << "do nothing because --part is missing\n";
        return 0;
    }

    if ( !Environment::vm().count("ifile") )
    {
        if ( Environment::isMasterRank() )
            std::cout << "do nothing because --ifile is missing\n";
        return 0;
    }

    fs::path pathInputMesh = fs::system_complete( soption("ifile") );
    if ( !fs::exists( pathInputMesh ) )
    {
        if ( Environment::isMasterRank() )
            std::cout << "do nothing because --ifile : " << pathInputMesh.string() << "not exist\n";
        return 0;
    }

    if ( dim == 1 )
    {
        run<Simplex<1>>( nParts );
    }
    else
    {
        if ( shape == "simplex" )
        {
            switch ( dim )
            {
            case 2 : run<Simplex<2>>( nParts ); break;
            case 3 : run<Simplex<3>>( nParts ); break;
            }
        }
        else if ( shape == "hypercube" )
        {
            switch ( dim )
            {
            case 2 : run<Hypercube<2>>( nParts ); break;
            case 3 : run<Hypercube<3>>( nParts ); break;
            }
        }
    }

    return 0;

}

