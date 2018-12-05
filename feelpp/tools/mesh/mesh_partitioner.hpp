//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
//!
//! This file is part of the Feel++ library
//!
//! Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! Date: 07 Jan 2017
//!
//! Copyright (C) 2017 Feel++ Consortium
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
#ifndef FEELPP_MESH_PARTITIONER_HPP
#define FEELPP_MESH_PARTITIONER_HPP 1

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelmesh/partitionmesh.hpp>
#include <feel/feelfilters/partitionio.hpp>
//#include <feel/feelfilters/savegmshmesh.hpp>

namespace Feel {


template <typename ShapeType>
void partition( std::vector<int> const& nParts)
{
    typedef Mesh<ShapeType> mesh_type;

    // only master rank because only metis is supported at this time
    if ( Environment::isMasterRank() )
    {
        std::cout << "run partioner with shape " << ShapeType::name() << "\n";
        fs::path inputPathMesh = fs::system_complete( soption("ifile") );

        tic();
        auto mesh = loadMesh(_mesh=new mesh_type(Environment::worldCommSeqPtr()), _savehdf5=0,
                             _filename=inputPathMesh.string(),
                             _update=size_type(MESH_UPDATE_ELEMENTS_ADJACENCY|MESH_NO_UPDATE_MEASURES|MESH_GEOMAP_NOT_CACHED),
                             _straighten=false );
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

        std::vector<elements_reference_wrapper_t<mesh_type>> partitionByRange;
        if ( Environment::vm().count("by-markers-desc") )
        {
            std::vector<std::string> inputMarkers = Environment::vm()["by-markers-desc"].template as<std::vector<std::string> >();
            std::string inputMarkersAsString;
            for ( std::string const& marker : inputMarkers )
                inputMarkersAsString += marker;

            boost::char_separator<char> sep(":");
            boost::char_separator<char> sep2(",");
            boost::tokenizer< boost::char_separator<char> > kvlist( inputMarkersAsString, sep );
            for( const auto& ikvl : kvlist )
            {
                boost::tokenizer< boost::char_separator<char> > kvlist2( ikvl, sep2);
                std::vector<std::string> markerList;
                for( const auto& ikvl2 : kvlist2 )
                {
                    markerList.push_back( ikvl2 );
                }
                if ( !markerList.empty() )
                    partitionByRange.push_back( markedelements(mesh,markerList) );
            }
        }
        else if ( Environment::vm().count("by-markers") )
        {
            for ( auto const& markPair : mesh->markerNames() )
            {
                std::string marker = markPair.first;
                if ( mesh->markerDim( marker ) == mesh_type::nDim ) // mesh->hasElementMarker( marker ) )
                    partitionByRange.push_back( markedelements(mesh,marker) );
            }
        }

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
            io.write( partitionMesh( mesh, nPartition, partitionByRange ) );
            toc("paritioning and save on disk done",FLAGS_v>0);
        }
    }

}



}
#endif
