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
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/partitionio.hpp>
#include <feel/feelmesh/remesh.hpp>
#include <feel/feelcore/json.hpp>

//#include <feel/feelfilters/savegmshmesh.hpp>
#include "dump.hpp"

namespace Feel {

namespace nl = nlohmann;

template <typename ShapeType>
void partition( std::vector<int> const& nParts, nl::json const& partconfig )
{
    typedef Mesh<ShapeType> mesh_type;

    // only master rank because only metis is supported at this time
    if ( Environment::isMasterRank() )
    {
        std::cout << "run partioner with shape " << ShapeType::name() << "\n";
        fs::path inputPathMesh = fs::system_complete( soption("ifile") );

        tic();
        size_type update_ = MESH_UPDATE_ELEMENTS_ADJACENCY|MESH_NO_UPDATE_MEASURES|MESH_GEOMAP_NOT_CACHED;
        if ( partconfig && partconfig["partitioner"].contains( "aggregates" ) ) //boption( "sc.ibc_partitioning" ) )
            update_ |= MESH_UPDATE_FACES_MINIMAL;
        auto mesh = loadMesh(_mesh=new mesh_type(Environment::worldCommSeqPtr()), _savehdf5=0,
                             _filename=inputPathMesh.string(),
                             _update=update_,
                             _straighten=false );
        //_update=size_type(MESH_UPDATE_FACES_MINIMAL|MESH_NO_UPDATE_MEASURES));
        //_update=size_type(MESH_UPDATE_FACES|MESH_UPDATE_EDGES));
        toc("loading mesh done",FLAGS_v>0);

        if constexpr ( is_simplex_v<ShapeType> )
        {
            if ( boption( "remesh" ) )
            {
                auto Xh = Pch<1>( mesh );
                auto met = Xh->element();
                if ( !soption( "remesh.metric" ).empty() )
                    met.on( _range=elements(mesh), _expr=expr(soption("remesh.metric")) );
                else
                {
                    auto [ havg, hmin, hmax ] = hMeasures( mesh );
                    std::map<std::string, double> hs { {"havg",havg },{"hmin",hmin },{"hmax",hmax } };
                    double h_ = hs.at("hmax");
                    if ( hs.count( soption("remesh.h") ) )
                        h_=hs.at( soption("remesh.h"));
                    met.on( _range=elements(mesh), _expr=cst(h_)  );
                }
                auto r =  remesher( mesh );
    
                r.setMetric( met );
                auto out = r.execute();
                out->updateForUse();
                mesh = out;
            }
        }
        dump(mesh);

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
            std::string outputFilenameWithoutExt = "";
            if( Environment::vm().count("ofile") )
            {
                std::string ofile = fs::path( soption("ofile")).stem().string();
                if( nParts.size() == 1 )
                    outputFilenameWithoutExt = ofile;
                else
                    outputFilenameWithoutExt = (boost::format("%1%_p%2%")%ofile %nPartition ).str();
            }
            else
                outputFilenameWithoutExt = (boost::format("%1%_p%2%")%inputFilenameWithoutExt %nPartition ).str();
            std::string outputFilenameWithExt = outputFilenameWithoutExt + ".json";
            fs::path outputDirPath;
            if ( Environment::vm().count("odir") )
                outputDirPath = fs::system_complete( soption("odir") );
            else
                outputDirPath = fs::current_path();
            std::string outputPathMesh = ( outputDirPath / fs::path(outputFilenameWithExt) ).string();

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
            io.write( partitionMesh( mesh, nPartition, partitionByRange, partconfig ) );
            toc("paritioning and save on disk done",FLAGS_v>0);
        }
    }

}



}
#endif
