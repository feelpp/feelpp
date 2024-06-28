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
#ifndef FEELPP_MESH_PARTITIONER_IMPL_HPP
#define FEELPP_MESH_PARTITIONER_IMPL_HPP 1

#include <mesh_partitioner.hpp>

#include <fmt/core.h>
#include <fmt/color.h>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelmesh/partitionmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/partitionio.hpp>
#if defined( FEELPP_HAS_MMG ) && defined( FEELPP_HAS_PARMMG )
#include <feel/feelmesh/remesh.hpp>
#endif
#include <feel/feelcore/json.hpp>
#include <feel/feelfilters/exporter.hpp>

#include <feel/feelmesh/dump.hpp>

namespace Feel {

template <typename ShapeType>
void partition( nl::json const& partconfig )
{
    typedef Mesh<ShapeType> mesh_type;

    // only master rank because only metis is supported at this time
    if ( Environment::isMasterRank() )
    {
        fmt::print( fmt::emphasis::bold,
                    "** Run partioner with shape {}...\n", ShapeType::name() );
        auto jInput = partconfig.at( "input" );
        fs::path inputPathMesh = jInput.at("filename").template get<std::string>();
        if ( !fs::exists( Environment::expand( inputPathMesh.string() ) ) )
        {
            std::cout << "input filename not exists : " << inputPathMesh.string() << std::endl;
            return;
        }

        auto jPartitioner = partconfig.at( "partitioner" );
        std::vector<int> nParts;
        if ( jPartitioner.contains("number-of-partition") )
        {
            auto jNumberOfPart = jPartitioner.at("number-of-partition");
            if ( jNumberOfPart.is_number_integer() )
                nParts.push_back( jNumberOfPart.template get<int>() );
            else if ( jNumberOfPart.is_array() )
                for ( auto const& [key,val] : jNumberOfPart.items() )
                    nParts.push_back( val.template get<int>() );
        }


        tic();
        size_type update_ = MESH_UPDATE_ELEMENTS_ADJACENCY|MESH_NO_UPDATE_MEASURES|MESH_GEOMAP_NOT_CACHED;
        if ( !partconfig.is_null() && !partconfig.empty() && partconfig["partitioner"].contains( "aggregates" ) ) //boption( "sc.ibc_partitioning" ) )
            update_ |= MESH_UPDATE_FACES_MINIMAL;
        auto mesh = loadMesh(_mesh=new mesh_type(Environment::worldCommSeqPtr()), _savehdf5=0,
                             _filename=inputPathMesh.string(),
                             _update=update_,
                             _straighten=false );
        //_update=size_type(MESH_UPDATE_FACES_MINIMAL|MESH_NO_UPDATE_MEASURES));
        //_update=size_type(MESH_UPDATE_FACES|MESH_UPDATE_EDGES));
        toc("loading mesh done",FLAGS_v>0);

#if 0
        if constexpr ( is_simplex_v<ShapeType> && ShapeType::nDim > 1 )
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
                nl::json j;
                if ( !soption( "remesh.json" ).empty() )
                {
                    std::ifstream ifs( soption( "remesh.json" ).c_str() );
                    ifs >> j;
                }

#if defined( FEELPP_HAS_MMG ) && defined( FEELPP_HAS_PARMMG )
                std::string metr = soption( "remesh.metric" );
                LOG(INFO) << fmt::format( "metric:{}\n",soption( "remesh.metric" ));
                LOG(INFO) << fmt::format( "h:{}\n",soption( "remesh.h" ));
                auto [out,c] = remesh( _mesh = mesh, _metric = metr, _params = j );
                mesh=out;
#else
                fmt::print( fg( fmt::color::crimson ) | fmt::emphasis::bold,
                            "mmg and parmmg are not configured with feelpp!\n"
                            " - remeshing is disabled\n"
                            " - result mesh is the initial mesh\n" );
#endif
            }
        }
        dump(mesh);
#endif

        fs::path inputDir = inputPathMesh.parent_path();
        std::string inputFilenameWithoutExt = inputPathMesh.stem().string();

        std::vector<Range<mesh_type,MESH_ELEMENTS>> partitionByRange;
        if ( jPartitioner.contains("splitting") )
        {
            auto jSplitting = jPartitioner.at("splitting");

            std::map<int,std::set<std::string>> splittingMarkers;
            int splitId = 0;

            for ( auto const& [key,val] : jSplitting.items() )
            {
                std::set<std::string> markerNames;
                if ( val.is_array() )
                {
                    for ( auto const& [subkey,subval] : val.items() )
                        markerNames.insert( subval.template get<std::string>() );
                }
                else if ( val.is_string() )
                    markerNames.insert( val.template get<std::string>() );
                if ( !markerNames.empty() )
                    splittingMarkers.emplace( splitId, std::move( markerNames) );
                ++splitId;
            }

            partitionByRange.resize( jSplitting.size() );
            auto collectionMarkerElts = collectionOfMarkedelements( mesh, splittingMarkers );
            for ( auto const& [splitId,rangeElt] : collectionMarkerElts )
                partitionByRange[splitId] = rangeElt;
        }

        auto jOutput = partconfig.at( "output" );

        fs::path outputDir = jOutput.at("directory").template get<std::string>();
        if ( !fs::exists( outputDir ) )
            fs::create_directories( outputDir );

        std::string outputBaseFilename = inputFilenameWithoutExt;
        if ( jOutput.contains("filename") )
            outputBaseFilename = jOutput.at("filename").template get<std::string>();

        auto rangeMeshElt = elements( mesh );
        std::shared_ptr<exporter_t<mesh_type,mesh_type::nOrder>> exporter;

        using pid_functionspace_type = Pdh_type<mesh_type,0>;
        using pid_element_type = Pdh_element_t<mesh_type,0>;
        std::shared_ptr<pid_functionspace_type> Vh;
        std::shared_ptr<pid_element_type> pidField;
        if ( partconfig.contains( "visualization-exporter" ) )
        {
            auto jVisuExporter = partconfig.at( "visualization-exporter" );
            bool doExportVisu = jVisuExporter.value( "enabled", false );
            if ( doExportVisu )
            {
                Vh = Pdh<0>( mesh );
                pidField = Vh->elementPtr();
                exporter = Feel::exporter( _mesh = mesh,_geo="static" );
            }
        }

        for ( int nPartition : nParts )
        {
            std::string outputFilenameWithoutExt = fmt::format("{}_p{}",outputBaseFilename,nPartition);
            std::string outputFilenameWithExt = outputFilenameWithoutExt + ".json";
            std::string outputPathMesh = ( outputDir / fs::path(outputFilenameWithExt) ).string();

            std::cout << "start mesh partitioning and save on disk : " << outputPathMesh << std::endl;

            // build a MeshPartitionSet based on a mesh partition that will feed a
            // partition io data structure to generate a parallel hdf5 file from which
            // the parallel mesh can be loaded
            tic();
            using io_t = PartitionIO<mesh_t<decltype(mesh)>>;
            io_t io( outputPathMesh );
            io.write( partitionMesh( mesh, nPartition, partitionByRange, partconfig ) );
            toc("paritioning and save on disk done",FLAGS_v>0);

            if ( exporter )
            {
                for ( rank_type partId=0; partId<nPartition; ++partId)
                {
                    for ( auto const& eltWrap : elements(mesh,partId) )
                    {
                        auto const& elt = unwrap_ref( eltWrap );
                        for( auto const& ldof : Vh->dof()->localDof( elt.id() ) )
                            pidField->set( ldof.second.index(), partId );
                    }
                }
            }
            // WARNING! partitionMesh should not modify mesh -> TODO!
            for ( auto const& eltWrap : allelements(mesh) )
                const_cast<typename mesh_type::element_type&>(unwrap_ref( eltWrap )).setProcessId( 0 );

            if ( exporter )
            {
                exporter->step( nPartition )->add( "pid", *pidField );
                exporter->save();
            }
        }
    }

}



}
#endif
