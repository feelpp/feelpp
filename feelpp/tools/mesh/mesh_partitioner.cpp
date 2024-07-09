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


#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <mesh_partitioner.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;
    using Feel::cout;

    po::options_description meshpartoptions( "Mesh Partitioner options" );
    meshpartoptions.add_options()
        ( "dim", po::value<int>()->default_value( 3 ), "mesh dimension" )
        ( "realdim", po::value<int>(), "real dimension of mesh nodes" )
        ( "shape", po::value<std::string>()->default_value( "simplex" ), "mesh basic unit" )
        ( "order", po::value<int>()->default_value( 1 ), "mesh geometric order" )
        ( "part", po::value<std::vector<int> >()->multitoken(), "number of partition" )
        ( "json", po::value<std::string>(), "json configuration file" )
        ( "splitting", po::value<std::string>(), "define splitting of partitioner, i.e. create partitioning on each split" )
        ( "ifile", po::value<std::string>(), "input mesh filename" )
        ( "odir", po::value<std::string>(), "output directory [optional]" )
        ( "ofile", po::value<std::string>(), "output base filename (name without extension and part number suffix) [optional]" )
        // ( "remesh", po::value<bool>()->default_value( 0 ), "remesh " )
        // ( "remesh.metric", po::value<std::string>()->default_value( "" ), "remesh metric expression" )
        // ( "remesh.h", po::value<std::string>()->default_value( "hmax" ), "remesh h size" )
        // ( "scalar_expr", po::value<std::vector<std::string>>()->default_value( {"g|sin(x):x|nodal|element"} ), "list of scalar expressions with name and representations" )
        // ( "vectorial_expr", po::value<std::vector<std::string>>()->default_value( {"gv|{sin(2*pi*x),sin(2*pi*x),sin(2*pi*x)}:x|nodal|element"} ), "list of vectorial  expressions with name and representations" )
        ( "export-visualization", po::value<bool>(), "export mesh with partitioning(s)" )
        ;

    auto initialCurrentPath = fs::current_path();

    Environment env( _argc=argc, _argv=argv,
                     _desc=meshpartoptions,
                     _about=about( _name="mesh_partitioner" ,
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ),
                     _directory=".",_subdir=false );

    nl::json partconfig;
    if ( Environment::vm().count("json") )
    {
        fs::path jsonFilename = Environment::expand( soption(_name="json") );
        if ( jsonFilename.is_relative() )
            jsonFilename = initialCurrentPath/jsonFilename;
        if ( fs::exists( jsonFilename ) )
        {
            Feel::cout << "reading partitioner configuration: " << jsonFilename << std::endl;
            std::ifstream ifs( jsonFilename );
            partconfig = nl::json::parse( ifs, nullptr, true, true );
        }
    }

    int dim = ioption(_name="dim");
    int realdim;
    if ( Environment::vm().count("realdim") )
    {
        realdim = ioption(_name="realdim");
        CHECK( realdim >= dim && realdim <=3 ) << "Realdim must be larger thank dim and at most 3";
    }
    else
    {
        realdim=dim;
    }
    std::string shape = soption(_name="shape");
    int order = ioption(_name="order");

    auto & jPartitioner = partconfig["partitioner"];

    if ( Environment::vm().count("part") )
    {
        std::vector<int> nParts = Environment::vm()["part"].as<std::vector<int> >();
        jPartitioner[ "number-of-partition" ] = nParts;
    }
    if ( Environment::vm().count("splitting") )
    {
        std::istringstream iss( soption(_name="splitting") );
        nl::json jSplittingVal = nl::json::parse(iss);
        jPartitioner["splitting"] = jSplittingVal;
    }

    if ( Environment::vm().count("ifile") )
    {
        fs::path ifileOption = Environment::expand( soption("ifile") );
        if ( ifileOption.is_relative() )
            ifileOption = initialCurrentPath/ifileOption;
        auto & jInput = partconfig["input"];
        jInput[ "filename" ] = fs::canonical( ifileOption ).string();
    }

    if ( Environment::vm().count("odir") )
    {
        fs::path odirOption = Environment::expand( soption("odir") );
        if ( odirOption.is_relative() )
            odirOption = initialCurrentPath/odirOption;
        auto & jOutput = partconfig["output"];
        jOutput[ "directory" ] = odirOption.string();
    }
    if ( Environment::vm().count("ofile") )
    {
        auto & jOutput = partconfig["output"];
        jOutput[ "filename" ] = soption(_name="ofile");
    }

    if ( Environment::vm().count("export-visualization") )
    {
        auto & jVisuExporter = partconfig["visualization-exporter"];
        jVisuExporter[ "enabled" ] = boption("export-visualization");
    }

    Feel::cout << "json config: " << partconfig.dump(1) << std::endl;

    if ( !partconfig.contains("input") )
    {
        if ( Environment::isMasterRank() )
            std::cout << "do nothing because no input file" << std::endl;
        return 0;
    }
    if ( !partconfig.contains("partitioner") )
    {
        if ( Environment::isMasterRank() )
            std::cout << "do nothing because partition setup is missing" << std::endl;
        return 0;
    }



    if ( dim == 1 )
    {
        partition<Simplex<1>>( partconfig );
    }
    else
    {
        if ( shape == "simplex" )
        {
            switch ( dim )
            {
            case 2 :
                if ( order == 1 && realdim==2)
                    partition<Simplex<2>>( partconfig );
                else if ( order == 1 && realdim==3)
                    partition<Simplex<2,1,3>>( partconfig );
                else if ( order==2 && realdim==2)
                    partition<Simplex<2,2>>( partconfig );
                else if ( order==2 && realdim==3)
                    partition<Simplex<2,2,3>>( partconfig );
                break;
            case 3 :
                if ( order == 1 )
                    partition<Simplex<3>>( partconfig );
                else if ( order==2 )
                    partition<Simplex<3,2>>( partconfig );
                break;
            }
        }
        else if ( shape == "hypercube" )
        {
            switch ( dim )
            {
            case 2 :
                if ( order == 1 )
                    partition<Hypercube<2>>( partconfig );
                else if ( order==2 )
                    partition<Hypercube<2,2>>( partconfig );
                break;
            case 3 :
                if ( order == 1 )
                    partition<Hypercube<3>>( partconfig );
                else if ( order==2 )
                    partition<Hypercube<3,2>>( partconfig );
                break;
            }
        }
    }

    return 0;

}
