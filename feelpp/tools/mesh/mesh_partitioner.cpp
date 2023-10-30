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
#include <mesh_partitioner.hpp>



namespace Feel {

extern template void partition<Simplex<1>>( std::vector<int> const& nParts, nl::json const& j );
extern template void partition<Simplex<2>>( std::vector<int> const& nParts, nl::json const& j );
extern template void partition<Simplex<2,2>>( std::vector<int> const& nParts, nl::json const& j );
extern template void partition<Simplex<3>>( std::vector<int> const& nParts, nl::json const& j );
extern template void partition<Simplex<3,2>>( std::vector<int> const& nParts, nl::json const& j );
extern template void partition<Hypercube<2>>( std::vector<int> const& nParts, nl::json const& j );
extern template void partition<Hypercube<2,2>>( std::vector<int> const& nParts, nl::json const& j );
extern template void partition<Hypercube<3>>( std::vector<int> const& nParts, nl::json const& j );
extern template void partition<Hypercube<3,3>>( std::vector<int> const& nParts, nl::json const& j );

}

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
        ( "by-markers", "partitioning by markers" )
        ( "by-markers-desc", po::value<std::vector<std::string> >()->multitoken(), "partitioning by markers. Example : --by-markers-desc=marker1:marker2,marker3" )
        ( "part", po::value<std::vector<int> >()->multitoken(), "number of partition" )
        ( "json", po::value<std::string>(), "json configuration file" )
        ( "ifile", po::value<std::string>(), "input mesh filename" )
        ( "odir", po::value<std::string>(), "output directory [optional]" )
        ( "ofile", po::value<std::string>(), "output mesh filename [optional]" )
        ( "remesh", po::value<bool>()->default_value( 0 ), "remesh " )
        ( "remesh.metric", po::value<std::string>()->default_value( "" ), "remesh metric expression" )
        ( "remesh.h", po::value<std::string>()->default_value( "hmax" ), "remesh h size" )
        ( "scalar_expr", po::value<std::vector<std::string>>()->default_value( {"g|sin(x):x|nodal|element"} ), "list of scalar expressions with name and representations" )
        ( "vectorial_expr", po::value<std::vector<std::string>>()->default_value( {"gv|{sin(2*pi*x),sin(2*pi*x),sin(2*pi*x)}:x|nodal|element"} ), "list of vectorial  expressions with name and representations" )

		;

    Environment env( _argc=argc, _argv=argv,
                     _desc=meshpartoptions,
                     _about=about( _name="mesh_partitioner" ,
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ),
                     _directory=".",_subdir=false );

    int dim = ioption(_name="dim");
    int realdim;
    if ( Environment::vm().count("realdim"))
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
    nl::json partconfig;
    if ( Environment::vm().count("json") && fs::exists( fs::current_path() / soption("json") ) )
    {
        Feel::cout << "reading partitioner configuration: " << ( fs::current_path() / soption("json") ) << std::endl;
        std::ifstream ifs( ( fs::current_path() / soption("json") ).string() );
        partconfig = nl::json::parse( ifs );
    }
    Feel::cout << "json config: " << partconfig.dump(1) << std::endl;

    fs::path pathInputMesh = fs::canonical( soption("ifile") );
    if ( !fs::exists( pathInputMesh ) )
    {
        if ( Environment::isMasterRank() )
            std::cout << "do nothing because --ifile : " << pathInputMesh.string() << " does not exist\n";
        return 0;
    }

    if ( dim == 1 )
    {
        partition<Simplex<1>>( nParts, partconfig );
    }
    else
    {
        if ( shape == "simplex" )
        {
            switch ( dim )
            {
            case 2 :
                if ( order == 1 && realdim==2)
                    partition<Simplex<2>>( nParts, partconfig );
                else if ( order == 1 && realdim==3)
                    partition<Simplex<2,1,3>>( nParts, partconfig );
                else if ( order==2 && realdim==2)
                    partition<Simplex<2,2>>( nParts, partconfig );
                else if ( order==2 && realdim==3)
                    partition<Simplex<2,2,3>>( nParts, partconfig );
                break;
            case 3 :
                if ( order == 1 )
                    partition<Simplex<3>>( nParts, partconfig );
                else if ( order==2 )
                    partition<Simplex<3,2>>( nParts, partconfig );
                break;
            }
        }
        else if ( shape == "hypercube" )
        {
            switch ( dim )
            {
            case 2 :
                if ( order == 1 )
                    partition<Hypercube<2>>( nParts, partconfig );
                else if ( order==2 )
                    partition<Hypercube<2,2>>( nParts, partconfig );
                break;
            case 3 :
                if ( order == 1 )
                    partition<Hypercube<3>>( nParts, partconfig );
                else if ( order==2 )
                    partition<Hypercube<3,2>>( nParts, partconfig );
                break;
            }
        }
    }
    return 0;

}
