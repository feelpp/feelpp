/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 22 Nov 2019

 Copyright (C) 2019 Feel++ Consortium

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
#include <map>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/ginac.hpp>
#include <feel/feelvf/geometricdata.hpp>
#include <boost/hana/tuple.hpp>
#include <boost/hana/pair.hpp>

template<typename MeshT>
void printMeshInfo( std::shared_ptr<MeshT> const& mesh )
{
    using namespace Feel;
    std::cout << "========================================" << std::endl;
    std::cout << "Information about the mesh:";
    std::cout << "      number of elements in memory : " << mesh->numGlobalElements() << std::endl;
    std::cout << "      number of faces in memory : " << mesh->numGlobalFaces() << std::endl;
    if ( dimension_v<MeshT> == 3 )
        std::cout << "      number of edges in memory : " << mesh->numGlobalEdges() << std::endl;
    std::cout << "      number of points  in memory : " << mesh->numGlobalPoints() << std::endl;
    for( auto marker: mesh->markerNames() )
    {
        auto [name,data] = marker;

        if ( data[1] == mesh->dimension() )
        {
            size_type nelts = nelements( markedelements(mesh, name ), true );
            std::cout << "      number of marked elements " << name << " with tag " << data[0] << " : " << nelts << std::endl;
        }
    }
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
    std::cout << "========================================" << std::endl;
}

template<typename mesh_t>
void
doExport()
{
    using namespace Feel;
    auto mesh = loadMesh( _mesh = new mesh_t );
    printMeshInfo( mesh );
    auto ex = exporter( _mesh=mesh );

    auto getfns = [=]( std::string const& sexpr ) {
                      std::vector<std::string> fields;
                      boost::split( fields, sexpr, boost::is_any_of("|"), boost::token_compress_on );
                      CHECK( fields.size()  > 1 ) << "bad expression format";
                      std::set<std::string> reps;
                      for( auto it=fields.begin()+2; it!=fields.end(); ++it )
                          reps.insert( *it );
                      return std::tuple{ fields[0], fields[1], reps };
                  };
    
    for( auto const& sexpr : vsoption( _name="scalar_expr" ) )
    {
        auto const& [strname,strexpr, reps ] = getfns( sexpr );
        ex->add( strname, expr( strexpr ), reps );
    }
    for( auto const& sexpr : vsoption( _name="vectorial_expr" ) )
    {
        auto const& [strname,strexpr, reps ] = getfns( sexpr );
        ex->add( strname, expr<3,1>( strexpr ), reps );
    } 
    ex->add( "P", P(), "nodal" );
    ex->add( "facesmarker", semarker(), faces(mesh), "nodal" );
    if constexpr ( dimension_v<mesh_t> > 2 )
        ex->add( "edgesmarker", semarker(), edges(mesh), "nodal" );
    ex->add( "pointmarker", semarker(), points(mesh), "nodal" );

    using namespace std::string_literals;
    auto exprs = hana::make_tuple( hana::pair{ "detJ"s, detJ() }, hana::pair{ "marker"s, emarker() }, hana::pair{ "pid"s, epid() }, hana::pair{ "h"s, h() }  );
    hana::for_each( exprs,
                    [&](const auto& x) { ex->add( hana::first(x), hana::second(x), "element" ); } );

    
    ex->save();

}

int main( int argc, char** argv )
{
    using namespace Feel;
    using Feel::cout;
    po::options_description meshpartoptions( "Mesh Partitioner options" );
	meshpartoptions.add_options()
        ( "dim", po::value<int>()->default_value( 3 ), "mesh dimension" )
        ( "shape", po::value<std::string>()->default_value( "simplex" ), "mesh dimension" )
        ( "order", po::value<int>()->default_value( 1 ), "mesh geometric order" )
        ( "scalar_expr", po::value<std::vector<std::string>>()->default_value( {"g|sin(x):x|nodal|element"} ), "list of scalar expressions with name and representations" )
        ( "vectorial_expr", po::value<std::vector<std::string>>()->default_value( {"gv|{sin(2*pi*x),sin(2*pi*x),sin(2*pi*x)}:x|nodal|element"} ), "list of vectorial  expressions with name and representations" )

        ;
    Environment env( _argc=argc, _argv=argv,
                     _desc=meshpartoptions,
                     _about=about( _name="mesh_exporter" ,
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" )
                     );


    if ( ioption( "dim" ) == 1 )
        doExport<Mesh<Simplex<1, 1>>>();
    if ( ioption( "dim" ) == 2 )
        doExport<Mesh<Simplex<2, 1>>>();
    if ( ioption( "dim" ) == 3 )
        doExport<Mesh<Simplex<3, 1>>>();
    return 0;

}

