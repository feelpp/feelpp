/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-04-04

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
/**
   \file geotool.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-04-04
 */

#include <feel/feel.hpp>

template<typename MeshType>
void myexport( std::string const& name, boost::shared_ptr<MeshType> mesh )
{
    using namespace Feel;

    typedef Mesh<Simplex<2> > mesh_type;
    typedef Exporter<mesh_type> exporter_type;
    auto exporter = exporter_type::New( "gmsh", name );
    exporter->step( 0 )->setMesh( mesh );
    exporter->save();
}

int main( int argc, char**argv )
{
    using namespace Feel;

    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="geotool",
                                  _author="Christophe Prud'homme",
                                  _email="christophe.prudhomme@feelpp.org") );


    typedef Mesh<Simplex<2> > mesh_type;
    typedef Exporter<mesh_type> exporter_type;

    GeoTool::Rectangle R1( 0.1,"R1",GeoTool::Node( 0,0 ),GeoTool::Node( 1,1 ) );

    GeoTool::Circle C1( 0.1,"C1",GeoTool::Node( 0.5,0.5 ),GeoTool::Node( 0.75,0.75 ) );

    auto R1mesh = R1.createMesh(_mesh=new mesh_type,_name="R1" );
    auto C1mesh = C1.createMesh(_mesh=new mesh_type,_name="C1" );
    auto R1mC1mesh = ( R1-C1 ).createMesh(_mesh=new mesh_type,_name="R1-C1" );
    auto R1pC1mesh = ( R1+C1 ).createMesh(_mesh=new mesh_type,_name="R1+C1" );


    myexport<mesh_type>( "R1",R1mesh );
    myexport<mesh_type>( "C1",C1mesh );
    myexport<mesh_type>( "R1-C1",R1mC1mesh );
    myexport<mesh_type>( "R1+C1",R1pC1mesh );



}
