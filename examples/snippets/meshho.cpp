/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2011-03-15

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
   \file dofpoints.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2011-03-15
 */
#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelvf/vf.hpp>

namespace Feel
{

}
int main( int argc, char** argv )
{
    double hsize = 2;
    // Declare the supported options.
    namespace po = boost::program_options;
    po::options_description desc( "Allowed options" );
    desc.add_options()
    ( "help", "produce help message" )
    ( "hsize", po::value<double>( &hsize )->default_value( 2 ), "h size" )
    ;

    po::variables_map vm;
    po::store( po::parse_command_line( argc, argv, desc ), vm );
    po::notify( vm );

    using namespace Feel;
    using namespace Feel::vf;
    Feel::Environment env( argc, argv );
    typedef Mesh<Simplex<2,4> > mesh_type;

    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _desc=domain( _name="ellipsoid-2",
                                        _usenames=true,
                                        _shape="ellipsoid",
                                        _dim=2,
                                        _order=4,
                                        _h=hsize ),
                                _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES );
    //straightenMesh( _mesh=mesh );

    //std::cout << "read mesh\n" << std::endl;

    std::cout << "ho  p1  opt " << std::endl
              << std::setprecision( 16 ) << std::scientific << integrate( _range=elements( mesh ), _quad=_Q<5>(), _expr=cst( 1. ), _geomap=GeomapStrategyType::GEOMAP_HO ).evaluate() << " "
              << std::setprecision( 16 ) << std::scientific << integrate( _range=elements( mesh ), _quad=_Q<5>(), _expr=cst( 1. ), _geomap=GeomapStrategyType::GEOMAP_O1 ).evaluate() << " "
              << std::setprecision( 16 ) << std::scientific << integrate( _range=elements( mesh ), _quad=_Q<5>(), _expr=cst( 1. ), _geomap=GeomapStrategyType::GEOMAP_OPT ).evaluate() << std::endl;


}
