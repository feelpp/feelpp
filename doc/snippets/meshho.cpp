/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-03-15
 */
#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/savegmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
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
    typedef Mesh<Simplex<3,2> > mesh_type;

    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _desc=domain( _name="ellipsoid-3",
                                        _usenames=true,
                                        _shape="ellipsoid",
                                        _h=hsize ),
                                _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES,
                                _straighten=0 );
    straightenMesh( mesh, Environment::worldComm(), false, true );

    //std::cout << "read mesh\n" << std::endl;

    std::cout << "ho  p1  opt " << std::endl
              << std::setprecision( 16 ) << std::scientific << integrate( _range=elements( mesh ), _quad=_Q<5>(), _expr=cst( 1. ), _geomap=GeomapStrategyType::GEOMAP_HO ).evaluate() << " "
              << std::setprecision( 16 ) << std::scientific << integrate( _range=elements( mesh ), _quad=_Q<5>(), _expr=cst( 1. ), _geomap=GeomapStrategyType::GEOMAP_O1 ).evaluate() << " "
              << std::setprecision( 16 ) << std::scientific << integrate( _range=elements( mesh ), _quad=_Q<5>(), _expr=cst( 1. ), _geomap=GeomapStrategyType::GEOMAP_OPT ).evaluate() << std::endl;


    std::cout << "ho  p1  opt " << std::endl
              << std::setprecision( 16 ) << std::scientific << integrate( _range=internalelements( mesh ), _quad=_Q<5>(), _expr=cst( 1. ), _geomap=GeomapStrategyType::GEOMAP_HO ).evaluate() << " "
              << std::setprecision( 16 ) << std::scientific << integrate( _range=internalelements( mesh ), _quad=_Q<5>(), _expr=cst( 1. ), _geomap=GeomapStrategyType::GEOMAP_O1 ).evaluate() << " "
              << std::setprecision( 16 ) << std::scientific << integrate( _range=internalelements( mesh ), _quad=_Q<5>(), _expr=cst( 1. ), _geomap=GeomapStrategyType::GEOMAP_OPT ).evaluate() << std::endl;

    std::cout << "ho  p1  opt " << std::endl
              << std::setprecision( 16 ) << std::scientific << integrate( _range=boundaryelements( mesh ), _quad=_Q<5>(), _expr=cst( 1. ), _geomap=GeomapStrategyType::GEOMAP_HO ).evaluate() << " "
              << std::setprecision( 16 ) << std::scientific << integrate( _range=boundaryelements( mesh ), _quad=_Q<5>(), _expr=cst( 1. ), _geomap=GeomapStrategyType::GEOMAP_O1 ).evaluate() << " "
              << std::setprecision( 16 ) << std::scientific << integrate( _range=boundaryelements( mesh ), _quad=_Q<5>(), _expr=cst( 1. ), _geomap=GeomapStrategyType::GEOMAP_OPT ).evaluate() << std::endl;


    auto range = internalelements(mesh);
    for( auto it = range.get<1>(), en = range.get<2>(); it != en; ++it )
    {
        double elt1 = integrate( _range=idedelements( mesh, it->id() ), _quad=_Q<5>(), _expr=cst( 1. ), _geomap=GeomapStrategyType::GEOMAP_HO ).evaluate()( 0, 0 );
        double elt2 = integrate( _range=idedelements( mesh, it->id() ), _quad=_Q<5>(), _expr=cst( 1. ), _geomap=GeomapStrategyType::GEOMAP_OPT ).evaluate()( 0, 0 );
        double elt3 = integrate( _range=idedelements( mesh, it->id() ), _quad=_Q<5>(), _expr=cst( 1. ), _geomap=GeomapStrategyType::GEOMAP_O1 ).evaluate()( 0, 0 );
        if ( math::abs( elt1-elt2 ) > 1e-13 )
        {
            std::cout << "problem with element: " << it->id() << " elt1= " << elt1 << " elt2 = " << elt2 << " elt3 = " << elt3 << "\n";
            auto meshe = createSubmesh( mesh, boost::make_tuple( mpl::int_<MESH_ELEMENTS>(), it, boost::next( it ) ) );
            std::ostringstream os;
            os << "elt-ho-" << it->id();
            saveGMSHMesh( _mesh=meshe, _filename=os.str() );
            //std::cout << "G=" << it->G() << "\n";
            //std::cout << "V=" << it->vertices() << "\n";
        }
    }

}
