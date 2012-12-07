/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-05-25

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file stvenant_kirchhoff_base.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-05-25
 */
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>

#include <feel/options.hpp>
#include <feel/feelmesh/geoentity.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>

#include <structurebase.hpp>
#include <stvenant_kirchhoff.hpp>

namespace Feel
{
StructureBase::structure_ptrtype
StructureBase::New( Feel::po::variables_map const& vm )
{
    LOG(INFO) << "Creating new structure model and solver\n";
    using namespace Feel;

    if ( vm["d"].as<int>() == 2 )
    {
        if ( vm["sorder"].as<int>() == 3 )
        {
            return structure_ptrtype( new StVenantKirchhoff<2,3>( vm ) );
        }

        else if ( vm["sorder"].as<int>() == 8 )
        {
            return structure_ptrtype( new StVenantKirchhoff<2,8>( vm ) );
        }
    }

    else if ( vm["d"].as<int>() == 3 )
    {

    }
}

Feel::po::options_description
StructureBase::makeOptions()
{
    Feel::po::options_description structureoptions( "Structure benchmark options" );
    structureoptions.add_options()
    ( "d", Feel::po::value<int>()->default_value( 2 ), "time step value" )
    ( "dt", Feel::po::value<double>()->default_value( 1 ), "time step value" )
    ( "nsubdt", Feel::po::value<int>()->default_value( 2 ), "number of sub time steps to save" )
    ( "ft", Feel::po::value<double>()->default_value( 1 ), "final time value" )

    ( "sorder", Feel::po::value<int>()->default_value( 8 ), "order of space discretisation for the displacement" )
    ( "torder", Feel::po::value<int>()->default_value( 2 ), "order of time discretisation" )

    ( "gammabc", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter for weak Dirichlet condition" )

    ( "hsize", Feel::po::value<double>()->default_value( 0.5 ), "first h value to start convergence" )
    ( "mesh-type", Feel::po::value<int>()->default_value( 1 ), "0 = oplagp1, 1 = Xh mesh" )
    ( "export", "export results(ensight, data file(1D)" )
    ( "export-matlab", "export matrix and vectors in matlab" )
    ;
    return structureoptions.add( Feel::feel_options() );
}

Feel::AboutData
StructureBase::makeAbout()
{
    Feel::AboutData about( "structure" ,
                           "structure" ,
                           "0.1",
                           "nD(n=2,3) structure  benchmark",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2008 Université Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}
StructureBase::StructureBase( int d )
    :
    M_dimension( d ),
    M_h( 0.1 ),
    M_T0( 0 ),
    M_T( 1 ),
    M_dt( 0.1 ),

    M_sorder( 2 ),
    M_torder( 2 ),

    M_dirichlet(),
    M_neumann()
{

    M_dirichlet.push_back( "left" );
    M_neumann.push_back( "bottom" );
    print();
}

StructureBase::StructureBase( StructureBase const& data )
    :
    M_dimension( data.M_dimension ),

    M_h( data.M_h ),

    M_T0( data.M_T0 ),
    M_T( data.M_T ),
    M_dt( data.M_dt ),

    M_sorder( data.M_sorder ),
    M_torder( data.M_torder ),

    M_dirichlet(),
    M_neumann()
{
    M_dirichlet.push_back( "left" );
    M_neumann.push_back( "bottom" );
    print();
}
StructureBase::StructureBase( int d, Feel::po::variables_map const& vm )
    :
    M_vm ( vm ),
    M_dirichlet(),
    M_neumann()
{

    M_dimension = d;//vm["dim"].as<int>();
    M_h = vm["hsize"].as<double>();
    M_dt = vm["dt"].as<double>();
    M_T = vm["ft"].as<double>();
    //M_T0 = vm["it"].as<double>();
    M_nsubsteps = vm["nsubdt"].as<int>();

    M_sorder = vm["sorder"].as<int>();
    M_torder = vm["torder"].as<int>();

    M_gammabc = vm["gammabc"].as<double>();

    M_dirichlet.push_back( "left" );
    M_neumann.push_back( "bottom" );
    print();
}
StructureBase::~StructureBase()
{}
void
StructureBase::print() const
{
    LOG(INFO) << "dt = " << this->dt() << "\n";
    LOG(INFO) << " T = " << this->T() << "\n";
    LOG(INFO) << "T0 = " << this->T0() << "\n";

    LOG(INFO) << "order in space = " << this->spaceOrder() << "\n";
    LOG(INFO) << "order in time = " << this->timeOrder() << "\n";

    LOG(INFO) << "gammabc = " << this->gammaBc() << "\n";



}

std::string
StructureBase::createMesh()
{
    std::ostringstream ostr;
    std::ostringstream nameStr;

    switch ( M_dimension )
    {
    case 2:
        ostr << "h=" << M_h << ";\n"
             << "\n"
             << "H=1;\n"
             << "L=10;\n"
             << "\n"
             << "Point(1) = {0,0,0,h};\n"
             << "Point(2) = {L,0,0,h};\n"
             << "Point(3) = {L,H,0,h};\n"
             << "Point(4) = {0,H,0,h};\n"
             << "\n"
             << "Line(1) = {1,2};\n"
             << "Line(2) = {2,3};\n"
             << "Line(3) = {3,4};\n"
             << "Line(4) = {4,1};\n"
             << "\n"
             << "Line Loop(1)={1,2,3,4};\n"
             << "\n"
             << "Plane Surface(20) = {1};\n"
             << "\n"
             << "Physical Line(\"left\") = {4};\n"
             << "Physical Line(\"bottom\") = {1};\n"
             << "Physical Line(\"top\") = {3};\n"
             << "Physical Line(\"right\") = {2};\n"
             << "\n"
             << "Physical Surface( 1 )={20};\n";

        nameStr << "structure";
        break;

    case 3:
        break;

    default:
        std::ostringstream os;
        os << "invalid dimension: " << M_dimension;
        throw std::logic_error( os.str() );
    }

    Gmsh gmsh;
    gmsh.setOrder( 1 );
    gmsh.setVersion( 2 );
    return gmsh.generate( nameStr.str(), ostr.str() );
}

#if 0
StructureBase::Inflow::Inflow( StructureBase const& data )
    :
    M_data( data )
{}


double
StructureBase::Inflow::operator()( uint16_type c1, uint16_type c2, node_type const& p, node_type const& /*n*/ ) const
{
    if ( M_data.d() == 2 )
    {
        switch ( c1 )
        {
        case 0:
            return 4*M_data.Um()*p[1]*( M_data.H()-p[1] )/pow( M_data.H(),2 );

        case 1:
        default:
            return 0;
        }
    }

    else
    {

        // 3D
        switch ( c1 )
        {
        case 0:
            return 16*M_data.Um()*p[1]*( M_data.H()-p[1] )*p[2]*( M_data.H()-p[2] )/pow( M_data.H(),4 );

        case 1:
        case 2:
        default:
            return 0;
        }
    }
}
#endif

} // Feel

