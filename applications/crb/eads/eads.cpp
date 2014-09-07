/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-06-03

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
   \file eads.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-06-03
 */
#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feelalg/solvereigen.hpp>
#include <eads.hpp>
#include <opusdata.hpp>



namespace Feel
{
AboutData
makeEadsAbout( std::string const& str )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "2D OPUS/EADS Benchmark",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2008-2010 Université de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Sylvain Vallaghe", "developer", "sylvain.vallaghe@ujf-grenoble.fr", "" );
    about.addAuthor( "Stéphane Veys", "developer", "stephane.veys@ujf-grenoble.fr", "" );

    return about;

}

Feel::po::options_description opusModelThermalOptions();
Feel::po::options_description opusModelFluidOptions();

Feel::po::options_description
opusModelFluidOptions()
{
    Feel::po::options_description fluidoptions( "Fluid Model Options" );
    fluidoptions.add_options()
    ( "fluid.hsize"   , Feel::po::value<double>()->default_value( 0.5 ),     "Fluid mesh size" )
    ( "fluid.flow-rate"   , Feel::po::value<double>()->default_value( 5e-3 ),     "Fluid flow rate" );
    fluidoptions
    .add( backend_options( "fluid" ) );
    return fluidoptions;
}

Feel::po::options_description
opusModelThermalOptions()
{
    Feel::po::options_description thermaloptions( "Thermal Model Options" );
    thermaloptions.add_options()
    ( "thermal.T0"   , Feel::po::value<double>()->default_value( 300 ),       "initial and at Gamma_4 boundary temperature" )
    ( "thermal.c"   , Feel::po::value<double>()->default_value( 100 ),       "thermal conductance (1e30 means no resistance)" )
    ( "thermal.gamma-temp", Feel::po::value<double>()->default_value( 2.5e-2 ), "stabilization parameter" )
    ( "thermal.stab", Feel::po::value<bool>()->default_value( 0 ),              "0=no stabilisation, 1=CIP stabilization" )
    ( "thermal.distance-disc", Feel::po::value<double>()->default_value( 1e-6 ), "distance to temperature discontinuity" );

    thermaloptions
    .add( backend_options( "thermal" ) );

#if defined( FEELPP_HAS_SLEPC )
    thermaloptions
    .add( solvereigen_options( "thermal-coerc" ) )
    .add( solvereigen_options( "thermal-cont" ) );
#endif

    return thermaloptions;
}

Feel::po::options_description makeEadsOptions()
{
    Feel::po::options_description options( "Eads Model Options" );
#if 1
    options
    .add( opusModelThermalOptions() )
    .add( opusModelFluidOptions() );

#endif
    return options.add( backend_options( "backend.crb.fem" ) )
           .add( backend_options( "backend.crb.norm" ) )
           //.add( Feel::bdf_options() )
           .add( Feel::OpusData::makeOptions() )
           .add( Feel::feel_options() );
}

}
