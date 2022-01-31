/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-07-20

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
   \file opuscomponent.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-07-20
 */
#include <string>
#include <feel/feelcore/feel.hpp>

#include "opuscomponent.hpp"

namespace Feel
{
po::options_description
makeComponentOptions()
{
    Feel::po::options_description pcboptions( "PCB Component" );
    pcboptions.add_options()
    ( "pcb.k"   , Feel::po::value<double>()->default_value( 0.2 ),     "Conductivity" )
    ( "pcb.rhoC", Feel::po::value<double>()->default_value( 2*1e6 ),  "Heat capacity" )
    ( "pcb.Q",    Feel::po::value<double>()->default_value( 0 ),       "Heat dissipation" )
    ( "pcb.h",    Feel::po::value<double>()->default_value( 13*1e-2 ), "height" )
    ( "pcb.e",    Feel::po::value<double>()->default_value( 2*1e-3 ),  "width" );
    Feel::po::options_description airoptions( "AIR Component" );
    airoptions.add_options()
    ( "air.k"   , Feel::po::value<double>()->default_value( 3*1e-2 ),     "Conductivity" )
    ( "air.rhoC", Feel::po::value<double>()->default_value( 1100 ),       "Heat capacity" )
    ( "air.Q",    Feel::po::value<double>()->default_value( 0 ),          "Heat dissipation" )
    ( "air.h",    Feel::po::value<double>()->default_value( 13*1e-2 ),    "height" )
    ( "air.e",    Feel::po::value<double>()->default_value( 5*1e-2 ),     "width" )
    ( "air.flow-rate",    Feel::po::value<double>()->default_value( 1e-3 ),     "Flow rate" );
    Feel::po::options_description ic1options( "IC1 Component" );
    ic1options.add_options()
    ( "ic1.k"   , Feel::po::value<double>()->default_value( 2 ),        "Conductivity" )
    ( "ic1.rhoC", Feel::po::value<double>()->default_value( 1.4*1e6 ),    "Heat capacity" )
    ( "ic1.Q",    Feel::po::value<double>()->default_value( 1e6 ),        "Heat dissipation" )
    ( "ic1.h",    Feel::po::value<double>()->default_value( 2*1e-2 ),     "height" )
    ( "ic1.e",    Feel::po::value<double>()->default_value( 2*1e-3 ),     "width" );

    Feel::po::options_description ic2options( "IC2 Component" );
    ic2options.add_options()
    ( "ic2.k"   , Feel::po::value<double>()->default_value( 2 ),         "Conductivity" )
    ( "ic2.rhoC", Feel::po::value<double>()->default_value( 1.4*1e6 ),    "Heat capacity" )
    ( "ic2.Q",    Feel::po::value<double>()->default_value( 1e6 ),        "Heat dissipation" )
    ( "ic2.h",    Feel::po::value<double>()->default_value( 2*1e-2 ),     "height" )
    ( "ic2.e",    Feel::po::value<double>()->default_value( 2*1e-3 ),     "width" );
    return pcboptions.add( airoptions ).add( ic1options ).add( ic2options );
}

#if 0
void
OpusComponent::save( boost::property_tree::ptree& pt )
{
    pt.put( "opus.eads.component." + M_name + ".name", M_name );
    pt.put( "opus.eads.component." + M_name + ".k", M_k );
    pt.put( "opus.eads.component." + M_name + ".rhoC", M_rhoC );
    pt.put( "opus.eads.component." + M_name + ".Q", M_Q );
    pt.put( "opus.eads.component." + M_name + ".h", M_h );
    pt.put( "opus.eads.component." + M_name + ".e", M_e );
}

void
OpusComponent::load( const boost::property_tree::ptree& pt )
{
    M_name = pt.get<std::string>( "opus.eads.component." + M_name + ".name", M_name );
    M_k = pt.get<double>( "opus.eads.component." + M_name + ".k", M_k );
    M_rhoC = pt.get<double>( "opus.eads.component." + M_name + ".rhoC", M_rhoC );
    M_Q = pt.get<double>( "opus.eads.component." + M_name + ".k", M_Q );
    M_h = pt.get<double>( "opus.eads.component." + M_name + ".h", M_h );
    M_e = pt.get<double>( "opus.eads.component." + M_name + ".e", M_e );
}
#endif // 0

std::ostream&
operator<<( std::ostream& os, OpusComponent const& oc )
{
    os << "Component : " << oc.name() << "\n"
       << "  - rhoC = " << oc.rhoC() << "\n"
       << "  - k    = " << oc.k() << "\n"
       << "  - Q    = " << oc.Q() << "\n"
       << "  - h    = " << oc.h() << "\n"
       << "  - e    = " << oc.e() << "\n";

    return os;
}
}
