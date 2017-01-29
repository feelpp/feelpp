/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 23 Feb 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
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


namespace Feel
{
po::options_description
stvenantkirchhoff_options(std::string prefix)
{
    po::options_description stvenantkirchhoffoptions( "StVenantKirchhoff problem options" );
    stvenantkirchhoffoptions.add_options()
        ( prefixvm( prefix, "young-modulus"), Feel::po::value<double>()->default_value( 1.4e6 ), "young-modulus" )
        ( prefixvm( prefix, "poisson-coeff"), Feel::po::value<double>()->default_value( 0.4 ), "poisson-coeff" )
        ( prefixvm( prefix, "rho"), Feel::po::value<double>()->default_value( 1000 ), "density [kg/m^3]" )
        ( prefixvm( prefix, "gravity"), Feel::po::value<std::string>()->default_value( "{0,0}" ), "gravity force expression" )
        ( prefixvm( prefix, "gravity-cst"), Feel::po::value<double>()->default_value( 2 ), "gravity-cst" )
        } 
}
