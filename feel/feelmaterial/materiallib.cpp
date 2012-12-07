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
   \file materiallib.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-05-25
 */
#include <feel/feelmaterial/materiallib.hpp>

#include <feel/feelmaterial/air.hpp>
#include <feel/feelmaterial/castiron.hpp>


namespace Feel
{
/**
 * \return the command lines options of the material library
 */
po::options_description material_options( std::string const& prefix )
{
    std::string _prefix = prefix;

    if ( !_prefix.empty() )
        _prefix += "-";

    po::options_description _options( "Material " + prefix + "  options" );
    _options.add_options()
    // material library options
    //((_prefix+"material-lib").c_str(), Feel::po::value<std::string>()->default_value( "stdmaterial.so" ), "Standard material library")
    ( ( _prefix+"material" ).c_str(), Feel::po::value<std::string>()->default_value( "Air" ), "material" )
    ;
    return _options;
}

MaterialLib::MaterialLib()
{}

MaterialLib::MaterialLib( po::variables_map const& vm )
{}
MaterialLib::~MaterialLib()
{}

material_ptrtype
MaterialLib::material( std::string const& name )
{
    return material_ptrtype( factory_type::instance().createObject( name ) );
}


const bool material_air = MaterialLib::factory_type::instance().registerProduct( "Air", &detail::createMaterial<Air> );
const bool material_castiron = MaterialLib::factory_type::instance().registerProduct( "CastIron", &detail::createMaterial<CastIron> );


}
