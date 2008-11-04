/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-05-25

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file materiallib.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-05-25
 */
#include <life/lifematerial/materiallib.hpp>

namespace Life
{
/**
 * \return the command lines options of the material library
 */
po::options_description material_options( std::string const& prefix )
{
    std::string _prefix = prefix;
    if ( !_prefix.empty() )
        _prefix += "-";

    po::options_description _options( "Material " + prefix + "  options");
    _options.add_options()
        // material library options
        ((_prefix+"material-lib").c_str(), Life::po::value<std::string>()->default_value( "stdmaterial.so" ), "Standard material library")
        ((_prefix+"material").c_str(), Life::po::value<std::string>()->default_value( "Air" ), "material")
        ;
    return _options;
}

MaterialLib::MaterialLib()
    :
    M_dll( "/scratch/life/opt2-gcc43/testsuite/lifecore/stdmaterial.so" ),
    M_plugin_factory( M_dll )
{}

MaterialLib::MaterialLib( po::variables_map const& vm )
    :
    M_dll( vm["material-lib"].as<std::string>() ),
    M_plugin_factory( M_dll )
{}
MaterialLib::~MaterialLib()
{}

material_ptrtype
MaterialLib::material( std::string const& name )
{
    return material_ptrtype( M_plugin_factory.create ( name ) );
}
}
