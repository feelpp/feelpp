/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-07-25

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
   \file feel.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-07-25
 */
#include <feel/feelcore/feel.hpp>
#include <boost/algorithm/string/trim_all.hpp>
namespace Feel
{

std::string
prefixvm( std::string const& prefix,
          std::string const& opt,
          std::string const& sep )
{
    std::string o = prefix;

    if ( !o.empty() )
        o += sep;

    return o+opt;
}

/**
 * @note Ensight variable description forbids ( [ + @ ! * $ ) ] - space # ^ /)
 * see ensight gold file format documentation
 */
std::string
sanitize( std::string const& s )
{
    
    
    return algorithm::trim_fill_copy_if(s, "_", algorithm::is_any_of(" ;*,:()[]@$/+-#^") );
}

std::vector<std::string>
sanitize( std::vector<std::string> const& s )
{
    std::vector<std::string> ss(s);
    for_each(ss.begin(), ss.end(), []( std::string& n ) { algorithm::trim_fill_if(n, "_", algorithm::is_any_of(" ;*,:()[]@$/+-#^") ); } );
    return ss;
}


}
