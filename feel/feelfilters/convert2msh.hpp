/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-12-24

  Copyright (C) 2013 Feel++ Consortium

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
   \file convert2msh.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#if !defined(FEELPP_CONVERT2MSH_HPP)
#define FEELPP_CONVERT2MSH_HPP 1

#include <feel/feelfilters/gmsh.hpp>

namespace Feel {

/**
 * \brief convert to gmsh format
 *
 * \arg filename
 * \arg dim (optional, default = 3)
 * \arg order (optional, default = 1)
 */
BOOST_PARAMETER_FUNCTION(
    ( gmsh_ptrtype ), // return type
    convert2msh,    // 2. function name
    tag,           // 3. namespace of tag types
    ( required
      ( filename,       *( boost::is_convertible<mpl::_,std::string> ) ) )
    ( optional
      ( dim,              *( boost::is_integral<mpl::_> ), 3 )
      ( order,              *( boost::is_integral<mpl::_> ), 1 ) )
    )
{
    gmsh_ptrtype gmsh_ptr( new Gmsh( 3, 1 ) );
#if BOOST_FILESYSTEM_VERSION == 3
    gmsh_ptr->setPrefix( fs::path( filename ).stem().string() );
#elif BOOST_FILESYSTEM_VERSION == 2
    gmsh_ptr->setPrefix( fs::path( filename ).stem() );
#endif

    // first try in the current path
    if ( fs::exists( filename ) )
        gmsh_ptr->setDescription( ( boost::format( "Merge \"%1%\";\n" ) % filename ).str() );

    else if ( fs::exists( fs::path( Environment::localGeoRepository() ) / filename ) )
        gmsh_ptr->setDescription( ( boost::format( "Merge \"%1%\";\n" ) % ( fs::path( Environment::localGeoRepository() ) / filename ).string() ).str() );

    else if ( Environment::systemGeoRepository().template get<1>()  &&
              fs::exists( fs::path( Environment::systemGeoRepository().get<0>() ) / filename ) )
        gmsh_ptr->setDescription( ( boost::format( "Merge \"%1%\";\n" ) % ( fs::path( Environment::systemGeoRepository().get<0>() ) / filename ).string() ).str() );

    else
    {
        std::ostringstream ostr;
        ostr << "File " << filename << " was not found neither in current directory or in " << Environment::localGeoRepository() << " or in " << Environment::systemGeoRepository();
        throw std::invalid_argument( ostr.str() );
    }

    return gmsh_ptr;
}

}

#endif /* FEELPP_CONVERT2MSH_HPP */
