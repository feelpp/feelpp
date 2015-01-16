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
   \file geo.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#if !defined(FEELPP_GEO_HPP)
#define FEELPP_GEO_HPP 1

#include <feel/feelfilters/gmsh.hpp>

namespace Feel {

/**
 * \brief geo return a gmsh_ptrtype of a .geo mesh
 *
 * \arg filename
 * \arg dimension
 * \arg order (optional, default = 1)
 * \arg h (optional, default = 0.1 )
 */
BOOST_PARAMETER_FUNCTION(
    ( gmsh_ptrtype ), // return type
    geo,    // 2. function name
    tag,           // 3. namespace of tag types
    ( required
      ( filename,       *( boost::is_convertible<mpl::_,std::string> ) ) )
    ( optional
      ( desc, *( boost::is_convertible<mpl::_,std::string> ), std::string() )
      ( h,              *( boost::is_arithmetic<mpl::_> ), doption(_name="gmsh.hsize") )
      ( geo_parameters,    *( boost::icl::is_map<mpl::_> ), Gmsh::gpstr2map( soption(_name="gmsh.geo-variables-list") ) )
      ( dim,              *( boost::is_integral<mpl::_> ), 3 )
      ( order,              *( boost::is_integral<mpl::_> ), 1 )
      ( files_path, *( boost::is_convertible<mpl::_,std::string> ), Environment::localGeoRepository() )
      ( depends, *( boost::is_convertible<mpl::_,std::string> ), soption(_name="gmsh.depends") )
      ( worldcomm,       (WorldComm), Environment::worldComm() ) )
    )

{
    gmsh_ptrtype gmsh_ptr( new Gmsh( 3, 1, worldcomm ) );

    gmsh_ptr->setCharacteristicLength( h );

#if BOOST_FILESYSTEM_VERSION == 3
        gmsh_ptr->setPrefix( fs::path( filename ).stem().string() );
#elif BOOST_FILESYSTEM_VERSION == 2
        gmsh_ptr->setPrefix( fs::path( filename ).stem() );
#endif

    if ( !desc.empty() )
    {
        gmsh_ptr->setDescription( desc );
    }
    else
    {
        std::string filename_with_path = Environment::findFile( filename );
        if ( filename_with_path.empty() )
        {
            std::vector<std::string> plist = Environment::geoPathList();
            std::ostringstream ostr;
            std::for_each( plist.begin(), plist.end(), [&ostr]( std::string s ) { ostr << " - " << s << "\n"; } );
            CHECK( !filename_with_path.empty() ) << "File " << filename << " cannot be found in the following paths list:\n " << ostr.str();
        }

        gmsh_ptr->setDescription( gmsh_ptr->getDescriptionFromFile( filename_with_path ) );

        if( worldcomm.globalRank() == worldcomm.masterRank() )
        {
            fs::path cp = fs::current_path();
            std::vector<std::string> depends_on_files;
            if ( !depends.empty() )
                algorithm::split( depends_on_files, depends, algorithm::is_any_of( ":,; " ), algorithm::token_compress_on );
            // copy include/merged files needed by geometry file
            boost::for_each( depends_on_files,
                             [&cp, &files_path]( std::string const& _filename )
                             {
                                 fs::path file_path( files_path );
                                 file_path /= _filename;

                                 try
                                 {
                                     boost::system::error_code ec;

                                     if ( !( fs::exists( file_path ) && fs::is_regular_file( file_path ) ) )
                                         std::cout << "File : " << file_path << " doesn't exist or is not a regular file" << std::endl;

                                     else if ( !fs::exists( cp / _filename )  )
                                         fs::copy_file( file_path, fs::path( _filename ), fs::copy_option::none );

                                 }

                                 catch ( const fs::filesystem_error& e )
                                 {
                                     std::cerr << "Error: " << e.what() << std::endl;
                                 }
                             } );
        }
    }
    gmsh_ptr->setGeoParameters( gmsh_ptr->retrieveGeoParameters( gmsh_ptr->description() ), 0 );
    gmsh_ptr->setGeoParameters( geo_parameters );

    return gmsh_ptr;

}


}

#endif /* FEELPP_GEO_HPP */
