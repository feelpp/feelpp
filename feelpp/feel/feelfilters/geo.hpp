/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-12-24

  Copyright (C) 2013-2016 Feel++ Consortium

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
#ifndef FEELPP_FILTERS_GEO_H
#define FEELPP_FILTERS_GEO_H

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
template <typename ... Ts>
gmsh_ptrtype geo( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    std::string const& filename = args.get(_filename );
    std::string const& prefix = args.get_else(_prefix,"" );
    po::variables_map const& vm = args.get_else( _vm, Environment::vm() );
    std::string const& desc = args.get_else( _desc, "" );
    double h = args.get_else_invocable( _h, [&prefix,&vm](){ return doption(_prefix=prefix,_name="gmsh.hsize",_vm=vm); } );
    auto && geo_parameters = args.get_else_invocable( _geo_parameters, [&prefix,&vm](){ return Gmsh::gpstr2map( soption(_prefix=prefix,_name="gmsh.geo-variables-list",_vm=vm) ); } );
    int dim = args.template get_else<std::is_integral>( _dim, 3 );
    int order = args.template get_else<std::is_integral>( _order, 1 );
    auto && files_path = args.get_else_invocable( _files_path, [](){ return Environment::localGeoRepository(); } );
    std::string const& depends = args.get_else_invocable(_depends, [&prefix,&vm](){ return soption(_prefix=prefix,_name="gmsh.depends",_vm=vm); } );
    worldcomm_ptr_t worldcomm = args.get_else(_worldcomm, Environment::worldCommPtr() );

    gmsh_ptrtype gmsh_ptr( new Gmsh( 3, 1, worldcomm ) );

    gmsh_ptr->setCharacteristicLength( h );

    gmsh_ptr->setPrefix( fs::path( filename ).stem().string() );

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

        if( worldcomm->globalRank() == worldcomm->masterRank() )
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
                                     if ( !( fs::exists( file_path ) && fs::is_regular_file( file_path ) ) )
                                         std::cout << "File : " << file_path << " doesn't exist or is not a regular file" << std::endl;

                                     else if ( !fs::exists( cp / _filename )  )
                                         fs::copy_file( file_path, fs::path( _filename ), fs::copy_options::none );
                                 }

                                 catch ( const fs::filesystem_error& e )
                                 {
                                     std::cerr << "Error: " << e.what() << std::endl;
                                 }
                             } );
        }
    }
    //gmsh_ptr->addGeoParameters( gmsh_ptr->retrieveGeoParameters( gmsh_ptr->description() ) );
    gmsh_ptr->addGeoParameters( geo_parameters );

    return gmsh_ptr;

}


}

#endif /* FEELPP_GEO_HPP */
