/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-07-21

  Copyright (C) 2007,2008,2009,2010 Université Joseph Fourier (Grenoble I)

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
   \file exporter.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-07-21
 */

#include <feel/feelfilters/exporter.hpp>

namespace Feel
{
//
// Exporter Options
//
/**
 * \return the command lines options for the exporter
 */
po::options_description exporter_options( std::string const& prefix )
{
    po::options_description _options( "Exporter " + prefix + " options" );
    _options.add_options()
        // do export
        ( prefixvm( prefix,"export" ).c_str(), Feel::po::value<bool>()->default_value( true ), "true if export, false otherwise" )
        // do export
        ( prefixvm( prefix,"exporter.export" ).c_str(), Feel::po::value<bool>()->default_value( true ), "true if export, false otherwise" )

        // exporter type
        ( prefixvm( prefix,"exporter.format" ).c_str(), Feel::po::value<std::string>()->default_value( "ensight" ), "type of exporter (ensight, gmsh, ensightgold or exodus)" )



        //  geometry
        ( prefixvm( prefix,"exporter.geometry" ).c_str(), Feel::po::value<int>()->default_value( (int)EXPORTER_GEOMETRY_CHANGE_COORDS_ONLY ), "type of geometry" )

        // prefix options
        ( prefixvm( prefix,"exporter.prefix" ).c_str(), Feel::po::value<std::string>()->default_value( prefix ), "prefix for exported files" )

        // directory options
        ( prefixvm( prefix,"exporter.directory" ).c_str(), Feel::po::value<std::string>()->default_value( "results" ), "directory for exported files" )

        // frequency options
        ( prefixvm( prefix,"exporter.freq" ).c_str(), Feel::po::value<int>()->default_value( 1 ), "frequency at which results are exported" )

        // file type options
        ( prefixvm( prefix,"exporter.file-type" ).c_str(), Feel::po::value<int>()->default_value( ASCII ), "file type in which the results are exported ('ascii' = 0 or 'binary' = 1)" )

        // matlab options
        ( prefixvm( prefix,"exporter.matlab" ).c_str(), Feel::po::value<bool>()->default_value( 0 ), "export matrices and vectors to matlab files" )

        //
        // ensightgold
        //
        ( prefixvm( prefix,"exporter.ensightgold.use-sos" ).c_str(), Feel::po::value<bool>()->default_value( true ), "use sos  (true) or first case file (false) for multiple case files" )
        //  single
        ( prefixvm( prefix,"exporter.fileset" ).c_str(), Feel::po::value<bool>()->default_value( false ), "use fileset for transient simulations" )

        ;
    return _options;
}

}
