/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-07-21

  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-07-21
 */

#include <life/lifefilters/exporter.hpp>

namespace Life
{
//
// Exporter Options
//
/**
 * \return the command lines options for the exporter
 */
po::options_description exporter_options( std::string const& prefix )
{
    std::string _prefix = prefix;
    if ( !_prefix.empty() )
        _prefix += "-";

    po::options_description _options( "Exporter " + prefix + " options");
    _options.add_options()
        // do export
        ((_prefix+"export").c_str(), Life::po::value<bool>()->default_value( true ), "true if export, false otherwise")

        // exporter type
        ((_prefix+"exporter-format").c_str(), Life::po::value<std::string>()->default_value( "ensight" ), "type of exporter")

        // prefix options
        ((_prefix+"exporter-prefix").c_str(), Life::po::value<std::string>()->default_value( prefix ), "prefix for exported files")

        // frequency options
        ((_prefix+"exporter-freq").c_str(), Life::po::value<int>()->default_value( 1 ), "frequency at which results are exported")

        // file type options
        ((_prefix+"exporter-file-type").c_str(), Life::po::value<int>()->default_value( ASCII ), "file type in which the results are exported ('ascii' = 0 or 'binary' = 1)")
        ;
    std::cout << "exporter options : " << _options << "\n";
    return _options;
}

}
