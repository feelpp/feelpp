/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 19 Jun 2016

 Copyright (C) 2016 Feel++ Consortium

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
#ifndef FEELPP_LOADCSV_HPP
#define FEELPP_LOADCSV_HPP 1

#include <map>
#include <string>

namespace Feel
{

//! load a XY dataset from a CSV file
//! \param filename name of the CSV file
//! \param abscissa name of the column to be used as X
//! \param ordinate name of the column to be used as Y
//! \return the map<key,value> {X,Y}
//! here is an example
//! \code
//!  // load the columns "time" and "ft" into a map<>
//!  auto m = loadXYFromCSV( "data.csv", "time", "ft" );
//! \endcode
std::map<std::vector<double>, double>
loadXYFromCSV( std::string const& filename,
               std::vector<std::string> const& abscissas,
               std::string const& ordinate );

inline std::map<double,double>
loadXYFromCSV( std::string const& filename,
               std::string const& abscissa,
               std::string const& ordinate )
{
    std::map<std::vector<double>, double> rs = loadXYFromCSV( filename, std::vector{ abscissa }, ordinate );
    std::map<double,double> m;
    for( auto const& [key,value] : rs )
    {
        m[key[0]] = value;
    }
    return m;
}


} // Feel
#endif
