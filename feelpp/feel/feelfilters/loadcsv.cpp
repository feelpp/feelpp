/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 20 Jun 2016

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
#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadcsv.hpp>

#include <iostream>     
#include <fstream>      
#include <vector>
#include <string>
#include <algorithm>    
#include <iterator>     
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

namespace Feel {

std::map<std::vector<double>,double>
loadXYFromCSV( std::string const& filename,
               std::vector<std::string> const& abscissas,
               std::string const& ordinate )
{
    std::ifstream in(filename.c_str());
    if (!in.is_open()) return std::map<std::vector<double>,double>();

    typedef boost::tokenizer< boost::escaped_list_separator<char> > Tokenizer;

    std::vector< std::string > vec;
    std::string line;
    std::map<std::vector<double>,double> data;
    std::vector<double> d_abscissas;
    int i_abs = 0;
    int i_ord = 1;
    int count = 0;
    while (std::getline(in,line))
    {
        count ++;
        
        vec.clear();
        d_abscissas.clear();

        Tokenizer tok(line);
        vec.assign(tok.begin(),tok.end());

        if (vec.size() < 2) continue;

        
        boost::trim(vec[0]);
        boost::trim(vec[1]);
        if ( count == 1 )
        {
           
            auto it_ord = std::find( vec.begin(), vec.end(), ordinate );
            if ( it_ord == vec.end() )
                throw std::logic_error( "Invalid ordinate data lookiup in CSV file " + filename + " (" + ordinate + ")" );
            i_ord = std::distance( vec.begin(), it_ord );

            for( auto abscissa : abscissas )
            {
                auto it_abs = std::find( vec.begin(), vec.end(), abscissa );
                if ( it_abs == vec.end()  )
                    throw std::logic_error( "Invalid abscissa data lookiup in CSV file " + filename + " (" + abscissa + ")" );
                i_abs = std::distance( vec.begin(), it_abs );
                d_abscissas.push_back( std::stod( vec[i_abs] ) );
            }
            continue;
        }
        data.insert( std::make_pair( d_abscissas, std::stod(vec[i_ord]) ) );
    }
    LOG(INFO) << "done reading CSV file " << filename;
    return data;
}


}
