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
#include <fmt/core.h>
#include <fmt/ranges.h>     
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

namespace Feel {

std::vector<Eigen::VectorXd>
loadXYFromCSV( std::string const& filename,
               std::vector<std::string> const& abscissas )
{
    std::ifstream in(filename.c_str());
    if (!in.is_open()) return std::vector<Eigen::VectorXd>();

    using Tokenizer = boost::tokenizer< boost::char_separator<char> >;

    std::vector< std::string > vec;
    std::string line;
    std::vector<Eigen::VectorXd> data;
    std::vector<int> i_abscissas;
    i_abscissas.reserve( abscissas.size() );
    std::vector<double> values( abscissas.size() );
    values.reserve( abscissas.size() );
    std::cout << fmt::format( "load {} from {} ", abscissas, filename ) << std::endl;
    boost::char_separator<char> els( " ;," );
    int count = 0;
    while (std::getline(in,line))
    {
        count ++;
        
        vec.clear();
        values.clear();

        Tokenizer tok(line, els);
        vec.assign(tok.begin(),tok.end());
        std::cout << fmt::format( "count : {} vec {} size {}", count, fmt::join( vec, " " ), vec.size() ) << std::endl;
        if (vec.size() < 2) continue;

        
        std::for_each(vec.begin(), vec.end(), [](std::string & s) { boost::trim(s); });

        if ( count == 1 )
        {
            for( auto abscissa : abscissas )
            {
                std::cout << fmt::format("Check abscissa _{}_ ", abscissa ) << std::endl;
                auto it_abs = std::find( vec.begin(), vec.end(), abscissa );
                if ( it_abs == vec.end()  )
                    throw std::logic_error( "Invalid abscissa data lookup in CSV file " + filename + " (" + abscissa + ")" );
                int i_abs = std::distance( vec.begin(), it_abs );
                std::cout << fmt::format( "Check abscissa _{}_: {} ", abscissa, i_abs ) << std::endl;
                i_abscissas.push_back( i_abs );
                
            }
            LOG(INFO) << fmt::format( "[loadcsv] indices {}",  i_abscissas ) << std::endl;
            continue;
        }
        std::cout << fmt::format( "count : {} i_abscissas {} ", count, fmt::join( i_abscissas, " " ) ) << std::endl;
        for( auto index : i_abscissas )
        {
            values.push_back( std::stod( vec[index] ) );
        }
        data.push_back( Eigen::Map<Eigen::VectorXd>( values.data(), values.size() ) );
    }
    LOG(INFO) << fmt::format( "[loadcsv] csv data loaded: {}",  data ) << std::endl;
    LOG(INFO) << "done reading CSV file " << filename;
    return data;
}


}
