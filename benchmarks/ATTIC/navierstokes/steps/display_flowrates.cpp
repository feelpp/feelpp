/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 03 juil. 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
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
#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <cmath>


bool
display_flowrates_header( std::ostream& os, std::map<std::string, double> const& flowrates )
{
    os << "| Time | ";
    for( auto & f : flowrates )
        os << "| Reynolds | " << std::setw(11) << std::right << f.first  << " | " << std::setw(11) << std::right << " Error | " ;
    os << "|\n";
    return true;
}

bool
display_flowrates( std::ostream& os, double t, std::map<std::string, double> const& flowrates, double Q )
{
    os << "| " << t << " | ";
    for( auto & f : flowrates )
        os << std::setw(11) << std::scientific << std::setprecision( 4 )<<std::abs(f.second)*0.012*1056/(0.0035*3.14*0.006*0.006)<< " | " << std::setw(11) << std::scientific << std::setprecision( 4 ) << std::abs(f.second)<< " | " << std::setw(11) << std::scientific << std::setprecision( 4 )<<((std::abs(f.second)-Q)/Q)*100 << "  | " ;
    os << "|\n";
    return true;
}
