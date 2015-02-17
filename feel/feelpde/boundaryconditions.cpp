/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 25 Jan 2015
 
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
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelpde/boundaryconditions.hpp>

namespace Feel
{

BoundaryConditions::BoundaryConditions()
    :
    BoundaryConditions( "" )
{}

BoundaryConditions::BoundaryConditions( std::string const& p )
    :
    super(),
    M_prefix( p )
{
    fs::path bc( Environment::expand( soption("bc-file") ) );
    
    if ( fs::exists( bc ) )
    {
        LOG(INFO) << "Loading Boundary Condition file " << bc.string();
        load( bc.string() );
    }
    else
    {
        LOG(WARNING) << "Boundary condition file " << bc.string() << " does not exist";
    }
}

void
BoundaryConditions::load(const std::string &filename)
{
    // Create an empty property tree object
    using boost::property_tree::ptree;
    ptree pt;

    read_json(filename, pt);

    for( auto const& v : pt )
    {
        
        //std::cout << "v.first:" << v.first  << "\n";
        std::string t = v.first; // field name
        for( auto const& f : v.second )
        {
            std::string k = t+"."+f.first; // condition type
            for( auto const& c : f.second ) // condition
            {
                try
                {
                    auto e= c.second.get<std::string>("expr");
                    LOG(INFO) << "adding boundary " << c.first << " with expression " << e << " to " << k;
                    this->operator[](t)[f.first].push_back( std::make_pair( c.first, e ) );
                }
                catch( ... )
                {
                    LOG(INFO) << "adding boundary " << c.first << " without expression" << " to " << k;
                    this->operator[]( t )[f.first].push_back( std::make_pair( c.first, std::string("") ) );
                }
            }
        }
        
    }
    if ( Environment::isMasterRank() )
    {
        for( auto const& s : *this )
        {
            LOG(INFO) << "field " << s.first << "\n";
            for( auto const& t : s.second )
            {
                LOG(INFO) << " - type " << t.first << "\n";
                for( auto const& c : t.second )
                {
                    LOG(INFO) << "  . boundary  " << c.first << " expr : " << c.second << "\n";
                }
                
            }
        }
    }
}
    


}//Feel
