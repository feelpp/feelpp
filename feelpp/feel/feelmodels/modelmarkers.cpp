/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 16 Mar 2015

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

#include <feel/feelmodels/modelmarkers.hpp>

namespace Feel {

ModelMarkers::ModelMarkers( std::string const& marker )
{
    this->insert(marker);
}

void
ModelMarkers::setPTree( pt::ptree const& p )
{
    M_p = p;
    this->clear();
    // markers = mark
    if( M_p.empty() )
        this->insert(M_p.template get_value<std::string>());
    else
    {
        for( auto const& item : M_p )
        {
            // markers = [mark1,mark2]
            if( item.first.empty() ) // array
                this->insert(item.second.template get_value<std::string>());
            else if( item.first == "name" )
            {
                // markers : { name = mark }
                if( item.second.empty() )
                    this->insert(item.second.template get_value<std::string>() );
                else
                {
                    // markers : { name = [mark1, mark2] }
                    for( auto const& item2 : item.second )
                        this->insert( item2.second.template get_value<std::string>() );
                }
            }
        }
        this->generateMarkersList();
    }
}

void ModelMarkers::generateMarkersList()
{
    for( int k = 0; k < 10; ++k )
    {
        std::vector<std::string> argK = this->getArguments(k);
        this->generateMarkersListForIndex(k, argK);
    }
}

std::vector<std::string>
ModelMarkers::getArguments( int k )
{
    std::string strK = (boost::format("index%1%")%k).str();
    std::vector<std::string> argK;
    if( auto indiceK = M_p.get_child_optional(strK) )
    {
        // indexN = [ A, B ]
        for( auto const& item : M_p.get_child(strK) )
            argK.push_back(item.second.template get_value<std::string>());
        // index N = start:end(:step)
        if( argK.empty() )
        {
            std::string rangeK = M_p.template get<std::string>(strK);
            boost::char_separator<char> sep(":");
            boost::tokenizer<boost::char_separator<char> > kvlist( rangeK, sep);
            int sizeRange = std::distance(kvlist.begin(),kvlist.end());
            if( sizeRange != 2 && sizeRange != 3 )
            {
                LOG(WARNING) << "range " << rangeK << " for markers has "
                             << sizeRange << " elements, should be 2 or 3";
                return argK;
            }

            std::vector<int> range;
            for( auto const& r : kvlist )
                range.push_back(std::stoi(r));
            if( sizeRange == 2 )
                range.push_back(1);
            auto argKInt = boost::irange(range[0],range[1],range[2]);
            for( auto const& i : argKInt )
                argK.push_back( std::to_string(i));
        }
    }
    return argK;
}

void ModelMarkers::generateMarkersListForIndex( int k, std::vector<std::string> argK )
{
    auto it = this->begin();
    while( it != this->end() )
    {
        std::string m = *it;
        auto pos = m.find("%"+std::to_string(k)+"%");
        if( pos != std::string::npos )
        {
            std::set<std::string> ms;
            for( auto const& s : argK )
            {
                std::string currentMarker = m;
                currentMarker.replace(pos, 3, s);
                ms.insert(currentMarker);
            }
            this->erase(it);
            this->insert(ms.begin(),ms.end());
            it = this->begin();
        }
        else
            std::advance(it,1);
    }
}

}
