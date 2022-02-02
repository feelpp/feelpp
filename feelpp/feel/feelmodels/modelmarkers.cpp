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
ModelMarkers::setPTree( pt::ptree const& p, ModelIndexes const& indexes )
{
    this->clear();
    // markers = mark
    if( p.empty() )
        this->insert( indexes.replace( p.template get_value<std::string>()) );
    else
    {
        auto indexesAllCases = ModelIndexes::generateAllCases( p, indexes );
        for ( auto const& indexesNew : indexesAllCases )
        {
            for( auto const& item : p )
            {
                // markers = [mark1,mark2]
                if( item.first.empty() ) // array
                    this->insert( indexesNew.replace( item.second.template get_value<std::string>() ) );
                else if( item.first == "name" )
                {
                    // markers : { name = mark }
                    if( item.second.empty() )
                        this->insert(indexesNew.replace( item.second.template get_value<std::string>() ) );
                    else
                    {
                        // markers : { name = [mark1, mark2] }
                        for( auto const& item2 : item.second )
                            this->insert( indexesNew.replace( item2.second.template get_value<std::string>() ) );
                    }
                }
            }
        }
    }
}

void
ModelMarkers::setup( nl::json const& jarg, ModelIndexes const& indexes )
{
    this->clear();

    if ( jarg.is_string() ) // markers = mark
        this->insert( indexes.replace( jarg.get<std::string>() ) );
    else
    {
        auto indexesAllCases = ModelIndexes::generateAllCases( jarg, indexes );
        for ( auto const& indexesNew : indexesAllCases )
        {
            if ( jarg.is_array() ) // markers = [mark1,mark2]
            {
                for ( auto const& el : jarg.items() )
                    this->insert( indexesNew.replace( el.value() ) );
            }
            else if ( jarg.is_object() )
            {
                if ( jarg.contains("name") )
                {
                    auto const& jName = jarg.at("name");
                     if ( jName.is_string() ) // markers : { name = mark }
                         this->insert( indexesNew.replace( jName.get<std::string>() ) );
                     else if ( jName.is_array() ) // markers : { name = [mark1, mark2] }
                     {
                         for ( auto const& el : jName.items() )
                             this->insert( indexesNew.replace( el.value() ) );
                     }
                }
            }
        }
    }
}

bool
ModelMarkers::empty() const
{
    if( this->size() > 1 )
        return false;
    if( this->size() == 0 )
        return true;
    return this->begin()->empty();
}

}
