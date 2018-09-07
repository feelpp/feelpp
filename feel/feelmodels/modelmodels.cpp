/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 9 Feb 2018

 Copyright (C) 2018 Feel++ Consortium

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

#include <feel/feelmodels/modelmodels.hpp>

namespace Feel {

ModelModel::ModelModel( pt::ptree const& p )
    :
    M_ptree( p )
{
    if ( auto eq = M_ptree.get_optional<std::string>("equations") )
        M_equations = *eq;
}

ModelModels::ModelModels()
    :
    M_useModelName( false )
{}

void
ModelModels::setup()
{
    if ( auto useModelName = M_p.get_optional<bool>("use-model-name") )
        M_useModelName = *useModelName;

    if ( M_useModelName )
    {
        for( auto const& p1 : M_p )
        {
            this->insert( std::make_pair( p1.first,ModelModel( p1.second ) ) );
        }
    }
    else
    {
        std::string name = "";
        this->insert( std::make_pair( name, ModelModel( M_p ) ) );
    }
}

bool
ModelModels::hasModel( std::string const& name ) const
{
    std::string nameUsed = (M_useModelName)? name : "";
    return this->find( nameUsed ) != this->end();
}

ModelModel const&
ModelModels::model( std::string const& name ) const
{
    std::string nameUsed = (M_useModelName)? name : "";
    auto itFindModel = this->find( nameUsed );
    if ( itFindModel != this->end() )
        return itFindModel->second;
    return M_emptyModel;
}

}
