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

ModelModel::ModelModel( std::string const& name, nl::json const& jarg, bool addVariants )
    :
    M_name( name ),
    M_ptree( jarg )
{
    if ( jarg.contains( "equations" ) )
    {
        auto const& j_eq = jarg.at( "equations" );
        if ( j_eq.is_string() )
            M_equations = j_eq.get<std::string>();
    }

    if ( jarg.contains( "submodels" ) )
    {
        auto const& j_submodels = jarg.at( "submodels" );
        if ( j_submodels.is_string() )
            M_submodels.insert( j_submodels.get<std::string>() );
        else if ( j_submodels.is_array() )
        {
            for ( auto const& [j_submodelskey,j_submodelsval]: j_submodels.items() )
                if ( j_submodelsval.is_string() )
                    M_submodels.insert( j_submodelsval.get<std::string>() );
        }
    }

#if 0 // TODO VINCENT
    if ( addVariants )
    {
        for( auto const& [key,ptreeVariant] : M_ptree )
        {
            if ( key == "equations" || key == "submodels" )
                continue;
            if ( ptreeVariant.empty() )
                continue;

            bool ptreeVariantIsArray = false;
            for ( auto const& item : ptreeVariant )
            {
                if ( item.first.empty() )
                {
                    ptreeVariantIsArray = true;
                    break;
                }
            }
            if ( ptreeVariantIsArray )
                continue;

            M_variants.emplace( key, ModelModel( key, ptreeVariant, false ) );
        }
    }
#endif
}

ModelModels::ModelModels()
    :
    M_useModelName( false )
{}

void
ModelModels::setup()
{
    auto const& jarg = M_p;
    if ( jarg.contains("use-model-name") )
    {
        auto const& j_useModelName = jarg.at("use-model-name");
        if ( j_useModelName.is_boolean() )
            M_useModelName = j_useModelName.template get<bool>();
        else if ( j_useModelName.is_string() )
            M_useModelName = boost::lexical_cast<bool>( j_useModelName.template get<std::string>() );
    }

    if ( M_useModelName )
    {
        for (auto const& [jargkey,jargval] : jarg.items())
        {
            if ( !jargval.is_object() )
                continue;
            std::string const& name = jargkey;
            this->insert( std::make_pair(name,ModelModel( name,jargval ) ) );
        }
    }
    else
    {
        std::string name = "";
        this->insert( std::make_pair( name, ModelModel(name, jarg ) ) );
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
