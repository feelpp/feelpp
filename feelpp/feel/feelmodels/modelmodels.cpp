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

ModelModel::ModelModel( std::string const& type, std::string const& name, nl::json const& jarg )
    :
    M_type( type ),
    M_name( name )
{
    if ( jarg.contains( "name" ) )
    {
        auto const& j_name = jarg.at( "name" );
        if ( j_name.is_string() )
            M_name = j_name.get<std::string>();
    }

    for ( std::string const& submodelName : { "submodels", "subphysics" } )
        if ( jarg.contains( submodelName ) )
        {
            auto const& j_submodels = jarg.at( submodelName );
            if ( j_submodels.is_object() )
            {
                for ( auto const& [j_submodelskey,j_submodelsval]: j_submodels.items() )
                {
                    if ( j_submodelsval.is_string() )
                        M_submodels[j_submodelskey].insert( j_submodelsval.get<std::string>() );
                    else if ( j_submodelsval.is_array() )
                        for ( auto const& [j_submodels_namekey,j_submodels_nameval] : j_submodelsval.items() )
                            if ( j_submodels_nameval.is_string() )
                                M_submodels[j_submodelskey].insert( j_submodels_nameval.get<std::string>() );
                }
            }
            else if ( j_submodels.is_array() )
            {
            }
    }

    for ( std::string const& setupName : { "equations", "equation", "setup" } )
        if ( jarg.contains( setupName ) )
            M_setup = jarg.at( setupName );


    if ( jarg.contains( "materials" ) )
    {
        auto const& j_materials = jarg.at( "materials" );
        if ( j_materials.is_string() )
            M_materials.insert( j_materials.get<std::string>() );
        else if ( j_materials.is_array() )
        {
            for ( auto const& [j_materialskey,j_materialsval]: j_materials.items() )
                if ( j_materialsval.is_string() )
                    M_materials.insert( j_materialsval.get<std::string>() );
        }
    }
}



void
ModelModelsSameType::setup( std::string type, nl::json const& jarg )
{
    ModelModel m( type, jarg );
    this->insert( std::make_pair( m.name(), std::move( m ) ) );
}

void
ModelModels::setup( nl::json const& jarg )
{
    for (auto const& [jargkey,jargval] : jarg.items())
    {
        std::string const& type = jargkey;
        if ( jargval.is_array() )
        {
            for (auto const& [j_modeltypekey,j_modeltypeval] : jargval.items())
            {
                this->operator[]( type ).setup( type, j_modeltypeval );
            }
        }
        else if ( jargval.is_object() )
        {
            if ( jargval.contains("common") )
            {
                auto const& j_common = jargval.at("common");
                CHECK( jargval.contains("models") ) << "required a models section";
                auto const& j_models = jargval.at("models");
                if ( j_models.is_array() )
                {
                    for (auto const& [j_modelskey,j_modelsval] : j_models.items())
                    {
                        nl::json j_data = j_common;
                        j_data.merge_patch( j_modelsval );
                        this->operator[]( type ).setup( type, j_data );
                    }
                }
                else if ( j_models.is_object() )
                {
                    //nl::json j_data = { { "setup", j_common } };
                    nl::json j_data = j_common;
                    j_data.merge_patch( j_models );
                    this->operator[]( type ).setup( type, j_data );
                }
            }
            else
                this->operator[]( type ).setup( type, jargval );
        }
    }
}

#if 0
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
#endif
}
