/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 22 Oct 2020

 Copyright (C) 2020 Feel++ Consortium

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
#pragma once

namespace Feel
{

template <typename MeshType>
class DataByMeshEntity
{
public :
    using mesh_type = MeshType;
    using index_type = typename mesh_type::index_type;
    using value_type = double;
    using mapping_id2value_type = std::map<index_type,value_type>;

    DataByMeshEntity( std::shared_ptr<mesh_type> const& mesh, ElementsType entityType, mapping_id2value_type const& entityIdToValue )
        :
        M_mesh( mesh ),
        M_entityType( entityType ),
        M_entityIdToValue( entityIdToValue )
        {}

    DataByMeshEntity( DataByMeshEntity const& ) = default;
    DataByMeshEntity( DataByMeshEntity && ) = default;

    std::shared_ptr<mesh_type> const& mesh() const { return M_mesh; }
    ElementsType entityType() const { return M_entityType; }
    mapping_id2value_type const& entityIdToValue() const { return M_entityIdToValue; }

    std::optional<value_type> valueAtEntityIdIfExists( index_type id ) const
        {
            auto itFindVal = M_entityIdToValue.find( id );
            if ( itFindVal != M_entityIdToValue.end() )
                return std::make_optional( itFindVal->second );
            else
                return std::optional<value_type>{};
        }
private :
    std::shared_ptr<mesh_type> M_mesh;
    ElementsType M_entityType;
    mapping_id2value_type M_entityIdToValue;
};

template <typename MeshType>
class CollectionOfDataByMeshEntity : public std::map<std::string,DataByMeshEntity<MeshType>>
{
public :
    using mesh_type = MeshType;
    using index_type = typename mesh_type::index_type;
    using data_by_mesh_entity_type = DataByMeshEntity<MeshType>;

    CollectionOfDataByMeshEntity( std::shared_ptr<mesh_type> const& mesh, std::string const& filename )
        :
        M_mesh( mesh )
        {}

    data_by_mesh_entity_type const& get( std::string const& field ) const { return this->find( field )->second; }

    void add( std::string const& field, ElementsType entityType, typename data_by_mesh_entity_type::mapping_id2value_type && data )
    {
        this->emplace( std::make_pair( field, data_by_mesh_entity_type( M_mesh, entityType, std::forward<typename data_by_mesh_entity_type::mapping_id2value_type>( data ) ) ) );
    }
private :
    std::shared_ptr<mesh_type> M_mesh;
};

} // namespace Feel
