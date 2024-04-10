/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

 This file is part of the Feel library

 Author(s): Vincent Chabannes <vincent.chabannes@cemosis.fr>
 Date: 2024-03-27

 Copyright (C) 2024 Feel++ Consortium

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
#ifndef FEELPP_MESH_MESHFRAGMENTATION_HPP
#define FEELPP_MESH_MESHFRAGMENTATION_HPP

#include <feel/feelmesh/filters.hpp>

namespace Feel
{

template<typename MeshType>
class MeshFragmentation
{
public:
    using mesh_type = MeshType;
    using range_element_type = Range<mesh_type,MESH_ELEMENTS>;
    using container_type = std::map<int,std::tuple<std::string,range_element_type>>; // partId -> ( partName, range )

    enum class Strategy{ None=0, AllMarkedElements, Custom };

    explicit MeshFragmentation( Strategy strategy = Strategy::AllMarkedElements ) : M_strategy( strategy ) {}
    MeshFragmentation( MeshFragmentation const& ) = default;
    MeshFragmentation( MeshFragmentation && ) = default;
    MeshFragmentation& operator=( MeshFragmentation const& ) = default;
    MeshFragmentation& operator=( MeshFragmentation && ) = default;

    //! build container of parts from a mesh
    container_type toContainer( mesh_type const& mesh ) const
        {
            container_type res;
            switch ( M_strategy )
            {
            case Strategy::None:
                res[0] = std::make_tuple("elements",elements(mesh));
                break;
            case Strategy::AllMarkedElements:
                for ( auto const& [fragmentId,fragmentData] : fragmentationMarkedElements( mesh ) )
                {
                    auto const& [range,mIds,fragmentName] = fragmentData;
                    res[fragmentId] = std::make_tuple(fragmentName,range);
                }
                break;
            case Strategy::Custom:
                return M_customParts;
            }
            return res;
        }

    //! add a part from fragmentId, fragmentName and range of element
    template <typename RT>
    auto addFragment( int fragmentId, std::string const& fragmentName, RT && r )
        {
            CHECK( M_strategy == Strategy::Custom ) << "add fragment only valid with custom strategy";
            return M_customParts.emplace( std::make_tuple( fragmentId, fragmentName, std::forward<RT>( r ) ) );
        }

private:
    Strategy M_strategy;
    container_type M_customParts;// M_partIdToRangeElement;
};

}

#endif
