/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@cemosis.fr>
 Date: 30 Sept 2024

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

#pragma once

#include <feel/feelmesh/filters.hpp>

namespace Feel {

//! This function takes a mesh and a range of entities defined on another mesh.
//! If a mesh relation exists between mesh and range (sub mesh or parent mesh),
//! we return a new range with the same entities but defined on this mesh.
//! If some entities can not be retrieved, we ignore them.
//! NOTE: For now, we support only faces entity and only relation with same kind of entity (i.e mesh dim == range mesh dim)
template<typename MeshType,typename RangeType>
FEELPP_EXPORT Range<typename MeshTraits<MeshType>::mesh_type,RangeType::mesh_entities>
migrate( MeshType const& mesh, RangeType const& range )
{
    using mesh_type = typename MeshTraits<MeshType>::mesh_type;
    using index_type = typename mesh_type::index_type;
    using range_mesh_type = typename RangeType::mesh_t;
    auto meshRange = range.mesh();
    Range<mesh_type,RangeType::mesh_entities> outputRange( mesh );

    int relationCase = -1;
    if ( meshRange->isSubMeshFrom( mesh ) )
        relationCase = 0;
    else if ( meshRange->isParentMeshOf( mesh ) )
        relationCase = 1;
    auto getRelatedEntityId = [&relationCase,&meshRange,&mesh]( auto const& entity ){
                                  switch ( relationCase )
                                  {
                                  case 0: return meshRange->subMeshToMesh( entity.id() );
                                  case 1: return mesh->meshToSubMesh( entity.id() );
                                  default: return entity.id();
                                  };
                              };

    if constexpr ( RangeType::mesh_entities == MESH_FACES )
    {
        index_type outputEltId = invalid_v<index_type>;
        for( auto const& faceWrap : range )
        {
            auto const& face = unwrap_ref( faceWrap );
            if ( !face.isConnectedTo0() )
                continue;
            if ( face.isGhostFace() ) continue;////

            auto const& elt0 = face.element0();
            outputEltId = getRelatedEntityId( elt0 );
            if ( outputEltId != invalid_v<index_type> )
            {
                auto const& outputElt = mesh->element( outputEltId );
                //if ( !outputElt.face( face.idInElement0() ).isGhostFace() )
                    outputRange.push_back( outputElt.face( face.idInElement0() ) );
            }
            else if ( face.isConnectedTo1() )
            {
                auto const& elt1 = face.element1();
                outputEltId = getRelatedEntityId( elt1 );
                if ( outputEltId != invalid_v<index_type> )
                {
                    auto const& outputElt = mesh->element( outputEltId );
                    //if ( !outputElt.face( face.idInElement1() ).isGhostFace() )
                        outputRange.push_back( outputElt.face( face.idInElement1() ) );
                }
            }
        }
    }
    else
        CHECK( false ) << "TODO";

    outputRange.shrink_to_fit();
    return outputRange;
}

}
