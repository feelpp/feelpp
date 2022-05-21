/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 Date: 05 Oct 2020

 Copyright (C) 2020 Université de Strasbourg

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
#include <feel/feeldiscr/traits.hpp>

namespace Feel
{

template <typename MeshT, typename = std::enable_if_t<is_mesh_v<MeshT>>>
void dump( std::shared_ptr<MeshT> const& mesh, std::string const& header = {} )
{
    //if ( mesh->worldComm().localRank() == 0 )
    {
        using namespace Feel;
        Feel::cout << fmt::format( "===================[ {} ]===================", header ) << std::endl;
        Feel::cout << "Information about the mesh:";
        Feel::cout << "      number of elements in memory : " << mesh->numGlobalElements() << std::endl;
        Feel::cout << "      number of faces in memory : " << mesh->numGlobalFaces() << std::endl;
        if ( dimension_v<MeshT> == 3 )
            Feel::cout << "      number of edges in memory : " << mesh->numGlobalEdges() << std::endl;
        Feel::cout << "      number of points  in memory : " << mesh->numGlobalPoints() << std::endl;
        for ( auto marker : mesh->markerNames() )
        {
            auto [name, data] = marker;

            if ( data[1] == mesh->dimension() )
            {
                size_type nelts = nelements( markedelements( mesh, name ), true );
                Feel::cout << "      number of marked elements " << name << " with tag " << data[0] << " : " << nelts << std::endl;
            }
        }
        for ( auto marker : mesh->markerNames() )
        {
            auto name = marker.first;
            auto data = marker.second;

            if ( data[1] == mesh->dimension() - 1 )
            {
                size_type nelts = nelements( markedfaces( mesh, name ), true );
                Feel::cout << "      number of marked faces " << name << " with tag " << data[0] << " : " << nelts << std::endl;
            }
        }
        Feel::cout << fmt::format( "===================[ {} ]===================", header ) << std::endl;
    }
}

} // namespace Feel
