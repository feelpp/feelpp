/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 09 Mar 2020

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
template<typename MeshType>
using value_t = typename MeshType::value_type;
template<typename MeshType>
using index_t = typename MeshType::index_type;

/**
 * compute h measure statistics over a range of elements or facets
 */
template<typename MeshType, typename RangeType, typename = std::enable_if_t<is_mesh_v<MeshType>>>
std::tuple<value_t<MeshType>, value_t<MeshType>, value_t<MeshType>>
hMeasures( std::shared_ptr<MeshType> const& m, RangeType r, typename MeshType::size_type nEntityInRange )
{
    using value_type=value_t<MeshType>;
    value_type h_avg = 0;
    value_type h_min = std::numeric_limits<value_type>::max();
    value_type h_max = 0;
    for ( auto const& eltWrap : r )
    {
        auto const& elt = unwrap_ref( eltWrap );
        h_avg += elt.h();
        h_min = std::min( h_min, elt.h() );
        h_max = std::max( h_max, elt.h() );
    }
    h_avg /= nEntityInRange;
    value_type reduction[3] = {h_avg, h_min, h_max};
    MPI_Op op;
    MPI_Op_create( (MPI_User_function*)(Functor::AvgMinMax<value_type, WorldComm::communicator_type>), 1, &op );
    MPI_Allreduce( MPI_IN_PLACE, reduction, 3, mpi::get_mpi_datatype<value_type>(), op, m->worldComm() );
    MPI_Op_free( &op );

    h_avg = reduction[0];
    h_min = reduction[1];
    h_max = reduction[2];
    return std::tuple{ h_avg, h_min, h_max };
}

template<typename MeshType, typename RangeType, typename = std::enable_if_t<is_mesh_v<MeshType>>>
std::tuple<value_t<MeshType>, value_t<MeshType>, value_t<MeshType>>
hMeasures( std::shared_ptr<MeshType> const& m, RangeType r )
{
    return hMeasures( m, elements( m ), nelements( r, Zone::GLOBAL ) );
}

/**
 * compute h measures statistics over all elements of a mesh
 */
template<typename MeshType, typename = std::enable_if_t<is_mesh_v<MeshType>>>
std::tuple<value_t<MeshType>, value_t<MeshType>, value_t<MeshType>>
hMeasures( std::shared_ptr<MeshType> const& m )
{
    return hMeasures( m, elements( m ), m->numGlobalElements() );
}
}
