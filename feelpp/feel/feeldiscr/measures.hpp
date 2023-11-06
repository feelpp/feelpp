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
template<typename RangeType>
using value_t = typename decay_type<RangeType>::value_t;

/**
 * compute h measure statistics over a range of elements or facets
 */
template<typename RangeType, typename = std::enable_if_t<is_range_v<RangeType>>>
std::tuple<value_t<RangeType>, value_t<RangeType>, value_t<RangeType>>
hMeasures( RangeType && r, size_t nEntityInRange )
{
    using value_type=value_t<RangeType>;
    value_type h_avg = 0;
    value_type h_min = std::numeric_limits<value_type>::max();
    value_type h_max = 0;
    std::for_each( std::forward<RangeType>(r).begin(), std::forward<RangeType>(r).end(), 
        [&]( element_t<RangeType> const& elt ) {
                h_avg += elt.h();
                h_min = std::min( h_min, elt.h() );
                h_max = std::max( h_max, elt.h() );
        } );

    value_type reduction[3] = {h_avg, h_min, h_max};
    MPI_Op op;
    MPI_Op_create( (MPI_User_function*)(Functor::AvgMinMax<value_type, WorldComm::communicator_type>), 1, &op );
    MPI_Allreduce( MPI_IN_PLACE, reduction, 3, mpi::get_mpi_datatype<value_type>(), op, std::forward<RangeType>(r).worldComm() );
    MPI_Op_free( &op );

    h_avg = reduction[0]/nEntityInRange;
    h_min = reduction[1];
    h_max = reduction[2];
    return std::tuple{ h_avg, h_min, h_max };
}

template<typename RangeType, typename = std::enable_if_t<is_mesh_v<RangeType>>>
std::tuple<value_t<RangeType>, value_t<RangeType>, value_t<RangeType>>
hMeasures( RangeType && r )
{
    return hMeasures( std::forward<RangeType>(r), nelements( std::forward<RangeType>(r), Zone::GLOBAL ) );
}

/**
 * compute h measures statistics over all elements of a mesh
 */
template<typename MeshType, typename = std::enable_if_t<is_mesh_v<MeshType>>>
auto
hMeasures( MeshType&& m )
{
    return hMeasures( elements( std::forward<MeshType>(m) ), std::forward<MeshType>(m).numGlobalElements() );
}

/**
 * @brief compute h measures statistics over all elements of a shared_ptr<mesh>
 * 
 * @tparam MeshType type of the mesh
 * @param m shared_ptr<mesh>
 * @return auto 
 */
template<typename MeshType, typename = std::enable_if_t<is_mesh_v<MeshType>>>
auto
hMeasures( std::shared_ptr<MeshType> const& m )
{
    return hMeasures( elements( m ), m->numGlobalElements() );
}
}
