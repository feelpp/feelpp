//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 22 Jul 2019
//! @copyright 2019 Feel++ Consortium
//!
#ifndef FEELPP_NULLSPACE_RIGIDBODY_HPP
#define FEELPP_NULLSPACE_RIGIDBODY_HPP 1

#include <type_traits>

namespace Feel
{
template <typename SpaceType>
NullSpace<double> qsNullSpace( SpaceType const& space )
{
    using space_type = std::remove_cvref_t<decltype( *space )>;
    constexpr auto dim = space_type::nRealDim;
    static_assert( dim == 2 || dim == 3, "qsNullSpace only supports 2D and 3D vector spaces" );

    auto mode1 = space->element( oneX() );
    auto mode2 = space->element( oneY() );

    if constexpr ( dim == 2 )
    {
        auto mode3 = space->element( vec( Py(), -Px() ) );
        return NullSpace<double>( { mode1, mode2, mode3 } );
    }
    else
    {
        auto mode3 = space->element( oneZ() );
        auto mode4 = space->element( vec( Py(), -Px(), cst( 0. ) ) );
        auto mode5 = space->element( vec( -Pz(), cst( 0. ), Px() ) );
        auto mode6 = space->element( vec( cst( 0. ), Pz(), -Py() ) );
        return NullSpace<double>( { mode1, mode2, mode3, mode4, mode5, mode6 } );
    }
}

}

#endif
