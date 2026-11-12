/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileContributor: Vincent Chabannes <vincent.chabannes@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2017-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
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
