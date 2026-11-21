/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#pragma once

#include <feel/feelcore/environment.hpp>
#include <feel/feelpoly/order.hpp>
#include <charconv>
#include <stdexcept>
#include <string>

namespace Feel::quickstart
{
inline RuntimeOrder checkedOrder( std::string const& option, int minimum = 1 )
{
    int const value = ioption( option );
    if ( value < minimum )
        throw std::invalid_argument( option + " must be >= " + std::to_string( minimum ) );
    return RuntimeOrder::checked( value );
}

inline RuntimeOrder geometryOrder()
{
    return checkedOrder( "geometry.order" );
}

template<typename MeshPtr>
void checkMeshOrder( MeshPtr const& mesh )
{
    auto const requested = geometryOrder();
    if ( !mesh || mesh->order() != requested.value )
        throw std::invalid_argument( "loaded mesh order differs from geometry.order" );
}

// Preserve existing case.discretization=Pk configurations while permitting
// an explicit integer --order for the same runtime Lagrange space.
inline RuntimeOrder lagrangeOrder()
{
    auto const& vm = Environment::vm();
    auto const discretization = soption( "case.discretization" );
    int value = 0;
    if ( discretization.size() < 2 || discretization.front() != 'P' )
        throw std::invalid_argument( "case.discretization must be P followed by a positive order" );
    auto const result = std::from_chars( discretization.data() + 1,
                                         discretization.data() + discretization.size(), value );
    if ( result.ec != std::errc{} || result.ptr != discretization.data() + discretization.size() || value < 1 )
        throw std::invalid_argument( "case.discretization must be P followed by a positive order" );
    if ( vm.count( "order" ) )
    {
        auto const order = checkedOrder( "order" );
        if ( !vm.at( "case.discretization" ).defaulted() && value != order.value )
            throw std::invalid_argument( "order and case.discretization specify different polynomial orders" );
        return order;
    }
    return RuntimeOrder::checked( value );
}
}
