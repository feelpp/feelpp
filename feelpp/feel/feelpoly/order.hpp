/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 2026-01-02

 Copyright (C) 2026 Feel++ Consortium

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
/**
 * @file order.hpp
 * @brief Order tags and DOF calculators for static/dynamic orders.
 * @author Christophe Prud'homme
 * @date 2026-01-02
 */
#ifndef FEELPP_FEELPOLY_ORDER_HPP
#define FEELPP_FEELPOLY_ORDER_HPP 1

#include <cstddef>

#include <feel/feelcore/feeltypes.hpp>
#include <feel/feelpoly/concepts.hpp>

namespace Feel
{

/**
 * @brief Order tag with static value.
 */
template <int N>
struct order_t
{
    static constexpr int value = N;
    static constexpr bool is_dynamic = (N == Dynamic);
    static constexpr bool is_static = !is_dynamic;
};

namespace detail
{
constexpr size_type pow_int( size_type base, int exp )
{
    size_type result = 1;
    for ( int i = 0; i < exp; ++i )
        result *= base;
    return result;
}

constexpr size_type binomial( int n, int k )
{
    if ( k < 0 || k > n )
        return 0;
    if ( k > n - k )
        k = n - k;
    size_type result = 1;
    for ( int i = 1; i <= k; ++i )
        result = result * static_cast<size_type>( n - k + i ) / static_cast<size_type>( i );
    return result;
}

template <int Dim, int Order>
constexpr size_type simplex_dof_static()
{
    static_assert( Dim >= 0, "Dimension must be non-negative" );
    static_assert( Order >= 0, "Order must be non-negative" );
    return binomial( Dim + Order, Dim );
}

template <int Dim, int Order>
constexpr size_type hypercube_dof_static()
{
    static_assert( Dim >= 0, "Dimension must be non-negative" );
    static_assert( Order >= 0, "Order must be non-negative" );
    return pow_int( static_cast<size_type>( Order + 1 ), Dim );
}

template <int Dim>
inline size_type simplex_dof_dynamic( int order )
{
    if ( order < 0 )
        return 0;
    return binomial( Dim + order, Dim );
}

template <int Dim>
inline size_type hypercube_dof_dynamic( int order )
{
    if ( order < 0 )
        return 0;
    return pow_int( static_cast<size_type>( order + 1 ), Dim );
}
} // namespace detail

/**
 * @brief DOF calculator for static or dynamic orders.
 */
template <int Dim, typename Order>
struct DofCalculator;

/**
 * @brief Static order specialization.
 */
template <int Dim, int Order>
struct DofCalculator<Dim, order_t<Order>>
{
    static constexpr size_type simplex()
    {
        return detail::simplex_dof_static<Dim, Order>();
    }

    static constexpr size_type hypercube()
    {
        return detail::hypercube_dof_static<Dim, Order>();
    }
};

/**
 * @brief Dynamic order specialization.
 */
template <int Dim>
struct DofCalculator<Dim, order_t<Dynamic>>
{
    static size_type simplex( int order )
    {
        return detail::simplex_dof_dynamic<Dim>( order );
    }

    static size_type hypercube( int order )
    {
        return detail::hypercube_dof_dynamic<Dim>( order );
    }
};

} // namespace Feel

#endif /* FEELPP_FEELPOLY_ORDER_HPP */
