/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
       Date: 2022-02-23

  Copyright (C) 2022 Université de Strasbourg

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
#pragma once

#include <fmt/core.h>
#include <fmt/compile.h>
#include <feel/feelcore/feel.hpp>

namespace Feel
{
/**
 * @brief Geometric dimension identifier collection generator
 *
 * @tparam Dmin minimal geometric dimension
 * @tparam Dmax maximal geometric dimension
 */
template <int Dmin = 2, int Dmax = 3, typename = std::enable_if_t<Dmin >= 1 && Dmax <= 3 && Dmin <= Dmax>>
auto dim_t = hana::unpack( hana::make_range( hana::int_c<Dmin>, hana::int_c<Dmax + 1> ), hana::make_tuple );

/**
 * @brief Polynomial order identifier collection generator
 *
 * @tparam Omin minimal polynomial order
 * @tparam Omax maximal polynomial order
 */
template <int Omin = 0, int Omax = FEELPP_INSTANTIATION_ORDER_MAX, typename T = std::enable_if_t<Omin >= 0 && Omax <= 3 && Omin <= Omax>>
auto order_t = hana::unpack( hana::make_range( hana::int_c<Omin>, hana::int_c<Omax + 1> ), hana::make_tuple );

/**
 * @brief Continuous Lagrange Function space identifier collection generator
 *
 * @tparam Dmin minimal geometric dimension
 * @tparam Dmax maximal geometric dimension
 * @tparam Omin minimal polynomial order
 * @tparam Omax maximal polynomial order
 */
template <int Dmin = 2, int Dmax = 3, int Omin = 1, int Omax = FEELPP_INSTANTIATION_ORDER_MAX>
auto Pc_t = hana::transform( hana::cartesian_product( hana::make_tuple( dim_t<Dmin, Dmax>, order_t<Omin, Omax> ) ), []( auto x )
                              { return hana::append( x, fmt::format( FMT_COMPILE( "P{}" ), std::decay_t<decltype( hana::at_c<1>( x ) )>::value ) ); } );

/**
 * @brief Discontinuous Lagrange Function space identifier collection generator
 *
 * @tparam Dmin minimal geometric dimension
 * @tparam Dmax maximal geometric dimension
 * @tparam Omin minimal polynomial order
 * @tparam Omax maximal polynomial order
 */
template <int Dmin = 2, int Dmax = 3, int Omin = 1, int Omax = FEELPP_INSTANTIATION_ORDER_MAX>
auto Pd_t = hana::transform( hana::cartesian_product( hana::make_tuple( dim_t<Dmin, Dmax>, order_t<Omin, Omax> ) ), []( auto x )
                              { return hana::append( x, fmt::format( FMT_COMPILE( "Pd{}" ), std::decay_t<decltype( hana::at_c<1>( x ) )>::value ) ); } );

} // namespace Feel
