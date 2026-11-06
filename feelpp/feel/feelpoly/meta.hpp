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
 * @file meta.hpp
 * @brief Metaprogramming utilities for feelpoly (mp11-based).
 * @author Christophe Prud'homme
 * @date 2026-01-02
 */
#ifndef FEELPP_FEELPOLY_META_HPP
#define FEELPP_FEELPOLY_META_HPP 1

#include <cstddef>
#include <type_traits>

#include <boost/mp11.hpp>

#include <feel/feelcore/traits.hpp>

namespace Feel
{
namespace mp = boost::mp11;

/**
 * @brief mp11-backed type list.
 */
template <typename... Ts>
using type_list = mp::mp_list<Ts...>;

/**
 * @brief Access the Ith type in a type list.
 */
template <typename L, std::size_t I>
using type_at = mp::mp_at_c<L, I>;

/**
 * @brief Find the index of T in a type list.
 */
template <typename L, typename T>
using type_find = mp::mp_find<L, T>;

/**
 * @brief Transform a type list with a unary metafunction.
 */
template <typename L, template <class...> class F>
using type_transform = mp::mp_transform<F, L>;

/**
 * @brief Conditional type selection.
 */
template <bool Cond, typename T, typename F>
using if_t = std::conditional_t<Cond, T, F>;

/**
 * @brief Integral constant helper.
 *
 * Note: defined in Feel::meta to avoid collision with ExtrapolationType_P1::constant.
 */
namespace meta
{
template <auto V>
using constant = std::integral_constant<decltype(V), V>;
} // namespace meta

/**
 * @brief int constant helper.
 */
template <int N>
using int_c = std::integral_constant<int, N>;

/**
 * @brief uint16_type constant helper.
 */
template <uint16_type N>
using uint16_c = std::integral_constant<uint16_type, N>;

/**
 * @brief bool constant helper.
 */
template <bool B>
using bool_c = std::integral_constant<bool, B>;

} // namespace Feel

#endif /* FEELPP_FEELPOLY_META_HPP */
