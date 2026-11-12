/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELSYSTEM_NAMEDARGUMENTS_HPP
#define FEELPP_FEELSYSTEM_NAMEDARGUMENTS_HPP 1

/**
 * @file namedarguments.hpp
 * @brief Small helpers for inspecting Feel++ named-argument packs.
 */

#include <concepts>
#include <type_traits>
#include <utility>

#include <feel/feelcore/parameter.hpp>

namespace Feel::systems::detail
{

/**
 * @brief True when `Arg` is a Feel++ named argument with the same tag as `NamedArg`.
 */
template <typename Arg, typename NamedArg>
inline constexpr bool isNamedArgumentWithTag =
    NA::is_named_argument_v<std::remove_cvref_t<Arg>> &&
    std::same_as<typename std::remove_cvref_t<Arg>::tag_type,
                 typename NamedArg::tag_type>;

/**
 * @brief True when an argument pack contains a named argument matching `NamedArg`.
 */
template <typename NamedArg, typename... Args>
inline constexpr bool hasNamedArgument =
    ( isNamedArgumentWithTag<Args, NamedArg> || ... || false );

/**
 * @brief Return the value carried by the first named argument matching `NamedArg`.
 *
 * This helper is used after `hasNamedArgument` has selected a valid overload.
 */
template <typename NamedArg, typename First, typename... Rest>
decltype(auto)
namedArgumentValue( First&& first, Rest&&... rest )
{
    if constexpr ( isNamedArgumentWithTag<First, NamedArg> )
        return std::forward<First>( first ).value();
    else
        return namedArgumentValue<NamedArg>( std::forward<Rest>( rest )... );
}

} // namespace Feel::systems::detail

#endif // FEELPP_FEELSYSTEM_NAMEDARGUMENTS_HPP
