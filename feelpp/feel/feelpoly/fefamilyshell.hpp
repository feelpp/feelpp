/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 * @file fefamilyshell.hpp
 * @brief Common metadata for public finite-element family descriptors.
 */

#ifndef FEELPP_FEELPOLY_FEFAMILYSHELL_HPP
#define FEELPP_FEELPOLY_FEFAMILYSHELL_HPP 1

#include <feel/feelpoly/order.hpp>

namespace Feel
{

/**
 * @brief Shared static/runtime-order metadata for FE family descriptors.
 * @tparam Order public family order or `Dynamic`
 * @tparam Tag compatibility tag carried by legacy descriptors
 */
template<int Order, uint16_type Tag>
struct FiniteElementFamilyDescriptor
{
    static_assert( Order == Dynamic || Order >= 0,
                   "finite-element family order must be non-negative or Dynamic" );

    static constexpr bool is_order_static = Order != Dynamic; /**< Whether order is known at compile time. */
    static constexpr bool is_order_dynamic = !is_order_static; /**< Whether order is supplied at construction. */
    static constexpr int order_spec = Order; /**< Unmodified public order specification. */
    static constexpr uint16_type nOrder = is_order_static ? static_cast<uint16_type>( Order ) : 0; /**< Legacy static-order value. */
    static constexpr uint16_type TAG = Tag; /**< Legacy family tag. */
    using order_policy_type = OrderBase<Order>; /**< Storage/evaluation order strategy. */
};

} // namespace Feel

#endif // FEELPP_FEELPOLY_FEFAMILYSHELL_HPP
