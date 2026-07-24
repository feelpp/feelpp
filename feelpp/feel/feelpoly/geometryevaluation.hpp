/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file geometryevaluation.hpp
 * @brief Evaluation strategies for geometric-mapping context kernels.
 */
#if !defined( FEELPP_FEELPOLY_GEOMETRYEVALUATION_HPP )
#define FEELPP_FEELPOLY_GEOMETRYEVALUATION_HPP 1

#include <concepts>

namespace Feel::geometry
{

/**
 * @enum geometry_evaluation_strategy
 * @brief Select the local kernel used to evaluate geometry derivatives.
 *
 * `automatic` selects the batched kernel when several non-affine evaluation
 * points are available. `pointwise` and `batched` are primarily useful for
 * validation and performance characterization; both implement the same
 * mathematical mapping.
 */
enum class geometry_evaluation_strategy
{
    automatic,
    pointwise,
    batched
};

/** @brief Compile-time strategy for an affine geometric mapping. */
struct affine_geometry_order
{
    /** Affine maps have constant first derivatives and zero second derivatives. */
    static constexpr bool is_affine = true;
};

/** @brief Compile-time strategy for a non-affine geometric mapping. */
struct high_order_geometry
{
    /** High-order maps may vary at every evaluation point. */
    static constexpr bool is_affine = false;
};

/**
 * @concept geometric_order_strategy
 * @brief Strategy accepted by a geometric context update kernel.
 */
template <typename Strategy>
concept geometric_order_strategy = requires {
    { Strategy::is_affine } -> std::convertible_to<bool>;
};

} // namespace Feel::geometry

#endif // FEELPP_FEELPOLY_GEOMETRYEVALUATION_HPP
