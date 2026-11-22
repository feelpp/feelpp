/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELSYSTEM_TAGS_HPP
#define FEELPP_FEELSYSTEM_TAGS_HPP 1

/**
 * @file tags.hpp
 * @brief Empty semantic tags used to classify FeelSystem problem descriptors.
 */

namespace Feel::systems
{

/** @brief Tag for problems assembled as linear algebraic systems. */
struct linear_t {};
/** @brief Tag for residual/Jacobian nonlinear problems. */
struct nonlinear_t {};
/** @brief Tag for Picard fixed-point nonlinear problems. */
struct picard_t {};

/** @brief Tag for stationary problems. */
struct stationary_t {};
/** @brief Tag for first-order evolution problems. */
struct first_order_t {};
/** @brief Tag for second-order evolution problems. */
struct second_order_t {};

/** @brief Tag for elliptic model problems. */
struct elliptic_t {};
/** @brief Tag for parabolic model problems. */
struct parabolic_t {};
/** @brief Tag for hyperbolic model problems. */
struct hyperbolic_t {};

/** @brief Tag for standard or generalized eigenvalue problems. */
struct generalized_eigen_t {};
/** @brief Tag for product-space block problems. */
struct block_t {};
/** @brief Tag for mixed product-space problems. */
struct mixed_t {};
/** @brief Tag for HDG product-space problems. */
struct hdg_t {};

/** @brief Value object for linear problem tagging. */
inline constexpr linear_t linear{};
/** @brief Value object for nonlinear problem tagging. */
inline constexpr nonlinear_t nonlinear{};
/** @brief Value object for Picard problem tagging. */
inline constexpr picard_t picard_tag{};

/** @brief Value object for stationary problem tagging. */
inline constexpr stationary_t stationary{};
/** @brief Value object for first-order evolution problem tagging. */
inline constexpr first_order_t first_order{};
/** @brief Value object for second-order evolution problem tagging. */
inline constexpr second_order_t second_order{};

/** @brief Value object for elliptic problem tagging. */
inline constexpr elliptic_t elliptic{};
/** @brief Value object for parabolic problem tagging. */
inline constexpr parabolic_t parabolic{};
/** @brief Value object for hyperbolic problem tagging. */
inline constexpr hyperbolic_t hyperbolic{};

/** @brief Value object for eigenvalue problem tagging. */
inline constexpr generalized_eigen_t generalized_eigen{};
/** @brief Value object for block problem tagging. */
inline constexpr block_t block{};
/** @brief Value object for mixed problem tagging. */
inline constexpr mixed_t mixed{};
/** @brief Value object for HDG problem tagging. */
inline constexpr hdg_t hdg{};

} // namespace Feel::systems

#endif // FEELPP_FEELSYSTEM_TAGS_HPP
