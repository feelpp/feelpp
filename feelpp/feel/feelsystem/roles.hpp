/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELSYSTEM_ROLES_HPP
#define FEELPP_FEELSYSTEM_ROLES_HPP 1

/**
 * @file roles.hpp
 * @brief Semantic role tags for block and HDG product-space components.
 */

#include <concepts>
#include <string_view>
#include <type_traits>

namespace Feel::systems
{

template <typename Role>
struct block_role_traits;

/**
 * @brief Concept for semantic block unknown/equation roles.
 *
 * A block role names one component of a product-space problem, for example a
 * flux, potential, trace, stress, velocity, or pressure block. It is distinct
 * from `term_role`, which names the algebraic purpose of a stored contribution
 * such as lhs, rhs, residual, or mass.
 */
template <typename Role>
concept BlockRole =
    requires {
        { block_role_traits<std::remove_cvref_t<Role>>::name } -> std::convertible_to<std::string_view>;
    };

/**
 * @brief Human-readable name for a semantic block role.
 */
template <BlockRole Role>
constexpr std::string_view
blockRoleName( Role ) noexcept
{
    return block_role_traits<std::remove_cvref_t<Role>>::name;
}

namespace roles
{

/** @brief Semantic role for flux unknowns/equations. */
struct flux_t {};
/** @brief Semantic role for scalar potential unknowns/equations. */
struct potential_t {};
/** @brief Semantic role for trace/skeleton unknowns/equations. */
struct trace_t {};
/** @brief Semantic role for stress unknowns/equations. */
struct stress_t {};
/** @brief Semantic role for velocity unknowns/equations. */
struct velocity_t {};
/** @brief Semantic role for pressure unknowns/equations. */
struct pressure_t {};
/** @brief Semantic role for displacement unknowns/equations. */
struct displacement_t {};
/** @brief Semantic role for constraint unknowns/equations. */
struct constraint_t {};

/** @brief Value object for the flux role. */
inline constexpr flux_t flux{};
/** @brief Value object for the potential role. */
inline constexpr potential_t potential{};
/** @brief Value object for the trace role. */
inline constexpr trace_t trace{};
/** @brief Value object for the stress role. */
inline constexpr stress_t stress{};
/** @brief Value object for the velocity role. */
inline constexpr velocity_t velocity{};
/** @brief Value object for the pressure role. */
inline constexpr pressure_t pressure{};
/** @brief Value object for the displacement role. */
inline constexpr displacement_t displacement{};
/** @brief Value object for the constraint role. */
inline constexpr constraint_t constraint{};

} // namespace roles

template <>
struct block_role_traits<roles::flux_t>
{
    /** @brief Human-readable role name. */
    inline static constexpr std::string_view name = "flux";
};

template <>
struct block_role_traits<roles::potential_t>
{
    /** @brief Human-readable role name. */
    inline static constexpr std::string_view name = "potential";
};

template <>
struct block_role_traits<roles::trace_t>
{
    /** @brief Human-readable role name. */
    inline static constexpr std::string_view name = "trace";
};

template <>
struct block_role_traits<roles::stress_t>
{
    /** @brief Human-readable role name. */
    inline static constexpr std::string_view name = "stress";
};

template <>
struct block_role_traits<roles::velocity_t>
{
    /** @brief Human-readable role name. */
    inline static constexpr std::string_view name = "velocity";
};

template <>
struct block_role_traits<roles::pressure_t>
{
    /** @brief Human-readable role name. */
    inline static constexpr std::string_view name = "pressure";
};

template <>
struct block_role_traits<roles::displacement_t>
{
    /** @brief Human-readable role name. */
    inline static constexpr std::string_view name = "displacement";
};

template <>
struct block_role_traits<roles::constraint_t>
{
    /** @brief Human-readable role name. */
    inline static constexpr std::string_view name = "constraint";
};

} // namespace Feel::systems

#endif // FEELPP_FEELSYSTEM_ROLES_HPP
