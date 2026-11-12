/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELSYSTEM_CONCEPTS_HPP
#define FEELPP_FEELSYSTEM_CONCEPTS_HPP 1

/**
 * @file concepts.hpp
 * @brief C++20 concepts selecting FeelSystem problem and system overloads.
 */

#include <concepts>
#include <cstddef>
#include <memory>
#include <string>
#include <type_traits>

#include <feel/feelsystem/tags.hpp>

namespace Feel::systems
{

/**
 * @brief Minimal interface expected from scalar Feel++ function spaces.
 */
template <typename T>
concept FunctionSpaceLike =
    requires( T space ) {
        typename std::remove_cvref_t<T>::element_type;
        typename std::remove_cvref_t<T>::value_type;
        { space.worldCommPtr() };
        { space.element() };
    };

/**
 * @brief True when a descriptor exposes a static `has_tag<Tag>` flag.
 */
template <typename T, typename Tag>
concept TaggedProblemLike =
    requires {
        std::remove_cvref_t<T>::template has_tag<Tag>;
    } && std::remove_cvref_t<T>::template has_tag<Tag>;

/**
 * @brief Common interface shared by scalar problem descriptors.
 */
template <typename T>
concept ProblemLike =
    requires( T const& p ) {
        typename std::remove_cvref_t<T>::space_type;
        typename std::remove_cvref_t<T>::space_ptrtype;
        typename std::remove_cvref_t<T>::value_type;
        { p.space() };
        { p.name() } -> std::convertible_to<std::string>;
        { p.nTerms() } -> std::convertible_to<std::size_t>;
    };

/**
 * @brief Problem descriptor that can be assembled by `LinearSystem`.
 */
template <typename T>
concept LinearProblemLike =
    ProblemLike<T> &&
    TaggedProblemLike<T, linear_t> &&
    requires( T const& p ) {
        { p.lhsTerms() };
        { p.rhsTerms() };
    };

/**
 * @brief Problem descriptor that provides residual and Jacobian terms.
 */
template <typename T>
concept NonlinearProblemLike =
    ProblemLike<T> &&
    TaggedProblemLike<T, nonlinear_t> &&
    requires( T const& p ) {
        { p.residualTerms() };
        { p.jacobianTerms() };
    };

/**
 * @brief Problem descriptor that provides Picard iteration callbacks.
 */
template <typename T>
concept PicardProblemLike =
    ProblemLike<T> &&
    TaggedProblemLike<T, picard_t> &&
    requires( T const& p ) {
        { p.iterationTerms() };
    };

/**
 * @brief Problem descriptor with a first-order time scheme such as BDF.
 */
template <typename T>
concept FirstOrderEvolutionProblemLike =
    ProblemLike<T> &&
    TaggedProblemLike<T, first_order_t> &&
    requires( T const& p ) {
        { p.massTerms() };
        { p.timeSchemeTerms() };
        { p.hasTimeScheme() } -> std::convertible_to<bool>;
        { p.bdfTimeScheme() };
    };

/**
 * @brief Problem descriptor with a second-order time scheme such as Newmark.
 */
template <typename T>
concept SecondOrderEvolutionProblemLike =
    ProblemLike<T> &&
    TaggedProblemLike<T, second_order_t> &&
    requires( T const& p ) {
        { p.massTerms() };
        { p.dampingTerms() };
        { p.stiffnessTerms() };
        { p.timeSchemeTerms() };
        { p.hasTimeScheme() } -> std::convertible_to<bool>;
        { p.newmarkTimeScheme() };
    };

/**
 * @brief First- or second-order evolution problem descriptor.
 */
template <typename T>
concept EvolutionProblemLike =
    FirstOrderEvolutionProblemLike<T> || SecondOrderEvolutionProblemLike<T>;

/**
 * @brief Problem descriptor containing eigen operators A and optional B.
 */
template <typename T>
concept GeneralizedEigenvalueProblemLike =
    ProblemLike<T> &&
    TaggedProblemLike<T, generalized_eigen_t> &&
    requires( T const& p ) {
        { p.operatorATerms() };
        { p.operatorBTerms() };
    };

/**
 * @brief Minimal interface expected from static Feel++ product spaces.
 */
template <typename T>
concept ProductSpaceLike =
    requires( T space ) {
        typename std::remove_cvref_t<T>::value_type;
        typename std::remove_cvref_t<T>::element_type;
        { space.numberOfSpaces() } -> std::convertible_to<int>;
        { space.element() };
    };

/**
 * @brief Product-space descriptor that can be assembled by `BlockSystem`.
 */
template <typename T>
concept BlockProblemLike =
    TaggedProblemLike<T, block_t> &&
    requires( T const& p ) {
        typename std::remove_cvref_t<T>::space_type;
        typename std::remove_cvref_t<T>::value_type;
        typename std::remove_cvref_t<T>::element_type;
        typename std::remove_cvref_t<T>::block_bilinear_form_type;
        typename std::remove_cvref_t<T>::block_linear_form_type;
        requires ProductSpaceLike<typename std::remove_cvref_t<T>::space_type>;
        { p.space() };
        { p.name() } -> std::convertible_to<std::string>;
        { p.nTerms() } -> std::convertible_to<std::size_t>;
        { p.lhsTerms() };
        { p.rhsTerms() };
        { p.constraintTerms() };
    };

/**
 * @brief Block problem descriptor with HDG role spaces and condensation policy.
 */
template <typename T>
concept HdgProblemLike =
    BlockProblemLike<T> &&
    TaggedProblemLike<T, hdg_t> &&
    requires( T const& p ) {
        { p.spaces() };
        { p.condensationPolicy() };
        { p.condensedBackendName() } -> std::convertible_to<std::string>;
        { p.postBackendName() } -> std::convertible_to<std::string>;
    };

/**
 * @brief Already assembled object exposing `assemble`, `assembled`, and `solve`.
 */
template <typename T>
concept SystemLike =
    requires( T system ) {
        { system.assemble() } -> std::same_as<T&>;
        { system.assembled() } -> std::convertible_to<bool>;
        { system.solve() };
    };

} // namespace Feel::systems

#endif // FEELPP_FEELSYSTEM_CONCEPTS_HPP
