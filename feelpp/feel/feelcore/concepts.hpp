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
 * @file concepts.hpp
 * @brief C++20 concepts for Feel++ core types
 * @author Christophe Prud'homme
 * @date 2026-01-02
 * 
 * This header provides modern C++20 concept definitions to replace SFINAE-based
 * type traits. These concepts provide:
 * - Clearer, more self-documenting code
 * - Better compiler error messages
 * - Faster compilation (concepts checked once, not repeatedly)
 * - Improved IDE support and code navigation
 * 
 * @note Requires C++20 or later.
 */
#ifndef FEELPP_CORE_CONCEPTS_HPP
#define FEELPP_CORE_CONCEPTS_HPP 1

// Include all standard library headers BEFORE entering Feel namespace
// to avoid conflicts with Feel::std namespace
#include <type_traits>
#include <concepts>
#include <iterator>
// Note: <ranges> not included to avoid conflicts with Feel::std namespace

// DO NOT wrap in namespace Feel {} - this file is included from traits.hpp
// which is already inside namespace Feel {}. Adding namespace here creates
// nested Feel::Feel:: and breaks qualified name lookups.

//
// Container Concepts
//

/**
 * @brief A type that can be iterated with begin() and end()
 * 
 * @details An Iterable type provides begin() and end() methods that return
 * iterators. This is commonly satisfied by standard containers, ranges, and
 * many Feel++ mesh entities.
 * 
 * @tparam T Type to check
 * 
 * @example
 * @code
 * template <Iterable Container>
 * void process(Container const& c) {
 *     for(auto const& elem : c) {
 *         // process elem
 *     }
 * }
 * @endcode
 */
template <typename T>
concept Iterable = requires(T t) {
    { t.begin() } -> std::input_or_output_iterator;
    { t.end() } -> std::sentinel_for<decltype(t.begin())>;
};

/**
 * @brief An Iterable type whose elements are convertible to a specific type
 * 
 * @tparam T Container type
 * @tparam V Element value type
 * 
 * @example
 * @code
 * template <IterableOf<double> Container>
 * double sum(Container const& c) {
 *     double total = 0;
 *     for(auto val : c) total += val;
 *     return total;
 * }
 * @endcode
 */
template <typename T, typename V>
concept IterableOf = Iterable<T> && requires(T t) {
    { *t.begin() } -> std::convertible_to<V>;
};

/**
 * @brief A standard vector type
 */
template <typename T>
concept StdVector = requires {
    typename T::value_type;
    requires std::same_as<T, std::vector<typename T::value_type>>;
};

//
// Numeric Concepts
//

/**
 * @brief A type that represents a scalar numeric value
 * 
 * @details Scalar types are fundamental numeric types (int, float, double, etc.)
 * or complex types. This concept is used to distinguish scalars from vectors,
 * tensors, and matrices in Feel++ expressions.
 * 
 * @note Named ScalarConcept to avoid collision with struct Scalar in feelpoly/policy.hpp
 */
template <typename T>
concept ScalarConcept = std::is_arithmetic_v<T> || 
                        (requires { typename T::value_type; } && std::is_arithmetic_v<typename T::value_type>);

/**
 * @brief A type that can be added
 */
template <typename T>
concept Addable = requires(T a, T b) {
    { a + b } -> std::convertible_to<T>;
};

/**
 * @brief A type that can be multiplied
 */
template <typename T>
concept Multipliable = requires(T a, T b) {
    { a * b } -> std::convertible_to<T>;
};

/**
 * @brief A type that forms a mathematical field (can add, subtract, multiply, divide)
 * 
 * @note Named FieldConcept to avoid collision with struct Field in feelpoly/policy.hpp
 */
template <typename T>
concept FieldConcept = requires(T a, T b) {
    { a + b } -> std::convertible_to<T>;
    { a - b } -> std::convertible_to<T>;
    { a * b } -> std::convertible_to<T>;
    { a / b } -> std::convertible_to<T>;
};

//
// Function Space Related Concepts
//

/**
 * @brief A type that has a value_type member
 * 
 * @details Most Feel++ types (function spaces, elements, expressions) define
 * a value_type indicating their scalar type (typically double or float).
 */
template <typename T>
concept HasValueType = requires {
    typename T::value_type;
};

/**
 * @brief A type that has a mesh_type member
 * 
 * @details Function spaces, elements, and many other Feel++ types are associated
 * with a mesh type.
 */
template <typename T>
concept HasMeshType = requires {
    typename T::mesh_type;
};

/**
 * @brief A type that provides world communicator access
 * 
 * @details Most Feel++ objects support MPI parallelism and provide worldComm()
 * or worldCommPtr() methods.
 */
template <typename T>
concept HasWorldComm = requires(T t) {
    { t.worldComm() };
} || requires(T t) {
    { t.worldCommPtr() };
};

//
// Smart Pointer Concepts
//

/**
 * @brief A shared_ptr to any type
 */
template <typename T>
concept SharedPtr = requires {
    typename T::element_type;
    requires std::same_as<T, std::shared_ptr<typename T::element_type>>;
};

/**
 * @brief Either a raw pointer or shared_ptr to a type
 */
template <typename T, typename U>
concept PtrTo = (std::is_pointer_v<T> && std::is_same_v<std::remove_pointer_t<T>, U>)
             || (SharedPtr<T> && std::is_same_v<typename T::element_type, U>);

//
// Expression Traits
//

/**
 * @brief A type that can be evaluated without a geometric context
 * 
 * @details Some Feel++ expressions (constants, parameters) can be evaluated
 * without needing mesh element information.
 */
template <typename T>
concept EvaluableWithoutContext = requires(T t) {
    { t.evaluate(true) };
};

//
// Metaprogramming Helpers
//

/**
 * @brief Check if a type is derived from a base (after decay)
 * 
 * @note Uses decay_type defined in traits.hpp
 */
template <typename Derived, typename Base>
concept DerivedFrom = std::derived_from<std::decay_t<Derived>, Base>;

/**
 * @brief Check if a type is the same as another (after decay)
 */
template <typename T, typename U>
concept SameAs = std::same_as<std::decay_t<T>, std::decay_t<U>>;

// Note: Backward compatibility bridges (is_iterable_v, is_iterable_of_v, etc.)
// are provided by traits.hpp which includes this file. No need to duplicate them here.

#endif /* FEELPP_CORE_CONCEPTS_HPP */
