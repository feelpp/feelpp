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
 * @brief C++20 concepts for Feel++ discretization (function spaces, elements, DOFs)
 * @author Christophe Prud'homme
 * @date 2026-01-02
 */
#ifndef FEELPP_FEELDISCR_CONCEPTS_HPP
#define FEELPP_FEELDISCR_CONCEPTS_HPP 1

#include <concepts>
#include <memory>
#include <feel/feelcore/concepts.hpp>

namespace Feel
{

//
// Function Space Concepts
//

/**
 * @brief A Feel++ function space
 * 
 * @details Function spaces are the fundamental discretization concept in Feel++.
 * They define the approximation space for finite element methods.
 * 
 * Requirements:
 * - value_type: The scalar type (typically double)
 * - mesh_type: The associated mesh type
 * - element_type: The finite element function type
 * - element(): Factory method to create elements
 * 
 * @example
 * @code
 * template <FunctionSpace SpaceT>
 * auto project(VfExpr auto&& expr, std::shared_ptr<SpaceT> space) {
 *     auto u = space->element();
 *     // project expr onto u
 *     return u;
 * }
 * @endcode
 */
template <typename T>
concept FunctionSpace = requires(T t) {
    typename T::value_type;
    typename T::mesh_type;
    typename T::element_type;
    { t.element() } -> std::same_as<typename T::element_type>;
};

/**
 * @brief A pointer (shared_ptr or raw) to a function space
 */
template <typename T>
concept FunctionSpacePtr = requires {
    requires (std::is_pointer_v<T> && FunctionSpace<std::remove_pointer_t<T>>)
          || (SharedPtr<T> && FunctionSpace<typename T::element_type>);
};

/**
 * @brief A function space element (finite element function)
 * 
 * @details Elements represent functions in the discretized space.
 * They have values at DOFs and can be evaluated at points.
 */
template <typename T>
concept FunctionSpaceElement = requires(T t) {
    typename T::functionspace_type;
    typename T::value_type;
    requires FunctionSpace<typename T::functionspace_type>;
};

/**
 * @brief A product of function spaces (for mixed formulations)
 * 
 * @details Product spaces combine multiple function spaces, used in
 * mixed finite element methods (e.g., Stokes: velocity × pressure).
 */
template <typename T>
concept ProductSpace = requires(T t) {
    typename T::spaces_tuple_type;
    { t.numberOfSpaces() } -> std::convertible_to<int>;
};

/**
 * @brief Multiple product spaces (product of product spaces)
 */
template <typename T>
concept ProductSpaces = ProductSpace<T> && requires {
    typename T::spaces_array_type;
};

//
// Mesh Concepts
//

/**
 * @brief A Feel++ mesh
 * 
 * @details Meshes represent the geometric domain discretization.
 */
template <typename T>
concept Mesh = requires(T t) {
    typename T::shape_type;
    typename T::element_type;
    typename T::face_type;
    { t.numElements() } -> std::convertible_to<size_t>;
};

/**
 * @brief A pointer to a mesh
 */
template <typename T>
concept MeshPtr = requires {
    requires (std::is_pointer_v<T> && Mesh<std::remove_pointer_t<T>>)
          || (SharedPtr<T> && Mesh<typename T::element_type>);
};

//
// DOF Concepts
//

/**
 * @brief A degree of freedom (DOF) identifier
 * 
 * @details DOFs represent discrete unknowns in the finite element system.
 */
template <typename T>
concept DofType = requires(T t) {
    { t.index() } -> std::convertible_to<size_t>;
    { t.sign() } -> std::convertible_to<int>;
};

/**
 * @brief A DOF table mapping elements to DOFs
 */
template <typename T>
concept DofTable = requires(T t) {
    typename T::dof_type;
    requires DofType<typename T::dof_type>;
};

//
// Geometric Mapping Concepts
//

/**
 * @brief A geometric mapping from reference to physical element
 */
template <typename T>
concept GeometricMapping = requires {
    typename T::gm_type;
    { T::nDim } -> std::convertible_to<int>;
    { T::nRealDim } -> std::convertible_to<int>;
};

/**
 * @brief A geometric mapping context (evaluation at quadrature points)
 */
template <typename T>
concept GeometricMappingContext = requires(T t) {
    typename T::gm_type;
    requires GeometricMapping<typename T::gm_type>;
};

//
// Basis Function Concepts
//

/**
 * @brief A finite element basis (shape functions)
 */
template <typename T>
concept Basis = requires {
    typename T::value_type;
    { T::nDof } -> std::convertible_to<int>;
    { T::nLocalDof } -> std::convertible_to<int>;
};

/**
 * @brief A basis context (evaluation at points)
 */
template <typename T>
concept BasisContext = requires(T t) {
    typename T::basis_type;
    requires Basis<typename T::basis_type>;
};

//
// Interpolation Concepts
//

/**
 * @brief A type that can be interpolated (projected) onto a function space
 */
template <typename T, typename SpaceT>
concept Interpolable = requires(T t, SpaceT space) {
    requires FunctionSpace<SpaceT>;
    // Can be evaluated to produce values
};

//
// Operator Concepts
//

/**
 * @brief A linear operator between function spaces
 */
template <typename T>
concept LinearOperator = requires(T t) {
    typename T::domain_space_type;
    typename T::dual_image_space_type;
    requires FunctionSpace<typename T::domain_space_type>;
    requires FunctionSpace<typename T::dual_image_space_type>;
};

/**
 * @brief A preconditioner for linear systems
 */
template <typename T>
concept Preconditioner = requires(T t) {
    typename T::operator_type;
    { t.setOperator(std::declval<typename T::operator_type>()) };
};

//
// Backward Compatibility Bridges
//

#ifdef FEELPP_ENABLE_CONCEPT_COMPATIBILITY

// Bridge old is_functionspace trait to new concept
template <typename T>
constexpr bool is_functionspace_v = FunctionSpace<T>;

template <typename T>
constexpr bool is_functionspace_element_v = FunctionSpaceElement<T>;

template <typename T>
constexpr bool is_product_space_v = ProductSpace<T>;

#endif // FEELPP_ENABLE_CONCEPT_COMPATIBILITY

} // namespace Feel

#endif /* FEELPP_FEELDISCR_CONCEPTS_HPP */
