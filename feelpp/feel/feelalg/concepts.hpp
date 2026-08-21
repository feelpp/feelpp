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
 * @brief C++20 concepts for Feel++ linear algebra (backends, vectors, matrices)
 * @author Christophe Prud'homme
 * @date 2026-01-02
 */
#ifndef FEELPP_FEELALG_CONCEPTS_HPP
#define FEELPP_FEELALG_CONCEPTS_HPP 1

#include <concepts>
#include <memory>
#include <boost/hana/concept/foldable.hpp>
#include <feel/feelcore/concepts.hpp>
#include <feel/feelcore/traits.hpp>
#include <feel/feelalg/condenser.hpp>
#include <feel/feelalg/productspaceconcepts.hpp>
#include <feel/feelalg/products.hpp>
#include <feel/feeldiscr/traits.hpp>

namespace Feel
{

//
// Backend Concepts
//

/**
 * @brief A linear algebra backend (PETSc, Eigen, Trilinos, etc.)
 * 
 * @details Backends provide vector and matrix types and solve operations.
 * They abstract different linear algebra libraries.
 * 
 * @example
 * @code
 * template <Backend BackendT>
 * void solve_system(std::shared_ptr<BackendT> backend, 
 *                   auto&& A, auto&& b, auto&& x) {
 *     backend->solve(A, b, x);
 * }
 * @endcode
 */
template <typename T>
concept Backend = requires(T t) {
    typename T::vector_type;
    typename T::sparse_matrix_type;
    typename T::graph_type;
    typename T::datamap_type;
};

/**
 * @brief A pointer to a backend
 */
template <typename T>
concept BackendPtr = requires {
    requires (std::is_pointer_v<T> && Backend<std::remove_pointer_t<T>>)
          || (SharedPtr<T> && Backend<typename T::element_type>);
};

//
// Vector Concepts
//

/**
 * @brief A distributed vector
 * 
 * @details Vectors store DOF values in finite element computations.
 * They support parallel distribution via datamaps.
 */
template <typename T>
concept Vector = requires(T t) {
    typename T::value_type;
    typename T::datamap_type;
    { t.size() } -> std::convertible_to<size_t>;
    { t.localSize() } -> std::convertible_to<size_t>;
};

/**
 * @brief A vector pointer (shared_ptr or raw)
 */
template <typename T>
concept VectorPtr = requires {
    requires (std::is_pointer_v<T> && Vector<std::remove_pointer_t<T>>)
          || (SharedPtr<T> && Vector<typename T::element_type>);
};

/**
 * @brief A block vector (concatenation of multiple vectors)
 */
template <typename T>
concept BlockVector = Vector<T> && requires(T t) {
    { t.nBlocks() } -> std::convertible_to<int>;
};

//
// Matrix Concepts
//

/**
 * @brief A sparse matrix
 * 
 * @details Sparse matrices represent finite element operators
 * (stiffness matrices, mass matrices, etc.)
 */
template <typename T>
concept SparseMatrix = requires(T t) {
    typename T::value_type;
    typename T::graph_type;
    { t.size1() } -> std::convertible_to<size_t>;
    { t.size2() } -> std::convertible_to<size_t>;
};

/**
 * @brief A matrix pointer
 */
template <typename T>
concept SparseMatrixPtr = requires {
    requires (std::is_pointer_v<T> && SparseMatrix<std::remove_pointer_t<T>>)
          || (SharedPtr<T> && SparseMatrix<typename T::element_type>);
};

/**
 * @brief A block matrix (block-structured sparse matrix)
 */
template <typename T>
concept BlockMatrix = SparseMatrix<T> && requires(T t) {
    { t.nBlockRows() } -> std::convertible_to<int>;
    { t.nBlockCols() } -> std::convertible_to<int>;
};

//
// Graph Concepts
//

/**
 * @brief A sparsity graph (defines matrix non-zero pattern)
 */
template <typename T>
concept Graph = requires(T t) {
    typename T::datamap_type;
    { t.nDof() } -> std::convertible_to<size_t>;
};

/**
 * @brief A datamap (DOF distribution in parallel)
 */
template <typename T>
concept DataMap = requires(T t) {
    { t.nDof() } -> std::convertible_to<size_t>;
    { t.nLocalDof() } -> std::convertible_to<size_t>;
    { t.nMyElements() } -> std::convertible_to<size_t>;
};

//
// Solver Concepts
//

/**
 * @brief A linear solver
 */
template <typename T>
concept LinearSolver = requires(T t) {
    typename T::backend_type;
    requires Backend<typename T::backend_type>;
};

/**
 * @brief A non-linear solver
 */
template <typename T>
concept NonLinearSolver = requires(T t) {
    typename T::backend_type;
    typename T::jacobian_type;
    typename T::residual_type;
    requires Backend<typename T::backend_type>;
};

/**
 * @brief An eigenvalue solver
 */
template <typename T>
concept EigenSolver = requires(T t) {
    typename T::backend_type;
    requires Backend<typename T::backend_type>;
};

//
// Preconditioner Concepts
//

/**
 * @brief A preconditioner/conditioner for linear systems
 */
template <typename T>
concept Conditioner = requires(T t) {
    typename T::backend_type;
    typename T::operator_type;
    requires Backend<typename T::backend_type>;
};

/**
 * @brief A block preconditioner
 */
template <typename T>
concept BlockConditioner = Conditioner<T> && requires(T t) {
    { t.nBlocks() } -> std::convertible_to<int>;
};

//
// Null Space Concepts
//

/**
 * @brief A null space (kernel) of an operator
 */
template <typename T>
concept NullSpace = requires(T t) {
    typename T::backend_type;
    typename T::vector_type;
    { t.size() } -> std::convertible_to<int>;
};

//
// Index Set Concepts
//

/**
 * @brief An index set (collection of DOF indices)
 */
template <typename T>
concept IndexSet = Iterable<T> && requires(T t) {
    { *t.begin() } -> std::convertible_to<size_t>;
};

//
// Backward Compatibility Bridges
//

#ifdef FEELPP_ENABLE_CONCEPT_COMPATIBILITY

template <typename T>
constexpr bool is_backend_v = Backend<T>;

template <typename T>
constexpr bool is_vector_v = Vector<T>;

template <typename T>
constexpr bool is_sparse_matrix_v = SparseMatrix<T>;

template <typename T>
constexpr bool is_graph_v = Graph<T>;

template <typename T>
constexpr bool is_datamap_v = DataMap<T>;

#endif // FEELPP_ENABLE_CONCEPT_COMPATIBILITY

} // namespace Feel

#endif /* FEELPP_FEELALG_CONCEPTS_HPP */
