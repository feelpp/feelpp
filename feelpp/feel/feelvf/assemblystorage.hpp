/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/** @file assemblystorage.hpp
 * @brief Checked local assembly shapes using the common Eigen storage policy.
 */
#pragma once

#include <algorithm>
#include <feel/feelpoly/eigenextent.hpp>

namespace Feel::detail
{
/** Existing assembly inline ceiling in bytes; Eigen's configured ceiling also applies. */
inline constexpr std::size_t assemblyInlineMatrixBytes = 128 * 128 * 8;

/** @brief Remove a mortar constraint and/or repeat a local extent without overflow.
 * Empty extents remain empty. Dynamic is a type-level sentinel only.
 */
[[nodiscard]] constexpr int assemblyExtent( int extent, int removed = 0, int cells = 1 )
{
    if ( removed < 0 || cells < 1 || extent < Dynamic )
        throw std::invalid_argument( "invalid local assembly extent adjustment" );
    if ( extent == Dynamic ) return Dynamic;
    int const remaining = std::max( extent - removed, 0 );
    if ( remaining > std::numeric_limits<int>::max() / cells )
        throw std::length_error( "local assembly extent exceeds int" );
    return cells * remaining;
}

/** @brief Check an FE's runtime local count before narrowing to Eigen's assembly extent. */
[[nodiscard]] constexpr int assemblyDofCount( std::size_t count )
{
    return staticEigenExtent( count );
}

/** @brief One shape vocabulary for ordinary, mortar, component and two-cell forms.
 * @tparam T Matrix scalar.
 * @tparam Test Fixed test count or Dynamic, independently of the trial side.
 * @tparam Trial Fixed trial count or Dynamic.
 * @tparam InlineBytes Existing assembly coefficient budget.
 */
template <typename T, int Test, int Trial, std::size_t InlineBytes = assemblyInlineMatrixBytes>
struct AssemblyStorage
{
    /** Matrix representation selected by the common polynomial/assembly byte policy. */
    template <int Rows, int Cols>
    using Matrix = typename MatrixStoragePolicy<T, Rows, Cols, InlineBytes>::matrix_type;
    /** Ordinary test-by-trial element contribution. */
    using Ordinary = Matrix<Test, Trial>;
    /** One eliminated test constraint, including an empty P0 mortar row space. */
    using MortarTest = Matrix<assemblyExtent( Test, 1 ), Trial>;
    /** One eliminated trial constraint, with legal Eigen vector storage order. */
    using MortarTrial = Matrix<Test, assemblyExtent( Trial, 1 )>;
    /** Four blocks for a pair of adjacent elements. */
    using TwoCell = Matrix<assemblyExtent( Test, 0, 2 ), assemblyExtent( Trial, 0, 2 )>;
    /** Test orientation signs, budgeted by their own coefficient count. */
    using RowSigns = typename MatrixStoragePolicy<int, Test, 1, InlineBytes>::matrix_type;
    /** Trial orientation signs, retaining a known trial extent in mixed forms. */
    using ColSigns = typename MatrixStoragePolicy<int, Trial, 1, InlineBytes>::matrix_type;
    /** Test signs for both cells. */
    using TwoCellRowSigns = typename MatrixStoragePolicy<int, assemblyExtent( Test, 0, 2 ), 1, InlineBytes>::matrix_type;
    /** Trial signs for both cells. */
    using TwoCellColSigns = typename MatrixStoragePolicy<int, assemblyExtent( Trial, 0, 2 ), 1, InlineBytes>::matrix_type;
};

/** @brief Validate allocation arithmetic and fixed sides before resizing any local matrix.
 * Eigen resize preserves storage for an unchanged shape. Validation also runs
 * in Release builds, where an Eigen assertion cannot enforce a fixed-axis contract.
 */
template <typename Matrix>
void resizeAssemblyMatrix( Matrix& matrix, Eigen::Index rows, Eigen::Index cols = 1 )
{
    static_cast<void>( matrixElementCount( rows, cols, sizeof( typename Matrix::Scalar ) ) );
    if constexpr ( Matrix::RowsAtCompileTime != Dynamic )
        if ( rows != Matrix::RowsAtCompileTime )
            throw std::invalid_argument( "local assembly static row mismatch" );
    if constexpr ( Matrix::ColsAtCompileTime != Dynamic )
        if ( cols != Matrix::ColsAtCompileTime )
            throw std::invalid_argument( "local assembly static column mismatch" );
    matrix.resize( rows, cols );
}
} // namespace Feel::detail
