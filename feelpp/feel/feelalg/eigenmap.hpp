/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file eigenmap.hpp
 * @brief Const-correct zero-copy Eigen views of contiguous uBLAS storage.
 */
#if !defined( FEELPP_FEELALG_EIGENMAP_HPP )
#define FEELPP_FEELALG_EIGENMAP_HPP 1

#include <Eigen/Core>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <concepts>
#include <memory>
#include <type_traits>

namespace Feel
{
namespace ublas = boost::numeric::ublas;

/** @brief Identify owning dense uBLAS matrices supported by `eigenMap()`. */
template <typename Type>
struct is_ublas_dense_matrix : std::false_type
{
};

/** @brief Dense uBLAS matrix specialization. */
template <typename T, typename Layout, typename Storage>
struct is_ublas_dense_matrix<ublas::matrix<T, Layout, Storage>>
    : std::true_type
{
};

/** @brief Identify owning dense uBLAS vectors supported by `eigenMap()`. */
template <typename Type>
struct is_ublas_dense_vector : std::false_type
{
};

/** @brief Dense uBLAS vector specialization. */
template <typename T, typename Storage>
struct is_ublas_dense_vector<ublas::vector<T, Storage>> : std::true_type
{
};

/**
 * @concept UblasDenseMatrix
 * @brief Owning contiguous uBLAS matrix accepted by the Eigen adapter.
 */
template <typename Type>
concept UblasDenseMatrix =
    is_ublas_dense_matrix<std::remove_cvref_t<Type>>::value;

/**
 * @concept UblasDenseVector
 * @brief Owning contiguous uBLAS vector accepted by the Eigen adapter.
 */
template <typename Type>
concept UblasDenseVector =
    is_ublas_dense_vector<std::remove_cvref_t<Type>>::value;

/**
 * @brief Create a mutable zero-copy Eigen view of a dense uBLAS matrix.
 *
 * The Eigen storage order follows the uBLAS layout. The returned unaligned
 * view is valid while the matrix storage remains alive and is not resized.
 *
 * @tparam T scalar type
 * @tparam Layout `ublas::row_major` or `ublas::column_major`
 * @tparam Storage contiguous uBLAS storage type
 * @param matrix owning uBLAS matrix
 * @return mutable Eigen view with matching dimensions and layout
 */
template <typename T, typename Layout, typename Storage>
    requires std::same_as<Layout, ublas::row_major> ||
             std::same_as<Layout, ublas::column_major>
[[nodiscard]] auto
eigenMap( ublas::matrix<T, Layout, Storage>& matrix ) noexcept
{
    constexpr int storageOrder =
        std::same_as<Layout, ublas::column_major> ? Eigen::ColMajor
                                                  : Eigen::RowMajor;
    using matrix_type =
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, storageOrder>;
    return Eigen::Map<matrix_type, Eigen::Unaligned>(
        std::to_address( matrix.data().begin() ), matrix.size1(),
        matrix.size2() );
}

/**
 * @brief Create a read-only zero-copy Eigen view of a dense uBLAS matrix.
 * @copydetails eigenMap(ublas::matrix<T, Layout, Storage>&)
 */
template <typename T, typename Layout, typename Storage>
    requires std::same_as<Layout, ublas::row_major> ||
             std::same_as<Layout, ublas::column_major>
[[nodiscard]] auto
eigenMap( ublas::matrix<T, Layout, Storage> const& matrix ) noexcept
{
    constexpr int storageOrder =
        std::same_as<Layout, ublas::column_major> ? Eigen::ColMajor
                                                  : Eigen::RowMajor;
    using matrix_type =
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, storageOrder>;
    return Eigen::Map<matrix_type const, Eigen::Unaligned>(
        std::to_address( matrix.data().begin() ), matrix.size1(),
        matrix.size2() );
}

/**
 * @brief Create a mutable zero-copy Eigen column-vector view of a uBLAS vector.
 *
 * @tparam T scalar type
 * @tparam Storage contiguous uBLAS storage type
 * @param vector owning uBLAS vector
 * @return mutable unaligned Eigen vector view
 */
template <typename T, typename Storage>
[[nodiscard]] auto eigenMap( ublas::vector<T, Storage>& vector ) noexcept
{
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    return Eigen::Map<vector_type, Eigen::Unaligned>(
        std::to_address( vector.data().begin() ), vector.size() );
}

/**
 * @brief Create a read-only zero-copy Eigen column-vector view of a uBLAS vector.
 * @copydetails eigenMap(ublas::vector<T, Storage>&)
 */
template <typename T, typename Storage>
[[nodiscard]] auto
eigenMap( ublas::vector<T, Storage> const& vector ) noexcept
{
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    return Eigen::Map<vector_type const, Eigen::Unaligned>(
        std::to_address( vector.data().begin() ), vector.size() );
}

} // namespace Feel

#endif // FEELPP_FEELALG_EIGENMAP_HPP
