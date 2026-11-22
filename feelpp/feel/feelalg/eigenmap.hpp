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

// clang-format off
#include <feel/feelcore/warnoff.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <feel/feelcore/warnon.hpp>
// clang-format on

#include <concepts>
#include <memory>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace Feel
{
namespace ublas = boost::numeric::ublas;

/** @brief Identify owning dense uBLAS matrices supported by `eigenMap()`. */
template <typename Type>
struct IsUblasDenseMatrix : std::false_type
{
};

/** @brief Dense uBLAS matrix specialization. */
template <typename T, typename Layout, typename Storage>
struct IsUblasDenseMatrix<ublas::matrix<T, Layout, Storage>>
    : std::true_type
{
};

/** @brief Identify owning dense uBLAS vectors supported by `eigenMap()`. */
template <typename Type>
struct IsUblasDenseVector : std::false_type
{
};

/** @brief Dense uBLAS vector specialization. */
template <typename T, typename Storage>
struct IsUblasDenseVector<ublas::vector<T, Storage>> : std::true_type
{
};

/** Compatibility spelling for existing matrix-trait users. */
template<typename Type> using is_ublas_dense_matrix = IsUblasDenseMatrix<Type>;
/** Compatibility spelling for existing vector-trait users. */
template<typename Type> using is_ublas_dense_vector = IsUblasDenseVector<Type>;

/**
 * @concept UblasDenseMatrix
 * @brief Owning contiguous uBLAS matrix accepted by the Eigen adapter.
 */
template <typename Type>
concept UblasDenseMatrix =
    IsUblasDenseMatrix<std::remove_cvref_t<Type>>::value;

/**
 * @concept UblasDenseVector
 * @brief Owning contiguous uBLAS vector accepted by the Eigen adapter.
 */
template <typename Type>
concept UblasDenseVector =
    IsUblasDenseVector<std::remove_cvref_t<Type>>::value;

namespace detail
{
/** @brief Construct a const-correct matrix view with checked independent Eigen axes.
 * @note This adapter has no polynomial-degree or allocation policy of its own.
 */
template<int Rows, int Cols, UblasDenseMatrix Matrix>
[[nodiscard]] auto makeEigenMatrixMap( Matrix& matrix )
{
    static_assert( Rows >= 0 || Rows == Eigen::Dynamic );
    static_assert( Cols >= 0 || Cols == Eigen::Dynamic );
    using owner_type = std::remove_cvref_t<Matrix>;
    using scalar_type = typename owner_type::value_type;
    if ( matrix.size1() > static_cast<std::size_t>( std::numeric_limits<Eigen::Index>::max() ) ||
         matrix.size2() > static_cast<std::size_t>( std::numeric_limits<Eigen::Index>::max() ) )
        throw std::length_error( "legacy matrix dimensions exceed Eigen::Index" );
    if ( ( Rows != Eigen::Dynamic && matrix.size1() != static_cast<std::size_t>( Rows ) ) ||
         ( Cols != Eigen::Dynamic && matrix.size2() != static_cast<std::size_t>( Cols ) ) )
        throw std::invalid_argument( "legacy matrix dimensions disagree with the requested Eigen view" );
    constexpr bool columnMajor = std::same_as<typename owner_type::orientation_category, ublas::column_major_tag>;
    constexpr int layout = Cols == 1 ? Eigen::ColMajor : Rows == 1 ? Eigen::RowMajor :
                           columnMajor ? Eigen::ColMajor : Eigen::RowMajor;
    using plain_type = Eigen::Matrix<scalar_type, Rows, Cols, layout>;
    using view_type = std::conditional_t<std::is_const_v<Matrix>, plain_type const, plain_type>;
    return Eigen::Map<view_type, Eigen::Unaligned>( std::to_address( matrix.data().begin() ),
        static_cast<Eigen::Index>( matrix.size1() ), static_cast<Eigen::Index>( matrix.size2() ) );
}
} // namespace detail

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
 * @throws std::length_error If a matrix dimension exceeds Eigen::Index.
 */
template <typename T, typename Layout, typename Storage>
    requires std::same_as<Layout, ublas::row_major> ||
             std::same_as<Layout, ublas::column_major>
[[nodiscard]] auto
eigenMap( ublas::matrix<T, Layout, Storage>& matrix )
{
    return ::Feel::detail::makeEigenMatrixMap<Eigen::Dynamic, Eigen::Dynamic>( matrix );
}

/**
 * @brief Create a read-only zero-copy Eigen view of a dense uBLAS matrix.
 * @copydetails eigenMap(ublas::matrix<T, Layout, Storage>&)
 */
template <typename T, typename Layout, typename Storage>
    requires std::same_as<Layout, ublas::row_major> ||
             std::same_as<Layout, ublas::column_major>
[[nodiscard]] auto
eigenMap( ublas::matrix<T, Layout, Storage> const& matrix )
{
    return ::Feel::detail::makeEigenMatrixMap<Eigen::Dynamic, Eigen::Dynamic>( matrix );
}

/** @brief Borrow a legacy matrix with checked fixed/mixed Eigen extents. */
template<int Rows, int Cols, typename T, typename Layout, typename Storage>
[[nodiscard]] auto eigenMap( ublas::matrix<T, Layout, Storage>& matrix )
{ return ::Feel::detail::makeEigenMatrixMap<Rows, Cols>( matrix ); }

/** @brief Borrow a const legacy matrix with checked fixed/mixed Eigen extents. */
template<int Rows, int Cols, typename T, typename Layout, typename Storage>
[[nodiscard]] auto eigenMap( ublas::matrix<T, Layout, Storage> const& matrix )
{ return ::Feel::detail::makeEigenMatrixMap<Rows, Cols>( matrix ); }

/** @brief Reject views of temporary legacy matrix owners. */
template<typename T, typename Layout, typename Storage>
auto eigenMap( ublas::matrix<T, Layout, Storage>&& ) = delete;
/** @brief Reject views of const temporary legacy matrix owners. */
template<typename T, typename Layout, typename Storage>
auto eigenMap( ublas::matrix<T, Layout, Storage> const&& ) = delete;
/** @brief Reject typed views of temporary legacy matrix owners. */
template<int Rows, int Cols, typename T, typename Layout, typename Storage>
auto eigenMap( ublas::matrix<T, Layout, Storage>&& ) = delete;
/** @brief Reject typed views of const temporary legacy matrix owners. */
template<int Rows, int Cols, typename T, typename Layout, typename Storage>
auto eigenMap( ublas::matrix<T, Layout, Storage> const&& ) = delete;

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

/** @brief Reject views of temporary legacy vector owners. */
template<typename T, typename Storage>
auto eigenMap( ublas::vector<T, Storage>&& ) = delete;
/** @brief Reject views of const temporary legacy vector owners. */
template<typename T, typename Storage>
auto eigenMap( ublas::vector<T, Storage> const&& ) = delete;

namespace detail
{
/** @brief Consume an Eigen expression without changing its declared extents.
 * @tparam T Scalar used only when the legacy overload needs a temporary owner.
 * @param matrix Expression borrowed for this call.
 * @param operation Callable returning void or an owning result; it must not
 * retain or return a view of the borrowed matrix.
 */
template<typename T, typename Derived, typename Operation>
auto withEigenMatrix( Eigen::MatrixBase<Derived> const& matrix, Operation&& operation )
{ return std::forward<Operation>( operation )( matrix.derived() ); }

/** @brief Borrow dense legacy input, materializing expressions or scalar conversions once.
 * @param matrix Legacy expression consumed before any temporary owner expires.
 * @param operation Callable returning void or an owning result, with no retained view.
 */
template<typename T, typename Derived, typename Operation>
auto withEigenMatrix( ublas::matrix_expression<Derived> const& matrix, Operation&& operation )
{
    if constexpr ( UblasDenseMatrix<Derived> && std::is_same_v<T, typename Derived::value_type> )
        return std::forward<Operation>( operation )( eigenMap( matrix() ) );
    else
    {
        ublas::matrix<T> owned( matrix );
        return std::forward<Operation>( operation )( eigenMap( owned ) );
    }
}
} // namespace detail

} // namespace Feel

#endif // FEELPP_FEELALG_EIGENMAP_HPP
