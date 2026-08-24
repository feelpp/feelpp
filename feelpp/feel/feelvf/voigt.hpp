/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme
       Date: 2026-03-25

  Copyright (C) 2026 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file voigt.hpp
   \author Christophe Prud'homme
   \date 2026-03-25
 */
#ifndef FEELPP_VF_VOIGT_HPP
#define FEELPP_VF_VOIGT_HPP 1

#include <array>
#include <concepts>
#include <numbers>
#include <type_traits>
#include <utility>

#include <feel/feelpoly/traits.hpp>
#include <feel/feelvf/concepts.hpp>
#include <feel/feelvf/matvec.hpp>
#include <feel/feelvf/tensorbasis.hpp>

namespace Feel
{
namespace vf
{

enum class SymmetricTensorNotation
{
    Mandel,
    Voigt
};

namespace detail
{

template <int Dim>
inline constexpr int symmetric_storage_size_v = Dim*( Dim+1 )/2;

template <int Dim, int StorageIndex>
consteval std::array<int, 2>
symmetricPairFromStorageIndex()
{
    static_assert( StorageIndex >= 0 && StorageIndex < symmetric_storage_size_v<Dim>,
                   "Voigt storage index is out of bounds" );

    for ( int i = 0; i < Dim; ++i )
        for ( int j = i; j < Dim; ++j )
            if ( Feel::detail::symmetricIndex( i, j, Dim ) == StorageIndex )
                return { i, j };

    return { -1, -1 };
}

template <int Dim, int StorageIndex>
inline constexpr int symmetric_storage_row_v = symmetricPairFromStorageIndex<Dim, StorageIndex>()[0];

template <int Dim, int StorageIndex>
inline constexpr int symmetric_storage_col_v = symmetricPairFromStorageIndex<Dim, StorageIndex>()[1];

template <typename T>
[[nodiscard]] inline auto
storageExpr( T const& value )
{
    if constexpr ( std::is_arithmetic_v<std::remove_cvref_t<T>> )
        return cst( static_cast<double>( value ) );
    else
        return value;
}

template <int StorageIndex, typename ExprT>
[[nodiscard]] inline auto
vectorComponent( ExprT const& expr )
{
    if constexpr ( expression_rows_v<ExprT> == 1 )
        return component<0, StorageIndex>( expr );
    else
        return component<StorageIndex, 0>( expr );
}

template <int Dim, int StorageIndex, typename ExprT>
[[nodiscard]] inline auto
voigtTerm( ExprT const& expr )
{
    constexpr int i = symmetric_storage_row_v<Dim, StorageIndex>;
    constexpr int j = symmetric_storage_col_v<Dim, StorageIndex>;
    return component<i, j>( expr );
}

template <int Dim, int StorageIndex, typename ExprT>
[[nodiscard]] inline auto
mandelTerm( ExprT const& expr )
{
    constexpr int i = symmetric_storage_row_v<Dim, StorageIndex>;
    constexpr int j = symmetric_storage_col_v<Dim, StorageIndex>;

    if constexpr ( i == j )
        return component<i, j>( expr );
    else
        return cst( std::numbers::sqrt2_v<double> )*component<i, j>( expr );
}

template <int Dim, typename ExprT, std::size_t... StorageIndices>
[[nodiscard]] inline auto
voigtImpl( ExprT const& expr, std::index_sequence<StorageIndices...> )
{
    return vec( voigtTerm<Dim, static_cast<int>( StorageIndices )>( expr )... );
}

template <int Dim, typename ExprT, std::size_t... StorageIndices>
[[nodiscard]] inline auto
mandelImpl( ExprT const& expr, std::index_sequence<StorageIndices...> )
{
    return vec( mandelTerm<Dim, static_cast<int>( StorageIndices )>( expr )... );
}

template <int Dim, int Row, int Col>
consteval int
symmetricStorageIndexFromComponent()
{
    constexpr int i = Row < Col ? Row : Col;
    constexpr int j = Row < Col ? Col : Row;
    return Feel::detail::symmetricIndex( i, j, Dim );
}

template <int Dim, int Row, int Col, bool Mandel, typename ExprT>
[[nodiscard]] inline auto
symmetricMatrixEntry( ExprT const& expr )
{
    constexpr int storageIndex = symmetricStorageIndexFromComponent<Dim, Row, Col>();

    if constexpr ( Row == Col )
        return vectorComponent<storageIndex>( expr );
    else if constexpr ( Mandel )
        return cst( 1.0/std::numbers::sqrt2_v<double> )*vectorComponent<storageIndex>( expr );
    else
        return vectorComponent<storageIndex>( expr );
}

template <int Dim, bool Mandel, typename ExprT, std::size_t... FlatIndices>
[[nodiscard]] inline auto
reconstructSymmetricMatrixImpl( ExprT const& expr, std::index_sequence<FlatIndices...> )
{
    return mat<Dim, Dim>( symmetricMatrixEntry<Dim,
                                               static_cast<int>( FlatIndices/Dim ),
                                               static_cast<int>( FlatIndices%Dim ),
                                               Mandel>( expr )... );
}

template <int Dim, int StorageIndex, typename LeftExprT, typename RightExprT>
[[nodiscard]] inline auto
voigtInnerTerm( LeftExprT const& left, RightExprT const& right )
{
    constexpr int i = symmetric_storage_row_v<Dim, StorageIndex>;
    constexpr int j = symmetric_storage_col_v<Dim, StorageIndex>;

    if constexpr ( i == j )
        return vectorComponent<StorageIndex>( left )*vectorComponent<StorageIndex>( right );
    else
        return cst( 2.0 )*vectorComponent<StorageIndex>( left )*vectorComponent<StorageIndex>( right );
}

template <int Dim, typename ExprT, std::size_t FirstStorageIndex, std::size_t... RemainingStorageIndices>
[[nodiscard]] inline auto
unvoigtImpl( ExprT const& expr, std::index_sequence<FirstStorageIndex, RemainingStorageIndices...> )
{
    return reconstructSymmetricMatrixImpl<Dim, false>( expr, std::make_index_sequence<Dim*Dim>{} );
}

template <int Dim, typename ExprT, std::size_t FirstStorageIndex, std::size_t... RemainingStorageIndices>
[[nodiscard]] inline auto
unmandelImpl( ExprT const& expr, std::index_sequence<FirstStorageIndex, RemainingStorageIndices...> )
{
    return reconstructSymmetricMatrixImpl<Dim, true>( expr, std::make_index_sequence<Dim*Dim>{} );
}

template <int Dim, typename LeftExprT, typename RightExprT, std::size_t FirstStorageIndex, std::size_t... RemainingStorageIndices>
[[nodiscard]] inline auto
voigtInnerImpl( LeftExprT const& left, RightExprT const& right,
                std::index_sequence<FirstStorageIndex, RemainingStorageIndices...> )
{
    return ( voigtInnerTerm<Dim, static_cast<int>( FirstStorageIndex )>( left, right ) + ... +
             voigtInnerTerm<Dim, static_cast<int>( RemainingStorageIndices )>( left, right ) );
}

template <int Dim>
consteval auto
naturalSymmetricPairs()
{
    std::array<std::array<int, 2>, symmetric_storage_size_v<Dim>> pairs{};
    int index = 0;

    for ( int i = 0; i < Dim; ++i )
        pairs[index++] = { i, i };

    for ( int i = 0; i < Dim; ++i )
        for ( int j = i+1; j < Dim; ++j )
            pairs[index++] = { i, j };

    return pairs;
}

template <int Dim, int I, int K>
consteval int
naturalSymmetricIndexFromComponent()
{
    constexpr auto pairs = naturalSymmetricPairs<Dim>();

    for ( int index = 0; index < symmetric_storage_size_v<Dim>; ++index )
        if ( pairs[index][0] == I && pairs[index][1] == K )
            return index;

    return -1;
}

template <int Dim, SymmetricTensorNotation Notation, int StorageIndex, typename TupleT>
[[nodiscard]] inline auto
symmetricStorageVectorTerm( TupleT const& components )
{
    constexpr int i = symmetric_storage_row_v<Dim, StorageIndex>;
    constexpr int j = symmetric_storage_col_v<Dim, StorageIndex>;
    constexpr int naturalIndex = naturalSymmetricIndexFromComponent<Dim, i, j>();
    static_assert( naturalIndex >= 0, "Invalid symmetric component requested" );

    auto value = storageExpr( std::get<naturalIndex>( components ) );

    if constexpr ( Notation == SymmetricTensorNotation::Mandel && i != j )
        return cst( std::numbers::sqrt2_v<double> )*value;
    else
        return value;
}

template <int Dim, SymmetricTensorNotation Notation, typename TupleT, std::size_t... StorageIndices>
[[nodiscard]] inline auto
symmetricStorageVectorImpl( TupleT const& components, std::index_sequence<StorageIndices...> )
{
    return vec( symmetricStorageVectorTerm<Dim, Notation, static_cast<int>( StorageIndices )>( components )... );
}

template <int Dim, SymmetricTensorNotation Notation, int I, int K, std::size_t... StorageIndices>
[[nodiscard]] inline auto
symmetricStorageBasisImpl( std::index_sequence<StorageIndices...> )
{
    constexpr int i = ( I < K ) ? I : K;
    constexpr int j = ( I < K ) ? K : I;
    constexpr int storageIndex = symmetricStorageIndexFromComponent<Dim, i, j>();
    constexpr bool isMandelOffDiagonal = ( Notation == SymmetricTensorNotation::Mandel ) && ( i != j );
    constexpr double value = isMandelOffDiagonal ? std::numbers::sqrt2_v<double> : 1.0;

    return vec( cst( static_cast<int>( StorageIndices ) == storageIndex ? value : 0.0 )... );
}

template <int Dim, SymmetricTensorNotation Notation, int StorageIndex, int I, int K, typename ExprT>
[[nodiscard]] inline auto
symmetricStorageComponentTerm( ExprT const& value )
{
    constexpr int row = symmetric_storage_row_v<Dim, StorageIndex>;
    constexpr int col = symmetric_storage_col_v<Dim, StorageIndex>;
    constexpr int i = ( I < K ) ? I : K;
    constexpr int j = ( I < K ) ? K : I;

    if constexpr ( row == i && col == j )
    {
        auto expr = storageExpr( value );
        if constexpr ( Notation == SymmetricTensorNotation::Mandel && i != j )
            return cst( std::numbers::sqrt2_v<double> )*expr;
        else
            return expr;
    }
    else
        return cst( 0.0 );
}

template <int Dim, SymmetricTensorNotation Notation, int I, int K, typename ExprT, std::size_t... StorageIndices>
[[nodiscard]] inline auto
symmetricStorageComponentImpl( ExprT const& value, std::index_sequence<StorageIndices...> )
{
    return vec( symmetricStorageComponentTerm<Dim,
                                             Notation,
                                             static_cast<int>( StorageIndices ),
                                             I,
                                             K>( value )... );
}

template <int Dim, SymmetricTensorNotation Notation, std::size_t... StorageIndices>
[[nodiscard]] inline auto
dynamicSymmetricStorageBasisImpl( int i, int k, std::index_sequence<StorageIndices...> )
{
    int const row = std::min( i, k );
    int const col = std::max( i, k );
    int const storageIndex = Feel::detail::symmetricIndex( row, col, Dim );
    double const value = ( Notation == SymmetricTensorNotation::Mandel && row != col ) ?
                         std::numbers::sqrt2_v<double> : 1.0;

    return vec( cst( static_cast<int>( StorageIndices ) == storageIndex ? value : 0.0 )... );
}

template <int Dim, SymmetricTensorNotation Notation, typename ExprT, std::size_t... StorageIndices>
[[nodiscard]] inline auto
dynamicSymmetricStorageComponentImpl( int i, int k, ExprT const& value,
                                      std::index_sequence<StorageIndices...> )
{
    int const row = std::min( i, k );
    int const col = std::max( i, k );
    auto expr = storageExpr( value );

    return vec(
        ( ( symmetric_storage_row_v<Dim, static_cast<int>( StorageIndices )> == row &&
            symmetric_storage_col_v<Dim, static_cast<int>( StorageIndices )> == col ) ?
          ( ( Notation == SymmetricTensorNotation::Mandel && row != col ) ?
            cst( std::numbers::sqrt2_v<double> )*expr :
            expr ) :
          cst( 0.0 ) )... );
}

template <typename ScaleExprT, typename StorageExprT, std::size_t... StorageIndices>
[[nodiscard]] inline auto
scaleSymmetricStorageImpl( ScaleExprT const& scale,
                           StorageExprT const& expr,
                           std::index_sequence<StorageIndices...> )
{
    auto scaleExpr = storageExpr( scale );
    return vec( scaleExpr*vectorComponent<static_cast<int>( StorageIndices )>( expr )... );
}

} // namespace detail

template <detail::StaticSquareMatrixExpression ExprT>
[[nodiscard]] inline auto
voigt( ExprT const& expr )
{
    constexpr int dim = detail::expression_rows_v<ExprT>;
    return detail::voigtImpl<dim>( expr, std::make_index_sequence<detail::symmetric_storage_size_v<dim>>{} );
}

template <detail::StaticSquareMatrixExpression ExprT>
[[nodiscard]] inline auto
mandel( ExprT const& expr )
{
    constexpr int dim = detail::expression_rows_v<ExprT>;
    return detail::mandelImpl<dim>( expr, std::make_index_sequence<detail::symmetric_storage_size_v<dim>>{} );
}

template <int Dim, SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel, typename... ComponentExprT>
requires ( sizeof...( ComponentExprT ) == detail::symmetric_storage_size_v<Dim> )
[[nodiscard]] inline auto
symm_storage_vec( ComponentExprT const&... components )
{
    auto tuple = std::forward_as_tuple( components... );
    return detail::symmetricStorageVectorImpl<Dim, Notation>(
        tuple, std::make_index_sequence<detail::symmetric_storage_size_v<Dim>>{} );
}

template <int Dim, typename... ComponentExprT>
requires ( sizeof...( ComponentExprT ) == detail::symmetric_storage_size_v<Dim> )
[[nodiscard]] inline auto
mandel_vec( ComponentExprT const&... components )
{
    return symm_storage_vec<Dim, SymmetricTensorNotation::Mandel>( components... );
}

template <int Dim, typename... ComponentExprT>
requires ( sizeof...( ComponentExprT ) == detail::symmetric_storage_size_v<Dim> )
[[nodiscard]] inline auto
voigt_vec( ComponentExprT const&... components )
{
    return symm_storage_vec<Dim, SymmetricTensorNotation::Voigt>( components... );
}

template <int Dim, SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel, int I, int K>
[[nodiscard]] inline auto
symm_storage_basis()
{
    static_assert( I >= 0 && I < Dim, "symm_storage_basis() row index out of bounds" );
    static_assert( K >= 0 && K < Dim, "symm_storage_basis() column index out of bounds" );

    return detail::symmetricStorageBasisImpl<Dim, Notation, I, K>(
        std::make_index_sequence<detail::symmetric_storage_size_v<Dim>>{} );
}

template <int Dim, SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel>
[[nodiscard]] inline auto
symm_storage_basis( int i, int k )
{
    CHECK( i >= 0 && i < Dim ) << "symm_storage_basis() row index " << i << " out of bounds for dimension " << Dim;
    CHECK( k >= 0 && k < Dim ) << "symm_storage_basis() column index " << k << " out of bounds for dimension " << Dim;

    return detail::dynamicSymmetricStorageBasisImpl<Dim, Notation>(
        i, k, std::make_index_sequence<detail::symmetric_storage_size_v<Dim>>{} );
}

template <int Dim, SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel, int I, int K, typename ExprT>
[[nodiscard]] inline auto
symm_storage_component( ExprT const& value )
{
    static_assert( I >= 0 && I < Dim, "symm_storage_component() row index out of bounds" );
    static_assert( K >= 0 && K < Dim, "symm_storage_component() column index out of bounds" );

    return detail::symmetricStorageComponentImpl<Dim, Notation, I, K>(
        value, std::make_index_sequence<detail::symmetric_storage_size_v<Dim>>{} );
}

template <int Dim, SymmetricTensorNotation Notation = SymmetricTensorNotation::Mandel, typename ExprT>
[[nodiscard]] inline auto
symm_storage_component( int i, int k, ExprT const& value )
{
    CHECK( i >= 0 && i < Dim ) << "symm_storage_component() row index " << i << " out of bounds for dimension " << Dim;
    CHECK( k >= 0 && k < Dim ) << "symm_storage_component() column index " << k << " out of bounds for dimension " << Dim;

    return detail::dynamicSymmetricStorageComponentImpl<Dim, Notation>(
        i, k, value, std::make_index_sequence<detail::symmetric_storage_size_v<Dim>>{} );
}

template <int Dim, int I, int K>
[[nodiscard]] inline auto
mandel_basis()
{
    return symm_storage_basis<Dim, SymmetricTensorNotation::Mandel, I, K>();
}

template <int Dim>
[[nodiscard]] inline auto
mandel_basis( int i, int k )
{
    return symm_storage_basis<Dim, SymmetricTensorNotation::Mandel>( i, k );
}

template <int Dim, int I, int K, typename ExprT>
[[nodiscard]] inline auto
mandel_component( ExprT const& value )
{
    return symm_storage_component<Dim, SymmetricTensorNotation::Mandel, I, K>( value );
}

template <int Dim, typename ExprT>
[[nodiscard]] inline auto
mandel_component( int i, int k, ExprT const& value )
{
    return symm_storage_component<Dim, SymmetricTensorNotation::Mandel>( i, k, value );
}

template <int Dim, int I, int K>
[[nodiscard]] inline auto
voigt_basis()
{
    return symm_storage_basis<Dim, SymmetricTensorNotation::Voigt, I, K>();
}

template <int Dim>
[[nodiscard]] inline auto
voigt_basis( int i, int k )
{
    return symm_storage_basis<Dim, SymmetricTensorNotation::Voigt>( i, k );
}

template <int Dim, int I, int K, typename ExprT>
[[nodiscard]] inline auto
voigt_component( ExprT const& value )
{
    return symm_storage_component<Dim, SymmetricTensorNotation::Voigt, I, K>( value );
}

template <int Dim, typename ExprT>
[[nodiscard]] inline auto
voigt_component( int i, int k, ExprT const& value )
{
    return symm_storage_component<Dim, SymmetricTensorNotation::Voigt>( i, k, value );
}

template <detail::StaticSymmetricVectorExpression StorageExprT, typename ScaleExprT>
[[nodiscard]] inline auto
scale_symm_storage( ScaleExprT const& scale, StorageExprT const& expr )
{
    return detail::scaleSymmetricStorageImpl( scale, expr,
                                              std::make_index_sequence<detail::vector_length_v<StorageExprT>>{} );
}

template <detail::StaticSymmetricVectorExpression ExprT>
[[nodiscard]] inline auto
unvoigt( ExprT const& expr )
{
    constexpr int dim = detail::symmetricDimensionFromStorageSize<detail::vector_length_v<ExprT>>();
    static_assert( dim > 0, "unvoigt() expects a vector whose size matches n(n+1)/2" );
    return detail::unvoigtImpl<dim>( expr, std::make_index_sequence<detail::symmetric_storage_size_v<dim>>{} );
}

template <detail::StaticSymmetricVectorExpression ExprT>
[[nodiscard]] inline auto
unmandel( ExprT const& expr )
{
    constexpr int dim = detail::symmetricDimensionFromStorageSize<detail::vector_length_v<ExprT>>();
    static_assert( dim > 0, "unmandel() expects a vector whose size matches n(n+1)/2" );
    return detail::unmandelImpl<dim>( expr, std::make_index_sequence<detail::symmetric_storage_size_v<dim>>{} );
}

template <detail::StaticSymmetricVectorExpression LeftExprT, detail::StaticSymmetricVectorExpression RightExprT>
[[nodiscard]] inline auto
voigt_inner( LeftExprT const& left, RightExprT const& right )
{
    constexpr int leftDim = detail::symmetricDimensionFromStorageSize<detail::vector_length_v<LeftExprT>>();
    constexpr int rightDim = detail::symmetricDimensionFromStorageSize<detail::vector_length_v<RightExprT>>();

    static_assert( leftDim == rightDim, "voigt_inner() expects vectors with matching symmetric storage size" );

    return detail::voigtInnerImpl<leftDim>( left, right,
                                            std::make_index_sequence<detail::symmetric_storage_size_v<leftDim>>{} );
}

} // namespace vf
} // namespace Feel

#endif /* FEELPP_VF_VOIGT_HPP */
