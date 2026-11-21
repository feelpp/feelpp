/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/** @file componentlayout.hpp
 * @brief Shared scalar identity and component packing at polynomial API boundaries.
 */
#ifndef FEELPP_FEELPOLY_COMPONENTLAYOUT_HPP
#define FEELPP_FEELPOLY_COMPONENTLAYOUT_HPP 1
#include <feel/feelalg/eigenmap.hpp>
#include <feel/feelpoly/eigenextent.hpp>
#include <feel/feelpoly/concepts.hpp>
#include <utility>

namespace Feel::detail
{
/** @brief Borrow lvalue identity conversions and own results from temporary expressions. */
struct IdentityCoefficientLayout
{
    /** @brief Preserve scalar coefficient extents and lvalue ownership. */
    template<typename Matrix>
        requires ( EigenMatrix<Matrix> || UblasDenseMatrix<Matrix> )
    static decltype(auto) toMatrix( Matrix&& matrix )
    {
        if constexpr ( std::is_lvalue_reference_v<Matrix&&> ) return std::forward<Matrix>( matrix );
        else if constexpr ( EigenMatrix<Matrix> )
        {
            using Owner = typename std::remove_cvref_t<Matrix>::PlainObject;
            return Owner( std::forward<Matrix>( matrix ) );
        }
        else return std::remove_cvref_t<Matrix>( std::forward<Matrix>( matrix ) );
    }

    /** @brief Scalar packing is the same identity conversion in either direction. */
    template<typename Matrix>
        requires ( EigenMatrix<Matrix> || UblasDenseMatrix<Matrix> )
    static decltype(auto) toType( Matrix&& matrix )
    { return toMatrix( std::forward<Matrix>( matrix ) ); }
};

/** @brief Pack interleaved component rows into component blocks of columns.
 * @tparam Components Number of stored field components, including full symmetric tensors.
 *
 * Interleaved coefficients C(f Components+c,j) map to A(f,c Np+j).
 * Both directions share this indexing rule for Eigen and legacy callers.
 */
template<int Components>
struct ComponentCoefficientLayout
{
    static_assert( Components > 0 );

    /** @brief Pack Eigen coefficients while preserving independent static extents. */
    template<typename Derived>
    static auto toMatrix( Eigen::MatrixBase<Derived> const& coefficients )
    { return reshape<false>( coefficients ); }

    /** @brief Unpack component blocks into interleaved Eigen coefficient rows. */
    template<typename Derived>
    static auto toType( Eigen::MatrixBase<Derived> const& coefficients )
    { return reshape<true>( coefficients ); }

    /** @brief Adapt a legacy matrix without another owning numerical implementation. */
    template<typename T, typename Layout, typename Storage>
    static ublas::matrix<T> toMatrix( ublas::matrix<T, Layout, Storage> const& coefficients )
    { return reshapeLegacy<false>( coefficients ); }

    /** @brief Materialize a legacy expression once at the compatibility boundary. */
    template<typename Derived>
    static auto toType( ublas::matrix_expression<Derived> const& expression )
    {
        if constexpr ( UblasDenseMatrix<Derived> ) return reshapeLegacy<true>( expression() );
        else
        {
            ublas::matrix<typename Derived::value_type> coefficients( expression );
            return reshapeLegacy<true>( coefficients );
        }
    }

private:
    /** @brief Derive a checked static component axis without erasing known extents. */
    template<int Extent, bool Multiply>
    static consteval int componentExtent()
    {
        if constexpr ( Extent == Dynamic ) return Dynamic;
        else if constexpr ( Multiply ) return staticEigenExtent( checkedProduct( Extent, Components ) );
        else
        {
            static_assert( Extent % Components == 0, "coefficient axis must be divisible by the component count" );
            return Extent / Components;
        }
    }

    /** @brief Validate runtime dimensions before allocation or Eigen block operations. */
    template<bool Unpack, typename Derived>
    static auto shape( Eigen::MatrixBase<Derived> const& input )
    {
        Eigen::Index const dividedAxis = Unpack ? input.cols() : input.rows();
        Eigen::Index const multipliedAxis = Unpack ? input.rows() : input.cols();
        if ( dividedAxis < 0 || dividedAxis % Components != 0 )
            throw std::invalid_argument( "coefficient axis must be divisible by the component count" );
        Eigen::Index const multiplied = matrixElementCount( multipliedAxis, Components, sizeof( typename Derived::Scalar ) );
        Eigen::Index const rows = Unpack ? multiplied : dividedAxis / Components;
        Eigen::Index const cols = Unpack ? dividedAxis / Components : multiplied;
        (void)matrixElementCount( rows, cols, sizeof( typename Derived::Scalar ) );
        return std::pair{ rows, cols };
    }

    /** @brief Apply the sole component indexing rule to an already allocated output. */
    template<bool Unpack, typename Target, typename Derived>
    static void reshapeInto( Target&& target, Eigen::MatrixBase<Derived> const& input )
    {
        Eigen::Index const functions = Unpack ? input.rows() : target.rows();
        Eigen::Index const modes = Unpack ? target.cols() : input.cols();
        for ( Eigen::Index function = 0; function < functions; ++function )
            for ( int component = 0; component < Components; ++component )
            {
                if constexpr ( Unpack )
                    target.row( function * Components + component ) = input.block( function, component * modes, 1, modes );
                else
                    target.block( function, component * modes, 1, modes ) = input.row( function * Components + component );
            }
    }

    /** @brief Allocate a budget-aware Eigen result and fill it through the shared rule. */
    template<bool Unpack, typename Derived>
    static auto reshape( Eigen::MatrixBase<Derived> const& input )
    {
        constexpr int rows = componentExtent<Derived::RowsAtCompileTime, Unpack>();
        constexpr int cols = componentExtent<Derived::ColsAtCompileTime, !Unpack>();
        using Result = typename EigenMatrixBuffer<typename Derived::Scalar, rows, cols>::matrix_type;
        auto const [outputRows, outputCols] = shape<Unpack>( input );
        Result result( outputRows, outputCols );
        reshapeInto<Unpack>( result, input );
        return result;
    }

    /** @brief Allocate only the requested legacy output, borrowing its storage for filling. */
    template<bool Unpack, UblasDenseMatrix Matrix>
    static auto reshapeLegacy( Matrix const& input )
    {
        auto source = eigenMap( input );
        auto const [rows, cols] = shape<Unpack>( source );
        ublas::matrix<typename Matrix::value_type> result( rows, cols );
        reshapeInto<Unpack>( eigenMap( result ), source );
        return result;
    }
};
} // namespace Feel::detail
#endif // FEELPP_FEELPOLY_COMPONENTLAYOUT_HPP
