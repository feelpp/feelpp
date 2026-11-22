/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/** @file polynomialcoefficients.hpp
 * @brief One Eigen coefficient owner and explicit legacy materialization boundaries.
 */
#ifndef FEELPP_FEELPOLY_POLYNOMIALCOEFFICIENTS_HPP
#define FEELPP_FEELPOLY_POLYNOMIALCOEFFICIENTS_HPP 1
#include <feel/feelalg/eigenmap.hpp>
#include <feel/feelpoly/componentlayout.hpp>
#include <type_traits>
#include <utility>

namespace Feel::detail
{

/** @brief Fill disjoint output with C Phi, retaining a contiguous single-sample specialization.
 * @note The caller owns output storage disjoint from both operands. Shape checks
 * precede Eigen evaluation; no borrowed expression survives this call.
 */
template <typename Target, typename Left, typename Right>
void multiplyPolynomialMatricesInto( Eigen::MatrixBase<Target>& output,
                                     Eigen::MatrixBase<Left> const& coefficients, Eigen::MatrixBase<Right> const& values )
{
    // Validate all product dimensions before evaluating into disjoint output storage.
    if ( coefficients.cols() != values.rows() || output.rows() != coefficients.rows() || output.cols() != values.cols() )
        throw std::invalid_argument( "polynomial product dimensions disagree" );

    auto multiply = [&]( auto const& factor )
    { output.derived().noalias() = coefficients * factor; };

    // Expose contiguous single-sample inputs as vectors for Eigen's product kernel.
    if constexpr ( Right::ColsAtCompileTime == Dynamic && ( Right::Flags & Eigen::DirectAccessBit ) )
    {
        if ( values.cols() == 1 && values.derived().rowStride() == 1 )
        {
            using Column = Eigen::Matrix<typename Right::Scalar, Right::RowsAtCompileTime, 1>;
            multiply( Eigen::Map<Column const>( values.derived().data(), values.rows() ) );
            return;
        }
    }

    // General point batches use the same coefficient product without a view refinement.
    multiply( values );
}

/** @brief Consume an Eigen coefficient expression without retaining a borrowed view. */
template<typename T, typename Derived, typename Operation>
void withPolynomialMatrix( Eigen::MatrixBase<Derived> const& coefficients, Operation&& operation )
{ withEigenMatrix<T>( coefficients, std::forward<Operation>( operation ) ); }

/** @brief Consume a dense legacy view, materializing a legacy expression only once.
 * @tparam T Scalar type for an expression that requires temporary ownership.
 * @param coefficients Legacy matrix or expression borrowed for this call.
 * @param operation Callable consuming the Eigen view before temporary storage expires.
 */
template<typename T, typename Derived, typename Operation>
void withPolynomialMatrix( ublas::matrix_expression<Derived> const& coefficients,
                           Operation&& operation )
{
    withEigenMatrix<T>( coefficients, std::forward<Operation>( operation ) );
}

/** @brief Materialize an Eigen expression once in a requested compatibility matrix. */
template<typename Matrix, typename Derived>
[[nodiscard]] Matrix materializePolynomialMatrix( Eigen::MatrixBase<Derived> const& expression )
{
    (void)matrixElementCount( expression.rows(), expression.cols(), sizeof( typename Derived::Scalar ) );
    Matrix result( expression.rows(), expression.cols() );
    if constexpr ( EigenMatrix<Matrix> ) result.noalias() = expression;
    else eigenMap( result ).noalias() = expression;
    return result;
}

/** @brief Derive coefficient columns from the expansion specification, independently of row count. */
template<typename Poly, int OrderSpec>
consteval int polynomialCoefficientExtent()
{
    if constexpr ( OrderSpec == Dynamic ) return Dynamic;
    else
    {
        constexpr auto kind = Poly::convex_type::is_simplex ? NodalCellKind::Simplex : NodalCellKind::Hypercube;
        return PolynomialStorage<typename Poly::value_type, Poly::nDim, OrderSpec, kind>::basisExtent;
    }
}

/** @brief Own coefficients C(component/function, mode) in one checked Eigen buffer.
 * @tparam T Scalar type.
 * @tparam Columns Static expansion cardinality or Dynamic; the row count is independent.
 *
 * Legacy matrices are materialized only when requested. No mirrored cache or
 * borrowed expression survives a call. Copies and moves retain value semantics.
 */
template<typename T, int Columns>
class PolynomialCoefficientStorage
{
public:
    /** Checked owner retains a known expansion axis with freely sized function rows. */
    using Buffer = EigenMatrixBuffer<T, Dynamic, Columns>;
    /** Scalar type used by coefficient expressions. */
    using value_type = T;
    /** Logical mode count available to compile-time consumers. */
    static constexpr int colsAtCompileTime = Columns;

    /** @brief Construct an empty coefficient set with its semantic column extent. */
    PolynomialCoefficientStorage() = default;

    /** @brief Own an Eigen expression without retaining any input view. */
    template<typename Derived>
    explicit PolynomialCoefficientStorage( Eigen::MatrixBase<Derived> const& coefficients )
    { assign( coefficients ); }

    /** @brief Copy legacy coefficients once into the canonical owner. */
    template<typename Derived>
    explicit PolynomialCoefficientStorage( ublas::matrix_expression<Derived> const& coefficients )
    { assign( coefficients ); }

    /** @brief Return the number of interleaved function/component rows. */
    [[nodiscard]] Eigen::Index rows() const noexcept { return M_coefficients.rows(); }
    /** @brief Return the effective number of expansion modes. */
    [[nodiscard]] Eigen::Index cols() const noexcept { return M_coefficients.cols(); }
    /** @brief Borrow mutable coefficients without exposing resizing. */
    [[nodiscard]] auto matrix() & { return M_coefficients.matrix(); }
    /** @brief Borrow immutable coefficients until the owner changes or is destroyed. */
    [[nodiscard]] auto matrix() const& { return M_coefficients.matrix(); }
    /** @brief Prevent borrowed views of temporary owners. */
    auto matrix() && = delete;
    /** @brief Prevent borrowed views of const temporary owners. */
    auto matrix() const&& = delete;

    /** @brief Replace coefficients, evaluating aliases before changing storage. */
    template<typename Derived>
    void assign( Eigen::MatrixBase<Derived> const& coefficients )
    { M_coefficients.assign( coefficients ); }

    /** @brief Borrow a dense legacy input or materialize a legacy expression once. */
    template<typename Derived>
    void assign( ublas::matrix_expression<Derived> const& coefficients )
    {
        withPolynomialMatrix<T>( coefficients, [this]( auto const& values ) { assign( values ); } );
    }

    /** @brief Apply one modal matrix product after checking its mathematical axes. */
    template<typename Derived>
    [[nodiscard]] auto multiply( Eigen::MatrixBase<Derived> const& rhs ) const
    {
        using Matrix = Eigen::Matrix<T, Dynamic, Derived::ColsAtCompileTime>;
        return multiplyAs<Matrix>( rhs );
    }

    /** @brief Write C Phi directly into a fresh result in the caller's representation.
     * @tparam Matrix Owning Eigen or legacy matrix; neither operand can alias it.
     * @throws std::invalid_argument If coefficient and modal dimensions disagree.
     */
    template<typename Matrix, typename Derived>
    [[nodiscard]] Matrix multiplyAs( Eigen::MatrixBase<Derived> const& rhs ) const
    {
        if ( cols() != rhs.rows() )
            throw std::invalid_argument( "polynomial coefficient and expansion dimensions are incompatible" );
        (void)matrixElementCount( rows(), rhs.cols(), sizeof( T ) );
        Matrix result( rows(), rhs.cols() );
        if constexpr ( EigenMatrix<Matrix> )
            multiplyPolynomialMatricesInto( result, matrix(), rhs );
        else
        {
            auto target = eigenMap( result );
            multiplyPolynomialMatricesInto( target, matrix(), rhs );
        }
        return result;
    }

    /** @brief Own C D directly without copying coefficients that the product replaces.
     * @param rhs Right factor, possibly borrowing this coefficient owner.
     * @throws std::invalid_argument If the inner or bound output dimensions disagree.
     * @note Source coefficients are unchanged; the checked owner evaluates every
     * borrowed operand before returning an independent result.
     */
    template<typename Derived>
    [[nodiscard]] PolynomialCoefficientStorage rightProduct( Eigen::MatrixBase<Derived> const& rhs ) const
    {
        if ( cols() != rhs.rows() )
            throw std::invalid_argument( "polynomial coefficient and operator dimensions are incompatible" );
        return PolynomialCoefficientStorage( matrix() * rhs );
    }

    /** @brief Compose functions in packed component space, or assign modal coefficients directly. */
    template<typename Layout, typename Derived>
    void set( Eigen::MatrixBase<Derived> const& coefficients, bool asIs )
    {
        if ( asIs ) { assign( coefficients ); return; }
        auto current = matrix();
        decltype(auto) packed = Layout::toMatrix( current );
        if ( coefficients.cols() != packed.rows() )
            throw std::invalid_argument( "polynomial change-of-basis dimensions are incompatible" );
        auto product = ( coefficients * packed ).eval();
        assign( Layout::toType( product ) );
    }

    /** @brief Adapt legacy coefficients to the same checked composition kernel. */
    template<typename Layout, typename Derived>
    void set( ublas::matrix_expression<Derived> const& coefficients, bool asIs )
    {
        withPolynomialMatrix<T>( coefficients,
            [this, asIs]( auto const& values ) { set<Layout>( values, asIs ); } );
    }

    /** @brief Return an owning legacy result; its lifetime is independent of this owner. */
    template<typename Matrix>
    [[nodiscard]] Matrix legacy() const
    { return materializePolynomialMatrix<Matrix>( matrix() ); }

private:
    /** Sole coefficient storage; legacy callers receive explicit owning copies. */
    Buffer M_coefficients;
};

/** @brief Copy native coefficients directly, adapting an unmigrated source only once. */
template<typename Storage, typename Source>
[[nodiscard]] Storage polynomialCoefficientStorageFrom( Source const& source )
{
    if constexpr ( requires { source.coefficients(); } )
    {
        if constexpr ( EigenMatrix<decltype( source.coefficients() )> )
            return Storage( source.coefficients() );
        else return Storage( source.coeff() );
    }
    else return Storage( source.coeff() );
}

/** @brief Own supplied modal coefficients directly, or compose source-relative coefficients.
 * @tparam Storage Canonical coefficient owner, including its logical expansion extent.
 * @tparam Layout Function/component packing policy.
 * @tparam PackInput Whether relative polynomial components also need packing.
 * @param source Basis and existing functions used by the relative interpretation.
 * @param coefficients Modal coefficients when asIs, otherwise weights on source functions.
 * @param asIs Whether the supplied matrix already uses the source's modal expansion.
 * @note Direct modal construction does not copy source coefficients that it would replace.
 */
template<typename Storage, typename Layout, bool PackInput = false, typename Source, typename Coefficients>
[[nodiscard]] Storage polynomialCoefficientStorageFrom( Source const& source,
                                                        Coefficients const& coefficients, bool asIs )
{
    if ( asIs ) return Storage( coefficients );
    auto result = polynomialCoefficientStorageFrom<Storage>( source );
    withPolynomialMatrix<typename Storage::value_type>( coefficients, [&]( auto const& values )
    {
        if constexpr ( PackInput )
        {
            decltype(auto) packed = Layout::toMatrix( values.derived() );
            result.template set<Layout>( packed, false );
        }
        else result.template set<Layout>( values, false );
    } );
    return result;
}
} // namespace Feel::detail
#endif // FEELPP_FEELPOLY_POLYNOMIALCOEFFICIENTS_HPP
