/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/** @file eigenextent.hpp
 * @brief Checked Eigen storage with independent polynomial, shape and allocation metadata.
 */
#ifndef FEELPP_FEELPOLY_EIGENEXTENT_HPP
#define FEELPP_FEELPOLY_EIGENEXTENT_HPP 1

#include <feel/feelpoly/order.hpp>
#include <Eigen/Core>
#include <array>
#include <concepts>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>

namespace Feel::detail
{
static_assert( Feel::Dynamic == Eigen::Dynamic );

/** Default per-buffer inline ceiling, following Eigen's configured limit.
 * A zero Eigen limit means unlimited; use the existing assembly ceiling then.
 * Callers divide this budget among simultaneously live workspace buffers.
 */
inline constexpr std::size_t defaultInlineMatrixBytes =
    EIGEN_STACK_ALLOCATION_LIMIT > 0 ? EIGEN_STACK_ALLOCATION_LIMIT : 128 * 1024;

/** @brief Check an extent used as an Eigen template argument before narrowing. */
[[nodiscard]] constexpr int staticEigenExtent( std::size_t count )
{
    if ( count > static_cast<std::size_t>( std::numeric_limits<int>::max() ) )
        throw std::length_error( "Eigen compile-time extent exceeds int" );
    return static_cast<int>( count );
}

/** @brief Validate signed dimensions, coefficient counts and allocation bytes.
 * @param rows Requested row count.
 * @param cols Requested column count.
 * @param scalarBytes Scalar size in bytes.
 * @return Checked coefficient count.
 * @throws std::invalid_argument For negative dimensions.
 * @throws std::length_error For index or allocation-size overflow.
 */
[[nodiscard]] constexpr Eigen::Index matrixElementCount(
    Eigen::Index rows, Eigen::Index cols, std::size_t scalarBytes )
{
    if ( rows < 0 || cols < 0 )
        throw std::invalid_argument( "matrix dimensions must be non-negative" );
    if ( cols != 0 && rows > std::numeric_limits<Eigen::Index>::max() / cols )
        throw std::length_error( "matrix coefficient count exceeds Eigen::Index" );
    auto const count = rows * cols;
    auto const maxBytes = static_cast<std::size_t>(
        std::numeric_limits<std::ptrdiff_t>::max() );
    if ( scalarBytes == 0 || static_cast<std::size_t>( count ) > maxBytes / scalarBytes )
        throw std::length_error( "matrix allocation byte count overflows" );
    return count;
}

/** @brief Empty static axis contract or one bound runtime extent.
 * @tparam Extent Compile-time extent or `Dynamic`.
 */
template<int Extent>
class MatrixAxis
{
public:
    /** @brief Validate an optional runtime binding against a static axis. */
    explicit constexpr MatrixAxis( Eigen::Index bound = Dynamic )
    {
        if ( bound != Dynamic && bound != Extent )
            throw std::invalid_argument( "axis binding disagrees with its static extent" );
    }
    /** @brief Whether a requested size respects this static axis. */
    [[nodiscard]] constexpr bool accepts( Eigen::Index count ) const noexcept { return count == Extent; }
    /** @brief Static axes need no runtime metadata on assignment. */
    constexpr void adopt( MatrixAxis const& ) noexcept {}
};

/** @brief Runtime axis that is either free or bound to one mathematical extent. */
template<>
class MatrixAxis<Dynamic>
{
public:
    /** @brief Bind a runtime axis, or retain a free axis with `Dynamic`. */
    explicit constexpr MatrixAxis( Eigen::Index bound = Dynamic ) : M_bound( bound )
    {
        if ( bound < 0 && bound != Dynamic )
            throw std::invalid_argument( "invalid runtime axis binding" );
    }
    /** @brief Whether a requested size respects a bound runtime axis. */
    [[nodiscard]] constexpr bool accepts( Eigen::Index count ) const noexcept
    { return M_bound == Dynamic || count == M_bound; }
    /** @brief An unbound destination adopts the source's mathematical contract. */
    constexpr void adopt( MatrixAxis const& source ) noexcept
    { if ( M_bound == Dynamic ) M_bound = source.M_bound; }
private:
    /** Bound extent or `Dynamic` for an independent, freely sized axis. */
    Eigen::Index M_bound;
};

/** @brief Shared logical shape, byte budget and Eigen allocation representation.
 * Polynomial owners and assembly matrices use this one allocation decision.
 * Mixed static/runtime axes retain the known side. Oversized fully fixed
 * matrices use heap storage while the policy retains their logical extents.
 */
template <typename T, int Rows, int Cols,
          std::size_t InlineBytes = defaultInlineMatrixBytes, int Layout = Eigen::RowMajor>
struct MatrixStoragePolicy
{
    static_assert( Rows >= 0 || Rows == Dynamic );
    static_assert( Cols >= 0 || Cols == Dynamic );
    /** Logical row extent before any allocation fallback. */
    static constexpr int rowsAtCompileTime = Rows;
    /** Logical column extent before any allocation fallback. */
    static constexpr int colsAtCompileTime = Cols;
    /** Vector-compatible storage order, independent of the allocation choice. */
    static constexpr int legalLayout = Cols == 1 ? Eigen::ColMajor :
                                       Rows == 1 ? Eigen::RowMajor : Layout;
    /** Whether the complete fixed shape fits the configured coefficient budget. */
    static constexpr bool fixedFits = []
    {
        if constexpr ( Rows == Dynamic || Cols == Dynamic ) return false;
        else
        {
            constexpr auto count = static_cast<std::size_t>( matrixElementCount( Rows, Cols, sizeof( T ) ) );
            constexpr std::size_t eigenLimit = EIGEN_STACK_ALLOCATION_LIMIT > 0 ?
                EIGEN_STACK_ALLOCATION_LIMIT : InlineBytes;
            constexpr auto ceiling = InlineBytes < eigenLimit ? InlineBytes : eigenLimit;
            if constexpr ( count > ceiling / sizeof( T ) ) return false;
            else return sizeof( Eigen::Matrix<T, Rows, Cols, legalLayout> ) <= InlineBytes;
        }
    }();
    /** Fully fixed shapes exceeding the budget use dynamic owning storage. */
    static constexpr bool heapFallback = Rows != Dynamic && Cols != Dynamic && !fixedFits;

    /** Whether coefficients reside inline in the owner. */
    static constexpr bool usesInlineStorage = fixedFits;
    /** Legal Eigen layout, including P0 and segment vector-shaped matrices. */
    static constexpr int storageOrder = legalLayout;
    /** Backing Eigen type; logical extents remain available on this owner. */
    using matrix_type = Eigen::Matrix<T,
        heapFallback && Cols == 1 ? Dynamic : Rows,
        heapFallback && Cols != 1 ? Dynamic : Cols, storageOrder>;
};

/** @brief Owning matrix whose logical extents survive a heap-allocation fallback.
 * @tparam T Scalar type.
 * @tparam Rows Static logical rows or `Dynamic`.
 * @tparam Cols Static logical columns or `Dynamic`.
 * @tparam InlineBytes Maximum inline coefficient bytes.
 * @tparam Layout Preferred layout, adjusted for Eigen row/column-vector rules.
 *
 * The owner validates all construction and assignment in Release. Borrowed
 * `matrix()` views support Eigen expressions without exposing an owning matrix
 * that could be resized behind the logical shape contract. Views expire when
 * the owner is resized, assigned, moved or destroyed. Rvalue views are forbidden.
 * A moved-from object supports destruction and assignment; its data are unspecified.
 */
template <typename T, int Rows, int Cols,
          std::size_t InlineBytes = defaultInlineMatrixBytes, int Layout = Eigen::RowMajor>
class EigenMatrixBuffer
{
    /** Select checked allocation for a kernel that overwrites every entry. */
    struct OverwriteStorage
    {
    };
    static_assert( Rows >= 0 || Rows == Dynamic );
    static_assert( Cols >= 0 || Cols == Dynamic );
    using StoragePolicy = MatrixStoragePolicy<T, Rows, Cols, InlineBytes, Layout>;

  public:
    /** Scalar type of coefficients and views. */
    using value_type = T;
    /** Logical row extent, independent of allocation. */
    static constexpr int rowsAtCompileTime = Rows;
    /** Logical column extent, independent of allocation. */
    static constexpr int colsAtCompileTime = Cols;
    /** Whether coefficients reside inline in the owner. */
    static constexpr bool usesInlineStorage = StoragePolicy::usesInlineStorage;
    /** Legal Eigen storage order, including row/column vectors. */
    static constexpr int storageOrder = StoragePolicy::storageOrder;
    /** Backing representation selected by the common byte budget. */
    using matrix_type = typename StoragePolicy::matrix_type;

    /** @brief Zero-initialize fixed axes, using zero for unspecified axes. */
    EigenMatrixBuffer() : EigenMatrixBuffer( Rows == Dynamic ? 0 : Rows,
                                             Cols == Dynamic ? 0 : Cols ) {}
    /** @brief Allocate and zero-initialize a validated logical shape. */
    EigenMatrixBuffer( Eigen::Index rows, Eigen::Index cols,
                       Eigen::Index boundRows = Dynamic, Eigen::Index boundCols = Dynamic )
        : EigenMatrixBuffer( OverwriteStorage{}, rows, cols, boundRows, boundCols ) { M_matrix.setZero(); }

    /** @brief Allocate checked storage for an internal full-write kernel.
     * @param rows Requested row count.
     * @param cols Requested sample/column count.
     * @param boundRows Optional semantic row binding.
     * @param boundCols Optional semantic column binding.
     * @note Every entry must be written before observation. Ordinary constructors
     * and polynomial storage factories retain their zero-initialization contract.
     */
    [[nodiscard]] static EigenMatrixBuffer forOverwrite( Eigen::Index rows, Eigen::Index cols,
        Eigen::Index boundRows = Dynamic, Eigen::Index boundCols = Dynamic )
    { return { OverwriteStorage{}, rows, cols, boundRows, boundCols }; }
    /** @brief Copy owned coefficients and shape. */
    EigenMatrixBuffer( EigenMatrixBuffer const& ) = default;
    /** @brief Transfer owned coefficients and shape. */
    EigenMatrixBuffer( EigenMatrixBuffer&& ) = default;
    /** @brief Copy owned coefficients and shape. */
    EigenMatrixBuffer& operator=( EigenMatrixBuffer const& source )
    {
        if ( this == &source ) return *this;
        validateShape( source.rows(), source.cols() );
        M_matrix = source.M_matrix;
        M_rows.adopt( source.M_rows );
        M_cols.adopt( source.M_cols );
        return *this;
    }
    /** @brief Transfer owned coefficients and shape. */
    EigenMatrixBuffer& operator=( EigenMatrixBuffer&& source )
    {
        if ( this == &source ) return *this;
        validateShape( source.rows(), source.cols() );
        M_matrix = std::move( source.M_matrix );
        M_rows.adopt( source.M_rows );
        M_cols.adopt( source.M_cols );
        return *this;
    }

    /** @brief Actual row count of the owned matrix. */
    [[nodiscard]] Eigen::Index rows() const noexcept { return M_matrix.rows(); }
    /** @brief Actual column count of the owned matrix. */
    [[nodiscard]] Eigen::Index cols() const noexcept { return M_matrix.cols(); }
    /** @brief Number of owned coefficients. */
    [[nodiscard]] Eigen::Index size() const noexcept { return M_matrix.size(); }

    /** @brief Borrow read-only backing storage for compatibility APIs requiring an Eigen owner reference.
     * @note Logical extents remain described by this buffer; no mutable or resizable
     * owner is exposed. The reference expires on mutation, movement or destruction.
     */
    [[nodiscard]] matrix_type const& matrixStorage() const& noexcept { return M_matrix; }

    /** Prevent a storage reference into a temporary owner. */
    matrix_type const&& matrixStorage() const&& = delete;

    /** @brief Resize free axes and zero-initialize every coefficient.
     * @throws std::invalid_argument For a negative or mismatched fixed extent.
     * @throws std::length_error For coefficient or allocation-size overflow.
     */
    void resize( Eigen::Index rows, Eigen::Index cols )
    {
        resizeForOverwrite( rows, cols );
        M_matrix.setZero();
    }

    /** @brief Prepare checked storage for a kernel that writes every coefficient.
     * @note Coefficients are unspecified until the caller completes its full write.
     * Point aliases must be snapshotted before this operation. Public factories
     * and resize() retain their zero-initialization contract. An unchanged shape
     * reuses the owner's established invariant; borrowed views cannot resize it.
     */
    void resizeForOverwrite( Eigen::Index rows, Eigen::Index cols )
    {
        if ( rows == M_matrix.rows() && cols == M_matrix.cols() ) return;
        prepareStorage( rows, cols );
    }

    /** @brief Copy an Eigen expression after validating shape, including self-aliasing.
     * @param expression Source, possibly a strided block or a view of this owner.
     * @note Evaluation completes before replacing the old storage.
     */
    template<typename Derived>
    void assign( Eigen::MatrixBase<Derived> const& expression )
    {
        validateShape( expression.rows(), expression.cols() );
        if constexpr ( ( Rows != Dynamic && Derived::RowsAtCompileTime != Dynamic && Rows != Derived::RowsAtCompileTime ) ||
                       ( Cols != Dynamic && Derived::ColsAtCompileTime != Dynamic && Cols != Derived::ColsAtCompileTime ) )
            throw std::invalid_argument( "Eigen expression disagrees with the logical matrix extents" );
        else
        {
            matrix_type evaluated( expression.rows(), expression.cols() );
            evaluated = expression;
            M_matrix.swap( evaluated );
        }
    }

    /** @brief Borrow a mutable Eigen view; resizing the view cannot resize the owner. */
    [[nodiscard]] auto matrix() & noexcept
    {
        return Eigen::Map<matrix_type, Eigen::Unaligned>(
            M_matrix.data(), M_matrix.rows(), M_matrix.cols() );
    }
    /** @brief Borrow an immutable Eigen view with the same lifetime contract. */
    [[nodiscard]] auto matrix() const& noexcept
    {
        return Eigen::Map<matrix_type const, Eigen::Unaligned>(
            M_matrix.data(), M_matrix.rows(), M_matrix.cols() );
    }
    /** @brief Prevent a borrowed view of a temporary owner. */
    auto matrix() && = delete;
    /** @brief Prevent a borrowed view of a const temporary owner. */
    auto matrix() const&& = delete;

private:
    /** @brief Prepare the same checked shape without a redundant initial clear. */
    EigenMatrixBuffer( OverwriteStorage, Eigen::Index rows, Eigen::Index cols,
                       Eigen::Index boundRows, Eigen::Index boundCols )
        : M_rows( boundRows ), M_cols( boundCols ) { prepareStorage( rows, cols ); }

    /** @brief Validate a new logical shape before changing the owned allocation.
     * @note Construction always takes this path, even when the backing Eigen
     * matrix already has the requested shape. Later resizes may reuse a shape
     * only after construction has established all logical bindings.
     */
    void prepareStorage( Eigen::Index rows, Eigen::Index cols )
    {
        validateShape( rows, cols );
        M_matrix.resize( rows, cols );
    }

    /** @brief Validate shape without relying on Eigen debug assertions. */
    void validateShape( Eigen::Index rows, Eigen::Index cols ) const
    {
        (void)matrixElementCount( rows, cols, sizeof( T ) );
        if ( !M_rows.accepts( rows ) || !M_cols.accepts( cols ) )
            throw std::invalid_argument( "matrix shape disagrees with its bound extents" );
    }
    /** Static axes occupy no runtime state; dynamic basis axes retain their binding. */
    [[no_unique_address]] MatrixAxis<Rows> M_rows;
    /** Column constraint is independent of the row constraint. */
    [[no_unique_address]] MatrixAxis<Cols> M_cols;
    /** Unique owning representation; borrowed Eigen views are never retained. */
    matrix_type M_matrix;
};

/** @brief Checked bridge owner exposing mutable, non-owning Eigen matrix views. */
template<typename Buffer>
concept EigenBuffer = requires( Buffer& buffer )
{
    typename Buffer::value_type;
    typename Buffer::matrix_type;
    { Buffer::rowsAtCompileTime } -> std::convertible_to<int>;
    { Buffer::colsAtCompileTime } -> std::convertible_to<int>;
    { buffer.matrix() } -> std::same_as<Eigen::Map<typename Buffer::matrix_type, Eigen::Unaligned>>;
};

/** @brief Borrow a fixed output collection without copying its coefficient owners.
 * @note Pointers remain valid only while the caller's collection stays alive and unmoved.
 */
template<typename Buffer, std::size_t Count>
[[nodiscard]] std::array<Buffer*, Count> matrixBufferPointers( std::array<Buffer, Count>& buffers )
{
    std::array<Buffer*, Count> result;
    for ( std::size_t i = 0; i < Count; ++i ) result[i] = &buffers[i];
    return result;
}

/** @brief Construct a fixed collection directly, without default heap buffers or rebinding.
 * @tparam Buffer Checked matrix owner.
 * @tparam Count Compile-time direction or workspace count, including zero.
 */
template<typename Buffer, std::size_t Count>
[[nodiscard]] auto makeMatrixBuffers( Eigen::Index rows, Eigen::Index cols,
                                      Eigen::Index boundRows = Dynamic, Eigen::Index boundCols = Dynamic )
{
    return [&]<std::size_t... I>( std::index_sequence<I...> )
    { return std::array<Buffer, Count>{ ( static_cast<void>( I ), Buffer{ rows, cols, boundRows, boundCols } )... }; }
        ( std::make_index_sequence<Count>{} );
}

/** @brief Order/cardinality policy for scalar polynomial storage.
 * @tparam T Scalar type.
 * @tparam Dim Fixed reference dimension, zero to three (zero denotes a point).
 * @tparam OrderSpec Static expansion degree or `Dynamic`.
 * @tparam CellKind Simplex total degree or hypercube tensor-product degree.
 * @tparam InlineBytes Aggregate budget for each direction collection.
 *
 * Polynomial rows, sample columns and allocation are independent. A runtime
 * degree equal to one keeps dynamic basis extents. Static degree remains in
 * the type even when an operator or nodal table uses heap-backed storage.
 */
template <typename T, int Dim, int OrderSpec, NodalCellKind CellKind,
          std::size_t InlineBytes = defaultInlineMatrixBytes,
          typename DegreePolicy = OrderBase<OrderSpec>>
class PolynomialStorage : public DegreePolicy
{
    static_assert( Dim >= 0 && Dim <= 3 );
    static_assert( OrderSpec == Dynamic || OrderSpec <= std::numeric_limits<uint16_type>::max() );
    using super = DegreePolicy;

  public:
    /** Unchanged expansion-degree specification. */
    static constexpr int order_spec = OrderSpec;
    /** Compile-time reference dimension. */
    static constexpr int dimension = Dim;
    /** Cell identity retained independently of matrix dimensions. */
    static constexpr NodalCellKind cellKind = CellKind;
    /** Checked static basis count, or the Eigen dynamic extent. */
    static constexpr int basisExtent = []
    {
        if constexpr ( OrderSpec == Dynamic ) return Dynamic;
        else if constexpr ( requires { DegreePolicy::cardinality; } )
            return staticEigenExtent( DegreePolicy::cardinality );
        else return staticEigenExtent( NodalCardinality<CellKind>::total( Dim, OrderSpec ) );
    }();
    /** Basis rows by independently sized sample columns. */
    template<int Samples = Dynamic>
    using ValuesFor = EigenMatrixBuffer<T, basisExtent, Samples, InlineBytes>;
    /** Modal values with a runtime sample count. */
    using Values = ValuesFor<>;
    /** Independently sized polynomial rows by modal coefficient columns. */
    template<int Polynomials = Dynamic>
    using CoefficientsFor = EigenMatrixBuffer<T, Polynomials, basisExtent, InlineBytes>;
    /** Modal coefficients with a runtime function count. */
    using Coefficients = CoefficientsFor<>;
    /** Fixed spatial rows by independently sized point columns. */
    template<int Samples = Dynamic>
    using PointsFor = EigenMatrixBuffer<T, Dim, Samples, InlineBytes>;
    /** Spatial coordinates with a runtime sample count. */
    using Points = PointsFor<>;
    /** Canonical nodal coordinates for point policies with exactly Np nodes. */
    using Nodes = EigenMatrixBuffer<T, Dim, basisExtent, InlineBytes>;
    /** One modal coefficient/value column. */
    using Vector = EigenMatrixBuffer<T, basisExtent, 1, InlineBytes>;
    /** A genuinely square basis operator. */
    using Operator = EigenMatrixBuffer<T, basisExtent, basisExtent, InlineBytes>;
    /** Gradient workspace, sharing its aggregate inline budget among directions. */
    template<int Samples = Dynamic>
    using GradientFor = std::array<EigenMatrixBuffer<T, basisExtent, Samples, ( Dim > 0 ? InlineBytes / Dim : 0 )>, Dim>;
    /** Gradient workspace with runtime sample count. */
    using Gradient = GradientFor<>;
    /** Hessian workspace indexed by i*Dim+j, with a bounded aggregate inline budget. */
    template<int Samples = Dynamic>
    using HessianFor = std::array<EigenMatrixBuffer<T, basisExtent, Samples, ( Dim > 0 ? InlineBytes / ( Dim * Dim ) : 0 )>, Dim * Dim>;
    /** Hessian workspace with runtime sample count. */
    using Hessian = HessianFor<>;

    /** @brief Construct a static policy or a runtime policy at its effective degree. */
    explicit constexpr PolynomialStorage( RuntimeOrder order =
        RuntimeOrder{ OrderSpec == Dynamic ? 0 : OrderSpec } ) : super( order ) {}

    /** @brief Preserve a richer degree configuration supplied by a tensor policy. */
    template <typename Degree>
        requires std::constructible_from<DegreePolicy, Degree>
    explicit constexpr PolynomialStorage( Degree degree )
        : super( degree )
    {
    }

    /** @brief Effective cardinality, constant-evaluated for a static order specification.
     * @note Runtime orders retain checked cardinality arithmetic without FE local-index limits.
     */
    [[nodiscard]] constexpr Eigen::Index basisCount() const
    {
        if constexpr ( OrderSpec != Dynamic )
        {
            static_assert( basisExtent <= std::numeric_limits<Eigen::Index>::max() );
            return basisExtent;
        }
        else
        {
            auto const count = [&]
            {
                if constexpr ( requires( DegreePolicy const& policy ) { policy.basisCount(); } )
                    return super::basisCount();
                else
                    return NodalCardinality<CellKind>::total( Dim, this->order() );
            }();
            if ( count > static_cast<std::size_t>( std::numeric_limits<Eigen::Index>::max() ) )
                throw std::length_error( "basis cardinality exceeds Eigen::Index" );
            return static_cast<Eigen::Index>( count );
        }
    }
    /** @brief Allocate initialized modal values, Np by sampleCount. */
    template<int Samples = Dynamic>
    [[nodiscard]] ValuesFor<Samples> makeValues( Eigen::Index sampleCount ) const
    {
        auto const count = basisCount();
        return { count, sampleCount, count, Dynamic };
    }
    /** @brief Allocate initialized coefficients, polynomialCount by Np. */
    template<int Polynomials = Dynamic>
    [[nodiscard]] CoefficientsFor<Polynomials> makeCoefficients( Eigen::Index polynomialCount ) const
    {
        auto const count = basisCount();
        return { polynomialCount, count, Dynamic, count };
    }
    /** @brief Allocate initialized spatial coordinates with independent point count. */
    template<int Samples = Dynamic>
    [[nodiscard]] PointsFor<Samples> makePoints( Eigen::Index pointCount ) const { return { Dim, pointCount }; }
    /** @brief Allocate canonical nodal coordinates for an Np-node policy. */
    [[nodiscard]] Nodes makeNodes() const
    {
        auto const count = basisCount();
        return { Dim, count, Dim, count };
    }
    /** @brief Allocate one initialized modal vector. */
    [[nodiscard]] Vector makeVector() const
    {
        auto const count = basisCount();
        return { count, 1, count, 1 };
    }
    /** @brief Allocate one initialized square modal operator. */
    [[nodiscard]] Operator makeOperator() const
    {
        auto const count = basisCount();
        return { count, count, count, count };
    }
    /** @brief Allocate gradient workspaces with a fixed direction count. */
    template<int Samples = Dynamic>
    [[nodiscard]] GradientFor<Samples> makeGradient( Eigen::Index sampleCount ) const
    { return makeDirections<GradientFor<Samples>>( sampleCount ); }
    /** @brief Allocate Hessian workspaces with a fixed direction-pair count. */
    template<int Samples = Dynamic>
    [[nodiscard]] HessianFor<Samples> makeHessian( Eigen::Index sampleCount ) const
    { return makeDirections<HessianFor<Samples>>( sampleCount ); }

    /** @brief Import a checked modal value expression without erasing static order. */
    template<typename Derived>
    [[nodiscard]] auto valuesFrom( Eigen::MatrixBase<Derived> const& expression ) const
    {
        if ( expression.rows() != basisCount() )
            throw std::invalid_argument( "modal value rows disagree with the effective degree" );
        auto result = makeValues<Derived::ColsAtCompileTime>( expression.cols() );
        result.assign( expression );
        return result;
    }
    /** @brief Convert typed modal values, preserving cell, dimension and effective degree.
     * @note Normalization/basis conventions are supplied by the owning basis layer.
     */
    template <int OtherSpec, std::size_t OtherBytes, typename OtherPolicy, EigenBuffer Buffer>
    [[nodiscard]] auto convertValues(
        PolynomialStorage<T, Dim, OtherSpec, CellKind, OtherBytes, OtherPolicy> const& source,
        Buffer const& values ) const
    {
        if ( source.order() != this->order() )
            throw std::invalid_argument( "modal conversion requires equal expansion degree" );
        if constexpr ( requires { source.degreeConfiguration(); this->degreeConfiguration(); } )
            if ( source.degreeConfiguration().template degrees<Dim>() != this->degreeConfiguration().template degrees<Dim>() )
                throw std::invalid_argument( "modal conversion requires equal coordinate degrees" );
        auto result = makeValues<Buffer::colsAtCompileTime>( values.cols() );
        result.assign( values.matrix() );
        return result;
    }

private:
    /** @brief Initialize a fixed collection using one shared shape rule. */
    template<typename Collection>
    [[nodiscard]] Collection makeDirections( Eigen::Index sampleCount ) const
    {
        auto const count = basisCount();
        return makeMatrixBuffers<typename Collection::value_type, std::tuple_size_v<Collection>>(
            count, sampleCount, count, Dynamic );
    }
};
} // namespace Feel::detail
#endif // FEELPP_FEELPOLY_EIGENEXTENT_HPP
