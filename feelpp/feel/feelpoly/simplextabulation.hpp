/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/** @file simplextabulation.hpp
 * @brief Singularity-free Dubiner tabulation in regular simplex coordinates.
 */
#ifndef FEELPP_FEELPOLY_SIMPLEXTABULATION_HPP
#define FEELPP_FEELPOLY_SIMPLEXTABULATION_HPP 1
#include <cmath>
#include <cstdint>
#include <feel/feelpoly/jacobirecurrence.hpp>
#include <feel/feelpoly/polynomialjet.hpp>
#include <feel/feelpoly/tabulation.hpp>
#include <utility>

namespace Feel::detail
{
/** @brief Legacy total-degree mode index for already validated multi-indices.
 * @note Triangle j varies fastest; tetrahedron k is outer and j is inner.
 * The 16-bit degree contract makes these 64-bit intermediates representable.
 */
template<int Dim>
[[nodiscard]] constexpr Eigen::Index simplexModeIndex( std::array<int, Dim> const& degrees )
{
    std::uint64_t total = 0;
    for ( int degree : degrees ) total += degree;
    if constexpr ( Dim == 1 ) return static_cast<Eigen::Index>( total );
    if constexpr ( Dim == 2 ) return static_cast<Eigen::Index>( total * ( total + 1 ) / 2 + degrees[1] );
    if constexpr ( Dim == 3 )
    {
        auto const k = static_cast<std::uint64_t>( degrees[2] );
        return static_cast<Eigen::Index>( total * ( total + 1 ) * ( total + 2 ) / 6 +
            k * ( 2 * total + 3 - k ) / 2 + degrees[1] );
    }
}
} // namespace Feel::detail

namespace Feel
{
/** @brief Dubiner basis on the biunit simplex, using homogeneous Jacobi recurrences.
 * @tparam Dim Reference dimension, one through three.
 * @tparam OrderSpec Static maximum total degree or Dynamic.
 * @tparam T Scalar type.
 * @tparam Normalized Whether to use the existing L2-orthonormal convention.
 *
 * H_n(u,t)=t^n P_n(u/t) is evaluated as a polynomial recurrence, never by
 * division. Each coordinate extension multiplies the preceding polynomial
 * factor by H_n; compact jets differentiate that same recurrence. Vertices,
 * edges and faces therefore use exactly the interior algorithm.
 * @see https://doi.org/10.1145/1644001.1644006
 * @see https://dlmf.nist.gov/18.9.E1
 */
template<int Dim, int OrderSpec, typename T = double, bool Normalized = true>
class SimplexDubinerBasis : public ::Feel::detail::PolynomialStorage<T, Dim, OrderSpec, ::Feel::detail::NodalCellKind::Simplex>
{
    /** Shared checked degree, cardinality and independent-axis storage policy. */
    using super = ::Feel::detail::PolynomialStorage<T, Dim, OrderSpec, ::Feel::detail::NodalCellKind::Simplex>;
public:
    /** Explicit modal multi-index; its entries sum to total polynomial degree. */
    using MultiIndex = std::array<int, Dim>;
    /** Modal rows by independently selected sample columns. */
    template<int Samples = Dynamic> using ValuesFor = typename super::template ValuesFor<Samples>;
    /** Values with a runtime sample axis. */
    using Values = ValuesFor<>;
    /** Fixed first-derivative direction collection. */
    template<int Samples = Dynamic> using GradientFor = typename super::template GradientFor<Samples>;
    /** Gradients with a runtime sample axis. */
    using Gradient = GradientFor<>;
    /** Ordered Hessian direction pairs, indexed by i*Dim+j. */
    template<int Samples = Dynamic> using HessianFor = typename super::template HessianFor<Samples>;
    /** Hessians with a runtime sample axis. */
    using Hessian = HessianFor<>;

    /** @brief Retain one point's jets and an alias-safe coordinate batch.
     * @tparam DerivativeOrder Zero for values, one for gradients, two for Hessians.
     * @tparam Samples Static sample count or Dynamic, independent of expansion degree.
     * @note Scratch is O(Np) in derivative records, independent of batch size.
     */
    template<unsigned DerivativeOrder, int Samples = Dynamic>
    class Workspace : public OrderBase<OrderSpec>
    {
        friend class SimplexDubinerBasis;
        /** Algebra selected at compile time; unused derivative slots do not exist. */
        using Jet = ::Feel::detail::PolynomialJet<Dim, DerivativeOrder, T>;
        /** One compact derivative record per mode, with budget-aware inline storage. */
        using Jets = ::Feel::detail::EigenMatrixBuffer<T, Jet::componentCount + 1, super::basisExtent,
            ::Feel::detail::defaultInlineMatrixBytes / 2, Eigen::ColMajor>;
        /** Snapshot protects all coordinate aliases before any output is changed. */
        using Coordinates = ::Feel::detail::EigenMatrixBuffer<T, Dim, Samples, ::Feel::detail::defaultInlineMatrixBytes / 2>;
    public:
        /** @brief Bind scratch to the effective degree and selected sample-axis contract. */
        explicit Workspace( RuntimeOrder degree = RuntimeOrder{ OrderSpec == Dynamic ? 0 : OrderSpec } )
            : OrderBase<OrderSpec>( degree ), M_jets( Jet::componentCount + 1,
                super{ degree }.basisCount(), Jet::componentCount + 1, super{ degree }.basisCount() ),
              M_coordinates( Dim, Samples == Dynamic ? 0 : Samples )
        {
            SimplexDubinerBasis{ degree }.forEachMode( [&]( Eigen::Index row, MultiIndex const& index )
            { M_jets.matrix()( Jet::componentCount, row ) = modeNormalization( index ); } );
        }
    private:
        /** @brief Read one already initialized modal derivative record. */
        [[nodiscard]] Jet jet( Eigen::Index row ) const
        {
            Jet result;
            auto const source = M_jets.matrix();
            for ( int c = 0; c < Jet::componentCount; ++c ) result.data[c] = source( c, row );
            return result;
        }
        /** @brief Write every scalar of one modal derivative record. */
        void setJet( Eigen::Index row, Jet const& value )
        {
            auto target = M_jets.matrix();
            for ( int c = 0; c < Jet::componentCount; ++c ) target( c, row ) = value.data[c];
        }
        /** Owned scratch plus one invariant normalization scalar per mode. */
        Jets M_jets;
        /** Owned coordinates; callers cannot alias or resize this storage. */
        Coordinates M_coordinates;
    };

    /** @brief Construct a checked static or runtime total-degree basis. */
    explicit SimplexDubinerBasis( RuntimeOrder degree = RuntimeOrder{ OrderSpec == Dynamic ? 0 : OrderSpec } )
        : super( degree ) {}

    /** @brief Map a valid modal multi-index into the preserved legacy ordering. */
    [[nodiscard]] Eigen::Index modeIndex( MultiIndex const& degrees ) const
    {
        std::int64_t total = 0;
        for ( int degree : degrees )
        {
            if ( degree < 0 ) throw std::out_of_range( "simplex mode degrees must be non-negative" );
            total += degree;
        }
        if ( total > this->order() ) throw std::out_of_range( "simplex mode exceeds the expansion degree" );
        return ::Feel::detail::simplexModeIndex<Dim>( degrees );
    }

    /** @brief Evaluate modal values, retaining both order and sample extents. */
    template<typename Derived>
    [[nodiscard]] auto values( Eigen::MatrixBase<Derived> const& points ) const
    {
        auto result = this->template makeValues<Derived::ColsAtCompileTime>( points.cols() );
        Workspace<0, Derived::ColsAtCompileTime> workspace{ RuntimeOrder{ this->order() } };
        valuesInto( result, points, workspace );
        return result;
    }
    /** @brief Evaluate reference gradients directly from the regular recurrence. */
    template<typename Derived>
    [[nodiscard]] auto gradients( Eigen::MatrixBase<Derived> const& points ) const
    {
        auto result = this->template makeGradient<Derived::ColsAtCompileTime>( points.cols() );
        Workspace<1, Derived::ColsAtCompileTime> workspace{ RuntimeOrder{ this->order() } };
        gradientsInto( result, points, workspace );
        return result;
    }
    /** @brief Evaluate actual reference Hessians, with symmetric pairs copied exactly. */
    template<typename Derived>
    [[nodiscard]] auto hessians( Eigen::MatrixBase<Derived> const& points ) const
    {
        auto result = this->template makeHessian<Derived::ColsAtCompileTime>( points.cols() );
        Workspace<2, Derived::ColsAtCompileTime> workspace{ RuntimeOrder{ this->order() } };
        hessiansInto( result, points, workspace );
        return result;
    }
    /** @brief Fill reusable value storage without constructing derivative operators. */
    template<int Samples, int WorkSamples, typename Derived>
    void valuesInto( ValuesFor<Samples>& output, Eigen::MatrixBase<Derived> const& points,
                     Workspace<0, WorkSamples>& workspace ) const
    {
        tabulate<Evaluation::Values>( std::array{ &output }, points, workspace );
    }
    /** @brief Fill reusable gradient storage through the same polynomial recurrence. */
    template<int Samples, int WorkSamples, typename Derived>
    void gradientsInto( GradientFor<Samples>& output, Eigen::MatrixBase<Derived> const& points,
                        Workspace<1, WorkSamples>& workspace ) const
    {
        tabulate<Evaluation::Gradients>( ::Feel::detail::matrixBufferPointers( output ), points, workspace );
    }
    /** @brief Fill reusable Hessian storage with fixed-dimensional derivative records. */
    template<int Samples, int WorkSamples, typename Derived>
    void hessiansInto( HessianFor<Samples>& output, Eigen::MatrixBase<Derived> const& points,
                       Workspace<2, WorkSamples>& workspace ) const
    {
        tabulate<Evaluation::Hessians>( ::Feel::detail::matrixBufferPointers( output ), points, workspace );
    }

    /** @brief Fill any requested combination through one regular-coordinate recurrence. */
    template <typename Tables, typename Derived, int WorkSamples>
    void tabulateInto( Tables& output, Eigen::MatrixBase<Derived> const& points,
                       Workspace<::Feel::detail::derivativeOrder( Tables::request ), WorkSamples>& workspace ) const
    {
        static_assert( Tables::dimension == Dim );
        tabulate<Tables::request>( output.pointers(), points, workspace );
    }

private:
    /** @brief Visit the preserved degree ordering once, independently of sample count. */
    template<typename Operation>
    void forEachMode( Operation&& operation ) const
    {
        Eigen::Index row = 0;
        for ( int total = 0; total <= this->order(); ++total )
            for ( int k = 0; k <= ( Dim == 3 ? total : 0 ); ++k )
                for ( int j = 0; j <= ( Dim >= 2 ? total - k : 0 ); ++j )
                {
                    MultiIndex index{};
                    index[0] = total - j - k;
                    if constexpr ( Dim >= 2 ) index[1] = j;
                    if constexpr ( Dim == 3 ) index[2] = k;
                    operation( row++, index );
                }
    }

    /** @brief Compute the invariant modal scale once per reusable workspace. */
    [[nodiscard]] static T modeNormalization( MultiIndex const& degrees )
    {
        if constexpr ( !Normalized ) return T( 1 );
        else
        {
            T scale = T( degrees[0] ) + T( .5 );
            if constexpr ( Dim >= 2 ) scale *= T( degrees[0] + degrees[1] + 1 );
            if constexpr ( Dim == 3 ) scale *= T( degrees[0] + degrees[1] + degrees[2] ) + T( 1.5 );
            using std::sqrt;
            return sqrt( scale );
        }
    }

    /** @brief Extend a lower-dimensional polynomial by the homogeneous Jacobi recurrence.
     * Hn=(constant*t+linear*u)Hn-1-previous*t*t*Hn-2.
     * All division is by degree/parameter constants, never by a coordinate.
     */
    template<int Axis, unsigned DerivativeOrder, int Samples>
    void extend( std::array<::Feel::detail::PolynomialJet<Dim, DerivativeOrder, T>, Dim> const& coordinates,
                 Workspace<DerivativeOrder, Samples>& workspace ) const
    {
        using Jet = ::Feel::detail::PolynomialJet<Dim, DerivativeOrder, T>;
        auto coordinate = [&]( int d ) { return d < Dim ? coordinates[d] : Jet( T( -1 ) ); };
        auto const t = T( -.5 ) * ( coordinate( Axis + 1 ) + coordinate( Axis + 2 ) );
        auto const u = coordinate( Axis ) + Jet( T( 1 ) ) - t;
        auto const square = t * t;
        int const degree = this->order();
        for ( int i = 0; i <= ( Axis >= 1 ? degree : 0 ); ++i )
            for ( int j = 0; j <= ( Axis >= 2 ? degree - i : 0 ); ++j )
            {
                MultiIndex index{};
                if constexpr ( Axis >= 1 ) index[0] = i;
                if constexpr ( Axis >= 2 ) index[1] = j;
                int const limit = degree - i - j;
                if ( limit == 0 ) continue;
                T const alpha = T( 2 * ( i + j ) + Axis );
                auto const first = ::Feel::detail::jacobiFirstCoefficients( alpha, T( 0 ) );
                auto previous = workspace.jet( ::Feel::detail::simplexModeIndex<Dim>( index ) );
                auto current = previous * ( first[0] * t + first[1] * u );
                index[Axis] = 1;
                workspace.setJet( ::Feel::detail::simplexModeIndex<Dim>( index ), current );
                for ( int n = 2; n <= limit; ++n )
                {
                    auto const c = ::Feel::detail::jacobiCoefficients( n, alpha, T( 0 ) );
                    auto next = ( c.constant * t + c.linear * u ) * current - c.previous * ( square * previous );
                    previous = current;
                    current = next;
                    index[Axis] = n;
                    workspace.setJet( ::Feel::detail::simplexModeIndex<Dim>( index ), current );
                }
            }
    }

    /** @brief Shared orientation, ownership, degree ordering and derivative-output path. */
    template <Evaluation Request, unsigned DerivativeOrder, typename Buffer, std::size_t Count, typename Derived, int WorkSamples>
    void tabulate( std::array<Buffer*, Count> output, Eigen::MatrixBase<Derived> const& points,
                   Workspace<DerivativeOrder, WorkSamples>& workspace ) const
    {
        static_assert( Count == ::Feel::detail::derivativeTableCount<Dim>( Request ) );
        static_assert( DerivativeOrder == ::Feel::detail::derivativeOrder( Request ) );

        // Validate the reference dimension and effective degree before evaluation.
        if ( points.rows() != Dim ) throw std::invalid_argument( "simplex point dimension mismatch" );
        if ( workspace.order() != this->order() ) throw std::invalid_argument( "simplex workspace degree mismatch" );

        // Snapshot coordinates before resizing outputs that may alias the input.
        auto const samples = points.cols();
        workspace.M_coordinates.resizeForOverwrite( Dim, samples );
        workspace.M_coordinates.matrix() = points.template cast<T>();

        // Retain the existing buffers and select only the necessary jet algebra.
        auto const count = this->basisCount();
        for ( auto* buffer : output ) buffer->resizeForOverwrite( count, samples );
        using Jet = ::Feel::detail::PolynomialJet<Dim, DerivativeOrder, T>;

        // Borrow stable output maps once before the point and mode loops.
        auto matrices = [&]<std::size_t... Slot>( std::index_sequence<Slot...> )
        { return std::array{ output[Slot]->matrix()... }; }( std::make_index_sequence<Count>{} );

        // Run one regular-coordinate recurrence for every sample and requested group.
        for ( Eigen::Index q = 0; q < samples; ++q )
        {
            std::array<Jet, Dim> coordinates;
            for ( int d = 0; d < Dim; ++d )
                coordinates[d] = Jet::coordinate( workspace.M_coordinates.matrix()( d, q ), d );
            workspace.setJet( 0, Jet( T( 1 ) ) );
            [&]<std::size_t... Axis>( std::index_sequence<Axis...> )
            { ( extend<Axis>( coordinates, workspace ), ... ); }( std::make_index_sequence<Dim>{} );

            // Extract actual reference derivatives from the shared normalized jets.
            for ( Eigen::Index row = 0; row < count; ++row )
            {
                auto const jet = workspace.jet( row );
                auto const scale = workspace.M_jets.matrix()( Jet::componentCount, row );
                constexpr std::size_t gradientStart = includes( Request, Evaluation::Values );
                constexpr std::size_t hessianStart = gradientStart + Dim * includes( Request, Evaluation::Gradients );

                // Keep requested derivative extraction in a flat numerical loop.
                for ( std::size_t slot = 0; slot < Count; ++slot )
                {
                    T value{};

                    // Values occupy the first slot when requested.
                    if constexpr ( includes( Request, Evaluation::Values ) )
                        if ( slot == 0 ) value = jet.value();

                    // Gradient slots carry one reference coordinate direction.
                    if constexpr ( includes( Request, Evaluation::Gradients ) )
                        if ( slot >= gradientStart && slot < hessianStart )
                            value = jet.gradient( slot - gradientStart );

                    // Hessian slots retain the established ordered direction pairs.
                    if constexpr ( includes( Request, Evaluation::Hessians ) )
                        if ( slot >= hessianStart )
                            value = jet.hessian( ( slot - hessianStart ) / Dim, ( slot - hessianStart ) % Dim );

                    matrices[slot]( row, q ) = scale * value;
                }
            }
        }
    }
};
} // namespace Feel
#endif // FEELPP_FEELPOLY_SIMPLEXTABULATION_HPP
