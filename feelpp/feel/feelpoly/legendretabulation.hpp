/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/** @file legendretabulation.hpp
 * @brief Shared Eigen tensor-product Legendre values, gradients and Hessians.
 */
#ifndef FEELPP_FEELPOLY_LEGENDRETABULATION_HPP
#define FEELPP_FEELPOLY_LEGENDRETABULATION_HPP 1
#include <array>
#include <cmath>
#include <feel/feelpoly/jacobirecurrence.hpp>
#include <feel/feelpoly/tabulation.hpp>
#include <utility>

namespace Feel
{
/** @brief Tensor-product Legendre basis with independent order and sample extents.
 * @tparam Dim Reference dimension, one to three.
 * @tparam OrderSpec Compile-time degree or `Dynamic`.
 * @tparam T Scalar type.
 * @tparam Normalized Whether each one-dimensional factor is L2-normalized.
 *
 * The last coordinate varies fastest in the basis index. All derivative orders
 * use the same tensor-product kernel and shifted Jacobi factors. Reusable
 * workspaces own their univariate tables; no borrowed sample views are retained.
 */
template<int Dim, int OrderSpec, typename T = double, bool Normalized = true>
class TensorLegendreBasis : public ::Feel::detail::PolynomialStorage<
    T, Dim, OrderSpec, ::Feel::detail::NodalCellKind::Hypercube>
{
    using super = ::Feel::detail::PolynomialStorage<T, Dim, OrderSpec, ::Feel::detail::NodalCellKind::Hypercube>;
    using Jacobi = JacobiPolynomial<OrderSpec, T>;
public:
    /** Modal basis rows by independent sample columns. */
    template<int Samples = Dynamic> using ValuesFor = typename super::template ValuesFor<Samples>;
    /** Values with a runtime sample count. */
    using Values = ValuesFor<>;
    /** One owning value table per reference derivative direction. */
    template<int Samples = Dynamic> using GradientFor = typename super::template GradientFor<Samples>;
    /** Gradient with a runtime sample count. */
    using Gradient = GradientFor<>;
    /** One owning value table per ordered direction pair. */
    template<int Samples = Dynamic> using HessianFor = typename super::template HessianFor<Samples>;
    /** Hessian with a runtime sample count. */
    using Hessian = HessianFor<>;

    /** @brief Univariate tables reused by one selected derivative operation.
     * @tparam DerivativeOrder Zero for values, one for gradients, two for Hessians.
     */
    template<unsigned DerivativeOrder, int Samples = Dynamic>
    class Workspace : public OrderBase<OrderSpec>
    {
        static_assert( DerivativeOrder <= 2 );
        friend class TensorLegendreBasis;
        /** Bound each univariate buffer under the aggregate workspace budget. */
        using UnivariateTable = typename Jacobi::template TableFor<Samples,
            ::Feel::detail::defaultInlineMatrixBytes / ( Dim * ( DerivativeOrder + 1 ) )>;
    public:
        /** @brief Bind workspace rows to the effective degree; samples may change. */
        explicit Workspace( RuntimeOrder degree = RuntimeOrder{ OrderSpec == Dynamic ? 0 : OrderSpec } )
            : OrderBase<OrderSpec>( degree ),
              M_values( ::Feel::detail::makeMatrixBuffers<UnivariateTable, Dim>(
                  degree.value + 1, Samples == Dynamic ? 0 : Samples, degree.value + 1 ) ),
              M_first( ::Feel::detail::makeMatrixBuffers<UnivariateTable, DerivativeOrder >= 1 ? Dim : 0>(
                  degree.value + 1, Samples == Dynamic ? 0 : Samples, degree.value + 1 ) ),
              M_second( ::Feel::detail::makeMatrixBuffers<UnivariateTable, DerivativeOrder >= 2 ? Dim : 0>(
                  degree.value + 1, Samples == Dynamic ? 0 : Samples, degree.value + 1 ) ) {}
    private:
        /** Value factors for each fixed spatial direction. */
        std::array<UnivariateTable, Dim> M_values;
        /** First derivatives, absent from value-only workspaces. */
        std::array<UnivariateTable, DerivativeOrder >= 1 ? Dim : 0> M_first;
        /** Second derivatives, absent unless a Hessian is requested. */
        std::array<UnivariateTable, DerivativeOrder >= 2 ? Dim : 0> M_second;
    };

    /** @brief Construct static or runtime degree through the same checked policy. */
    explicit TensorLegendreBasis( RuntimeOrder degree = RuntimeOrder{ OrderSpec == Dynamic ? 0 : OrderSpec } )
        : super( degree ) {}

    /** @brief Evaluate tensor-product modal values at column-oriented points. */
    template<typename Derived>
    [[nodiscard]] auto values( Eigen::MatrixBase<Derived> const& points ) const
    {
        auto const modeCount = this->basisCount();
        auto result = ValuesFor<Derived::ColsAtCompileTime>::forOverwrite(
            modeCount, points.cols(), modeCount, Dynamic );
        Workspace<0, Derived::ColsAtCompileTime> workspace{ RuntimeOrder{ this->order() } };
        valuesInto( result, points, workspace );
        return result;
    }
    /** @brief Evaluate every reference gradient direction. */
    template<typename Derived>
    [[nodiscard]] auto gradients( Eigen::MatrixBase<Derived> const& points ) const
    {
        auto result = this->template makeGradient<Derived::ColsAtCompileTime>( points.cols() );
        Workspace<1, Derived::ColsAtCompileTime> workspace{ RuntimeOrder{ this->order() } };
        gradientsInto( result, points, workspace );
        return result;
    }
    /** @brief Evaluate ordered Hessian directions, indexed by i*Dim+j. */
    template<typename Derived>
    [[nodiscard]] auto hessians( Eigen::MatrixBase<Derived> const& points ) const
    {
        auto result = this->template makeHessian<Derived::ColsAtCompileTime>( points.cols() );
        Workspace<2, Derived::ColsAtCompileTime> workspace{ RuntimeOrder{ this->order() } };
        hessiansInto( result, points, workspace );
        return result;
    }
    /** @brief Evaluate values using reusable output and univariate storage. */
    template<int Samples, int WorkSamples, typename Derived>
    void valuesInto( ValuesFor<Samples>& output, Eigen::MatrixBase<Derived> const& points,
                     Workspace<0, WorkSamples>& workspace ) const
    {
        tabulate<Evaluation::Values>( std::array{ &output }, points, workspace );
    }
    /** @brief Evaluate gradients using reusable output and univariate storage. */
    template<int Samples, int WorkSamples, typename Derived>
    void gradientsInto( GradientFor<Samples>& output, Eigen::MatrixBase<Derived> const& points,
                        Workspace<1, WorkSamples>& workspace ) const
    {
        tabulate<Evaluation::Gradients>( ::Feel::detail::matrixBufferPointers( output ), points, workspace );
    }
    /** @brief Evaluate Hessians using reusable output and univariate storage. */
    template<int Samples, int WorkSamples, typename Derived>
    void hessiansInto( HessianFor<Samples>& output, Eigen::MatrixBase<Derived> const& points,
                       Workspace<2, WorkSamples>& workspace ) const
    {
        tabulate<Evaluation::Hessians>( ::Feel::detail::matrixBufferPointers( output ), points, workspace );
    }

    /** @brief Fill requested values and derivatives from one set of univariate factors. */
    template <typename Tables, typename Derived, int WorkSamples>
    void tabulateInto( Tables& output, Eigen::MatrixBase<Derived> const& points,
                       Workspace<::Feel::detail::derivativeOrder( Tables::request ), WorkSamples>& workspace ) const
    {
        static_assert( Tables::dimension == Dim );
        tabulate<Tables::request>( output.pointers(), points, workspace );
    }

  private:
    /** @brief Shared tensor indexing, normalization and product-rule kernel. */
    template <Evaluation Request, unsigned DerivativeOrder, typename Buffer, std::size_t Count, typename Derived, int WorkSamples>
    void tabulate( std::array<Buffer*, Count> output, Eigen::MatrixBase<Derived> const& points,
                   Workspace<DerivativeOrder, WorkSamples>& workspace ) const
    {
        static_assert( Count == ::Feel::detail::derivativeTableCount<Dim>( Request ) );
        static_assert( DerivativeOrder == ::Feel::detail::derivativeOrder( Request ) );

        // Validate coordinates and retained scratch before touching any output.
        if ( points.rows() != Dim ) throw std::invalid_argument( "Legendre point dimension mismatch" );
        if ( workspace.order() != this->order() ) throw std::invalid_argument( "Legendre workspace degree mismatch" );

        auto const samples = points.cols();

        // Fill one Jacobi derivative table while retaining expression alias protection.
        auto fill = [&]<unsigned Derivative>( auto& buffer, auto const& coordinates )
        {
            if constexpr ( Derived::Flags & Eigen::DirectAccessBit )
            {
                // Workspace cells are private and cannot alias caller coordinates.
                buffer.resizeForOverwrite( this->order() + 1, samples );
                auto target = buffer.matrix();
                ::Feel::detail::jacobiTableRows<Derivative, OrderSpec>(
                    this->order(), T( 0 ), T( 0 ), target, coordinates.template cast<T>().array() );
            }
            else
            {
                // Keep snapshot handling for arbitrary coordinate expressions.
                Jacobi polynomial{ RuntimeOrder{ this->order() } };
                if constexpr ( Derivative == 0 )
                    polynomial.valuesInto( buffer, coordinates );
                else
                    polynomial.template derivativesInto<Derivative>( buffer, coordinates );
            }
        };

        // Read every coordinate before resizing outputs, which may own those coordinates.
        if ( samples != 0 )
            for ( int d = 0; d < Dim; ++d )
            {
                fill.template operator()<0>( workspace.M_values[d], points.row( d ) );
                if constexpr ( DerivativeOrder >= 1 )
                    fill.template operator()<1>( workspace.M_first[d], points.row( d ) );
                if constexpr ( DerivativeOrder >= 2 )
                    fill.template operator()<2>( workspace.M_second[d], points.row( d ) );
            }

        // Reuse result buffers now that all coordinate factors are in private scratch.
        auto const modeCount = this->basisCount();
        for ( auto* buffer : output )
            buffer->resizeForOverwrite( modeCount, samples );

        // Empty outputs have their final shape without forming any coordinate or result row view.
        if ( samples == 0 ) return;

        // Tensor products use one numerical loop after factor setup and shape validation.
        combineFactors<Request>( output, samples, workspace );
    }

    /** @brief Form requested tensor products from private one-dimensional factors. */
    template <Evaluation Request, unsigned DerivativeOrder, typename Buffer, std::size_t Count, int WorkSamples>
    EIGEN_ALWAYS_INLINE void combineFactors( std::array<Buffer*, Count> const& output, Eigen::Index samples,
                                             Workspace<DerivativeOrder, WorkSamples> const& workspace ) const
    {
        auto const modeCount = this->basisCount();
        int const radix = this->order() + 1;

        // Decode tensor modes in the established ordering and apply their normalization.
        for ( Eigen::Index row = 0; row < modeCount; ++row )
        {
            std::array<int, Dim> degrees;
            Eigen::Index remaining = row;
            T normalization = T( 1 );
            for ( int d = Dim - 1; d >= 0; --d )
            {
                degrees[d] = static_cast<int>( remaining % radix );
                remaining /= radix;
                if constexpr ( Normalized ) normalization *= T( degrees[d] ) + T( .5 );
            }
            using std::sqrt;
            if constexpr ( Normalized ) normalization = sqrt( normalization );

            // Keep a flat numerical loop over the requested output tables.
            for ( std::size_t slot = 0; slot < Count; ++slot )
            {
                auto target = output[slot]->matrix();
                std::array<decltype( &workspace.M_values[0] ), Dim> factors{};
                bool vanishes = false;

                // Select one private univariate factor for each reference coordinate.
                for ( int d = 0; d < Dim; ++d )
                {
                    unsigned const count = ::Feel::detail::derivativeExponent<Dim, Request>( slot, d );

                    // A derivative beyond an axis mode's degree vanishes identically.
                    if constexpr ( includes( Request, Evaluation::Hessians ) )
                        if ( degrees[d] < count )
                        {
                            vanishes = true;
                            break;
                        }

                    factors[d] = &workspace.M_values[d];
                    if constexpr ( DerivativeOrder >= 1 )
                        if ( count == 1 ) factors[d] = &workspace.M_first[d];
                    if constexpr ( DerivativeOrder >= 2 )
                        if ( count == 2 ) factors[d] = &workspace.M_second[d];
                }

                // A vanishing mode derivative needs one zero fill and no factor reads.
                if ( vanishes )
                {
                    target.row( row ).setZero();
                    continue;
                }

                // Fuse normalization and all axis factors into one output assignment.
                auto assignProduct = [&]<std::size_t... Axis>( std::index_sequence<Axis...> )
                {
                    target.row( row ).array() = ( normalization * ... * factors[Axis]->matrix().row( degrees[Axis] ).array() );
                };
                assignProduct( std::make_index_sequence<Dim>{} );
            }
        }
    }
};
} // namespace Feel
#endif // FEELPP_FEELPOLY_LEGENDRETABULATION_HPP
