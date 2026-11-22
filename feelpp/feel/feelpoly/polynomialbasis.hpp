// SPDX-FileCopyrightText: 2026 University of Strasbourg
// SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <array>
#include <feel/feelalg/eigenmap.hpp>
#include <feel/feelpoly/order.hpp>
#include <feel/feelpoly/tabulation.hpp>
#include <memory>
#include <mutex>
#include <stdexcept>

namespace Feel::detail
{
/** @brief Modal descriptors exposing immutable all-direction requested-order ownership. */
template <typename Basis>
concept SharedReferenceOperators = requires( Basis const& basis ) {
    basis.referenceOperators( RuntimeOrder{ 0 } );
};

/** @brief Storage-free derivative policy for a matching static basis. */
template <class Basis, bool Cached, bool Shared = SharedReferenceOperators<Basis>>
struct PolynomialDerivativeCache
{
};

/** @brief Shared, immutable-after-initialization operators for an effective order. */
template <class Basis>
struct PolynomialDerivativeCache<Basis, true, false>
{
    /** @brief Operators initialized once per coordinate direction. */
    struct Operators
    {
        /** Initialization guards make concurrent const evaluation safe. */
        std::array<std::once_flag, Basis::nDim> ready;
        /** Modal differentiation matrices at the effective degree. */
        std::array<typename Basis::matrix_type, Basis::nDim> matrices;
    };
    /** Copies share operators; changing order detaches this state. */
    std::shared_ptr<Operators> operators = std::make_shared<Operators>();
};

/** @brief Retain one immutable all-direction handle without copying its coefficient matrices. */
template <typename Basis>
struct PolynomialDerivativeCache<Basis, true, true>
{
    /** Lazy ownership copied by polynomial value semantics and detached on order changes. */
    struct Operators
    {
        /** Publish one complete reference handle exactly once for this effective order. */
        std::once_flag ready;
        /** Family-owned immutable matrices with lifetime independent of other requested orders. */
        decltype( std::declval<Basis const&>().referenceOperators( RuntimeOrder{ 0 } ) ) reference;
    };
    /** Copies share first access and ownership; an order change creates a new state. */
    std::shared_ptr<Operators> operators = std::make_shared<Operators>();
};

/**
 * @brief Common order-aware evaluation and operator policy for polynomials.
 *
 * Matching static bases retain their native evaluation and fixed operators.
 * Runtime orders and converted static extents use the same requested-order
 * kernels. Operator references remain valid until the owning state is assigned,
 * changes order or is destroyed; copies share immutable operators.
 */
template<class Basis, int OrderSpec>
class PolynomialBasisState : public OrderBase<OrderSpec>,
                             private PolynomialDerivativeCache<Basis, OrderSpec == Dynamic || OrderSpec != Basis::nOrder>
{
    /** Whether operators must be owned independently of the seed basis. */
    static constexpr bool cached = OrderSpec == Dynamic || OrderSpec != Basis::nOrder;
    /** Order-policy base. */
    using order_base = OrderBase<OrderSpec>;
    /** Optional operator storage. */
    using cache_base = PolynomialDerivativeCache<Basis, cached>;
public:
    /** @brief Construct at the static degree or the dynamic seed degree. */
    PolynomialBasisState() : order_base( RuntimeOrder{OrderSpec == Dynamic ? Basis::nOrder : OrderSpec} ) {}

    /** @brief Update dynamic order or validate a static degree. */
    void setOrder( uint16_type order )
    {
        if constexpr ( OrderSpec == Dynamic )
        {
            if ( order != this->order() )
            {
                this->operators = std::make_shared<typename cache_base::Operators>();
                order_base::setOrder( order );
            }
        }
        else
            (void)order_base{RuntimeOrder{order}};
    }

    /** @brief Evaluate modal values at a batch of points. */
    template<class Points>
    auto evaluate( Basis const& basis, Points const& points ) const
    {
        if constexpr ( cached && requires { basis.evaluate( points, this->order() ); } )
            return basis.evaluate( points, this->order() );
        else
        {
            validateSeed();
            return basis.evaluate( points );
        }
    }

    /** @brief Evaluate first derivatives at a batch of points. */
    template<class Points>
    auto derivate( Basis const& basis, Points const& points ) const
    {
        if constexpr ( cached && requires { basis.derivate( points, this->order() ); } )
            return basis.derivate( points, this->order() );
        else
        {
            validateSeed();
            return basis.derivate( points );
        }
    }

    /** @brief Consume modal values while their native or legacy owner is alive.
     * @param basis Modal evaluator, possibly seeded at a different degree.
     * @param points Column-oriented sample coordinates.
     * @param operation Callable returning an owning result from a borrowed Eigen matrix.
     * @note Static orders retain their extent even when the seed degree differs.
     */
    template<class Points, class Operation>
    auto withValues( Basis const& basis, Points const& points, Operation&& operation ) const
    {
        if constexpr ( requires { basis.template values<OrderSpec>( points, RuntimeOrder{ this->order() } ); } )
        {
            auto values = basis.template values<OrderSpec>( points, RuntimeOrder{ this->order() } );
            return operation( values.matrix() );
        }
        else
        {
            auto values = evaluate( basis, points );
            return operation( eigenMap( values ) );
        }
    }

    /** @brief Consume reference gradients without materializing native tables as legacy matrices.
     * @param basis Modal evaluator, possibly seeded at a different degree.
     * @param points Column-oriented sample coordinates.
     * @param operation Callable receiving a direction and its borrowed Eigen table.
     */
    template<class Points, class Operation>
    void forEachGradient( Basis const& basis, Points const& points, Operation&& operation ) const
    {
        if constexpr ( requires { basis.template gradients<OrderSpec>( points, RuntimeOrder{ this->order() } ); } )
        {
            auto gradients = basis.template gradients<OrderSpec>( points, RuntimeOrder{ this->order() } );
            for ( int d = 0; d < Basis::nDim; ++d ) operation( d, gradients[d].matrix() );
        }
        else
        {
            auto gradients = derivate( basis, points );
            for ( int d = 0; d < Basis::nDim; ++d ) operation( d, eigenMap( gradients[d] ) );
        }
    }

    /** Whether this modal facade supplies direct, operator-free reference Hessians. */
    template<class Points>
    static constexpr bool hasDirectHessians = requires( Basis const& basis, Points const& points )
    { basis.template hessians<OrderSpec>( points, RuntimeOrder{ 0 } ); };

    /** @brief Select direct sampled Hessians except for measured low-order tensor operators. */
    [[nodiscard]] bool useDirectHessians() const noexcept
    {
        // Linear tensor expansions favor the existing coefficient operators on warm calls.
        if constexpr ( requires { Basis::convex_type::is_simplex; } )
            return Basis::convex_type::is_simplex || this->order() != 1;
        else
            return true;
    }

    /** @brief Consume native Hessians while their owning tables remain alive.
     * @param basis Modal facade supplying a direct second-derivative recurrence.
     * @param points Column-oriented sample coordinates.
     * @param operation Callable consuming (i,j,table); no coefficient operators are requested.
     */
    template<class Points, class Operation>
        requires hasDirectHessians<Points>
    void forEachHessian( Basis const& basis, Points const& points, Operation&& operation ) const
    {
        auto hessians = basis.template hessians<OrderSpec>( points, RuntimeOrder{ this->order() } );
        for ( int i = 0; i < Basis::nDim; ++i )
            for ( int j = 0; j < Basis::nDim; ++j )
                operation( i, j, hessians[i * Basis::nDim + j].matrix() );
    }

    /** @brief Return the owned differentiation operator in a coordinate direction. */
    typename Basis::matrix_type const& d( Basis const& basis, uint16_type direction ) const
    {
        if ( direction >= Basis::nDim )
            throw std::out_of_range("polynomial derivative direction is outside the reference dimension");
        if constexpr ( cached )
        {
            if constexpr ( SharedReferenceOperators<Basis> )
            {
                auto initialize = [&]
                {
                    this->operators->reference = basis.referenceOperators( RuntimeOrder{ this->order() } );
                };
                std::call_once( this->operators->ready, initialize );
                return this->operators->reference->matrix( direction );
            }
            else
            {
                std::call_once( this->operators->ready[direction], [&] {
                    if constexpr ( requires { basis.d( direction, this->order() ); } )
                        this->operators->matrices[direction] = basis.d( direction, this->order() );
                    else
                    {
                        validateSeed();
                        this->operators->matrices[direction] = basis.d( direction );
                    }
                } );
                return this->operators->matrices[direction];
            }
        }
        else
            return basis.d( direction );
    }
private:
    /** @brief Reject unavailable requested-order kernels before matrix algebra. */
    void validateSeed() const
    {
        if constexpr ( cached )
            if ( this->order() != Basis::nOrder )
                throw std::invalid_argument("basis does not support the requested polynomial order");
    }
};
} // namespace Feel::detail
