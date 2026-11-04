// SPDX-FileCopyrightText: 2026 University of Strasbourg
// SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <array>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <feel/feelpoly/order.hpp>

namespace Feel::detail
{
/** @brief Obtain the expansion degree, which may differ from a FE family index. */
template<class Source>
[[nodiscard]] uint16_type polynomialOrder( Source const& source )
{
    if constexpr ( requires { source.polynomialDegree(); } )
        return source.polynomialDegree();
    else if constexpr ( requires { source.order(); } )
        return source.order();
    else
        return staticOrderValue<Source>();
}

/** @brief Storage-free derivative policy for a matching static basis. */
template<class Basis, bool Cached>
struct PolynomialDerivativeCache {};

/** @brief Shared, immutable-after-initialization operators for an effective order. */
template<class Basis>
struct PolynomialDerivativeCache<Basis, true>
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

    /** @brief Return the owned differentiation operator in a coordinate direction. */
    typename Basis::matrix_type const& d( Basis const& basis, uint16_type direction ) const
    {
        if ( direction >= Basis::nDim )
            throw std::out_of_range("polynomial derivative direction is outside the reference dimension");
        if constexpr ( cached )
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
}
