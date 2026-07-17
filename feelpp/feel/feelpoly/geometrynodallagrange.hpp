/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file geometrynodallagrange.hpp
 * @brief Nodal Lagrange representation of polynomial geometric mappings.
 */
#if !defined( FEELPP_FEELPOLY_GEOMETRYNODALLAGRANGE_HPP )
#define FEELPP_FEELPOLY_GEOMETRYNODALLAGRANGE_HPP 1

#include <feel/feelpoly/bernsteinbezier.hpp>

#include <Eigen/Core>
#include <Eigen/LU>

#include <array>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace Feel::geometry
{
//!
//! \class NodalLagrange
//! \brief Nodal interpolation representation backed by the Bernstein kernel.
//!
//! The inverse Bernstein Vandermonde matrix defines the nodal cardinal basis.
//! It also provides explicit and mathematically correct conversion between
//! interpolation values and Bezier control points.
//!
//! \tparam Dim reference dimension
//! \tparam Order static polynomial order or `Dynamic`
//! \tparam Entity reference-entity family
//! \tparam T scalar type
template <int Dim, int Order = 1,
          template <int, int, int> class Entity = Simplex, typename T = double>
class NodalLagrange : public OrderBase<Order>
{
    static_assert( Order == Dynamic || Order >= 1,
                   "A geometry representation order must be at least one" );
    using order_base_type = OrderBase<Order>;

  public:
    /** Scalar value type. */
    using value_type = T;
    /** Reference convex used by the representation. */
    using convex_type = Entity<Dim, ( Order == Dynamic ? 1 : Order ), Dim>;
    /** Matrix whose columns are reference interpolation nodes. */
    using points_type = Eigen::Matrix<T, Dim, Eigen::Dynamic>;
    /** Matrix whose rows are basis functions and columns are points. */
    using basis_matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    /** @brief Return the mathematical coefficient semantics. */
    [[nodiscard]] static consteval CoefficientSemantics coefficientSemantics() noexcept
    {
        return CoefficientSemantics::InterpolationPoint;
    }
    /** @brief Return the reference dimension. */
    [[nodiscard]] static consteval int dimension() noexcept { return Dim; }
    /** @brief Return whether the polynomial order is selected at runtime. */
    [[nodiscard]] static consteval bool isOrderDynamic() noexcept { return Order == Dynamic; }
    /** @brief Return whether the polynomial order is known at compile time. */
    [[nodiscard]] static consteval bool isOrderStatic() noexcept { return !isOrderDynamic(); }
    /** @brief Return the template order, including `Dynamic`. */
    [[nodiscard]] static consteval int staticOrder() noexcept { return Order; }
    /** @brief Return whether exact subdivision is supported. */
    [[nodiscard]] static consteval bool supportsSubdivision() noexcept { return false; }
    /** @brief Return whether coefficient convex-hull bounds are available. */
    [[nodiscard]] static consteval bool supportsConvexHull() noexcept { return false; }

    /**
     * @brief Construct from the interpolation nodes of a static representation.
     * @param referenceNodes reference nodes stored by column
     */
    explicit NodalLagrange( points_type referenceNodes )
        requires( Order != Dynamic )
        : M_bernstein(), M_referenceNodes( std::move( referenceNodes ) )
    {
        this->initializeVandermonde();
    }

    /**
     * @brief Construct from runtime order and interpolation nodes.
     * @param order polynomial order
     * @param referenceNodes reference nodes stored by column
     */
    NodalLagrange( RuntimeOrder order, points_type referenceNodes )
        : order_base_type( order ), M_bernstein( order ),
          M_referenceNodes( std::move( referenceNodes ) )
    {
        if constexpr ( Order != Dynamic )
            if ( order.value != Order )
                throw std::invalid_argument( "Runtime order disagrees with static nodal representation" );
        this->initializeVandermonde();
    }

    /** @brief Return the effective polynomial order. */
    [[nodiscard]] constexpr uint16_type order() const noexcept
    {
        return order_base_type::order();
    }

    /** @brief Return the number of interpolation coefficients. */
    [[nodiscard]] size_type coefficientCount() const noexcept
    {
        return M_referenceNodes.cols();
    }

    /** @brief Return a stable human-readable representation name. */
    [[nodiscard]] static constexpr std::string_view name() noexcept
    {
        return "NodalLagrange";
    }

    /** @brief Return the reference interpolation nodes. */
    [[nodiscard]] points_type const& referenceNodes() const noexcept
    {
        return M_referenceNodes;
    }

    /** @brief Evaluate the nodal cardinal basis at reference points. */
    [[nodiscard]] basis_matrix_type basisValues( points_type const& points ) const
    {
        return M_inverseVandermonde * M_bernstein.basisValues( points );
    }

    /** @brief Evaluate reference gradients of the nodal cardinal basis. */
    [[nodiscard]] std::array<basis_matrix_type, Dim>
    basisGradients( points_type const& points ) const
    {
        auto result = M_bernstein.basisGradients( points );
        for ( auto& derivative : result )
            derivative = M_inverseVandermonde * derivative;
        return result;
    }

    /** @brief Evaluate reference Hessians of the nodal cardinal basis. */
    [[nodiscard]] std::array<basis_matrix_type, Dim * Dim>
    basisHessians( points_type const& points ) const
    {
        auto result = M_bernstein.basisHessians( points );
        for ( auto& derivative : result )
            derivative = M_inverseVandermonde * derivative;
        return result;
    }

    /**
     * @brief Convert physical nodal values to Bezier control points.
     * @param nodalValues component-by-node coefficient matrix
     * @return component-by-control-point coefficient matrix
     */
    template <typename Derived>
    [[nodiscard]] basis_matrix_type
    toBezier( Eigen::MatrixBase<Derived> const& nodalValues ) const
    {
        this->validateCoefficientColumns( nodalValues.cols() );
        return nodalValues * M_inverseVandermonde;
    }

    /**
     * @brief Convert Bezier control points to physical nodal values.
     * @param controlPoints component-by-control-point coefficient matrix
     * @return component-by-node coefficient matrix
     */
    template <typename Derived>
    [[nodiscard]] basis_matrix_type
    fromBezier( Eigen::MatrixBase<Derived> const& controlPoints ) const
    {
        this->validateCoefficientColumns( controlPoints.cols() );
        return controlPoints * M_vandermonde;
    }

    /** @brief Return the Bernstein evaluation matrix at the interpolation nodes. */
    [[nodiscard]] basis_matrix_type const& bernsteinVandermonde() const noexcept
    {
        return M_vandermonde;
    }

  private:
    /** Build and factor the Bernstein interpolation matrix. */
    void initializeVandermonde()
    {
        if ( M_referenceNodes.cols() != static_cast<Eigen::Index>( M_bernstein.coefficientCount() ) )
            throw std::invalid_argument( "Nodal geometry has the wrong number of reference nodes" );
        M_vandermonde = M_bernstein.basisValues( M_referenceNodes );
        Eigen::FullPivLU<basis_matrix_type> factorization( M_vandermonde );
        if ( !factorization.isInvertible() )
            throw std::invalid_argument( "Geometry interpolation nodes are not unisolvent" );
        M_inverseVandermonde = factorization.inverse();
    }

    /** Validate a component-by-coefficient table cardinality. */
    void validateCoefficientColumns( Eigen::Index columns ) const
    {
        if ( columns != static_cast<Eigen::Index>( this->coefficientCount() ) )
            throw std::invalid_argument( "Geometry coefficient table has an invalid column count" );
    }

    /** Bernstein representation supplying the common evaluation kernel. */
    BernsteinBezier<Dim, Order, Entity, T> M_bernstein;
    /** Reference interpolation nodes stored by column. */
    points_type M_referenceNodes;
    /** Bernstein basis evaluated at the interpolation nodes. */
    basis_matrix_type M_vandermonde;
    /** Inverse interpolation matrix defining the nodal cardinal basis. */
    basis_matrix_type M_inverseVandermonde;
};

} // namespace Feel::geometry

#endif // FEELPP_FEELPOLY_GEOMETRYNODALLAGRANGE_HPP
