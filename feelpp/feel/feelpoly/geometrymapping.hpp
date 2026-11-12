/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file geometrymapping.hpp
 * @brief Representation-independent polynomial geometric mapping.
 */
#if !defined( FEELPP_FEELPOLY_GEOMETRYMAPPING_HPP )
#define FEELPP_FEELPOLY_GEOMETRYMAPPING_HPP 1

#include <feel/feelpoly/geometrycoefficients.hpp>

#include <array>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Feel::geometry
{
//! \class Mapping
//! \brief Representation-independent polynomial geometric mapping.
//!
//! The class composes a mathematical representation with its typed coefficient
//! table. Point, Jacobian, and Hessian assembly are shared by nodal Lagrange,
//! Bernstein/Bezier, static order, and runtime order.
//!
//! \tparam RealDim ambient physical dimension
//! \tparam Representation geometry representation policy
template <int RealDim, typename Representation>
class Mapping
{
  public:
    /** Geometry representation policy. */
    using representation_type = Representation;
    /** Scalar value type. */
    using value_type = typename representation_type::value_type;
    /** Reference convex type. */
    using convex_type = typename representation_type::convex_type;
    /** Typed coefficient table. */
    using coefficients_type = CoefficientTable<representation_type>;
    /** Reference points stored by column. */
    using points_type = typename representation_type::points_type;
    /** Dynamic physical matrix type. */
    using matrix_type = typename coefficients_type::matrix_type;
    /** One physical Jacobian per evaluation point. */
    using jacobian_type = Eigen::Matrix<value_type, RealDim, representation_type::dimension()>;
    /** Collection of pointwise physical Jacobians. */
    using jacobians_type = std::vector<jacobian_type>;
    /** Reference dimension. */
    static constexpr int S_dim = representation_type::dimension();

    /**
     * @brief Construct a geometric mapping from typed coefficients.
     * @param coefficients representation-bound physical coefficients
     * @throws std::invalid_argument if the coefficient row count is not `RealDim`
     */
    explicit Mapping( coefficients_type coefficients )
        : M_coefficients( std::move( coefficients ) )
    {
        if ( M_coefficients.values().rows() != RealDim )
            throw std::invalid_argument( "Geometry coefficient rows must equal RealDim" );
    }

    /** @brief Return the effective polynomial order. */
    [[nodiscard]] uint16_type order() const noexcept
    {
        return M_coefficients.representation().order();
    }

    /** @brief Return whether the effective mapping order is affine. */
    [[nodiscard]] bool isLinear() const noexcept { return this->order() == 1; }

    /** @brief Return whether the representation order is known at compile time. */
    [[nodiscard]] static consteval bool isOrderStatic() noexcept
    {
        return representation_type::isOrderStatic();
    }

    /** @brief Return whether the compile-time mapping is affine. */
    [[nodiscard]] static consteval bool isStaticallyLinear() noexcept
    {
        if constexpr ( representation_type::isOrderStatic() )
            return representation_type::staticOrder() == 1;
        else
            return false;
    }

    /** @brief Return the reference dimension. */
    [[nodiscard]] static consteval int dim() noexcept { return S_dim; }

    /** @brief Return the ambient physical dimension. */
    [[nodiscard]] static consteval int realDim() noexcept { return RealDim; }

    /** @brief Return the representation-bound coefficient table. */
    [[nodiscard]] coefficients_type const& geometryCoefficients() const noexcept
    {
        return M_coefficients;
    }

    /**
     * @brief Map reference points to physical points.
     * @param points reference points stored by column
     * @return `RealDim`-by-point physical coordinate matrix
     */
    [[nodiscard]] matrix_type transform( points_type const& points ) const
    {
        return M_coefficients.values() *
               M_coefficients.representation().basisValues( points );
    }

    /**
     * @brief Assemble pointwise physical Jacobians.
     * @param points reference points stored by column
     * @return one `RealDim`-by-`nDim` Jacobian per point
     */
    [[nodiscard]] jacobians_type jacobians( points_type const& points ) const
    {
        auto const derivatives = M_coefficients.representation().basisGradients( points );
        std::array<matrix_type, S_dim> physicalDerivatives;
        for ( int d = 0; d < S_dim; ++d )
            physicalDerivatives[d] = M_coefficients.values() * derivatives[d];
        jacobians_type result( points.cols() );
        for ( Eigen::Index q = 0; q < points.cols(); ++q )
            for ( int d = 0; d < S_dim; ++d )
                result[q].col( d ) = physicalDerivatives[d].col( q );
        return result;
    }

    /**
     * @brief Assemble physical second derivatives of the mapping.
     * @param points reference points stored by column
     * @return row-major `(i,j)` array of `RealDim`-by-point matrices
     */
    [[nodiscard]] std::array<matrix_type, S_dim * S_dim>
    hessians( points_type const& points ) const
    {
        auto const derivatives = M_coefficients.representation().basisHessians( points );
        std::array<matrix_type, S_dim * S_dim> result;
        for ( int derivative = 0; derivative < S_dim * S_dim; ++derivative )
            result[derivative] = M_coefficients.values() * derivatives[derivative];
        return result;
    }

  private:
    /** Representation-bound physical coefficients. */
    coefficients_type M_coefficients;
};

//!
//! \brief Construct a representation-independent geometric mapping.
//! \tparam RealDim ambient physical dimension
//! \tparam Representation geometry representation policy
//! \param representation representation object
//! \param values component-by-coefficient physical values
//! \return composed geometric mapping
template <int RealDim, typename Representation>
[[nodiscard]] auto makeMapping(
    Representation representation,
    typename CoefficientTable<Representation>::matrix_type values )
{
    return Mapping<RealDim, Representation>(
        CoefficientTable<Representation>( std::move( representation ), std::move( values ) ) );
}

} // namespace Feel::geometry

#endif // FEELPP_FEELPOLY_GEOMETRYMAPPING_HPP
