/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 2026-01-02

 Copyright (C) 2026 Feel++ Consortium

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
/**
 * @file concepts.hpp
 * @brief C++20 concepts for Feel++ polynomial spaces and basis functions
 * @author Christophe Prud'homme
 * @date 2026-01-02
 */
#ifndef FEELPP_FEELPOLY_CONCEPTS_HPP
#define FEELPP_FEELPOLY_CONCEPTS_HPP 1

#include <concepts>
#include <feel/feelcore/concepts.hpp>

namespace Feel
{

//
// Polynomial Space Concepts
//

/**
 * @brief A polynomial space (Lagrange, Nedelec, Raviart-Thomas, etc.)
 * 
 * @details Polynomial spaces define the shape functions used in
 * finite element discretizations.
 */
template <typename T>
concept PolynomialSet = requires {
    typename T::value_type;
    typename T::points_type;
    { T::nDim } -> std::convertible_to<int>;
    { T::nOrder } -> std::convertible_to<int>;
    { T::nComponents } -> std::convertible_to<int>;
};

/**
 * @brief A scalar polynomial space (e.g., P1, P2, Pk)
 */
template <typename T>
concept ScalarPolynomialSet = PolynomialSet<T> && requires {
    requires T::nComponents == 1;
};

/**
 * @brief A vectorial polynomial space (e.g., P1^d, Nedelec, RT)
 */
template <typename T>
concept VectorialPolynomialSet = PolynomialSet<T> && requires {
    requires T::nComponents > 1;
};

//
// Basis Concepts
//

/**
 * @brief A finite element basis
 */
template <typename T>
concept Basis = requires {
    typename T::value_type;
    typename T::polyset_type;
    { T::nDof } -> std::convertible_to<int>;
    { T::nLocalDof } -> std::convertible_to<int>;
    requires PolynomialSet<typename T::polyset_type>;
};

/**
 * @brief A continuous basis (C0 continuity)
 */
template <typename T>
concept ContinuousBasis = Basis<T> && requires {
    requires T::is_continuous;
};

/**
 * @brief A discontinuous basis (DG)
 */
template <typename T>
concept DiscontinuousBasis = Basis<T> && requires {
    requires T::is_discontinuous;
};

/**
 * @brief An H(div)-conforming basis (Raviart-Thomas, BDM)
 */
template <typename T>
concept HDivBasis = VectorialPolynomialSet<typename T::polyset_type> && Basis<T>;

/**
 * @brief An H(curl)-conforming basis (Nedelec)
 */
template <typename T>
concept HCurlBasis = VectorialPolynomialSet<typename T::polyset_type> && Basis<T>;

//
// Point Set Concepts
//

/**
 * @brief A set of points (for interpolation, quadrature)
 */
template <typename T>
concept PointSet = requires(T t) {
    typename T::value_type;
    { t.nPoints() } -> std::convertible_to<int>;
};

/**
 * @brief An equidistributed point set
 */
template <typename T>
concept EquidistributedPointSet = PointSet<T>;

/**
 * @brief A Fekete point set (optimal for high-order)
 */
template <typename T>
concept FeketePointSet = PointSet<T>;

/**
 * @brief A Gauss-Lobatto point set
 */
template <typename T>
concept GaussLobattoPointSet = PointSet<T>;

//
// Quadrature Concepts
//

/**
 * @brief A quadrature rule (integration)
 */
template <typename T>
concept Quadrature = requires {
    typename T::value_type;
    typename T::node_type;
    typename T::weights_type;
    { T::Degree } -> std::convertible_to<int>;
};

/**
 * @brief A Gauss quadrature rule
 */
template <typename T>
concept GaussQuadrature = Quadrature<T>;

/**
 * @brief A Gauss-Lobatto quadrature rule
 */
template <typename T>
concept GaussLobattoQuadrature = Quadrature<T>;

//
// Polynomial Order Concepts
//

/**
 * @brief A compile-time polynomial order
 */
template <int Order>
concept PolynomialOrder = (Order >= 0);

/**
 * @brief Low-order polynomial (P0, P1, P2)
 */
template <int Order>
concept LowOrder = PolynomialOrder<Order> && (Order <= 2);

/**
 * @brief High-order polynomial (P3+)
 */
template <int Order>
concept HighOrder = PolynomialOrder<Order> && (Order > 2);

//
// Continuity Concepts
//

/**
 * @brief C0 continuous finite element
 */
template <typename T>
concept C0Continuous = requires {
    requires T::continuity == 0;
};

/**
 * @brief C1 continuous finite element
 */
template <typename T>
concept C1Continuous = requires {
    requires T::continuity == 1;
};

/**
 * @brief Discontinuous finite element
 */
template <typename T>
concept Discontinuous = requires {
    requires T::continuity == -1;
};

//
// Modal vs Nodal Basis
//

/**
 * @brief A nodal basis (values at nodes)
 */
template <typename T>
concept NodalBasis = Basis<T> && requires {
    requires T::is_nodal;
};

/**
 * @brief A modal basis (hierarchical)
 */
template <typename T>
concept ModalBasis = Basis<T> && requires {
    requires T::is_modal;
};

//
// Backward Compatibility Bridges
//

#ifdef FEELPP_ENABLE_CONCEPT_COMPATIBILITY

template <typename T>
constexpr bool is_polynomial_set_v = PolynomialSet<T>;

template <typename T>
constexpr bool is_basis_v = Basis<T>;

template <typename T>
constexpr bool is_continuous_v = ContinuousBasis<T>;

template <typename T>
constexpr bool is_discontinuous_v = DiscontinuousBasis<T>;

#endif // FEELPP_ENABLE_CONCEPT_COMPATIBILITY

} // namespace Feel

#endif /* FEELPP_FEELPOLY_CONCEPTS_HPP */
