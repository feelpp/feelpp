/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
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
#include <cstddef>
#include <string_view>
#include <type_traits>

#include <Eigen/Core>

// clang-format off
#include <feel/feelcore/warnoff.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <feel/feelcore/warnon.hpp>
// clang-format on

#include <feel/feelpoly/traits.hpp>
#include <feel/feelpoly/order.hpp>
#include <feel/feelpoly/geometrycoefficients.hpp>

namespace Feel
{

// Dynamic is now defined in order.hpp to avoid circular dependencies

//
// Linear algebra concepts (Eigen + ublas)
//

/**
 * @brief An Eigen dense matrix or expression.
 */
template <typename T>
concept EigenMatrix = requires {
    typename std::remove_cvref_t<T>::Scalar;
    { std::remove_cvref_t<T>::RowsAtCompileTime } -> std::convertible_to<int>;
    { std::remove_cvref_t<T>::ColsAtCompileTime } -> std::convertible_to<int>;
} && std::is_base_of_v<Eigen::MatrixBase<std::remove_cvref_t<T>>, std::remove_cvref_t<T>>;

/**
 * @brief An Eigen vector (column or row).
 */
template <typename T>
concept EigenVector = EigenMatrix<T> &&
                      (std::remove_cvref_t<T>::RowsAtCompileTime == 1 ||
                       std::remove_cvref_t<T>::ColsAtCompileTime == 1);

/**
 * @brief A ublas matrix or matrix expression.
 */
template <typename T>
concept UBlasMatrix = requires(std::remove_reference_t<T> m) {
    typename std::remove_cvref_t<T>::value_type;
    { m.size1() } -> std::convertible_to<std::size_t>;
    { m.size2() } -> std::convertible_to<std::size_t>;
};

/**
 * @brief A ublas vector or vector expression.
 */
template <typename T>
concept UBlasVector = requires(std::remove_reference_t<T> v) {
    typename std::remove_cvref_t<T>::value_type;
    { v.size() } -> std::convertible_to<std::size_t>;
};

/**
 * @brief Matrix-like type (Eigen or ublas).
 */
template <typename T>
concept MatrixLike = EigenMatrix<T> || UBlasMatrix<T>;

/**
 * @brief Vector-like type (Eigen or ublas).
 */
template <typename T>
concept VectorLike = EigenVector<T> || UBlasVector<T>;

/**
 * @brief Storage that is safe for Kokkos kernels.
 */
template <typename T>
concept KokkosCompatibleStorage =
    std::is_standard_layout_v<std::remove_cvref_t<T>> &&
    std::is_trivially_copyable_v<std::remove_cvref_t<T>> &&
    std::is_trivially_destructible_v<std::remove_cvref_t<T>>;

//
// Polynomial order concepts
//

/**
 * @brief A compile-time static order tag.
 */
template <typename T>
concept StaticOrder = requires {
    { std::remove_cvref_t<T>::value } -> std::convertible_to<int>;
} && (std::remove_cvref_t<T>::value >= 0);

/**
 * @brief A compile-time dynamic order tag.
 */
template <typename T>
concept DynamicOrder = requires {
    { std::remove_cvref_t<T>::value } -> std::convertible_to<int>;
} && (std::remove_cvref_t<T>::value == Dynamic);

/**
 * @brief A type that exposes an order (static or runtime).
 */
template <typename T>
concept HasOrder = StaticOrder<T> || DynamicOrder<T> ||
                   requires {
                       { std::remove_cvref_t<T>::nOrder } -> std::convertible_to<int>;
                   } ||
                   requires(std::remove_reference_t<T> t) {
                       { t.order() } -> std::convertible_to<int>;
                   };

/**
 * @brief Detect types exposing static and dynamic-order flags.
 */
template <typename T>
concept HasOrderSupport = requires {
    { std::remove_cvref_t<T>::is_order_dynamic } -> std::convertible_to<bool>;
    { std::remove_cvref_t<T>::is_order_static } -> std::convertible_to<bool>;
};

/**
 * @brief Detect types with dynamic-order support.
 */
template <typename T>
concept HasDynamicOrder = HasOrderSupport<T> && std::remove_cvref_t<T>::is_order_dynamic;

/**
 * @brief Detect types with static-order support.
 */
template <typename T>
concept HasStaticOrder = HasOrderSupport<T> && std::remove_cvref_t<T>::is_order_static;

/**
 * @brief Detect types exposing per-entity dof topology APIs.
 */
template <typename T>
concept HasDofTopology = HasOrderSupport<T> &&
                         requires(std::remove_reference_t<T> const& t) {
                             { t.dofPerVertex() } -> std::convertible_to<uint16_type>;
                             { t.dofPerEdge() } -> std::convertible_to<uint16_type>;
                             { t.dofPerFace() } -> std::convertible_to<uint16_type>;
                             { t.dofPerVolume() } -> std::convertible_to<uint16_type>;
                             { t.localDof() } -> std::convertible_to<uint16_type>;
                             { t.dofPerEntity( uint16_type{}, uint16_type{} ) } -> std::convertible_to<uint16_type>;
                         };

//
// Convex concepts
//

/**
 * @brief A convex reference shape.
 */
template <typename T>
concept ConvexConcept = requires {
    { std::remove_cvref_t<T>::nDim } -> std::convertible_to<int>;
    { std::remove_cvref_t<T>::nOrder } -> std::convertible_to<int>;
    { std::remove_cvref_t<T>::nRealDim } -> std::convertible_to<int>;
} && is_convex<std::remove_cvref_t<T>>::value;

/**
 * @brief A simplex convex (segment/triangle/tetra).
 */
template <typename T>
concept SimplexConvex = ConvexConcept<T> && is_simplex_v<std::remove_cvref_t<T>>;

/**
 * @brief A hypercube convex (line/quad/hex).
 */
template <typename T>
concept HypercubeConvex = ConvexConcept<T> && is_hypercube_v<std::remove_cvref_t<T>>;

//
// Field concepts (scalar/vector/tensor policies)
//

/**
 * @brief Scalar field policy.
 */
template <typename T>
concept ScalarFieldConcept = std::derived_from<std::remove_cvref_t<T>, ScalarBase> ||
                             (requires { requires std::remove_cvref_t<T>::is_scalar; });

/**
 * @brief Vector field policy.
 */
template <typename T>
concept VectorFieldConcept = std::derived_from<std::remove_cvref_t<T>, VectorialBase> ||
                             (requires { requires std::remove_cvref_t<T>::is_vectorial; });

/**
 * @brief Tensor2 field policy.
 */
template <typename T>
concept Tensor2FieldConcept = std::derived_from<std::remove_cvref_t<T>, Tensor2Base> ||
                              (requires { requires std::remove_cvref_t<T>::is_tensor2; });

//
// Geometric mapping concepts
//

/**
 * @brief Polynomial representation of a geometric mapping.
 *
 * A representation defines coefficient semantics and evaluates a scalar
 * reference basis and its first two derivatives. It is independent of the
 * physical-space dimension and of consumers of the resulting map.
 *
 * @tparam T candidate representation type
 */
template <typename T>
concept GeometryRepresentation =
    requires( std::remove_cvref_t<T> const& representation,
              typename std::remove_cvref_t<T>::points_type const& points ) {
        typename std::remove_cvref_t<T>::value_type;
        typename std::remove_cvref_t<T>::points_type;
        typename std::remove_cvref_t<T>::basis_matrix_type;
        { std::remove_cvref_t<T>::coefficientSemantics() } -> std::convertible_to<geometry::CoefficientSemantics>;
        { std::remove_cvref_t<T>::dimension() } -> std::convertible_to<int>;
        { std::remove_cvref_t<T>::supportsSubdivision() } -> std::convertible_to<bool>;
        { std::remove_cvref_t<T>::supportsConvexHull() } -> std::convertible_to<bool>;
        { std::remove_cvref_t<T>::name() } -> std::convertible_to<std::string_view>;
        { representation.order() } -> std::convertible_to<uint16_type>;
        { representation.coefficientCount() } -> std::convertible_to<size_type>;
        { representation.basisValues( points ) } -> std::same_as<typename std::remove_cvref_t<T>::basis_matrix_type>;
        representation.basisGradients( points );
        representation.basisHessians( points );
    };

/**
 * @brief Geometry coefficient table bound to a representation.
 * @tparam T candidate coefficient-table type
 */
template <typename T>
concept GeometryCoefficients =
    requires( std::remove_cvref_t<T> const& coefficients ) {
        typename std::remove_cvref_t<T>::representation_type;
        typename std::remove_cvref_t<T>::value_type;
        typename std::remove_cvref_t<T>::matrix_type;
        requires GeometryRepresentation<
            typename std::remove_cvref_t<T>::representation_type>;
        { coefficients.representation() } -> std::same_as<typename std::remove_cvref_t<T>::representation_type const&>;
        { coefficients.values() } -> std::same_as<typename std::remove_cvref_t<T>::matrix_type const&>;
        { coefficients.size() } -> std::convertible_to<size_type>;
    };

/**
 * @brief A geometric mapping (GeoMap-like type).
 */
template <typename T>
concept GeometricMappingConcept =
    requires( std::remove_cvref_t<T> const& mapping ) {
        typename std::remove_cvref_t<T>::value_type;
        typename std::remove_cvref_t<T>::convex_type;
        { mapping.dim() } -> std::convertible_to<int>;
        { mapping.realDim() } -> std::convertible_to<int>;
        { mapping.order() } -> std::convertible_to<uint16_type>;
        { mapping.isLinear() } -> std::convertible_to<bool>;
    } && ConvexConcept<typename std::remove_cvref_t<T>::convex_type>;

namespace detail
{
/** Return whether a geometric mapping is statically affine. */
template <typename T>
[[nodiscard]] consteval bool isStaticallyLinearGeometricMapping()
{
    using mapping_type = std::remove_cvref_t<T>;
    if constexpr ( requires {
                       mapping_type::is_order_static;
                       mapping_type::is_linear;
                   } )
        return mapping_type::is_order_static && mapping_type::is_linear;
    else if constexpr ( requires {
                            mapping_type::isOrderStatic();
                            mapping_type::isStaticallyLinear();
                        } )
        return mapping_type::isOrderStatic() && mapping_type::isStaticallyLinear();
    else
        return false;
}

/** Return whether a geometric mapping can be nonlinear. */
template <typename T>
[[nodiscard]] consteval bool isNonLinearGeometricMapping()
{
    using mapping_type = std::remove_cvref_t<T>;
    if constexpr ( requires {
                       mapping_type::is_order_dynamic;
                       mapping_type::is_linear;
                   } )
        return mapping_type::is_order_dynamic || !mapping_type::is_linear;
    else if constexpr ( requires {
                            mapping_type::isOrderStatic();
                            mapping_type::isStaticallyLinear();
                        } )
        return !mapping_type::isOrderStatic() || !mapping_type::isStaticallyLinear();
    else
        return false;
}
} // namespace detail

/**
 * @brief A linear geometric mapping.
 */
template <typename T>
concept LinearGeometricMappingConcept =
    GeometricMappingConcept<T> && detail::isStaticallyLinearGeometricMapping<T>();

/**
 * @brief A nonlinear geometric mapping.
 */
template <typename T>
concept NonLinearGeometricMappingConcept =
    GeometricMappingConcept<T> && detail::isNonLinearGeometricMapping<T>();

/**
 * @brief A pointwise geometric-mapping evaluation context.
 *
 * The contract exposes the physical point, Jacobian matrix, inverse Jacobian
 * transpose, and integration Jacobian required by finite-element pullbacks.
 * It also requires opt-in quality diagnostics and a scale-aware validity
 * predicate so algorithms can reject singular or inverted mappings without
 * depending on a concrete context implementation.
 *
 * @tparam T geometric-mapping context type
 */
template <typename T>
concept GeometricMappingContextConcept =
    requires( std::remove_cvref_t<T> const& context, int q ) {
        typename std::remove_cvref_t<T>::gm_type;
        typename std::remove_cvref_t<T>::value_type;
        requires GeometricMappingConcept<typename std::remove_cvref_t<T>::gm_type>;
        { context.nPoints() } -> std::convertible_to<int>;
        { context.J( q ) } -> std::convertible_to<typename std::remove_cvref_t<T>::value_type>;
        { context.integrationJacobian( q ) } -> std::convertible_to<typename std::remove_cvref_t<T>::value_type>;
        context.K( q );
        context.jacobianMatrix( q );
        context.B( q );
        context.inverseJacobianTranspose( q );
        context.xReal( q );
        context.mappingQuality( q );
        { context.isValid( q ) } -> std::convertible_to<bool>;
    };

/**
 * @brief A geometric context for a cell embedded in a higher-dimensional space.
 *
 * This extends the common mapping contract with the intrinsic metric and the
 * orthogonal tangent/normal decomposition. It intentionally does not require
 * a single normal, which is not canonical beyond codimension one.
 *
 * @tparam T embedded geometric-mapping context type
 */
template <typename T>
concept EmbeddedGeometricMappingContextConcept =
    GeometricMappingContextConcept<T> &&
    ( std::remove_cvref_t<T>::PDim < std::remove_cvref_t<T>::NDim ) &&
    requires( std::remove_cvref_t<T> const& context, int q ) {
        context.metricTensor( q );
        context.inverseMetricTensor( q );
        { context.manifoldJacobian( q ) } ->
            std::convertible_to<typename std::remove_cvref_t<T>::value_type>;
        context.tangentProjector( q );
        context.normalProjector( q );
        context.normalFrame( q );
    };

//
// Polynomial space concepts
//

/**
 * @brief A polynomial basis (e.g., Dubiner, Legendre).
 */
template <typename T>
concept PolynomialBasis = requires {
    typename std::remove_cvref_t<T>::value_type;
    typename std::remove_cvref_t<T>::points_type;
    typename std::remove_cvref_t<T>::matrix_type;
    typename std::remove_cvref_t<T>::convex_type;
    { std::remove_cvref_t<T>::nDim } -> std::convertible_to<int>;
    { std::remove_cvref_t<T>::nRealDim } -> std::convertible_to<int>;
    { std::remove_cvref_t<T>::nOrder } -> std::convertible_to<int>;
} && ConvexConcept<typename std::remove_cvref_t<T>::convex_type>;

/**
 * @brief A polynomial set (collection of basis polynomials).
 */
template <typename T>
concept PolynomialSetConcept = requires {
    typename std::remove_cvref_t<T>::value_type;
    typename std::remove_cvref_t<T>::points_type;
    typename std::remove_cvref_t<T>::basis_type;
    typename std::remove_cvref_t<T>::convex_type;
    { std::remove_cvref_t<T>::nDim } -> std::convertible_to<int>;
    { std::remove_cvref_t<T>::nRealDim } -> std::convertible_to<int>;
    { std::remove_cvref_t<T>::nOrder } -> std::convertible_to<int>;
    { std::remove_cvref_t<T>::nComponents } -> std::convertible_to<int>;
} && ConvexConcept<typename std::remove_cvref_t<T>::convex_type>;

/**
 * @brief A scalar polynomial set (single component).
 */
template <typename T>
concept ScalarPolynomialSet = PolynomialSetConcept<T> && requires {
    requires std::remove_cvref_t<T>::nComponents == 1;
};

/**
 * @brief A vectorial polynomial set (multiple components).
 */
template <typename T>
concept VectorialPolynomialSet = PolynomialSetConcept<T> && requires {
    requires std::remove_cvref_t<T>::nComponents > 1;
};

//
// Basis Concepts
//

/**
 * @brief A finite element basis
 */
template <typename T>
concept BasisConcept = requires {
    typename T::value_type;
    typename T::polyset_type;
    { T::nDof } -> std::convertible_to<int>;
    { T::nLocalDof } -> std::convertible_to<int>;
    requires PolynomialSetConcept<typename T::polyset_type>;
};

/**
 * @brief A continuous basis (C0 continuity)
 */
template <typename T>
concept ContinuousBasis = BasisConcept<T> && requires {
    requires T::is_continuous;
};

/**
 * @brief A discontinuous basis (DG)
 */
template <typename T>
concept DiscontinuousBasis = BasisConcept<T> && requires {
    requires T::is_discontinuous;
};

/**
 * @brief An H(div)-conforming basis (Raviart-Thomas, BDM)
 */
template <typename T>
concept HDivBasis = VectorialPolynomialSet<typename T::polyset_type> && BasisConcept<T>;

/**
 * @brief An H(curl)-conforming basis (Nedelec)
 */
template <typename T>
concept HCurlBasis = VectorialPolynomialSet<typename T::polyset_type> && BasisConcept<T>;

//
// Point Set Concepts
//

/**
 * @brief A set of points (for interpolation, quadrature)
 */
template <typename T>
concept PointSetConcept = requires(T t) {
    typename T::value_type;
    { t.nPoints() } -> std::convertible_to<int>;
};

/**
 * @brief An equidistributed point set
 */
template <typename T>
concept EquidistributedPointSet = PointSetConcept<T>;

/**
 * @brief A Fekete point set (optimal for high-order)
 */
template <typename T>
concept FeketePointSet = PointSetConcept<T>;

/**
 * @brief A Gauss-Lobatto point set
 */
template <typename T>
concept GaussLobattoPointSet = PointSetConcept<T>;

//
// Quadrature Concepts
//

/**
 * @brief A quadrature rule (integration)
 */
template <typename T>
concept QuadratureConcept = requires {
    typename T::value_type;
    typename T::node_type;
    typename T::weights_type;
    { T::Degree } -> std::convertible_to<int>;
};

/**
 * @brief A Gauss quadrature rule
 */
template <typename T>
concept GaussQuadrature = QuadratureConcept<T>;

/**
 * @brief A Gauss-Lobatto quadrature rule
 */
template <typename T>
concept GaussLobattoQuadrature = QuadratureConcept<T>;

//
// Polynomial Order Concepts
//

/**
 * @brief A compile-time polynomial order
 */
template <int Order>
concept PolynomialOrder = (Order >= 0);

/**
 * @brief Static polynomial order (Order >= 0)
 *
 * Use with requires clauses to provide constexpr accessors:
 * @code
 * constexpr uint16_type order() const requires is_static_order<nOrder> { return nOrder; }
 * @endcode
 */
template <int Order>
concept is_static_order = (Order >= 0);

/**
 * @brief Dynamic polynomial order (Order < 0, i.e., Order == Dynamic)
 *
 * Use with requires clauses to provide runtime accessors:
 * @code
 * uint16_type order() const requires is_dynamic_order<nOrder> { return M_runtime_order; }
 * @endcode
 */
template <int Order>
concept is_dynamic_order = (Order < 0);

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
concept DiscontinuousConcept = requires {
    requires T::continuity == -1;
};

//
// Modal vs Nodal Basis
//

/**
 * @brief A nodal basis (values at nodes)
 */
template <typename T>
concept NodalBasis = BasisConcept<T> && requires {
    requires T::is_nodal;
};

/**
 * @brief A modal basis (hierarchical)
 */
template <typename T>
concept ModalBasis = BasisConcept<T> && requires {
    requires T::is_modal;
};

//
// Backward Compatibility Bridges
//

#ifdef FEELPP_ENABLE_CONCEPT_COMPATIBILITY

template <typename T>
constexpr bool is_polynomial_set_v = PolynomialSetConcept<T>;

template <typename T>
constexpr bool is_basis_v = BasisConcept<T>;

template <typename T>
constexpr bool is_continuous_v = ContinuousBasis<T>;

template <typename T>
constexpr bool is_discontinuous_v = DiscontinuousBasis<T>;

#endif // FEELPP_ENABLE_CONCEPT_COMPATIBILITY

} // namespace Feel

#endif /* FEELPP_FEELPOLY_CONCEPTS_HPP */
