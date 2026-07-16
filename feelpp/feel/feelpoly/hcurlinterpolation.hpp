/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

/** @file hcurlinterpolation.hpp @brief Reusable H(curl) interpolation kernels. */

#ifndef FEELPP_FEELPOLY_HCURLINTERPOLATION_HPP
#define FEELPP_FEELPOLY_HCURLINTERPOLATION_HPP 1

#include <type_traits>
#include <cmath>

#include <Eigen/Core>

#include <feel/feelcore/feeltypes.hpp>
#include <feel/feelpoly/pointsetquadrature.hpp>
#include <feel/feelpoly/hcurlmomentpolicy.hpp>

namespace Feel::detail
{

/**
 * Apply edge-tangential point DoFs to an expression.
 *
 * Geometry/tangent acquisition is a policy supplied by the caller, keeping
 * this local kernel common to element and facet interpolation and suitable for
 * static or runtime cardinalities.
 *
 * @tparam Expr expression evaluated at functional support points
 * @tparam LocalVector mutable local coefficient vector
 * @tparam Tangent tangent vector storage
 * @tparam TangentProvider callable providing the physical edge tangent
 * @param expression expression to interpolate
 * @param localValues destination local coefficients
 * @param edgeCount number of edges in the traversed entity
 * @param dofPerEdge number of ordered point DoFs per edge
 * @param tangent reusable tangent workspace
 * @param tangentProvider geometry policy called for every edge point
 */
template<typename Expr, typename LocalVector, typename Tangent, typename TangentProvider>
void
interpolateHCurlEdgePointDofs( Expr const& expression,
                              LocalVector& localValues,
                              uint16_type edgeCount,
                              uint16_type dofPerEdge,
                              Tangent& tangent,
                              TangentProvider&& tangentProvider )
{
    using shape_type = typename std::remove_cvref_t<Expr>::shape;
    for ( uint16_type edge = 0; edge < edgeCount; ++edge )
    {
        for ( uint16_type ordinal = 0; ordinal < dofPerEdge; ++ordinal )
        {
            const uint16_type point = static_cast<uint16_type>( edge*dofPerEdge + ordinal );
            tangentProvider( edge, point, tangent );
            for ( uint16_type component = 0; component < shape_type::M; ++component )
                localValues( point ) += expression.evalq( component, 0, point ) * tangent( component );
        }
    }
}

/**
 * Evaluate the reference-cell integrand of an H(curl) interior moment.
 *
 * For the covariant Piola map `v = K^{-T} v_hat`, the inverse pullback is
 * `v_hat = K^T v`. Hence `v_hat . q_hat = v . (K q_hat)`. Reference-only
 * test expressions have no geometric context and use the identity branch.
 *
 * @tparam Expr physical or reference expression type
 * @tparam PolynomialValues tabulated reference vector-polynomial values
 * @param expression expression being interpolated
 * @param polynomialValues reference test-polynomial values
 * @param polynomialIndex selected vector polynomial
 * @param quadraturePoint quadrature point in the polynomial table
 * @param expressionPoint corresponding expression evaluation point
 * @param componentCount vector dimension
 * @return covariant-Piola-correct reference moment integrand
 */
template<typename Expr, typename PolynomialValues>
auto
hcurlInteriorMomentIntegrand( Expr const& expression,
                              PolynomialValues const& polynomialValues,
                              int polynomialIndex,
                              int quadraturePoint,
                              int expressionPoint,
                              int componentCount )
{
    using value_type = std::decay_t<decltype( expression.evalq( 0, 0, expressionPoint ) )>;
    value_type result = 0;

    if constexpr ( requires { expression.geom()->K( expressionPoint ); } )
    {
        auto const& jacobian = expression.geom()->K( expressionPoint );
        for ( int physicalComponent = 0; physicalComponent < componentCount; ++physicalComponent )
        {
            value_type transformedTest = 0;
            for ( int referenceComponent = 0; referenceComponent < componentCount; ++referenceComponent )
                transformedTest += jacobian( physicalComponent, referenceComponent ) *
                                   polynomialValues( componentCount * polynomialIndex + referenceComponent,
                                                     quadraturePoint );
            result += expression.evalq( physicalComponent, 0, expressionPoint ) * transformedTest;
        }
    }
    else
    {
        for ( int component = 0; component < componentCount; ++component )
            result += expression.evalq( component, 0, expressionPoint ) *
                      polynomialValues( componentCount * polynomialIndex + component, quadraturePoint );
    }
    return result;
}

/**
 * Apply one triangular-face block of tangential NED1 moments.
 *
 * The reference moment is evaluated through the inverse covariant pullback
 * identity `u_hat.t_hat = u.(K t_hat)`. The functional ordering is identical
 * to `appendHCurlTriangularFaceMoments`: tangent, total degree, exponent of the
 * second affine face coordinate.
 *
 * @tparam Expr expression evaluated at the face quadrature samples
 * @tparam LocalVector mutable local interpolation vector
 * @tparam ReferenceConvexType reference tetrahedron type
 * @param expression expression being interpolated
 * @param localValues destination local coefficients
 * @param referenceConvex reference tetrahedron
 * @param face reference face identifier
 * @param polynomialDegree scalar face test-space degree
 * @param firstDof first destination face-moment DoF
 * @param expressionPointOffset first face sample in the expression point table
 * @param expectedMomentCount expected number of face moments
 */
template<typename Expr, typename LocalVector, typename ReferenceConvexType>
void
interpolateHCurlTriangularFaceMoments(
    Expr const& expression,
    LocalVector& localValues,
    ReferenceConvexType const& referenceConvex,
    uint16_type face,
    int polynomialDegree,
    uint16_type firstDof,
    uint32_type expressionPointOffset,
    uint16_type expectedMomentCount )
{
    if ( polynomialDegree < 0 )
        return;
    using value_type = typename ReferenceConvexType::value_type;
    IMGeneral<3, value_type, Simplex> quadrature(
        hcurlTriangularFaceMomentQuadratureOrder( polynomialDegree ) );
    const uint32_type pointCount = quadrature.nPointsOnFace( face );
    if constexpr ( requires { expression.nPoints(); } )
    {
        CHECK( expression.nPoints() >= expressionPointOffset+pointCount )
            << "H(curl) face-moment interpolation requires expression values at appended face quadrature points: got "
            << expression.nPoints() << " points, need " << expressionPointOffset+pointCount;
    }

    auto const vertices = referenceConvex.faceVertices( face );
    Eigen::Matrix<value_type,3,1> origin, tangent0, tangent1;
    for ( int component = 0; component < 3; ++component )
    {
        origin( component ) = vertices( component, 0 );
        tangent0( component ) = vertices( component, 1 )-vertices( component, 0 );
        tangent1( component ) = vertices( component, 2 )-vertices( component, 0 );
    }
    const value_type g00 = tangent0.dot( tangent0 );
    const value_type g01 = tangent0.dot( tangent1 );
    const value_type g11 = tangent1.dot( tangent1 );
    const value_type metricDeterminant = g00*g11-g01*g01;
    const value_type inverseDeterminant = value_type( 1 )/metricDeterminant;
    const value_type surfaceJacobian = std::sqrt( metricDeterminant );
    auto const& points = quadrature.fpoints( face );

    uint16_type moment = 0;
    for ( auto const& referenceTangent : { tangent0, tangent1 } )
    {
        for ( int degree = 0; degree <= polynomialDegree; ++degree )
        {
            for ( int exponentB = 0; exponentB <= degree; ++exponentB, ++moment )
            {
                const int exponentA = degree-exponentB;
                const uint16_type dof = static_cast<uint16_type>( firstDof+moment );
                for ( uint32_type q = 0; q < pointCount; ++q )
                {
                    Eigen::Matrix<value_type,3,1> delta;
                    for ( int component = 0; component < 3; ++component )
                        delta( component ) = points( component, q )-origin( component );
                    const value_type rhs0 = tangent0.dot( delta );
                    const value_type rhs1 = tangent1.dot( delta );
                    const value_type a = ( g11*rhs0-g01*rhs1 )*inverseDeterminant;
                    const value_type b = ( g00*rhs1-g01*rhs0 )*inverseDeterminant;
                    const uint32_type expressionPoint = expressionPointOffset+q;

                    Eigen::Matrix<value_type,3,1> testTangent = referenceTangent;
                    if constexpr ( requires { expression.geom()->K( expressionPoint ); } )
                    {
                        auto const& jacobian = expression.geom()->K( expressionPoint );
                        testTangent.setZero();
                        for ( int physical = 0; physical < 3; ++physical )
                            for ( int reference = 0; reference < 3; ++reference )
                                testTangent( physical ) +=
                                    jacobian( physical, reference )*referenceTangent( reference );
                    }
                    value_type tangentialValue = 0;
                    for ( int component = 0; component < 3; ++component )
                        tangentialValue += expression.evalq(
                            component, 0, expressionPoint )*testTangent( component );
                    localValues( dof ) += quadrature.weight( face, q )/surfaceJacobian*
                        std::pow( a, exponentA )*std::pow( b, exponentB )*tangentialValue;
                }
            }
        }
    }
    CHECK_EQ( moment, expectedMomentCount )
        << "Invalid H(curl) triangular-face moment count";
}

/**
 * Apply one NED2 triangular-face block of RT tangential moments.
 *
 * This is the interpolation counterpart of
 * `appendHCurlTriangularFaceRTMoments`; both use the same intrinsic face RT
 * basis, quadrature, ordering, and covariant-Piola inverse pullback.
 */
template<typename Expr, typename LocalVector, typename ReferenceConvexType,
         typename FaceMomentSpace>
void
interpolateHCurlTriangularFaceRTMoments(
    Expr const& expression,
    LocalVector& localValues,
    ReferenceConvexType const& referenceConvex,
    FaceMomentSpace const& moments,
    uint16_type face,
    uint16_type firstDof,
    uint32_type expressionPointOffset )
{
    using value_type = typename ReferenceConvexType::value_type;
    const uint16_type quadratureOrder = static_cast<uint16_type>(
        2*( moments.publicOrder()+4 ) );
    IMGeneral<3, value_type, Simplex> quadrature( quadratureOrder );
    const uint32_type pointCount = quadrature.nPointsOnFace( face );
    if constexpr ( requires { expression.nPoints(); } )
        CHECK( expression.nPoints() >= expressionPointOffset+pointCount );

    auto const vertices = referenceConvex.faceVertices( face );
    Eigen::Matrix<value_type,3,1> origin, tangent0, tangent1;
    for ( int component = 0; component < 3; ++component )
    {
        origin( component ) = vertices( component, 0 );
        tangent0( component ) = vertices( component, 1 )-vertices( component, 0 );
        tangent1( component ) = vertices( component, 2 )-vertices( component, 0 );
    }
    const value_type g00 = tangent0.dot( tangent0 );
    const value_type g01 = tangent0.dot( tangent1 );
    const value_type g11 = tangent1.dot( tangent1 );
    const value_type metricDeterminant = g00*g11-g01*g01;
    const value_type inverseDeterminant = value_type( 1 )/metricDeterminant;
    const value_type surfaceJacobian = std::sqrt( metricDeterminant );
    auto const& points = quadrature.fpoints( face );
    typename FaceMomentSpace::points_type intrinsicPoints( 2, pointCount );
    for ( uint16_type q = 0; q < pointCount; ++q )
    {
        Eigen::Matrix<value_type,3,1> delta;
        for ( int component = 0; component < 3; ++component )
            delta( component ) = points( component, q )-origin( component );
        const value_type rhs0 = tangent0.dot( delta );
        const value_type rhs1 = tangent1.dot( delta );
        const value_type a = ( g11*rhs0-g01*rhs1 )*inverseDeterminant;
        const value_type b = ( g00*rhs1-g01*rhs0 )*inverseDeterminant;
        // Keep interpolation in the same Feel++ triangle coordinates as the
        // dual RT moment builder.
        intrinsicPoints( 0, q ) = value_type( -1 )+value_type( 2 )*a;
        intrinsicPoints( 1, q ) = value_type( -1 )+value_type( 2 )*b;
    }
    auto const momentValues = moments.evaluate( intrinsicPoints );

    for ( uint16_type moment = 0; moment < moments.polynomialDimension(); ++moment )
    {
        const uint16_type dof = static_cast<uint16_type>( firstDof+moment );
        for ( uint16_type q = 0; q < pointCount; ++q )
        {
            Eigen::Matrix<value_type,3,1> referenceMoment =
                tangent0*momentValues( 2*moment, q )+
                tangent1*momentValues( 2*moment+1, q );
            Eigen::Matrix<value_type,3,1> physicalMoment = referenceMoment;
            const uint32_type expressionPoint = expressionPointOffset+q;
            if constexpr ( requires { expression.geom()->K( expressionPoint ); } )
            {
                auto const& jacobian = expression.geom()->K( expressionPoint );
                physicalMoment.setZero();
                for ( int physical = 0; physical < 3; ++physical )
                    for ( int reference = 0; reference < 3; ++reference )
                        physicalMoment( physical ) +=
                            jacobian( physical, reference )*referenceMoment( reference );
            }
            value_type value = 0;
            for ( int component = 0; component < 3; ++component )
                value += expression.evalq( component, 0, expressionPoint )*
                         physicalMoment( component );
            localValues( dof ) += quadrature.weight( face, q )/surfaceJacobian*value;
        }
    }
}

/**
 * Apply vector cell-moment DoFs to an expression.
 *
 * The caller supplies the already selected vector polynomial test space, so
 * static and runtime NED1 differ only in how that space stores its order.
 * Quadrature selection, covariant-Piola handling, point offsets, and moment
 * accumulation are shared.
 *
 * @tparam Dim reference-cell dimension
 * @tparam Expr expression evaluated at the FE support/quadrature points
 * @tparam LocalVector mutable local interpolation vector
 * @tparam MomentPolynomialSet vector polynomial moment space
 * @param expression expression being interpolated
 * @param localValues destination local coefficients
 * @param moments vector polynomial test space
 * @param internalOrder ambient polynomial degree
 * @param firstDof first cell-interior local DoF
 * @param expressionPointOffset first cell quadrature expression sample
 * @param expectedMomentCount expected number of cell moments
 * @param componentCount vector dimension
 */
template<uint16_type Dim, typename Expr, typename LocalVector, typename MomentPolynomialSet>
void
interpolateHCurlCellMoments( Expr const& expression,
                             LocalVector& localValues,
                             MomentPolynomialSet const& moments,
                             uint16_type internalOrder,
                             uint16_type firstDof,
                             uint32_type expressionPointOffset,
                             uint16_type expectedMomentCount,
                             uint16_type componentCount )
{
    using value_type = typename MomentPolynomialSet::value_type;
    CHECK_EQ( moments.polynomialDimension(), expectedMomentCount )
        << "Invalid H(curl) cell-moment dimension";
    IMGeneral<Dim, value_type, Simplex> quadrature(
        static_cast<uint16_type>( 2*( internalOrder+1 ) ) );
    auto const momentValues = moments.evaluate( quadrature.points() );
    if constexpr ( requires { expression.nPoints(); } )
    {
        CHECK( expression.nPoints() >= expressionPointOffset+quadrature.nPoints() )
            << "H(curl) cell-moment interpolation requires expression values at appended quadrature points: got "
            << expression.nPoints() << " points, need " << expressionPointOffset+quadrature.nPoints();
    }
    for ( uint16_type moment = 0; moment < expectedMomentCount; ++moment )
    {
        const uint16_type dof = static_cast<uint16_type>( firstDof+moment );
        for ( uint16_type q = 0; q < quadrature.nPoints(); ++q )
        {
            const uint32_type expressionPoint = expressionPointOffset+q;
            localValues( dof ) += quadrature.weight( q )*
                hcurlInteriorMomentIntegrand(
                    expression, momentValues, moment, q, expressionPoint, componentCount );
        }
    }
}

} // namespace Feel::detail

#endif // FEELPP_FEELPOLY_HCURLINTERPOLATION_HPP
