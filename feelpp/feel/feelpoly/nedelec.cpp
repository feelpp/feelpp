/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2006-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 * @file nedelec.cpp
 * @brief Cached triangular-face moment transformations shared by Nedelec elements.
 */
#include <feel/feelpoly/nedelec.hpp>

#include <map>
#include <mutex>

namespace Feel
{
namespace detail
{
/**
 * @details Construct the NED1 tangent/monomial basis action on the unit triangle
 * using a full-pivot LU solve, then transpose it to act on functional values.
 * Cache access is synchronized; construction runs outside the mutex, and
 * concurrent insertions retain the first stored result for each permutation
 * and degree. The cache is owned by the core library.
 */
std::shared_ptr<DofTransform::MatrixStorage const>
nedelecTriangularFaceMomentTransform(
    std::array<uint16_type,3> const& vertexPermutation,
    uint16_type degree )
{
    using matrix_type = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
    using vector2_type = Eigen::Vector2d;

    using storage_ptrtype = std::shared_ptr<DofTransform::MatrixStorage const>;
    static std::mutex cacheMutex;
    static std::map<std::pair<uint16_type,uint16_type>, storage_ptrtype> cache;
    const uint16_type permutationCode = static_cast<uint16_type>(
        9*vertexPermutation[0]+3*vertexPermutation[1]+vertexPermutation[2] );
    auto const cacheKey = std::make_pair( permutationCode, degree );
    {
        std::scoped_lock lock( cacheMutex );
        if ( auto const it = cache.find( cacheKey ); it != cache.end() )
            return it->second;
    }
    auto const& p = vertexPermutation;

    std::array<vector2_type,3> const vertices{{ vector2_type( 0., 0. ),
                                                vector2_type( 1., 0. ),
                                                vector2_type( 0., 1. ) }};
    vector2_type const origin = vertices[p[0]];
    Eigen::Matrix2d jacobian;
    jacobian.col( 0 ) = vertices[p[1]] - origin;
    jacobian.col( 1 ) = vertices[p[2]] - origin;
    Eigen::Matrix2d const inverse = jacobian.inverse();

    uint16_type const scalarDimension = static_cast<uint16_type>( ( degree + 1 ) * ( degree + 2 ) / 2 );
    uint16_type const dimension = static_cast<uint16_type>( 2 * scalarDimension );
    matrix_type canonical = matrix_type::Zero( dimension, dimension );
    matrix_type permuted = matrix_type::Zero( dimension, dimension );

    std::vector<vector2_type> points;
    points.reserve( scalarDimension );
    if ( degree == 0 )
        points.emplace_back( 1. / 3., 1. / 3. );
    else
        for ( uint16_type i = 0; i <= degree; ++i )
            for ( uint16_type j = 0; j <= degree - i; ++j )
                points.emplace_back( double( i ) / degree, double( j ) / degree );

    auto monomial = []( vector2_type const& x, uint16_type d, uint16_type b )
    {
        return std::pow( x[0], d - b ) * std::pow( x[1], b );
    };
    for ( uint16_type q = 0; q < scalarDimension; ++q )
    {
        vector2_type const local = inverse * ( points[q] - origin );
        uint16_type scalarBasis = 0;
        for ( uint16_type d = 0; d <= degree; ++d )
            for ( uint16_type b = 0; b <= d; ++b, ++scalarBasis )
                for ( uint16_type tangent = 0; tangent < 2; ++tangent )
                {
                    uint16_type const column = static_cast<uint16_type>( tangent * scalarDimension + scalarBasis );
                    canonical( 2*q + tangent, column ) = monomial( points[q], d, b );
                    vector2_type const valueAtPoint = jacobian.col( tangent ) * monomial( local, d, b );
                    permuted( 2*q, column ) = valueAtPoint[0];
                    permuted( 2*q + 1, column ) = valueAtPoint[1];
                }
    }

    // `facePermutation()` maps local face coordinates to the face-entity
    // coordinates. Columns express the local test basis in that canonical
    // basis, hence functionals transform with the transpose action.
    matrix_type const testBasisTransform = canonical.fullPivLu().solve( permuted );
    matrix_type const functionalTransform = testBasisTransform.transpose();
    auto storage = std::make_shared<DofTransform::MatrixStorage>();
    storage->rows = dimension;
    storage->cols = dimension;
    storage->values.resize( static_cast<std::size_t>( dimension ) * dimension );
    for ( uint16_type row = 0; row < dimension; ++row )
        for ( uint16_type col = 0; col < dimension; ++col )
            storage->values[static_cast<std::size_t>( row ) * dimension + col] = functionalTransform( row, col );
    {
        std::scoped_lock lock( cacheMutex );
        return cache.try_emplace( cacheKey, std::move( storage ) ).first->second;
    }
}

/**
 * @details Tabulate RT test polynomials at reference-triangle lattice points,
 * apply the contravariant Piola map, and recover the basis action with a
 * column-pivoted QR solve. Its transpose acts on the face functional values.
 * Cache access is synchronized; construction runs outside the mutex, and
 * concurrent insertions retain the first stored result for each permutation
 * and public order. The cache is owned by the core library.
 */
std::shared_ptr<DofTransform::MatrixStorage const>
nedelecTriangularFaceRTMomentTransform(
    std::array<uint16_type,3> const& vertexPermutation,
    uint16_type publicOrder )
{
    CHECK_GT( publicOrder, 0 );
    using matrix_type = Eigen::MatrixXd;
    using vector2_type = Eigen::Vector2d;
    using storage_ptrtype = std::shared_ptr<DofTransform::MatrixStorage const>;
    static std::mutex cacheMutex;
    static std::map<std::pair<uint16_type,uint16_type>, storage_ptrtype> cache;
    const uint16_type permutationCode = static_cast<uint16_type>(
        9*vertexPermutation[0]+3*vertexPermutation[1]+vertexPermutation[2] );
    auto const cacheKey = std::make_pair( permutationCode, publicOrder );
    {
        std::scoped_lock lock( cacheMutex );
        if ( auto const it = cache.find( cacheKey ); it != cache.end() )
            return it->second;
    }

    // RaviartThomasDynamicPolynomialSet uses the Feel++ reference triangle,
    // whose vertices are (-1,-1), (1,-1), and (-1,1).  The face functional
    // and interpolation kernels convert their unit barycentric coordinates
    // to this frame before tabulating the same RT basis.
    std::array<vector2_type,3> const vertices{{ vector2_type( -1., -1. ),
                                                vector2_type(  1., -1. ),
                                                vector2_type( -1.,  1. ) }};
    vector2_type const referenceOrigin = vertices[0];
    vector2_type const origin = vertices[vertexPermutation[0]];
    Eigen::Matrix2d jacobian;
    jacobian.col( 0 ) = ( vertices[vertexPermutation[1]]-origin )/2.;
    jacobian.col( 1 ) = ( vertices[vertexPermutation[2]]-origin )/2.;
    const Eigen::Matrix2d inverse = jacobian.inverse();
    const double determinant = std::abs( jacobian.determinant() );

    RaviartThomasDynamicPolynomialSet<2> moments(
        RuntimeOrder{ static_cast<uint16_type>( publicOrder-1 ) } );
    const uint16_type latticeDegree = Feel::detail::checkedLocalIndexValue(
        static_cast<size_type>( publicOrder ) + 2 );
    const uint16_type pointCount = Feel::detail::checkedLocalIndexValue(
        Feel::detail::SimplexCardinality::total( 2, latticeDegree ) );
    typename decltype( moments )::points_type canonicalPoints( 2, pointCount );
    typename decltype( moments )::points_type localPoints( 2, pointCount );
    uint16_type point = 0;
    for ( uint16_type i = 0; i <= latticeDegree; ++i )
        for ( uint16_type j = 0; j <= latticeDegree-i; ++j, ++point )
        {
            vector2_type const x( -1.+2.*double( i )/latticeDegree,
                                  -1.+2.*double( j )/latticeDegree );
            vector2_type const local = referenceOrigin+inverse*( x-origin );
            canonicalPoints( 0, point ) = x[0];
            canonicalPoints( 1, point ) = x[1];
            localPoints( 0, point ) = local[0];
            localPoints( 1, point ) = local[1];
        }
    auto const canonicalValues = moments.evaluate( canonicalPoints );
    auto const localValues = moments.evaluate( localPoints );
    const uint16_type dimension = moments.polynomialDimension();
    matrix_type canonical( 2*pointCount, dimension );
    matrix_type permuted( 2*pointCount, dimension );
    for ( uint16_type q = 0; q < pointCount; ++q )
        for ( uint16_type basis = 0; basis < dimension; ++basis )
        {
            vector2_type localValue( localValues( 2*basis, q ),
                                     localValues( 2*basis+1, q ) );
            vector2_type mappedValue = jacobian*localValue/determinant;
            for ( uint16_type component = 0; component < 2; ++component )
            {
                canonical( 2*q+component, basis ) =
                    canonicalValues( 2*basis+component, q );
                permuted( 2*q+component, basis ) = mappedValue[component];
            }
        }
    matrix_type const testBasisTransform =
        canonical.colPivHouseholderQr().solve( permuted );
    matrix_type const functionalTransform = testBasisTransform.transpose();
    auto storage = std::make_shared<DofTransform::MatrixStorage>();
    storage->rows = dimension;
    storage->cols = dimension;
    storage->values.resize( static_cast<std::size_t>( dimension )*dimension );
    for ( uint16_type row = 0; row < dimension; ++row )
        for ( uint16_type column = 0; column < dimension; ++column )
            storage->values[static_cast<std::size_t>( row )*dimension+column] =
                functionalTransform( row, column );
    {
        std::scoped_lock lock( cacheMutex );
        return cache.try_emplace( cacheKey, std::move( storage ) ).first->second;
    }
}
} // namespace detail
} // namespace Feel
