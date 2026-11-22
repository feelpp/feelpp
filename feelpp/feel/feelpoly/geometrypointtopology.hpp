/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 * @file geometrypointtopology.hpp
 * @brief Order-independent reference-cell point restrictions.
 */
#if !defined( FEELPP_FEELPOLY_GEOMETRYPOINTTOPOLOGY_HPP )
#define FEELPP_FEELPOLY_GEOMETRYPOINTTOPOLOGY_HPP 1

#include <feel/feelmesh/hypercube.hpp>
#include <feel/feelmesh/simplex.hpp>
#include <feel/feelpoly/equispaced.hpp>
#include <feel/feelpoly/order.hpp>

#include <Eigen/Core>
#include <Eigen/LU>

#include <cmath>
#include <array>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace Feel::geometry
{
/**
 * @class GeometryPointTopology
 * @brief Reference-node restrictions from a cell to its faces and edges.
 *
 * The restriction tables are derived once from the mathematical equispaced
 * point set. This makes the reference-node layout identical for compile-time
 * and runtime orders and removes the order ceiling of legacy hand-written
 * `fToP` tables.
 *
 * @tparam Shape simplex or hypercube geometry shape
 */
template <typename Shape>
class GeometryPointTopology
{
    static_assert( Shape::nDim >= 0 && Shape::nDim <= 3,
                   "GeometryPointTopology supports dimensions zero to three" );

    static constexpr int nDim = Shape::nDim;
    using parent_shape_type = std::conditional_t<
        Shape::is_simplex, Simplex<nDim, 1, nDim>,
        Hypercube<nDim, 1, nDim>>;
    using parent_pointset_type =
        PointSetEquiSpaced<parent_shape_type, Dynamic, double>;
    using line_shape_type = Simplex<1, 1, 1>;
    using line_pointset_type =
        PointSetEquiSpaced<line_shape_type, Dynamic, double>;
    using face_shape_type = std::conditional_t<
        Shape::is_simplex, Simplex<2, 1, 2>,
        Hypercube<2, 1, 2>>;
    using face_pointset_type =
        PointSetEquiSpaced<face_shape_type, Dynamic, double>;
    using index_table_type = std::vector<std::vector<uint16_type>>;

  public:
    /**
     * @brief Return the process-wide shared topology for an effective order.
     *
     * Tables are immutable after construction and cached per shape and order.
     * Meshes therefore share the setup cost without storing duplicate
     * connectivity metadata in every mesh or geometric entity.
     *
     * @param order strictly positive geometric polynomial order
     */
    [[nodiscard]] static GeometryPointTopology const&
    get( RuntimeOrder order )
    {
        static std::mutex cacheMutex;
        static std::unordered_map<
            uint16_type, std::unique_ptr<GeometryPointTopology>>
            cache;

        std::scoped_lock lock( cacheMutex );
        auto const effectiveOrder = validatedOrder( order );
        auto [iterator, inserted] = cache.try_emplace(
            effectiveOrder.value, nullptr );
        if ( inserted )
            iterator->second =
                std::make_unique<GeometryPointTopology>(
                    effectiveOrder );
        return *iterator->second;
    }

    /**
     * @brief Construct the restriction tables for an effective order.
     * @param order strictly positive geometric polynomial order
     */
    explicit GeometryPointTopology( RuntimeOrder order )
        : M_order( validatedOrder( order ) )
    {
        if constexpr ( nDim > 0 )
            this->build();
    }

    /** @brief Return the effective geometric polynomial order. */
    [[nodiscard]] uint16_type order() const noexcept
    {
        return M_order.value;
    }

    /**
     * @brief Map one face-local point number to the cell-local point number.
     * @param face local face number
     * @param point face-local point number
     */
    [[nodiscard]] uint16_type
    faceToPoint( uint16_type face, uint16_type point ) const
    {
        return checkedEntry( M_faceToPoint, face, point, "face" );
    }

    /**
     * @brief Map one edge-local point number to the cell-local point number.
     * @param edge local edge number
     * @param point edge-local point number
     */
    [[nodiscard]] uint16_type
    edgeToPoint( uint16_type edge, uint16_type point ) const
    {
        return checkedEntry( M_edgeToPoint, edge, point, "edge" );
    }

    /** @brief Return the number of points on one cell face. */
    [[nodiscard]] size_type pointsOnFace() const noexcept
    {
        return M_faceToPoint.empty() ? 0 : M_faceToPoint.front().size();
    }

    /** @brief Return the number of points on one cell edge. */
    [[nodiscard]] size_type pointsOnEdge() const noexcept
    {
        return M_edgeToPoint.empty() ? 0 : M_edgeToPoint.front().size();
    }

  private:
    /** Validate the positive-order invariant. */
    [[nodiscard]] static RuntimeOrder
    validatedOrder( RuntimeOrder order )
    {
        if ( order.value < 1 )
            throw std::invalid_argument(
                "Geometry point topology requires an order of at least one" );
        if constexpr ( Shape::is_simplex )
            return ::Feel::detail::validatedNodalOrder<
                ::Feel::detail::SimplexCardinality>( nDim, order );
        else
            return ::Feel::detail::validatedNodalOrder<
                ::Feel::detail::HypercubeCardinality>( nDim, order );
    }

    /** Return a bounds-checked restriction entry. */
    [[nodiscard]] static uint16_type
    checkedEntry( index_table_type const& table, uint16_type entity,
                  uint16_type point, char const* entityName )
    {
        if ( entity >= table.size() || point >= table[entity].size() )
            throw std::out_of_range(
                std::string( "Invalid " ) + entityName +
                " point restriction index" );
        return table[entity][point];
    }

    /** Build the complete edge and face restriction tables. */
    void build()
    {
        parent_pointset_type parentPoints( M_order );
        this->buildEdges( parentPoints );
        this->buildFaces( parentPoints );
    }

    /** Build all cell-to-edge point restrictions. */
    void buildEdges( parent_pointset_type const& parentPoints )
    {
        line_pointset_type edgePoints( M_order );
        M_edgeToPoint.resize( parent_shape_type::numEdges );

        for ( uint16_type edge = 0;
              edge < parent_shape_type::numEdges; ++edge )
        {
            auto const first =
                parent_shape_type::e2p( edge, 0 );
            auto const second =
                parent_shape_type::e2p( edge, 1 );
            auto& restriction = M_edgeToPoint[edge];
            restriction.reserve( edgePoints.nPoints() );

            double const x0 = edgePoints.points()( 0, 0 );
            double const x1 = edgePoints.points()( 0, 1 );
            for ( size_type q = 0;
                  q < edgePoints.nPoints(); ++q )
            {
                double const t =
                    ( edgePoints.points()( 0, q ) - x0 ) / ( x1 - x0 );
                Eigen::VectorXd mapped( nDim );
                for ( int d = 0; d < nDim; ++d )
                    mapped( d ) =
                        ( 1.0 - t ) * parentPoints.points()( d, first ) +
                        t * parentPoints.points()( d, second );
                restriction.push_back(
                    this->findParentPoint( parentPoints, mapped ) );
            }
        }
    }

    /** Build all cell-to-face point restrictions. */
    void buildFaces( parent_pointset_type const& parentPoints )
    {
        if constexpr ( nDim == 1 )
        {
            M_faceToPoint.resize(
                parent_shape_type::numTopologicalFaces );
            for ( uint16_type face = 0;
                  face < parent_shape_type::numTopologicalFaces; ++face )
                M_faceToPoint[face].push_back( face );
        }
        else if constexpr ( nDim == 2 )
        {
            M_faceToPoint = M_edgeToPoint;
        }
        else if constexpr ( nDim == 3 )
        {
            face_pointset_type facePoints( M_order );
            M_faceToPoint.resize(
                parent_shape_type::numTopologicalFaces );

            constexpr int basisVertex2 =
                Shape::is_simplex ? 2 : 3;
            std::array<int, 3> const basisVertices{
                0, 1, basisVertex2 };

            Eigen::Matrix2d localBasis;
            for ( int j = 0; j < 2; ++j )
                for ( int d = 0; d < 2; ++d )
                    localBasis( d, j ) =
                        facePoints.points()(
                            d, basisVertices[j + 1] ) -
                        facePoints.points()( d, basisVertices[0] );
            Eigen::Matrix2d const localBasisInverse =
                localBasis.inverse();

            for ( uint16_type face = 0;
                  face < parent_shape_type::numTopologicalFaces; ++face )
            {
                std::array<uint16_type, 3> parentVertices{};
                for ( int j = 0; j < 3; ++j )
                    parentVertices[j] =
                        parent_shape_type::f2p(
                            face, basisVertices[j] );

                auto& restriction = M_faceToPoint[face];
                restriction.reserve( facePoints.nPoints() );
                for ( size_type q = 0;
                      q < facePoints.nPoints(); ++q )
                {
                    Eigen::Vector2d local;
                    for ( int d = 0; d < 2; ++d )
                        local( d ) =
                            facePoints.points()( d, q ) -
                            facePoints.points()(
                                d, basisVertices[0] );
                    Eigen::Vector2d const coordinates =
                        localBasisInverse * local;

                    Eigen::VectorXd mapped( nDim );
                    for ( int d = 0; d < nDim; ++d )
                    {
                        double const origin =
                            parentPoints.points()(
                                d, parentVertices[0] );
                        mapped( d ) = origin;
                        for ( int j = 0; j < 2; ++j )
                            mapped( d ) +=
                                coordinates( j ) *
                                ( parentPoints.points()(
                                      d, parentVertices[j + 1] ) -
                                  origin );
                    }
                    restriction.push_back(
                        this->findParentPoint( parentPoints, mapped ) );
                }
            }
        }
    }

    /** Match one restricted reference node to the cell point table. */
    [[nodiscard]] uint16_type
    findParentPoint( parent_pointset_type const& parentPoints,
                     Eigen::VectorXd const& point ) const
    {
        double nearestDistanceSquared =
            std::numeric_limits<double>::max();
        size_type nearest = parentPoints.nPoints();
        for ( size_type candidate = 0;
              candidate < parentPoints.nPoints(); ++candidate )
        {
            double distanceSquared = 0.0;
            for ( int d = 0; d < nDim; ++d )
            {
                double const delta =
                    parentPoints.points()( d, candidate ) - point( d );
                distanceSquared += delta * delta;
            }
            if ( distanceSquared < nearestDistanceSquared )
            {
                nearestDistanceSquared = distanceSquared;
                nearest = candidate;
            }
        }

        constexpr double tolerance = 1e-24;
        if ( nearest == parentPoints.nPoints() ||
             nearestDistanceSquared > tolerance )
            throw std::runtime_error(
                "Unable to match a restricted geometry point to the "
                "reference-cell point table" );
        if ( nearest >
             static_cast<size_type>(
                 std::numeric_limits<uint16_type>::max() ) )
            throw std::overflow_error(
                "Geometry point index exceeds the supported index range" );
        return static_cast<uint16_type>( nearest );
    }

    RuntimeOrder M_order;
    index_table_type M_faceToPoint;
    index_table_type M_edgeToPoint;
};

} // namespace Feel::geometry

#endif // FEELPP_FEELPOLY_GEOMETRYPOINTTOPOLOGY_HPP
