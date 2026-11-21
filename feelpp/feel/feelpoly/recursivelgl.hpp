/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/** @file recursivelgl.hpp
 * @brief Isaac recursive LGL nodes with shared immutable simplex/entity data.
 */
#pragma once

#include <Eigen/LU>
#include <algorithm>
#include <array>
#include <feel/feelalg/eigenmap.hpp>
#include <feel/feelpoly/equispaced.hpp>
#include <feel/feelpoly/jacobi.hpp>
#include <memory>
#include <mutex>
#include <numeric>
#include <span>

namespace Feel
{

/** @brief Canonical geometric owner of one interpolation node. */
struct SimplexNodeAttachment
{
    /** Dimension of the owning entity, excluding its boundary. */
    uint16_type dimension = 0;
    /** Entity number in the Feel++ reference simplex. */
    uint16_type entity = 0;
    /** Node ordinal within its entity. */
    size_type ordinal = 0;
};

namespace detail
{
/** Largest validated degree of the recursive LGL point policy. */
inline constexpr uint16_type recursiveLGLMaxOrder = 32;

/** @brief Immutable point, integer-index and entity data for one simplex degree. */
template <int Dim, typename T>
struct RecursiveLGLData
{
    /** Barycentric integer index of a node of total degree p. */
    using MultiIndex = std::array<uint16_type, Dim + 1>;
    /** Fixed coordinate rows, independently sized node columns. */
    using Coordinates = Eigen::Matrix<T, Dim, Eigen::Dynamic>;
    /** Per-entity node lists, stored in legacy local ordering. */
    using EntityNodes = std::array<std::vector<std::vector<size_type>>, Dim + 1>;

    /** Reference coordinates in the established entity ordering. */
    Coordinates points;
    /** Barycentric coordinates in the same ordering as points. */
    Eigen::Matrix<T, Dim + 1, Eigen::Dynamic> barycentric;
    /** Lexicographically ordered barycentric multi-indices. */
    std::vector<MultiIndex> canonicalIndices;
    /** Canonical multi-index number to legacy local-node number. */
    std::vector<size_type> canonicalToLegacy;
    /** Inverse of canonicalToLegacy. */
    std::vector<size_type> legacyToCanonical;
    /** Owning entity and ordinal of each legacy local node. */
    std::vector<SimplexNodeAttachment> attachments;
    /** Interior nodes of each entity. */
    EntityNodes interiors;
    /** Entity closures in the existing local subentity ordering. */
    EntityNodes closures;
};

/** @brief Enumerate barycentric integer indices once for any simplex dimension. */
template <std::size_t Size, typename F>
void forEachSimplexMultiIndex( uint16_type degree, F&& consume )
{
    std::array<uint16_type, Size> alpha{};
    auto visit = [&]( auto&& self, std::size_t axis, uint16_type remaining ) -> void
    {
        if ( axis + 1 == Size )
        {
            alpha[axis] = remaining;
            consume( alpha );
            return;
        }
        for ( uint16_type value = 0; value <= remaining; ++value )
        {
            alpha[axis] = value;
            self( self, axis + 1, remaining - value );
        }
    };
    visit( visit, 0, degree );
}

/** @brief Build shifted LGL sets using the existing shared Jacobi root algorithm. */
template <typename T>
auto recursiveLGLFactors( uint16_type degree )
{
    std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1>> factors( degree + 1 );
    factors[0] = Eigen::Matrix<T, 1, 1>::Constant( T( .5 ) );

    // Every recursive subproblem uses a lower degree of the same one-dimensional family.
    for ( uint16_type p = 1; p <= degree; ++p )
    {
        auto& nodes = factors[p];
        nodes.resize( p + 1 );
        nodes[0] = T( 0 );
        nodes[p] = T( 1 );
        if ( p > 1 )
        {
            Eigen::Matrix<T, Eigen::Dynamic, 1> roots( p - 1 );
            dyna::Jacobi<T> polynomial( p - 1, T( 1 ), T( 1 ) );
            details::roots( polynomial, roots );
            nodes.segment( 1, p - 1 ) = ( roots.array() + T( 1 ) ) / T( 2 );
        }

        // Enforce the exact reflection symmetry required by the recursive construction.
        for ( uint16_type i = 1; i <= p / 2; ++i )
        {
            T const left = ( nodes[i] + T( 1 ) - nodes[p - i] ) / T( 2 );
            nodes[i] = left;
            nodes[p - i] = T( 1 ) - left;
        }
        if ( p % 2 == 0 ) nodes[p / 2] = T( .5 );
        for ( uint16_type i = 1; i <= p; ++i )
            if ( !std::isfinite( nodes[i] ) || nodes[i] <= nodes[i - 1] )
                throw std::runtime_error( "recursive LGL factor construction did not converge" );
    }
    return factors;
}

/** @brief Isaac's barycentric recurrence, shared by all degree and dimension policies.
 * @see https://doi.org/10.1137/20M1321802
 * @see https://tisaac.gitlab.io/recursivenodes/ (equations 1 and 2).
 */
template <std::size_t Size, typename Factors>
auto recursiveLGLPoint( std::array<uint16_type, Size> const& alpha, Factors const& factors )
{
    using T = typename Factors::value_type::Scalar;
    Eigen::Matrix<T, Size, 1> point = Eigen::Matrix<T, Size, 1>::Zero();
    auto const degree = std::accumulate( alpha.begin(), alpha.end(), uint16_type( 0 ) );

    // The degree-zero space has one barycenter, including recursive zero-degree faces.
    if ( degree == 0 ) return Eigen::Matrix<T, Size, 1>::Constant( T( 1 ) / Size ).eval();
    if constexpr ( Size == 2 )
    {
        point << factors[degree][alpha[0]], factors[degree][alpha[1]];
    }
    else
    {
        T weight = 0;
        for ( std::size_t omitted = 0; omitted < Size; ++omitted )
        {
            std::array<uint16_type, Size - 1> face{};
            for ( std::size_t i = 0, j = 0; i < Size; ++i )
                if ( i != omitted ) face[j++] = alpha[i];
            T const contribution = factors[degree][degree - alpha[omitted]];
            if ( contribution == T( 0 ) ) continue;
            auto const lower = recursiveLGLPoint( face, factors );
            for ( std::size_t i = 0, j = 0; i < Size; ++i )
                if ( i != omitted ) point[i] += contribution * lower[j++];
            weight += contribution;
        }
        point /= weight;
    }
    return point;
}

/** @brief Construct mathematical nodes while adapting the existing integer/entity layout once. */
template <int Dim, typename T>
auto makeRecursiveLGLData( uint16_type degree )
{
    using Data = RecursiveLGLData<Dim, T>;
    auto data = std::make_shared<Data>();
    PointSetEquiSpaced<Simplex<Dim>, Dynamic, T> layout{ RuntimeOrder{ degree } };
    Reference<Simplex<Dim>, Dim, 1, Dim, T> reference;
    auto const vertices = eigenMap( reference.vertices() );
    auto const count = layout.nPoints();
    data->points.resize( Dim, count );
    data->barycentric.resize( Dim + 1, count );
    data->canonicalToLegacy.resize( count );
    data->legacyToCanonical.resize( count );
    data->attachments.resize( count );
    forEachSimplexMultiIndex<Dim + 1>( degree, [&]( auto const& alpha )
                                       { data->canonicalIndices.push_back( alpha ); } );

    // Decode the established rational lattice, without inferring ownership from generated LGL coordinates.
    Eigen::Matrix<T, Dim, Dim> edges;
    for ( int d = 0; d < Dim; ++d )
        edges.col( d ) = vertices.col( d + 1 ) - vertices.col( 0 );
    auto const transform = edges.partialPivLu();
    for ( int d = 0; d <= Dim; ++d )
    {
        auto const entities = reference.entityRange( d );
        data->interiors[d].resize( entities.end() );
        data->closures[d].resize( entities.end() );
    }
    for ( size_type legacy = 0; legacy < count; ++legacy )
    {
        typename Data::MultiIndex alpha{};
        if ( degree > 0 )
        {
            Eigen::Matrix<T, Dim, 1> const local = transform.solve( eigenMap( layout.points() ).col( legacy ) - vertices.col( 0 ) );
            int first = degree;
            for ( int d = 0; d < Dim; ++d )
            {
                alpha[d + 1] = static_cast<uint16_type>( std::llround( degree * local[d] ) );
                first -= alpha[d + 1];
            }
            if ( first < 0 ) throw std::logic_error( "invalid simplex lattice index" );
            alpha[0] = first;
        }
        auto const found = std::lower_bound( data->canonicalIndices.begin(), data->canonicalIndices.end(), alpha );
        if ( found == data->canonicalIndices.end() || *found != alpha ) throw std::logic_error( "invalid simplex lattice permutation" );
        auto const canonical = std::distance( data->canonicalIndices.begin(), found );
        data->canonicalToLegacy[canonical] = legacy;
        data->legacyToCanonical[legacy] = canonical;
        auto const [dimension, entity] = layout.pointToEntity( legacy );
        auto& nodes = data->interiors[dimension][entity];
        data->attachments[legacy] = { dimension, entity, nodes.size() };
        nodes.push_back( legacy );
    }

    // Reuse existing local subentity ordering for CG/DG and orientation consumers.
    for ( int d = 0; d <= Dim; ++d )
        for ( size_type entity = 0; entity < data->interiors[d].size(); ++entity )
        {
            if ( d == 0 )
                data->closures[d][entity] = data->interiors[d][entity];
            else
            {
                auto const entities = layout.entityToLocal( d, entity, true );
                for ( int lower = 0; lower <= d; ++lower )
                    for ( auto id : entities[lower] )
                    {
                        auto const& nodes = data->interiors[lower][id];
                        data->closures[d][entity].insert( data->closures[d][entity].end(), nodes.begin(), nodes.end() );
                    }
            }
        }

    // Apply the same recursive rule to every canonical index, then map once to reference coordinates.
    auto const factors = recursiveLGLFactors<T>( degree );
    for ( size_type canonical = 0; canonical < count; ++canonical )
    {
        auto const legacy = data->canonicalToLegacy[canonical];
        auto const point = recursiveLGLPoint( data->canonicalIndices[canonical], factors );
        data->barycentric.col( legacy ) = point;
        data->points.col( legacy ).noalias() = vertices * point;
    }
    return std::shared_ptr<Data const>( std::move( data ) );
}

/** @brief Bounded weak reuse; live point sets retain ownership, never mutable shared scratch. */
template <int Dim, typename T>
auto recursiveLGLData( uint16_type degree )
{
    static std::mutex mutex;
    static std::array<std::weak_ptr<RecursiveLGLData<Dim, T> const>, recursiveLGLMaxOrder + 1> cache;
    if ( degree > recursiveLGLMaxOrder ) throw std::out_of_range( "recursive LGL degrees are supported in [0,32]" );
    std::lock_guard lock( mutex );
    auto data = cache[degree].lock();
    if ( !data )
    {
        data = makeRecursiveLGLData<Dim, T>( degree );
        cache[degree] = data;
    }
    return data;
}
} // namespace detail

/** @brief Recursive LGL interpolation nodes on a simplex, for static or runtime degree.
 * @tparam Convex Simplex of dimension one to three; coordinates are reference coordinates.
 * @tparam Order Compile-time degree or Dynamic, independently of geometric order.
 * @tparam T Floating-point coordinate scalar.
 * @note Full nodes have recursively matching traces. Interior-only selections are DG data.
 * @note Eigen data are canonical owners; inherited coordinates are an explicit legacy copy.
 */
template <class Convex, int Order, typename T = double>
class PointSetRecursiveLGL : public PointSet<Convex, T>, public OrderBase<Order>
{
    static_assert( Convex::is_simplex && Convex::nDim >= 1 && Convex::nDim <= 3 );
    static_assert( std::is_floating_point_v<T>, "recursive LGL coordinates require a floating-point scalar" );
    static_assert( Order == Dynamic || ( Order >= 0 && Order <= ::Feel::detail::recursiveLGLMaxOrder ) );
    using Base = PointSet<Convex, T>;

  public:
    /** Spatial dimension independent of polynomial order. */
    static constexpr int Dim = Convex::nDim;
    /** Immutable node-data representation, shared across static/runtime views. */
    using Data = ::Feel::detail::RecursiveLGLData<Dim, T>;
    /** Legacy point-matrix return type retained at FE compatibility boundaries. */
    using points_type = typename Base::nodes_type;
    /** Coordinate scalar type. */
    using value_type = T;
    /** Full-node compile-time extent, independent of a runtime interior selection. */
    static constexpr int pointExtent = Order == Dynamic ? Dynamic : ::Feel::detail::SimplexCardinality::total( Dim, Order );
    using OrderBase<Order>::order;
    using Base::nPoints;

    /** @brief Construct static-order nodes, optionally restricting to the cell interior. */
    explicit PointSetRecursiveLGL( int interior = 0 )
        requires( Order != Dynamic )
        : PointSetRecursiveLGL( RuntimeOrder{ Order }, interior )
    {
    }

    /** @brief Construct checked static/runtime nodes and retain immutable reference data. */
    explicit PointSetRecursiveLGL( RuntimeOrder degree, int interior = 0 )
        : Base( 0, Dim ), OrderBase<Order>( degree ), M_data( ::Feel::detail::recursiveLGLData<Dim, T>( degree.value ) )
    {
        if ( interior != 0 && interior != 1 ) throw std::invalid_argument( "recursive LGL interior selector must be 0 or 1" );
        M_begin = interior && !M_data->interiors[Dim][0].empty() ? M_data->interiors[Dim][0].front() : 0;
        M_count = interior ? M_data->interiors[Dim][0].size() : M_data->points.cols();
        points_type legacy( Dim, M_count );
        eigenMap( legacy ) = pointsEigen();
        this->setPoints( legacy );
        this->setName( "recursive-lgl", degree.value );
    }

    /** @brief Retain an immutable owner independently of this object's lifetime. */
    [[nodiscard]] auto referenceData() const noexcept { return M_data; }
    /** @brief Borrow selected reference points until this point set is destroyed. */
    [[nodiscard]] auto pointsEigen() const& { return M_data->points.middleCols( M_begin, M_count ); }
    /** Prevent a borrowed view into a destroyed temporary point set. */
    auto pointsEigen() const&& = delete;
    /** @brief Borrow the complete node matrix with a static column extent when available. */
    [[nodiscard]] auto fullPointsEigen() const&
    {
        return Eigen::Map<Eigen::Matrix<T, Dim, pointExtent> const>( M_data->points.data(), Dim, M_data->points.cols() );
    }
    /** Prevent a full-node view into a temporary. */
    auto fullPointsEigen() const&& = delete;
    /** @brief Borrow selected barycentric coordinates in legacy local ordering. */
    [[nodiscard]] auto barycentricPoints() const& { return M_data->barycentric.middleCols( M_begin, M_count ); }
    /** Prevent a barycentric view into a temporary. */
    auto barycentricPoints() const&& = delete;
    /** @brief Entity owner of a selected node. */
    [[nodiscard]] SimplexNodeAttachment attachment( size_type point ) const
    {
        if ( point >= M_count ) throw std::out_of_range( "recursive LGL node index" );
        return M_data->attachments[M_begin + point];
    }
    /** @brief Canonical-to-legacy permutation of the full point set. */
    [[nodiscard]] std::span<size_type const> canonicalToLegacyPermutation() const& { return M_data->canonicalToLegacy; }
    /** Prevent a permutation view into a temporary. */
    auto canonicalToLegacyPermutation() const&& = delete;
    /** @brief Reference coordinates on an entity, optionally including its boundary. */
    [[nodiscard]] points_type pointsBySubEntity( uint16_type dimension, uint16_type entity, bool boundary = false ) const
    {
        if ( dimension > Dim || entity >= M_data->interiors[dimension].size() ) throw std::out_of_range( "recursive LGL entity index" );
        auto const& indices = boundary ? M_data->closures[dimension][entity] : M_data->interiors[dimension][entity];
        size_type count = 0;
        for ( auto index : indices )
            if ( index >= M_begin && index < M_begin + M_count ) ++count;
        points_type result( Dim, count );
        for ( size_type q = 0; auto index : indices )
            if ( index >= M_begin && index < M_begin + M_count ) eigenMap( result ).col( q++ ) = M_data->points.col( index );
        return result;
    }
    /** @brief Local-node permutation induced by a bijection of simplex vertices. */
    [[nodiscard]] std::vector<size_type> nodePermutation( std::array<uint16_type, Dim + 1> const& vertices ) const
    {
        auto sorted = vertices;
        std::sort( sorted.begin(), sorted.end() );
        for ( int i = 0; i <= Dim; ++i )
            if ( sorted[i] != i ) throw std::invalid_argument( "simplex vertex permutation is not a bijection" );
        std::vector<size_type> result( M_count );
        for ( size_type i = 0; i < M_count; ++i )
        {
            auto const& alpha = M_data->canonicalIndices[M_data->legacyToCanonical[M_begin + i]];
            typename Data::MultiIndex mapped{};
            for ( int v = 0; v <= Dim; ++v )
                mapped[vertices[v]] = alpha[v];
            auto const found = std::lower_bound( M_data->canonicalIndices.begin(), M_data->canonicalIndices.end(), mapped );
            auto const target = M_data->canonicalToLegacy.at( std::distance( M_data->canonicalIndices.begin(), found ) );
            if ( target < M_begin || target >= M_begin + M_count ) throw std::logic_error( "simplex permutation changed the interior selection" );
            result[i] = target - M_begin;
        }
        return result;
    }

  private:
    /** Owner shared by static/runtime point sets at the same mathematical configuration. */
    std::shared_ptr<Data const> M_data;
    /** First selected node in the immutable full-node ordering. */
    size_type M_begin = 0;
    /** Number of nodes selected by the full/interior construction. */
    size_type M_count = 0;
};
} // namespace Feel
