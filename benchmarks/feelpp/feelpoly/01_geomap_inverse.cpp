/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file 01_geomap_inverse.cpp
 * @brief Microbenchmarks for full-dimensional inverse maps and embedded closest points.
 */

#include <benchmark/benchmark.h>
#include <feel/feelpoly/geomap.hpp>

#include <algorithm>
#include <cmath>
#include <memory>
#include <type_traits>

using namespace Feel;

namespace
{
/**
 * @struct inverse_benchmark_geometry
 * @brief Minimal owner of physical geometric coefficients for an inverse map.
 * @tparam GeomapType geometric mapping type
 */
template <typename GeomapType>
struct inverse_benchmark_geometry
{
    using matrix_type = typename GeomapType::matrix_node_t_type;

    /** Physical interpolation coefficients stored by columns. */
    matrix_type coefficients;

    /** @brief Return all physical geometric coefficients. */
    [[nodiscard]] matrix_type const& G() const noexcept { return coefficients; }

    /** @brief Return the affine compatibility coefficient view. */
    [[nodiscard]] matrix_type const& vertices() const noexcept
    {
        return coefficients;
    }
};

/**
 * @brief Select a static- or runtime-order geometric mapping.
 * @tparam Dim reference dimension
 * @tparam Order effective geometric polynomial order
 * @tparam RealDim ambient physical dimension
 * @tparam RuntimeOrderV whether the order is selected at runtime
 */
template <int Dim, int Order, int RealDim, bool RuntimeOrderV>
using inverse_benchmark_geomap_type = std::conditional_t<
    RuntimeOrderV,
    GeoMap<Dim, Dynamic, RealDim, double, Simplex>,
    GeoMap<Dim, Order, RealDim, double, Simplex>>;

/**
 * @brief Construct the selected static- or runtime-order geometric mapping.
 * @tparam Dim reference dimension
 * @tparam Order effective geometric polynomial order
 * @tparam RealDim ambient physical dimension
 * @tparam RuntimeOrderV whether the order is selected at runtime
 */
template <int Dim, int Order, int RealDim, bool RuntimeOrderV>
[[nodiscard]] auto makeInverseBenchmarkGeomap()
{
    using geomap_type =
        inverse_benchmark_geomap_type<Dim, Order, RealDim, RuntimeOrderV>;
    if constexpr ( RuntimeOrderV )
        return std::make_shared<geomap_type>( RuntimeOrder{ Order } );
    else
        return std::make_shared<geomap_type>();
}

/**
 * @brief Evaluate a polynomial benchmark mapping exactly representable at `Order`.
 *
 * Supported configurations are a 2D cell, a 3D cell, and a two-dimensional
 * surface embedded in three dimensions.
 *
 * @tparam Dim reference dimension
 * @tparam Order geometric polynomial order
 * @tparam RealDim ambient physical dimension
 * @param point reference coordinate
 * @return physical coordinate
 */
template <int Dim, int Order, int RealDim>
[[nodiscard]] Eigen::Matrix<double, RealDim, 1>
inverseBenchmarkMap( Eigen::Matrix<double, Dim, 1> const& point )
{
    static_assert( ( Dim == 2 && RealDim == 2 ) ||
                       ( Dim == 3 && RealDim == 3 ) ||
                       ( Dim == 2 && RealDim == 3 ),
                   "unsupported inverse benchmark dimensions" );

    Eigen::Matrix<double, RealDim, 1> result;
    if constexpr ( Dim == 2 && RealDim == 2 )
    {
        double const x = point( 0 );
        double const y = point( 1 );
        result << 0.15 + 1.08 * x + 0.09 * y,
            -0.2 + 0.05 * x + 1.04 * y;
        if constexpr ( Order >= 2 )
        {
            result( 0 ) += 0.06 * x * y + 0.025 * y * y;
            result( 1 ) += 0.045 * x * x - 0.02 * x * y;
        }
        if constexpr ( Order >= 3 )
        {
            result( 0 ) += 0.012 * x * x * x - 0.008 * x * y * y;
            result( 1 ) += 0.01 * y * y * y + 0.006 * x * x * y;
        }
    }
    else if constexpr ( Dim == 3 && RealDim == 3 )
    {
        double const x = point( 0 );
        double const y = point( 1 );
        double const z = point( 2 );
        result << 0.1 + 1.03 * x + 0.05 * y + 0.025 * z,
            -0.16 + 0.035 * x + 1.02 * y + 0.04 * z,
            0.08 + 0.025 * x + 0.03 * y + 1.05 * z;
        if constexpr ( Order >= 2 )
        {
            result( 0 ) += 0.025 * x * y + 0.018 * z * z;
            result( 1 ) += 0.022 * y * z + 0.014 * x * x;
            result( 2 ) += 0.02 * x * z - 0.012 * y * y;
        }
        if constexpr ( Order >= 3 )
        {
            result( 0 ) += 0.005 * x * y * z;
            result( 1 ) += 0.004 * x * x * z;
            result( 2 ) += 0.0045 * x * y * y;
        }
    }
    else
    {
        double const x = point( 0 );
        double const y = point( 1 );
        result << 0.1 + 1.02 * x + 0.04 * y,
            -0.12 + 0.03 * x + 1.01 * y,
            0.18;
        if constexpr ( Order >= 2 )
        {
            result( 0 ) += 0.035 * x * y;
            result( 1 ) += 0.025 * x * x;
            result( 2 ) += 0.2 * x * x + 0.17 * y * y + 0.045 * x * y;
        }
        if constexpr ( Order >= 3 )
            result( 2 ) += 0.018 * x * x * y - 0.012 * x * y * y;
    }
    return result;
}

/**
 * @brief Sample physical interpolation coefficients at the GeoMap nodes.
 * @tparam Dim reference dimension
 * @tparam Order effective geometric order
 * @tparam RealDim ambient physical dimension
 * @tparam GeomapType geometric mapping type
 * @param geomap geometric mapping
 * @return benchmark geometry coefficients
 */
template <int Dim, int Order, int RealDim, typename GeomapType>
[[nodiscard]] inverse_benchmark_geometry<GeomapType>
makeInverseBenchmarkGeometry( std::shared_ptr<GeomapType> const& geomap )
{
    inverse_benchmark_geometry<GeomapType> geometry{
        typename GeomapType::matrix_node_t_type(
            RealDim, geomap->nbPoints() ) };
    for ( uint16_type q = 0; q < geomap->nbPoints(); ++q )
    {
        Eigen::Matrix<double, Dim, 1> referencePoint;
        for ( int component = 0; component < Dim; ++component )
            referencePoint( component ) = geomap->refNode( q )( component );
        auto const physicalPoint =
            inverseBenchmarkMap<Dim, Order, RealDim>( referencePoint );
        for ( int component = 0; component < RealDim; ++component )
            geometry.coefficients( component, q ) = physicalPoint( component );
    }
    return geometry;
}

/**
 * @brief Return an interior reference point with a controllable batch offset.
 * @tparam Dim reference dimension
 * @param parameter value in `[-1,1]` used to vary neighboring batch points
 * @return admissible simplex reference coordinate
 */
template <int Dim>
[[nodiscard]] Eigen::Matrix<double, Dim, 1>
inverseBenchmarkReferencePoint( double parameter = 0.0 )
{
    Eigen::Matrix<double, Dim, 1> result;
    if constexpr ( Dim == 2 )
        result << -0.58 + 0.2 * parameter,
            -0.52 + 0.07 * parameter;
    else
        result << -0.62 + 0.08 * parameter,
            -0.57 + 0.04 * parameter,
            -0.52 - 0.02 * parameter;
    return result;
}

/**
 * @brief Return the simplex reference barycenter used as a cold initial guess.
 * @tparam Dim reference dimension
 */
template <int Dim>
[[nodiscard]] Eigen::Matrix<double, Dim, 1>
inverseBenchmarkInitialGuess()
{
    return Eigen::Matrix<double, Dim, 1>::Constant(
        double( 1 - Dim ) / double( Dim + 1 ) );
}

/**
 * @brief Convert a fixed Eigen vector to a Feel++ node.
 * @tparam VectorType fixed Eigen vector type
 * @param point Eigen coordinate
 * @return Feel++ node
 */
template <typename VectorType>
[[nodiscard]] node_t<double> inverseBenchmarkNode( VectorType const& point )
{
    node_t<double> result( point.rows() );
    for ( Eigen::Index component = 0; component < point.rows(); ++component )
        result( component ) = point( component );
    return result;
}

/**
 * @brief Return a physical exact-inverse or embedded closest-point target.
 * @tparam Dim reference dimension
 * @tparam Order geometric polynomial order
 * @tparam RealDim ambient physical dimension
 * @tparam Projected whether to offset an embedded target in the normal direction
 * @param referencePoint reference coordinate
 * @return physical target
 */
template <int Dim, int Order, int RealDim, bool Projected>
[[nodiscard]] Eigen::Matrix<double, RealDim, 1>
inverseBenchmarkTarget(
    Eigen::Matrix<double, Dim, 1> const& referencePoint )
{
    auto target = inverseBenchmarkMap<Dim, Order, RealDim>( referencePoint );
    if constexpr ( Projected )
    {
        static_assert( Dim == 2 && RealDim == 3,
                       "projected benchmark requires a surface in 3D" );
        constexpr double epsilon = 1e-6;
        Eigen::Matrix<double, 3, 2> jacobian;
        for ( int direction = 0; direction < 2; ++direction )
        {
            auto plus = referencePoint;
            auto minus = referencePoint;
            plus( direction ) += epsilon;
            minus( direction ) -= epsilon;
            jacobian.col( direction ) =
                ( inverseBenchmarkMap<Dim, Order, RealDim>( plus ) -
                  inverseBenchmarkMap<Dim, Order, RealDim>( minus ) ) /
                ( 2.0 * epsilon );
        }
        Eigen::Vector3d normal =
            jacobian.col( 0 ).cross( jacobian.col( 1 ) ).normalized();
        target += 0.075 * normal;
    }
    return target;
}

/**
 * @struct inverse_benchmark_data
 * @brief Reusable untimed mapping, coefficients, and inverse facade.
 * @tparam Dim reference dimension
 * @tparam Order effective geometric order
 * @tparam RealDim ambient physical dimension
 * @tparam RuntimeOrderV whether the order is selected at runtime
 */
template <int Dim, int Order, int RealDim, bool RuntimeOrderV>
struct inverse_benchmark_data
{
    using geomap_type =
        inverse_benchmark_geomap_type<Dim, Order, RealDim, RuntimeOrderV>;
    using geometry_type = inverse_benchmark_geometry<geomap_type>;
    using inverse_type = typename geomap_type::Inverse;

    /** Geometric mapping shared with the inverse facade. */
    std::shared_ptr<geomap_type> geomap;
    /** Physical geometric coefficients. */
    geometry_type geometry;
    /** Stateful compatibility facade using the shared numerical kernel. */
    inverse_type inverse;

    /** @brief Construct all benchmark data outside the timed region. */
    inverse_benchmark_data()
        : geomap( makeInverseBenchmarkGeomap<Dim, Order, RealDim,
                                             RuntimeOrderV>() ),
          geometry( makeInverseBenchmarkGeometry<Dim, Order, RealDim>(
              geomap ) ),
          inverse( geomap, geometry )
    {
    }
};

/**
 * @brief Benchmark construction of the inverse facade and its cached factors.
 * @tparam Dim reference dimension
 * @tparam Order effective geometric order
 * @tparam RealDim ambient physical dimension
 * @tparam RuntimeOrderV whether the order is selected at runtime
 * @param state Google Benchmark state
 */
template <int Dim, int Order, int RealDim, bool RuntimeOrderV>
void BM_GeoMapInverseConstruction( benchmark::State& state )
{
    using geomap_type =
        inverse_benchmark_geomap_type<Dim, Order, RealDim, RuntimeOrderV>;
    using inverse_type = typename geomap_type::Inverse;
    auto geomap =
        makeInverseBenchmarkGeomap<Dim, Order, RealDim, RuntimeOrderV>();
    auto geometry =
        makeInverseBenchmarkGeometry<Dim, Order, RealDim>( geomap );

    for ( [[maybe_unused]] auto _ : state )
    {
        inverse_type inverse( geomap, geometry );
        benchmark::DoNotOptimize( inverse.J() );
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed( state.iterations() );
}

/**
 * @brief Benchmark a cold single-point inverse or closest-point operation.
 * @tparam Dim reference dimension
 * @tparam Order effective geometric order
 * @tparam RealDim ambient physical dimension
 * @tparam RuntimeOrderV whether the order is selected at runtime
 * @tparam Projected whether the target is offset from an embedded surface
 * @param state Google Benchmark state
 */
template <int Dim, int Order, int RealDim, bool RuntimeOrderV, bool Projected>
void BM_GeoMapInverseSingle( benchmark::State& state )
{
    inverse_benchmark_data<Dim, Order, RealDim, RuntimeOrderV> data;
    auto const referencePoint = inverseBenchmarkReferencePoint<Dim>();
    auto const target = inverseBenchmarkNode(
        inverseBenchmarkTarget<Dim, Order, RealDim, Projected>(
            referencePoint ) );
    auto const initialGuess = inverseBenchmarkInitialGuess<Dim>();
    typename decltype( data )::inverse_type::inverse_options_type options;
    options.warm_start = false;
    if constexpr ( Projected )
        options.domain_policy =
            Feel::geometry::InverseDomainPolicy::constrained;

    auto const validation = data.inverse.inverse(
        target, options, initialGuess );
    if ( !validation.converged() )
    {
        state.SkipWithError( "inverse benchmark setup did not converge" );
        return;
    }

    for ( [[maybe_unused]] auto _ : state )
    {
        auto result = data.inverse.inverse( target, options, initialGuess );
        benchmark::DoNotOptimize( result.reference_point.data() );
        benchmark::DoNotOptimize( result.status );
        benchmark::ClobberMemory();
    }
    state.counters["iterations"] = validation.iterations;
    state.SetItemsProcessed( state.iterations() );
}

/**
 * @brief Benchmark batched inverse or embedded closest-point operations.
 * @tparam Dim reference dimension
 * @tparam Order effective geometric order
 * @tparam RealDim ambient physical dimension
 * @tparam RuntimeOrderV whether the order is selected at runtime
 * @tparam Projected whether targets are offset from an embedded surface
 * @tparam WarmStart whether consecutive targets reuse the preceding result
 * @tparam EigenInput whether targets use the Eigen-native batch overload
 * @tparam PointwiseAffine force the compatibility pointwise affine path
 * @param state Google Benchmark state; argument is the target count
 */
template <int Dim, int Order, int RealDim, bool RuntimeOrderV,
          bool Projected, bool WarmStart, bool EigenInput = false,
          bool PointwiseAffine = false>
void BM_GeoMapInverseBatch( benchmark::State& state )
{
    inverse_benchmark_data<Dim, Order, RealDim, RuntimeOrderV> data;
    if constexpr ( PointwiseAffine )
        data.inverse.setGeometryEvaluationStrategy(
            geometry::geometry_evaluation_strategy::pointwise );
    auto const targetCount = static_cast<std::size_t>( state.range( 0 ) );
    using targets_type = std::conditional_t<
        EigenInput,
        typename decltype( data )::inverse_type::eigen_physical_points_type,
        typename decltype( data )::geomap_type::matrix_node_t_type>;
    targets_type targets( RealDim, targetCount );
    for ( std::size_t q = 0; q < targetCount; ++q )
    {
        double const parameter = targetCount > 1
                                     ? -1.0 + 2.0 * double( q ) /
                                                  double( targetCount - 1 )
                                     : 0.0;
        auto const target = inverseBenchmarkTarget<Dim, Order, RealDim,
                                                   Projected>(
            inverseBenchmarkReferencePoint<Dim>( parameter ) );
        for ( int component = 0; component < RealDim; ++component )
            targets( component, q ) = target( component );
    }

    typename decltype( data )::inverse_type::inverse_options_type options;
    options.warm_start = WarmStart;
    if constexpr ( Projected )
        options.domain_policy =
            Feel::geometry::InverseDomainPolicy::constrained;
    auto const validation = data.inverse.inverse( targets, options );
    if ( !std::ranges::all_of( validation,
                               []( auto const& result )
                               { return result.converged(); } ) )
    {
        state.SkipWithError( "batched inverse benchmark setup did not converge" );
        return;
    }

    double totalIterations = 0.0;
    for ( auto const& result : validation )
        totalIterations += result.iterations;
    for ( [[maybe_unused]] auto _ : state )
    {
        auto results = data.inverse.inverse( targets, options );
        benchmark::DoNotOptimize( results.data() );
        benchmark::ClobberMemory();
    }
    state.counters["iterations/point"] =
        totalIterations / double( targetCount );
    state.SetItemsProcessed( state.iterations() * targetCount );
}
} // namespace

BENCHMARK_TEMPLATE( BM_GeoMapInverseConstruction, 2, 2, 2, false )
    ->Name( "geomap_inverse/construction/2d/g2/static" );
BENCHMARK_TEMPLATE( BM_GeoMapInverseConstruction, 2, 2, 2, true )
    ->Name( "geomap_inverse/construction/2d/g2/runtime" );
BENCHMARK_TEMPLATE( BM_GeoMapInverseConstruction, 3, 2, 3, false )
    ->Name( "geomap_inverse/construction/3d/g2/static" );
BENCHMARK_TEMPLATE( BM_GeoMapInverseConstruction, 3, 2, 3, true )
    ->Name( "geomap_inverse/construction/3d/g2/runtime" );
BENCHMARK_TEMPLATE( BM_GeoMapInverseConstruction, 2, 2, 3, false )
    ->Name( "geomap_inverse/construction/surface2d-in-3d/g2/static" );
BENCHMARK_TEMPLATE( BM_GeoMapInverseConstruction, 2, 2, 3, true )
    ->Name( "geomap_inverse/construction/surface2d-in-3d/g2/runtime" );

BENCHMARK_TEMPLATE( BM_GeoMapInverseSingle, 2, 1, 2, false, false )
    ->Name( "geomap_inverse/single/2d/g1/static/exact" );
BENCHMARK_TEMPLATE( BM_GeoMapInverseSingle, 2, 1, 2, true, false )
    ->Name( "geomap_inverse/single/2d/g1/runtime/exact" );
BENCHMARK_TEMPLATE( BM_GeoMapInverseSingle, 2, 2, 2, false, false )
    ->Name( "geomap_inverse/single/2d/g2/static/exact" );
BENCHMARK_TEMPLATE( BM_GeoMapInverseSingle, 2, 2, 2, true, false )
    ->Name( "geomap_inverse/single/2d/g2/runtime/exact" );
BENCHMARK_TEMPLATE( BM_GeoMapInverseSingle, 2, 3, 2, false, false )
    ->Name( "geomap_inverse/single/2d/g3/static/exact" );

BENCHMARK_TEMPLATE( BM_GeoMapInverseSingle, 3, 1, 3, false, false )
    ->Name( "geomap_inverse/single/3d/g1/static/exact" );
BENCHMARK_TEMPLATE( BM_GeoMapInverseSingle, 3, 1, 3, true, false )
    ->Name( "geomap_inverse/single/3d/g1/runtime/exact" );
BENCHMARK_TEMPLATE( BM_GeoMapInverseSingle, 3, 2, 3, false, false )
    ->Name( "geomap_inverse/single/3d/g2/static/exact" );
BENCHMARK_TEMPLATE( BM_GeoMapInverseSingle, 3, 2, 3, true, false )
    ->Name( "geomap_inverse/single/3d/g2/runtime/exact" );

BENCHMARK_TEMPLATE( BM_GeoMapInverseSingle, 2, 2, 3, false, false )
    ->Name( "geomap_inverse/single/surface2d-in-3d/g2/static/exact" );
BENCHMARK_TEMPLATE( BM_GeoMapInverseSingle, 2, 2, 3, true, false )
    ->Name( "geomap_inverse/single/surface2d-in-3d/g2/runtime/exact" );
BENCHMARK_TEMPLATE( BM_GeoMapInverseSingle, 2, 2, 3, false, true )
    ->Name( "geomap_inverse/single/surface2d-in-3d/g2/static/closest" );
BENCHMARK_TEMPLATE( BM_GeoMapInverseSingle, 2, 2, 3, true, true )
    ->Name( "geomap_inverse/single/surface2d-in-3d/g2/runtime/closest" );

BENCHMARK_TEMPLATE( BM_GeoMapInverseBatch, 2, 1, 2, false, false, false,
                    true )
    ->Name( "geomap_inverse/batch/2d/g1/static/exact/eigen-input" )
    ->Arg( 64 )
    ->Arg( 256 );
BENCHMARK_TEMPLATE( BM_GeoMapInverseBatch, 2, 1, 2, true, false, false,
                    true )
    ->Name( "geomap_inverse/batch/2d/g1/runtime/exact/eigen-input" )
    ->Arg( 64 )
    ->Arg( 256 );
BENCHMARK_TEMPLATE( BM_GeoMapInverseBatch, 3, 1, 3, false, false, false,
                    true )
    ->Name( "geomap_inverse/batch/3d/g1/static/exact/eigen-input" )
    ->Arg( 64 );
BENCHMARK_TEMPLATE( BM_GeoMapInverseBatch, 3, 1, 3, true, false, false,
                    true )
    ->Name( "geomap_inverse/batch/3d/g1/runtime/exact/eigen-input" )
    ->Arg( 64 );
BENCHMARK_TEMPLATE( BM_GeoMapInverseBatch, 2, 1, 2, true, false, false,
                    true, true )
    ->Name(
        "geomap_inverse/batch/2d/g1/runtime/exact/eigen-input/pointwise" )
    ->Arg( 64 )
    ->Arg( 256 );
BENCHMARK_TEMPLATE( BM_GeoMapInverseBatch, 3, 1, 3, true, false, false,
                    true, true )
    ->Name(
        "geomap_inverse/batch/3d/g1/runtime/exact/eigen-input/pointwise" )
    ->Arg( 64 );

BENCHMARK_TEMPLATE( BM_GeoMapInverseBatch, 2, 2, 2, false, false, false )
    ->Name( "geomap_inverse/batch/2d/g2/static/exact/cold" )
    ->Arg( 8 )
    ->Arg( 64 )
    ->Arg( 256 );
BENCHMARK_TEMPLATE( BM_GeoMapInverseBatch, 2, 2, 2, false, false, true )
    ->Name( "geomap_inverse/batch/2d/g2/static/exact/warm" )
    ->Arg( 8 )
    ->Arg( 64 )
    ->Arg( 256 );
BENCHMARK_TEMPLATE( BM_GeoMapInverseBatch, 2, 2, 2, true, false, true )
    ->Name( "geomap_inverse/batch/2d/g2/runtime/exact/warm" )
    ->Arg( 8 )
    ->Arg( 64 )
    ->Arg( 256 );
BENCHMARK_TEMPLATE( BM_GeoMapInverseBatch, 2, 2, 2, true, false, true,
                    true )
    ->Name( "geomap_inverse/batch/2d/g2/runtime/exact/warm/eigen-input" )
    ->Arg( 256 );

BENCHMARK_TEMPLATE( BM_GeoMapInverseBatch, 3, 2, 3, false, false, false )
    ->Name( "geomap_inverse/batch/3d/g2/static/exact/cold" )
    ->Arg( 8 )
    ->Arg( 64 );
BENCHMARK_TEMPLATE( BM_GeoMapInverseBatch, 3, 2, 3, false, false, true )
    ->Name( "geomap_inverse/batch/3d/g2/static/exact/warm" )
    ->Arg( 8 )
    ->Arg( 64 );
BENCHMARK_TEMPLATE( BM_GeoMapInverseBatch, 3, 2, 3, true, false, true )
    ->Name( "geomap_inverse/batch/3d/g2/runtime/exact/warm" )
    ->Arg( 8 )
    ->Arg( 64 );
BENCHMARK_TEMPLATE( BM_GeoMapInverseBatch, 3, 2, 3, true, false, true,
                    true )
    ->Name( "geomap_inverse/batch/3d/g2/runtime/exact/warm/eigen-input" )
    ->Arg( 64 );

BENCHMARK_TEMPLATE( BM_GeoMapInverseBatch, 2, 2, 3, false, true, false )
    ->Name(
        "geomap_inverse/batch/surface2d-in-3d/g2/static/closest/cold" )
    ->Arg( 8 )
    ->Arg( 64 )
    ->Arg( 256 );
BENCHMARK_TEMPLATE( BM_GeoMapInverseBatch, 2, 2, 3, false, true, true )
    ->Name(
        "geomap_inverse/batch/surface2d-in-3d/g2/static/closest/warm" )
    ->Arg( 8 )
    ->Arg( 64 )
    ->Arg( 256 );
BENCHMARK_TEMPLATE( BM_GeoMapInverseBatch, 2, 2, 3, true, true, true )
    ->Name(
        "geomap_inverse/batch/surface2d-in-3d/g2/runtime/closest/warm" )
    ->Arg( 8 )
    ->Arg( 64 )
    ->Arg( 256 );
BENCHMARK_TEMPLATE( BM_GeoMapInverseBatch, 2, 2, 3, true, true, true,
                    true )
    ->Name( "geomap_inverse/batch/surface2d-in-3d/g2/runtime/closest/warm/eigen-input" )
    ->Arg( 256 );

/** @brief Initialize Feel++ and execute the registered inverse benchmarks. */
int main( int argc, char** argv )
{
    Environment environment(
        _argc = argc,
        _argv = argv,
        _about = about( _name = "feelpp_bench_geomap_inverse",
                        _author = "Feel++ Consortium",
                        _email = "feelpp-devel@feelpp.org" ) );

    benchmark::Initialize( &argc, argv );
    benchmark::RunSpecifiedBenchmarks();
}
