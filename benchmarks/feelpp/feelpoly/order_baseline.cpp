/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#include <benchmark/benchmark.h>

#include <feel/feelcore/environment.hpp>
#include <feel/feelpoly/lagrange.hpp>

#include "allocation_probe.hpp"

#include <array>
#include <initializer_list>
#include <string>
#include <utility>

namespace
{
namespace allocations = Feel::benchmarking::allocation_probe;

/** @brief Scalar equispaced triangular element used by both order modes.
 * @tparam OrderSpec Static degree or `Feel::Dynamic`.
 */
template<int OrderSpec>
using element_type = typename Feel::Lagrange<OrderSpec, Feel::Scalar, Feel::Continuous,
    Feel::PointSetEquiSpaced>::template apply<2, 2>::type;

/** @brief Construct the same element family using its selected order policy.
 * @tparam OrderSpec Static degree or `Feel::Dynamic`.
 * @param order Effective degree for runtime order; static cases use OrderSpec.
 * @return Constructed reference finite element.
 */
template<int OrderSpec>
auto makeElement( int order )
{
    if constexpr ( OrderSpec == Feel::Dynamic )
        return element_type<OrderSpec>{ Feel::RuntimeOrder{ order } };
    else
        return element_type<OrderSpec>{};
}

/** @brief Time an operation and report its mean C++ allocation calls/bytes.
 * @tparam Operation Nullary callable performing one measured operation.
 * @param state Google Benchmark iteration state receiving allocation counters.
 * @param operation Callable whose setup is outside the measured iteration loop.
 * @note Framework/setup allocations are excluded; probe overhead remains timed.
 */
template<typename Operation>
void benchmarkOperation( benchmark::State& state, Operation&& operation )
{
    allocations::AllocationSnapshot total;
    for ( [[maybe_unused]] auto iteration : state )
    {
        allocations::begin();
        try
        {
            operation();
        }
        catch ( ... )
        {
            (void)allocations::end();
            throw;
        }
        auto const snapshot = allocations::end();
        total.calls += snapshot.calls;
        total.bytes += snapshot.bytes;
    }
    state.counters["cpp-allocation-calls"] = benchmark::Counter(
        static_cast<double>( total.calls ), benchmark::Counter::kAvgIterations );
    state.counters["cpp-requested-bytes"] = benchmark::Counter(
        static_cast<double>( total.bytes ), benchmark::Counter::kAvgIterations );
}

/** @brief Expose every output buffer to the optimizer's memory barrier.
 * @tparam Matrix Contiguous uBLAS matrix returned by the polynomial API.
 * @param matrix Computed matrix whose contents must remain observable.
 */
template<typename Matrix>
void consumeMatrix( Matrix const& matrix )
{
    benchmark::DoNotOptimize( matrix.data().begin() );
    benchmark::ClobberMemory();
}

/** @brief Measure reference-element construction and destruction.
 * @tparam OrderSpec Static degree or `Feel::Dynamic`.
 * @param state Benchmark state; argument zero is the effective degree.
 */
template<int OrderSpec>
void benchmarkConstruction( benchmark::State& state )
{
    benchmarkOperation( state, [&]()
    {
        auto element = makeElement<OrderSpec>( static_cast<int>( state.range( 0 ) ) );
        consumeMatrix( element.coeff() );
    } );
    state.SetItemsProcessed( state.iterations() );
}

/** @brief Measure values, reference gradients or reference Hessians at fixed points.
 * @tparam OrderSpec Static degree or `Feel::Dynamic`.
 * @tparam DerivativeOrder Zero for values, one for gradients, two for Hessians.
 * @param state Benchmark state; arguments are effective degree and point count.
 * @note Element and coordinate construction precede the timed iteration loop.
 */
template<int OrderSpec, unsigned DerivativeOrder>
void benchmarkTabulation( benchmark::State& state )
{
    auto element = makeElement<OrderSpec>( static_cast<int>( state.range( 0 ) ) );
    typename decltype( element )::points_type points( 2, state.range( 1 ) );
    for ( unsigned q = 0; q < points.size2(); ++q )
    {
        points( 0, q ) = -0.8 + q * 0.1;
        points( 1, q ) = -0.7 + q * 0.05;
    }
    benchmarkOperation( state, [&]()
    {
        if constexpr ( DerivativeOrder == 0 )
            consumeMatrix( element.evaluate( points ) );
        else if constexpr ( DerivativeOrder == 1 )
        {
            auto const gradients = element.derivate( points );
            for ( auto const& gradient : gradients )
                consumeMatrix( gradient );
        }
        else
        {
            auto const hessians = element.derivate2( points );
            for ( unsigned i = 0; i < 2; ++i )
                for ( unsigned j = 0; j < 2; ++j )
                    consumeMatrix( hessians( i, j ) );
        }
    } );
    state.SetItemsProcessed( state.iterations() * state.range( 1 ) );
}

/** @brief Measure extraction of one edge from a runtime Fekete point set.
 * @param state Benchmark state; argument zero is the point-set order.
 */
void benchmarkFeketeEdgeQuery( benchmark::State& state )
{
    Feel::PointSetFekete<Feel::Simplex<2>, Feel::Dynamic, double> points{
        Feel::RuntimeOrder{ static_cast<int>( state.range( 0 ) ) } };
    benchmarkOperation( state, [&]()
    {
        consumeMatrix( points.pointsBySubEntity( 1, 0 ) );
    } );
    state.SetItemsProcessed( state.iterations() );
}

/** @brief Register construction and one/four-point tabulation datasets.
 * @tparam OrderSpec Static degree or `Feel::Dynamic`.
 * @param name Descriptive benchmark prefix identifying the cell and order policy.
 * @param orders Effective degrees; static registrations use their template degree.
 */
template<int OrderSpec>
void registerElementBenchmarks( std::string const& name, std::initializer_list<int> orders )
{
    auto* construction = benchmark::RegisterBenchmark(
        ( name + "/construct" ).c_str(), &benchmarkConstruction<OrderSpec> )
        ->Unit( benchmark::kMicrosecond );
    for ( int order : orders )
        construction->Arg( order );

    for ( auto const& [operation, function] : std::array{
        std::pair{ "values", &benchmarkTabulation<OrderSpec, 0> },
        std::pair{ "gradients", &benchmarkTabulation<OrderSpec, 1> },
        std::pair{ "hessians", &benchmarkTabulation<OrderSpec, 2> } } )
    {
        auto* tabulation = benchmark::RegisterBenchmark(
            ( name + "/" + operation ).c_str(), function )->Unit( benchmark::kMicrosecond );
        for ( int order : orders )
            for ( int pointCount : { 1, 4 } )
                tabulation->Args( { order, pointCount } );
    }
}
} // namespace

/** @brief Initialize Feel++ and execute the polynomial-order benchmarks.
 * @param argc Command-line argument count.
 * @param argv Feel++ and Google Benchmark command-line options.
 * @return Zero after the selected benchmark run or on a non-master MPI rank.
 */
int main( int argc, char** argv )
{
    Feel::Environment environment( argc, argv );
    if ( !Feel::Environment::isMasterRank() )
        return 0;
    benchmark::Initialize( &argc, argv );
    registerElementBenchmarks<3>( "simplex2/static", { 3 } );
    registerElementBenchmarks<Feel::Dynamic>( "simplex2/runtime", { 2, 3, 5 } );
    benchmark::RegisterBenchmark( "fekete2/runtime/edge_query", &benchmarkFeketeEdgeQuery )
        ->Arg( 3 )->Unit( benchmark::kMicrosecond );
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
}
