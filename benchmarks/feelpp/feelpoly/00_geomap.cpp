/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file 00_geomap.cpp
 * @brief Microbenchmarks for mapping contexts and coefficient representations.
 */

#include <benchmark/benchmark.h>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelpoly/geometryoperations.hpp>

#include <map>
#include <memory>

using namespace Feel;

/**
 * @brief Benchmark a geometric-mapping context update over a complete mesh.
 *
 * Mesh construction and context precomputation are excluded from the timed
 * region. The first benchmark argument controls the mesh resolution and the
 * second controls the integration order used to build the point table.
 *
 * @tparam MeshType mesh type supplying the geometric mapping
 * @tparam ContextFlags compile-time geometric context mask
 * @param state Google Benchmark state and integer arguments
 */
template <typename MeshType, size_type ContextFlags>
void BM_GeoMapContextUpdate( benchmark::State& state )
{
    using mesh_type = MeshType;
    using integration_method_type = im_t<typename mesh_type::element_type, double>;
    using context_type =
        typename mesh_type::gm_type::template Context<typename mesh_type::element_type>;

    std::map<int, std::shared_ptr<mesh_type>> meshes;
    std::map<int, std::shared_ptr<context_type>> contexts;
    size_type itemCount = 0;

    for ( [[maybe_unused]] auto _ : state )
    {
        state.PauseTiming();
        if ( !meshes.contains( state.range( 0 ) ) )
        {
            meshes[state.range( 0 )] = loadMesh(
                _mesh = new mesh_type,
                _savehdf5 = false,
                _verbose = false,
                _update = 1,
                _h = 1. / state.range( 0 ) );
        }
        auto const& mesh = meshes.at( state.range( 0 ) );

        if ( !contexts.contains( state.range( 1 ) ) )
        {
            auto precompute = mesh->gm()->preCompute(
                integration_method_type( state.range( 1 ) ).points() );
            auto const& initialElement = mesh->beginElement()->second;
            contexts[state.range( 1 )] =
                mesh->gm()->template context<ContextFlags>( initialElement, precompute );
        }
        auto const& context = contexts.at( state.range( 1 ) );
        state.ResumeTiming();

        for ( auto const& element : mesh->elements() )
            context->template update<ContextFlags>( element.second );
        itemCount = mesh->numElements() * context->pc()->nPoints();
    }
    state.SetItemsProcessed( state.iterations() * itemCount );
}

/**
 * @brief Benchmark construction of the nodal representation compatibility view.
 * @tparam Order static geometric order
 * @param state Google Benchmark state
 */
template <int Order>
void BM_NodalRepresentationConstruction( benchmark::State& state )
{
    GeoMap<2, Order> geomap;
    for ( [[maybe_unused]] auto _ : state )
    {
        auto representation = geomap.geometryRepresentation();
        benchmark::DoNotOptimize( representation.bernsteinVandermonde().data() );
    }
    state.SetItemsProcessed( state.iterations() );
}

/**
 * @brief Benchmark explicit nodal-to-Bezier coefficient conversion.
 * @tparam Order static geometric order
 * @param state Google Benchmark state
 */
template <int Order>
void BM_NodalToBezierConversion( benchmark::State& state )
{
    GeoMap<2, Order> geomap;
    auto nodal = geomap.geometryRepresentation();
    geometry::BernsteinBezier<2, Order> bezier;
    Eigen::MatrixXd values( 2, nodal.coefficientCount() );
    values.setRandom();
    geometry::CoefficientTable coefficientTable( nodal, values );
    for ( [[maybe_unused]] auto _ : state )
    {
        auto converted = geometry::convert( coefficientTable, bezier );
        benchmark::DoNotOptimize( converted.values().data() );
    }
    state.SetItemsProcessed( state.iterations() * nodal.coefficientCount() );
}

/**
 * @brief Benchmark direct Bernstein/Bezier geometry evaluation.
 * @tparam Order static geometric order
 * @param state Google Benchmark state; argument is the point count
 */
template <int Order>
void BM_BezierEvaluation( benchmark::State& state )
{
    geometry::BernsteinBezier<2, Order> representation;
    Eigen::MatrixXd values( 2, representation.coefficientCount() );
    values.setRandom();
    geometry::CoefficientTable coefficientTable( representation, values );
    typename decltype( representation )::points_type points( 2, state.range( 0 ) );
    points.setRandom();
    for ( [[maybe_unused]] auto _ : state )
    {
        auto physicalPoints = geometry::evaluate( coefficientTable, points );
        benchmark::DoNotOptimize( physicalPoints.data() );
    }
    state.SetItemsProcessed( state.iterations() * state.range( 0 ) );
}

BENCHMARK_TEMPLATE( BM_GeoMapContextUpdate, Mesh<Simplex<2>>, 0 )
    ->Unit( benchmark::kMillisecond )
    ->Args( { 4, 1 } )
    ->Args( { 4, 2 } )
    ->Args( { 4, 4 } )
    ->Args( { 8, 1 } )
    ->Args( { 8, 2 } )
    ->Args( { 8, 4 } );
BENCHMARK_TEMPLATE( BM_GeoMapContextUpdate, Mesh<Simplex<2>>, vm::JACOBIAN )
    ->Unit( benchmark::kMillisecond )
    ->Args( { 4, 1 } )
    ->Args( { 4, 2 } )
    ->Args( { 4, 4 } )
    ->Args( { 8, 1 } )
    ->Args( { 8, 2 } )
    ->Args( { 8, 4 } );
BENCHMARK_TEMPLATE( BM_GeoMapContextUpdate, Mesh<Simplex<2>>, vm::JACOBIAN | vm::KB )
    ->Unit( benchmark::kMillisecond )
    ->Args( { 4, 1 } )
    ->Args( { 4, 2 } )
    ->Args( { 4, 4 } )
    ->Args( { 8, 1 } )
    ->Args( { 8, 2 } )
    ->Args( { 8, 4 } );
BENCHMARK_TEMPLATE( BM_GeoMapContextUpdate, Mesh<Simplex<2>>,
                    vm::JACOBIAN | vm::KB | vm::SECOND_DERIVATIVE )
    ->Unit( benchmark::kMillisecond )
    ->Args( { 4, 1 } )
    ->Args( { 4, 2 } )
    ->Args( { 4, 4 } )
    ->Args( { 8, 1 } )
    ->Args( { 8, 2 } )
    ->Args( { 8, 4 } );

BENCHMARK_TEMPLATE( BM_NodalRepresentationConstruction, 1 );
BENCHMARK_TEMPLATE( BM_NodalRepresentationConstruction, 2 );
BENCHMARK_TEMPLATE( BM_NodalRepresentationConstruction, 3 );
BENCHMARK_TEMPLATE( BM_NodalToBezierConversion, 1 );
BENCHMARK_TEMPLATE( BM_NodalToBezierConversion, 2 );
BENCHMARK_TEMPLATE( BM_NodalToBezierConversion, 3 );
BENCHMARK_TEMPLATE( BM_BezierEvaluation, 1 )->Arg( 4 )->Arg( 16 )->Arg( 64 );
BENCHMARK_TEMPLATE( BM_BezierEvaluation, 2 )->Arg( 4 )->Arg( 16 )->Arg( 64 );
BENCHMARK_TEMPLATE( BM_BezierEvaluation, 3 )->Arg( 4 )->Arg( 16 )->Arg( 64 );

int main( int argc, char** argv )
{
    Environment environment(
        _argc = argc,
        _argv = argv,
        _about = about( _name = "feelpp_bench_geomap",
                        _author = "Feel++ Consortium",
                        _email = "feelpp-devel@feelpp.org" ) );

    benchmark::Initialize( &argc, argv );
    benchmark::RunSpecifiedBenchmarks();
}
