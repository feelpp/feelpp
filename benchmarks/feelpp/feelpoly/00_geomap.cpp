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
#include <string>

using namespace Feel;

using simplex2_g1_mesh_type = Mesh<Simplex<2>>;
using simplex2_g2_mesh_type = Mesh<Simplex<2, 2, 2>>;
using simplex2_static_g1_type = GeoMap<2, 1>;
using simplex2_runtime_type = GeoMap<2, Dynamic>;
using simplex2_static_g2_type = GeoMap<2, 2>;

/**
 * @brief Benchmark a geometric-mapping context update over a complete mesh.
 *
 * Mesh construction and context precomputation are excluded from the timed
 * region. The first benchmark argument controls the mesh resolution and the
 * second controls the integration order used to build the point table.
 *
 * @tparam MeshType mesh type supplying the geometric mapping
 * @tparam ContextFlags compile-time geometric context mask
 * @tparam GeomapType geometric mapping type used by the context
 * @tparam EffectiveOrder runtime order when `GeomapType` is dynamic
 * @tparam EvaluationStrategy pointwise, batched, or automatic derivative kernel
 * @param state Google Benchmark state and integer arguments
 */
template <typename MeshType, size_type ContextFlags,
          typename GeomapType = typename MeshType::gm_type,
          int EffectiveOrder = GeomapType::nOrder,
          geometry::geometry_evaluation_strategy EvaluationStrategy =
              geometry::geometry_evaluation_strategy::automatic>
void BM_GeoMapContextUpdate( benchmark::State& state )
{
    using mesh_type = MeshType;
    using geomap_type = GeomapType;
    using integration_method_type = im_t<typename mesh_type::element_type, double>;
    using context_type =
        typename geomap_type::template Context<typename mesh_type::element_type>;

    auto geomap = []()
    {
        if constexpr ( geomap_type::is_order_dynamic )
            return std::make_shared<geomap_type>(
                RuntimeOrder{ EffectiveOrder } );
        else
            return std::make_shared<geomap_type>();
    }();

    std::map<int, std::shared_ptr<mesh_type>> meshes;
    std::map<int, std::shared_ptr<context_type>> contexts;
    size_type itemCount = 0;

    for ( [[maybe_unused]] auto _ : state )
    {
        state.PauseTiming();
        if ( !meshes.contains( state.range( 0 ) ) )
        {
            auto const meshSize = 1. / state.range( 0 );
            meshes[state.range( 0 )] = createGMSHMesh(
                _mesh = new mesh_type,
                _desc = domain(
                    _name = "geomap-context-simplex2-g" +
                            std::to_string( mesh_type::shape_type::nOrder ),
                    _shape = "simplex", _dim = mesh_type::nDim,
                    _order = mesh_type::shape_type::nOrder, _h = meshSize ),
                _verbose = false,
                _straighten = false,
                _update = 1,
                _h = meshSize );
        }
        auto const& mesh = meshes.at( state.range( 0 ) );

        if ( !contexts.contains( state.range( 1 ) ) )
        {
            auto precompute = geomap->preCompute(
                integration_method_type( state.range( 1 ) ).points() );
            auto const& initialElement = mesh->beginElement()->second;
            contexts[state.range( 1 )] =
                geomap->template context<ContextFlags>( initialElement,
                                                        precompute );
            contexts[state.range( 1 )]->setGeometryEvaluationStrategy(
                EvaluationStrategy );
        }
        auto const& context = contexts.at( state.range( 1 ) );
        state.ResumeTiming();

        for ( auto const& element : mesh->elements() )
            context->template update<ContextFlags>( element.second );
        itemCount = mesh->numElements() * context->pc()->nPoints();
    }
    if ( !contexts.empty() )
        state.counters["computed-points"] =
            contexts.begin()->second->nComputedPoints();
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

BENCHMARK_TEMPLATE(
    BM_GeoMapContextUpdate, simplex2_g1_mesh_type,
    vm::JACOBIAN | vm::KB, simplex2_static_g1_type, 1,
    geometry::geometry_evaluation_strategy::automatic )
    ->Name( "geomap/context/2d/g1/static/automatic" )
    ->Unit( benchmark::kMillisecond )
    ->Args( { 8, 4 } );
BENCHMARK_TEMPLATE(
    BM_GeoMapContextUpdate, simplex2_g1_mesh_type,
    vm::JACOBIAN | vm::KB, simplex2_runtime_type, 1,
    geometry::geometry_evaluation_strategy::automatic )
    ->Name( "geomap/context/2d/g1/runtime/automatic" )
    ->Unit( benchmark::kMillisecond )
    ->Args( { 8, 4 } );
BENCHMARK_TEMPLATE(
    BM_GeoMapContextUpdate, simplex2_g2_mesh_type,
    vm::JACOBIAN | vm::KB | vm::SECOND_DERIVATIVE,
    simplex2_static_g2_type, 2,
    geometry::geometry_evaluation_strategy::automatic )
    ->Name( "geomap/context/2d/g2/static/automatic" )
    ->Unit( benchmark::kMillisecond )
    ->Args( { 4, 4 } )
    ->Args( { 8, 4 } );
BENCHMARK_TEMPLATE(
    BM_GeoMapContextUpdate, simplex2_g2_mesh_type,
    vm::JACOBIAN | vm::KB | vm::SECOND_DERIVATIVE,
    simplex2_runtime_type, 2,
    geometry::geometry_evaluation_strategy::automatic )
    ->Name( "geomap/context/2d/g2/runtime/automatic" )
    ->Unit( benchmark::kMillisecond )
    ->Args( { 4, 4 } )
    ->Args( { 8, 4 } );
BENCHMARK_TEMPLATE(
    BM_GeoMapContextUpdate, simplex2_g2_mesh_type,
    vm::JACOBIAN | vm::KB | vm::SECOND_DERIVATIVE,
    simplex2_runtime_type, 2,
    geometry::geometry_evaluation_strategy::pointwise )
    ->Name( "geomap/context/2d/g2/runtime/pointwise" )
    ->Unit( benchmark::kMillisecond )
    ->Args( { 4, 4 } )
    ->Args( { 8, 4 } );
BENCHMARK_TEMPLATE(
    BM_GeoMapContextUpdate, simplex2_g2_mesh_type,
    vm::JACOBIAN | vm::KB | vm::SECOND_DERIVATIVE,
    simplex2_runtime_type, 2,
    geometry::geometry_evaluation_strategy::batched )
    ->Name( "geomap/context/2d/g2/runtime/batched" )
    ->Unit( benchmark::kMillisecond )
    ->Args( { 4, 4 } )
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
