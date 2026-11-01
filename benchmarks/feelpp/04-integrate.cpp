/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2019-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 * @file 04-integrate.cpp
 * @brief End-to-end Feel++ DSEL assembly benchmarks on curved geometry.
 */

#include <benchmark/benchmark.h>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelpoly/pointsetquadrature.hpp>
#include <feel/feelpoly/quadratureevaluation.hpp>
#include <feel/feelvf/vf.hpp>

#include <algorithm>
#include <cmath>
#include <string>

using namespace Feel;
using namespace Feel::vf;

namespace
{

/** @brief Local bilinear form assembled by a benchmark instance. */
enum class dsel_form_kind
{
    mass,
    stiffness
};

/** @brief Local assembly route selected by a benchmark instance. */
enum class dsel_contraction_path
{
    entrywise,
    local_tensor
};

/** @brief Unit geometric weight used by the expression-only benchmark. */
struct unit_quadrature_geometry
{
    /** @return unit volume Jacobian. */
    double J( std::uint32_t /*q*/ ) const noexcept { return 1.; }

    /** @return unit boundary Jacobian. */
    double normalNorm( std::uint32_t /*q*/ ) const noexcept { return 1.; }
};

/**
 * @brief Scalar mass-like tensor used to measure the pointwise contract.
 *
 * Storage is identical in the pointwise and batch fixtures so the benchmark
 * isolates expression traversal rather than table construction.
 */
struct pointwise_mass_tensor
{
    using table_type = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                     Eigen::RowMajor>;

    pointwise_mass_tensor( Eigen::Index dofs, Eigen::Index points )
        : coefficient( points ), testValues( dofs, points ),
          trialValues( dofs, points )
    {
        for ( Eigen::Index q = 0; q < points; ++q )
        {
            coefficient[q] = static_cast<double>( q ) /
                             std::max<Eigen::Index>( points - 1, 1 );
            for ( Eigen::Index i = 0; i < dofs; ++i )
            {
                testValues( i, q ) = 0.2 + 0.03*i + 0.002*q;
                trialValues( i, q ) = 0.3 - 0.01*i + 0.001*q;
            }
        }
    }

    /** @brief Start evaluation on a new element. */
    void beginElement() const noexcept {}

    /** @return representative composed variable coefficient at one point. */
    double coefficientValue( std::uint32_t q ) const noexcept
    {
        auto const x = coefficient[q];
        auto const indicator = x > 0.5 ? 1. : 0.;
        return ( 1. + x )*( 1. + x ) + indicator*( 0.25 + x*x );
    }

    /** @return one coefficient/test/trial product. */
    double evalijq( std::uint32_t i, std::uint32_t j,
                    std::uint16_t /*c1*/, std::uint16_t /*c2*/,
                    std::uint32_t q ) const noexcept
    {
        return coefficientValue( q )*testValues( i, q )*trialValues( j, q );
    }

    Eigen::VectorXd coefficient;
    table_type testValues;
    table_type trialValues;
};

/** @brief Eigen quadrature-block realization of the same mass-like tensor. */
struct batched_mass_tensor : pointwise_mass_tensor
{
    using pointwise_mass_tensor::pointwise_mass_tensor;

    /** @brief Invalidate the coefficient cache for the next element. */
    void beginElement() const noexcept { coefficientBatchValid = false; }

    /** @brief Evaluate one local matrix entry over every quadrature point. */
    void evaluateQuadratureBatch(
        quadrature::evaluation_request const& request,
        Eigen::Ref<Eigen::VectorXd> values ) const
    {
        if ( !coefficientBatchValid )
        {
            coefficientBatch.resize( coefficient.size() );
            coefficientBatch.array() =
                ( 1. + coefficient.array() ).square() +
                ( coefficient.array() > 0.5 )
                    .select( 0.25 + coefficient.array().square(), 0. );
            coefficientBatchValid = true;
        }
        values.array() = coefficientBatch.array() *
                         testValues.row( request.test ).transpose().array() *
                         trialValues.row( request.trial ).transpose().array();
    }

    mutable Eigen::VectorXd coefficientBatch;
    mutable bool coefficientBatchValid = false;
};

/**
 * @brief Measure only expression evaluation and quadrature reduction.
 * @tparam Tensor pointwise or Eigen-batched tensor contract
 * @param state arguments are local DoF count and quadrature order
 */
template <typename Tensor>
void BM_QuadratureExpression( benchmark::State& state )
{
    auto const dofs = static_cast<Eigen::Index>( state.range( 0 ) );
    auto const quadratureOrder = static_cast<int>( state.range( 1 ) );
    PointSetQuadrature<Simplex<2>, double> quadratureRule( quadratureOrder );
    quadratureRule.update( unit_quadrature_geometry{} );
    Tensor tensor( dofs, quadratureRule.nPoints() );

    for ( [[maybe_unused]] auto _ : state )
    {
        tensor.beginElement();
        double localSum = 0.;
        for ( Eigen::Index i = 0; i < dofs; ++i )
            for ( Eigen::Index j = 0; j < dofs; ++j )
                localSum += quadratureRule( tensor, i, j, 0, 0 );
        benchmark::DoNotOptimize( localSum );
    }

    auto const valuesPerIteration = static_cast<std::int64_t>(
        dofs*dofs*quadratureRule.nPoints() );
    state.SetItemsProcessed( state.iterations()*valuesPerIteration );
    state.counters["local-dofs"] = static_cast<double>( dofs );
    state.counters["quadrature-points"] = quadratureRule.nPoints();
}

/** @brief Register representative local sizes and quadrature orders. */
void configureExpressionBenchmark(
    benchmark::internal::Benchmark* registration )
{
    registration->Unit( benchmark::kMicrosecond )
        ->Args( { 6, 4 } )
        ->Args( { 6, 8 } )
        ->Args( { 15, 8 } )
        ->Args( { 15, 12 } );
}

using mesh_type = Mesh<Simplex<2, 2>>;
using space_type = Pch_type<mesh_type, 2>;

/**
 * @brief Build a genuinely curved second-order triangular mesh of the disk.
 * @param resolution reciprocal target mesh size
 * @return distributed curved mesh
 */
std::shared_ptr<mesh_type> makeCurvedMesh( int resolution )
{
    auto const meshSize = 1. / static_cast<double>( resolution );
    return createGMSHMesh(
        _mesh = new mesh_type,
        _desc = domain( _name = "dsel-batching-circle-g2-" +
                                std::to_string( resolution ),
                        _shape = "ellipsoid", _dim = 2, _order = 2,
                        _xmin = -1., _xmax = 1., _ymin = -1., _ymax = 1.,
                        _h = meshSize, _straighten = false ),
        _h = meshSize, _straighten = false, _verbose = false );
}

/**
 * @brief Assemble one DSEL mass or stiffness matrix.
 *
 * The mathematical expression is independent of the geometry execution
 * strategy. `GEOMAP_HO` ensures that all G2 cells retain their high-order
 * mapping, while `_geometry_evaluation` selects only the equivalent local
 * geometry kernel.
 *
 * @tparam Kind mass or stiffness form
 * @tparam ContractionPath entrywise quadrature batches or one local tensor
 * @tparam EvaluationStrategy pointwise, batched, or automatic geometry kernel
 * @tparam FormType bilinear-form type
 * @tparam ElementType finite-element field type
 * @param form matrix-backed bilinear form to overwrite
 * @param mesh curved integration mesh
 * @param trial trial field used by the DSEL expression
 * @param test test field used by the DSEL expression
 * @param quadratureOrder integration order
 */
template <dsel_form_kind Kind, dsel_contraction_path ContractionPath,
          geometry::geometry_evaluation_strategy EvaluationStrategy,
          typename FormType, typename ElementType>
void assembleDselForm( FormType& form, std::shared_ptr<mesh_type> const& mesh,
                       ElementType const& trial, ElementType const& test,
                       int quadratureOrder )
{
    constexpr auto lowering =
        ContractionPath == dsel_contraction_path::local_tensor ?
            vf::local_form_lowering::enabled : vf::local_form_lowering::disabled;
    if constexpr ( Kind == dsel_form_kind::mass )
    {
        auto mass = idt( trial )*id( test );
        form = integrate(
            _range = elements( mesh ), _expr = mass,
            _local_form_lowering = lowering,
            _quad = quadratureOrder,
            _geomap = GeomapStrategyType::GEOMAP_HO,
            _geometry_evaluation = EvaluationStrategy );
    }
    else
    {
        auto stiffness = inner( gradt( trial ), grad( test ) );
        form = integrate(
            _range = elements( mesh ), _expr = stiffness,
            _local_form_lowering = lowering,
            _quad = quadratureOrder,
            _geomap = GeomapStrategyType::GEOMAP_HO,
            _geometry_evaluation = EvaluationStrategy );
    }
    form.close();
}

/**
 * @brief Check geometry and contraction paths through complete matrix norms.
 * @tparam Kind mass or stiffness form
 * @tparam ContractionPath entrywise quadrature batches or one local tensor
 * @param space finite-element space
 * @param mesh curved integration mesh
 * @param trial trial field
 * @param test test field
 * @param quadratureOrder integration order
 */
template <dsel_form_kind Kind, dsel_contraction_path ContractionPath,
          typename ElementType>
void checkAssemblyParity( std::shared_ptr<space_type> const& space,
                          std::shared_ptr<mesh_type> const& mesh,
                          ElementType const& trial, ElementType const& test,
                          int quadratureOrder )
{
    auto pointwise = form2( _test = space, _trial = space );
    auto batched = form2( _test = space, _trial = space );
    assembleDselForm<Kind, ContractionPath,
                     geometry::geometry_evaluation_strategy::pointwise>(
        pointwise, mesh, trial, test, quadratureOrder );
    assembleDselForm<Kind, ContractionPath,
                     geometry::geometry_evaluation_strategy::batched>(
        batched, mesh, trial, test, quadratureOrder );

    auto difference = batched.matrixPtr()->clone();
    difference->addMatrix( -1., pointwise.matrixPtr() );
    difference->close();
    auto const scale = std::max( 1., std::abs( pointwise.matrixPtr()->l1Norm() ) );
    CHECK_LE( difference->l1Norm() / scale, 1e-12 )
        << "pointwise and batched DSEL matrices disagree in l1 norm";
    CHECK_LE( difference->linftyNorm() / scale, 1e-12 )
        << "pointwise and batched DSEL matrices disagree in infinity norm";

    auto otherPath = form2( _test = space, _trial = space );
    constexpr auto otherContraction =
        ContractionPath == dsel_contraction_path::local_tensor ?
            dsel_contraction_path::entrywise : dsel_contraction_path::local_tensor;
    assembleDselForm<Kind, otherContraction,
                     geometry::geometry_evaluation_strategy::pointwise>(
        otherPath, mesh, trial, test, quadratureOrder );
    otherPath.matrixPtr()->addMatrix( -1., pointwise.matrixPtr() );
    otherPath.matrixPtr()->close();
    CHECK_LE( otherPath.matrixPtr()->l1Norm() / scale, 1e-12 )
        << "entrywise and local-tensor DSEL matrices disagree in l1 norm";
    CHECK_LE( otherPath.matrixPtr()->linftyNorm() / scale, 1e-12 )
        << "entrywise and local-tensor DSEL matrices disagree in infinity norm";
}

/**
 * @brief Benchmark complete DSEL assembly, including sparse insertion/close.
 *
 * Mesh, space, expression fields, and sparse graph construction are outside
 * the timed loop. Each iteration measures matrix zeroing, element geometry,
 * FE/expression quadrature, local tensor formation, sparse insertion, and
 * backend close—the end-to-end operation visible to a Feel++ application.
 *
 * @tparam Kind mass or stiffness form
 * @tparam ContractionPath entrywise quadrature batches or one local tensor
 * @tparam EvaluationStrategy pointwise, batched, or automatic geometry kernel
 * @param state Google Benchmark state; arguments are resolution and quadrature
 * order
 */
template <dsel_form_kind Kind, dsel_contraction_path ContractionPath,
          geometry::geometry_evaluation_strategy EvaluationStrategy>
void BM_DselAssembly( benchmark::State& state )
{
    auto mesh = makeCurvedMesh( state.range( 0 ) );
    auto space = space_type::New( _mesh = mesh );
    auto trial = space->element();
    auto test = space->element();
    auto form = form2( _test = space, _trial = space );
    auto const quadratureOrder = static_cast<int>( state.range( 1 ) );

    checkAssemblyParity<Kind, ContractionPath>(
        space, mesh, trial, test, quadratureOrder );

    for ( [[maybe_unused]] auto _ : state )
        assembleDselForm<Kind, ContractionPath, EvaluationStrategy>(
            form, mesh, trial, test, quadratureOrder );

    state.counters["elements"] = static_cast<double>( mesh->numElements() );
    state.counters["dofs"] = static_cast<double>( space->nLocalDof() );
    state.counters["quadrature-order"] = quadratureOrder;
    state.SetItemsProcessed( state.iterations() * mesh->numElements() );
}

/** @brief Register the common curved-mesh benchmark sizes. */
void configureDselBenchmark( benchmark::internal::Benchmark* registration )
{
    registration->Unit( benchmark::kMillisecond )
        ->Args( { 4, 4 } )
        ->Args( { 8, 4 } )
        ->Args( { 8, 8 } )
        ->Args( { 16, 8 } );
}

} // namespace

BENCHMARK_TEMPLATE( BM_QuadratureExpression, pointwise_mass_tensor )
    ->Name( "quadrature-expression/mass/pointwise" )
    ->Apply( configureExpressionBenchmark );
BENCHMARK_TEMPLATE( BM_QuadratureExpression, batched_mass_tensor )
    ->Name( "quadrature-expression/mass/batched" )
    ->Apply( configureExpressionBenchmark );

BENCHMARK_TEMPLATE(
    BM_DselAssembly, dsel_form_kind::mass,
    dsel_contraction_path::local_tensor,
    geometry::geometry_evaluation_strategy::pointwise )
    ->Name( "dsel/mass/2d/g2/p2/pointwise" )
    ->Apply( configureDselBenchmark );
BENCHMARK_TEMPLATE(
    BM_DselAssembly, dsel_form_kind::mass,
    dsel_contraction_path::local_tensor,
    geometry::geometry_evaluation_strategy::batched )
    ->Name( "dsel/mass/2d/g2/p2/batched" )
    ->Apply( configureDselBenchmark );
BENCHMARK_TEMPLATE(
    BM_DselAssembly, dsel_form_kind::mass,
    dsel_contraction_path::local_tensor,
    geometry::geometry_evaluation_strategy::automatic )
    ->Name( "dsel/mass/2d/g2/p2/automatic" )
    ->Apply( configureDselBenchmark );

BENCHMARK_TEMPLATE(
    BM_DselAssembly, dsel_form_kind::stiffness,
    dsel_contraction_path::local_tensor,
    geometry::geometry_evaluation_strategy::pointwise )
    ->Name( "dsel/stiffness/2d/g2/p2/pointwise" )
    ->Apply( configureDselBenchmark );
BENCHMARK_TEMPLATE(
    BM_DselAssembly, dsel_form_kind::stiffness,
    dsel_contraction_path::local_tensor,
    geometry::geometry_evaluation_strategy::batched )
    ->Name( "dsel/stiffness/2d/g2/p2/batched" )
    ->Apply( configureDselBenchmark );
BENCHMARK_TEMPLATE(
    BM_DselAssembly, dsel_form_kind::stiffness,
    dsel_contraction_path::local_tensor,
    geometry::geometry_evaluation_strategy::automatic )
    ->Name( "dsel/stiffness/2d/g2/p2/automatic" )
    ->Apply( configureDselBenchmark );

BENCHMARK_TEMPLATE(
    BM_DselAssembly, dsel_form_kind::mass,
    dsel_contraction_path::entrywise,
    geometry::geometry_evaluation_strategy::automatic )
    ->Name( "dsel/mass/2d/g2/p2/entrywise-batch" )
    ->Apply( configureDselBenchmark );
BENCHMARK_TEMPLATE(
    BM_DselAssembly, dsel_form_kind::mass,
    dsel_contraction_path::local_tensor,
    geometry::geometry_evaluation_strategy::automatic )
    ->Name( "dsel/mass/2d/g2/p2/local-tensor" )
    ->Apply( configureDselBenchmark );
BENCHMARK_TEMPLATE(
    BM_DselAssembly, dsel_form_kind::stiffness,
    dsel_contraction_path::entrywise,
    geometry::geometry_evaluation_strategy::automatic )
    ->Name( "dsel/stiffness/2d/g2/p2/entrywise-batch" )
    ->Apply( configureDselBenchmark );
BENCHMARK_TEMPLATE(
    BM_DselAssembly, dsel_form_kind::stiffness,
    dsel_contraction_path::local_tensor,
    geometry::geometry_evaluation_strategy::automatic )
    ->Name( "dsel/stiffness/2d/g2/p2/local-tensor" )
    ->Apply( configureDselBenchmark );

/** @brief Run the Feel++ DSEL integration benchmarks. */
int main( int argc, char** argv )
{
    Environment environment(
        _argc = argc, _argv = argv,
        _about = about( _name = "feelpp_bench_integrate",
                        _author = "Feel++ Consortium",
                        _email = "feelpp-devel@feelpp.org" ) );

    benchmark::Initialize( &argc, argv );
    benchmark::RunSpecifiedBenchmarks();
}
