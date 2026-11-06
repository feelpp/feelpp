//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme
//! @date 25 Mar 2026
//!

#include <benchmark/benchmark.h>

#include <numbers>

#include <feel/feeldiscr/meshstructured.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;

namespace
{

enum class SelectorKind
{
    Canonical,
    Symmetric,
    Mandel
};

template <SelectorKind Kind>
const char*
selectorName()
{
    if constexpr ( Kind == SelectorKind::Canonical )
        return "delta";
    else if constexpr ( Kind == SelectorKind::Symmetric )
        return "symm_delta";
    else
        return "mandel_delta";
}

template <bool UseStructured, SelectorKind Kind>
void
BM_TensorBasisContraction( benchmark::State& state )
{
    using mesh_t = MeshStructured<Hypercube<2>>;
    constexpr auto invSqrt2 = 1.0/std::numbers::sqrt2_v<double>;

    auto mesh = std::make_shared<mesh_t>( nl::json{ { "Discretisation", { { "n_points", { state.range( 0 ), state.range( 0 ) } } } } } );
    mesh->components().set( size_type( MESH_UPDATE_FACES|MESH_UPDATE_EDGES ) );
    mesh->updateForUse();

    auto A = mat<2, 2>( Px() + cst( 1.0 ),
                        Py(),
                        Px()*Py(),
                        cst( 2.0 ) + Px() );

    auto denseDelta = mat<2, 2>( cst( 0.0 ), cst( 1.0 ),
                                 cst( 0.0 ), cst( 0.0 ) );
    auto denseSymm = mat<2, 2>( cst( 0.0 ), cst( 1.0 ),
                                cst( 1.0 ), cst( 0.0 ) );
    auto denseMandel = mat<2, 2>( cst( 0.0 ), cst( invSqrt2 ),
                                  cst( invSqrt2 ), cst( 0.0 ) );

    for ( auto _ : state )
    {
        double value = 0.0;

        if constexpr ( UseStructured && Kind == SelectorKind::Canonical )
            value = integrate( _range=elements( mesh ), _expr=inner( A, delta<2, 0, 1>() ) ).evaluate()( 0, 0 );
        else if constexpr ( UseStructured && Kind == SelectorKind::Symmetric )
            value = integrate( _range=elements( mesh ), _expr=inner( A, symm_delta<2, 0, 1>() ) ).evaluate()( 0, 0 );
        else if constexpr ( UseStructured && Kind == SelectorKind::Mandel )
            value = integrate( _range=elements( mesh ), _expr=inner( A, mandel_delta<2, 0, 1>() ) ).evaluate()( 0, 0 );
        else if constexpr ( Kind == SelectorKind::Canonical )
            value = integrate( _range=elements( mesh ), _expr=inner( A, denseDelta ) ).evaluate()( 0, 0 );
        else if constexpr ( Kind == SelectorKind::Symmetric )
            value = integrate( _range=elements( mesh ), _expr=inner( A, denseSymm ) ).evaluate()( 0, 0 );
        else
            value = integrate( _range=elements( mesh ), _expr=inner( A, denseMandel ) ).evaluate()( 0, 0 );

        benchmark::DoNotOptimize( value );
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed( state.iterations()*mesh->numElements() );
    state.SetLabel( std::string( UseStructured ? "structured" : "dense" ) + ":" + selectorName<Kind>() );
}

template <SelectorKind Kind>
void
BM_RuntimeTensorBasisContraction( benchmark::State& state )
{
    using mesh_t = MeshStructured<Hypercube<2>>;

    auto mesh = std::make_shared<mesh_t>( nl::json{ { "Discretisation", { { "n_points", { state.range( 0 ), state.range( 0 ) } } } } } );
    mesh->components().set( size_type( MESH_UPDATE_FACES|MESH_UPDATE_EDGES ) );
    mesh->updateForUse();

    auto A = mat<2, 2>( Px() + cst( 1.0 ),
                        Py(),
                        Px()*Py(),
                        cst( 2.0 ) + Px() );

    for ( auto _ : state )
    {
        double value = 0.0;

        if constexpr ( Kind == SelectorKind::Canonical )
            value = integrate( _range=elements( mesh ), _expr=inner( A, delta<2>( 0, 1 ) ) ).evaluate()( 0, 0 );
        else if constexpr ( Kind == SelectorKind::Symmetric )
            value = integrate( _range=elements( mesh ), _expr=inner( A, symm_delta<2>( 0, 1 ) ) ).evaluate()( 0, 0 );
        else
            value = integrate( _range=elements( mesh ), _expr=inner( A, mandel_delta<2>( 0, 1 ) ) ).evaluate()( 0, 0 );

        benchmark::DoNotOptimize( value );
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed( state.iterations()*mesh->numElements() );
    state.SetLabel( std::string( "runtime:" ) + selectorName<Kind>() );
}

template <bool UseStructured, SelectorKind Kind>
void
BM_TensorBasisVectorProduct( benchmark::State& state )
{
    using mesh_t = MeshStructured<Hypercube<2>>;
    constexpr auto invSqrt2 = 1.0/std::numbers::sqrt2_v<double>;

    auto mesh = std::make_shared<mesh_t>( nl::json{ { "Discretisation", { { "n_points", { state.range( 0 ), state.range( 0 ) } } } } } );
    mesh->components().set( size_type( MESH_UPDATE_FACES|MESH_UPDATE_EDGES ) );
    mesh->updateForUse();

    auto v = vec( Px() + 2.0*Py(),
                  cst( 1.0 ) - Py() );

    auto denseDelta = mat<2, 2>( cst( 0.0 ), cst( 1.0 ),
                                 cst( 0.0 ), cst( 0.0 ) );
    auto denseSymm = mat<2, 2>( cst( 0.0 ), cst( 1.0 ),
                                cst( 1.0 ), cst( 0.0 ) );
    auto denseMandel = mat<2, 2>( cst( 0.0 ), cst( invSqrt2 ),
                                  cst( invSqrt2 ), cst( 0.0 ) );

    for ( auto _ : state )
    {
        double value = 0.0;

        if constexpr ( UseStructured && Kind == SelectorKind::Canonical )
            value = integrate( _range=elements( mesh ), _expr=inner( delta<2, 0, 1>()*v, delta<2, 0, 1>()*v ) ).evaluate()( 0, 0 );
        else if constexpr ( UseStructured && Kind == SelectorKind::Symmetric )
            value = integrate( _range=elements( mesh ), _expr=inner( symm_delta<2, 0, 1>()*v, symm_delta<2, 0, 1>()*v ) ).evaluate()( 0, 0 );
        else if constexpr ( UseStructured && Kind == SelectorKind::Mandel )
            value = integrate( _range=elements( mesh ), _expr=inner( mandel_delta<2, 0, 1>()*v, mandel_delta<2, 0, 1>()*v ) ).evaluate()( 0, 0 );
        else if constexpr ( Kind == SelectorKind::Canonical )
            value = integrate( _range=elements( mesh ), _expr=inner( denseDelta*v, denseDelta*v ) ).evaluate()( 0, 0 );
        else if constexpr ( Kind == SelectorKind::Symmetric )
            value = integrate( _range=elements( mesh ), _expr=inner( denseSymm*v, denseSymm*v ) ).evaluate()( 0, 0 );
        else
            value = integrate( _range=elements( mesh ), _expr=inner( denseMandel*v, denseMandel*v ) ).evaluate()( 0, 0 );

        benchmark::DoNotOptimize( value );
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed( state.iterations()*mesh->numElements() );
    state.SetLabel( std::string( UseStructured ? "structured" : "dense" ) + ":vector:" + selectorName<Kind>() );
}

template <bool UseStructured, SelectorKind Kind>
void
BM_TensorBasisRightMatrixProduct( benchmark::State& state )
{
    using mesh_t = MeshStructured<Hypercube<2>>;
    constexpr auto invSqrt2 = 1.0/std::numbers::sqrt2_v<double>;

    auto mesh = std::make_shared<mesh_t>( nl::json{ { "Discretisation", { { "n_points", { state.range( 0 ), state.range( 0 ) } } } } } );
    mesh->components().set( size_type( MESH_UPDATE_FACES|MESH_UPDATE_EDGES ) );
    mesh->updateForUse();

    auto A = mat<2, 2>( Px() + cst( 1.0 ),
                        Py(),
                        Px()*Py(),
                        cst( 2.0 ) + Px() );

    auto denseDelta = mat<2, 2>( cst( 0.0 ), cst( 1.0 ),
                                 cst( 0.0 ), cst( 0.0 ) );
    auto denseSymm = mat<2, 2>( cst( 0.0 ), cst( 1.0 ),
                                cst( 1.0 ), cst( 0.0 ) );
    auto denseMandel = mat<2, 2>( cst( 0.0 ), cst( invSqrt2 ),
                                  cst( invSqrt2 ), cst( 0.0 ) );

    for ( auto _ : state )
    {
        double value = 0.0;

        if constexpr ( UseStructured && Kind == SelectorKind::Canonical )
            value = integrate( _range=elements( mesh ), _expr=inner( A*delta<2, 0, 1>(), A*delta<2, 0, 1>() ) ).evaluate()( 0, 0 );
        else if constexpr ( UseStructured && Kind == SelectorKind::Symmetric )
            value = integrate( _range=elements( mesh ), _expr=inner( A*symm_delta<2, 0, 1>(), A*symm_delta<2, 0, 1>() ) ).evaluate()( 0, 0 );
        else if constexpr ( UseStructured && Kind == SelectorKind::Mandel )
            value = integrate( _range=elements( mesh ), _expr=inner( A*mandel_delta<2, 0, 1>(), A*mandel_delta<2, 0, 1>() ) ).evaluate()( 0, 0 );
        else if constexpr ( Kind == SelectorKind::Canonical )
            value = integrate( _range=elements( mesh ), _expr=inner( A*denseDelta, A*denseDelta ) ).evaluate()( 0, 0 );
        else if constexpr ( Kind == SelectorKind::Symmetric )
            value = integrate( _range=elements( mesh ), _expr=inner( A*denseSymm, A*denseSymm ) ).evaluate()( 0, 0 );
        else
            value = integrate( _range=elements( mesh ), _expr=inner( A*denseMandel, A*denseMandel ) ).evaluate()( 0, 0 );

        benchmark::DoNotOptimize( value );
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed( state.iterations()*mesh->numElements() );
    state.SetLabel( std::string( UseStructured ? "structured" : "dense" ) + ":right-mat:" + selectorName<Kind>() );
}

template <bool UseStructured, SelectorKind Kind>
void
BM_TensorBasisLeftMatrixProduct( benchmark::State& state )
{
    using mesh_t = MeshStructured<Hypercube<2>>;
    constexpr auto invSqrt2 = 1.0/std::numbers::sqrt2_v<double>;

    auto mesh = std::make_shared<mesh_t>( nl::json{ { "Discretisation", { { "n_points", { state.range( 0 ), state.range( 0 ) } } } } } );
    mesh->components().set( size_type( MESH_UPDATE_FACES|MESH_UPDATE_EDGES ) );
    mesh->updateForUse();

    auto A = mat<2, 2>( Px() + cst( 1.0 ),
                        Py(),
                        Px()*Py(),
                        cst( 2.0 ) + Px() );

    auto denseDelta = mat<2, 2>( cst( 0.0 ), cst( 1.0 ),
                                 cst( 0.0 ), cst( 0.0 ) );
    auto denseSymm = mat<2, 2>( cst( 0.0 ), cst( 1.0 ),
                                cst( 1.0 ), cst( 0.0 ) );
    auto denseMandel = mat<2, 2>( cst( 0.0 ), cst( invSqrt2 ),
                                  cst( invSqrt2 ), cst( 0.0 ) );

    for ( auto _ : state )
    {
        double value = 0.0;

        if constexpr ( UseStructured && Kind == SelectorKind::Canonical )
            value = integrate( _range=elements( mesh ), _expr=inner( delta<2, 0, 1>()*A, delta<2, 0, 1>()*A ) ).evaluate()( 0, 0 );
        else if constexpr ( UseStructured && Kind == SelectorKind::Symmetric )
            value = integrate( _range=elements( mesh ), _expr=inner( symm_delta<2, 0, 1>()*A, symm_delta<2, 0, 1>()*A ) ).evaluate()( 0, 0 );
        else if constexpr ( UseStructured && Kind == SelectorKind::Mandel )
            value = integrate( _range=elements( mesh ), _expr=inner( mandel_delta<2, 0, 1>()*A, mandel_delta<2, 0, 1>()*A ) ).evaluate()( 0, 0 );
        else if constexpr ( Kind == SelectorKind::Canonical )
            value = integrate( _range=elements( mesh ), _expr=inner( denseDelta*A, denseDelta*A ) ).evaluate()( 0, 0 );
        else if constexpr ( Kind == SelectorKind::Symmetric )
            value = integrate( _range=elements( mesh ), _expr=inner( denseSymm*A, denseSymm*A ) ).evaluate()( 0, 0 );
        else
            value = integrate( _range=elements( mesh ), _expr=inner( denseMandel*A, denseMandel*A ) ).evaluate()( 0, 0 );

        benchmark::DoNotOptimize( value );
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed( state.iterations()*mesh->numElements() );
    state.SetLabel( std::string( UseStructured ? "structured" : "dense" ) + ":left-mat:" + selectorName<Kind>() );
}

} // namespace

BENCHMARK_TEMPLATE( BM_TensorBasisContraction, false, SelectorKind::Canonical )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_TensorBasisContraction, true, SelectorKind::Canonical )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_RuntimeTensorBasisContraction, SelectorKind::Canonical )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_TensorBasisContraction, false, SelectorKind::Symmetric )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_TensorBasisContraction, true, SelectorKind::Symmetric )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_RuntimeTensorBasisContraction, SelectorKind::Symmetric )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_TensorBasisContraction, false, SelectorKind::Mandel )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_TensorBasisContraction, true, SelectorKind::Mandel )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_RuntimeTensorBasisContraction, SelectorKind::Mandel )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_TensorBasisVectorProduct, false, SelectorKind::Canonical )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_TensorBasisVectorProduct, true, SelectorKind::Canonical )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_TensorBasisVectorProduct, false, SelectorKind::Symmetric )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_TensorBasisVectorProduct, true, SelectorKind::Symmetric )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_TensorBasisVectorProduct, false, SelectorKind::Mandel )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_TensorBasisVectorProduct, true, SelectorKind::Mandel )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_TensorBasisRightMatrixProduct, false, SelectorKind::Canonical )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_TensorBasisRightMatrixProduct, true, SelectorKind::Canonical )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_TensorBasisRightMatrixProduct, false, SelectorKind::Symmetric )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_TensorBasisRightMatrixProduct, true, SelectorKind::Symmetric )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_TensorBasisRightMatrixProduct, false, SelectorKind::Mandel )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_TensorBasisRightMatrixProduct, true, SelectorKind::Mandel )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_TensorBasisLeftMatrixProduct, false, SelectorKind::Canonical )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_TensorBasisLeftMatrixProduct, true, SelectorKind::Canonical )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_TensorBasisLeftMatrixProduct, false, SelectorKind::Symmetric )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_TensorBasisLeftMatrixProduct, true, SelectorKind::Symmetric )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_TensorBasisLeftMatrixProduct, false, SelectorKind::Mandel )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_TensorBasisLeftMatrixProduct, true, SelectorKind::Mandel )->Unit( benchmark::kMillisecond )->Arg( 64 )->Arg( 128 );

int main( int argc, char** argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _about=about( _name="feelpp_bench_tensor_basis",
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ) );

    benchmark::Initialize( &argc, argv );
    benchmark::RunSpecifiedBenchmarks();
}
