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

#include <feel/feeldiscr/meshstructured.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;

namespace
{

template <bool UseManualExpansion>
void
BM_MatVecProduct( benchmark::State& state )
{
    using mesh_t = MeshStructured<Hypercube<2>>;

    auto mesh = std::make_shared<mesh_t>( nl::json{ { "Discretisation", { { "n_points", { state.range( 0 ), state.range( 0 ) } } } } } );
    mesh->components().set( size_type( MESH_UPDATE_FACES|MESH_UPDATE_EDGES ) );
    mesh->updateForUse();

    auto a00 = Px() + cst( 1.0 );
    auto a01 = Py();
    auto a10 = Px()*Py();
    auto a11 = cst( 2.0 ) + Px();
    auto A = mat<2, 2>( a00, a01,
                        a10, a11 );

    auto v0 = Px() + cst( 2.0 )*Py();
    auto v1 = cst( 1.0 ) - Py();
    auto v = vec( v0, v1 );

    auto expanded = vec( a00*v0 + a01*v1,
                         a10*v0 + a11*v1 );
    auto generic = A*v;

    for ( auto _ : state )
    {
        auto value = UseManualExpansion ?
            integrate( _range=elements( mesh ), _expr=inner( expanded, expanded ) ).evaluate()( 0, 0 ) :
            integrate( _range=elements( mesh ), _expr=inner( generic, generic ) ).evaluate()( 0, 0 );

        benchmark::DoNotOptimize( value );
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed( state.iterations()*mesh->numElements() );
    state.SetLabel( UseManualExpansion ? "manual:matvec" : "generic:matvec" );
}

template <bool UseManualExpansion>
void
BM_MatMatProduct( benchmark::State& state )
{
    using mesh_t = MeshStructured<Hypercube<2>>;

    auto mesh = std::make_shared<mesh_t>( nl::json{ { "Discretisation", { { "n_points", { state.range( 0 ), state.range( 0 ) } } } } } );
    mesh->components().set( size_type( MESH_UPDATE_FACES|MESH_UPDATE_EDGES ) );
    mesh->updateForUse();

    auto a00 = Px() + cst( 1.0 );
    auto a01 = Py();
    auto a10 = Px()*Py();
    auto a11 = cst( 2.0 ) + Px();
    auto A = mat<2, 2>( a00, a01,
                        a10, a11 );

    auto b00 = cst( 2.0 ) - Py();
    auto b01 = Px();
    auto b10 = cst( 1.0 ) + Py();
    auto b11 = cst( 3.0 ) + Px()*Py();
    auto B = mat<2, 2>( b00, b01,
                        b10, b11 );

    auto expanded = mat<2, 2>( a00*b00 + a01*b10, a00*b01 + a01*b11,
                               a10*b00 + a11*b10, a10*b01 + a11*b11 );
    auto generic = A*B;

    for ( auto _ : state )
    {
        auto value = UseManualExpansion ?
            integrate( _range=elements( mesh ), _expr=inner( expanded, expanded ) ).evaluate()( 0, 0 ) :
            integrate( _range=elements( mesh ), _expr=inner( generic, generic ) ).evaluate()( 0, 0 );

        benchmark::DoNotOptimize( value );
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed( state.iterations()*mesh->numElements() );
    state.SetLabel( UseManualExpansion ? "manual:matmat" : "generic:matmat" );
}

template <bool UseManualExpansion>
void
BM_MatHadamardProduct( benchmark::State& state )
{
    using mesh_t = MeshStructured<Hypercube<2>>;

    auto mesh = std::make_shared<mesh_t>( nl::json{ { "Discretisation", { { "n_points", { state.range( 0 ), state.range( 0 ) } } } } } );
    mesh->components().set( size_type( MESH_UPDATE_FACES|MESH_UPDATE_EDGES ) );
    mesh->updateForUse();

    auto a00 = Px() + cst( 1.0 );
    auto a01 = Py();
    auto a10 = Px()*Py();
    auto a11 = cst( 2.0 ) + Px();
    auto A = mat<2, 2>( a00, a01,
                        a10, a11 );

    auto b00 = cst( 2.0 ) - Py();
    auto b01 = Px();
    auto b10 = cst( 1.0 ) + Py();
    auto b11 = cst( 3.0 ) + Px()*Py();
    auto B = mat<2, 2>( b00, b01,
                        b10, b11 );

    auto expanded = mat<2, 2>( a00*b00, a01*b01,
                               a10*b10, a11*b11 );
    auto generic = hadamard( A, B );

    for ( auto _ : state )
    {
        auto value = UseManualExpansion ?
            integrate( _range=elements( mesh ), _expr=inner( expanded, expanded ) ).evaluate()( 0, 0 ) :
            integrate( _range=elements( mesh ), _expr=inner( generic, generic ) ).evaluate()( 0, 0 );

        benchmark::DoNotOptimize( value );
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed( state.iterations()*mesh->numElements() );
    state.SetLabel( UseManualExpansion ? "manual:hadamard-mat" : "explicit:hadamard-mat" );
}

template <bool UseManualExpansion>
void
BM_VecHadamardProduct( benchmark::State& state )
{
    using mesh_t = MeshStructured<Hypercube<2>>;

    auto mesh = std::make_shared<mesh_t>( nl::json{ { "Discretisation", { { "n_points", { state.range( 0 ), state.range( 0 ) } } } } } );
    mesh->components().set( size_type( MESH_UPDATE_FACES|MESH_UPDATE_EDGES ) );
    mesh->updateForUse();

    auto v0 = Px() + cst( 2.0 )*Py();
    auto v1 = cst( 1.0 ) - Py();
    auto v = vec( v0, v1 );

    auto w0 = cst( 3.0 ) - Px();
    auto w1 = cst( 2.0 ) + Px()*Py();
    auto w = vec( w0, w1 );

    auto expanded = vec( v0*w0,
                         v1*w1 );
    auto generic = hadamard( v, w );

    for ( auto _ : state )
    {
        auto value = UseManualExpansion ?
            integrate( _range=elements( mesh ), _expr=inner( expanded, expanded ) ).evaluate()( 0, 0 ) :
            integrate( _range=elements( mesh ), _expr=inner( generic, generic ) ).evaluate()( 0, 0 );

        benchmark::DoNotOptimize( value );
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed( state.iterations()*mesh->numElements() );
    state.SetLabel( UseManualExpansion ? "manual:hadamard-vec" : "explicit:hadamard-vec" );
}

template <int Mode>
void
BM_MatComponentAccess( benchmark::State& state )
{
    using mesh_t = MeshStructured<Hypercube<2>>;

    auto mesh = std::make_shared<mesh_t>( nl::json{ { "Discretisation", { { "n_points", { state.range( 0 ), state.range( 0 ) } } } } } );
    mesh->components().set( size_type( MESH_UPDATE_FACES|MESH_UPDATE_EDGES ) );
    mesh->updateForUse();

    auto a00 = Px() + cst( 1.0 );
    auto a01 = Py();
    auto a10 = Px()*Py();
    auto a11 = cst( 2.0 ) + Px();
    auto A = mat<2, 2>( a00, a01,
                        a10, a11 );

    for ( auto _ : state )
    {
        auto value =
            [&]() -> double
            {
                if constexpr ( Mode == 0 )
                    return integrate( _range=elements( mesh ), _expr=A( 0, 1 )*A( 0, 1 ) ).evaluate()( 0, 0 );
                else if constexpr ( Mode == 1 )
                    return integrate( _range=elements( mesh ), _expr=component<0,1>( A )*component<0,1>( A ) ).evaluate()( 0, 0 );
                else
                    return integrate( _range=elements( mesh ), _expr=a01*a01 ).evaluate()( 0, 0 );
            }();

        benchmark::DoNotOptimize( value );
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed( state.iterations()*mesh->numElements() );
    if constexpr ( Mode == 0 )
        state.SetLabel( "runtime:mat-component" );
    else if constexpr ( Mode == 1 )
        state.SetLabel( "static:mat-component" );
    else
        state.SetLabel( "manual:mat-entry" );
}

template <int Mode>
void
BM_VecComponentAccess( benchmark::State& state )
{
    using mesh_t = MeshStructured<Hypercube<2>>;

    auto mesh = std::make_shared<mesh_t>( nl::json{ { "Discretisation", { { "n_points", { state.range( 0 ), state.range( 0 ) } } } } } );
    mesh->components().set( size_type( MESH_UPDATE_FACES|MESH_UPDATE_EDGES ) );
    mesh->updateForUse();

    auto v0 = Px() + cst( 2.0 )*Py();
    auto v1 = cst( 1.0 ) - Py();
    auto v = vec( v0, v1 );

    for ( auto _ : state )
    {
        auto value =
            [&]() -> double
            {
                if constexpr ( Mode == 0 )
                    return integrate( _range=elements( mesh ), _expr=v( 1, 0 )*v( 1, 0 ) ).evaluate()( 0, 0 );
                else if constexpr ( Mode == 1 )
                    return integrate( _range=elements( mesh ), _expr=component<1,0>( v )*component<1,0>( v ) ).evaluate()( 0, 0 );
                else
                    return integrate( _range=elements( mesh ), _expr=v1*v1 ).evaluate()( 0, 0 );
            }();

        benchmark::DoNotOptimize( value );
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed( state.iterations()*mesh->numElements() );
    if constexpr ( Mode == 0 )
        state.SetLabel( "runtime:vec-component" );
    else if constexpr ( Mode == 1 )
        state.SetLabel( "static:vec-component" );
    else
        state.SetLabel( "manual:vec-entry" );
}

BENCHMARK_TEMPLATE( BM_MatVecProduct, false )->Unit( benchmark::kMillisecond )->Arg( 32 )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_MatVecProduct, true )->Unit( benchmark::kMillisecond )->Arg( 32 )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_MatMatProduct, false )->Unit( benchmark::kMillisecond )->Arg( 32 )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_MatMatProduct, true )->Unit( benchmark::kMillisecond )->Arg( 32 )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_MatHadamardProduct, false )->Unit( benchmark::kMillisecond )->Arg( 32 )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_MatHadamardProduct, true )->Unit( benchmark::kMillisecond )->Arg( 32 )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_VecHadamardProduct, false )->Unit( benchmark::kMillisecond )->Arg( 32 )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_VecHadamardProduct, true )->Unit( benchmark::kMillisecond )->Arg( 32 )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_MatComponentAccess, 0 )->Unit( benchmark::kMillisecond )->Arg( 32 )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_MatComponentAccess, 1 )->Unit( benchmark::kMillisecond )->Arg( 32 )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_MatComponentAccess, 2 )->Unit( benchmark::kMillisecond )->Arg( 32 )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_VecComponentAccess, 0 )->Unit( benchmark::kMillisecond )->Arg( 32 )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_VecComponentAccess, 1 )->Unit( benchmark::kMillisecond )->Arg( 32 )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_VecComponentAccess, 2 )->Unit( benchmark::kMillisecond )->Arg( 32 )->Arg( 64 )->Arg( 128 );

} // namespace

int main( int argc, char** argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _about=about( _name="feelpp_bench_matvec",
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ) );

    benchmark::Initialize( &argc, argv );
    benchmark::RunSpecifiedBenchmarks();
}
