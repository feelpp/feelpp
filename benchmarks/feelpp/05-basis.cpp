/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#include <benchmark/benchmark.h>

#include <feel/feeldiscr/meshstructured.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;

template <bool UseBasisProxy>
void
BM_ElasticityLikeAssembly( benchmark::State& state )
{
    using mesh_t = MeshStructured<Hypercube<2>>;

    auto mesh = std::make_shared<mesh_t>( nl::json{ { "Discretisation", { { "n_points", { state.range( 0 ), state.range( 0 ) } } } } } );
    mesh->components().set( size_type( MESH_UPDATE_FACES|MESH_UPDATE_EDGES ) );
    mesh->updateForUse();

    auto Xh = Pchv<1>( mesh );
    auto ndofs = Xh->nLocalDof();

    for ( auto _ : state )
    {
        auto a = form2( _test=Xh, _trial=Xh );

        if constexpr ( UseBasisProxy )
        {
            auto u = trial( Xh );
            auto v = test( Xh );
            a = integrate( _range=elements( mesh ),
                           _expr=inner( symm_grad( u ), symm_grad( v ) ) + inner( u, v ) );
        }
        else
        {
            auto u = Xh->element( "u" );
            auto v = Xh->element( "v" );
            a = integrate( _range=elements( mesh ),
                           _expr=inner( sym( gradt( u ) ), sym( grad( v ) ) ) + inner( idt( u ), id( v ) ) );
        }

        a.close();
        benchmark::DoNotOptimize( a.matrixPtr() );
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed( state.iterations()*ndofs );
    state.SetLabel( std::string( UseBasisProxy? "basis":"legacy" ) +
                    " n=" + std::to_string( state.range( 0 ) ) +
                    " ndofs=" + std::to_string( ndofs ) );
}

BENCHMARK_TEMPLATE( BM_ElasticityLikeAssembly, false )->Unit( benchmark::kMillisecond )->Arg( 32 )->Arg( 64 )->Arg( 128 );
BENCHMARK_TEMPLATE( BM_ElasticityLikeAssembly, true )->Unit( benchmark::kMillisecond )->Arg( 32 )->Arg( 64 )->Arg( 128 );

int main( int argc, char** argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _about=about( _name="feelpp_bench_basis",
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ) );

    benchmark::Initialize( &argc, argv );
    benchmark::RunSpecifiedBenchmarks();
}
