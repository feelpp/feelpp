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
