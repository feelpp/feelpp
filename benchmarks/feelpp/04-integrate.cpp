/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2019-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 07 Apr 2019

#include <benchmark/benchmark.h>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/meshstructured.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;
template<typename SpaceT>
void BM_Space( benchmark::State& state )
{
    using namespace Feel;
    size_type nitems = 0;
    for (auto _ : state)
    {
        using mesh_t = typename SpaceT::mesh_type;
        using mesh_ptr_t = typename SpaceT::mesh_ptrtype;
        const int Dim = mesh_t::nDim;

        state.PauseTiming();

        mesh_ptr_t m;
        if constexpr ( is_mesh_structured_v<mesh_t> )
        {
            typename mesh_t::index_type nx_ = state.range( 0 );
            typename mesh_t::index_type ny_ = state.range( 1 );
            //m = std::make_shared<mesh_t>( {nx_, ny_} );
            m = std::make_shared<mesh_t>( nl::json{ { "Discretisation", {  { "n_points", { nx_,ny_ } } } } } );
            m->components().set( size_type( MESH_UPDATE_FACES|MESH_UPDATE_EDGES ) );
            m->updateForUse();
        }
        else
        {
            m = loadMesh( _mesh = new mesh_t,
                          _filename = "feelpp-" + std::to_string( Dim ) + "d-" + std::to_string( state.range( 0 ) ) + ".msh",
                          _savehdf5 = false, _verbose = false,
                          _h = 1. / state.range( 0 ) );
        }

        auto Xh = SpaceT::New( _mesh = m );
        nitems = Xh->nLocalDof();
        auto a = form2( _test=Xh, _trial=Xh );
        auto v = Xh->element();
        state.ResumeTiming();
        a=integrate( _range=elements(support(Xh)), _expr=inner(grad(v),gradt(v)));
        a.close();
    }
    state.SetItemsProcessed( state.iterations() * nitems );
    state.SetLabel( std::to_string( SpaceT::nDim ) + "D h=1/" + std::to_string( state.range( 0 ) ) + " ndofs=" + std::to_string( nitems ) );
}

#if 0
BENCHMARK_TEMPLATE(BM_Space,Pch_type<Mesh<Simplex<2>>,1>)->Unit(benchmark::kMillisecond)->Args({4,0})->Args({8,0})->Args({16,0})->Args({32,0})->Args({64,0});
BENCHMARK_TEMPLATE(BM_Space,Pch_type<Mesh<Simplex<2>>,2>)->Unit(benchmark::kMillisecond)->Args({4,0})->Args({8,0})->Args({16,0})->Args({32,0})->Args({64,0});
BENCHMARK_TEMPLATE(BM_Space,Pch_type<Mesh<Simplex<2>>,3>)->Unit(benchmark::kMillisecond)->Args({4,0})->Args({8,0})->Args({16,0})->Args({32,0})->Args({64,0});
BENCHMARK_TEMPLATE(BM_Space,Pch_type<Mesh<Simplex<3>>,1>)->Unit(benchmark::kMillisecond)->Args({4,0})->Args({8,0})->Args({16,0})->Args({32,0});
BENCHMARK_TEMPLATE(BM_Space,Pch_type<Mesh<Simplex<3>>,2>)->Unit(benchmark::kMillisecond)->Args({4,0})->Args({8,0})->Args({16,0})->Args({32,0});
BENCHMARK_TEMPLATE(BM_Space,Pch_type<Mesh<Simplex<3>>,3>)->Unit(benchmark::kMillisecond)->Args({4,0})->Args({8,0})->Args({16,0})->Args({32,0});

BENCHMARK_TEMPLATE(BM_Space,Pdh_type<Mesh<Simplex<2>>,1>)->Unit(benchmark::kMillisecond)->Args({4,0})->Args({8,0})->Args({16,0})->Args({32,0});
BENCHMARK_TEMPLATE(BM_Space,Pdh_type<Mesh<Simplex<3>>,1>)->Unit(benchmark::kMillisecond)->Args({4,0})->Args({8,0})->Args({16,0})->Args({32,0});
#endif
BENCHMARK_TEMPLATE(BM_Space,Pch_type<MeshStructured<Hypercube<2>>,1>)->Unit(benchmark::kMillisecond)
                    ->Args({4,4})->Args({8,8})->Args({16,16})->Args({32,32})->Args({64,64})
                    ->Args({128,128})->Args({256,256})->Args({512,512})->Args({1024,1024});
//                    ->Args({1920,1024})->Args({2048,2048});

int main(int argc, char** argv)
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="feelpp_bench_space",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    //benchmark::FLAGS_benchmark_out = "feelpp_bench_mesh.json";
    benchmark::Initialize(&argc, argv);
    //benchmark::ConsoleReporter cr;
    //benchmark::JSONReporter jr;
    benchmark::RunSpecifiedBenchmarks();
}
