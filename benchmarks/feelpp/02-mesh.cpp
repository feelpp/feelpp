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
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 07 Apr 2019
//! @copyright 2019 Feel++ Consortium
//!

#include <benchmark/benchmark.h>
#include <feel/feelfilters/loadmesh.hpp>

using namespace Feel;
// loadmesh
void BM_Loadmesh3D( benchmark::State& state )
{
    using namespace Feel;
    size_type nitems;
    for (auto _ : state)
    {
        using mesh_t = Mesh<Simplex<3>>;
        auto m = loadMesh( _mesh=new mesh_t,
                           _filename="hypercube-"+std::to_string(state.range(0))+".msh",
                           _savehdf5=false,_verbose=false,
                           _update=0, _h=1./state.range(0) );
        nitems=m->numElements();
    }
    state.SetItemsProcessed( nitems );
}
BENCHMARK(BM_Loadmesh3D)->Unit(benchmark::kMillisecond)->Arg(4)->Arg(8)->Arg(16);//->Arg(32);

// loadmesh
template<typename MeshT>
void BM_updateForUse( benchmark::State& state )
{
    using namespace Feel;
    using mesh_t = MeshT;
    const int nDim = mesh_t::nDim;

    size_type nitems;
    for (auto _ : state)
    {
        state.PauseTiming();
        auto m = loadMesh( _mesh=new mesh_t,
                           _filename="feel-"+std::to_string(nDim)+"d-"+std::to_string(state.range(0))+".msh",
                           _savehdf5=false,_verbose=false,
                           _update=0, _h=1./state.range(0) );
        state.ResumeTiming();
        m->components().reset();
        m->components().set( state.range(1) );
        m->updateForUse();
        nitems=m->numElements();
    }
    state.SetItemsProcessed( state.iterations()*nitems );
    state.SetLabel( std::to_string(nDim)+"D h=1/"+std::to_string(state.range(0))+" update="+std::to_string(state.range(1)) + " nelts=" + std::to_string(nitems) );
}
BENCHMARK_TEMPLATE(BM_updateForUse,Mesh<Simplex<2>>)->Unit(benchmark::kMillisecond)
->Args({4,0})
->Args({4,MESH_NO_UPDATE_MEASURES})
->Args({4,MESH_UPDATE_FACES})
->Args({4,MESH_UPDATE_EDGES})
->Args({4,MESH_UPDATE_ELEMENTS_ADJACENCY})
->Args({4,MESH_UPDATE_FACES_MINIMAL})
->Args({8,0})
->Args({8,MESH_NO_UPDATE_MEASURES})
->Args({8,MESH_UPDATE_FACES})
->Args({8,MESH_UPDATE_EDGES})
->Args({8,MESH_UPDATE_ELEMENTS_ADJACENCY})
->Args({8,MESH_UPDATE_FACES_MINIMAL})
->Args({16,0})
->Args({16,MESH_NO_UPDATE_MEASURES})
->Args({16,MESH_UPDATE_FACES})
->Args({16,MESH_UPDATE_EDGES})
->Args({16,MESH_UPDATE_ELEMENTS_ADJACENCY})
->Args({16,MESH_UPDATE_FACES_MINIMAL})
->Args({32,0})
->Args({32,MESH_NO_UPDATE_MEASURES})
->Args({32,MESH_UPDATE_FACES})
->Args({32,MESH_UPDATE_EDGES})
->Args({32,MESH_UPDATE_ELEMENTS_ADJACENCY})
->Args({32,MESH_UPDATE_FACES_MINIMAL});
BENCHMARK_TEMPLATE(BM_updateForUse,Mesh<Simplex<3>>)->Unit(benchmark::kMillisecond)
->Args({4,0})
->Args({4,MESH_NO_UPDATE_MEASURES})
->Args({4,MESH_UPDATE_FACES})
->Args({4,MESH_UPDATE_EDGES})
->Args({4,MESH_UPDATE_ELEMENTS_ADJACENCY})
->Args({4,MESH_UPDATE_FACES_MINIMAL})
->Args({8,0})
->Args({8,MESH_NO_UPDATE_MEASURES})
->Args({8,MESH_UPDATE_FACES})
->Args({8,MESH_UPDATE_EDGES})
->Args({8,MESH_UPDATE_ELEMENTS_ADJACENCY})
->Args({8,MESH_UPDATE_FACES_MINIMAL})
->Args({16,0})
->Args({16,MESH_NO_UPDATE_MEASURES})
->Args({16,MESH_UPDATE_FACES})
->Args({16,MESH_UPDATE_EDGES})
->Args({16,MESH_UPDATE_ELEMENTS_ADJACENCY})
->Args({16,MESH_UPDATE_FACES_MINIMAL});


int main(int argc, char** argv)
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="feelpp_bench_mesh",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    //benchmark::FLAGS_benchmark_out = "feelpp_bench_mesh.json";
    benchmark::Initialize(&argc, argv);
    //benchmark::ConsoleReporter cr;
    //benchmark::JSONReporter jr;
    benchmark::RunSpecifiedBenchmarks();
}
