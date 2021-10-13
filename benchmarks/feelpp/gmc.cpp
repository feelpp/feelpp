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

// gmc update
template<typename MeshT, size_type gmc_v>
void BM_gmc( benchmark::State& state )
{
    using namespace Feel;
    using mesh_t = MeshT;
    const int nDim = mesh_t::nDim;
    using im_type = im_t<typename mesh_t::element_type, double>;
    std::map<int,std::shared_ptr<mesh_t>> meshes;

    using gmc_t = typename mesh_t::gm_type::template Context<typename mesh_t::element_type>;
    std::map<int,std::shared_ptr<gmc_t>> gmcs;

    size_type nitems;
    for (auto _ : state)
    {
        state.PauseTiming();
        if ( meshes.find( state.range(0) ) == meshes.end() )
        {
            meshes[ state.range(0) ] = loadMesh( _mesh=new mesh_t,
                                                 //_filename="feel-"+std::to_string(nDim)+"d-"+std::to_string(state.range(0))+".msh",
                                                 _savehdf5=false,_verbose=false,
                                                 _update=1, _h=1./state.range(0) );
        }
        auto m = meshes[ state.range(0) ];

        if ( gmcs.find( state.range(1) ) == gmcs.end() )
        {
            auto geopc = m->gm()->preCompute( im_type(state.range(1)).points() );
            auto const& eltInit = m->beginElement()->second;
            gmcs[ state.range(1)] = m->gm()->template context<gmc_v>(eltInit,geopc);
        }
        auto gmc = gmcs[ state.range(1) ];
        state.ResumeTiming();
        for ( auto const& [eltId,elt] : m->elements() )
        {
            gmc->template update<gmc_v>( elt );
        }
        nitems=m->numElements()*gmc->pc()->nPoints();
    }
    state.SetItemsProcessed( state.iterations()*nitems );
    //state.SetLabel( std::to_string(nDim)+"D h=1/"+std::to_string(state.range(0))+" update="+std::to_string(state.range(1)) + " nelts=" + std::to_string(nitems) );
}
BENCHMARK_TEMPLATE(BM_gmc,Mesh<Simplex<2>>,0)->Unit(benchmark::kMillisecond)
->Args({4,1})->Args({4,2})->Args({4,4})
->Args({8,1})->Args({8,2})->Args({8,4})
;
BENCHMARK_TEMPLATE(BM_gmc,Mesh<Simplex<2>>,vm::JACOBIAN)->Unit(benchmark::kMillisecond)
->Args({4,1})->Args({4,2})->Args({4,4})
->Args({8,1})->Args({8,2})->Args({8,4})
;
BENCHMARK_TEMPLATE(BM_gmc,Mesh<Simplex<2>>,vm::JACOBIAN|vm::KB)->Unit(benchmark::kMillisecond)
->Args({4,1})->Args({4,2})->Args({4,4})
->Args({8,1})->Args({8,2})->Args({8,4})
;
BENCHMARK_TEMPLATE(BM_gmc,Mesh<Simplex<2>>,vm::JACOBIAN|vm::KB|vm::SECOND_DERIVATIVE)->Unit(benchmark::kMillisecond)
->Args({4,1})->Args({4,2})->Args({4,4})
->Args({8,1})->Args({8,2})->Args({8,4})
;


int main(int argc, char** argv)
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="feelpp_bench_gmc",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    //benchmark::FLAGS_benchmark_out = "feelpp_bench_mesh.json";
    benchmark::Initialize(&argc, argv);
    //benchmark::ConsoleReporter cr;
    //benchmark::JSONReporter jr;
    benchmark::RunSpecifiedBenchmarks();
}
