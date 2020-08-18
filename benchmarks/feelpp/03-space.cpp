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
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
using namespace Feel;

template<typename SpaceT>
void BM_Space( benchmark::State& state )
{
    using namespace Feel;
    size_type nitems = 0;
    for (auto _ : state)
    {
        using mesh_t = typename SpaceT::mesh_type;
        const int Dim = mesh_t::nDim;

        state.PauseTiming();
        auto m = loadMesh( _mesh=new mesh_t,
                           _filename="feelpp-"+std::to_string(Dim)+"d-"+std::to_string(state.range(0))+".msh",
                           _savehdf5=false,_verbose=false,
                           _h=1./state.range(0) );
        state.ResumeTiming();
        auto Xh = SpaceT::New( _mesh = m,
                               _worldscomm = makeWorldsComm( 1, m->worldComm() ) );
        nitems = Xh->nLocalDof();
    }
    state.SetItemsProcessed( state.iterations()*nitems );
    state.SetLabel( std::to_string(SpaceT::nDim)+"D h=1/"+std::to_string(state.range(0)) + " ndofs=" + std::to_string(nitems) );
}

BENCHMARK_TEMPLATE(BM_Space,Pch_type<Mesh<Simplex<2>>,1>)->Unit(benchmark::kMillisecond)->Arg(4)->Arg(8)->Arg(16)->Arg(32)->Arg(64);
BENCHMARK_TEMPLATE(BM_Space,Pch_type<Mesh<Simplex<2>>,2>)->Unit(benchmark::kMillisecond)->Arg(4)->Arg(8)->Arg(16)->Arg(32)->Arg(64);
BENCHMARK_TEMPLATE(BM_Space,Pch_type<Mesh<Simplex<2>>,3>)->Unit(benchmark::kMillisecond)->Arg(4)->Arg(8)->Arg(16)->Arg(32)->Arg(64);
BENCHMARK_TEMPLATE(BM_Space,Pch_type<Mesh<Simplex<3>>,1>)->Unit(benchmark::kMillisecond)->Arg(4)->Arg(8)->Arg(16)->Arg(32);
BENCHMARK_TEMPLATE(BM_Space,Pch_type<Mesh<Simplex<3>>,2>)->Unit(benchmark::kMillisecond)->Arg(4)->Arg(8)->Arg(16)->Arg(32);
BENCHMARK_TEMPLATE(BM_Space,Pch_type<Mesh<Simplex<3>>,3>)->Unit(benchmark::kMillisecond)->Arg(4)->Arg(8)->Arg(16)->Arg(32);

BENCHMARK_TEMPLATE(BM_Space,Pdh_type<Mesh<Simplex<2>>,1>)->Unit(benchmark::kMillisecond)->Arg(4)->Arg(8)->Arg(16)->Arg(32);
BENCHMARK_TEMPLATE(BM_Space,Pdh_type<Mesh<Simplex<3>>,1>)->Unit(benchmark::kMillisecond)->Arg(4)->Arg(8)->Arg(16)->Arg(32);

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
