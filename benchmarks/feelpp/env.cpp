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
#include <feel/feelcore/environment.hpp>

// loadmesh

static void BM_env(benchmark::State& st)
{
    using namespace Feel;
    for (auto _ : st)
    {
    }
};
BENCHMARK(BM_env)->Arg(8)->Arg(64)->Arg(512)->Arg(1<<10)->Arg(8<<10);
//benchmark::internal::Benchmark*  dummy1 = benchmark::RegisterBenchmark(
//    "BM_env_registration", BM_env);

int main(int argc, char** argv)
{
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
}
