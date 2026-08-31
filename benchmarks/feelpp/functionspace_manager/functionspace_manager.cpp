/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-

   SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

   SPDX-FileCopyrightText: 2026 University of Strasbourg

   SPDX-License-Identifier: LGPL-3.0-or-later
*/

#include <cstdint>
#include <set>
#include <string>
#include <vector>

#include <benchmark/benchmark.h>

#include <feel/feeldiscr/functionspacebuildinstrumentation.hpp>
#include <feel/feeldiscr/functionspacemanager.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/unitcube.hpp>
#include <feel/feelfilters/unitsquare.hpp>

namespace Feel
{
namespace
{

template <int Dim>
auto makeUnitMesh( std::int64_t subdivisions )
{
    auto const h = 1.0 / static_cast<double>( subdivisions );
    if constexpr ( Dim == 2 )
        return unitSquare( h );
    else
        return unitCube( h );
}

template <int Dim, typename Factory>
void benchmarkReusedWholeMeshFactory( benchmark::State& state, std::string const& spaceName,
                                      Factory&& factory )
{
    auto mesh = makeUnitMesh<Dim>( state.range( 0 ) );
    auto& manager = FunctionSpaceManager::instance();
    manager.clear( mesh );
    auto const memoryBefore = Environment::logMemoryUsage( spaceName + " before reused construction" );

    bool const instrumentationWasEnabled = FunctionSpaceBuildInstrumentation::enabled();
    FunctionSpaceBuildInstrumentation::reset();
    FunctionSpaceBuildInstrumentation::setEnabled( true );

    using space_ptrtype = decltype( factory( mesh ) );
    std::vector<space_ptrtype> spaces;
    std::set<void const*> distinctSpaces;
    double lastResidentMemory = memoryBefore.memory_usage;
    std::size_t request = 0;

    for ( auto _ : state )
    {
        auto space = factory( mesh );
        benchmark::DoNotOptimize( space.get() );

        state.PauseTiming();
        spaces.push_back( space );
        distinctSpaces.insert( space.get() );
        ++request;
        auto const memory = Environment::logMemoryUsage(
            spaceName + " after reused request " + std::to_string( request ) );
        lastResidentMemory = memory.memory_usage;
        state.counters["rss_request_" + std::to_string( request ) + "_bytes"] =
            memory.memory_usage;
        state.ResumeTiming();
    }

    auto const counts = FunctionSpaceBuildInstrumentation::counts();
    auto const stats = manager.stats( mesh );
    FunctionSpaceBuildInstrumentation::setEnabled( instrumentationWasEnabled );

    state.counters["function_space_constructions"] = counts.functionSpaceConstructions;
    state.counters["dof_table_builds"] = counts.dofTableBuilds;
    state.counters["distinct_space_pointers"] = distinctSpaces.size();
    state.counters["manager_hits"] = stats.hits;
    state.counters["manager_misses"] = stats.misses;
    state.counters["manager_retained_entries"] = stats.retainedEntries;
    state.counters["manager_retained_dofs"] = stats.retainedDofs;
    state.counters["rss_delta_bytes"] = lastResidentMemory - memoryBefore.memory_usage;
    state.counters["local_dofs"] = spaces.empty() ? 0 : spaces.front()->nLocalDof();
    state.counters["local_mesh_elements"] = mesh->numElements();
    state.counters["global_mesh_elements"] = mesh->numGlobalElements();
    state.counters["global_mesh_points"] = mesh->numGlobalPoints();
    state.SetLabel( spaceName + " reused " + std::to_string( Dim ) + "D h=1/" +
                    std::to_string( state.range( 0 ) ) );

    auto const expectedHits = static_cast<std::uint64_t>( state.iterations() ) - 1;
    if ( counts.functionSpaceConstructions != 1 || counts.dofTableBuilds != 1 ||
         distinctSpaces.size() != 1 || stats.misses != 1 || stats.hits != expectedHits )
        state.SkipWithError( "function-space manager did not canonicalize repeated requests" );

    manager.clear( mesh );
}

void BM_Pch2_2D_Reuse( benchmark::State& state )
{
    benchmarkReusedWholeMeshFactory<2>(
        state, "Pch2", []( auto const& mesh )
        {
            return Pch<2>(
                _mesh=mesh,
                _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
        } );
}

void BM_Pdh1_2D_Reuse( benchmark::State& state )
{
    benchmarkReusedWholeMeshFactory<2>(
        state, "Pdh1", []( auto const& mesh )
        {
            return Pdh<1>(
                _mesh=mesh,
                _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
        } );
}

void BM_Pch2_3D_Reuse( benchmark::State& state )
{
    benchmarkReusedWholeMeshFactory<3>(
        state, "Pch2", []( auto const& mesh )
        {
            return Pch<2>(
                _mesh=mesh,
                _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
        } );
}

void BM_Pdh1_3D_Reuse( benchmark::State& state )
{
    benchmarkReusedWholeMeshFactory<3>(
        state, "Pdh1", []( auto const& mesh )
        {
            return Pdh<1>(
                _mesh=mesh,
                _fspace_reuse_policy=FunctionSpaceReusePolicy::reuse );
        } );
}

BENCHMARK( BM_Pch2_2D_Reuse )->Unit( benchmark::kMillisecond )->Arg( 16 )->Iterations( 8 );
BENCHMARK( BM_Pdh1_2D_Reuse )->Unit( benchmark::kMillisecond )->Arg( 16 )->Iterations( 8 );
BENCHMARK( BM_Pch2_3D_Reuse )->Unit( benchmark::kMillisecond )->Arg( 8 )->Arg( 12 )->Iterations( 8 );
BENCHMARK( BM_Pdh1_3D_Reuse )->Unit( benchmark::kMillisecond )->Arg( 8 )->Arg( 12 )->Iterations( 8 );
BENCHMARK( BM_Pch2_3D_Reuse )
    ->Name( "BM_Pch2_3D_Large_Reuse" )
    ->Unit( benchmark::kMillisecond )
    ->Arg( 50 )
    ->Arg( 100 )
    ->Iterations( 8 );

} // namespace
} // namespace Feel

int main( int argc, char** argv )
{
    Feel::Environment env(
        Feel::_argc = argc,
        Feel::_argv = argv,
        Feel::_about = Feel::about(
            Feel::_name = "feelpp_bench_functionspace_manager",
            Feel::_author = "Feel++ Consortium",
            Feel::_email = "feelpp-devel@feelpp.org" ) );

    benchmark::Initialize( &argc, argv );
    benchmark::RunSpecifiedBenchmarks();
}
