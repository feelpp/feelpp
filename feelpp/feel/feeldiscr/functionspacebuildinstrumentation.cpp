/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

   SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

   SPDX-FileCopyrightText: 2026 University of Strasbourg

   SPDX-License-Identifier: LGPL-3.0-or-later
*/

#include <atomic>

#include <feel/feeldiscr/functionspacebuildinstrumentation.hpp>

namespace Feel
{
namespace
{
std::atomic<bool> S_enabled{ false };
std::atomic<std::uint64_t> S_functionSpaceConstructions{ 0 };
std::atomic<std::uint64_t> S_dofTableBuilds{ 0 };
} // namespace

void FunctionSpaceBuildInstrumentation::setEnabled( bool enabled ) noexcept
{
    S_enabled.store( enabled, std::memory_order_relaxed );
}

bool FunctionSpaceBuildInstrumentation::enabled() noexcept
{
    return S_enabled.load( std::memory_order_relaxed );
}

void FunctionSpaceBuildInstrumentation::reset() noexcept
{
    S_functionSpaceConstructions.store( 0, std::memory_order_relaxed );
    S_dofTableBuilds.store( 0, std::memory_order_relaxed );
}

FunctionSpaceBuildCounts
FunctionSpaceBuildInstrumentation::counts() noexcept
{
    return {
        S_functionSpaceConstructions.load( std::memory_order_relaxed ),
        S_dofTableBuilds.load( std::memory_order_relaxed ) };
}

void FunctionSpaceBuildInstrumentation::recordFunctionSpaceConstruction() noexcept
{
    if ( S_enabled.load( std::memory_order_relaxed ) )
        S_functionSpaceConstructions.fetch_add( 1, std::memory_order_relaxed );
}

void FunctionSpaceBuildInstrumentation::recordDofTableBuild() noexcept
{
    if ( S_enabled.load( std::memory_order_relaxed ) )
        S_dofTableBuilds.fetch_add( 1, std::memory_order_relaxed );
}

} // namespace Feel
