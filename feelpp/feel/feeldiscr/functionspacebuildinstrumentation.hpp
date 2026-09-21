/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

   SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

   SPDX-FileCopyrightText: 2026 University of Strasbourg

   SPDX-License-Identifier: LGPL-3.0-or-later
*/

/**
 * @file functionspacebuildinstrumentation.hpp
 * @brief Opt-in construction counters for function-space tests and benchmarks.
 */

#ifndef FEELPP_DISCR_FUNCTIONSPACEBUILDINSTRUMENTATION_HPP
#define FEELPP_DISCR_FUNCTIONSPACEBUILDINSTRUMENTATION_HPP 1

#include <cstdint>

#include <feel/feelcore/feelmacros.hpp>

namespace Feel
{

/**
 * @brief Snapshot of function-space construction instrumentation counters.
 */
struct FunctionSpaceBuildCounts
{
    std::uint64_t functionSpaceConstructions = 0; ///< FunctionSpaceBase constructions observed while enabled.
    std::uint64_t dofTableBuilds = 0;             ///< DOF-table builds observed while enabled.
};

/**
 * @brief Runtime-disabled counters used by focused tests and benchmarks.
 *
 * Reset and enable/disable operations are expected to be performed while no
 * function-space construction is running. Recording is thread-safe.
 */
class FEELPP_EXPORT FunctionSpaceBuildInstrumentation
{
  public:
    /**
     * @brief Enable or disable construction event recording.
     * @param enabled whether subsequent events are recorded
     */
    static void setEnabled( bool enabled ) noexcept;

    /**
     * @brief Report whether construction event recording is active.
     * @return @c true when recording is enabled
     */
    static bool enabled() noexcept;

    /**
     * @brief Reset both counters to zero.
     *
     * Call while no function-space or DOF-table construction is running.
     */
    static void reset() noexcept;

    /**
     * @brief Read both construction counters.
     * @return an atomic snapshot of both process-local counters
     */
    static FunctionSpaceBuildCounts counts() noexcept;

    /** @brief Record one FunctionSpaceBase construction when enabled. */
    static void recordFunctionSpaceConstruction() noexcept;

    /** @brief Record one DOF-table build when enabled. */
    static void recordDofTableBuild() noexcept;
};

} // namespace Feel

#endif // FEELPP_DISCR_FUNCTIONSPACEBUILDINSTRUMENTATION_HPP
