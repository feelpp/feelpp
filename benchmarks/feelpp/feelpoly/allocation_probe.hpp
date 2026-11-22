/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#pragma once

#include <cstdint>

namespace Feel::benchmarking::allocation_probe
{
/** C++ allocation requests; direct malloc calls inside libraries are excluded. */
struct AllocationSnapshot
{
    /** Number of instrumented C++ allocation requests. */
    std::uint64_t calls = 0;
    /** Sum of requested sizes, rather than retained or peak heap memory. */
    std::uint64_t bytes = 0;
};

/** Reset the counters and enable recording in this benchmark executable. */
void begin() noexcept;

/** @brief Disable recording and collect allocation activity.
 * @return Requested calls and bytes accumulated since the last begin().
 */
[[nodiscard]] AllocationSnapshot end() noexcept;
} // namespace Feel::benchmarking::allocation_probe
