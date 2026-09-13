// SPDX-License-Identifier: LGPL-3.0-or-later
//! \file exporterio.hpp
//! \brief Minimal exporter policy for safe shared-file writes.
#ifndef FEELPP_FILTERS_EXPORTERIO_HPP
#define FEELPP_FILTERS_EXPORTERIO_HPP 1

namespace Feel
{
//! \brief Choose the owner of distributed payload writes.
//! This is a correctness policy, not an MPI tuning interface. Automatic selects
//! one bounded-memory writer on detected Linux NFS; other filesystems retain
//! independent writes. Independent requires a validated filesystem/backend.
enum class ExporterIOPolicy
{
    Automatic,   //!< Detect Linux NFS; use independent writes otherwise.
    Independent, //!< Disjoint per-rank writes on a validated filesystem.
    Root         //!< One bounded-memory streaming writer per file.
};
} // namespace Feel

#endif // FEELPP_FILTERS_EXPORTERIO_HPP
