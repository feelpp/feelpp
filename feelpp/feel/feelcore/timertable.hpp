/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-

   SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   SPDX-FileContributor: Guillaume Dollé <dolle.guillaume@gmail.com>

   SPDX-FileCopyrightText: 2012-2026 University of Strasbourg
   SPDX-License-Identifier: LGPL-2.1-or-later
*/
#ifndef FEELPP_TIMERTABLE_HPP
#define FEELPP_TIMERTABLE_HPP 1

#include <cstddef>
#include <cstdint>
#include <map>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include <feel/feelcore/journalwatcher.hpp>
#include <feel/feelstats/iterativestatistics.hpp>
#include <feel/feelstats/kll.hpp>

namespace Feel
{

/** @brief Serializable local moments and quantile sketch for one label. */
struct TimerProfileSnapshot
{
    /// Mergeable local count, sum, extrema, mean, and variance state.
    Stats::IterativeStatistics<double> statistics;
    /// Mergeable quantile summary of local durations.
    Stats::KllSketch<double> quantiles{};

    /** @brief Serialize the profile for the MPI gather. */
    template <typename Archive>
    void serialize( Archive& archive, unsigned int )
    {
        archive & statistics & quantiles;
    }
};

/**
 * @brief Streaming local statistics for one tic()/toc() label.
 *
 * TimerTable owns these records. Summary statistics are updated on each sample.
 * Moments use constant storage; the equal-capacity quantile sketch retains
 * O(k log N) values for N calls and capacity k per level. Raw samples are not
 * retained. Recording is intended for a single thread per process.
 */
class TimerData
{
public:
    /** @brief Create an empty timer record. */
    TimerData() = default;
    /** @brief Copy the timer label, journal name, and accumulated statistics. */
    TimerData( TimerData const& ) = default;

    /** @brief Create a labeled timer record with an optional journal instance name. */
    explicit TimerData( std::string const& message, std::string const& instanceName = {} );

    /** @brief Accumulate a finite, nonnegative duration; ignore invalid samples. */
    void add( std::pair<double, int> const& sample );

    /** @brief Number of accumulated samples. */
    std::size_t size() const { return M_Statistics.count(); }
    /** @brief Whether any sample has been accumulated. */
    bool empty() const { return M_Statistics.empty(); }
    /** @brief Latest elapsed duration; valid when the record is nonempty. */
    double back() const { return M_Last; }
    /** @brief Sum of elapsed durations in seconds. */
    double total() const { return M_Statistics.total(); }
    /** @brief Smallest elapsed duration in seconds. */
    double min() const { return M_Statistics.min(); }
    /** @brief Largest elapsed duration in seconds. */
    double max() const { return M_Statistics.max(); }
    /** @brief Arithmetic mean elapsed duration in seconds. */
    double mean() const { return M_Statistics.mean(); }
    /** @brief Population standard deviation in seconds. */
    double stdDev() const;
    /** @brief Sum of squared deviations, for merging rank statistics. */
    double m2() const { return M_Statistics.m2(); }
    /** @brief Approximate quantile of all samples for this label. */
    double quantile( double probability ) const { return M_Quantiles.quantile( probability ); }
    /** @brief Copy the local statistics needed for collective profiling. */
    TimerProfileSnapshot profileSnapshot() const;

    /// Human-readable timer label.
    std::string msg;
    /// Name of this timer record in the environment journal.
    std::string instance_name;
    /// Nesting level of the most recent sample.
    int level = 0;

private:
    /// Mergeable moments and extrema of durations for this label.
    Stats::IterativeStatistics<double> M_Statistics;
    /// Latest accumulated duration in seconds.
    double M_Last = 0;
    /// Quantile summary with bounded retention per level.
    Stats::KllSketch<double> M_Quantiles;
};

/**
 * @brief Local timer statistics and their text, Markdown, and journal reports.
 *
 * Each tic()/toc() call updates this table locally. Gather profiles only after
 * the timed work; every rank in the communicator must call the gather with the
 * same labels and rank selection. Rank IDs are local to that communicator.
 * The pooled profile combines every call from the selected ranks for each label.
 *
 * @par Record repeated calls under one or more labels
 * @code{.cpp}
 * #include <feel/feeltiming/tic.hpp>
 *
 * auto const& comm = Feel::Environment::worldComm().localComm();
 * Feel::tic();
 * // Assemble a system.
 * Feel::toc( "assembly", Feel::no_display );
 * Feel::tic();
 * // Solve the system.
 * Feel::toc( "solve", Feel::no_display );
 * @endcode
 *
 * @par Default: all labels recorded on communicator rank 0
 * @code{.cpp}
 * auto const& master = Feel::Environment::gatherTimerProfileStatistics( comm );
 * // One report per label on rank 0; each report.ranks contains only rank 0.
 * @endcode
 *
 * @par Specific labels, with rank 0 selected by default
 * @code{.cpp}
 * auto const& chosenLabels = Feel::Environment::gatherTimerProfileStatistics(
 *     { "assembly", "solve" }, comm );
 * // Reports follow the requested label order, even if a label has no calls.
 * @endcode
 *
 * @par All labels from all ranks
 * @code{.cpp}
 * auto const& all = Feel::Environment::gatherTimerProfileStatistics( comm, {} );
 * // Labels are the sorted union across ranks; each report has all rank profiles.
 * @endcode
 *
 * @par All labels from one rank or a group of ranks
 * @code{.cpp}
 * int targetRank = comm.size() - 1;
 * auto const& one = Feel::Environment::gatherTimerProfileStatistics( comm, { targetRank } );
 * auto const& group = Feel::Environment::gatherTimerProfileStatistics( comm, { 0, targetRank } );
 * // group[0].pooled covers only rank 0 and targetRank (once if they coincide).
 * @endcode
 *
 * @par Specific labels from all ranks or a group
 * @code{.cpp}
 * auto const& allAssembly = Feel::Environment::gatherTimerProfileStatistics(
 *     { "assembly" }, comm, {} );
 * auto const& groupSolve = Feel::Environment::gatherTimerProfileStatistics(
 *     { "solve" }, comm, { 0, targetRank } );
 * @endcode
 *
 * Each gather replaces the cached reports, so consume a returned reference
 * before the next gather. Quantiles are approximate; counts and moments are
 * accumulated from all recorded calls.
 */
class TimerTable : private std::map<std::string, TimerData>, public JournalWatcher
{
public:
    /** @brief Register a local timer table with the environment journal. */
    TimerTable();
    /** @brief Finalize journal registration before destroying samples. */
    ~TimerTable() override;

    /** @brief Store a valid duration. Empty labels and invalid durations are ignored. */
    void add( std::string const& message, std::pair<double, int> const& sample,
              std::string const& instanceName = {} );
    /** @brief Return all recorded local labels in lexical order. */
    std::vector<std::string> labels() const;
    /** @brief Return a local profile, or an empty profile for an absent label. */
    TimerProfileSnapshot profileSnapshot( std::string const& message ) const;
    /** @brief Print local timer statistics on the master rank when display is true. */
    void save( bool display );
    /** @brief Write local timer statistics as Markdown. */
    void saveMD( std::ostream& output );
    /** @brief Retained compatibility hook; currently performs no work. */
    void test();

private:
    /** @brief Refresh the timer section of the Environment journal. */
    void updateInformationObject( nl::json& output ) const override;

    /// Width of the longest indented timer label in text output.
    std::size_t M_max_len = 0;
};

} // namespace Feel

#endif
