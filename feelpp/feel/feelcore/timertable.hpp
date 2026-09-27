/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-

   SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   SPDX-FileContributor: Guillaume Dollé <dolle.guillaume@gmail.com>

   SPDX-FileCopyrightText: 2012-2026 University of Strasbourg
   SPDX-License-Identifier: LGPL-2.1-or-later
*/
#ifndef FEELPP_TIMERTABLE_HPP
#define FEELPP_TIMERTABLE_HPP 1

#include <cstddef>
#include <map>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include <feel/feelcore/journalwatcher.hpp>

namespace Feel
{

/**
 * @brief All local elapsed-time samples recorded for one tic()/toc() label.
 *
 * TimerTable owns these records. Samples are stored in insertion order.
 */
class TimerData : public std::vector<double>
{
public:
    /** @brief Create an empty timer record. */
    TimerData() = default;
    /** @brief Copy the timer label, journal name, and recorded samples. */
    TimerData( TimerData const& ) = default;

    /** @brief Create a labeled timer record with an optional journal instance name. */
    explicit TimerData( std::string const& message, std::string const& instanceName = {} );

    /** @brief Append one elapsed duration and its nesting level. */
    void add( std::pair<double, int> const& sample );

    /// Human-readable timer label.
    std::string msg;
    /// Name of this timer record in the environment journal.
    std::string instance_name;
    /// Nesting level of the most recent sample.
    int level = 0;
};

/** @brief Local timer samples and their text, Markdown, and journal reports. */
class TimerTable : private std::map<std::string, TimerData>, public JournalWatcher
{
public:
    /** @brief Register a local timer table with the environment journal. */
    TimerTable();
    /** @brief Finalize journal registration before destroying samples. */
    ~TimerTable() override;

    /** @brief Store a local sample. Empty labels are ignored. */
    void add( std::string const& message, std::pair<double, int> const& sample,
              std::string const& instanceName = {} );
    /** @brief Return the latest local duration or NaN when the label is absent. */
    double lastOrNaN( std::string const& message ) const;
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
