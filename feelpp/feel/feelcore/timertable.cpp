/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-

   SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

   SPDX-FileCopyrightText: 2012-2026 University of Strasbourg
   SPDX-License-Identifier: LGPL-2.1-or-later
*/

#include <feel/feelcore/timertable.hpp>

#include <feel/feelcore/environment.hpp>

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>

namespace Feel
{
namespace
{
/** @brief Derived statistics for one local timer record. */
struct TimerSummary
{
    /// Number of recorded samples.
    std::size_t count = 0;
    /// Sum of recorded durations in seconds.
    double total = 0;
    /// Minimum recorded duration in seconds.
    double min = 0;
    /// Maximum recorded duration in seconds.
    double max = 0;
    /// Arithmetic mean duration in seconds.
    double mean = 0;
    /// Population standard deviation in seconds.
    double stdDev = 0;
};

/** @brief Compute descriptive statistics for a timer's local samples. */
TimerSummary summarize( TimerData const& timer )
{
    TimerSummary result;
    result.count = timer.size();
    if ( result.count == 0 )
        return result;

    auto [minimum, maximum] = std::minmax_element( timer.begin(), timer.end() );
    result.min = *minimum;
    result.max = *maximum;
    result.total = std::accumulate( timer.begin(), timer.end(), 0.0 );
    result.mean = result.total / result.count;
    double variance = 0;
    for ( double value : timer )
        variance += ( value - result.mean ) * ( value - result.mean );
    result.stdDev = std::sqrt( variance / result.count );
    return result;
}

/** @brief Allocate a process-local unique journal timer name. */
std::string nextTimerInstanceName()
{
    static std::atomic<std::uint64_t> nextId{0};
    return "timer-" + std::to_string( nextId.fetch_add( 1, std::memory_order_relaxed ) );
}
} // namespace

TimerData::TimerData( std::string const& message, std::string const& instanceName )
    : msg( message ), instance_name( instanceName.empty() ? nextTimerInstanceName() : instanceName )
{}

void TimerData::add( std::pair<double, int> const& sample )
{
    this->push_back( sample.first );
    level = sample.second;
}

TimerTable::TimerTable()
    : JournalWatcher( "TimerTable", false )
{}

TimerTable::~TimerTable()
{
    this->journalFinalize();
}

void TimerTable::add( std::string const& message, std::pair<double, int> const& sample,
                      std::string const& instanceName )
{
    if ( message.empty() )
        return;
    auto it = this->find( message );
    if ( it == this->end() )
        it = this->emplace( message, TimerData( message, instanceName ) ).first;
    it->second.add( sample );
    M_max_len = std::max( M_max_len, message.size() + 2 * static_cast<std::size_t>( std::max( sample.second, 0 ) ) );
}

double TimerTable::lastOrNaN( std::string const& message ) const
{
    auto it = this->find( message );
    if ( it == this->end() || it->second.empty() )
        return std::numeric_limits<double>::quiet_NaN();
    return it->second.back();
}

void TimerTable::save( bool display )
{
    if ( !display || !Environment::isMasterRank() )
        return;

    std::ostringstream output;
    auto printHeader = [&]()
    {
        output << std::setw( M_max_len ) << std::left << "Timer" << ' '
               << std::setw( 7 ) << std::right << "Count" << ' '
               << std::setw( 11 ) << "Total(s)" << ' '
               << std::setw( 11 ) << "Max(s)" << ' '
               << std::setw( 11 ) << "Min(s)" << ' '
               << std::setw( 11 ) << "Mean(s)" << ' '
               << std::setw( 11 ) << "StdDev(s)" << '\n';
    };
    auto printRow = [&]( TimerData const& timer )
    {
        auto summary = summarize( timer );
        auto indent = 2 * static_cast<std::size_t>( std::max( timer.level, 0 ) );
        output << std::setw( indent ) << " "
               << std::setw( M_max_len > indent ? M_max_len - indent : 0 ) << std::left << timer.msg << ' '
               << std::setw( 7 ) << std::right << summary.count << ' '
               << std::setw( 11 ) << std::scientific << std::setprecision( 2 ) << summary.total << ' '
               << std::setw( 11 ) << summary.max << ' '
               << std::setw( 11 ) << summary.min << ' '
               << std::setw( 11 ) << summary.mean << ' '
               << std::setw( 11 ) << summary.stdDev << '\n';
    };

    std::vector<std::pair<double, TimerData const*>> byTotal;
    byTotal.reserve( this->size() );
    printHeader();
    for ( auto const& [label, timer] : *this )
    {
        printRow( timer );
        byTotal.emplace_back( summarize( timer ).total, &timer );
    }
    std::stable_sort( byTotal.begin(), byTotal.end(), []( auto const& left, auto const& right )
                      { return left.first > right.first; } );
    output << "--------------------------------------------------------------------------------\n";
    printHeader();
    for ( auto const& [total, timer] : byTotal )
        printRow( *timer );
    std::cout << output.str() << std::endl;
}

void TimerTable::saveMD( std::ostream& output )
{
    auto printHeader = [&]()
    {
        output << "| Timer | Count | Total(s) | Max(s) | Min(s) | Mean(s) | StdDev(s) |\n"
               << "| --- | ---: | ---: | ---: | ---: | ---: | ---: |\n";
    };
    auto printRow = [&]( TimerData const& timer )
    {
        auto summary = summarize( timer );
        output << "| " << timer.msg << " | " << summary.count << " | "
               << std::scientific << std::setprecision( 2 ) << summary.total << " | "
               << summary.max << " | " << summary.min << " | " << summary.mean << " | "
               << summary.stdDev << " |\n";
    };

    std::vector<std::pair<double, TimerData const*>> byTotal;
    byTotal.reserve( this->size() );
    printHeader();
    for ( auto const& [label, timer] : *this )
    {
        printRow( timer );
        byTotal.emplace_back( summarize( timer ).total, &timer );
    }
    std::stable_sort( byTotal.begin(), byTotal.end(), []( auto const& left, auto const& right )
                      { return left.first > right.first; } );
    output << '\n';
    printHeader();
    for ( auto const& [total, timer] : byTotal )
        printRow( *timer );
    output << '\n' << std::defaultfloat;
}

void TimerTable::test()
{}

void TimerTable::updateInformationObject( nl::json& output ) const
{
    output.clear();
    for ( auto const& [label, timer] : *this )
    {
        auto summary = summarize( timer );
        auto& value = output[timer.instance_name];
        value["message"] = timer.msg;
        value["count"] = summary.count;
        value["total"] = summary.total;
        value["max"] = summary.max;
        value["min"] = summary.min;
        value["mean"] = summary.mean;
        value["std_dev"] = summary.stdDev;
    }
}

} // namespace Feel
