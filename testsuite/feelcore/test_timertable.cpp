/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-

   SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   SPDX-FileCopyrightText: 2026 University of Strasbourg
   SPDX-License-Identifier: LGPL-2.1-or-later
*/

#define BOOST_TEST_MODULE timertable
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/timertable.hpp>
#include <feel/feeltiming/tic.hpp>

#include <cmath>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

FEELPP_ENVIRONMENT_NO_OPTIONS

/** @brief Invalid durations neither alter accumulated statistics nor create labels. */
BOOST_AUTO_TEST_CASE( invalid_timer_samples )
{
    using namespace Feel;
    TimerData timer( "valid-durations" );
    TimerTable table;
    for ( double value : { -1.0, std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::infinity(),
                           -std::numeric_limits<double>::infinity() } )
    {
        timer.add( { value, -1 } );
        table.add( "invalid-only", { value, -1 } );
    }
    BOOST_CHECK( timer.empty() );
    BOOST_CHECK( table.labels().empty() );
    BOOST_CHECK_EQUAL( timer.profileSnapshot().quantiles.count(), 0 );

    timer.add( { 0.0, 1 } );
    timer.add( { 2.0, 1 } );
    timer.add( { -1.0, -1 } );
    BOOST_CHECK_EQUAL( timer.size(), 2 );
    BOOST_CHECK_EQUAL( timer.total(), 2.0 );
    BOOST_CHECK_EQUAL( timer.mean(), 1.0 );
    BOOST_CHECK_EQUAL( timer.back(), 2.0 );
    BOOST_CHECK_EQUAL( timer.level, 1 );
    BOOST_CHECK_EQUAL( timer.profileSnapshot().quantiles.count(), 2 );

    table.add( "valid-zero", { 0.0, 0 } );
    BOOST_REQUIRE_EQUAL( table.labels().size(), 1 );
    BOOST_CHECK_EQUAL( table.profileSnapshot( "valid-zero" ).statistics.count(), 1 );
}

/** @brief An unmatched toc() keeps its error result out of timer reports. */
BOOST_AUTO_TEST_CASE( unmatched_toc_is_not_recorded )
{
    using namespace Feel;
    BOOST_REQUIRE( details::sec_timer.times().empty() );
    BOOST_CHECK_EQUAL( toc( "timertable-unmatched-toc", no_display ), -1.0 );
    BOOST_CHECK( details::sec_timer.times().empty() );

    auto const& comm = Environment::worldComm().localComm();
    auto const& reports = Environment::gatherTimerProfileStatistics( comm, {} );
    BOOST_CHECK( std::none_of( reports.begin(), reports.end(), []( auto const& report )
                               { return report.label == "timertable-unmatched-toc"; } ) );
}

/** @brief Verify local Welford statistics over many repeated durations. */
BOOST_AUTO_TEST_CASE( online_statistics )
{
    Feel::TimerData timer( "repeated" );
    BOOST_CHECK( timer.empty() );
    BOOST_CHECK_EQUAL( timer.stdDev(), 0.0 );

    for ( int i = 0; i < 100000; ++i )
        timer.add( { 1.0 + i % 4, 2 } );

    BOOST_CHECK_EQUAL( timer.size(), 100000 );
    BOOST_CHECK_EQUAL( timer.back(), 4.0 );
    BOOST_CHECK_EQUAL( timer.total(), 250000.0 );
    BOOST_CHECK_EQUAL( timer.min(), 1.0 );
    BOOST_CHECK_EQUAL( timer.max(), 4.0 );
    BOOST_CHECK_CLOSE( timer.mean(), 2.5, 1e-10 );
    BOOST_CHECK_CLOSE( timer.stdDev(), std::sqrt( 1.25 ), 1e-10 );
    BOOST_CHECK_EQUAL( timer.level, 2 );
}

/** @brief Verify one-line phase laps and disabled profiling leave tic()/toc() balanced. */
BOOST_AUTO_TEST_CASE( conditional_profile_laps )
{
    using namespace Feel;
    /** @brief Phase totals updated by the conditional lap timer. */
    struct PhaseTotals
    {
        /// Accumulated metadata phase durations.
        double metadata = 0;
        /// Accumulated marker phase durations.
        double marker = 0;
        /// Accumulated final phase durations.
        double container = 0;
    } totals;

    auto const& comm = Environment::worldComm().localComm();
    auto stackSize = details::sec_timer.times().size();
    {
        ProfileLapTimer timer{ &totals };
        BOOST_CHECK_EQUAL( details::sec_timer.times().size(), stackSize );
        timer.lap( "timertable-lap-metadata", &PhaseTotals::metadata );
        timer.lap( "timertable-lap-marker", &PhaseTotals::marker );
        timer.lap( "timertable-lap-metadata", &PhaseTotals::metadata );
        timer.finish( "timertable-lap-container", &PhaseTotals::container );
    }
    BOOST_CHECK_EQUAL( details::sec_timer.times().size(), stackSize );
    BOOST_CHECK( totals.metadata >= 0 && totals.marker >= 0 && totals.container >= 0 );

    auto recordedTotals = totals;
    ProfileLapTimer disabled{ static_cast<PhaseTotals*>( nullptr ) };
    disabled.lap( "timertable-lap-disabled", &PhaseTotals::metadata );
    disabled.finish( "timertable-lap-disabled", &PhaseTotals::container );
    BOOST_CHECK_EQUAL( details::sec_timer.times().size(), stackSize );
    BOOST_CHECK_EQUAL( totals.metadata, recordedTotals.metadata );
    BOOST_CHECK_EQUAL( totals.container, recordedTotals.container );

    auto const& reports = Environment::gatherTimerProfileStatistics(
        { "timertable-lap-metadata", "timertable-lap-marker",
          "timertable-lap-container", "timertable-lap-disabled" }, comm );
    BOOST_REQUIRE_EQUAL( reports.size(), 4 );
    BOOST_CHECK_EQUAL( reports[0].pooled.count, 2 );
    BOOST_CHECK_EQUAL( reports[1].pooled.count, 1 );
    BOOST_CHECK_EQUAL( reports[2].pooled.count, 1 );
    BOOST_CHECK_EQUAL( reports[3].pooled.count, 0 );
}

/** @brief Gather every call from each rank in one pooled profile. */
BOOST_AUTO_TEST_CASE( profile_includes_all_calls )
{
    using namespace Feel;
    auto const& comm = Environment::worldComm().localComm();
    double rank = static_cast<double>( comm.rank() );
    Environment::addTimer( "timertable-all-calls", { rank + 1.0, 0 }, "" );
    Environment::addTimer( "timertable-all-calls", { rank + 2.0, 0 }, "" );

    auto const& reports = Environment::gatherTimerProfileStatistics( { "timertable-all-calls" }, comm, {} );
    BOOST_REQUIRE_EQUAL( reports.size(), 1 );
    BOOST_REQUIRE_EQUAL( reports.front().ranks.size(), comm.size() );
    BOOST_CHECK_EQUAL( reports.front().pooled.count, 2 * comm.size() );
    for ( int sampleRank = 0; sampleRank < comm.size(); ++sampleRank )
    {
        BOOST_CHECK_EQUAL( reports.front().ranks.at( sampleRank ).count, 2 );
        BOOST_CHECK_EQUAL( reports.front().ranks.at( sampleRank ).total, 2 * sampleRank + 3 );
    }
}

/** @brief Verify balanced timer calls and the profile's sample counts. */
BOOST_AUTO_TEST_CASE( balanced_tic_toc_calls )
{
    using namespace Feel;
    bool allFinite = true;
    for ( int i = 0; i < 10000; ++i )
    {
        tic();
        double elapsed = toc( "timertable-balanced-calls", no_display );
        allFinite = allFinite && std::isfinite( elapsed ) && elapsed >= 0;
    }
    BOOST_CHECK( allFinite );
    BOOST_CHECK( details::sec_timer.times().empty() );

    auto const& comm = Environment::worldComm().localComm();
    auto const& profiles = Environment::gatherTimerProfileStatistics( { "timertable-balanced-calls" }, comm, {} );
    BOOST_REQUIRE_EQUAL( profiles.size(), 1 );
    BOOST_REQUIRE_EQUAL( profiles.front().ranks.size(), comm.size() );
    BOOST_CHECK_EQUAL( profiles.front().ranks.at( comm.rank() ).count, 10000 );
    BOOST_CHECK_EQUAL( profiles.front().pooled.count, 10000 * comm.size() );
}

/** @brief Verify per-rank profiles and count-weighted pooled statistics. */
BOOST_AUTO_TEST_CASE( per_rank_and_pooled_profiles )
{
    using namespace Feel;
    auto const& comm = Environment::worldComm().localComm();
    double base = 10.0 * comm.rank();
    for ( double value : { 1.0, 3.0, 5.0 } )
        Environment::addTimer( "timertable-profiles", { base + value, 0 }, "" );
    if ( comm.rank() == 0 )
        Environment::addTimer( "timertable-profiles", { 7.0, 0 }, "" );

    auto const& reports = Environment::gatherTimerProfileStatistics( { "timertable-profiles", "timertable-absent" }, comm, {} );
    BOOST_REQUIRE_EQUAL( reports.size(), 2 );
    BOOST_REQUIRE_EQUAL( reports[0].ranks.size(), comm.size() );
    BOOST_CHECK_EQUAL( reports[0].ranks.at( comm.rank() ).count, comm.rank() == 0 ? 4 : 3 );
    BOOST_CHECK_EQUAL( reports[0].ranks.at( comm.rank() ).minimum, base + 1 );
    BOOST_CHECK_EQUAL( reports[0].ranks.at( comm.rank() ).maximum, comm.rank() == 0 ? 7 : base + 5 );

    std::vector<double> allSamples;
    for ( int rank = 0; rank < comm.size(); ++rank )
    {
        double rankBase = 10.0 * rank;
        allSamples.push_back( rankBase + 1 );
        allSamples.push_back( rankBase + 3 );
        allSamples.push_back( rankBase + 5 );
    }
    allSamples.push_back( 7 );
    std::sort( allSamples.begin(), allSamples.end() );
    double total = 0;
    for ( double value : allSamples )
        total += value;
    double mean = total / allSamples.size();
    double variance = 0;
    for ( double value : allSamples )
        variance += ( value - mean ) * ( value - mean );
    variance /= allSamples.size();

    auto const& pooled = reports[0].pooled;
    BOOST_CHECK_EQUAL( pooled.count, allSamples.size() );
    BOOST_CHECK_CLOSE( pooled.total, total, 1e-10 );
    BOOST_CHECK_CLOSE( pooled.mean, mean, 1e-10 );
    BOOST_CHECK_CLOSE( pooled.stdDev, std::sqrt( variance ), 1e-10 );
    BOOST_CHECK_EQUAL( pooled.minimum, allSamples.front() );
    BOOST_CHECK_EQUAL( pooled.maximum, allSamples.back() );
    BOOST_CHECK_EQUAL( pooled.median, allSamples[static_cast<std::size_t>( std::ceil( 0.5 * allSamples.size() ) - 1 )] );
    BOOST_CHECK_EQUAL( pooled.p90, allSamples[static_cast<std::size_t>( std::ceil( 0.9 * allSamples.size() ) - 1 )] );

    BOOST_CHECK_EQUAL( reports[1].pooled.count, 0 );
    for ( auto const& [rank, statistics] : reports[1].ranks )
        BOOST_CHECK_EQUAL( statistics.count, 0 );

    auto json = Environment::timerProfileReportsJson();
    BOOST_CHECK_EQUAL( json.at( 0 ).at( "pooled" ).at( "count" ).get<std::size_t>(), allSamples.size() );
    BOOST_CHECK_EQUAL( json.at( 0 ).at( "ranks" ).size(), comm.size() );
}

/** @brief Select the master, one rank, and a deduplicated rank group. */
BOOST_AUTO_TEST_CASE( selected_rank_profiles )
{
    using namespace Feel;
    auto const& comm = Environment::worldComm().localComm();
    auto rank = comm.rank();
    auto lastRank = comm.size() - 1;
    Environment::addTimer( "timertable-selected", { rank + 1.0, 0 }, "" );
    Environment::addTimer( "timertable-selected", { rank + 3.0, 0 }, "" );

    auto const& master = Environment::gatherTimerProfileStatistics( { "timertable-selected" }, comm );
    BOOST_REQUIRE_EQUAL( master.front().ranks.size(), 1 );
    BOOST_CHECK_EQUAL( master.front().ranks.at( 0 ).count, 2 );
    BOOST_CHECK_EQUAL( master.front().pooled.total, 4.0 );
    BOOST_CHECK_EQUAL( Environment::timerProfileReportsJson().at( 0 ).at( "ranks" ).size(), 1 );

    auto const& specific = Environment::gatherTimerProfileStatistics( { "timertable-selected" }, comm, { lastRank } );
    BOOST_REQUIRE_EQUAL( specific.front().ranks.size(), 1 );
    BOOST_CHECK_EQUAL( specific.front().ranks.at( lastRank ).total, 2 * lastRank + 4 );
    BOOST_CHECK_EQUAL( specific.front().pooled.count, 2 );

    auto const& group = Environment::gatherTimerProfileStatistics( { "timertable-selected" }, comm,
                                                                   { lastRank, 0, lastRank } );
    BOOST_REQUIRE_EQUAL( group.front().ranks.size(), comm.size() == 1 ? 1 : 2 );
    BOOST_CHECK_EQUAL( group.front().pooled.count, comm.size() == 1 ? 2 : 4 );
    BOOST_CHECK_EQUAL( group.front().pooled.total, comm.size() == 1 ? 4 : 2 * lastRank + 8 );
    auto json = Environment::timerProfileReportsJson();
    BOOST_CHECK_EQUAL( json.at( 0 ).at( "ranks" ).at( 0 ).at( "rank" ).get<int>(), 0 );
    BOOST_CHECK_THROW( Environment::gatherTimerProfileStatistics( { "timertable-selected" }, comm,
                                                                 { comm.size() } ), std::out_of_range );
}

/** @brief Discover the union of labels on the selected communicator ranks. */
BOOST_AUTO_TEST_CASE( discovered_label_profiles )
{
    using namespace Feel;
    auto const& comm = Environment::worldComm().localComm();
    int lastRank = comm.size() - 1;
    auto rankLabel = []( int rank ) { return "timertable-discovered-rank-" + std::to_string( rank ); };
    Environment::addTimer( "timertable-discovered-shared", { comm.rank() + 1.0, 0 }, "" );
    Environment::addTimer( rankLabel( comm.rank() ), { comm.rank() + 2.0, 0 }, "" );

    auto findLabel = []( auto const& reports, std::string const& label )
    {
        return std::find_if( reports.begin(), reports.end(), [&]( auto const& report ) { return report.label == label; } );
    };

    auto const& master = Environment::gatherTimerProfileStatistics( comm );
    auto masterLabel = findLabel( master, rankLabel( 0 ) );
    BOOST_REQUIRE( masterLabel != master.end() );
    BOOST_CHECK_EQUAL( masterLabel->ranks.size(), 1 );
    BOOST_CHECK_EQUAL( masterLabel->pooled.count, 1 );
    BOOST_CHECK( findLabel( master, "timertable-discovered-shared" ) != master.end() );
    if ( comm.size() > 1 )
        BOOST_CHECK( findLabel( master, rankLabel( lastRank ) ) == master.end() );

    auto const& all = Environment::gatherTimerProfileStatistics( comm, {} );
    auto shared = findLabel( all, "timertable-discovered-shared" );
    BOOST_REQUIRE( shared != all.end() );
    BOOST_CHECK_EQUAL( shared->ranks.size(), comm.size() );
    BOOST_CHECK_EQUAL( shared->pooled.count, comm.size() );
    for ( int rank = 0; rank < comm.size(); ++rank )
    {
        auto found = findLabel( all, rankLabel( rank ) );
        BOOST_REQUIRE( found != all.end() );
        BOOST_CHECK_EQUAL( found->pooled.count, 1 );
        BOOST_CHECK_EQUAL( found->ranks.at( rank ).count, 1 );
    }

    auto const& group = Environment::gatherTimerProfileStatistics( comm, { 0, lastRank } );
    BOOST_CHECK( findLabel( group, rankLabel( 0 ) ) != group.end() );
    BOOST_CHECK( findLabel( group, rankLabel( lastRank ) ) != group.end() );
    if ( comm.size() > 2 )
        BOOST_CHECK( findLabel( group, rankLabel( 1 ) ) == group.end() );

    auto const& chosen = Environment::gatherTimerProfileStatistics(
        { "timertable-discovered-shared" }, comm, { lastRank } );
    BOOST_REQUIRE_EQUAL( chosen.size(), 1 );
    BOOST_CHECK_EQUAL( chosen.front().pooled.count, 1 );
    BOOST_CHECK_EQUAL( chosen.front().ranks.at( lastRank ).total, lastRank + 1.0 );
    BOOST_CHECK_THROW( Environment::gatherTimerProfileStatistics( comm, { comm.size() } ), std::out_of_range );
}
