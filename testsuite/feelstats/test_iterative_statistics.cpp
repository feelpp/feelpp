/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-

   SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   SPDX-FileCopyrightText: 2026 University of Strasbourg
   SPDX-License-Identifier: LGPL-2.1-or-later
*/

#define BOOST_TEST_MODULE iterative_statistics
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelstats/iterativestatistics.hpp>

#include <cmath>

FEELPP_ENVIRONMENT_NO_OPTIONS

/** @brief Verify online moments and Chan's merge against exact values. */
BOOST_AUTO_TEST_CASE( online_and_merged_moments )
{
    Feel::Stats::IterativeStatistics<double> left;
    Feel::Stats::IterativeStatistics<double> right;
    for ( int value : { 1, 2, 3, 4 } )
        left.add( value );
    for ( int value : { 5, 6, 7, 8 } )
        right.add( value );
    left.merge( right );

    BOOST_CHECK_EQUAL( left.count(), 8 );
    BOOST_CHECK_EQUAL( left.total(), 36 );
    BOOST_CHECK_EQUAL( left.min(), 1 );
    BOOST_CHECK_EQUAL( left.max(), 8 );
    BOOST_CHECK_CLOSE( left.mean(), 4.5, 1e-10 );
    BOOST_CHECK_CLOSE( left.variance(), 5.25, 1e-10 );
    BOOST_CHECK_CLOSE( left.stdDev(), std::sqrt( 5.25 ), 1e-10 );
}

/** @brief Verify compensated summation when tiny durations follow a large one. */
BOOST_AUTO_TEST_CASE( mixed_scale_total_and_empty_merge )
{
    Feel::Stats::IterativeStatistics<double> statistics;
    Feel::Stats::IterativeStatistics<double> empty;
    BOOST_CHECK( statistics.empty() );
    statistics.merge( empty );
    statistics.add( 1e9 );
    for ( int i = 0; i < 1000000; ++i )
        statistics.add( 1e-9 );
    BOOST_CHECK_SMALL( statistics.total() - ( 1e9 + 1e-3 ), 1e-6 );
    BOOST_CHECK_EQUAL( statistics.count(), 1000001 );
    BOOST_CHECK_EQUAL( statistics.min(), 1e-9 );
    BOOST_CHECK_EQUAL( statistics.max(), 1e9 );
}
