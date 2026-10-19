/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-

   SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   SPDX-FileCopyrightText: 2026 University of Strasbourg
   SPDX-License-Identifier: LGPL-2.1-or-later
*/

#define BOOST_TEST_MODULE kll
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelstats/kll.hpp>

#include <cmath>
#include <concepts>

FEELPP_ENVIRONMENT_NO_OPTIONS

/** @brief Verify exact rank selection while every observation is retained. */
BOOST_AUTO_TEST_CASE( exact_before_compaction )
{
    static_assert( std::same_as<Feel::Stats::KllSketch<double>::value_type, double> );
    Feel::Stats::KllSketch<double> sketch;
    for ( double value : { 4.0, 1.0, 3.0, 2.0 } )
        sketch.add( value );
    BOOST_CHECK_EQUAL( sketch.count(), 4 );
    BOOST_CHECK_EQUAL( sketch.retainedCount(), 4 );
    BOOST_CHECK_EQUAL( sketch.quantile( 0 ), 1 );
    BOOST_CHECK_EQUAL( sketch.quantile( 0.5 ), 2 );
    BOOST_CHECK_EQUAL( sketch.quantile( 0.99 ), 4 );
    BOOST_CHECK_EQUAL( sketch.quantile( 1 ), 4 );
}

/** @brief Verify bounded retention and approximate ranks after merging. */
BOOST_AUTO_TEST_CASE( bounded_sketch_and_merge )
{
    Feel::Stats::KllSketch<double> left( 200 );
    Feel::Stats::KllSketch<double> right( 200 );
    for ( int value = 0; value < 100000; ++value )
    {
        if ( value % 2 == 0 )
            left.add( value );
        else
            right.add( value );
    }
    left.merge( right );

    BOOST_CHECK_EQUAL( left.count(), 100000 );
    BOOST_CHECK_EQUAL( left.min(), 0 );
    BOOST_CHECK_EQUAL( left.max(), 99999 );
    BOOST_CHECK_LT( left.retainedCount(), 4000 );
    BOOST_CHECK_SMALL( left.quantile( 0.5 ) - 50000, 3000.0 );
    BOOST_CHECK_SMALL( left.quantile( 0.9 ) - 90000, 3000.0 );
    BOOST_CHECK_SMALL( left.quantile( 0.99 ) - 99000, 3000.0 );
}

/** @brief Verify quantiles of a distribution with widely separated timings. */
BOOST_AUTO_TEST_CASE( mixed_scale_quantiles )
{
    Feel::Stats::KllSketch<double> sketch;
    for ( int i = 0; i < 10000; ++i )
        sketch.add( i % 10 == 0 ? 1e6 : 1e-9 );
    BOOST_CHECK_EQUAL( sketch.quantile( 0.5 ), 1e-9 );
    BOOST_CHECK_EQUAL( sketch.quantile( 0.99 ), 1e6 );
}
