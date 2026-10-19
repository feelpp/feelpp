/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-

   SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   SPDX-FileCopyrightText: 2026 University of Strasbourg
   SPDX-License-Identifier: LGPL-2.1-or-later
*/

#define BOOST_TEST_MODULE otherstatistics
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelstats/iterativecovariance.hpp>
#include <feel/feelstats/iterativeextrema.hpp>
#include <feel/feelstats/iterativemoments.hpp>

#include <concepts>
#include <string>

FEELPP_ENVIRONMENT_NO_OPTIONS

/** @brief Verify mergeable extrema for a nonnumeric value type. */
BOOST_AUTO_TEST_CASE( ordered_extrema )
{
    using Extrema = Feel::Stats::IterativeExtrema<std::string>;
    static_assert( std::same_as<Extrema::value_type, std::string> );
    Extrema left;
    Extrema right;
    left.add( "pear" );
    left.add( "apple" );
    right.add( "orange" );
    left.merge( right );
    BOOST_CHECK_EQUAL( left.count(), 3 );
    BOOST_CHECK_EQUAL( left.min(), "apple" );
    BOOST_CHECK_EQUAL( left.max(), "pear" );
}

/** @brief Verify centered covariance and Chan's merge on paired samples. */
BOOST_AUTO_TEST_CASE( paired_covariance )
{
    Feel::Stats::IterativeCovariance<double> left;
    Feel::Stats::IterativeCovariance<double> right;
    for ( double x : { 1.0, 2.0 } )
        left.add( x, 2 * x + 1 );
    for ( double x : { 3.0, 4.0 } )
        right.add( x, 2 * x + 1 );
    left.merge( right );
    BOOST_CHECK_EQUAL( left.count(), 4 );
    BOOST_CHECK_CLOSE( left.meanX(), 2.5, 1e-10 );
    BOOST_CHECK_CLOSE( left.meanY(), 6.0, 1e-10 );
    BOOST_CHECK_CLOSE( left.covariance(), 2.5, 1e-10 );
    BOOST_CHECK_CLOSE( left.sampleCovariance(), 10.0 / 3.0, 1e-10 );
    BOOST_CHECK_CLOSE( left.correlation(), 1.0, 1e-10 );
}

/** @brief Verify second through fourth moments before and after a merge. */
BOOST_AUTO_TEST_CASE( fourth_order_moments )
{
    Feel::Stats::IterativeMoments<double> left;
    Feel::Stats::IterativeMoments<double> right;
    for ( double value : { -2.0, -1.0, 0.0 } )
        left.add( value );
    for ( double value : { 1.0, 2.0 } )
        right.add( value );
    left.merge( right );
    BOOST_CHECK_EQUAL( left.count(), 5 );
    BOOST_CHECK_SMALL( left.mean(), 1e-12 );
    BOOST_CHECK_CLOSE( left.variance(), 2.0, 1e-10 );
    BOOST_CHECK_SMALL( left.skewness(), 1e-12 );
    BOOST_CHECK_CLOSE( left.kurtosis(), 1.7, 1e-10 );
    BOOST_CHECK_CLOSE( left.excessKurtosis(), -1.3, 1e-10 );
}
