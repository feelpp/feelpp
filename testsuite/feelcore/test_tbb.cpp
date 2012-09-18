/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-10-26

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file test_tbb.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-10-26
 */
#define USE_BOOST_TEST 1
// Boost.Test

// make sure that the init_unit_test function is defined by UTF
#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE TBB testsuite
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN


#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;
#include <boost/test/floating_point_comparison.hpp>

#include <Eigen/Core>

#include <feel/feelcore/feel.hpp>
#include <feel/options.hpp>
#include <feel/feelcore/environment.hpp>

#if defined (FEELPP_HAS_TBB)
int N=100000000;
class A
{
public:
    A(  )  {}
    ~A() {}
    void operator() ( const tbb::blocked_range<std::vector<double>::iterator >& r ) const
    {
        for ( auto i = r.begin(); i != r.end(); ++i )
        {
            *i += std::cos( *i )+std::sin( *i );
        }
    }
};
class B
{
public:
    B(  ): m( Eigen::Matrix<double,3,3>::Zero() )
    {}
    B( B const&b  ): m( b.m )
    {}
    B( B& b, tbb::split )
        :
        m( Eigen::Matrix<double,3,3>::Zero() )
    {
    }

    ~B() {}
    void operator() ( const tbb::blocked_range<std::vector<double>::iterator >& r )
    {
        for ( auto i = r.begin(); i != r.end(); ++i )
        {
            m += Eigen::Matrix<double,3,3>::Constant( *i );
        }
    }
    void join( B const& other )
    {
        m += other.m;
    }
    Eigen::Matrix<double,3,3> m;
};

BOOST_AUTO_TEST_CASE( test_tbb )
{

    BOOST_TEST_MESSAGE( "Test TBB" );
    std::vector<double> vs( N,M_PI ),v( N,M_PI );
    tbb::blocked_range<std::vector<double>::iterator> r( v.begin(), v.end() );
    BOOST_CHECK_EQUAL( r.size(), N );
    tbb::tick_count serial_t0 = tbb::tick_count::now();

    for ( auto i = vs.begin(); i != vs.end(); ++i )
    {
        *i += std::cos( *i )+std::sin( *i );
    }

    tbb::tick_count serial_t1 = tbb::tick_count::now();

    int n = tbb::task_scheduler_init::default_num_threads();

    for ( int p=1; p<=n; ++p )
    {
        BOOST_TEST_MESSAGE( "[test_tbb] start tests with nthreads = " << p );
        tbb::task_scheduler_init init( p );

        tbb::tick_count parallel_t0 = tbb::tick_count::now();
        tbb::parallel_for( r, A() );
        tbb::tick_count parallel_t1 = tbb::tick_count::now();

        BOOST_TEST_MESSAGE( "Serial version ran in " << ( serial_t1 - serial_t0 ).seconds() << " seconds" << "\n"
                            << "Parallel version ran in " <<  ( parallel_t1 - parallel_t0 ).seconds() << " seconds" << "\n"
                            << "Resulting in a speedup of " << ( serial_t1 - serial_t0 ).seconds() / ( parallel_t1 - parallel_t0 ).seconds() << "\n" );
    }

    BOOST_TEST_MESSAGE( "Test TBB done" );

}
BOOST_AUTO_TEST_CASE( test_tbb_reduce )
{

    BOOST_TEST_MESSAGE( "Test TBB Reduce" );
    std::vector<double> vs( N,M_PI ),v( N,M_PI );
    tbb::blocked_range<std::vector<double>::iterator> r( v.begin(), v.end() );
    BOOST_CHECK_EQUAL( r.size(), N );
    Eigen::Matrix<double,3,3> m;
    tbb::tick_count serial_t0 = tbb::tick_count::now();

    for ( auto i = vs.begin(); i != vs.end(); ++i )
    {
        m += Eigen::Matrix<double,3,3>::Constant( *i );
    }

    tbb::tick_count serial_t1 = tbb::tick_count::now();
    std::cout << m << "\n";

    int n = tbb::task_scheduler_init::default_num_threads();

    for ( int p=1; p<=n; ++p )
    {
        BOOST_TEST_MESSAGE( "[test_tbb_reduce] start tests with nthreads = " << p );
        tbb::task_scheduler_init init( p );

        B b;

        tbb::tick_count parallel_t0 = tbb::tick_count::now();
        tbb::parallel_reduce( r, b );
        tbb::tick_count parallel_t1 = tbb::tick_count::now();

        std::cout << b.m << "\n";
        BOOST_TEST_MESSAGE( "Serial version ran in " << ( serial_t1 - serial_t0 ).seconds() << " seconds" << "\n"
                            << "Parallel version ran in " <<  ( parallel_t1 - parallel_t0 ).seconds() << " seconds" << "\n"
                            << "Resulting in a speedup of " << ( serial_t1 - serial_t0 ).seconds() / ( parallel_t1 - parallel_t0 ).seconds() << "\n" );
    }

    BOOST_TEST_MESSAGE( "Test TBB reduce done" );

}

#endif // FEELPP_HAS_TBB
