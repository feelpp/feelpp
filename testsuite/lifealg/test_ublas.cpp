/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-06-24

  Copyright (C) 2005,2006 EPFL

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
   \file test_ublas.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-06-24
 */
#include <iostream>
#include <boost/timer.hpp>

#include <boost/log/log.hpp>
#include <boost/log/functions.hpp>
#include <boost/numeric/ublas/operation.hpp>
BOOST_DECLARE_LOG(info);
BOOST_DEFINE_LOG(info, "info");

void write_to_cout(const std::string &, const std::string &msg) { std::cout << msg; }

void init_logs() {
    using namespace boost::logging;
    // all logs prefix the message by time
    //add_modifier("*", prepend_time("$hh:$mm:$ss "), DEFAULT_INDEX + 1 );
    manipulate_logs("*").add_modifier(prepend_time("$hh:$mm:$ss "), DEFAULT_INDEX + 1 );
    // all log' messages are prefixed by the log name ("app", or "DEBUG" or "info")
    manipulate_logs("*").add_modifier(prepend_prefix);

    //manipulate_logs("*").add_appender(&write_to_cout);
    manipulate_logs("*").add_appender(write_to_file( "test_ublas.log") );
    flush_log_cache();
}

#include <life/lifealg/glas.hpp>
//#include <boost/numeric/bindings/traits/traits.hpp>
//#include <boost/numeric/bindings/traits/ublas_vector.hpp>
//#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
//#include <boost/numeric/bindings/atlas/cblas.hpp>
//#include <boost/numeric/bindings/blas/blas.hpp>

//namespace blas = boost::numeric::bindings::blas;
//namespace traits = boost::numeric::bindings::traits;


int main( int argc, char** argv )
{
    init_logs();

    const int S = 10;
    int N;
    if ( argc == 2 )
        N = std::atoi( argv[1] );
    else
        N = 100000;

    using namespace Life;
    typedef ublas::matrix<double, ublas::column_major> dcm_t;
    typedef ublas::matrix<double, ublas::row_major> drm_t;

    dcm_t mc1( 1, S );
    mc1 = ublas::scalar_matrix<double>( mc1.size1(), mc1.size2(), 1 );
    dcm_t mc2( 1, S );
    mc2 = ublas::scalar_matrix<double>( mc2.size1(), mc2.size2(), 2 );

    BOOST_LOG(info) << "mc1 = " << mc1 << std::endl;
    BOOST_LOG(info) << "mc2 = " << mc2 << std::endl;

    dcm_t mc3( S, S );

    BOOST_LOG(info) << "mc3 = trans(mc1)*mc2 = " << ublas::prod( ublas::trans( mc1 ), mc2 ) << std::endl;
    BOOST_LOG(info) << "scalar = mc1*trans(mc2) = " << ublas::prod( mc1, ublas::trans( mc2 ) ) << std::endl;

#if 0
    blas::gemm( traits::TRANSPOSE, traits::NO_TRANSPOSE,
                1.0, mc1, mc2,
                0.0, mc3 );
    BOOST_LOG(info) << "mc3 = trans(mc1)*mc2 using gemm = " << mc3 << std::endl;
#endif
    boost::timer timer;

    for ( int i = 0; i < N; ++i )
    {
        mc3.assign( ublas::scalar_matrix<double>( mc3.size1(), mc3.size2(), 1 ) );
        mc3.assign( ublas::prod( ublas::trans( mc1 ), mc2 ) );
    }
    BOOST_LOG(info) << "ublas::prod : " << timer.elapsed() << std::endl;

    timer.restart();
    dcm_t mc6( S, S );
    for ( int i = 0; i < N; ++i )
    {
        mc6.assign( ublas::scalar_matrix<double>( mc3.size1(), mc3.size2(), 1 ) );
        ublas::opb_prod( ublas::trans( mc1 ), mc2, mc6 );
    }
    BOOST_LOG(info) << "ublas::obp_prod : " << timer.elapsed() << std::endl;

    timer.restart();
    dcm_t mc4( S, S );
#if 0
    for ( int i = 0; i < N; ++i )
    {
        mc4.assign( ublas::scalar_matrix<double>( mc3.size1(), mc3.size2(), 1 ) );
        blas::gemm( traits::TRANSPOSE, traits::NO_TRANSPOSE,
                    1.0, mc1, mc2,
                    0.0, mc4 );
    }
    BOOST_LOG(info) << "blas::gemm : " << timer.elapsed() << std::endl;
#endif
    BOOST_LOG(info) << "||mc3-mc4|| : " << ublas::norm_frobenius( mc3-mc4 ) << std::endl;
    BOOST_LOG(info) << "||mc3-mc6|| : " << ublas::norm_frobenius( mc3-mc6 ) << std::endl;

    timer.restart();

    dcm_t mc5( S, S );
#if 0
    for ( int i = 0; i < N; ++i )
    {
        mc5.assign( ublas::scalar_matrix<double>( mc5.size1(), mc5.size2(), 1 ) );
        blas::gemm( traits::NO_TRANSPOSE, traits::NO_TRANSPOSE,
                    1.0, mc3, mc4,
                    0.0, mc5 );
    }
#endif
    BOOST_LOG(info) << "blas::gemm mc3*mc4 : " << timer.elapsed() << std::endl;

    timer.restart();

    for ( int i = 0; i < N; ++i )
    {
        mc5.assign( ublas::scalar_matrix<double>( mc5.size1(), mc5.size2(), 1 ) );
        mc5.assign( prod( mc3, mc4 ) );
    }
    BOOST_LOG(info) << "ublas::prod mc3*mc4 : " << timer.elapsed() << std::endl;

    timer.restart();

    for ( int i = 0; i < N; ++i )
    {
        mc5.assign( ublas::scalar_matrix<double>( mc5.size1(), mc5.size2(), 1 ) );
        ublas::axpy_prod( mc3, mc4,mc5 );
    }
    BOOST_LOG(info) << "ublas::axpy_prod mc3*mc4 : " << timer.elapsed() << std::endl;

}
