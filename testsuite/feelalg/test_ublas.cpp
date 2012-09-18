/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-06-24

  Copyright (C) 2005,2006 EPFL

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-06-24
 */
#include <iostream>
#include <boost/timer.hpp>

#include <boost/numeric/ublas/operation.hpp>

#include <feel/feelalg/glas.hpp>
//#include <boost/numeric/bindings/traits/traits.hpp>
//#include <boost/numeric/bindings/traits/ublas_vector.hpp>
//#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
//#include <boost/numeric/bindings/atlas/cblas.hpp>
//#include <boost/numeric/bindings/blas/blas.hpp>

//namespace blas = boost::numeric::bindings::blas;
//namespace traits = boost::numeric::bindings::traits;


int main( int argc, char** argv )
{
    const int S = 10;
    int N;

    if ( argc == 2 )
        N = std::atoi( argv[1] );

    else
        N = 100000;

    using namespace Feel;
    typedef ublas::matrix<double, ublas::column_major> dcm_t;
    typedef ublas::matrix<double, ublas::row_major> drm_t;

    dcm_t mc1( 1, S );
    mc1 = ublas::scalar_matrix<double>( mc1.size1(), mc1.size2(), 1 );
    dcm_t mc2( 1, S );
    mc2 = ublas::scalar_matrix<double>( mc2.size1(), mc2.size2(), 2 );

    LOG(INFO) << "mc1 = " << mc1 << "\n";
    LOG(INFO) << "mc2 = " << mc2 << "\n";

    dcm_t mc3( S, S );

    drm_t r1 = ublas::prod( ublas::trans( mc1 ), mc2 );
    drm_t r2 = ublas::prod( mc1, ublas::trans( mc2 ) );
    LOG(INFO) << "mc3 = trans(mc1)*mc2 = " << r1 << "\n";
    LOG(INFO) << "scalar = mc1*trans(mc2) = " << r2 << "\n";

#if 0
    blas::gemm( traits::TRANSPOSE, traits::NO_TRANSPOSE,
                1.0, mc1, mc2,
                0.0, mc3 );
    LOG(INFO) << "mc3 = trans(mc1)*mc2 using gemm = " << mc3 << "\n";
#endif
    boost::timer timer;

    for ( int i = 0; i < N; ++i )
    {
        mc3.assign( ublas::scalar_matrix<double>( mc3.size1(), mc3.size2(), 1 ) );
        mc3.assign( ublas::prod( ublas::trans( mc1 ), mc2 ) );
    }

    LOG(INFO) << "ublas::prod : " << timer.elapsed() << "\n";

    timer.restart();
    dcm_t mc6( S, S );

    for ( int i = 0; i < N; ++i )
    {
        mc6.assign( ublas::scalar_matrix<double>( mc3.size1(), mc3.size2(), 1 ) );
        ublas::opb_prod( ublas::trans( mc1 ), mc2, mc6 );
    }

    LOG(INFO) << "ublas::obp_prod : " << timer.elapsed() << "\n";

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

    LOG(INFO) << "blas::gemm : " << timer.elapsed() << "\n";
#endif
    LOG(INFO) << "||mc3-mc4|| : " << ublas::norm_frobenius( mc3-mc4 ) << "\n";
    LOG(INFO) << "||mc3-mc6|| : " << ublas::norm_frobenius( mc3-mc6 ) << "\n";

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
    LOG(INFO) << "blas::gemm mc3*mc4 : " << timer.elapsed() << "\n";

    timer.restart();

    for ( int i = 0; i < N; ++i )
    {
        mc5.assign( ublas::scalar_matrix<double>( mc5.size1(), mc5.size2(), 1 ) );
        mc5.assign( prod( mc3, mc4 ) );
    }

    LOG(INFO) << "ublas::prod mc3*mc4 : " << timer.elapsed() << "\n";

    timer.restart();

    for ( int i = 0; i < N; ++i )
    {
        mc5.assign( ublas::scalar_matrix<double>( mc5.size1(), mc5.size2(), 1 ) );
        ublas::axpy_prod( mc3, mc4,mc5 );
    }

    LOG(INFO) << "ublas::axpy_prod mc3*mc4 : " << timer.elapsed() << "\n";

}
