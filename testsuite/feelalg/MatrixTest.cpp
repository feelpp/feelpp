/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2004-08-29

  Copyright (C) 2004 EPFL

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
   \file MatrixTest.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2004-08-29
*/
#include <MatrixTest.hpp>

namespace Feel
{
//
// Mass
//
MatrixMass::MatrixMass( int n )
    :
    M_mat( 0 ),
    M_pattern( 0 ),
    M_val()
{
    // Defining constants.

    value_type sub  = 1.0/value_type( n+1 );
    value_type diag = 4.0/value_type( n+1 );

    // Defining the number of nonzero matrix elements.

    int nnz = 3*n-2;

    std::vector<uint> __ia( n+1 ), __ja( nnz );
    M_val.resize( nnz );

    __ia[0] = 0;
    int __j = 0;

    for ( int __i = 0; __i < n; ++__i )
    {
        if ( __i != 0 )
        {
            M_val[__j] =  sub;
            __ja[__j++] = __i-1;
        }

        M_val[__j] =  diag;
        __ja[__j++] = __i;

        if ( __i != ( n-1 ) )
        {
            M_val[__j] =  sub;
            __ja[__j++] = __i+1;
        }

        __ia[__i+1] = __j;
    }

    M_pattern = new CSRPatt( nnz, n, n, __ia, __ja );
    M_mat = new CSRMatr<CSRPatt, double>( *M_pattern, M_val );

}


//
// Convection Diffusion
//

MatrixConvectionDiffusion::MatrixConvectionDiffusion( int nx, value_type __rho )
    :
    M_rho( __rho ),
    M_mat( 0 ),
    M_pattern( 0 ),
    M_val()
{

    int N = nx*nx;
    int nnz = 5*N-4*nx;

    std::vector<uint> __ia( N+1 ), __ja( nnz );
    M_val.resize( nnz );

    __ia[0] = 0;
    int i = 0;

    value_type h  = 1.0/value_type( nx+1 );
    value_type h2 = h*h;
    value_type dd = 4.0/h2;
    value_type df = -1.0/h2;
    value_type dl = df - 0.5*M_rho/h;
    value_type du = df + 0.5*M_rho/h;

    for ( int j = 0; j < N; j++ )
    {
        if ( j >= nx )
        {
            M_val[i] = df;
            __ja[i++]=j-nx;
        }

        if ( ( j%nx ) != 0 )
        {
            M_val[i] = du;
            __ja[i++] = j-1;
        }

        M_val[i] = dd;
        __ja[i++] = j;

        if ( ( ( j+1 )%nx ) != 0 )
        {
            M_val[i] = dl;
            __ja[i++] = j+1;
        }

        if ( j < N-nx )
        {
            M_val[i] = df;
            __ja[i++] = j+nx;
        }

        __ia[j+1]=i;
    }

    M_pattern = new CSRPatt( nnz, N, N, __ia, __ja );
    M_mat = new CSRMatr<CSRPatt, double>( *M_pattern, M_val );
    assert( M_mat != 0 );

    std::cerr << __PRETTY_FUNCTION__ << " matrix constructed" << std::endl;
}
}
