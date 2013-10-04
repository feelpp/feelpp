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
   \file MatrixTest.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2004-08-29
*/
#ifndef __MatrixTest_H
#define __MatrixTest_H 1

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

#include <feel/feelalg/glas.hpp>


namespace Feel
{
template<typename ordering>
struct matrix
{
    typedef ublas::compressed_matrix<double,ordering> type;
};

template<typename ordering>
class MatrixMass
{
public:
    typedef double value_type;
    typedef typename matrix<ordering>::type matrix_type;

    MatrixMass( int n )
        :
        M_mat( n, n )
    {
        // Defining constants.

        value_type sub  = 1.0/value_type( n+1 );
        value_type diag = 4.0/value_type( n+1 );

        for ( int __i = 0; __i < n; ++__i )
        {
            if ( __i != 0 )
                M_mat( __i, __i-1 ) = sub;

            M_mat( __i, __i ) = diag;

            if ( __i != n-1 )
                M_mat( __i, __i+1 ) = sub;
        }
    }

    ~MatrixMass()
    {
    }
    uint const * ia()const
    {
        return M_mat.index1_data();
    }
    uint const * jaData()const
    {
        return M_mat.index2_data();
    }
    double* valueData()
    {
        return M_mat.value_data();
    }

    matrix_type const& matrix() const
    {
        return M_mat;
    }
    matrix_type &      matrix()
    {
        return M_mat;
    }

private:
    matrix_type M_mat;
};
#if 0
/*!
      \class MatrixConvectionDiffusion


    | T -I          |
    |-I  T -I       |
A = |   -I  T       |
    |        ...  -I|
    |           -I T|

    derived from the standard central difference discretization of the
     2-dimensional convection-diffusion operator (Laplacian u) + rho*(du/dx)
    on a unit square with zero Dirichlet boundary conditions.
    When rho*h/2 <= 1, the discrete convection-diffusion operator has real
    eigenvalues.  When rho*h/2 > 1, it has COMPLEX eigenvalues.

    */
class MatrixConvectionDiffusion
{
public:
    typedef double value_type;
    typedef CSRMatr<CSRPatt,value_type> matrix_type;

    MatrixConvectionDiffusion( int nx, value_type __rho = 0.0 );

    ~MatrixConvectionDiffusion()
    {
        delete M_mat;
        delete M_pattern;
    }
    uint const * iaData()const
    {
        return M_mat->Patt()->giveRawCSR_ia();
    }
    uint const * jaData()const
    {
        return M_mat->Patt()->giveRawCSR_ja();
    }
    double* valueData()
    {
        return M_mat->giveRawCSR_value();
    }

    matrix_type const& matrix() const
    {
        return *M_mat;
    }
    matrix_type &      matrix()
    {
        return *M_mat;
    }

private:
    value_type M_rho;
    matrix_type* M_mat;
    CSRPatt* M_pattern;
    std::vector<double> M_val;
};
#endif
}
#endif /* __MatrixTest_H */
