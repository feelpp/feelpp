/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-06-12

  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file test_openmp.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-06-12
 */
#include <iostream>

#include <boost/timer.hpp>

#include <life/lifecore/life.hpp>

#if defined( HAVE_OPENMP )

#include <omp.h>

int main()
{
    const int size = 100000000;
    double* x = new double[size];
    double* y = new double[size];
    double start;
    double end;
    start = omp_get_wtime();


    //omp_set_dynamic(1);
    omp_set_num_threads(2);

#pragma omp parallel
    {
#pragma omp single
        printf("Num threads in dynamic region is = %d\n",
               omp_get_num_threads());

        int i;
#pragma omp for nowait
        for( i = 1; i < size; ++i)
            {
                x[i] = (y[i-1] + y[i+1])/2;
                x[i]++;
            }
    }
    end = omp_get_wtime();
    printf("Work took %f seconds\n", end - start);

    start = omp_get_wtime();


    //omp_set_dynamic();
    omp_set_num_threads(1);

#pragma omp parallel
    {
#pragma omp single
        printf("Num threads in dynamic region is = %d\n",
               omp_get_num_threads());

#pragma omp for
        for(int i = 1; i < size; ++i)
            {
                x[i] = (y[i-1] + y[i+1])/2;
                x[i]++;
            }
    }
    end = omp_get_wtime();
    printf("Work took %f seconds\n", end - start);

    printf("prcision: %f seconds\n", omp_get_wtick());

}
#else
int main()
{
    std::cout << "This compiler does not support OpenMP.\n"
              << "Exiting.\n";
    return 1;
}
#endif
