/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2010-08-05

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
   \file residualestimator.oct
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2010-08-05
 */
#include <octave/oct.h>

#include <residualestimator.hpp>

#define OCTNAME(name,dim, order) BOOST_PP_CAT(BOOST_PP_CAT(BOOST_PP_CAT(name,dim),_),order)

#define RESIDUAL_ESTIMATOR( dim, order )                                \
    DEFUN_DLD (OCTNAME(residualestimator_,dim,order), args, nargout, "Residual Estimator for the Laplacian") \
    {                                                                   \
        int nargin = args.length ();                                    \
        if (nargin != 1)                                                \
            print_usage ();                                             \
        else                                                            \
        {                                                               \
            NDArray A = args(0).array_value ();                         \
            dim_vector dims = A.dims();                                 \
            dim_vector ydims(2);                                        \
            ydims(0)=dims(0);                                           \
            ydims(1)=4;                                                 \
                                                                        \
            NDArray Y(ydims);                                           \
            {                                                           \
                static bool is_init = false;                            \
                if ( !is_init )                                         \
                {                                                       \
                    if (!Life::Environment::initialized() )             \
                        new Life::Environment();                        \
                    is_init = true;                                     \
                }                                                       \
                boost::shared_ptr<ResidualEstimator<dim,order> > OCTNAME(app_,dim,order)( new ResidualEstimator<dim,order>( makeAbout() ) ); \
                for( int i = 0; i < dims(0); ++i )                      \
                {                                                       \
                    std::vector<double> x( dims(1) );                   \
                    for( int j = 0; j < dims(1); ++j )                  \
                    {                                                   \
                        x[j] = A(i,j);                                  \
                        std::cout << "x["<< j << "]=" << x[j] << "\n";  \
                    }                                                   \
                    std::vector<double> y( 4 );                         \
                    OCTNAME(app_,dim,order)->run( x.data(), dims(1), y.data(), 4 ); \
                    Y(i,0)=y[0];                                        \
                    Y(i,1)=y[1];                                        \
                    Y(i,2)=y[2];                                        \
                    Y(i,3)=y[3];                                        \
                }                                                       \
            }                                                           \
            if (! error_state)                                          \
                return octave_value (Y);                                \
        }                                                               \
        return octave_value_list ();                                    \
    }
/**/

