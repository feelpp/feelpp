/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-06-20

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file octave_wrapper.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-06-20
 */
#include <octave/oct.h>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/options.hpp>

#include <@CRB_MODEL_SHORT_NAME@.hpp>
#include <feel/feelcrb/opusapp.hpp>

static boost::shared_ptr<Feel::OpusApp<Feel::@CRB_MODEL_LONG_NAME@ > > app;

DEFUN_DLD ( @CRB_MODEL_WRAPPER_NAME@, args, nargout, "Oct file for @CRB_MODEL_LONG_NAME@" )
{
    int nargin = args.length ();

    if ( nargin != 1 )
        print_usage ();

    else
    {
        NDArray A = args( 0 ).array_value ();
        dim_vector dims = A.dims();
        dim_vector ydims( 2 );
        ydims( 0 )=dims( 0 );
        ydims( 1 )=2;

        NDArray Y( ydims );
        {
            static bool is_init = false;

            if ( !is_init )
            {
                app = boost::shared_ptr<Feel::OpusApp<Feel::@CRB_MODEL_LONG_NAME@ > >( new Feel::OpusApp<Feel::@CRB_MODEL_LONG_NAME@ >( Feel::make@CRB_MODEL_LONG_NAME@About( "@CRB_MODEL_SHORT_NAME@" ),
                        Feel::make@CRB_MODEL_LONG_NAME@Options() ) );
                app->setMode( @CRB_MODEL_WRAPPER_TYPE@ );
                is_init = true;
            }


            for ( int i = 0; i < dims( 0 ); ++i )
            {
                std::vector<double> x( dims( 1 ) );

                for ( int j = 0; j < dims( 1 ); ++j )
                    x[j] = A( i,j );

                std::vector<double> y( 2 );
                app->run( x.data(), dims( 1 ), y.data(), 2 );
                Y( i,0 )=y[0];
                Y( i,1 )=y[1];
            }
        }

        //MPI_Finalize();
        if ( ! error_state )
            return octave_value ( Y );
    }

    return octave_value_list ();
}


