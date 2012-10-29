/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-11-13

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
   \file opuseadspfem.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-11-13
 */
#include <octave/oct.h>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/options.hpp>


#include <opusmodelrb.hpp>
#include <feel/feelcrb/pfemapp.hpp>

static boost::shared_ptr<Feel::PFemApp<Feel::OpusModelRB<2,1,2> > > app;

DEFUN_DLD ( opuseadspfem, args, nargout, "Opus EADS Test Case (pfem)" )
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
                app = boost::shared_ptr<Feel::PFemApp<Feel::OpusModelRB<2,1,2> > >( new Feel::PFemApp<Feel::OpusModelRB<2,1,2> >( Feel::makeEadsAbout( "eadspfem" ),
                        Feel::OpusData::makeOptions() ) );
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

