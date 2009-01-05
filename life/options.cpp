/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-07-21

  Copyright (C) 2007-2008 Université Joseph Fourier (Grenoble I)

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
   \file options.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-07-21
 */
#include <life/options.hpp>
#include <life/lifealg/backend.hpp>
#include <life/lifealg/backendgmm.hpp>
#include <life/lifealg/backendpetsc.hpp>
#include <life/lifealg/solvereigenslepc.hpp>
#include <life/lifealg/backendtrilinos.hpp>
#include <life/lifediscr/oseendata.hpp>
#include <life/lifediscr/bdf2.hpp>
#include <life/lifefilters/exporter.hpp>
#include <life/lifematerial/materiallib.hpp>

namespace Life
{

po::options_description
life_options( std::string const& prefix  )
{
    return
        /* alg options */
        backend_options()
        .add( backendgmm_options( prefix ) )
#if defined(HAVE_PETSC_H)
        .add( backendpetsc_options( prefix ) )
#endif
#if defined(HAVE_SLEPC) && defined(HAVE_PETSC)
        .add( solvereigenslepc_options( prefix ) )
#endif // HAVE_SLEPC && HAVE_PETSC
#if defined( HAVE_TRILINOS_EPETRA )
        .add( backendtrilinos_options( prefix ) )
#endif
        /* nonlinear solver options */
        .add( nlsolver_options() )

        /* discr options */
        .add( oseen_options( prefix ) )
        .add( bdf_options( prefix ) )

        /* exporter options */
        .add( exporter_options( prefix ) )

        /* material options */
        .add( material_options( prefix ) );


}
}
