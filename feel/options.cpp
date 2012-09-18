/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-07-21

  Copyright (C) 2007-2012 Universite Joseph Fourier (Grenoble I)

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
   \file options.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-07-21
 */
#include <feel/options.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/backendpetsc.hpp>
#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelalg/backendtrilinos.hpp>
#include <feel/feelmaterial/materiallib.hpp>

#include <feel/feeldiscr/bdf2.hpp>


#include <feel/feelfilters/exporter.hpp>

namespace Feel
{

po::options_description
feel_options( std::string const& prefix  )
{
    return
        /* alg options */
        backend_options()
#if defined(FEELPP_HAS_PETSC_H)
        .add( backendpetsc_options( prefix ) )
#endif
        .add( solvereigen_options( prefix ) )
#if defined( FEELPP_HAS_TRILINOS_EPETRA )
        .add( backendtrilinos_options( prefix ) )
#endif
        /* nonlinear solver options */
        .add( nlsolver_options() )

        /* discr options */
        .add( bdf_options( prefix ) )

        /* exporter options */
        .add( exporter_options( prefix ) )

        /* material options */
        .add( material_options( prefix ) );


}
}

