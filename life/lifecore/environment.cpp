/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2010-04-14

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
   \file environment.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2010-04-14
 */
#include <life/lifecore/environment.hpp>

#if defined( HAVE_PETSC_H )
extern "C"
{
#include <petsc.h>
#include <petscerror.h>
}
#if defined( HAVE_SLEPC )
# include <slepc/slepc.h>
#endif /* HAVE_SLEPC */

#endif /* HAVE_PETSC_H */

namespace Life
{
Environment::Environment()
    :
    M_env(),
    M_comm()
{
#if defined ( HAVE_PETSC_H )
    PetscTruth is_petsc_initialized;
    PetscInitialized( &is_petsc_initialized );

    if ( !is_petsc_initialized )
    {
        i_initialized = true;
        int ierr = PetscInitializeNoArguments();

        boost::ignore_unused_variable_warning(ierr);
        CHKERRABORT(M_comm,ierr);
    }
#endif // HAVE_PETSC_H
}
Environment::Environment( int& argc, char**& argv )
    :
    M_env( argc, argv ),
    M_comm()
{

#if defined ( HAVE_PETSC_H )
    PetscTruth is_petsc_initialized;
    PetscInitialized( &is_petsc_initialized );

    if ( !is_petsc_initialized )
    {
        i_initialized = true;
#if defined( HAVE_SLEPC )
        int ierr = SlepcInitialize(&argc,&argv, PETSC_NULL, PETSC_NULL );
#else
        int ierr = PetscInitialize( &argc, &argv, PETSC_NULL, PETSC_NULL );
#endif

        boost::ignore_unused_variable_warning(ierr);
        CHKERRABORT(M_comm,ierr);
    }
#endif // HAVE_PETSC_H

}

Environment::~Environment()
{
    if ( i_initialized )
    {
#if defined ( HAVE_PETSC_H )
#if defined( HAVE_SLEPC )
        SlepcFinalize();
#else
        PetscFinalize();
#endif
#endif
    }
}

bool
Environment::initialized()
{
    return mpi::environment::initialized();
}

bool
Environment::finalized()
{
    return mpi::environment::finalized();
}

}
