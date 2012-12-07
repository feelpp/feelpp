/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-10-18

  Copyright (C) 2005,2006,2009 EPFL
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
   \file application.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-10-18
 */
#include <boost/concept_check.hpp>

#include <feel/feelcore/application.hpp>

#if defined(FEELPP_HAS_SLEPC)
#include <slepc/slepceps.h>
#endif /* FEELPP_HAS_SLEPC */

namespace Feel
{
#if defined( FEELPP_HAS_PETSC_H )

FEELPP_NO_EXPORT
po::options_description
petscOptions()
{
    po::options_description petsc( "PETSC options" );
    petsc.add_options()
    ( "disable-petsc", "disable petsc" )
    ;
    return petsc;
}

Application::Application( int argc,
                          char** argv,
                          AboutData const& ad,
                          MPI_Comm comm )
    :
    super( argc, argv, ad, petscOptions(), comm )
{
    PETSC_COMM_WORLD = super::COMM_WORLD;
    int __argc = this->unknownArgc();
    char** __argv = this->unknownArgv();
#if defined( FEELPP_HAS_SLEPC )
    int ierr = SlepcInitialize( &__argc,&__argv, PETSC_NULL, PETSC_NULL );
#else
    int ierr = PetscInitialize( &__argc, &__argv, PETSC_NULL, PETSC_NULL );
#endif
    std::cerr << "[Application] argc " << __argc << "\n";

    for ( int i = 0; i < argc; ++i )
        std::cerr << "[Application] argv[" << i << "]="<< __argv[i] << "\n";

    boost::ignore_unused_variable_warning( ierr );
    CHKERRABORT( super::COMM_WORLD,ierr );
}


Application::Application( int argc,
                          char** argv,
                          AboutData const& ad,
                          po::options_description const& od,
                          MPI_Comm comm )
    :
    super( argc, argv, ad, petscOptions().add( od ), comm )
{
    PETSC_COMM_WORLD = super::COMM_WORLD;
    int __argc = this->unknownArgc();
    char** __argv = this->unknownArgv();
#if defined( FEELPP_HAS_SLEPC )
    int ierr = SlepcInitialize( &__argc,&__argv, PETSC_NULL, PETSC_NULL );
#else
    int ierr = PetscInitialize( &__argc, &__argv, PETSC_NULL, PETSC_NULL );
#endif
#if 0
    std::cerr << "[Application] argc " << __argc << "\n";
    --__argc;

    for ( int i = 0; i < argc; ++i )
        if ( __argv[i] )
            std::cerr << "[Application] argv[" << i << "]="<< __argv[i] << "\n";

#endif
    //int ierr = PetscInitializeNoArguments();
    boost::ignore_unused_variable_warning( ierr );
    CHKERRABORT( super::COMM_WORLD,ierr );
}

Application::~Application()
{
#if defined( FEELPP_HAS_SLEPC )
    SlepcFinalize();
#else
    PetscFinalize();
#endif
}

#endif // PETSC

}
