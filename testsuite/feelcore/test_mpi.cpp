/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-10-18

  Copyright (C) 2005,2006 EPFL

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
   \file test_mpi.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-10-18
 */
#include <feel/feelcore/application.hpp>
#include <feel/feelcore/application.hpp>

Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "test_mpi" ,
                           "test_mpi" ,
                           "0.1",
                           "MPI test",
                           Feel::AboutData::License_LGPL,
                           "Copyright (c) 2005,2006 EPFL" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "" );
    return about;

}
int main( int argc, char** argv )
{
    using namespace Feel;

    Feel::po::options_description test( "test options" );
    test.add_options()
    ( "testall", "run all test cases" )
    ;
#if defined( FEELPP_HAS_PETSC_H )
    //MPI_Init( &argc, &argv );
    //PetscInitialize (&argc, &argv, PETSC_NULL, PETSC_NULL);
    Application app( argc, argv, makeAbout(), test );

    std::cout << "[Application] N process: " << Application::nProcess() << "\n"
              << "[Application] Id : " << Application::processId() << "\n";
#else
    Application app( argc, argv, makeAbout(), test );
    std::cout << "[Application] N process: " << Application::nProcess() << "\n"
              << "[Application] Id : " << Application::processId() << "\n";
#endif
    //PetscFinalize();
}
