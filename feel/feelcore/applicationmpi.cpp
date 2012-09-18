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
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>

#if defined(FEELPP_HAS_TAU)
#include <Profile/Profiler.h>
#endif /* FEELPP_HAS_TAU */


namespace Feel
{
FEELPP_NO_EXPORT
po::options_description
mpiOptions()
{
    po::options_description mpi( "MPI options" );
    mpi.add_options()
    ( "disable-mpi", "Disable MPI" )
    ;
    return mpi;
}

bool Application::_S_is_mpi_initialized = false;

#if defined( FEELPP_HAS_MPI )
MPI_Comm Application::COMM_WORLD = MPI_COMM_NULL;

Application::Application( int argc,
                          char** argv,
                          AboutData const& ad,
                          MPI_Comm comm )
#else
Application::Application( int argc,
                          char** argv,
                          AboutData const& ad )
#endif // FEELPP_HAS_MPI
    :
    super( argc, argv, ad, mpiOptions(), false )
#if defined( FEELPP_HAS_MPI )
    ,
    M_env()
#endif
{
#if defined( FEELPP_HAS_MPI )
    int is_mpi_initialized;
    MPI_Initialized ( &is_mpi_initialized );

    //std::cout << "is_mpi_initialized = " << is_mpi_initialized << "\n";
    if ( !is_mpi_initialized )
    {
        M_env = boost::shared_ptr<mpi::environment>( new mpi::environment( argc, argv ) );
        ///std::cout << "processor name = " << M_env->processor_name() << "\n";
#if 0
        //int __argc = this->unknownArgc();
        //char** __argv = this->unknownArgv();
        MPI_Init ( &argc, &argv );
#endif
        _S_is_mpi_initialized = true;
    }

    MPI_Comm_dup ( comm, &Application::COMM_WORLD );
    //MPI_Comm_dup ( comm, (MPI_Comm*)&S_world );
#if 0
    MPI_Comm_rank ( Application::COMM_WORLD, &_S_process_id );
    MPI_Comm_size ( Application::COMM_WORLD, &_S_n_process );
#else
    _S_process_id = S_world.rank();
    _S_n_process = S_world.size();
#endif
#endif // FEELPP_HAS_MPI

    this->doOptions( argc, argv );

#if defined( FEELPP_HAS_MPI )
    char * __env = getenv( "DEBUG" );
    std::string env_str;

    if ( __env )
        env_str = __env;

    mpi::broadcast( S_world, env_str, 0 );

    if ( _S_process_id != 0 )
    {
        setenv( "DEBUG", env_str.c_str(), 1 );
        //VLOG(1) << "DEBUG is set to " << env_str << "\n";
        //std::cout << "DEBUG is set to " << env_str << "\n";
    }

#endif // MPI

#if defined(FEELPP_HAS_TAU)
    TAU_PROFILE_SET_NODE( _S_process_id );
#endif /* FEELPP_HAS_TAU */
}


#if defined( FEELPP_HAS_MPI )
Application::Application( int argc,
                          char** argv,
                          AboutData const& ad,
                          po::options_description const& od,
                          MPI_Comm comm )
#else
Application::Application( int argc,
                          char** argv,
                          AboutData const& ad,
                          po::options_description const& od )
#endif // FEELPP_HAS_MPI
    :
    super( argc, argv, ad, mpiOptions().add( od ), false )
#if defined( FEELPP_HAS_MPI )
    ,
    M_env()
#endif
{
#if defined( FEELPP_HAS_MPI )
    int is_mpi_initialized;
    MPI_Initialized ( &is_mpi_initialized );

    //std::cout << "is_mpi_initialized = " << is_mpi_initialized << "\n";
    if ( !is_mpi_initialized )
    {
        M_env = boost::shared_ptr<mpi::environment>( new mpi::environment( argc, argv ) );
        //std::cout << "processor name = " << M_env->processor_name() << "\n";
#if 0
        //int __argc = this->unknownArgc();
        //char** __argv = this->unknownArgv();
        MPI_Init ( &argc, &argv );
#endif
        _S_is_mpi_initialized = true;
    }

    //Application::COMM_WORLD = comm;
    MPI_Comm_dup ( comm, &Application::COMM_WORLD );
    //MPI_Comm_dup ( comm, (MPI_Comm*)&S_world );
#if 0
    MPI_Comm_rank ( Application::COMM_WORLD, &_S_process_id );
    MPI_Comm_size ( Application::COMM_WORLD, &_S_n_process );
#else
    _S_process_id = S_world.rank();
    _S_n_process = S_world.size();
#endif
#endif

    this->doOptions( argc, argv );

#if defined( FEELPP_HAS_MPI )
    char * __env = getenv( "DEBUG" );
    std::string env_str;

    if ( __env )
        env_str = __env;

    //Application::Broadcast( env_str, 0 );
    mpi::broadcast( S_world, env_str, 0 );

    if ( _S_process_id != 0 )
    {
        setenv( "DEBUG", env_str.c_str(), 1 );
        //VLOG(1) << "DEBUG is set to " << env_str << "\n";
        //std::cout << "DEBUG is set to " << env_str << "\n";
    }

#endif // MPI

#if defined(FEELPP_HAS_TAU)
    TAU_PROFILE_SET_NODE( _S_process_id );
#endif /* FEELPP_HAS_TAU */
}

Application::~Application()
{
#if 0 // defined ( FEELPP_HAS_MPI )
    MPI_Comm_free ( &Application::COMM_WORLD );

    if ( _S_is_mpi_initialized )
        MPI_Finalize();

#endif // MPI
}

mpi::communicator Application::S_world;

}
