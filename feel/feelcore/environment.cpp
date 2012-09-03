/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2010-04-14

  Copyright (C) 2010,2011 Université Joseph Fourier (Grenoble I)

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
#include <feel/feelconfig.h>

#include <boost/program_options.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

#include <feel/feelcore/environment.hpp>

#include <feel/feelcore/feelpetsc.hpp>


namespace google
{
namespace glog_internal_namespace_
{
bool IsGoogleLoggingInitialized();
}
}
namespace Feel
{
Environment::Environment()
    :
    M_env(false)
{
    // Initialize Google's logging library.
    if ( !google::glog_internal_namespace_::IsGoogleLoggingInitialized() )
        google::InitGoogleLogging("feel++");
    google::InstallFailureSignalHandler();
#if defined( FEELPP_HAS_TBB )
    int n = tbb::task_scheduler_init::default_num_threads();
    Log() << "[Feel++] TBB running with " << n << " threads\n";
#else
    int n = 1 ;
#endif

    //tbb::task_scheduler_init init(1);

    mpi::communicator world;
#if defined ( FEELPP_HAS_PETSC_H )
    PetscTruth is_petsc_initialized;
    PetscInitialized( &is_petsc_initialized );

    if ( !is_petsc_initialized )
    {
        i_initialized = true;
        int ierr = PetscInitializeNoArguments();

        boost::ignore_unused_variable_warning( ierr );
        CHKERRABORT( world,ierr );
    }

    // make sure that petsc do not catch signals and hence do not print long
    //and often unuseful messages
    PetscPopSignalHandler();
#endif // FEELPP_HAS_PETSC_H

    S_worldcomm = worldcomm_type::New( world );
    FEELPP_ASSERT( S_worldcomm ).error( "creating worldcomm failed" );
}
Environment::Environment( int& argc, char**& argv )
    :
    M_env( argc, argv, false )
{
    // Initialize Google's logging library.
    if ( !google::glog_internal_namespace_::IsGoogleLoggingInitialized() )
    {
        if ( argc > 0 )
            google::InitGoogleLogging(argv[0]);
        else
            google::InitGoogleLogging("feel++");
    }
    google::InstallFailureSignalHandler();

#if defined( FEELPP_HAS_TBB )
    int n = tbb::task_scheduler_init::default_num_threads();
    //int n = 2;
    //Log() << "[Feel++] TBB running with " << n << " threads\n";
    //tbb::task_scheduler_init init(2);
#endif

    mpi::communicator world;
#if defined ( FEELPP_HAS_PETSC_H )
    PetscTruth is_petsc_initialized;
    PetscInitialized( &is_petsc_initialized );

    if ( !is_petsc_initialized )
    {
        i_initialized = true;
#if defined( FEELPP_HAS_SLEPC )
        int ierr = SlepcInitialize( &argc,&argv, PETSC_NULL, PETSC_NULL );
#else
        int ierr = PetscInitialize( &argc, &argv, PETSC_NULL, PETSC_NULL );
#endif
        boost::ignore_unused_variable_warning( ierr );
        CHKERRABORT( world,ierr );
    }

    // make sure that petsc do not catch signals and hence do not print long
    //and often unuseful messages
    PetscPopSignalHandler();
#endif // FEELPP_HAS_PETSC_H


    if ( argc >= 1 )
    {
        std::ostringstream ostr;
        ostr << argv[0] << ".assertions";
        Assert::setLog( ostr.str().c_str() );
    }

    S_worldcomm = worldcomm_type::New( world );
    FEELPP_ASSERT( S_worldcomm ).error( "creating worldcomm failed" );
}

Environment::~Environment()
{
    Debug(900) << "[~Environment] sending delete to all deleters" << "\n";

    // send signal to all deleters
    S_deleteObservers();
    Debug(900) << "[~Environment] delete signal sent" << "\n";

    if ( i_initialized )
    {
        Debug(900) << "[~Environment] finalizing slepc,petsc and mpi\n";
#if defined ( FEELPP_HAS_PETSC_H )
        PetscTruth is_petsc_initialized;
        PetscInitialized( &is_petsc_initialized );

        if ( is_petsc_initialized )
        {
#if defined( FEELPP_HAS_SLEPC )
            SlepcFinalize();
#else
            PetscFinalize();
#endif // FEELPP_HAS_SLEPC
        }

#endif // FEELPP_HAS_PETSC_H

        google::ShutdownGoogleLogging();
    }


}



bool
Environment::initialized()
{

#if defined( FEELPP_HAS_PETSC_H )
    PetscTruth is_petsc_initialized;
    PetscInitialized( &is_petsc_initialized );
    return mpi::environment::initialized() && is_petsc_initialized ;
#else
    return mpi::environment::initialized() ;
#endif
}

bool
Environment::finalized()
{
#if defined( FEELPP_HAS_PETSC_H )
    PetscTruth is_petsc_initialized;
    PetscInitialized( &is_petsc_initialized );
    return mpi::environment::finalized() && !is_petsc_initialized;
#else
    return mpi::environment::finalized();
#endif
}


std::string
Environment::rootRepository()
{
    std::string env;

    if ( ::getenv( "FEELPP_REPOSITORY" ) )
    {
        env = ::getenv( "FEELPP_REPOSITORY" );
    }

    else
    {
        // by default create $HOME/feel
        env = ::getenv( "HOME" );
        env += "/feel";
    }

    return env;
}
std::string
Environment::localGeoRepository()
{
    fs::path rep_path;

    rep_path = Environment::rootRepository();
    rep_path /= "geo";

    if ( !fs::exists( rep_path ) )
        fs::create_directory( rep_path );

    return rep_path.string();
}
boost::tuple<std::string,bool>
Environment::systemGeoRepository()
{
    fs::path rep_path;

    rep_path = BOOST_PP_STRINGIZE( INSTALL_PREFIX );
    rep_path /= "share/feel/geo";
    return boost::make_tuple( rep_path.string(), fs::exists( rep_path ) );
}

std::string
Environment::localConfigRepository()
{
    fs::path rep_path;

    rep_path = Environment::rootRepository();
    rep_path /= "config";

    if ( !fs::exists( rep_path ) )
        fs::create_directory( rep_path );

    return rep_path.string();
}
boost::tuple<std::string,bool>
Environment::systemConfigRepository()
{
    fs::path rep_path;

    rep_path = BOOST_PP_STRINGIZE( INSTALL_PREFIX );
    rep_path /= "share/feel/config";
    return boost::make_tuple( rep_path.string(), fs::exists( rep_path ) );
}

void
Environment::changeRepository( boost::format fmt, std::string const& logfilename )
{
    fs::path rep_path;

    rep_path = Environment::rootRepository();

    if ( !fs::exists( rep_path ) )
        fs::create_directory( rep_path );

    typedef std::vector< std::string > split_vector_type;

    split_vector_type dirs; // #2: Search for tokens
    std::string fmtstr = fmt.str();
    boost::split( dirs, fmtstr, boost::is_any_of( "/" ) );

    BOOST_FOREACH( std::string const& dir, dirs )
    {
        //Debug( 1000 ) << "[Application::Application] option: " << s << "\n";
        rep_path = rep_path / dir;

        if ( !fs::exists( rep_path ) )
            fs::create_directory( rep_path );
    }

    ::chdir( rep_path.string().c_str() );

    setLogs( logfilename );
}

po::variables_map
Environment::vm( po::options_description const& desc )
{
    po::variables_map vm;
    po::store( po::parse_command_line( 0, ( char** )0, desc ), vm );
    po::notify( vm );

    return vm;
}

void
Environment::setLogs( std::string const& prefix )
{
    mpi::communicator world;
    Log().detachAll();
    std::ostringstream ostr;
    ostr << prefix << "-" << world.size()  << "." << world.rank();
    Log().attach( ostr.str() );

    std::ostringstream ostr_assert;
    ostr_assert << prefix  << "-" << world.size()  << "." << world.rank() << ".assertions";
    Assert::setLog( ostr_assert.str().c_str() );
}

std::vector<WorldComm> const&
Environment::worldsComm( int n )
{
    if ( !S_worldcomm )
    {
        mpi::communicator world;
        S_worldcomm = worldcomm_type::New( world );
        FEELPP_ASSERT( S_worldcomm ).error( "worldcomm not allocated" );
    }

    return S_worldcomm->subWorlds(n);
}

std::vector<WorldComm> const&
Environment::worldsCommGroupBySubspace( int n )
{
#if 0
    std::cout << "n=" << n << "\n";
    S_worldcomm->showMe();
    S_worldcomm->masterWorld(n).showMe();
    std::cout << "size=" << S_worldcomm->subWorlds(n).size() <<  "\n";
    S_worldcomm->subWorlds(n).begin()->showMe();
#endif
    return S_worldcomm->subWorldsGroupBySubspace(n);
}


WorldComm const&
Environment::masterWorldComm( int n )
{
    return S_worldcomm->masterWorld(n);
}

boost::signals2::signal<void ()> Environment::S_deleteObservers;

boost::shared_ptr<WorldComm> Environment::S_worldcomm;

}
