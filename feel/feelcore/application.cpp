/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  Date: 2005-03-17

  Copyright (C) 2007-2012 Universite de Grenoble 1
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
   \file application.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-03-17
 */
#include <cstdlib>
#include <locale>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <boost/assign/list_of.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/application.hpp>

#include <feel/feelcore/feelpetsc.hpp>

#if defined(FEELPP_HAS_TRILINOS_EPETRA)
#if defined(FEELPP_HAS_MPI_H)
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif /* FEELPP_HAS_MPI_H */
#endif /* FEELPP_HAS_TRILINOS_EPETRA */
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION

#include <feel/feelcore/mpicompat.hpp>

namespace google
{
namespace glog_internal_namespace_
{
bool IsGoogleLoggingInitialized();
}
}

namespace Feel
{
namespace fs = boost::filesystem;
namespace ptree = boost::property_tree;

namespace detail
{
const int spaces = 30;
}

FEELPP_NO_EXPORT
std::pair<std::string, std::string>
at_option_parser( std::string const&s )
{
    if ( '@' == s[0] )
        return std::make_pair( std::string( "response-file" ), s.substr( 1 ) );

    else
        return std::pair<std::string, std::string>();
}


void
Application::initPETSc()
{
#if defined ( FEELPP_HAS_PETSC_H )
    //if ( M_vm["backend"].as<std::string>() == "petsc" )
    {
        PETSC_COMM_WORLD = COMM_WORLD;
        int __argc = this->unknownArgc();
        char** __argv = this->unknownArgv();
#if defined( FEELPP_HAS_SLEPC )
        int ierr = SlepcInitialize( &__argc,&__argv, PETSC_NULL, PETSC_NULL );
#else
        int ierr = PetscInitialize( &__argc, &__argv, PETSC_NULL, PETSC_NULL );
#endif
        // make sure that petsc do not catch signals and hence do not print long
        //and often unuseful messages
        PetscPopSignalHandler();
#if 0
        std::cerr << "[Application] argc " << __argc << "\n";
        --__argc;

        for ( int i = 0; i < argc; ++i )
            if ( __argv[i] )
                std::cerr << "[Application] argv[" << i << "]="<< __argv[i] << "\n";

#endif
        //int ierr = PetscInitializeNoArguments();
        boost::ignore_unused_variable_warning( ierr );
        CHKERRABORT( COMM_WORLD,ierr );
    }
#endif // FEELPP_HAS_PETSC_H

}
void
Application::initTrilinos()
{
#if defined( FEELPP_HAS_TRILINOS_EPETRA )
    //if ( M_vm["backend"].as<std::string>() == "trilinos" )
    {

    }
#endif // FEELPP_HAS_TRILINOS_EPETRA
}
void
Application::initMPI( int argc, char** argv, MPI_Comm comm )
{
#if defined( FEELPP_HAS_TBB )
    int n = tbb::task_scheduler_init::default_num_threads();
#else
    int n = 1 ;
#endif
    DVLOG(2) << "[Feel++] TBB running with " << n << " threads\n";
    //tbb::task_scheduler_init init(1);

#if defined( FEELPP_HAS_MPI_H )

    if ( !mpi::environment::initialized() )
    {
        MPI_Init ( &argc, &argv );
#if MPI_VERSION >= 2
        MPI_Comm_set_errhandler( MPI_COMM_WORLD, MPI_ERRORS_RETURN );
#else
        MPI_Errhandler_set( MPI_COMM_WORLD, MPI_ERRORS_RETURN );
#endif
    }

    M_comm = Environment::worldComm();
    //MPI_Comm_dup ( comm, &COMM_WORLD);
    //MPI_Comm_dup ( comm, (MPI_Comm*)&S_world );
#if 0
    MPI_Comm_rank ( COMM_WORLD, &_S_process_id );
    MPI_Comm_size ( COMM_WORLD, &_S_n_process );
#else
    //std::cout << "rank : " << M_comm->rank() << " " << "size : " << M_comm->size() << "\n";
    //_S_process_id = S_world.rank();
    //_S_n_process = S_world.size();
#endif
#endif // FEELPP_HAS_MPI_H

}

Application::Application()
    :
    M_about( Environment::about() ),
    M_desc( Environment::optionsDescription() ),
    M_vm( Environment::vm() )
{

}
#if defined( FEELPP_HAS_MPI_H )
MPI_Comm Application::COMM_WORLD = MPI_COMM_WORLD;

Application::Application( int argc,
                          char** argv,
                          AboutData const& ad,
                          MPI_Comm comm )
#else
Application::Application( int argc,
                          char** argv,
                          AboutData const& ad )
#endif // FEELPP_HAS_MPI_H
    :
    M_about( ad ),
    M_desc( Environment::optionsDescription() ),
    M_vm( Environment::vm() ),
    M_to_pass_further()
#if defined( FEELPP_HAS_MPI_H )
    ,
    M_env()
#endif
{
    //M_desc.add( Feel::feel_options() );
    if ( !google::glog_internal_namespace_::IsGoogleLoggingInitialized() )
    {
        // Initialize Google's logging library.
        google::InitGoogleLogging(M_about.appName().c_str());
    }

    initMPI( argc, argv, comm );

    doOptions( argc, argv );

#if defined( FEELPP_HAS_MPI_H )
    char * __env = getenv( "DEBUG" );
    std::string env_str;

    if ( __env )
        env_str = __env;

    mpi::broadcast( M_comm, env_str, 0 );

    if ( processId() != 0 )
    {
        setenv( "DEBUG", env_str.c_str(), 1 );
        //VLOG(1) << "DEBUG is set to " << env_str << "\n";
        //std::cout << "DEBUG is set to " << env_str << "\n";
    }

#endif // MPI

    initPETSc();
    initTrilinos();

#if defined(FEELPP_HAS_TAU)
    TAU_PROFILE( "Application", "Application::Application( int, char**, AboutData const&, bool)", TAU_DEFAULT );
    TAU_PROFILE_INIT( argc,argv );
    TAU_PROFILE_SET_NODE( 0 );
#endif /* FEELPP_HAS_TAU */
}

#if defined( FEELPP_HAS_MPI_H )
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
#endif // FEELPP_HAS_MPI_H
    :
    M_about( ad ),
    M_desc( Environment::optionsDescription() ),
    M_vm( Environment::vm() ),
    M_to_pass_further()
#if defined( FEELPP_HAS_MPI_H )
    ,
    M_env()
#endif

{
    if ( !google::glog_internal_namespace_::IsGoogleLoggingInitialized() )
    {
        // Initialize Google's logging library.
        google::InitGoogleLogging(M_about.appName().c_str());
    }

    //M_desc.add( Feel::feel_options() ).add( od );
    M_desc.add( od );


    initMPI( argc, argv, comm );

    //doOptions( argc, argv );

#if defined( FEELPP_HAS_MPI_H )
    char * __env = getenv( "DEBUG" );
    std::string env_str;

    if ( __env )
        env_str = __env;

    mpi::broadcast( M_comm, env_str, 0 );

    if ( processId() != 0 )
    {
        setenv( "DEBUG", env_str.c_str(), 1 );
        //VLOG(1) << "DEBUG is set to " << env_str << "\n";
        //std::cout << "DEBUG is set to " << env_str << "\n";
    }

#endif // MPI

    initPETSc();
    initTrilinos();

#if defined(FEELPP_HAS_TAU)
    TAU_PROFILE( "Application", "Application::Application( int, char**, AboutData const&, po::options_description const&, bool)", TAU_DEFAULT );
    TAU_PROFILE_INIT( argc,argv );
    TAU_PROFILE_SET_NODE( 0 );
#endif /* FEELPP_HAS_TAU */




}

#if defined( FEELPP_HAS_MPI_H )
Application::Application( AboutData const& ad,
                          po::options_description const& od,
                          MPI_Comm comm )
#else
Application::Application( AboutData const& ad,
                          po::options_description const& od )
#endif // FEELPP_HAS_MPI_H
    :
    M_about( ad ),
    M_desc( Environment::optionsDescription() ),
    M_vm( Environment::vm() ),
    M_to_pass_further()
#if defined( FEELPP_HAS_MPI_H )
    ,
    M_env()
#endif

{
    if ( !google::glog_internal_namespace_::IsGoogleLoggingInitialized() )
    {
        // Initialize Google's logging library.
        google::InitGoogleLogging(M_about.appName().c_str());
    }

    //M_desc.add( Feel::feel_options() ).add( od );
    M_desc.add( od );

    //
    // if we are using openmpi, then we need to dlopen mpi with some special flags
    // to avoid some undefined symbol when using Application classes in Python code
    // for example.
    //
#if OPEN_MPI
    OPENMPI_dlopen_libmpi();
#endif

    int argc = 1;
    char** argv = new char*[argc];
    argv[0] = new char[M_about.appName().size()+1];
    ::strcpy( argv[0], M_about.appName().c_str() );


    initMPI( argc, argv, comm );

    //doOptions( argc, argv );

#if defined( FEELPP_HAS_MPI_H )
    char * __env = getenv( "DEBUG" );
    std::string env_str;

    if ( __env )
        env_str = __env;

    mpi::broadcast( M_comm, env_str, 0 );

    if ( processId() != 0 )
    {
        setenv( "DEBUG", env_str.c_str(), 1 );
        //VLOG(1) << "DEBUG is set to " << env_str << "\n";
        //std::cout << "DEBUG is set to " << env_str << "\n";
    }

#endif // MPI

    initPETSc();
    initTrilinos();

#if defined(FEELPP_HAS_TAU)
    TAU_PROFILE( "Application", "Application::Application( int, char**, AboutData const&, po::options_description const&, bool)", TAU_DEFAULT );
    TAU_PROFILE_INIT( argc,argv );
    TAU_PROFILE_SET_NODE( 0 );
#endif /* FEELPP_HAS_TAU */




}

#if defined( FEELPP_HAS_MPI_H )
Application::Application( AboutData const& ad,
                          MPI_Comm comm )
#else
Application::Application( AboutData const& ad )

#endif // FEELPP_HAS_MPI_H
    :
    M_about( ad ),
    M_desc( Environment::optionsDescription() ),
    M_vm( Environment::vm() ),
    M_to_pass_further()
#if defined( FEELPP_HAS_MPI_H )
    ,
    M_env()
#endif

{
    if ( !google::glog_internal_namespace_::IsGoogleLoggingInitialized() )
    {
        // Initialize Google's logging library.
        google::InitGoogleLogging(M_about.appName().c_str());
    }
#if 1

    //
    // if we are using openmpi, then we need to dlopen mpi with some special flags
    // to avoid some undefined symbol when using Application classes in Python code
    // for example.
    //
#if OPEN_MPI
    OPENMPI_dlopen_libmpi();
#endif

    int argc = 1;
    char** argv = new char*[argc];
    argv[0] = new char[M_about.appName().size()+1];
    ::strcpy( argv[0], M_about.appName().c_str() );


    initMPI( argc, argv, comm );
#endif

#if defined( FEELPP_HAS_MPI_H )
    char * __env = getenv( "DEBUG" );
    std::string env_str;

    if ( __env )
        env_str = __env;

    mpi::broadcast( M_comm, env_str, 0 );

    if ( processId() != 0 )
    {
        setenv( "DEBUG", env_str.c_str(), 1 );
        //VLOG(1) << "DEBUG is set to " << env_str << "\n";
        //std::cout << "DEBUG is set to " << env_str << "\n";
    }

#endif // MPI

    initPETSc();
    initTrilinos();

#if defined(FEELPP_HAS_TAU)
    TAU_PROFILE( "Application", "Application::Application( int, char**, AboutData const&, po::options_description const&, bool)", TAU_DEFAULT );
    TAU_PROFILE_INIT( argc,argv );
    TAU_PROFILE_SET_NODE( 0 );
#endif /* FEELPP_HAS_TAU */

}
Application::Application( Application const& __app )
    :
    M_about( __app.M_about ),
    M_desc( __app.M_desc ),
    M_vm( __app.M_vm ),
    M_to_pass_further( __app.M_to_pass_further )
{
}
Application::~Application()
{
#if 0
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

#endif // 0
}
po::options_description
benchmark_options( std::string const& prefix  )
{
    po::options_description benchopt( "Benchmarking options" );
    benchopt.add_options()
    ( prefixvm( prefix,"benchmark.nlevels" ).c_str(), po::value<int>()->default_value( 1 ), "number of mesh levels to benchmark" )
    ( prefixvm( prefix,"benchmark.hsize" ).c_str(), po::value<double>()->default_value( 0.1 ), "default mesh size" )
    ( prefixvm( prefix,"benchmark.refine" ).c_str(), po::value<double>()->default_value( 2 ), "refine ratio for meshes" )
    ( prefixvm( prefix,"benchmark.prepare" ).c_str(), po::value<bool>()->default_value( false ), "prepare data for the benchmark (e.g. precompute meshes/partitioning)" )
    ( prefixvm( prefix,"benchmark.partitions" ).c_str(), po::value<int>()->default_value( 1 ), "number of partitions used for the benchmark" )
    ;

    // make sense only for global benchmark options
    if ( prefix.empty() )
        benchopt.add_options()
        ( "benchmark.only", po::value<std::string>()->default_value( "" ), "benchmarks to run, empty means all" )
        ;

    return benchopt;
}
void
Application::doOptions( int argc, char** argv )
{
    //std::locale::global(std::locale(""));
    try
    {
        po::options_description generic( "Generic options" );
        generic.add_options()
        ( "authors", "prints the authors list" )
        ( "copyright", "prints the copyright statement" )
        ( "help", "prints this help message" )
        ( "license", "prints the license text" )
        ( "version,v", "prints the version" )
        ( "feelinfo", "prints feel libraries information" )
        ( "nochdir", "Don't change repository directory even though it is called" )
        ( "config-file", po::value<std::string>(), "specify .cfg file" )
        ( "result-file", po::value<std::string>()->default_value(this->about().appName()), "specify .res file" )
        ( "response-file", po::value<std::string>(), "can be specified with '@name', too" )
        ;
        po::options_description debug( "Debugging options" );
        debug.add_options()
        ( "debug", po::value<std::string>()->default_value( "" ), "specify a debugging area list" );
        M_desc.add( generic ).add( debug ).add( benchmark_options() );

        this->parseAndStoreOptions( po::command_line_parser( argc, argv ), true );
        processGenericOptions();

        /**
         * parse config file if given to command line
         */
        if ( M_vm.count( "config-file" ) )
        {
            DVLOG(2) << "[Application] parsing " << M_vm["config-file"].as<std::string>() << "\n";

            if ( fs::exists(  M_vm["config-file"].as<std::string>() ) )
            {

                std::ifstream ifs( M_vm["config-file"].as<std::string>().c_str() );
                po::store( parse_config_file( ifs, M_desc, true ), M_vm );
                po::notify( M_vm );
            }
        }

        std::vector<fs::path> prefixes = boost::assign::list_of( fs::current_path() )
                                         ( fs::path ( Environment::localConfigRepository() ) )
                                         ( fs::path ( Environment::systemConfigRepository().get<0>() ) )
                                         ( fs::path ( "/usr/share/feel/config" ) )
                                         ( fs::path ( "/usr/local/share/feel/config" ) )
                                         ( fs::path ( "/opt/local/share/feel/config" ) );

        BOOST_FOREACH( auto prefix, prefixes )
        {
            std::string config_name = ( boost::format( "%1%/%2%.cfg" ) % prefix.string() % this->about().appName() ).str();
            DVLOG(2) << "[Application] Looking for " << config_name << "\n";
            DVLOG(2) << "[Application] Looking for " << config_name << "\n";

            if ( fs::exists( config_name ) )
            {
                DVLOG(2) << "[Application] parsing " << config_name << "\n";
                std::ifstream ifs( config_name.c_str() );
                store( parse_config_file( ifs, M_desc, true ), M_vm );
                break;
            }

            else
            {
                // try with a prefix feel_
                std::string config_name = ( boost::format( "%1%/feel_%2%.cfg" ) % prefix.string() % this->about().appName() ).str();
                DVLOG(2) << "[Application] Looking for " << config_name << "\n";

                if ( fs::exists( config_name ) )
                {
                    DVLOG(2) << "[Application] loading configuration file " << config_name << "...\n";
                    std::ifstream ifs( config_name.c_str() );
                    store( parse_config_file( ifs, M_desc, true ), M_vm );
                    break;
                }
            }
        }
        //po::store(po::parse_command_line(argc, argv, M_desc), M_vm);
        po::notify( M_vm );

    }

    // catches program_options exceptions
    catch ( boost::program_options::multiple_occurrences const& e )
    {
        std::cout << "Command line or config file option parsing error: " << e.what() << "\n"
                  << "  o faulty option: " << e.get_option_name() << "\n"
                  << "Warning: the .cfg file or some options may not have been read properly\n";
    }

    catch ( boost::program_options::ambiguous_option const& e )
    {
        std::cout << "Command line or config file option parsing error: " << e.what() << "\n"
                  << "  o faulty option: " << e.get_option_name() << "\n"
                  << "  o possible alternatives: " ;
        std::for_each( e.alternatives().begin(), e.alternatives().end(), []( std::string const& s )
        {
            std::cout << s << " ";
        } );
        std::cout << "\n"
                  << "Warning: the .cfg file or some options may not have been read properly\n";
    }

    // catches program_options exceptions
    catch ( std::exception& e )
    {
        std::cout << "Application option parsing: unknown option:" << e.what() << " (the .cfg file or some options may not have been read properly)\n";
    }

    catch ( ... )
    {
        std::cout << "Application option parsing: unknown exception triggered  (the .cfg file or some options may not have been read properly)\n";
    }
}
char**
Application::unknownArgv() const
{
    char** argv = new char*[ M_to_pass_further.size()+1 ];
    argv[0] = new char[ std::strlen( this->about().appName().c_str() )+1 ];
    strcpy( argv[0], this->about().appName().c_str() );
    argv[0][std::strlen( this->about().appName().c_str() )] = '\0';
    int n_a = 0;
    DVLOG(2) << "argv[ " << n_a << " ]=" << argv[0] << "\n";
    ++n_a;
    BOOST_FOREACH( std::string const& s, M_to_pass_further )
    {
        size_type ssize=s.size();
        DVLOG(2) << "new arg " << s << " size = " << ssize << "\n";
        argv[n_a] = new char[ s.size()+1 ];
        strncpy( argv[n_a], s.c_str(), s.size() );
        argv[n_a][s.size()] = '\0';
        ++n_a;

        DVLOG(2) << "argv[ " << n_a-1 << " ]=" << argv[n_a-1] << "\n";
    }
    return argv;
}
void
Application::setName1( std::string const& name1 )
{
    M_name1 = name1;
}

void
Application::setName2( std::string const& name2 )
{
    M_name2 = name2;
}

void
Application::setH( double h, int precision )
{
    M_h = std::make_pair( h, precision );
}

void
Application::setDimension( int dim )
{
    M_dim = dim;
}
std::string
Application::resultFileName() const
{
    return this->vm()["result-file"].as<std::string>();
}
void
Application::setLogs()
{
}
void
Application::processGenericOptions()
{
    // leave this to subclasses or users
#if 0
    if ( M_vm.count( "help" ) )
        std::cout << M_desc << "\n";

#endif


    if ( M_vm.count( "response-file" ) )
    {
        using namespace std;
        // Load the file and tokenize it
        ifstream ifs( M_vm["response-file"].as<std::string>().c_str() );

        if ( !ifs )
        {
            cout << "Could not open the response file\n";
            return ;
        }

        // Read the whole file into a string
        stringstream ss;
        ss << ifs.rdbuf();
        // Split the file content
        boost::char_separator<char> sep( " \n\r" );
        boost::tokenizer<boost::char_separator<char> > tok( ss.str(), sep );
        vector<string> args;
        copy( tok.begin(), tok.end(), back_inserter( args ) );

        this->parseAndStoreOptions( po::command_line_parser( args ) );
    }

    if ( this->comm().rank() == 0 )
    {

        if ( M_vm.count( "feelinfo" ) )
            std::cout << std::setw( 15 ) << std::right << "Feel Version : " << Info::versionString() << "\n"
                      << std::setw( 15 ) << std::right << "Major : " << Info::versionMajor() << "\n"
                      << std::setw( 15 ) << std::right << "Minor : " << Info::versionMinor() << "\n"
                      << std::setw( 15 ) << std::right << "Micro : " << Info::versionMicro() << "\n"
                      << std::setw( 15 ) << std::right << "Revision : " << Info::revision() << "\n"
                      << std::setw( 15 ) << std::right << "BuildId : " << Info::buildId() << "\n"
                      << std::setw( 15 ) << std::right << "Feel Prefix : " << Info::prefix() << "\n"
                      << std::setw( 15 ) << std::right << "Feel DataDir : " << Info::datadir() << "\n";

        if ( M_vm.count( "help" ) ||
             M_vm.count( "version" ) ||
             M_vm.count( "copyright" ) ||
             M_vm.count( "license" ) ||
             M_vm.count( "authors" ) )
        {
            std::cout << M_about.appName() << ": " << M_about.shortDescription() <<  "\n";
        }

        if ( M_vm.count( "version" ) )
        {
            std::cout << " version : " << M_about.version() << "\n";
        }

        if ( M_vm.count( "copyright" ) )
        {
            std::cout << " copyright : " << M_about.copyrightStatement() << "\n";
        }

        if ( M_vm.count( "license" ) )
        {
            std::cout << " license : " << M_about.license() << "\n";
        }

        if ( M_vm.count( "authors" ) )
        {
            std::cout << std::setw( 30 )
                      << "Author Name"
                      << " " << std::setw( 15 )
                      << "Task"
                      << " " << std::setw( 40 )
                      << "Email Address"
                      << "\n";
            std::cout << std::setw( 85+3 ) << std::setfill( '-' ) << "\n" << std::setfill( ' ' );
            std::for_each( M_about.authors().begin(),
                           M_about.authors().end(),
                           std::cout
                           << std::setw( 30 )
                           << lambda::bind( &AboutPerson::name,
                                            lambda::_1 )
                           << " " << std::setw( 15 )
                           << lambda::bind( &AboutPerson::task,
                                            lambda::_1 )
                           << " " << std::setw( 40 )
                           << lambda::bind( &AboutPerson::emailAddress,
                                            lambda::_1 )
                           << "\n" );
        }

        if ( M_vm.count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
        }


    }
    if ( M_vm.count( "help" ) ||
         M_vm.count( "version" ) ||
         M_vm.count( "copyright" ) ||
         M_vm.count( "license" ) ||
         M_vm.count( "authors" ) )
    {
        this->comm().barrier();
        MPI_Finalize();
        exit(0);
    }
#if 0
    std::cout << "count = " << M_vm.count( "debug" ) << "\n"
              << "string = " << M_vm["debug"].as<std::string>() << "\n";
#endif

    DVLOG(2) << "[Application::processGenericOptions] done\n";
}
std::string
Application::rootRepository() const
{
    return Environment::rootRepository();
}

Application&
Application::changeRepository( boost::format fmt )
{
    if ( M_vm.count( "nochdir" ) )
    {
        this->setLogs();
        return *this;
    }

    Environment::changeRepository( fmt, M_about.appName() );
    this->setLogs();
    return *this;
}
void
Application::parseAndStoreOptions( po::command_line_parser parser, bool extra_parser )
{
    DVLOG(2) << "[Application::Application] parsing options...\n";

    boost::shared_ptr<po::parsed_options> parsed;

    if ( extra_parser )
    {
        parsed = boost::shared_ptr<po::parsed_options>( new po::parsed_options( parser
                 .options( M_desc )
                 .extra_parser( at_option_parser )
                 .allow_unregistered()
                 .run() ) );
    }

    else
    {
        parsed = boost::shared_ptr<po::parsed_options>( new po::parsed_options( parser
                 .options( M_desc )
                 .allow_unregistered()
                 .run() ) );
    }

    DVLOG(2) << "[Application::Application] parsing options done\n";

    M_to_pass_further = po::collect_unrecognized( parsed->options, po::include_positional );
    DVLOG(2) << "[Application::Application] number of unrecognized options: " << ( M_to_pass_further.size() ) << "\n";

    BOOST_FOREACH( std::string const& s, M_to_pass_further )
    {
        DVLOG(2) << "[Application::Application] option: " << s << "\n";
    }
    std::vector<po::basic_option<char> >::iterator it = parsed->options.begin();
    std::vector<po::basic_option<char> >::iterator en  = parsed->options.end();

    for ( ; it != en ; ++it )
        if ( it->unregistered )
        {
            DVLOG(2) << "[Application::Application] remove from vector " << it->string_key << "\n";
            parsed->options.erase( it );
        }

    po::store( *parsed, M_vm );
}

void
Application::add( Simget* simget )
{
    M_simgets.push_back( simget );
}
std::vector<ptree::ptree>&
updateStats( std::vector<ptree::ptree>& stats )
{
    try
    {
        for ( auto it = stats.begin(), en =stats.end(); it!=en; ++it )
        {
            double h  = it->get<double>( "h" );
            int l  = it->get<int>( "level" );
            double hp = h;
            if ( l > 1 )
                hp = boost::prior( it )->get<double>( "h" );

            auto ept = it->get_child( "e" );
            // loop on norms
            for( auto itchild = ept.begin(), enchild = ept.end(); itchild != enchild; ++itchild )
            {

                // for each norm compute the rate of convergence (roc)
                for( auto itv = itchild->second.begin(), env = itchild->second.end(); itv != env; ++itv )
                {
                    std::ostringstream keystr;
                    keystr << "e." << itchild->first << "." << itv->first;
                    std::string key = keystr.str();
                    //std::cout << "update key " <<key << " ...\n";

                    double u  = 1;
                    try {
                        u = it->get<double>( key );
                    }
                    catch( ptree::ptree_bad_data const& e )
                    {
                        std::cout << "Application::updateStats: (ptree::bad_data) : "  << e.what() << "\n";
                        u = 1;
                    }

                    double up = u;
                    double roc = 1;

                    if ( l > 1 )
                    {
                        up  = boost::prior( it )->get<double>( key );
                        roc = std::log10( up/u )/std::log10( hp/h );
                        //std::cout << "u = "  << u << " up = " << up << " roc=" << roc <<" \n";
                    }
                    // create roc entry in the last statistics
                    stats.back().put( key+".roc", roc );
                }
            }
        }
    }
    catch( ptree::ptree_bad_data const& e )
    {
        std::cout << "ptree::bad_data : "  << e.what() << "\n";
    }
    catch( ptree::ptree_bad_path const& e )
    {
        std::cout << "ptree::bad_path : "  << e.what() << "\n";
    }
    catch( ... )
    {

    }
    return stats;
}
void
Application::run()
{
    std::string runonly = M_vm["benchmark.only"].as<std::string>();
    bool prepare = M_vm["benchmark.prepare"].as<bool>();

    // get current work directory before a simget eventually change the working
    // directory to store its results. the statistics files will be stored in the
    // same directory as the executable for now
    fs::path cp = fs::current_path();

    for ( auto i = M_simgets.begin(), end = M_simgets.end(); i != end; ++i )
    {
        if ( ( runonly.empty() == false  ) &&
                runonly.find( i->name() ) == std::string::npos )
            continue;

        std::string s1 = prefixvm( i->name(),"benchmark.nlevels" );
        int nlevels = M_vm.count( s1 )?M_vm[s1].as<int>():M_vm["benchmark.nlevels"].as<int>();
        std::string s2 = prefixvm( i->name(),"benchmark.hsize" );
        double hsize = M_vm.count( s2 )?M_vm[s2].as<double>():M_vm["benchmark.hsize"].as<double>();
        std::string s3 = prefixvm( i->name(),"benchmark.refine" );
        double refine = M_vm.count( s3 )?M_vm[s3].as<double>():M_vm["benchmark.refine"].as<double>();
        i->setMeshSizeInit( hsize );
        bool has_stats = false;
        for ( int l = 0; l < nlevels; ++l )
        {
            double meshSize= hsize/std::pow( refine,l );
            i->setMeshSize( meshSize );
            i->setLevel( l+1 );
            i->run();
            if ( !prepare && !i->stats().empty() )
                {
                    has_stats = true;

                    i->stats().put( "h",i->meshSize() );
                    i->stats().put( "level", i->level() );

                    M_stats[i->name()].push_back( i->stats() );
                    updateStats( M_stats[i->name()] );
                    this->printStats( std::cout, Application::ALL );
                }

        }
        if ( !prepare && has_stats == true )
        {
            std::string fname = (boost::format( "%1%-%2%.tsv" )% i->name()% Environment::numberOfProcessors() ).str();
            fs::ofstream ofs( cp / fname );
            std::string fnameall = (boost::format( "%1%-%2%-all.tsv" )% i->name()% Environment::numberOfProcessors() ).str();
            fs::ofstream ofsall( cp / fnameall );
            std::string fnameerrors = (boost::format( "%1%-%2%-errors.tsv" )% i->name()% Environment::numberOfProcessors() ).str();
            fs::ofstream ofserrors( cp / fnameerrors );
            std::string fnametime = (boost::format( "%1%-%2%-timings.tsv" )% i->name()% Environment::numberOfProcessors() ).str();
            fs::ofstream ofstime( cp / fnametime );
            std::string fnamedata = (boost::format( "%1%-%2%-data.tsv" )% i->name()% Environment::numberOfProcessors() ).str();
            fs::ofstream ofsdata( cp / fnamedata );

            this->printStats( ofs, Application::ALL );
            this->printStats( ofsall, Application::ALL|Application::FLAT );
            this->printStats( ofserrors, Application::ERRORS|Application::FLAT );
            this->printStats( ofstime, Application::TIME|Application::FLAT );
            this->printStats( ofsdata, Application::DATA|Application::NUMBERS|Application::FLAT );
        }
    }
}

void
Application::run( const double* X, unsigned long P, double* Y, unsigned long N )
{
    auto it = M_simgets.begin(), en = M_simgets.end();
    int i = 0;

    for ( ; it != en; ++i, ++it )
    {
        it->run( &X[i*P], P, &Y[i*N], N );
    }
}

template<typename StatsIterator>
void
printErrors( std::ostream& out,
             StatsIterator statsit,
             StatsIterator statsen,
             std::string const& key,
             size_type stats = Application::ALL|Application::HEADER )
{
    if ( !( stats & Application::ERRORS ) ) return;
    bool header = stats & Application::HEADER;
    bool flat = stats & Application::FLAT;
    if ( header )
    {
        if ( !flat )
            out << std::setw( 10 ) << std::right << "levels"
                << std::setw( 10 ) << std::right << "h";
        BOOST_FOREACH( auto v, statsit->get_child( key ) )
        {
            if ( !flat )
                out << std::setw( detail::spaces ) << std::right << v.first
                    << std::setw( detail::spaces ) << std::right << v.first+".roc";
            else
                out << std::setw( detail::spaces ) << std::right << key+"."+v.first
                    << std::setw( detail::spaces ) << std::right << key+"."+v.first+".roc";
        }
        if ( !flat )
            out << "\n";
    }

    if ( ! ( (header == true) && ( flat == true ) ) )
        for ( auto it = statsit, en =statsen; it!=en; ++it )
        {
            //std::for_each( it->begin(),it->end(), []( std::pair<std::string,boost::any> const& o ) { std::cout << o.first << "\n"; } );
            //std::map<std::string,boost::any> data = *it;
            //std::map<std::string,boost::any> datap;
            double h  = it->template get<double>( "h" );
            int l  = it->template get<double>( "level" );
            double hp = h;
            if ( !flat )
                out << std::right << std::setw( 10 ) << l
                    << std::right << std::setw( 10 ) << std::fixed  << std::setprecision( 4 ) << h;
            BOOST_FOREACH( auto v, it->get_child( key ) )
            {
                double u  = it->template get<double>( key+"."+v.first );
                double roc  = it->template get<double>( key+"."+v.first+".roc" );

                out << std::right << std::setw( detail::spaces ) << std::scientific << std::setprecision( 2 ) << u
                    << std::right << std::setw( detail::spaces ) << std::fixed << std::setprecision( 2 ) << roc;
            }
            if ( !flat )
                out << "\n";
        }
}
template<typename StatsIterator>
void
printNumbers( std::ostream& out,
              StatsIterator statsit,
              StatsIterator statsen,
              std::string const& key,
              size_type stats = Application::ALL|Application::HEADER )
{
    if ( !( stats & Application::DATA ) ) return;
    bool header = stats & Application::HEADER;
    bool flat = stats & Application::FLAT;
    if ( header )
    {
        if ( !flat )
            out << std::setw( 10 ) << std::right << "levels"
                << std::setw( 10 ) << std::right << "h";
        BOOST_FOREACH( auto v, statsit->get_child( key ) )
        {
            if ( !flat )
                out << std::setw( detail::spaces ) << std::right << v.first;
            else
                out << std::setw( detail::spaces ) << std::right << key+"."+v.first;
        }
        if ( !flat )
            out << "\n";
    }
    if ( ! ( (header == true) && ( flat == true ) ) )
        for ( auto it = statsit, en =statsen; it!=en; ++it )
        {
            double h  = it->template get<double>( "h" );
            int l  = it->template get<double>( "level" );
            if ( !flat )
                out << std::right << std::setw( 10 ) << l
                    << std::right << std::setw( 10 ) << std::fixed  << std::setprecision( 4 ) << h;
            BOOST_FOREACH( auto v, it->get_child( key ) )
            {
                size_type u  = it->template get<size_type>( key+"."+v.first );
                out << std::right << std::setw( detail::spaces )  << u;
            }
            if ( !flat )
                out << "\n";
        }
}
template<typename StatsIterator>
void
printData( std::ostream& out,
           StatsIterator statsit,
           StatsIterator statsen,
           std::string const& key,
           size_type stats = Application::ALL|Application::HEADER )
{
    try
    {
        if ( !( stats & Application::DATA ) ) return;
        bool header = stats & Application::HEADER;
        bool flat = stats & Application::FLAT;
        if ( header )
        {
            if ( !flat )
                out << std::setw( 10 ) << std::right << "levels"
                    << std::setw( 10 ) << std::right << "h";
            try {
                BOOST_FOREACH( auto v, statsit->get_child( key+".bool" ) )
                {
                    out << std::setw( detail::spaces ) << std::right << (flat?key+"."+v.first:v.first);
                }
            }
            catch(...)
            {}
            try {
                BOOST_FOREACH( auto v, statsit->get_child( key+".int" ) )
                {
                    out << std::setw( detail::spaces ) << std::right << (flat?key+"."+v.first:v.first);
                }
            }
            catch(...)
            {}
            try{
                BOOST_FOREACH( auto v, statsit->get_child( key+".double" ) )
                {
                    out << std::setw( detail::spaces ) << std::right << (flat?key+"."+v.first:v.first);
                }
            }
            catch(...)
            {}
            if ( !flat )
                out << "\n";
        }
        int l=1;
        if ( ! ( (header == true) && ( flat == true ) ) )
            for ( auto it = statsit, en =statsen; it!=en; ++it,++l )
            {
                double h  = it->template get<double>( "h" );
                if ( !flat )
                    out << std::right << std::setw( 10 ) << l
                        << std::right << std::setw( 10 ) << std::fixed  << std::setprecision( 4 ) << h;
                try {
                    BOOST_FOREACH( auto v, it->get_child( key+".bool" ) )
                    {
                        bool u  = it->template get<bool>( key+".bool."+v.first );
                        out << std::right << std::setw( detail::spaces )  << u;
                    }
                }
                catch(...){}
                try {
                    BOOST_FOREACH( auto v, it->get_child( key+".int" ) )
                    {
                        size_type u  = it->template get<size_type>( key+".int."+v.first );
                        out << std::right << std::setw( detail::spaces )  << u;
                    }
                }
                catch(...){}
                try {
                    BOOST_FOREACH( auto v, it->get_child( key+".double" ) )
                    {
                        double u  = it->template get<double>( key+".double."+v.first );
                        out << std::right << std::setw( detail::spaces ) << std::scientific << std::setprecision( 2 ) << u;
                    }
                }
                catch(...){}
                if ( !flat )
                    out << "\n";
            }
    }
    catch( ptree::ptree_bad_data const& e )
    {
        LOG(WARNING) << "Invalid property tree data : " << e.what();
    }
}
template<typename StatsIterator>
void
printTime( std::ostream& out,
           StatsIterator statsit,
           StatsIterator statsen,
           std::string const& key,
           size_type stats = Application::ALL|Application::HEADER )
{
    try {

        if ( !( stats & Application::TIME ) ) return;
        bool header = stats & Application::HEADER;
        bool flat = stats & Application::FLAT;
        if ( header )
        {
            if ( !flat )
                out << std::setw( 10 ) << std::right << "levels"
                    << std::setw( 10 ) << std::right << "h";
            BOOST_FOREACH( auto v, statsit->get_child( key ) )
            {
                int len1 = std::max( detail::spaces,( int )v.first.size() );
                int len2 = std::max( detail::spaces,( int )( std::string( " normalized" ).size() ) );
                out << std::setw( len1 ) << std::right << (flat?key+"."+v.first:v.first) << " "
                    << std::setw( len2 ) << std::right << (flat?key+"."+v.first+".n":v.first+".n");
            }
            if (!flat)
                out << "\n";
        }
        int l=1;
        std::map<std::string,double> t0;
        if ( ! ( (header == true) && ( flat == true ) ) )
            for ( auto it = statsit,  en =statsen; it!=en; ++it,++l )
            {
                double h  = it->template get<double>( "h" );
                if ( !flat )
                    out << std::right << std::setw( 10 ) << l
                        << std::right << std::setw( 10 ) << std::fixed  << std::setprecision( 4 ) << h;
                BOOST_FOREACH( auto v, it->get_child( key ) )
                {
                    std::string thekey = key+"."+v.first;
                    int len1 = std::max( detail::spaces,( int )(flat?key+"."+v.first:v.first).size() );
                    int len2 = std::max( detail::spaces,(int)std::string( " normalized" ).size()-2);
                    double u  = it->template get<double>( thekey );

                    if ( l == 1 )
                        t0[thekey] = u;

                    out << std::right << std::setw( len1 ) << std::scientific << std::setprecision( 2 ) << u << " "
                        << std::right << std::setw( len2 ) << std::scientific << std::setprecision( 2 ) << u/t0[thekey];
                }
                if ( !flat )
                    out << "\n";
            }
    }
    catch( ptree::ptree_bad_data const& e )
    {
        LOG(WARNING) << "Invalid property tree data : " << e.what();
    }
}

void
Application::setStats( std::vector<std::string> const& keys )
{
    M_keys = keys;
}
void
Application::storeStats( std::string const&  n, ptree::ptree const& s )
{
    M_stats[n].push_back( s );
}
void
Application::printStats( std::ostream& out, size_type stats ) const
{
    printStats( out, M_keys, stats );
}
void
Application::printStats( std::ostream& out,
                         std::vector<std::string> const& keys,
                         size_type stats ) const
{
    bool header = stats & Application::HEADER;
    bool flat = stats & Application::FLAT;
    if ( keys.empty() ) return;
    if ( M_comm.rank() != 0 ) return ;
    std::string runonly = M_vm["benchmark.only"].as<std::string>();
    bool prepare = M_vm["benchmark.prepare"].as<bool>();
    if ( prepare ) return;
    for ( auto i = M_simgets.begin(), end = M_simgets.end(); i != end; ++i )
    {
        if ( ( runonly.empty() == false  ) &&
                runonly.find( i->name() ) == std::string::npos )
            continue;

        if ( !flat )
        {
            out << "================================================================================\n";
            out << "Simulation " << i->name() << "\n";
        }
        BOOST_FOREACH( auto key, keys )
        {
            if ( flat == false )
            {
                out << "------------------------------------------------------------\n";
                out << "Key: " << key << "\n";
            }
            if ( M_stats.find( i->name() ) != M_stats.end() )
            {
                if ( !flat )
                {
                    auto it = M_stats.find( i->name() )->second.begin();
                    auto en = M_stats.find( i->name() )->second.end();
                    if ( key.find( "e." ) != std::string::npos )
                        printErrors( out, it, en, key, stats|Application::HEADER );
                    if ( key.find( "n." ) != std::string::npos )
                        printNumbers( out, it, en, key, stats|Application::HEADER );
                    if ( key.find( "t." ) != std::string::npos )
                        printTime( out, it, en, key, stats|Application::HEADER );
                    if ( key.find( "d." ) != std::string::npos )
                        printData( out, it, en, key, stats|Application::HEADER );
                }
            }
        }
        if ( flat && M_stats.find( i->name() ) != M_stats.end() )
        {
            // first print header
            auto it = M_stats.find( i->name() )->second.begin();
            auto en = M_stats.find( i->name() )->second.end();
            out << std::setw( 10 ) << std::right << "levels"
                << std::setw( 10 ) << std::right << "h";
            for( auto const& key : keys )
            {
                if ( key.find( "e." ) != std::string::npos )
                    printErrors( out, it, en, key, stats|Application::HEADER );
                if ( key.find( "n." ) != std::string::npos )
                    printNumbers( out, it, en, key, stats|Application::HEADER );
                if ( key.find( "t." ) != std::string::npos )
                    printTime( out, it, en, key, stats|Application::HEADER );
                if ( key.find( "d." ) != std::string::npos )
                    printData( out, it, en, key, stats|Application::HEADER );
            }
            out << "\n";
            // then print data
            it = M_stats.find( i->name() )->second.begin();
            en = M_stats.find( i->name() )->second.end();
            for ( ; it != en; ++it )
            {
                double h  = it->get<double>( "h" );
                int l  = it->get<int>( "level" );
                out << std::right << std::setw( 10 ) << l
                    << std::right << std::setw( 10 ) << std::fixed  << std::setprecision( 4 ) << h;
                for(auto const& key: keys )
                {
                    if ( key.find( "e." ) != std::string::npos )
                        printErrors( out, it, boost::next(it), key, stats&(~Application::HEADER) );
                    if ( key.find( "n." ) != std::string::npos )
                        printNumbers( out, it, boost::next(it), key, stats&(~Application::HEADER) );
                    if ( key.find( "t." ) != std::string::npos )
                        printTime( out, it, boost::next(it), key, stats&(~Application::HEADER) );
                    if ( key.find( "d." ) != std::string::npos )
                        printData( out, it, boost::next(it), key, stats&(~Application::HEADER) );
                }
                out << "\n";
            }
        }
    }
}
}
