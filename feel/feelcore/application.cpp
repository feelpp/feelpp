/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
  Date: 2005-03-17

  Copyright (C) 2007-2010 Universite de Grenoble 1
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-03-17
 */
#include <cstdlib>

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

#if defined(HAVE_TRILINOS_EPETRA)
#if defined(HAVE_MPI_H)
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif /* HAVE_MPI_H */
#endif /* HAVE_TRILINOS_EPETRA */
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION

#include <feel/feelcore/mpicompat.hpp>


namespace Feel
{
namespace fs = boost::filesystem;
namespace ptree = boost::property_tree;


FEEL_NO_EXPORT
std::pair<std::string, std::string>
at_option_parser(std::string const&s)
{
    if ('@' == s[0])
        return std::make_pair(std::string("response-file"), s.substr(1));
    else
        return std::pair<std::string, std::string>();
}


void
Application::initPETSc()
{
#if defined ( HAVE_PETSC_H )
    //if ( _M_vm["backend"].as<std::string>() == "petsc" )
        {
            PETSC_COMM_WORLD = COMM_WORLD;
            int __argc = this->unknownArgc();
            char** __argv = this->unknownArgv();
#if defined( FEEL_HAVE_SLEPC )
            int ierr = SlepcInitialize(&__argc,&__argv, PETSC_NULL, PETSC_NULL );
#else
            int ierr = PetscInitialize( &__argc, &__argv, PETSC_NULL, PETSC_NULL );
#endif
            // make sure that petsc do not catch signals and hence do not print long
            //and often unuseful messages
            PetscPopSignalHandler();
#if 0
            std::cerr << "[Application] argc " << __argc << "\n";
            --__argc;
            for( int i = 0; i < argc; ++i )
                if ( __argv[i] )
                    std::cerr << "[Application] argv[" << i << "]="<< __argv[i] << "\n";
#endif
            //int ierr = PetscInitializeNoArguments();
            boost::ignore_unused_variable_warning(ierr);
            CHKERRABORT(COMM_WORLD,ierr);
        }
#endif // HAVE_PETSC_H

}
void
Application::initTrilinos()
{
#if defined( HAVE_TRILINOS_EPETRA )
    //if ( _M_vm["backend"].as<std::string>() == "trilinos" )
        {

        }
#endif // HAVE_TRILINOS_EPETRA
}
void
Application::initMPI( int argc, char** argv, MPI_Comm comm )
{
#if defined( HAVE_TBB )
    int n = tbb::task_scheduler_init::default_num_threads();
#else
    int n = 1 ;
#endif
    Debug( 1000 ) << "[Feel++] TBB running with " << n << " threads\n";
    //tbb::task_scheduler_init init(1);

#if defined( HAVE_MPI_H )
    if (!mpi::environment::initialized())
    {
        MPI_Init (&argc, &argv);
        MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
    }
    M_comm = boost::shared_ptr<mpi::communicator>( new mpi::communicator() );
    //MPI_Comm_dup ( comm, &COMM_WORLD);
    //MPI_Comm_dup ( comm, (MPI_Comm*)&S_world );
#if 0
    MPI_Comm_rank (COMM_WORLD, &_S_process_id);
    MPI_Comm_size (COMM_WORLD, &_S_n_process);
#else
    //std::cout << "rank : " << M_comm->rank() << " " << "size : " << M_comm->size() << "\n";
    //_S_process_id = S_world.rank();
    //_S_n_process = S_world.size();
#endif
#endif // HAVE_MPI_H

}

#if defined( HAVE_MPI_H )
MPI_Comm Application::COMM_WORLD = MPI_COMM_WORLD;

Application::Application( int argc,
                          char** argv,
                          AboutData const& ad,
                          MPI_Comm comm )
#else
Application::Application( int argc,
                          char** argv,
                          AboutData const& ad )
#endif // HAVE_MPI_H
    :
_M_about( ad ),
_M_desc( "Allowed options" ),
_M_vm(),
_M_to_pass_further()
#if defined( HAVE_MPI_H )
    ,
    M_env()
#endif
{
    //_M_desc.add( Feel::feel_options() );

    initMPI( argc, argv, comm );

    doOptions( argc, argv );

#if defined( HAVE_MPI_H )
    char * __env = getenv("DEBUG");
    std::string env_str;
    if ( __env )
        env_str = __env;
    mpi::broadcast( *M_comm, env_str, 0 );
    if ( processId() != 0 )
        {
            setenv( "DEBUG", env_str.c_str(), 1 );
            //Debug() << "DEBUG is set to " << env_str << "\n";
            //std::cout << "DEBUG is set to " << env_str << "\n";
        }
#endif // MPI

    initPETSc();
    initTrilinos();

#if defined(HAVE_TAU)
    TAU_PROFILE("Application", "Application::Application( int, char**, AboutData const&, bool)", TAU_DEFAULT);
    TAU_PROFILE_INIT(argc,argv);
    TAU_PROFILE_SET_NODE(0);
#endif /* HAVE_TAU */
}

#if defined( HAVE_MPI_H )
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
#endif // HAVE_MPI_H
    :
    _M_about( ad ),
    _M_desc( "Allowed options" ),
    _M_vm(),
    _M_to_pass_further()
#if defined( HAVE_MPI_H )
    ,
    M_env()
#endif

{
    //_M_desc.add( Feel::feel_options() ).add( od );
    _M_desc.add( od );


    initMPI( argc, argv, comm );

    doOptions( argc, argv );

#if defined( HAVE_MPI_H )
    char * __env = getenv("DEBUG");
    std::string env_str;
    if ( __env )
        env_str = __env;
    mpi::broadcast( *M_comm, env_str, 0 );
    if ( processId() != 0 )
        {
            setenv( "DEBUG", env_str.c_str(), 1 );
            //Debug() << "DEBUG is set to " << env_str << "\n";
            //std::cout << "DEBUG is set to " << env_str << "\n";
        }
#endif // MPI

    initPETSc();
    initTrilinos();

#if defined(HAVE_TAU)
    TAU_PROFILE("Application", "Application::Application( int, char**, AboutData const&, po::options_description const&, bool)", TAU_DEFAULT);
    TAU_PROFILE_INIT(argc,argv);
    TAU_PROFILE_SET_NODE(0);
#endif /* HAVE_TAU */




}

#if defined( HAVE_MPI_H )
Application::Application( AboutData const& ad,
                          po::options_description const& od,
                          MPI_Comm comm )
#else
Application::Application( AboutData const& ad,
                          po::options_description const& od )
#endif // HAVE_MPI_H
    :
    _M_about( ad ),
    _M_desc( "Allowed options" ),
    _M_vm(),
    _M_to_pass_further()
#if defined( HAVE_MPI_H )
    ,
    M_env()
#endif

{
    //_M_desc.add( Feel::feel_options() ).add( od );
    _M_desc.add( od );

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
    argv[0] = new char[_M_about.appName().size()+1];
    ::strcpy( argv[0], _M_about.appName().c_str() );


    initMPI( argc, argv, comm );

    doOptions( argc, argv );

#if defined( HAVE_MPI_H )
    char * __env = getenv("DEBUG");
    std::string env_str;
    if ( __env )
        env_str = __env;
    mpi::broadcast( *M_comm, env_str, 0 );
    if ( processId() != 0 )
        {
            setenv( "DEBUG", env_str.c_str(), 1 );
            //Debug() << "DEBUG is set to " << env_str << "\n";
            //std::cout << "DEBUG is set to " << env_str << "\n";
        }
#endif // MPI

    initPETSc();
    initTrilinos();

#if defined(HAVE_TAU)
    TAU_PROFILE("Application", "Application::Application( int, char**, AboutData const&, po::options_description const&, bool)", TAU_DEFAULT);
    TAU_PROFILE_INIT(argc,argv);
    TAU_PROFILE_SET_NODE(0);
#endif /* HAVE_TAU */




}

#if defined( HAVE_MPI_H )
Application::Application( AboutData const& ad,
                          MPI_Comm comm )
#else
Application::Application( AboutData const& ad )

#endif // HAVE_MPI_H
    :
    _M_about( ad ),
    _M_desc( "Allowed options" ),
    _M_vm(),
    _M_to_pass_further()
#if defined( HAVE_MPI_H )
    ,
    M_env()
#endif

{
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
    argv[0] = new char[_M_about.appName().size()+1];
    ::strcpy( argv[0], _M_about.appName().c_str() );


    initMPI( argc, argv, comm );
#endif

#if defined( HAVE_MPI_H )
    char * __env = getenv("DEBUG");
    std::string env_str;
    if ( __env )
        env_str = __env;
    mpi::broadcast( *M_comm, env_str, 0 );
    if ( processId() != 0 )
        {
            setenv( "DEBUG", env_str.c_str(), 1 );
            //Debug() << "DEBUG is set to " << env_str << "\n";
            //std::cout << "DEBUG is set to " << env_str << "\n";
        }
#endif // MPI

    initPETSc();
    initTrilinos();

#if defined(HAVE_TAU)
    TAU_PROFILE("Application", "Application::Application( int, char**, AboutData const&, po::options_description const&, bool)", TAU_DEFAULT);
    TAU_PROFILE_INIT(argc,argv);
    TAU_PROFILE_SET_NODE(0);
#endif /* HAVE_TAU */

}
Application::Application( Application const& __app )
    :
    _M_about( __app._M_about ),
    _M_desc( __app._M_desc ),
    _M_vm( __app._M_vm ),
    _M_to_pass_further( __app._M_to_pass_further)
{
}
Application::~Application()
{
#if 0
#if defined ( HAVE_PETSC_H )
    PetscTruth is_petsc_initialized;
    PetscInitialized( &is_petsc_initialized );
    if ( is_petsc_initialized )
    {
#if defined( FEEL_HAVE_SLEPC )
        SlepcFinalize();
#else
        PetscFinalize();
#endif // FEEL_HAVE_SLEPC
    }
#endif // HAVE_PETSC_H

#endif // 0
}
po::options_description
benchmark_options( std::string const& prefix  )
{
    po::options_description benchopt( "Benchmarking options" );
    benchopt.add_options()
        (prefixvm(prefix,"benchmark.nlevels").c_str(), po::value<int>()->default_value(1), "number of mesh levels to benchmark")
        (prefixvm(prefix,"benchmark.hsize").c_str(), po::value<double>()->default_value(0.1), "default mesh size")
        (prefixvm(prefix,"benchmark.refine").c_str(), po::value<double>()->default_value(2), "refine ratio for meshes")
        ;
    // make sense only for global benchmark options
    if ( prefix.empty() )
        benchopt.add_options()
            ("benchmark.only", po::value<std::string>()->default_value(""), "benchmarks to run, empty means all")
            ;
    return benchopt;
}
void
Application::doOptions( int argc, char** argv )
{
    try{
        po::options_description generic( "Generic options" );
        generic.add_options()
            ("authors", "prints the authors list")
            ("copyright", "prints the copyright statement")
            ("help", "prints this help message")
            ("license", "prints the license text")
            ("version,v", "prints the version")
            ("feelinfo", "prints feel libraries information")
            ("verbose,V", "verbose mode")
            ("nochdir", "Don't change repository directory even though it is called")
            ("config-file", po::value<std::string>(), "specify .cfg file")
            ("response-file", po::value<std::string>(), "can be specified with '@name', too")
            ;
        po::options_description debug( "Debugging options" );
        debug.add_options()
            ("debug", po::value<std::string>()->default_value( "" ), "specify a debugging area list");
        _M_desc.add( generic ).add( debug ).add( benchmark_options() );

        this->parseAndStoreOptions( po::command_line_parser(argc, argv), true );
        processGenericOptions();
        /**
         * parse config file if given to command line
         */
        if ( _M_vm.count("config-file") )
        {
            Debug( 1000 ) << "[Application] parsing " << _M_vm["config-file"].as<std::string>() << "\n";
            if ( fs::exists(  _M_vm["config-file"].as<std::string>() ) )
            {

                std::ifstream ifs( _M_vm["config-file"].as<std::string>().c_str() );
                po::store(parse_config_file(ifs, _M_desc, true), _M_vm);
                po::notify(_M_vm);
            }
        }
        std::vector<fs::path> prefixes = boost::assign::list_of( fs::current_path() )
            ( fs::path (Environment::localConfigRepository() ) )
            ( fs::path (Environment::systemConfigRepository().get<0>() ) )
            ( fs::path ("/usr/share/feel/config") )
            ( fs::path ("/usr/local/share/feel/config") )
            ( fs::path ("/opt/local/share/feel/config") );

        BOOST_FOREACH( auto prefix, prefixes )
        {
            std::string config_name = (boost::format( "%1%/%2%.cfg" ) % prefix.string() % this->about().appName()).str();
            Debug( 1000 ) << "[Application] Looking for " << config_name << "\n";
            Debug( 1000 ) << "[Application] Looking for " << config_name << "\n";
            if ( fs::exists( config_name ) )
            {
                Debug( 1000 ) << "[Application] parsing " << config_name << "\n";
                std::ifstream ifs( config_name.c_str() );
                store(parse_config_file(ifs, _M_desc, true), _M_vm);
                break;
            }
            else
            {
                // try with a prefix feel_
                std::string config_name = (boost::format( "%1%/feel_%2%.cfg" ) % prefix.string() % this->about().appName()).str();
                Debug( 1000 ) << "[Application] Looking for " << config_name << "\n";
                if ( fs::exists( config_name ) )
                {
                    Debug( 1000 ) << "[Application] loading configuration file " << config_name << "...\n";
                    std::ifstream ifs( config_name.c_str() );
                    store(parse_config_file(ifs, _M_desc, true), _M_vm);
                    break;
                }
            }
        }
        //po::store(po::parse_command_line(argc, argv, _M_desc), _M_vm);
        po::notify(_M_vm);

    }
    catch( boost::program_options::unknown_option const& e )
    {
        std::cout << "[Application::Application] unknown option:" << e.what() << "\n";
    }
    catch( ... )
    {
        std::cout << "[Application::Application] unknown exception\n";
    }
}
char**
Application::unknownArgv() const
{
    char** argv = new char*[ _M_to_pass_further.size()+1 ];
    argv[0] = new char[ std::strlen(this->about().appName().c_str())+1 ];
    strcpy( argv[0], this->about().appName().c_str() );
    argv[0][std::strlen(this->about().appName().c_str())] = '\0';
    int n_a = 0;
    Debug( 1000 ) << "argv[ " << n_a << " ]=" << argv[0] << "\n";
    ++n_a;
    BOOST_FOREACH( std::string const& s, _M_to_pass_further )
        {
            size_type ssize=s.size();
            Debug( 1000 ) << "new arg " << s << " size = " << ssize << "\n";
            argv[n_a] = new char[ s.size()+1 ];
            strncpy( argv[n_a], s.c_str(), s.size() );
            argv[n_a][s.size()] = '\0';
            ++n_a;

            Debug( 1000 ) << "argv[ " << n_a-1 << " ]=" << argv[n_a-1] << "\n";
        }
    return argv;
}
void
Application::setName1( std::string const& name1 )
{
    _M_name1 = name1;
}

void
Application::setName2( std::string const& name2 )
{
    _M_name2 = name2;
}

void
Application::setH( double h, int precision )
{
    _M_h = std::make_pair( h, precision );
}

void
Application::setDimension( int dim )
{
    _M_dim = dim;
}
void
Application::setLogs()
{
    Debug( 1000 ).detachAll();
    std::ostringstream ostr;
    ostr << this->about().appName() << "-" << nProcess()  << "." << processId();
    Debug( 1000 ).attach( ostr.str() );

    std::ostringstream ostr_assert;
    ostr_assert << this->about().appName() << "_assertions" << "-" << nProcess()  << "." << processId();
    Assert::setLog( ostr_assert.str().c_str() );
}
void
Application::processGenericOptions()
{
    // leave this to subclasses or users
#if 0
    if ( _M_vm.count( "help" ) )
        std::cout << _M_desc << "\n";
#endif


    if ( _M_vm.count("response-file") )
        {
            using namespace std;
            // Load the file and tokenize it
            ifstream ifs( _M_vm["response-file"].as<std::string>().c_str() );
            if (!ifs) {
                cout << "Could not open the response file\n";
                return ;
            }
            // Read the whole file into a string
            stringstream ss;
            ss << ifs.rdbuf();
            // Split the file content
            boost::char_separator<char> sep(" \n\r");
            boost::tokenizer<boost::char_separator<char> > tok(ss.str(), sep);
            vector<string> args;
            copy(tok.begin(), tok.end(), back_inserter(args));

            this->parseAndStoreOptions( po::command_line_parser( args ) );
        }

    if ( _M_vm.count( "feelinfo" ) )
        std::cout << std::setw( 15 ) << std::right << "Feel Version : " << Info::versionString() << "\n"
                  << std::setw( 15 ) << std::right << "Major : " << Info::versionMajor() << "\n"
                  << std::setw( 15 ) << std::right << "Minor : " << Info::versionMinor() << "\n"
                  << std::setw( 15 ) << std::right << "Micro : " << Info::versionMicro() << "\n"
                  << std::setw( 15 ) << std::right << "Revision : " << Info::revision() << "\n"
                  << std::setw( 15 ) << std::right << "BuildId : " << Info::buildId() << "\n"
                  << std::setw( 15 ) << std::right << "Feel Prefix : " << Info::prefix() << "\n"
                  << std::setw( 15 ) << std::right << "Feel DataDir : " << Info::datadir() << "\n";

    if ( _M_vm.count( "verbose" ) ||
         _M_vm.count( "help" ) ||
         _M_vm.count( "version" ) ||
         _M_vm.count( "copyright" ) ||
         _M_vm.count( "license" ) ||
         _M_vm.count( "authors" ) )
        std::cout << _M_about.appName() << ": " << _M_about.shortDescription() <<  "\n";

    if ( _M_vm.count( "version" ) )
        std::cout << " version : " << _M_about.version() << "\n";
    if ( _M_vm.count( "copyright" ) )
        std::cout << " copyright : " << _M_about.copyrightStatement() << "\n";
    if ( _M_vm.count( "license" ) )
        std::cout << " license : " << _M_about.license() << "\n";
    if ( _M_vm.count( "authors" ) )
        {
            std::cout << std::setw( 30 )
                      << "Author Name"
                      << " " << std::setw( 15 )
                      << "Task"
                      << " " << std::setw( 40 )
                      << "Email Address"
                      << "\n";
            std::cout << std::setw( 85+3 ) << std::setfill( '-' ) << "\n" << std::setfill( ' ' );
            std::for_each( _M_about.authors().begin(),
                           _M_about.authors().end(),
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
                           << "\n");
        }

#if 0
    std::cout << "count = " << _M_vm.count( "debug" ) << "\n"
              << "string = " << _M_vm["debug"].as<std::string>() << "\n";
#endif
    if ( _M_vm.count( "debug" ) && !_M_vm["debug"].as<std::string>().empty() )
        DebugStream::showDebugAreas( _M_vm["debug"].as<std::string>() );

    Debug( 1000 ) << "[Application::processGenericOptions] done\n";
}
std::string
Application::rootRepository() const
{
    return Environment::rootRepository();
}

Application&
Application::changeRepository( boost::format fmt )
{
    if ( _M_vm.count( "nochdir" ) )
    {
        this->setLogs();
        return *this;
    }
    Environment::changeRepository( fmt, _M_about.appName() );
    this->setLogs();
    return *this;
}
void
Application::parseAndStoreOptions( po::command_line_parser parser, bool extra_parser )
{
    Debug( 1000 ) << "[Application::Application] parsing options...\n";

    boost::shared_ptr<po::parsed_options> parsed;
    if ( extra_parser )
        {
            parsed = boost::shared_ptr<po::parsed_options>( new po::parsed_options( parser
                                                                                    .options(_M_desc)
                                                                                    .extra_parser(at_option_parser)
                                                                                    .allow_unregistered()
                                                                                    .run() ) );
        }
    else
        {
            parsed = boost::shared_ptr<po::parsed_options>( new po::parsed_options( parser
                                                                                    .options(_M_desc)
                                                                                    .allow_unregistered()
                                                                                    .run() ) );
        }

    Debug( 1000 ) << "[Application::Application] parsing options done\n";

    _M_to_pass_further = po::collect_unrecognized( parsed->options, po::include_positional );
    Debug( 1000 ) << "[Application::Application] number of unrecognized options: " << (_M_to_pass_further.size()) << "\n";

    BOOST_FOREACH( std::string const& s, _M_to_pass_further )
        {
            Debug( 1000 ) << "[Application::Application] option: " << s << "\n";
        }
    std::vector<po::basic_option<char> >::iterator it = parsed->options.begin();
    std::vector<po::basic_option<char> >::iterator en  = parsed->options.end();
    for ( ; it != en ; ++it )
        if ( it->unregistered )
            {
                Debug( 1000 ) << "[Application::Application] remove from vector " << it->string_key << "\n";
                parsed->options.erase( it );
            }

    po::store(*parsed, _M_vm );
}

void
Application::add( Simget* simget )
{
    M_simgets.push_back( simget );
}

void
Application::run()
{
    std::string runonly = _M_vm["benchmark.only"].as<std::string>();
    for( auto i = M_simgets.begin(), end = M_simgets.end(); i != end; ++i )
    {
        if ( ( runonly.empty() == false  ) &&
             runonly.find( i->name() ) == std::string::npos )
            continue;


        std::string s1 = prefixvm(i->name(),"benchmark.nlevels");
        int nlevels = _M_vm.count( s1 )?_M_vm[s1].as<int>():_M_vm["benchmark.nlevels"].as<int>();
        std::string s2 = prefixvm(i->name(),"benchmark.hsize");
        double hsize = _M_vm.count( s2 )?_M_vm[s2].as<double>():_M_vm["benchmark.hsize"].as<double>();
        std::string s3 = prefixvm(i->name(),"benchmark.refine");
        double refine = _M_vm.count( s3 )?_M_vm[s3].as<double>():_M_vm["benchmark.refine"].as<double>();

        for(int l = 0; l < nlevels; ++l )
        {
            double meshSize= hsize/std::pow(refine,l);
            i->setMeshSize( meshSize );
            i->run();
            M_stats[i->name()].push_back( i->stats() );
        }

    }
}

void
Application::run( const double* X, unsigned long P, double* Y, unsigned long N )
{
    auto it = M_simgets.begin(), en = M_simgets.end();
    int i = 0;
    for( ; it != en; ++i, ++it )
    {
        it->run( &X[i*P], P, &Y[i*N], N );
    }
}

void
printErrors( std::ostream& out, std::vector<ptree::ptree> const& stats, std::string const& key )
{
    out << std::setw(10) << std::right << "levels"
        << std::setw(10) << std::right << "h";
    BOOST_FOREACH(auto v, stats.front().get_child(key))
        {
            out << std::setw(15) << std::right << v.first
                << std::setw(15) << std::right << "ROC";
        }
    out << "\n";
    int l=1;
    for( auto it = stats.begin(), en =stats.end(); it!=en; ++it,++l )
    {
        //std::for_each( it->begin(),it->end(), []( std::pair<std::string,boost::any> const& o ) { std::cout << o.first << "\n"; } );
        //std::map<std::string,boost::any> data = *it;
        //std::map<std::string,boost::any> datap;
        double rocu = 1, rocp=1;
        double h  = it->get<double>("h");
        double hp = h;
        out << std::right << std::setw(10) << l
            << std::right << std::setw(10) << std::fixed  << std::setprecision( 4 ) << h;
        if ( l > 1 )
            hp = boost::prior(it)->get<double>("h");
        BOOST_FOREACH(auto v, it->get_child(key))
        {
            double u  = it->get<double>( key+"."+v.first );
            double up = u;
            double roc = 1;
            if ( l > 1 )
            {
                up  = boost::prior(it)->get<double>( key+"."+v.first );
                roc = std::log10( up/u )/std::log10( hp/h );
            }
            out << std::right << std::setw(15) << std::scientific << std::setprecision( 2 ) << u
                << std::right << std::setw(15) << std::fixed << std::setprecision( 2 ) << roc;
        }
        out << "\n";
    }
}
void
printNumbers( std::ostream& out, std::vector<ptree::ptree> const& stats, std::string const& key )
{
    out << std::setw(10) << std::right << "levels"
        << std::setw(10) << std::right << "h";
    BOOST_FOREACH(auto v, stats.front().get_child(key))
        {
            out << std::setw(10) << std::right << v.first;
        }
    out << "\n";
    int l=1;
    for( auto it = stats.begin(), en =stats.end(); it!=en; ++it,++l )
    {
        double h  = it->get<double>("h");
        out << std::right << std::setw(10) << l
            << std::right << std::setw(10) << std::fixed  << std::setprecision( 4 ) << h;
        BOOST_FOREACH(auto v, it->get_child(key))
        {
            size_type u  = it->get<size_type>( key+"."+v.first );
            out << std::right << std::setw(10)  << u;
        }
        out << "\n";
    }
}

void
printTime( std::ostream& out, std::vector<ptree::ptree> const& stats, std::string const& key )
{
    out << std::setw(10) << std::right << "levels"
        << std::setw(10) << std::right << "h";
    BOOST_FOREACH(auto v, stats.front().get_child(key))
        {
            out << std::setw(15) << std::right << v.first;
        }
    out << "\n";
    int l=1;
    for( auto it = stats.begin(), en =stats.end(); it!=en; ++it,++l )
    {
        double h  = it->get<double>("h");
        out << std::right << std::setw(10) << l
            << std::right << std::setw(10) << std::fixed  << std::setprecision( 4 ) << h;
        BOOST_FOREACH(auto v, it->get_child(key))
        {
            double u  = it->get<double>( key+"."+v.first );
            out << std::right << std::setw(15) << std::scientific << std::setprecision( 2 ) << u;
        }
        out << "\n";
    }
}

void
Application::printStats( std::ostream& out, std::vector<std::string> const& keys ) const
{
    std::string runonly = _M_vm["benchmark.only"].as<std::string>();
    for( auto i = M_simgets.begin(), end = M_simgets.end(); i != end; ++i )
    {
        if ( ( runonly.empty() == false  ) &&
             runonly.find( i->name() ) == std::string::npos )
            continue;
        std::cout << "================================================================================\n";
        std::cout << "Simulation " << i->name() << "\n";
        BOOST_FOREACH( auto key, keys )
        {
            std::cout << "------------------------------------------------------------\n";
            std::cout << "Key: " << key << "\n";
            if ( key.find("e.") != std::string::npos )
            {
                if ( M_stats.find(i->name()) != M_stats.end() )
                    printErrors( out, M_stats.find(i->name())->second, key );
            }
            if ( key.find("n.") != std::string::npos )
            {
                if ( M_stats.find(i->name()) != M_stats.end() )
                    printNumbers( out, M_stats.find(i->name())->second, key );
            }
            if ( key.find("t.") != std::string::npos )
            {
                if ( M_stats.find(i->name()) != M_stats.end() )
                    printTime( out, M_stats.find(i->name())->second, key );
            }
        }
    }
}
}
