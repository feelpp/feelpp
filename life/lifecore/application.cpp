/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
  Date: 2005-03-17

  Copyright (C) 2007,2008,2009 Université de Grenoble 1
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

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifecore/application.hpp>

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

#include <life/lifecore/mpicompat.hpp>


namespace Life
{
namespace fs = boost::filesystem;


int Application::_S_n_process = 1;
int Application::_S_process_id = 0;


LIFE_NO_EXPORT
std::pair<std::string, std::string>
at_option_parser(std::string const&s)
{
    //std::cout << "string =" << s << "\n";
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
#if defined( HAVE_SLEPC )
            int ierr = SlepcInitialize(&__argc,&__argv, PETSC_NULL, PETSC_NULL );
#else
            int ierr = PetscInitialize( &__argc, &__argv, PETSC_NULL, PETSC_NULL );
#endif
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
#if defined( HAVE_MPI_H )
    int is_mpi_initialized;
    MPI_Initialized (&is_mpi_initialized);
    //std::cout << "is_mpi_initialized = " << is_mpi_initialized << "\n";
    if (!is_mpi_initialized)
    {
        M_env = boost::shared_ptr<mpi::environment>( new mpi::environment( argc, argv ) );
        ///std::cout << "processor name = " << M_env->processor_name() << "\n";
#if 0
        //int __argc = this->unknownArgc();
        //char** __argv = this->unknownArgv();
        MPI_Init (&argc, &argv);
#endif
        _S_is_mpi_initialized = true;
    }
    MPI_Comm_dup ( comm, &COMM_WORLD);
    //MPI_Comm_dup ( comm, (MPI_Comm*)&S_world );
#if 0
    MPI_Comm_rank (COMM_WORLD, &_S_process_id);
    MPI_Comm_size (COMM_WORLD, &_S_n_process);
#else
    _S_process_id = S_world.rank();
    _S_n_process = S_world.size();
#endif
#endif // HAVE_MPI_H

}

#if defined( HAVE_MPI_H )
MPI_Comm Application::COMM_WORLD = MPI_COMM_NULL;

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
    //_M_desc.add( Life::life_options() );

    initMPI( argc, argv, comm );

    doOptions( argc, argv );

#if defined( HAVE_MPI_H )
    char * __env = getenv("DEBUG");
    std::string env_str;
    if ( __env )
        env_str = __env;
    mpi::broadcast( S_world, env_str, 0 );
    if ( _S_process_id != 0 )
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
    //_M_desc.add( Life::life_options() ).add( od );
    _M_desc.add( od );


    initMPI( argc, argv, comm );

    doOptions( argc, argv );

#if defined( HAVE_MPI_H )
    char * __env = getenv("DEBUG");
    std::string env_str;
    if ( __env )
        env_str = __env;
    mpi::broadcast( S_world, env_str, 0 );
    if ( _S_process_id != 0 )
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
    mpi::broadcast( S_world, env_str, 0 );
    if ( _S_process_id != 0 )
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
#if defined ( HAVE_PETSC_H )
#if defined( HAVE_SLEPC )
    SlepcFinalize();
#else
    PetscFinalize();
#endif
#endif
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
            ("lifeinfo", "prints life libraries information")
            ("verbose,V", "verbose mode")
            ("response-file", po::value<std::string>(), "can be specified with '@name', too")
            ;
        po::options_description debug( "Debugging options" );
        debug.add_options()
            ("debug", po::value<std::string>()->default_value( "" ), "specify a debugging area list");
        _M_desc.add( generic ).add( debug );


        this->parseAndStoreOptions( po::command_line_parser(argc, argv), true );

        std::string config_name = (boost::format( "%1%.cfg" ) % this->about().appName()).str();
        Debug( 1000 ) << "[Application] Looking for " << config_name << "\n";
        if ( fs::exists( config_name ) )
        {

            std::ifstream ifs( config_name.c_str() );
            store(parse_config_file(ifs, _M_desc), _M_vm);
        }
        else
        {
            // try with a prefix life_
            std::string config_name = (boost::format( "life_%1%.cfg" ) % this->about().appName()).str();
            Debug( 1000 ) << "[Application] Looking for " << config_name << "\n";

            if ( fs::exists( config_name ) )
            {
                std::ifstream ifs( config_name.c_str() );
                store(parse_config_file(ifs, _M_desc), _M_vm);
            }
        }

        //po::store(po::parse_command_line(argc, argv, _M_desc), _M_vm);
        po::notify(_M_vm);

        processGenericOptions();
    }
    catch( boost::program_options::unknown_option const& e )
        {
            std::cout << "[Application::Application] unknown option\n";
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
    Log().detachAll();
    std::ostringstream ostr;
    ostr << this->about().appName() << "-" << _S_n_process  << "." << _S_process_id;
    Log().attach( ostr.str() );

    std::ostringstream ostr_assert;
    ostr_assert << this->about().appName() << "_assertions" << "-" << _S_n_process  << "." << _S_process_id;
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
            ifstream ifs( _M_vm["response-file"].as<string>().c_str() );
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

    if ( _M_vm.count( "lifeinfo" ) )
        std::cout << std::setw( 15 ) << std::right << "Life Version : " << Info::versionString() << "\n"
                  << std::setw( 15 ) << std::right << "Major : " << Info::versionMajor() << "\n"
                  << std::setw( 15 ) << std::right << "Minor : " << Info::versionMinor() << "\n"
                  << std::setw( 15 ) << std::right << "Micro : " << Info::versionMicro() << "\n"
                  << std::setw( 15 ) << std::right << "Revision : " << Info::revision() << "\n"
                  << std::setw( 15 ) << std::right << "BuildId : " << Info::buildId() << "\n"
                  << std::setw( 15 ) << std::right << "Life Prefix : " << Info::prefix() << "\n"
                  << std::setw( 15 ) << std::right << "Life DataDir : " << Info::datadir() << "\n";

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
    std::string env;
    if ( ::getenv( "LIFE_REPOSITORY" ) )
        {
            env = ::getenv( "LIFE_REPOSITORY" );
        }
    else
        {
            // by default create $HOME/life
            env = ::getenv( "HOME" );
            env += "/life";
        }
    return env;
}

Application&
Application::changeRepository( boost::format fmt )
{
    fs::path rep_path;

    rep_path = rootRepository();
    if ( !fs::exists( rep_path ) )
        fs::create_directory( rep_path );

    typedef std::vector< std::string > split_vector_type;

    split_vector_type dirs; // #2: Search for tokens
    std::string fmtstr = fmt.str();
    boost::split( dirs, fmtstr, boost::is_any_of("/") );

    BOOST_FOREACH( std::string const& dir, dirs )
        {
            //Debug( 1000 ) << "[Application::Application] option: " << s << "\n";
            rep_path = rep_path / dir;
            if (!fs::exists( rep_path ) )
                fs::create_directory( rep_path );
        }

    ::chdir( rep_path.string().c_str() );
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

                                                                                    // new at least boost 1.33.1 for this
#if BOOST_VERSION >= 103301

                                                                                    .allow_unregistered()

#endif
                                                                                    .run() ) );
        }
    else
        {
            parsed = boost::shared_ptr<po::parsed_options>( new po::parsed_options( parser
                                                                                    .options(_M_desc)

                                                                                    // new at least boost 1.33.1 for this
#if BOOST_VERSION >= 103301

                                                                                    .allow_unregistered()

#endif
                                                                                    .run() ) );
        }

    Debug( 1000 ) << "[Application::Application] parsing options done\n";

#if BOOST_VERSION >= 103301

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
#endif

    po::store(*parsed, _M_vm );

}

mpi::communicator Application::S_world;
bool Application::_S_is_mpi_initialized = false;
}
