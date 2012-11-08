/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-04-14
 */

#include <boost/program_options.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

#include <gflags/gflags.h>

#include <feel/feelconfig.h>
#include <feel/feelcore/feel.hpp>


#include <feel/feelcore/environment.hpp>

#include <feel/feelcore/feelpetsc.hpp>
#include <feel/options.hpp>


namespace google
{
namespace glog_internal_namespace_
{
bool IsGoogleLoggingInitialized();
}
}
namespace Feel
{

namespace detail
{
FEELPP_NO_EXPORT
std::pair<std::string, std::string>
at_option_parser_2( std::string const&s )
{
    if ( '@' == s[0] )
        return std::make_pair( std::string( "response-file" ), s.substr( 1 ) );

    else
        return std::pair<std::string, std::string>();
}

/**
   \fn freeargv -- free an argument vector

   void freeargv (char** vector)

   Free an argument vector that was built using dupargv.  Simply scans
   through the vector, freeing the memory for each argument until the
   terminating NULL is found, and then frees the vector itself.

*/

void freeargv (char** vector)
{
    char **scan;

    if (vector != NULL)
    {
        for (scan = vector; *scan != NULL; scan++)
        {
            free (*scan);
        }
        free (vector);
    }
}


/**
   \fn dupargv -- duplicate an argument vector

   char **dupargv (char** vector)

   Duplicate an argument vector.  Simply scans through the
   vector, duplicating each argument until the
   terminating NULL is found.

   \return a pointer to the argument vector if
   successful. Returns NULL if there is insufficient memory to
   complete building the argument vector.
*/
char **
dupargv (char** argv)
{
  int argc;
  char **copy;

  if (argv == NULL)
    return NULL;

  /* the vector */
  for (argc = 0; argv[argc] != NULL; argc++);
  copy = (char **) malloc ((argc + 1) * sizeof (char *));
  if (copy == NULL)
      return NULL;

  /* the strings */
  for (argc = 0; argv[argc] != NULL; argc++)
  {
      int len = strlen (argv[argc]);
      copy[argc] = (char*)malloc (sizeof (char *) * (len + 1));
      if (copy[argc] == NULL)
      {
          freeargv (copy);
          return NULL;
      }
      strcpy (copy[argc], argv[argc]);
  }
  copy[argc] = NULL;
  return copy;
}
void
Environment::processGenericOptions()
{
//     // leave this to subclasses or users
// #if 0
//     if ( S_vm.count( "help" ) )
//         std::cout << S_desc << "\n";

// #endif


//     if ( S_vm.count( "response-file" ) )
//     {
//         using namespace std;
//         // Load the file and tokenize it
//         ifstream ifs( S_vm["response-file"].as<std::string>().c_str() );

//         if ( !ifs )
//         {
//             cout << "Could not open the response file\n";
//             return ;
//         }

//         // Read the whole file into a string
//         stringstream ss;
//         ss << ifs.rdbuf();
//         // Split the file content
//         boost::char_separator<char> sep( " \n\r" );
//         boost::tokenizer<boost::char_separator<char> > tok( ss.str(), sep );
//         vector<string> args;
//         copy( tok.begin(), tok.end(), back_inserter( args ) );

//         parseAndStoreOptions( po::command_line_parser( args ) );
//     }

    //if ( worldComm().rank() == 0 )
    {

        if ( S_vm.count( "feelinfo" ) )
            std::cout << std::setw( 15 ) << std::right << "Feel Version : " << Info::versionString() << "\n"
                      << std::setw( 15 ) << std::right << "Major : " << Info::versionMajor() << "\n"
                      << std::setw( 15 ) << std::right << "Minor : " << Info::versionMinor() << "\n"
                      << std::setw( 15 ) << std::right << "Micro : " << Info::versionMicro() << "\n"
                      << std::setw( 15 ) << std::right << "Revision : " << Info::revision() << "\n"
                      << std::setw( 15 ) << std::right << "BuildId : " << Info::buildId() << "\n"
                      << std::setw( 15 ) << std::right << "Feel Prefix : " << Info::prefix() << "\n"
                      << std::setw( 15 ) << std::right << "Feel DataDir : " << Info::datadir() << "\n";

        if ( S_vm.count( "verbose" ) ||
             S_vm.count( "help" ) ||
             S_vm.count( "version" ) ||
             S_vm.count( "copyright" ) ||
             S_vm.count( "license" ) ||
             S_vm.count( "authors" ) )
        {
            std::cout << S_about.appName() << ": " << S_about.shortDescription() <<  "\n";
        }

        if ( S_vm.count( "version" ) )
        {
            std::cout << " version : " << S_about.version() << "\n";
        }

        if ( S_vm.count( "copyright" ) )
        {
            std::cout << " copyright : " << S_about.copyrightStatement() << "\n";
        }

        if ( S_vm.count( "license" ) )
        {
            std::cout << " license : " << S_about.license() << "\n";
        }

        if ( S_vm.count( "authors" ) )
        {
#if 0
            std::cout << std::setw( 30 )
                      << "Author Name"
                      << " " << std::setw( 15 )
                      << "Task"
                      << " " << std::setw( 40 )
                      << "Email Address"
                      << "\n";
            std::cout << std::setw( 85+3 ) << std::setfill( '-' ) << "\n" << std::setfill( ' ' );
            std::for_each( S_about.authors().begin(),
                           S_about.authors().end(),
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
#endif
        }

        if ( S_vm.count( "help" ) )
        {
            std::cout << optionsDescription() << "\n";
        }


    }
    if ( S_vm.count( "verbose" ) ||
         S_vm.count( "help" ) ||
         S_vm.count( "version" ) ||
         S_vm.count( "copyright" ) ||
         S_vm.count( "license" ) ||
         S_vm.count( "authors" ) )
    {
        if ( Environment::initialized() )
        {
            worldComm().barrier();
            MPI_Finalize();
        }
        exit(0);
    }
#if 0
    std::cout << "count = " << S_vm.count( "debug" ) << "\n"
              << "string = " << S_vm["debug"].as<std::string>() << "\n";
#endif

    if ( S_vm.count( "debug" ) && !S_vm["debug"].as<std::string>().empty() )
        DebugStream::showDebugAreas( S_vm["debug"].as<std::string>() );

    VLOG(2) << "[processGenericOptions] done\n";
}

void
Environment::parseAndStoreOptions( po::command_line_parser parser, bool extra_parser )
{
    VLOG(2) << " parsing options...\n";

    boost::shared_ptr<po::parsed_options> parsed;

    if ( extra_parser )
    {
        parsed = boost::shared_ptr<po::parsed_options>( new po::parsed_options( parser
                                                                                .options( *S_desc )
                                                                                .extra_parser( at_option_parser_2 )
                                                                                .allow_unregistered()
                                                                                .run() ) );
    }

    else
    {
        parsed = boost::shared_ptr<po::parsed_options>( new po::parsed_options( parser
                                                                                .options( *S_desc )
                                                                                .allow_unregistered()
                                                                                .run() ) );
    }

    VLOG(2) << "[parseAndStoreOptions] parsing options done\n";

    S_to_pass_further = po::collect_unrecognized( parsed->options, po::include_positional );
    VLOG(2)<< " number of unrecognized options: " << ( S_to_pass_further.size() ) << "\n";

    BOOST_FOREACH( std::string const& s, S_to_pass_further )
    {
        VLOG(2)<< " option: " << s << "\n";
    }
    std::vector<po::basic_option<char> >::iterator it = parsed->options.begin();
    std::vector<po::basic_option<char> >::iterator en  = parsed->options.end();

    for ( ; it != en ; ++it )
        if ( it->unregistered )
        {
            VLOG(2)<< " remove from vector " << it->string_key << "\n";
            parsed->options.erase( it );
        }

    po::store( *parsed, S_vm );
}


void
Environment::doOptions( int argc, char** argv, po::options_description const& desc, std::string const& appName )
{
    //std::locale::global(std::locale(""));
    try
    {
        parseAndStoreOptions( po::command_line_parser( argc, argv ), true );
        processGenericOptions();

        /**
         * parse config file if given to command line
         */
        if ( S_vm.count( "config-file" ) )
        {
            VLOG(2)<< " parsing " << S_vm["config-file"].as<std::string>() << "\n";

            if ( fs::exists(  S_vm["config-file"].as<std::string>() ) )
            {

                std::ifstream ifs( S_vm["config-file"].as<std::string>().c_str() );
                po::store( parse_config_file( ifs, desc, true ), S_vm );
                po::notify( S_vm );
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
            std::string config_name = ( boost::format( "%1%/%2%.cfg" ) % prefix.string() % appName ).str();
            VLOG(2)<< " Looking for " << config_name << "\n";
            VLOG(2)<< " Looking for " << config_name << "\n";

            if ( fs::exists( config_name ) )
            {
                VLOG(2)<< " parsing " << config_name << "\n";
                std::ifstream ifs( config_name.c_str() );
                store( parse_config_file( ifs, desc, true ), S_vm );
                break;
            }

            else
            {
                // try with a prefix feel_
                std::string config_name = ( boost::format( "%1%/feel_%2%.cfg" ) % prefix.string() % appName ).str();
                VLOG(2)<< " Looking for " << config_name << "\n";

                if ( fs::exists( config_name ) )
                {
                    VLOG(2)<< " loading configuration file " << config_name << "...\n";
                    std::ifstream ifs( config_name.c_str() );
                    store( parse_config_file( ifs, desc, true ), S_vm );
                    break;
                }
            }
        }
        //po::store(po::parse_command_line(argc, argv, desc), S_vm);
        po::notify( S_vm );

    }

    // catches program_options exceptions
    catch ( boost::program_options::multiple_occurrences const& e )
    {
        LOG(WARNING) << "Command line or config file option parsing error: " << e.what() << "\n"
                     << "  o faulty option: " << e.get_option_name() << "\n"
                     << "Warning: the .cfg file or some options may not have been read properly\n";
    }

    catch ( boost::program_options::ambiguous_option const& e )
    {
        LOG(WARNING) << "Command line or config file option parsing error: " << e.what() << "\n"
                     << "  o faulty option: " << e.get_option_name() << "\n"
                     << "  o possible alternatives: " ;
        std::for_each( e.alternatives().begin(), e.alternatives().end(), []( std::string const& s )
                       {
                           LOG(WARNING) << s << " ";
                       } );
        LOG(WARNING) << "\n"
                     << "Warning: the .cfg file or some options may not have been read properly\n";
    }

    // catches program_options exceptions
    catch ( std::exception& e )
    {
        LOG(WARNING) << "Application option parsing: unknown option:" << e.what() << " (the .cfg file or some options may not have been read properly)\n";
    }

    catch ( ... )
    {
        LOG(WARNING) << "Application option parsing: unknown exception triggered  (the .cfg file or some options may not have been read properly)\n";
    }
}

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
    LOG(INFO) << "[Feel++] TBB running with " << n << " threads\n";
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
    // if scratsch dir not defined, define it
    const char* env;
    env = getenv("FEELPP_SCRATCHDIR");
    if (env != NULL && env[0] != '\0')
    {
        env = getenv("SCRATCHDIR");
        if (env != NULL && env[0] != '\0')
        {
            std::string value = (boost::format("%1%/feelpp/") % env).str();
            ::setenv("FEELPP_SCRATCHDIR",value.c_str(), 0 );
        }
        else
        {
            std::string value = (boost::format("/tmp/feelpp/") % env).str();
            ::setenv("FEELPP_SCRATCHDIR",value.c_str(), 0 );
        }
    }
    env = getenv("FEELPP_SCRATCHDIR");
    S_scratchdir = fs::path( env );


    google::AllowCommandLineReparsing();
    google::ParseCommandLineFlags(&argc, &argv, false);

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
    //LOG(INFO) << "[Feel++] TBB running with " << n << " threads\n";
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
    CHECK( S_worldcomm ) << "Feel++ Environment: creang worldcomm failed!";
}

void
Environment::init( int argc, char** argv, po::options_description const& desc, AboutData const& about )
{
    // duplicate argv before passing to gflags because gflags is going to
    // rearrange them and it screws badly the flags for PETSc/SLEPc
    char** envargv = dupargv( argv );

    google::AllowCommandLineReparsing();
    google::ParseCommandLineFlags(&argc, &argv, false);

#if 0
    std::cout << "argc=" << argc << "\n";
    for(int i = 0; i < argc; ++i )
    {
        std::cout << "argv[" << i << "]=" << argv[i] << "\n";
    }
#endif
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
    //LOG(INFO) << "[Feel++] TBB running with " << n << " threads\n";
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
        int ierr = SlepcInitialize( &argc,&envargv, PETSC_NULL, PETSC_NULL );
#else
        int ierr = PetscInitialize( &argc, &envargv, PETSC_NULL, PETSC_NULL );
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
    CHECK( S_worldcomm ) << "Feel++ Environment: creang worldcomm failed!";

    S_about = about;
    doOptions( argc, envargv, *S_desc, about.appName() );

    freeargv( envargv );

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
Environment::changeRepositoryImpl( boost::format fmt, std::string const& logfilename, bool add_subdir_np )
{
    if ( Environment::vm().count( "nochdir" ) )
        return;

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
        //VLOG(2)<< " option: " << s << "\n";
        rep_path = rep_path / dir;

        if ( !fs::exists( rep_path ) )
            fs::create_directory( rep_path );
    }
    if ( add_subdir_np )
    {
        rep_path = rep_path / (boost::format( "np_%1%" ) % Environment::numberOfProcessors() ).str();
        if ( !fs::exists( rep_path ) )
            fs::create_directory( rep_path );
        LOG(INFO) << "rep_path=" << rep_path << "\n";
    }
    ::chdir( rep_path.string().c_str() );

    setLogs( logfilename );
}

#if 0
po::variables_map
Environment::vm( po::options_description const& desc )
{
    po::variables_map vm;
    po::store( po::parse_command_line( 0, ( char** )0, desc ), vm );
    po::notify( vm );

    return vm;
}
#endif

void
Environment::setLogs( std::string const& prefix )
{

    mpi::communicator world;
#if 0
    LOG(INFO).detachAll();
    std::ostringstream ostr;
    ostr << prefix << "-" << world.size()  << "." << world.rank();
    LOG(INFO).attach( ostr.str() );
#endif

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

AboutData Environment::S_about;
po::variables_map Environment::S_vm;
boost::shared_ptr<po::options_description> Environment::S_desc;
std::vector<std::string> Environment::S_to_pass_further;

boost::signals2::signal<void ()> Environment::S_deleteObservers;

boost::shared_ptr<WorldComm> Environment::S_worldcomm;

fs::path Environment::S_scratchdir;

} // detail

}
