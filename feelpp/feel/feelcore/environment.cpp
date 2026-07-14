/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
    SPDX-FileContributor: Vincent Chabannes <vincent.chabannes@feelpp.org>

    SPDX-FileCopyrightText: 2007-2011 Joseph Fourier University
    SPDX-FileCopyrightText: 2011-2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

#include <cstdlib>
#include <pwd.h>
#include <utility>
#ifdef __cplusplus
extern "C"
{
#endif
#include <sys/stat.h>
#ifdef __cplusplus
}
#endif
#if defined(FEELPP_HAS_PYTHON)
#include <feel/feelpython/pybind11/pybind11.h>
#include <feel/feelpython/pybind11/embed.h>
#endif

#include <boost/program_options.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/smart_ptr/make_shared.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/local_time_adjustor.hpp>
#include <boost/date_time/c_local_time_adjustor.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <range/v3/view/take_exactly.hpp>

#include <feel/feelcore/mongocxx.hpp>

#include <fmt/chrono.h>
//#include <gflags/gflags.h>

#if defined(FEELPP_HAS_SPDLOG)
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/null_sink.h>
#include <feel/feelcore/logger.hpp>
#else
#include <glog/logging.h>
#endif

#include <feel/feelinfo.h>
#include <feel/feelconfig.h>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/table.hpp>

#if defined ( FEELPP_HAS_PETSC_H )
#include <petscsys.h>
#endif
#if defined( FEELPP_HAS_GMSH_H )
#if defined( FEELPP_HAS_GMSH_API )
#include <gmsh.h>
#else
#include <Gmsh.h>
#endif
#endif

#include <feel/feelcore/environment.hpp>

#if FEELPP_HAS_PETSC
#include <feel/feelcore/feelpetsc.hpp>
#endif
#include <feel/feelcore/timertable.hpp>
#include <feel/feelcore/utility.hpp>
#include <feel/feeltiming/tic.hpp>
#include <feel/options.hpp>
#include <feel/feelcore/remotedata.hpp>

#define stringize2(x) #x
#define stringize(x) stringize2(x)



namespace GiNaC
{
extern void cleanup_ex( bool verbose );
}
namespace detail
{
void DebugWait(int rank)
{
    char    a;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    std::cout << "PID " << getpid() << " on " << hostname <<  " ready for attach" << std::endl;
    if(rank == 0) {
        std::cin >> a;
        std::cout << rank << ": Starting now\n";
    }

    MPI_Bcast(&a, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
    std::cout << rank << ": Starting now\n";
}
class Env
{
public:
    static std::string getUserName()
    {
        struct passwd *pw;
        uid_t uid;
        int c;

        uid = geteuid ();
        pw = getpwuid ( uid );

        if ( pw )
        {
            return std::string( pw->pw_name );
        }

        return std::string( "" );
    }
};
}
namespace google
{
namespace glog_internal_namespace_
{
bool IsGoogleLoggingInitialized();
}
}
namespace Feel
{
namespace pt =  boost::property_tree;
#if defined(FEELPP_HAS_PYTHON)
namespace py = pybind11;
using namespace py::literals;
#endif
//namespace detail
//{
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

void freeargv ( char** vector )
{
    char **scan;

    if ( vector != NULL )
    {
        for ( scan = vector; *scan != NULL; scan++ )
        {
            free ( *scan );
        }

        free ( vector );
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
dupargv ( char** argv )
{
    int argc;
    char **copy;

    if ( argv == NULL )
        return NULL;

    /* the vector */
    for ( argc = 0; argv[argc] != NULL; argc++ );

    copy = ( char ** ) malloc ( ( argc + 1 ) * sizeof ( char * ) );

    if ( copy == NULL )
        return NULL;

    /* the strings */
    for ( argc = 0; argv[argc] != NULL; argc++ )
    {
        int len = strlen ( argv[argc] );
        copy[argc] = ( char* )malloc ( sizeof ( char * ) * ( len + 1 ) );

        if ( copy[argc] == NULL )
        {
            freeargv ( copy );
            return NULL;
        }

        strcpy ( copy[argc], argv[argc] );
    }

    copy[argc] = NULL;
    return copy;
}

AboutData makeAboutDefault( std::string name )
{
    AboutData about( name,
                     name,
                     "0.1",
                     name,
                     AboutData::License_GPL,
                     "Copyright (c) 2012-2020 Feel++ Consortium" );

    about.addAuthor( "Feel++ Consortium",
                     "",
                     "feelpp-devel@feelpp.org", "" );
    return about;
}

fs::path scratchdir()
{
    const char* env;
    // if scratsch dir not defined, define it
    env = getenv( "FEELPP_SCRATCHDIR" );

    if ( env == NULL || env[0] == '\0' )
    {
        env = getenv( "SCRATCHDIR" );

        if ( env != NULL && env[0] != '\0' )
        {
            std::string value = ( boost::format( "%1%/%2%/feelpp/" ) % env % ::detail::Env::getUserName() ).str();
            setenv( "FEELPP_SCRATCHDIR", ( boost::format( "%1%/%2%/feelpp/" ) % env % ::detail::Env::getUserName() ).str().c_str(),0 );
        }

        else
        {
            env = getenv( "SCRATCH" );

            if ( env != NULL && env[0] != '\0' )
            {
                std::string value = ( boost::format( "%1%/%2%/feelpp/" ) % env % ::detail::Env::getUserName() ).str();
                setenv( "FEELPP_SCRATCHDIR", ( boost::format( "%1%/%2%/feelpp/" ) % env % ::detail::Env::getUserName() ).str().c_str(),0 );
            }

            else
            {
                std::string value = ( boost::format( "/tmp/%1%/feelpp/" ) % ::detail::Env::getUserName() ).str();
                setenv( "FEELPP_SCRATCHDIR", value.c_str(),0 );
            }
        }
    }

    env = getenv( "FEELPP_SCRATCHDIR" );

    if ( env != NULL && env[0] != '\0' )
    {
        return fs::path( env );
    }

    std::string value = ( boost::format( "/tmp/%1%/feelpp/" ) % ::detail::Env::getUserName() ).str();
    return fs::path( value );
}



//! Default constructor.
Environment::Environment()
    :
#if BOOST_VERSION >= 105500
    Environment( 0, nullptr, mpi::threading::single, feel_nooptions(), feel_options(), makeAboutDefault("feelpp"), globalRepository(makeAboutDefault("feelpp").appName()) )
#else
    Environment( 0, nullptr, feel_nooptions(), feel_options(), makeAboutDefault("feelpp"), makeAboutDefault("feelpp").appName() )
#endif

{
}



//! Constructor
Environment::Environment( int& argc, char**& argv )
    :
    Environment( argc, argv, mpi::threading::single, feel_nooptions(), feel_options(),
                 makeAboutDefault(argv[0]), globalRepository( makeAboutDefault(argv[0]).appName() ) )
{
}



#if defined(FEELPP_HAS_PYTHON)
struct PythonArgs
{
#if defined(FEELPP_HAS_BOOST_PYTHON)
    PythonArgs( boost::python::list arg )
        {
            if ( argv == nullptr )
            {
                /* Convert python options into argc/argv format */

                argc = boost::python::len( arg );

                argv =new char* [argc+1];
                boost::python::stl_input_iterator<std::string> begin( arg ), end;
                int i=0;

                while ( begin != end )
                {
                    //std::cout << *begin << std::endl ;
                    argv[i] =strdup( ( *begin ).c_str() );
                    begin++;
                    i++;
                }

                argv[argc]=nullptr;
            }
        }
#endif
    PythonArgs( pybind11::list arg )
        {
            argc = 0;
            for (auto item : arg)
            {
                ++argc;
            }
            argv =new char* [argc+1];
            int i = 0;
            for (auto item : arg)
            {
                argv[i++] = strdup( std::string(pybind11::str(item) ).c_str() );
                argv[argc]=nullptr;
            }

        }
    static int argc;
    static char** argv;
};
int PythonArgs::argc = 1;
char** PythonArgs::argv = nullptr;

// Constructor
Environment::Environment( pybind11::list arg )
    :
    Environment( PythonArgs(arg).argc, PythonArgs::argv, mpi::threading::single, feel_nooptions(), feel_options(), makeAboutDefault(PythonArgs::argv[0]), globalRepository( makeAboutDefault(PythonArgs::argv[0]).appName() ) )
{
}

// Constructor
Environment::Environment( pybind11::list arg, po::options_description const& desc, Repository::Config const& config )
    :
    Environment( PythonArgs(arg).argc, PythonArgs::argv, mpi::threading::single, desc, feel_options(), makeAboutDefault(PythonArgs::argv[0]), config )
{
}

#if 0
Environment::Environment( boost::python::list arg )
    :
#if BOOST_VERSION >= 105500
    Environment( PythonArgs(arg).argc, PythonArgs::argv, mpi::threading::single, feel_nooptions(), feel_options(), makeAboutDefault(PythonArgs::argv[0]), makeAboutDefault(PythonArgs::argv[0]).appName() )
#else
    Environment( PythonArgs(arg).argc, PythonArgs::argv, feel_nooptions(), feel_options(), makeAboutDefault(PythonArgs::argv[0]), makeAboutDefault(PythonArgs::argv[0]).appName() )
#endif
{
}
#endif // 0

#endif // FEELPP_HAS_PYTHON

#if defined ( FEELPP_HAS_PETSC_H )
void
Environment::initPetsc( int * argc, char *** argv )
{
    PetscTruth is_petsc_initialized;
    PetscInitialized( &is_petsc_initialized );

    if ( !is_petsc_initialized )
    {
        int ierr;
        if( (*argc > 0) && ( argv != nullptr ) )
            ierr = PetscInitialize( argc, argv, PETSC_IGNORE, PETSC_IGNORE );
        else
            ierr = PetscInitializeNoArguments();
        CHKERRABORT( *S_worldcomm,ierr );
    }

    // make sure that petsc do not catch signals and hence do not print long
    //and often unuseful messages
    PetscPopSignalHandler();

#if defined( FEELPP_HAS_SLEPC )
    PetscBool is_slepc_initialized;
    SlepcInitialized( &is_slepc_initialized );
    if ( !is_slepc_initialized )
    {
        int ierr;
        if( (*argc > 0) && ( argv != nullptr ) )
            ierr = SlepcInitialize( argc, argv, PETSC_IGNORE, PETSC_IGNORE );
        else
            ierr = SlepcInitializeNoArguments();
        CHKERRABORT( *S_worldcomm,ierr );
    }
#endif
}
#endif // FEELPP_HAS_PETSC_H

// Constructor
Environment::Environment( int argc, char** argv,
                          mpi::threading::level lvl,
                          po::options_description const& desc,
                          po::options_description const& desc_lib,
                          AboutData const& about,
                          Repository::Config const& config )
{
    // we can initialized only once
    if (S_initialized)
        return;
    S_initialized = true;

    if ( argc == 0 )
    {
#if BOOST_VERSION >= 105500
        M_env = std::make_unique<boost::mpi::environment>(lvl, false);
#else
        M_env = std::make_unique<boost::mpi::environment>(false);
#endif
    }
    else
    {
#if BOOST_VERSION >= 105500
        M_env = std::make_unique<boost::mpi::environment>(argc, argv, lvl, false);
#else
        M_env = std::make_unique<boost::mpi::environment>(argc, argv, false);
#endif
    }
    CHECK( M_env->initialized()) << "MPI environment failed to initialize properly.";
    S_argc = argc;
    S_argv = argv;

#if defined(FEELPP_HAS_SPDLOG)
    // Immediately disable console logging by setting up a null sink as default
    // This prevents early LOG() calls from appearing on console
    // The proper logger (with file/console based on options) will be set up later in startLogging()
    Logger::setDefaultLogger(Logger::createNullLogger("feelpp_early"));
#endif

    //
    // setup worldcomm
    //
    S_worldcomm = worldcomm_type::New();
    CHECK( S_worldcomm ) << "Feel++ Environment: creating worldcomm failed!";
    S_worldcommSeq.reset( new WorldComm( S_worldcomm->subWorldCommSeq() ) );
    cout.attachWorldComm( S_worldcomm );
    cerr.attachWorldComm( S_worldcomm );
    clog.attachWorldComm( S_worldcomm );

    //
    // setup repository 
    //
    S_repository = Repository( config );
#if 0    
    if ( !S_repository.isCustom() )
    {
        //S_repository.configure();
        S_rootdir = S_repository.root();
        S_appdir = S_repository.directory();
        S_appdirWithoutNumProc = S_repository.directoryWithoutAppenders();
    }
    else
    {
        S_rootdir.clear();
        S_appdir.clear();
        S_appdirWithoutNumProc.clear();
    }
#endif
    //
    // setup options
    //
    S_desc_app = std::make_shared<po::options_description>( desc );
    S_desc_lib = std::make_shared<po::options_description>( desc_lib );
    S_desc = std::make_shared<po::options_description>();
    S_desc->add( *S_desc_app );


    // try to see if the feel++ lib options are already in S_desc_app, if yes then we do not add S_desc_lib
    // otherwise we will have duplicated options
    std::vector<boost::shared_ptr<po::option_description>> opts = Environment::optionsDescriptionApplication().options();
    auto it = std::find_if( opts.begin(), opts.end(),
                            []( boost::shared_ptr<po::option_description> const&o )
                            {
                                return o->format_name().erase( 0,2 ) == "backend";
                            } );

    if   ( it == opts.end() )
        S_desc->add( *S_desc_lib );

    S_desc->add( file_options( about.appName() ) );
    S_desc->add( generic_options() );
    S_about = about;

    // duplicate argv before passing to gflags because gflags is going to
    // rearrange them and it screws badly the flags for PETSc/SLEPc
    char** envargv = dupargv( argv );

    //
    // Initialize PETSc
    //
#if defined ( FEELPP_HAS_PETSC_H )
    initPetsc( &argc, &envargv );
#endif

    //
    // Initialize Gmsh
    //
#if defined( FEELPP_HAS_GMSH_H )
#if defined( FEELPP_HAS_GMSH_API )
    gmsh::initialize();
#else
    GmshInitialize();
#endif
#endif
#if defined(FEELPP_HAS_PYTHON)
    if ( !Py_IsInitialized() )
    {
        py::initialize_interpreter();
        S_init_python = true;
    }
    else
        S_init_python = false;
#else
    S_init_python = false;
#endif

    cout << "[ Feel++ ] "
         << "application " << about.appName()
         << " version " << about.version()
         << " initializing..." << std::endl;
    //
    // parse options
    //
    doOptions( argc, envargv, *S_desc, *S_desc_lib, about.appName() );

    // Enable auto mode for all observers.
    Environment::setJournalEnable( boption(_name="journal") );
    Environment::setJournalAutoMode( boption(_name="journal.auto") );

    // Force environment to connect to the journal.
    S_informationObject = std::make_unique<JournalWatcher>( std::bind( &Environment::updateInformationObject, this, std::placeholders::_1 ), "Environment", "", false );

    S_timers = std::make_unique<TimerTable>();

    auto today = std::chrono::system_clock::now();
    tic();
    Logger::console()->info("[ Starting Feel++ ] application {} version {} date {:%Y-%m-%d}", about.appName(), about.version(), today);
    Logger::console()->flush();

    //
    // setup work directory
    //
    fs::path directory;
    if ( S_vm.count( "directory" ) )
        directory = expand(S_vm["directory"].as<std::string>());
    if ( S_vm.count( "repository.prefix" ) )
        directory = expand(S_vm["repository.prefix"].as<std::string>());
    if ( S_vm.count( "repository.case" ) )
    {
        fs::path d{ directory };
        d /= expand(S_vm["repository.case"].as<std::string>());
        directory = d.string();
    }
    // For custom repositories, changeRepository() will invoke the callback
    // The directory parameter is only used for non-custom repositories
    changeRepository( _directory = boost::format{ directory.string() } );
    //
    // use --dirs to check the directories of Feel++ environment
    //
    if ( S_vm.count( "dirs" ) == 1 )
    {
        cout<< "- root: " << rootRepository() << std::endl
            << "-  app: " << appRepository() << std::endl
            // << "-  geo: " << geoRepository() << std::endl
            << "- home: " << expand("$home") << std::endl
            << "-  cfg: " << expand("$cfgdir") << std::endl
            << "- data: " << expand("$datadir") << std::endl;
        if ( Environment::initialized() )
        {
            worldComm().barrier();
            MPI_Finalize();
        }
#if defined(FEELPP_HAS_MONGOCXX )
        MongoCxx::reset();
#endif
        exit( 0 );
    }

    //
    // setup journal
    //
    if( S_vm.count( "journal.filename" ) )
    {
        // TODO relative or absolute path
        Environment::setJournalFilename( (fs::path( Environment::appRepository() )/fs::path(soption(_name="journal.filename"))).string() );
    }
    else
        Environment::setJournalFilename( (fs::path( Environment::appRepository() )/fs::path("journal.json")).string() );

    //
    // setup MongoDB
    //
#if defined(FEELPP_HAS_MONGOCXX )
    MongoConfig journaldbconf;
    if( S_vm.count( "journal.database.name" ) )
       journaldbconf.name = S_vm["journal.database.name"].as<std::string>();
    if( S_vm.count( "journal.database.host" ) )
       journaldbconf.host = S_vm["journal.database.host"].as<std::string>();
    if( S_vm.count( "journal.database.port" ) )
       journaldbconf.port = S_vm["journal.database.port"].as<std::string>();
    if( S_vm.count( "journal.database.user" ) )
    {
       journaldbconf.user = S_vm["journal.database.user"].as<std::string>();
    }
    else
    {
        if( auto user = std::getenv( "FEELPP_DB_JOURNAL_USER" ) )
            journaldbconf.user = user;
    }

    if( S_vm.count( "journal.database.password" ) )
    {
        std::string password = S_vm["journal.database.password"].as<std::string>();
        if( S_vm["journal.database"].as<bool>() )
        {
            // TODO Fix in parallel
            if( password == "?" )
                password = askPassword("Enter your mongodb password:");
            journaldbconf.password = password;
        }
    }
    // Environment variable.
    else
    {
        if( auto password = std::getenv( "FEELPP_DB_JOURNAL_PASSWORD" ) )
            journaldbconf.password = password;
    }
    if( S_vm.count( "journal.database.authsrc" ) )
        journaldbconf.authsrc = S_vm["journal.database.authsrc"].as<std::string>();
    if( S_vm.count( "journal.database.collection" ) )
        journaldbconf.collection = S_vm["journal.database.collection"].as<std::string>();
    Environment::journalDBConfig( journaldbconf );
    Feel::MongoCxx::instance();
#endif

    if ( not S_hwSysInstance )
    {
#if defined(FEELPP_HAS_KWSYS )
        // Use kwsys library.
        S_hwSysInstance = std::make_unique<Sys::KWSys>();
#else
        S_hwSysInstance = std::make_unique<Sys::HWSys>();
#endif
    }

#if defined( FEELPP_HAS_TBB )
    int n = tbb::task_scheduler_init::default_num_threads();
    //int n = 2;
    //VLOG(2) << "[Feel++] TBB running with " << n << " threads\n";
    //tbb::task_scheduler_init init(2);
#endif

    // Note: verbosity is now handled in Environment::startLogging() for spdlog
    // or by glog directly when FEELPP_HAS_SPDLOG is not defined
#if !defined(FEELPP_HAS_SPDLOG)
    if ( S_vm.count( "v" ) )
        Environment::logVerbosityLevel() = S_vm["v"].as<int>();
    if ( S_vm.count( "vmodule" ) )
    {
        //Environment::logVerbosityLevel()module = S_vm["vmodule"].as<std::string>();
        //google::SetVLOGLevel( "*btpcd", 2 );
    }
#endif

#if 0
    if ( S_vm.count( "nochdir" ) == 0 )
    {
        if ( S_vm.count( "directory" ) )
            directory = S_vm["directory"].as<std::string>();

        LOG( INFO ) << "change directory to " << directory << "\n";
        boost::format f( directory );
        bool createSubdir = add_subdir_np && S_vm["npdir"].as<bool>();
        changeRepository( _directory=f,_subdir=createSubdir );
    }
#endif //0

    freeargv( envargv );

    /* Initialize hwloc topology */
    /* to extract info about architecture */
#if defined(FEELPP_HAS_HARTS)
    Environment::initHwlocTopology();
#endif

    Environment::journalCheckpoint();

    //::detail::DebugWait( worldComm().globalRank() );
}
void
Environment::clearSomeMemory()
{
    Environment::logMemoryUsage( "Environment::clearSomeMemory before:" );

    // send signal to all deleters
    S_deleteObservers();
#if defined(FEELPP_HAS_SPDLOG)
    Logger::flushOn(0);
#else
    google::FlushLogFiles( google::GLOG_INFO );
#endif
    VLOG( 2 ) << "clearSomeMemory: delete signal sent" << "\n";

    Environment::logMemoryUsage( "Environment::clearSomeMemory after:" );
}
bool
Environment::shouldLog()
{
    if ( !Environment::initialized() )
        return true;

    std::string mode = Environment::logMpiMode();
    int rank = Environment::rank();

    if ( mode == "none" )
        return false;
    else if ( mode == "master" )
        return ( rank == 0 );
    // else mode == "all": all ranks log
    return true;
}
// Destructor.
Environment::~Environment()
{
    if ( boption( _name="display-stats" ) )
        Environment::saveTimers( true );

    double t = toc("env",no_display);
    Table summary;
    summary.add_row( { S_about.appName() } );
    summary( 0, 0 ).format().setFontAlign( Font::Align::center );
    Table data;
    data.add_row( { "logs", Environment::logsRepository() } );
    if ( Environment::journalEnabled()  )
        data.add_row( { "journal", Environment::journalFilename() } );
    Table paths;
    paths.add_row( { Environment::appRepository() } );
    for ( auto p : S_paths | ranges::views::take_exactly( S_paths.size() - 3 ) )
    {
        paths.add_row({p.string()});
    }
    data.add_row( { "directories", paths } );
    summary.add_row({data});
    cout << summary << std::endl;
    cout << "[ Stopping Feel++ ] " << tc::green << "application " << S_about.appName()
         << " execution time " << t << "s" << tc::reset << std::endl;

#if defined(FEELPP_HAS_HARTS)
    /* if we used hwloc, we free topology data */
    Environment::destroyHwlocTopology();
#endif

    /* if we were using onelab */
    /* we write the file containing the filename marked for automatic loading in Gmsh */
    /* we serialize the writing of the size by the different MPI processes */


    //std::cout << S_vm["onelab.enable"].as<int>() << std::endl;

    // Only execute onelab cleanup if MPI is still active (not finalized)
    if ( ioption( _name="onelab.enable" ) == 2 && initialized() && !finalized() )
    {
        for ( int i = 0; i < worldComm().size(); i++ )
        {
            /* only one process at a time */
            if ( i == worldComm().globalRank() )
            {
                std::cout << Environment::olAppPath << std::endl;
                int i;
                std::ofstream ool;

                /* Generate a file containing the name of the outputs for the current dataset */
                /* eother truncate the file if we are process 0 or complete it if we are an other process */
                if ( worldComm().globalRank() == 0 )
                {
                    ool.open( Environment::olAppPath + ".onelab.out", std::ofstream::out | std::ofstream::trunc );
                }

                else
                {
                    ool.open( Environment::olAppPath + ".onelab.out", std::ofstream::out | std::ofstream::app );
                }

                fs::path p( Environment::olAppPath );

                /* If we have dataset to load */
                /* we add each of them to the file containing the files to load */
                if ( Environment::olAutoloadFiles.size() > 0 )
                {
                    // Files marked for autoloading
                    ool << "#";

                    for ( i = 0; i < Environment::olAutoloadFiles.size(); i++ )
                    {
                        ool << " ";

                        if ( S_vm.count( "onelab.remote" ) && S_vm["onelab.remote"].as<std::string>() != "" )
                        {
                            ool << S_vm["onelab.remote"].as<std::string>() << ":";
                        }

                        ool << p.parent_path().string() << "/" << Environment::olAutoloadFiles[i];
                    }

                    ool << std::endl;

                    i = 0;
                    ool << "FeelApp.merge(" << Environment::olAutoloadFiles[i];

                    for ( i = 1; i < Environment::olAutoloadFiles.size(); i++ )
                    {
                        ool << ", " << Environment::olAutoloadFiles[i];
                    }

                    ool << ");" << std::endl;
                }

                else
                {
                    std::cout << worldComm().globalRank() << " No files to load" << std::endl;
                }

                ool.close();
            }

            /* wait for the current process to finish */
            Environment::worldComm().barrier();
        }
    }

    VLOG( 2 ) << "[~Environment] sending delete to all deleters" << "\n";

    Environment::clearSomeMemory();

#if defined(FEELPP_HAS_PYTHON)
    if ( S_init_python )
        py::finalize_interpreter();
#endif
#if defined(FEELPP_HAS_MONGOCXX )
    VLOG( 2 ) << "cleaning mongocxxInstance";
    MongoCxx::reset();
#endif

#if 0
#if defined( FEELPP_HAS_GMSH_H )
#if defined( FEELPP_HAS_GMSH_API )
    gmsh::finalize();
#else
    GmshFinalize();
#endif
#endif
#endif

    VLOG( 2 ) << "clearing known paths\n";
    S_paths.clear();

    VLOG( 2 ) << "[~Environment] cleaning up global excompiler\n";
    GiNaC::cleanup_ex( false );

    VLOG( 2 ) << "[~Environment] finalizing slepc,petsc and mpi\n";
#if defined ( FEELPP_HAS_PETSC_H )
    PetscTruth is_petsc_initialized;
    PetscInitialized( &is_petsc_initialized );
    if ( is_petsc_initialized && !Environment::aborted() )
    {
#if defined( FEELPP_HAS_SLEPC )
        SlepcFinalize();
#else
        PetscFinalize();
#endif // FEELPP_HAS_SLEPC
    }
#endif // FEELPP_HAS_PETSC_H

    stopLogging();

    JournalManager::journalFinalize();
    S_timers.reset();
    S_hwSysInstance.reset(); // call deleter
    S_informationObject.reset();

    // make sure everybody is here (only if MPI is still active)
    if ( !Environment::aborted() && initialized() && !finalized() )
        Environment::worldComm().barrier();
    // Handle --rm cleanup: only master rank should remove files to avoid race conditions
    // We check initialized() to ensure MPI is still valid before checking rank
    if ( Environment::isMasterRank() && S_vm.count("rm") )
    {
        // Simplified cleanup without MPI rank check (unsafe after PetscFinalize)
        // This will run on all ranks but that's safer than crashing
        Logger::console()->info("Removing files (--rm) in {}...", appRepository());
        Logger::console()->flush();
        fs::remove_all( S_appdir );
        // should remove expression dir
        fs::remove_all( S_repository.exprs() );
        if ( fs::exists( S_repository.root()/"crbdb"/S_about.appName()))
        {
            Logger::console()->info("Removing files (--rm) in {}", S_repository.root()/"crbdb"/S_about.appName());
            Logger::console()->flush();
            fs::remove_all( S_repository.root()/"crbdb"/S_about.appName() );
        }
    }

    if ( !Environment::aborted() )
    {
    // call gmsh::finalize() at the end because if gmsh is compiled with MPI support, the gmsh lib call MPI_Finalize
#if defined( FEELPP_HAS_GMSH_H )
#if defined( FEELPP_HAS_GMSH_API )
        gmsh::finalize();
#else
        GmshFinalize();
#endif
#endif
    }
}


void
Environment::generateOLFiles( int argc, char** argv, std::string const& appName )
{
    //Application path
    int i;
    bool isNum = false;
    fs::path p( argv[0] );
    std::ostringstream appPath;

    /* get app name */
    appPath.str( "" );
    appPath << fs::absolute( p ).string();

    std::ostringstream optionPath;

    std::ofstream ol;
    ol.open( appPath.str() + ".ol", std::ofstream::out | std::ofstream::trunc ); //.ol file
    std::ofstream cfgol;
    cfgol.open( appPath.str() + ".onelab.cfg.ol", std::ofstream::out | std::ofstream::trunc ); //.cfg.ol file

    /* map from feel option name to onelab option path */
    std::map<std::string, std::string> mOptToOptPath;

    std::map<std::string, std::vector<boost::shared_ptr<po::option_description> >> moptions
    {
        {"Feelpp", Environment::optionsDescriptionLibrary().options() },
        {S_about.appName(), Environment::optionsDescriptionApplication().options() },
    };

    for ( auto o : moptions )
    {
        for ( boost::shared_ptr<po::option_description> option : o.second )
        {
            //Information about the option
            std::string optName = option->format_name().erase( 0,2 ); //Putting the option name in a variable for easier manipulations
            std::string defVal = ""; //option->format_parameter(); //Putting the option default value in a variable for easier manipulations
            std::string desc=option->description(); // Option description

            // reset option path
            optionPath.str( "" );
            // reset type
            isNum = false;

            //std::cout << optName << ";" << defVal << ";" << desc << std::endl;

            std::string ens,funcName;

            std::vector<std::string> strings; //Vector of the split name
            boost::split( strings,optName,boost::is_any_of( "." ) ); //Spliting option name
            ens = "";

            if ( strings.size() > 1 )
            {
                ens = strings[0] + "/"; //Getting the first split element for the option set in the .cfg.ol file

                for ( size_t i = 1; i < strings.size() - 1; i++ ) //Getting the option set
                {
                    ens += strings[i] + "/";
                }
            }

            funcName = strings[strings.size() - 1]; //Raw option name

            /* skip some options */
            /* they won't be displayed in Gmsh */
            if ( funcName == "config-file" )
            {
                continue;
            }

            /* if an option has been set either through command line */
            /* or through the initial config file */
            /* we use its configuration */
            //if(S_vm.count(optName) && !(S_vm[optName].defaulted()))
            if ( S_vm.count( optName ) )
            {
                std::ostringstream oss;
                oss.str( "" );
                //std::cout << defVal;

                //std::cout << "Entry for " << optName << ": ";
                // if the option if defaulted and soesn't starts with onelab,
                // we put it in the end of Gmsh options
                if ( S_vm[optName].defaulted() && optName.find( "onelab." ) == std::string::npos )
                {
                    //std::cout << "defaulted ";
                    optionPath << "GeneralParameters/" << o.first << "/" << ens;
                }

                // if we have a user defined option or a onelab option
                // we want them to be on top of the list for easier access
                else
                {
                    optionPath << "DefinedParameters/" << o.first << "/" << ens;
                }

                //optionPath << "Parameters/" << ens;

                if ( optName == "licence" )
                {
                    const std::type_info & ti = S_vm[optName].value().type();
                    std::cout << ti.name() << " " << std::endl;
                }

                if ( S_vm[optName].empty() )
                {
                    //std::cout << "empty ";
                }
                else
                {
                    const std::type_info & ti = S_vm[optName].value().type();

                    //std::cout << ti.name() << " ";
                    if ( ti == typeid( bool ) )
                    {
                        oss.str( "" );
                        oss << ( S_vm[optName].as<bool>() ? "1" : "0" );
                        isNum = false;
                    }

                    else if ( ti == typeid( int ) )
                    {
                        oss.str( "" );
                        oss << S_vm[optName].as<int>();
                        isNum = true;
                    }

                    else if ( ti == typeid( size_type ) )
                    {
                        oss.str( "" );
                        oss << S_vm[optName].as<size_type>();
                        isNum = true;
                    }

                    else if ( ti == typeid( float ) )
                    {
                        oss.str( "" );
                        oss << S_vm[optName].as<float>();
                        isNum = true;
                    }

                    else if ( ti == typeid( double ) )
                    {
                        oss.str( "" );
                        oss << S_vm[optName].as<double>();
                        isNum = true;
                    }

                    else if ( ti == typeid( std::string ) )
                    {
                        oss.str( "" );
                        oss <<  S_vm[optName].as<std::string>();
                        isNum = false;
                    }

                    else
                    {
                        std::cout << "Unknown type for parameter " << optName << "(" << typeid( void ).name() << ")" << std::endl;
                        isNum = false;
                    }
                }

                //std::cout << oss.str() << std::endl;
                defVal = oss.str();

                /* Force Gmsh as a the default exporter */
                /* as we are using OneLab */
                if ( ens == "exporter/" && funcName == "format" )
                {
                    defVal = "gmsh";
                }

                if ( optName == "onelab.enable" )
                {
                    defVal = "2";
                }

                /*
                   if(defVal.size() != 0) //Excluding options without a default value
                   {
                 */
                if ( isNum )
                {
                    ol << funcName << ".number(" << defVal << ", " << optionPath.str() << ");" << " # "<< desc << std::endl;
                    cfgol << optName << "=OL.get(" << optionPath.str() << funcName << ")" << std::endl;
                }

                else
                {
                    ol << funcName << ".string(" << defVal << ", " << optionPath.str() << ");" << " # "<< desc << std::endl;
                    cfgol << optName << "=OL.get(" << optionPath.str() << funcName << ")" << std::endl;
                }

                //}

                /* Hide some options from users */
                if ( optName == "onelab.enable"
                        || optName == "onelab.remote"
                        || optName == "onelab.sync.script" )
                {
                    ol << funcName << ".setVisible(0);" << std::endl;
                }

                ol << funcName << ".setReadOnly(0);" << std::endl;

                /* store some option paths for building ol script */
                if ( optName == "onelab.chroot"
                        || optName == "onelab.remote"
                        || optName == "onelab.np"
                        || optName == "onelab.sync.script" )
                {
                    mOptToOptPath[optName] = optionPath.str() + funcName;
                }

            }

        }
    }

    ol << "" << std::endl;

    /* Mesher instructions */
    ol << "Mesher.register(native," << stringize( GMSH_EXECUTABLE ) << ");" << std::endl;
    ol << "OL.if(OL.get(Parameters/gmsh/filename) == untitled.geo)" << std::endl;
    ol << "OL.msg(No geo file specified. Using a default one);" << std::endl;
    ol << "OL.endif" << std::endl;

    if ( S_vm.count( "onelab.remote" )
            && S_vm["onelab.remote"].as<std::string>() != ""
            && S_vm["onelab.remote"].as<std::string>() != "localhost" )
    {
        ol << "FeelApp.remote(" << "OL.get(" + mOptToOptPath["onelab.remote"] << "), " << p.parent_path().string() << "/" << ");" << std::endl;

        ol << "FeelApp.register(interfaced, ./OL.get(Arguments/FileName).onelab.py);" << std::endl;

        ol << "FeelApp.in(OL.get(Arguments/FileName).onelab.cfg.ol);" << std::endl;

        /* test for chroots */
        ol << "OL.if(OL.get(" << mOptToOptPath["onelab.chroot"] << "))" << std::endl;
        ol << "FeelApp.run( schroot -c OL.get(" << mOptToOptPath["onelab.chroot"] << ") -- ";
        ol << stringize( MPIEXEC ) << " " << stringize( MPIEXEC_NUMPROC_FLAG ) << " OL.get(" << mOptToOptPath["onelab.np"] << ") " << appPath.str();
        ol << " --config-file OL.get(Arguments/FileName).onelab.cfg --nochdir );" << std::endl;
        ol << "OL.else" << std::endl;
        ol << "FeelApp.run(" << stringize( MPIEXEC ) << " " << stringize( MPIEXEC_NUMPROC_FLAG ) << " OL.get(" << mOptToOptPath["onelab.np"] << ") " << appPath.str();
        ol << " --config-file OL.get(Arguments/FileName).onelab.cfg --nochdir );" << std::endl;
        ol << "OL.endif" << std::endl;

        ol << "FeelApp.out(OL.get(Arguments/FileName).onelab.out);" << std::endl;

        ol << "SyncData.register(interfaced, OL.get(" << mOptToOptPath["onelab.sync.script"] << "));" << std::endl;
        ol << "SyncData.in(OL.get(Arguments/FileName).onelab.out);" << std::endl;
        ol << "SyncData.run(OL.get(Arguments/FileName).onelab.out);" << std::endl;

        ol << "OL.include(OL.get(Arguments/FileName).onelab.out);" << std::endl;
    }

    else
    {
        /* setup chroot */
        ol << "FeelApp.register(interfaced, " << appPath.str() << ".onelab.py);" << std::endl;

        std::string cpath = "";

        if ( S_vm.count( "onelab.remote" )
                && ( S_vm["onelab.remote"].as<std::string>() == ""
                     || S_vm["onelab.remote"].as<std::string>() == "localhost" ) )
        {
            cpath = fs::current_path().string();
            size_t n = std::count( cpath.begin(), cpath.end(), '/' );
            cpath = "";

            for ( int i = 0; i < n; i++ )
            {
                cpath = cpath + "../";
            }
        }

        ol << "FeelApp.in(" << cpath << appPath.str() << ".onelab.cfg.ol);" << std::endl;

        /* test for chroots */
        ol << "OL.if(OL.get(" << mOptToOptPath["onelab.chroot"] << "))" << std::endl;
        ol << "FeelApp.run( schroot -c OL.get(" << mOptToOptPath["onelab.chroot"] << ") -- ";
        ol << stringize( MPIEXEC ) << " " << stringize( MPIEXEC_NUMPROC_FLAG ) << " OL.get(" << mOptToOptPath["onelab.np"] << ") " << appPath.str();
        ol << " --config-file " << appPath.str() << ".onelab.cfg --nochdir );" << std::endl;
        ol << "OL.else" << std::endl;
        ol << "FeelApp.run(" << stringize( MPIEXEC ) << " " << stringize( MPIEXEC_NUMPROC_FLAG ) << " OL.get(" << mOptToOptPath["onelab.np"] << ") " << appPath.str();
        ol << " --config-file " << appPath.str() << ".onelab.cfg --nochdir );" << std::endl;
        ol << "OL.endif" << std::endl;

        ol << "FeelApp.out( " << cpath << appPath.str() << ".onelab.out);" << std::endl;

        ol << "OL.include(" << cpath << appPath.str() << ".onelab.out);" << std::endl;
    }

    ol.close();
    cfgol.close();

    /* generate a script for executing the application */
    /* to avoid patching Gmsh */
    std::string pyscript = appPath.str() + ".onelab.py";
    std::ofstream shs;
    shs.open( pyscript, std::ofstream::out | std::ofstream::trunc );

    shs << "#!/usr/bin/python" << std::endl;
    shs << "import sys, subprocess" << std::endl << std::endl;

    shs << "def main():" << std::endl;

    shs << "  cmd = sys.argv[1:]" << std::endl;
    shs << "  print cmd" << std::endl;
    shs << "  retval = subprocess.call(cmd)" << std::endl;
    shs << "  return retval" << std::endl;

    shs << "main()" << std::endl;

    shs.close();

    chmod( pyscript.c_str(), S_IRWXU|S_IRGRP|S_IROTH );

}

void
Environment::setLogVerbosityLevel( int v )
{
    LOG(INFO) << fmt::format( "set log verbosity level to {}, previously {}", v, Environment::logVerbosityLevel() );
#if defined(FEELPP_HAS_SPDLOG)
    Logger::verbosity() = v;
    Logger::setLevel(v);
#else
    Environment::logVerbosityLevel() = v;
#endif
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

    if ( worldComm().isMasterRank() )
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
                S_vm.count( "help-lib" ) ||
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
            std::cout << optionsDescriptionApplication() << "\n";
            std::cout << file_options( S_about.appName() ) << "\n";
            std::cout << generic_options() << "\n";
        }

        if ( S_vm.count( "help-lib" ) )
        {
            std::cout << optionsDescriptionLibrary() << "\n";
            std::cout << file_options( S_about.appName() ) << "\n";
            std::cout << generic_options() << "\n";
        }
    }

    if ( S_vm.count( "verbose" ) ||
            S_vm.count( "help" ) ||
            S_vm.count( "help-lib" ) ||
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
#if defined(FEELPP_HAS_MONGOCXX )
        MongoCxx::reset();
#endif
        exit( 0 );
    }

#if 0
    std::cout << "count = " << S_vm.count( "debug" ) << "\n"
              << "string = " << S_vm["debug"].as<std::string>() << "\n";
#endif

    VLOG( 2 ) << "[processGenericOptions] done\n";
}

void
Environment::parseAndStoreOptions( po::command_line_parser parser, bool extra_parser )
{
    VLOG( 2 ) << " parsing options...\n";

    std::shared_ptr<po::parsed_options> parsed;

    if ( extra_parser )
    {
        parsed = std::shared_ptr<po::parsed_options>( new po::parsed_options( parser
                 .options( *S_desc )
                 .extra_parser( at_option_parser_2 )
                 .allow_unregistered()
                 .run() ) );
    }


    else
    {
        parsed = std::shared_ptr<po::parsed_options>( new po::parsed_options( parser
                 .options( *S_desc )
                 .allow_unregistered()
                 .run() ) );
    }

    VLOG( 2 ) << "[parseAndStoreOptions] parsing options done\n";

    S_to_pass_further = po::collect_unrecognized( parsed->options, po::include_positional );

    if ( Environment::isMasterRank() && S_to_pass_further.size() )
    {
        LOG( ERROR ) << "Some options (" << ( S_to_pass_further.size() ) << ") were not recognized.";
        LOG( ERROR ) << "We remove them from Feel++ options management system and pass them to PETSc/SLEPc";
        LOG( ERROR ) << "and other third party libraries";

        for ( std::string const& s: S_to_pass_further )
        {
            LOG( ERROR ) << "  |- unrecognized option: " << s << "\n";
        }
    }
    std::vector<po::basic_option<char> >::iterator it = parsed->options.begin();
    //std::vector<po::basic_option<char> >::iterator en  = parsed->options.end();
    for ( ; it != parsed->options.end() ; )
    {
        if ( it->unregistered )
        {
            if ( Environment::isMasterRank() )
                LOG( ERROR ) << "  |- remove " << it->string_key << " from Feel++ options management system"  << "\n";
            it = parsed->options.erase( it );
        }
        else
            ++it;
    }

    po::store( *parsed, S_vm );

    if ( boption( _name="fail-on-unknown-option" ) && S_to_pass_further.size() )
    {
        std::stringstream ostr;

        for ( std::string const& s: S_to_pass_further )
        {
            ostr << s << " ";
        }

        if ( Environment::isMasterRank() )
            LOG( ERROR ) << "Unknown options [" << ostr.str() << "] passed to Feel++. Quitting application...";

        //MPI_Barrier( S_worldcomm->comm() );
        MPI_Abort( S_worldcomm->comm(), 1 );
    }
}


std::string
Environment::findFileRemotely( std::string const& fname, std::string const& subdir )
{
    RemoteData rdTool( fname, worldCommPtr());
    if ( rdTool.canDownload() )
    {
        // Use temp directory if rootRepository hasn't been configured yet
        fs::path downloadRoot = S_rootdir.empty() ? fs::temp_directory_path() / "feelpp" : S_rootdir;
        auto downloadedFolder = rdTool.download( (downloadRoot/fs::path("downloads")/fs::path(Environment::about().appName())/fs::path(subdir)).string() );
        for( auto dl : downloadedFolder )
            std::cout << dl << std::endl;

        //CHECK( downloadedFolder.size() == 1 ) << "download only one folder";
        return downloadedFolder[0];
    }
    return fname;
}


void
Environment::doOptions( int argc, char** argv,
                        po::options_description const& desc,
                        po::options_description const& desc_lib,
                        std::string const& appName )
{
    //std::locale::global(std::locale(""));
    try
    {
        S_commandLineParser = std::shared_ptr<po::command_line_parser>( new po::command_line_parser( argc, argv ) );
        parseAndStoreOptions( po::command_line_parser( argc, argv ), true );
        processGenericOptions();

        VLOG( 2 ) << "options parsed and stored in database";
        S_log_mpi_mode = S_vm["log.mpi"].as<std::string>();
        std::vector<std::string> configFiles;

        if ( S_vm.count( "case" ) )
        {
            std::vector<std::string> cfgsInCaseDir;
            std::string caseDir = S_vm["case"].as<std::string>();
            RemoteData rdTool( caseDir, worldCommPtr());
            if ( rdTool.canDownload() )
            {
                // Use temp directory if rootRepository hasn't been configured yet
                fs::path downloadRoot = S_rootdir.empty() ? fs::temp_directory_path() / "feelpp" : S_rootdir;
                auto downloadedFolder = rdTool.download( (downloadRoot/fs::path("downloads")/fs::path(appName)/fs::path("cases")).string() );
                CHECK( downloadedFolder.size() == 1 ) << "download only one folder";
                caseDir = downloadedFolder[0];
            }

            fs::path fscaseDir( caseDir );
            CHECK( fs::is_directory( fscaseDir ) ) << "case must be a directory";
            std::string dirName = fscaseDir.filename().string();
            if ( Feel::filename_is_dot( fscaseDir.filename() ) )
                dirName = fscaseDir.parent_path().filename().string();
            std::string caseConfigFile = dirName + ".cfg";
            if ( S_vm.count( "case.config-file" ) )
                caseConfigFile = fs::path(S_vm["case.config-file"].as<std::string>()).filename().string();

            fs::directory_iterator end_itr;
            for ( fs::directory_iterator itr( fscaseDir ); itr != end_itr; ++itr )
            {
                if ( fs::is_regular_file( itr->path() ) && ( itr->path().extension() == ".cfg" ) )
                {
                    cfgsInCaseDir.push_back( itr->path().string() );
                }
            }
            if ( cfgsInCaseDir.size() == 1 )
                configFiles.push_back( cfgsInCaseDir.front() );
            else
            {
                for ( std::string const& cfgFile : cfgsInCaseDir )
                {
                    if ( fs::path( cfgFile ).filename().string() == caseConfigFile )
                    {
                        configFiles.push_back( cfgFile );
                        break;
                    }
                }
            }
        }
         // parse config file if given to command line
        if ( S_vm.count( "config-file" ) || S_vm.count( "config-files" ) )
        {
            std::vector<std::string> configFilesFromCmd;
            if ( S_vm.count( "config-file" ) )
                configFilesFromCmd.push_back( S_vm["config-file"].as<std::string>() );
            if ( S_vm.count( "config-files" ) )
            {
                std::vector<std::string> configFilesOptVec = S_vm["config-files"].as<std::vector<std::string> >();
                configFilesFromCmd.insert( configFilesFromCmd.end(), configFilesOptVec.begin(), configFilesOptVec.end() );
            }
            for ( std::string const& cfgFile : configFilesFromCmd )
            {
                RemoteData rdTool( cfgFile, worldCommPtr());
                if ( rdTool.canDownload() )
                {
                    // Use temp directory if rootRepository hasn't been configured yet
                    fs::path downloadRoot = S_rootdir.empty() ? fs::temp_directory_path() / "feelpp" : S_rootdir;
                    auto dowloadedData = rdTool.download( (downloadRoot/fs::path("downloads")/fs::path(appName)/fs::path("cfgs")).string() );
                    for ( std::string const& data : dowloadedData )
                        configFiles.push_back( data );
                }
                else
                    configFiles.push_back( cfgFile );
            }
        }
#if 0
        std::cout << "CONFIG-FILES\n";
        for ( std::string const& cfgFile : configFiles )
            std::cout << cfgFile << "\n";
#endif
        // Use setConfigFiles for consistent config file handling
        if ( !configFiles.empty() )
        {
            setConfigFiles( configFiles );
        }
        else
        {
            po::notify( S_vm );
        }



        /* handle the generation of onelab files after having processed */
        /* the regular config file, so we have parsed user defined parameters */
        /* or restored a previous configuration */

        /* We store the application path for further use */
        fs::path p( argv[0] );
        Environment::olAppPath = fs::absolute( p ).string();

        if ( worldComm().isMasterRank() )
        {
            if ( S_vm.count( "onelab.enable" ) )
            {

                if ( S_vm["onelab.enable"].as<int>() == 1 )
                {
                    Environment::generateOLFiles( argc, argv, appName );

                    if ( Environment::initialized() )
                    {
                        worldComm().barrier();
                        MPI_Finalize();
                    }

                    exit( 0 );
                }
            }
        }
    }

    // catches program_options exceptions
    catch ( boost::program_options::multiple_occurrences const& e )
    {
        LOG( WARNING ) << "Command line or config file option parsing error: " << e.what() << "\n"
                       << "  o faulty option: " << e.get_option_name() << "\n"
                       << "Warning: the .cfg file or some options may not have been read properly\n";

    }

    catch ( boost::program_options::ambiguous_option const& e )
    {
        LOG( WARNING ) << "Command line or config file option parsing error: " << e.what() << "\n"
                       << "  o faulty option: " << e.get_option_name() << "\n"
                       << "  o possible alternatives: " ;
        std::for_each( e.alternatives().begin(), e.alternatives().end(), []( std::string const& s )
        {
            LOG( WARNING ) << s << " ";
        } );
        LOG( WARNING ) << "\n"
                       << "Warning: the .cfg file or some options may not have been read properly\n";
    }

    catch ( pt::ptree_error & e )
    {
        LOG(ERROR) << "Error parsing the JSON file. Please check the JSON for syntax errors, missing commas..." << std::endl;
        LOG(ERROR) << "We suggest using 'yamllint' or 'jsonlint' available in docker or singularity images to check json files." << std::endl;
        throw;
    }

    // catches program_options exceptions
    catch ( std::exception& e )
    {
        //LOG( WARNING ) << "Application option parsing: unknown option:" << e.what() << " (the .cfg file or some options may not have been read properly)\n";
        throw;
    }

    catch ( ... )
    {
        //LOG( WARNING ) << "Application option parsing: unknown exception triggered  (the .cfg file or some options may not have been read properly)\n";
        throw;
    }
}

void
Environment::setConfigFiles( std::vector<std::string> const& cfgfiles )
{
    std::vector<fs::path> cfgAbsolutePaths;
    cfgAbsolutePaths.reserve( cfgfiles.size() );

    for ( std::string const& cfgfile : cfgfiles )
    {
        if ( cfgfile.empty() )
            continue;

        // Check if file already exists (might be from doOptions with absolute path)
        fs::path cfgPath( cfgfile );
        if ( fs::exists( cfgPath ) )
        {
            cfgAbsolutePaths.push_back( fs::absolute( cfgPath ) );
            continue;
        }

        // Try to locate the file using findFile
        std::string locatedFile = findFile( cfgfile, {} );
        if ( locatedFile.empty() )
            continue;

        fs::path cfgAbsolutePath = fs::absolute( locatedFile );
        if ( !fs::exists( cfgAbsolutePath ) )
            continue;

        cfgAbsolutePaths.push_back( cfgAbsolutePath );
    }

    if ( cfgAbsolutePaths.empty() )
        return;

    // Clear config files list but not S_vm (it may have command-line options)
    S_configFiles.clear();

    // reverse order (priority for the last)
    std::reverse( cfgAbsolutePaths.begin(), cfgAbsolutePaths.end() );

    for ( fs::path const& cfgAbsolutePath : cfgAbsolutePaths )
    {
        cout << tc::green << "Reading " << cfgAbsolutePath.string() << "..." << tc::reset << std::endl;
        S_cfgdir = cfgAbsolutePath.parent_path();
        std::ifstream ifs( cfgAbsolutePath.string().c_str() );
        std::istringstream iss( readFromFile( cfgAbsolutePath.string() ) );
        po::store( parse_config_file( ifs, *S_desc, true ), S_vm );
        S_configFiles.emplace_back( cfgAbsolutePath.string(), std::move( iss ) );
    }

    po::notify( S_vm );
}

void
Environment::setConfigFile( std::string const& cfgfile )
{
    setConfigFiles( std::vector<std::string>{ cfgfile } );
}

void
Environment::addOptions( po::options_description const& desc )
{
    if ( !S_desc )
        throw std::logic_error(
            "Environment::addOptions requires an initialized Environment" );

    bool anyRegistered = false;
    bool allRegistered = true;
    for ( auto const& option : desc.options() )
    {
        bool const registered =
            S_desc->find_nothrow( option->long_name(), false ) != nullptr;
        anyRegistered = anyRegistered || registered;
        allRegistered = allRegistered && registered;
    }

    if ( allRegistered )
        return;
    if ( anyRegistered )
        throw std::invalid_argument(
            "Environment::addOptions received a partially registered option description" );

    S_desc->add( desc );
    po::store(
        po::command_line_parser( S_argc, S_argv )
            .options( desc )
            .allow_unregistered()
            .run(),
        S_vm );

    for ( auto& config : S_configFiles )
    {
        auto& stream = std::get<1>( config );
        stream.clear();
        stream.seekg( 0 );
        po::store( po::parse_config_file( stream, desc, true ), S_vm );
    }
    po::notify( S_vm );
}

bool
Environment::initialized()
{
    return S_initialized;
}

bool Environment::aborted()
{
    return S_aborted;
}
bool
Environment::finalized()
{
    return mpi::environment::finalized();
}

mpi::threading::level
Environment::threadLevel()
{
    return mpi::environment::thread_level();
}
bool
Environment::isMainThread()
{
    return mpi::environment::is_main_thread();
}
void
Environment::abort( int error_code )
{
    S_aborted = true;
    return mpi::environment::abort( error_code );
}

/**
 * @brief Create a bootstrap directory path with MPI synchronization
 * 
 * Helper function to avoid code duplication. Creates a directory on rank 0
 * and synchronizes all ranks before returning the path.
 * 
 * @param subpath Path relative to scratchdir/bootstrap/appName
 * @return fs::path The full bootstrap path
 */
fs::path
Environment::createBootstrapPath( fs::path const& subpath )
{
    static std::map<std::string, fs::path> bootstrap_cache;
    
    std::string key = subpath.string();
    auto it = bootstrap_cache.find( key );
    if ( it != bootstrap_cache.end() )
        return it->second;
    
    fs::path bootstrap_path = scratchdir() / "bootstrap" / S_about.appName() / subpath;
    
    // Only use MPI if it's initialized and not yet finalized
    if ( initialized() && !finalized() )
    {
        if ( isMasterRank() && !fs::exists( bootstrap_path ) )
            fs::create_directories( bootstrap_path );
        worldComm().barrier();
    }
    else
    {
        // Before MPI init or after MPI finalize, just create the directories without synchronization
        if ( !fs::exists( bootstrap_path ) )
            fs::create_directories( bootstrap_path );
    }
    
    bootstrap_cache[key] = bootstrap_path;
    return bootstrap_path;
}

fs::path const&
Environment::rootRepository()
{
    if ( !repositoryConfigured() )
    {
        static fs::path bootstrap_root = createBootstrapPath( "" );
        return bootstrap_root;
    }
    return S_rootdir;
}

std::string
Environment::findFile( std::string const& filename, std::vector<std::string> paths )
{
    fs::path cp = fs::current_path();

    fs::path p( filename );

    if ( p.is_absolute() && fs::exists( p ) )
    {
        LOG( INFO ) << "File " << filename << " found";
        return filename;
    }

#if 0

    // first try in the current path
    if ( fs::exists( cp / filename ) )
    {
        LOG( INFO ) << "File " << ( cp/filename ) << " found";
        return ( cp/filename ).string();
    }

#endif

    auto filename_ = fs::path(filename).filename().string();
    for( auto const& ps : paths )
    {

        fs::path p( Environment::expand(ps) );
        if ( fs::exists( p / filename_ ) )

            return ( p / filename_ ).string();
    }
    // look in to paths list from end-1 to begin
    auto it = std::find_if( S_paths.rbegin(), S_paths.rend(),
                            [&filename] ( fs::path const& p ) -> bool
    {
        LOG(INFO) << " looking for " << p/filename << std::endl;

        if ( fs::exists( p/filename ) )
            return true;
        return false;
    } );

    if ( it != S_paths.rend() )
    {
        LOG( INFO ) << "File " << ( *it/filename ) << " found";
        return ( *it / filename ).string();
    }

    if ( fs::exists( cp / filename ) )
    {
        LOG( INFO ) << "File " << ( cp/filename ) << " found";
        return ( cp/filename ).string();
    }

    if ( fs::path( filename ).extension() == ".geo" ||
         fs::path( filename ).extension() == ".msh" ||
         fs::path( filename ).extension() == ".mesh" ||
         fs::path( filename ).extension() == ".med" )
    {
        auto filename_only = fs::path( filename ).filename();
        
        // Helper lambda to recursively search for a file in a directory
        auto search_recursive = []( fs::path const& root, fs::path const& target_filename ) -> std::string
        {
            if ( !fs::exists( root ) || !fs::is_directory( root ) )
                return std::string();
            
            try
            {
                for ( auto const& entry : fs::recursive_directory_iterator( root, fs::directory_options::follow_directory_symlink ) )
                {
                    if ( fs::is_regular_file( entry ) && entry.path().filename() == target_filename )
                    {
                        LOG( INFO ) << "File " << entry.path() << " found recursively";
                        return entry.path().string();
                    }
                }
            }
            catch ( fs::filesystem_error const& e )
            {
                LOG( WARNING ) << "Error during recursive search in " << root << ": " << e.what();
            }
            
            return std::string();
        };
        
        // First try exact path
        if ( fs::exists( fs::path( Environment::localGeoRepository() ) / filename ) )
        {
            LOG( INFO ) << "File " << ( fs::path( Environment::localGeoRepository() ) / filename ) << " found";
            return ( fs::path( Environment::localGeoRepository() ) / filename ).string();
        }

        if ( Environment::systemGeoRepository().get<1>()  &&
                fs::exists( fs::path( Environment::systemGeoRepository().get<0>() ) / filename ) )
        {
            LOG( INFO ) << "File " << ( fs::path( Environment::systemGeoRepository().get<0>() ) / filename ) << " found";
            return ( fs::path( Environment::systemGeoRepository().get<0>() ) / filename ).string();
        }
        
        // If not found, try recursive search with just the filename
        if ( filename != filename_only.string() )
        {
            // Already tried with a relative path, skip recursive search
        }
        else
        {
            // Search recursively in localGeoRepository
            if ( auto found = search_recursive( fs::path( Environment::localGeoRepository() ), filename_only ); !found.empty() )
                return found;
            
            // Search recursively in systemGeoRepository
            if ( Environment::systemGeoRepository().get<1>() )
            {
                if ( auto found = search_recursive( fs::path( Environment::systemGeoRepository().get<0>() ), filename_only ); !found.empty() )
                    return found;
            }
        }
    }

    LOG( INFO ) << "File " << filename << " not found";
    return std::string();
}
std::vector<std::string>
Environment::geoPathList()
{
    std::vector<std::string> plist;
    plist.push_back( fs::current_path().string() );
    std::for_each( S_paths.rbegin(), S_paths.rend(),
                   [&plist] ( fs::path const& p )
                   {
                       plist.push_back( p.string() );
                   } );

    if ( fs::exists( Environment::localGeoRepository() ) )
        plist.push_back( Environment::localGeoRepository() );

    if ( Environment::systemGeoRepository().get<1>()  &&
            fs::exists( Environment::systemGeoRepository().get<0>() ) )
        plist.push_back( Environment::systemGeoRepository().get<0>() );

    return plist;
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
    fs::path rep_path = Info::datadir();
    rep_path /= "geo";
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

    rep_path = Info::prefix();
    rep_path /= "share/feel/config";
    return boost::make_tuple( rep_path.string(), fs::exists( rep_path ) );
}
std::string
Environment::appRepository()
{
    if ( !repositoryConfigured() )
    {
        static fs::path bootstrap_app = createBootstrapPath( "app" );
        return bootstrap_app.string();
    }
    return S_appdir.string();
}
std::string
Environment::appRepositoryWithoutNumProc()
{
    if ( !repositoryConfigured() )
    {
        static fs::path bootstrap_app = createBootstrapPath( "app" );
        return bootstrap_app.string();
    }
    return S_appdirWithoutNumProc.string();
}
std::string
Environment::exprRepository()
{
    if ( !repositoryConfigured() )
    {
        static fs::path bootstrap_exprs = createBootstrapPath( "app/exprs" );
        return bootstrap_exprs.string();
    }
    return S_repository.exprs().string();
}

std::string
Environment::logsRepository()
{
    if ( !repositoryConfigured() )
    {
        static fs::path bootstrap_logs = createBootstrapPath( "app/logs" );
        return bootstrap_logs.string();
    }
    return (S_appdir / "logs").string();
}

std::string
Environment::exportsRepository()
{
    return (S_appdir / "exports").string();
}

std::string
Environment::downloadsRepository()
{
    if ( !repositoryConfigured() )
    {
        static fs::path bootstrap_downloads = scratchdir() / "downloads" / S_about.appName();
        if ( isMasterRank() && !fs::exists( bootstrap_downloads ) )
            fs::create_directories( bootstrap_downloads );
        worldComm().barrier();
        return bootstrap_downloads.string();
    }
    return (S_appdirWithoutNumProc / "downloads").string();
}

uuids::uuid
Environment::randomUUID( bool parallel )
{
    auto uuid = S_generator();
    if ( parallel )
    {
        // overwrite uuid with the uuid from master rank process
        std::string suuid = uuids::to_string(uuid);
        mpi::broadcast( Environment::worldComm().globalComm(), suuid, 0 );
        uuid = boost::lexical_cast<uuids::uuid>( suuid );
    }
    return uuid;
}

uuids::uuid
Environment::nameUUID( uuids::uuid const& dns_namespace_uuid, std::string const& name )
{
    boost::uuids::name_generator gen( dns_namespace_uuid );
    return gen( name );
}

void
Environment::changeRepositoryImpl( boost::format fmt, std::string const& logfilename, Location location, bool add_subdir_np, WorldComm const& worldcomm, bool remove )
{
    stopLogging( remove );
    fs::path rep_path;
    S_paths.push_back( S_appdir );
    std::string directory = fmt.str();

    // Update repository configuration from options (if available)
    if ( Environment::vm().count( "repository.append.np" ) )
        S_repository.config().append_np = boption( "repository.append.np" );
    if ( Environment::vm().count( "repository.append.date" ) )
        S_repository.config().append_date = boption( "repository.append.date" );
    
    // if we are in relative mode then first go back to the initial current path
    if ( location == Location::relative )
        ::chdir( S_paths.front().string().c_str() );
    
    // Configure repository
    // For custom locations: callback determines the path (directory parameter is ignored)
    // For other locations: use directory parameter if provided
    if ( !directory.empty() && !S_repository.isCustom() )
    {
        S_repository.configure( directory, location );
    }
    else
    {
        S_repository.configure();
    }
    S_repository.cd();

    S_rootdir = S_repository.root();
    S_appdir = S_repository.directory();
    S_appdirWithoutNumProc = S_repository.directoryWithoutAppenders();

    startLogging( Environment::about().appName() );

    Logger::console()->info("[feelpp::{}] files are stored in {}", Environment::about().appName(), Environment::appRepository());
    Logger::console()->info("[feelpp::{}] logs are stored in {}", Environment::about().appName(), Environment::logsRepository());
    Logger::console()->flush();

    // eventually cleanup
    if ( remove )
    {
        LOG(INFO) << "removing " << S_paths.back() << " after changing repository";
        if ( worldcomm.isMasterRank() )
            fs::remove_all( S_paths.back() );
    }
}

void
Environment::updateInformationObject( nl::json & p ) const
{
    if ( !p.contains( "application" ) )
    {
        p["/application/name"_json_pointer] = Environment::about().appName();

        std::string commandLineUsed;
        for ( int i = 0; i < S_argc; ++i )
        {
            commandLineUsed += S_argv[i];
            if (i < S_argc-1 )
                commandLineUsed += " ";
        }

        auto ptCfg = nl::json::array({});
        for ( auto const& cfg : S_configFiles )
            ptCfg.push_back( std::get<0>( cfg ) );

        p.emplace( "run", nl::json( {
                    { "uuid", uuids::to_string( Environment::randomUUID() ) },
                    { "directories", {
                            { "app", Environment::appRepository() },
                            { "export", Environment::exportsRepository() },
                            { "logs", Environment::logsRepository() },
                            { "exprs", Environment::exprRepository() },
                            { "downloads", Environment::downloadsRepository() }
                        }
                    },
                    { "command-line", commandLineUsed },
                    { "config-files", ptCfg },
                    { "number_of_processors", Environment::numberOfProcessors() }
                } ) );

        // Software
        p["/software/boost/version"_json_pointer] = BOOST_LIB_VERSION;
#if BOOST_VERSION >= 106700
        std::stringstream mpi_version;
        mpi_version << M_env->version().first << "."
                    << M_env->version().second;
        p["/software/mpi/version"_json_pointer] = mpi_version.str();
#else
        p["/software/mpi/version"_json_pointer] = "unknown";
#endif
#if defined( OMPI_MPI_H )
        p["/software/mpi/library/openmpi"_json_pointer] = nl::json( {
                { "version", { { "major", OMPI_MAJOR_VERSION }, { "minor", OMPI_MINOR_VERSION }, { "release", OMPI_RELEASE_VERSION } } }
            } );
#endif
#if defined ( FEELPP_HAS_PETSC_H )
        p["/software/petsc"_json_pointer] = nl::json( {
                { "version", {
                        { "major", PETSC_VERSION_MAJOR },
                        { "minor", PETSC_VERSION_MINOR },
                        { "subminor", PETSC_VERSION_SUBMINOR },
                        { "release", PETSC_VERSION_RELEASE }
                    } },
                { "date", {
                        { "release", PETSC_RELEASE_DATE },
                        { "version", PETSC_VERSION_DATE }
                    } }
            });
#endif
    }

    if ( S_hwSysInstance )
        S_hwSysInstance->updateInformationObject( p["hardware"] );
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


}

void
Environment::startLogging( std::string decorate )
{
#if defined(FEELPP_HAS_SPDLOG)
    // Determine MPI logging mode from options
    std::string log_mpi_mode = "master"; // default
    if (S_vm.count("log.mpi"))
        log_mpi_mode = S_vm["log.mpi"].as<std::string>();
    
    // Check if this rank should log based on log.mpi setting
    bool should_log = true;
    int rank = S_worldcomm->rank();
    
    if (log_mpi_mode == "none")
    {
        should_log = false;
    }
    else if (log_mpi_mode == "master")
    {
        should_log = (rank == 0);
    }
    else if (log_mpi_mode == "all")
    {
        should_log = true;
    }
    else
    {
        // Invalid value, warn and default to master
        if (rank == 0)
            std::cerr << "Warning: Invalid log.mpi value '" << log_mpi_mode 
                      << "'. Using 'master' mode. Valid values: none, master, all\n";
        should_log = (rank == 0);
    }
    
    // Use spdlog for logging
    fs::path dir = logsRepository();
    const int Nproc = 200;

    if ( S_worldcomm->size() > Nproc )
    {
        std::string smin = boost::lexical_cast<std::string>( Nproc*std::floor( S_worldcomm->rank()/Nproc ) );
        std::string smax = boost::lexical_cast<std::string>( Nproc*std::ceil( double( S_worldcomm->rank()+1 )/Nproc )-1 );
        std::string replog = smin + "-" + smax;
        dir /= replog;
    }

    // only one processor every Nproc creates the corresponding log directory
    if ( S_worldcomm->rank() % Nproc == 0 )
    {
        if ( !fs::exists( dir ) )
            fs::create_directories( dir );
    }

    // Wait for directory creation - ALL ranks must reach this barrier
    worldComm().barrier();

    // Now handle logging setup based on should_log
    if (!should_log)
    {
        // This rank should not log - set up null sink
        Logger::setDefaultLogger(Logger::createNullLogger("feelpp"));
        Logger::disable();
        Logger::verbosity() = 0;
        return;
    }

    // Build per-rank log filename
    auto file = (dir / fmt::format("rank-{:06d}.log", S_worldcomm->rank())).string();

    // Ensure the directory exists before creating file sink (safety check for NFS delays)
    if (!fs::exists(dir))
    {
        try {
            fs::create_directories(dir);
        } catch (const fs::filesystem_error& e) {
            // Directory might have been created by another rank between the check and creation
            // This is safe to ignore if the directory now exists
            if (!fs::exists(dir)) {
                std::cerr << fmt::format("[feelpp] Failed to create log directory {}: {}\n", 
                                        dir.string(), e.what());
                throw;
            }
        }
    }

    // Create file sink
    bool log_to_console = false;
    if (S_vm.count("log.console"))
        log_to_console = S_vm["log.console"].as<bool>();
    
    // Create logger with file sink, optionally add stderr sink
    if (log_to_console)
    {
        Logger::setDefaultLogger(Logger::createMultiLogger("feelpp", file, false));
    }
    else
    {
        Logger::setDefaultLogger(Logger::createFileLogger("feelpp", file, false));
    }
    Logger::setPattern("%Y-%m-%d %H:%M:%S.%e [%l] %n %P:%t %v");

    // Map Environment::logVerbosityLevel() to spdlog level and VLOG verbosity
    int v = Environment::logVerbosityLevel();
    Logger::verbosity() = v;
    Logger::setLevel(v);
#else
    // Use glog for logging
    fs::path a0 = logsRepository();
    const int Nproc = 200;

    if ( S_worldcomm->size() > Nproc )
    {
        std::string smin = boost::lexical_cast<std::string>( Nproc*std::floor( S_worldcomm->rank()/Nproc ) );
        std::string smax = boost::lexical_cast<std::string>( Nproc*std::ceil( double( S_worldcomm->rank()+1 )/Nproc )-1 );
        std::string replog = smin + "-" + smax;
        a0 /= replog;
    }

    // only one processor every Nproc creates the corresponding log directory
    if ( S_worldcomm->rank() % Nproc == 0 )
    {
        if ( !fs::exists( a0 ) )
            fs::create_directories( a0 );
    }

    FLAGS_log_dir=a0.string();

    google::AllowCommandLineReparsing();

    // duplicate argv before passing to gflags because gflags is going to rearrange them
    char** envargv = dupargv( S_argv );
    int envargc = S_argc;
    google::ParseCommandLineFlags( &envargc, &envargv/*S_argv*/, false );
    freeargv( envargv );

    // Initialize Google's logging library.
    if ( !google::glog_internal_namespace_::IsGoogleLoggingInitialized() )
    {
        if ( FLAGS_no_log )
        {
            if ( S_worldcomm->rank() == 0 && FLAGS_no_log == 1 )
                FLAGS_no_log = 0;
        }
        google::InitGoogleLogging( S_argv[0] );
    }
    google::InstallFailureSignalHandler();
#endif
}

void
Environment::stopLogging( bool remove )
{
#if defined(FEELPP_HAS_SPDLOG)
    // Shutdown spdlog
    Logger::shutdown();
    
    // Only attempt cleanup if MPI is still active and we can safely check rank
    // After PetscFinalize/SlepcFinalize, MPI may be finalized making isMasterRank() unsafe
    bool can_check_rank = initialized() && !finalized() && S_worldcomm;
    bool is_master = can_check_rank ? S_worldcomm->isMasterRank() : true;
    
    // Determine who should perform cleanup based on log.mpi mode
    bool should_cleanup = false;
    if (remove || Environment::vm().count( "rmlogs" ))
    {
        std::string log_mpi_mode = "master"; // default
        if (S_vm.count("log.mpi"))
            log_mpi_mode = S_vm["log.mpi"].as<std::string>();
        
        if (log_mpi_mode == "none")
        {
            // If no rank was logging, only master cleans up (to be safe)
            should_cleanup = is_master;
        }
        else if (log_mpi_mode == "master")
        {
            // Only master was logging, only master cleans up
            should_cleanup = is_master;
        }
        else if (log_mpi_mode == "all")
        {
            // All ranks were logging, each rank cleans up its own logs
            // But to avoid race conditions on shared directories, only master removes the whole tree
            should_cleanup = is_master;
        }
        else
        {
            // Default to master cleanup
            should_cleanup = is_master;
        }
    }
    
    if ( should_cleanup )
    {
        std::cout << tc::red << "Removing log files (--rmlogs) in " << Environment::logsRepository() << tc::reset << std::endl;
        fs::remove_all( Environment::logsRepository() );
    }
#else
    // Use glog shutdown
    if ( google::glog_internal_namespace_::IsGoogleLoggingInitialized() )
    {
        google::ShutdownGoogleLogging();
        
        // Only attempt cleanup if MPI is still active and we can safely check rank
        // After PetscFinalize/SlepcFinalize, MPI may be finalized making isMasterRank() unsafe
        bool can_check_rank = initialized() && !finalized() && S_worldcomm;
        bool is_master = can_check_rank ? S_worldcomm->isMasterRank() : true;
        
        if ( (remove || Environment::vm().count( "rmlogs" )) && is_master )
        {
            std::cout  << tc::red << "Removing log files (--rmlogs) in " << Environment::logsRepository() << tc::reset << std::endl;
            fs::remove_all( Environment::logsRepository() );
        }
    }
#endif
}

worldscomm_ptr_t &
Environment::worldsComm( int n )
{
    CHECK( S_worldcomm ) << "Environment: worldcomm not allocated\n";
    return S_worldcomm->subWorlds( n );
}

worldscomm_ptr_t &
Environment::worldsCommSeq( int n )
{
    CHECK( S_worldcommSeq ) << "Environment: worldcomm not allocated\n";
    return S_worldcommSeq->subWorlds( n );
}

worldscomm_ptr_t &
Environment::worldsCommGroupBySubspace( int n )
{
#if 0
    std::cout << "n=" << n << "\n";
    S_worldcomm->showMe();
    S_worldcomm->masterWorld( n ).showMe();
    std::cout << "size=" << S_worldcomm->subWorlds( n ).size() <<  "\n";
    S_worldcomm->subWorlds( n ).begin()->showMe();
#endif
    return S_worldcomm->subWorldsGroupBySubspace( n );
}


worldcomm_t &
Environment::masterWorldComm( int n )
{
    return S_worldcomm->masterWorld( n );
}

#if defined(FEELPP_HAS_HARTS)

void Environment::initHwlocTopology()
{
    /* init and load hwloc topology for the current node */
    if ( !( Environment::S_hwlocTopology ) )
    {
        hwloc_topology_init( &( Environment::S_hwlocTopology ) );
        hwloc_topology_load( Environment::S_hwlocTopology );
    }
}

void Environment::destroyHwlocTopology()
{
    if ( Environment::S_hwlocTopology )
    {
        hwloc_topology_destroy( Environment::S_hwlocTopology );
    }
}

void Environment::bindToCore( unsigned int id )
{
    int err;
    hwloc_cpuset_t set;
    hwloc_obj_t coren;

    /* get the nth core object */
    coren = hwloc_get_obj_by_type( Environment::S_hwlocTopology, HWLOC_OBJ_CORE, id );
    /* get the cpu mask of the nth core */
    set = hwloc_bitmap_dup( coren->cpuset );
    /* bind the process thread to this core */
    err = hwloc_set_cpubind( Environment::S_hwlocTopology, set, 0 );

    /* free memory */
    hwloc_bitmap_free( set );
}

int Environment::getNumberOfCores(bool logical)
{
    int nCores = -1;
    int depth = HWLOC_TYPE_DEPTH_UNKNOWN;
    if(logical)
    { depth = hwloc_get_type_depth( Environment::S_hwlocTopology, HWLOC_OBJ_PU ); }
    else
    { depth = hwloc_get_type_depth( Environment::S_hwlocTopology, HWLOC_OBJ_CORE ); }

    if(depth != HWLOC_TYPE_DEPTH_UNKNOWN)
    {
        nCores = hwloc_get_nbobjs_by_depth(Environment::S_hwlocTopology, depth);
    }

    return nCores;
}

int Environment::countCoresInSubtree( hwloc_obj_t node, bool logical )
{
    int res = 0;

    /* get the number of cores in the subtree */
    for ( int i = 0; i < node->arity; i++ )
    {
        res += Environment::countCoresInSubtree( node->children[i] );
    }

    /* if we are a core node, we increment the counter */
    /* count the number of real cores or logical cores */
    /* according to the logical parameter */
    if ( (logical && node->type == HWLOC_OBJ_PU)
    || (!logical && node->type == HWLOC_OBJ_CORE) )
    {
        res++;
    }

    return res;
}

void Environment::bindNumaRoundRobin( int lazy )
{
    int err, depth;
    int nbCoresPerNuma = 0, nbCoresTotal = 0, nbNumaNodesTotal = 0;
    hwloc_cpuset_t set;
    hwloc_obj_t numaNode;

    std::cout << "Round Robin Numa" << std::endl;

    /* get the first numa node */
    numaNode = hwloc_get_obj_by_type( Environment::S_hwlocTopology, HWLOC_OBJ_NODE, 0 );
    nbCoresPerNuma = Environment::countCoresInSubtree( numaNode );

    /* count the number of numaNodes */
    depth = hwloc_get_type_depth( Environment::S_hwlocTopology, HWLOC_OBJ_NODE );

    if ( depth != HWLOC_TYPE_DEPTH_UNKNOWN )
    {
        nbNumaNodesTotal = hwloc_get_nbobjs_by_depth( Environment::S_hwlocTopology, depth );
    }

    /* count the number of cores on the current server */
    depth = hwloc_get_type_depth( Environment::S_hwlocTopology, HWLOC_OBJ_CORE );

    if ( depth != HWLOC_TYPE_DEPTH_UNKNOWN )
    {
        nbCoresTotal = hwloc_get_nbobjs_by_depth( Environment::S_hwlocTopology, depth );
    }

    /* compute the virtual core index of the first core of the numa node to use */
    int vcoreid = Environment::worldComm().rank() * nbCoresPerNuma;
    /* compute the rank of the Numa processor for the current process */
    int numaRank = Environment::worldComm().rank() % nbCoresPerNuma;

    /* get the numa node where to place the process */
    numaNode = hwloc_get_obj_by_type( Environment::S_hwlocTopology, HWLOC_OBJ_NODE, numaRank );

    /* duplicate the node set of the Numa node */
    set = hwloc_bitmap_dup( numaNode->cpuset );
    /*
    char * a;
    hwloc_bitmap_asprintf(&a, set);
    std::cout << Environment::worldComm().rank() << " " << a << ";" << std::endl;
    free(a);
    */

    /* if we do not want to bind lazily, i.e. to generally bind on the numa node */
    /* we select the specific core */
    int bid = -1;

    if ( !lazy )
    {
        /* get the cpuset corresponding to the core we want to bind to */
        /* compute the core number that we want to bind to on the current Numa node */
        int tid = ( vcoreid / nbCoresTotal ) % nbCoresPerNuma;
        /* get the id of the first core */
        bid = hwloc_bitmap_first( set );

        /* iterate to find the core we want to bind to */
        for ( int i = 0; i < tid; i++ )
        {
            bid = hwloc_bitmap_next( set, bid );
        }

        hwloc_bitmap_only( set, bid );
        /*
           hwloc_bitmap_asprintf(&a, set);
           std::cout << Environment::worldComm().rank() << " " << a << ";" << std::endl;
           free(a);
           */
    }

    int coreid = vcoreid % nbCoresTotal + vcoreid / nbCoresTotal;
    std::cout << Environment::worldComm().rank() << " nbCoresNuma:" << nbCoresPerNuma << " Total:" << nbCoresTotal << " "
              << " coreid=" << coreid << " "
              << " nbCoresPerNuma=" << nbCoresPerNuma
              << " idOnNuma=" << bid
              << " numaRank=" << numaRank
              << std::endl;

    /* bind the process thread to this core */
    err = hwloc_set_cpubind( Environment::S_hwlocTopology, set, 0 );

    /* free memory */
    hwloc_bitmap_free( set );
}

void Environment::getLastBoundCPU( std::vector<int> * lastCPU, std::vector<int> * cpuAffinity )
{
    int cid;
    hwloc_cpuset_t set;

    /* get a cpuset object */
    set = hwloc_bitmap_alloc();

    if(cpuAffinity)
    {
        /* Get the cpu thread affinity info of the current process/thread */
        hwloc_get_cpubind( Environment::S_hwlocTopology, set, 0 );

        /* write the corresponding processor indexes */
        cid = hwloc_bitmap_first( set );

        while ( cid != -1 )
        {
            cpuAffinity->push_back(cid);
            cid = hwloc_bitmap_next( set, cid );
        }
    }

    hwloc_bitmap_zero(set);

    if(lastCPU)
    {
        /* Get the latest core location of the current process/thread */
        hwloc_get_last_cpu_location( Environment::S_hwlocTopology, set, 0 );

        /* write the corresponding processor indexes */
        cid = hwloc_bitmap_first( set );

        while ( cid != -1 )
        {
            lastCPU->push_back(cid);
            cid = hwloc_bitmap_next( set, cid );
        }
    }

    /* free memory */
    hwloc_bitmap_free( set );
}

void Environment::writeCPUData( std::string fname )
{
    hwloc_cpuset_t set;
    int cid;
    char * a;
    char buf[256];
    unsigned int depth;

    std::ostringstream oss;

    /* get a cpuset object */
    set = hwloc_bitmap_alloc();

    /* Get the cpu thread affinity info of the current process/thread */
    hwloc_get_cpubind( Environment::S_hwlocTopology, set, 0 );
    hwloc_bitmap_asprintf( &a, set );
    oss << a;
    free( a );

    /* write the corresponding processor indexes */
    cid = hwloc_bitmap_first( set );
    oss << " (";

    while ( cid != -1 )
    {
        oss << cid << " ";
        cid = hwloc_bitmap_next( set, cid );
    }

    oss << ")|";

    /* Get the latest core location of the current process/thread */
    hwloc_get_last_cpu_location( Environment::S_hwlocTopology, set, 0 );
    hwloc_bitmap_asprintf( &a, set );
    oss << a;
    free( a );

    /* write the corresponding processor indexes */
    cid = hwloc_bitmap_first( set );
    oss << " (";

    while ( cid != -1 )
    {
        oss << cid << " ";
        cid = hwloc_bitmap_next( set, cid );
    }

    oss << ");";

    /* free memory */
    hwloc_bitmap_free( set );

    /* if filename is empty, we write to stdout */
    if ( fname == "" )
    {
        std::cout << Environment::worldComm().rank() << " " << oss.str() << std::endl;
    }

    else
    {
        /* Write the gathered information with MPIIO */
        MPI_File fh;
        MPI_Status status;

        if ( fs::exists( fname ) )
        {
            MPI_File_delete( const_cast<char *>( fname.c_str() ), MPI_INFO_NULL );
        }

        MPI_File_open( Environment::worldComm().comm(), const_cast<char *>( fname.c_str() ), MPI_MODE_RDWR | MPI_MODE_CREATE | MPI_MODE_APPEND , MPI_INFO_NULL, &fh );
        MPI_File_write_ordered( fh, const_cast<char *>( oss.str().c_str() ), oss.str().size(), MPI_CHAR, &status );

        MPI_File_close( &fh );
    }
}

#endif

MemoryUsage
Environment::logMemoryUsage( std::string const& message )
{
    MemoryUsage mem;
#if defined ( FEELPP_HAS_PETSC_H )
    PetscMemoryGetCurrentUsage( &mem.memory_usage );
    LOG( INFO ) << message << " PETSC get current memory usage (resident memory): "  << mem.memory_usage/1e3 << "  KBytes "  << mem.memory_usage/1e6 << "  MBytes " << mem.memory_usage/1e9 << " GBytes" ;
    //PetscMemoryGetMaximumUsage( &mem );
    //LOG(INFO) << logMessage << " PETSC get maximum memory usag (resident memory): " << mem/1e6 << "  MBytes " << mem/1e9 << " GBytes" ;

    PetscMallocGetCurrentUsage( &mem.petsc_malloc_usage );
    LOG( INFO ) << message << " PETSC get current PETSC Malloc usage: "  << mem.petsc_malloc_usage/1e3 << "  KBytes " << mem.petsc_malloc_usage/1e6 << " MBytes " << mem.petsc_malloc_usage/1e9 << " GBytes" ;
    PetscMallocGetMaximumUsage( &mem.petsc_malloc_maximum_usage );
    LOG( INFO ) << message << " PETSC get maximum PETSC Malloc usage(largest memory ever used so far): "  << mem.petsc_malloc_maximum_usage/1e3 << "  KBytes " << mem.petsc_malloc_maximum_usage/1e6 << " MBytes " << mem.petsc_malloc_maximum_usage/1e9 << " GBytes" ;
#endif
    return mem;
}

std::string
Environment::expand( std::string const& expr )
{
    std::string topSrcDir = BOOST_PP_STRINGIZE( FEELPP_SOURCE_DIR );
    std::string topBuildDir = BOOST_PP_STRINGIZE( FEELPP_BUILD_DIR );
    std::string cfgDir = S_cfgdir.string();
    std::string homeDir = ::getenv( "HOME" );
    
    // Prefer build data directory if it exists (for testing without install)
    std::string buildDataDir = topBuildDir + "/share/feelpp/data";
    std::string installDataDir = BOOST_PP_STRINGIZE( FEELPP_DATADIR );
    std::string dataDir = fs::exists(fs::path(buildDataDir)) ? buildDataDir : installDataDir;
    
    std::string exprdbDir = ( fs::path( Environment::rootRepository() )/fs::path( "exprDB" ) ).string();

    VLOG( 2 ) << "topSrcDir " << topSrcDir << "\n"
              << "topBuildDir " << topBuildDir << "\n"
              << "cfgDir " << cfgDir << "\n"
              << "HOME " << homeDir << "\n"
              << "Environment::rootRepository() " << Environment::rootRepository()
              << "dataDir " << dataDir << "\n"
              << "exprdbdir " << exprdbDir << "\n"
              << "\n";

    std::string res=expr;

    boost::replace_all( res, "${feelpp_srcdir}", topSrcDir );
    boost::replace_all( res, "${feelpp_builddir}", topBuildDir );
    boost::replace_all( res, "${feelpp_databasesdir}", topSrcDir + "/databases/" );
    boost::replace_all( res, "${top_srcdir}", topSrcDir );
    boost::replace_all( res, "${toolboxes_srcdir}", topSrcDir + "/toolboxes/" );
    boost::replace_all( res, "${top_builddir}", topBuildDir );
    boost::replace_all( res, "${cfgdir}", cfgDir );
    boost::replace_all( res, "${home}", homeDir );
    boost::replace_all( res, "${repository}", Environment::rootRepository().string() );
    boost::replace_all( res, "${appdir}", Environment::appRepository() );
    boost::replace_all( res, "${datadir}", dataDir );
    boost::replace_all( res, "${exprdbdir}", exprdbDir );
    boost::replace_all( res, "${h}", std::to_string(doption(_name="gmsh.hsize") ) );

    boost::replace_all( res, "$feelpp_srcdir", topSrcDir );
    boost::replace_all( res, "$feelpp_builddir", topBuildDir );
    boost::replace_all( res, "$feelpp_databasesdir", topSrcDir + "/databases/" );
    boost::replace_all( res, "$top_srcdir", topSrcDir );
    boost::replace_all( res, "$toolboxes_srcdir", topSrcDir + "/toolboxes/" );
    boost::replace_all( res, "$top_builddir", topBuildDir );
    boost::replace_all( res, "$cfgdir", cfgDir );
    boost::replace_all( res, "$home", homeDir );
    boost::replace_all( res, "$repository", Environment::rootRepository().string() );
    boost::replace_all( res, "$appdir", Environment::appRepository() );
    boost::replace_all( res, "$datadir", dataDir );
    boost::replace_all( res, "$exprdbdir", exprdbDir );
    boost::replace_all( res, "$h", std::to_string(doption(_name="gmsh.hsize") ) );
    boost::replace_all( res, "$np", std::to_string(Environment::numberOfProcessors()) );

    typedef std::vector< std::string > split_vector_type;

#if defined(FEELPP_ENABLED_PROJECTS)
    split_vector_type SplitVec; // #2: Search for tokens
    boost::split( SplitVec, FEELPP_ENABLED_PROJECTS, boost::is_any_of(" "), boost::token_compress_on );
    for( auto const& s : SplitVec )
    {
        std::ostringstream oo1,oo2,oo3;
        oo1 << "${" << s << "_srcdir}";
        oo2 << "${" << s << "_builddir}";
        oo3 << "${" << s << "_databasesdir}";

        boost::replace_all( res, oo1.str(), topSrcDir + "/research/" + s );
        VLOG(2) << oo1.str() << " : " << topSrcDir + "/research/" + s;
        boost::replace_all( res, oo2.str(),  topBuildDir + "/research/" + s );
        VLOG(2) << oo2.str() << " : " << topBuildDir + "/research/" + s;
        boost::replace_all( res, oo3.str(),  topSrcDir + "/research/" + s + "/databases/" );
        VLOG(2) << oo3.str() << " : " << topSrcDir + "/research/" + s + "/databases/";;

        std::ostringstream o1,o2,o3;
        o1 << "$" << s << "_srcdir";
        o2 << "$" << s << "_builddir";
        o3 << "$" << s << "_databasesdir";

        boost::replace_all( res, o1.str(), topSrcDir + "/research/" + s );
        VLOG(2) << o1.str() << " : " << topSrcDir + "/research/" + s;
        boost::replace_all( res, o2.str(),  topBuildDir + "/research/" + s );
        VLOG(2) << o2.str() << " : " << topBuildDir + "/research/" + s;
        boost::replace_all( res, o3.str(),  topSrcDir + "/research/" + s + "/databases/" );
        VLOG(2) << o3.str() << " : " << topSrcDir + "/research/" + s + "/databases/";;
    }
#endif

    VLOG(1) << "Expand " << expr << " to "  << res;
    return res;
}

//std::unique_ptr<TimerTable>
//Environment::timers()
//{
//    return S_timers;
//}

void
Environment::addTimer( std::string const& msg,
                       std::pair<double,int> const& t,
                       std::string const& uiname = "" )
{
    S_timers->add( msg, t, uiname );
}

void
Environment::saveTimers( bool display )
{
    //S_timers.save( Environment::about().appName(), display );
    S_timers->save( display );
}

void
Environment::saveTimersMD( std::ostream &os )
{
    //S_timers.save( Environment::about().appName(), display );
    S_timers->saveMD( os );
}

int Environment::S_argc = 0;
char** Environment::S_argv = 0;

AboutData Environment::S_about;
std::shared_ptr<po::command_line_parser> Environment::S_commandLineParser;
std::vector<std::tuple<std::string,std::istringstream> > Environment::S_configFiles;
po::variables_map Environment::S_vm;
std::shared_ptr<po::options_description> Environment::S_desc;
std::shared_ptr<po::options_description> Environment::S_desc_app;
std::shared_ptr<po::options_description> Environment::S_desc_lib;
std::vector<std::string> Environment::S_to_pass_further;

boost::signals2::signal<void()> Environment::S_deleteObservers;

std::string Environment::S_log_mpi_mode = "master";
std::shared_ptr<WorldComm> Environment::S_worldcomm;
std::shared_ptr<WorldComm> Environment::S_worldcommSeq;
boost::uuids::random_generator Environment::S_generator;

std::vector<fs::path> Environment::S_paths = { fs::current_path(),
                                               Environment::systemConfigRepository().get<0>(),
                                               Environment::systemGeoRepository().get<0>()
                                             };
fs::path Environment::S_rootdir = fs::current_path();
fs::path Environment::S_appdir = fs::current_path();
fs::path Environment::S_appdirWithoutNumProc;
fs::path Environment::S_scratchdir;
fs::path Environment::S_cfgdir;

std::string Environment::olAppPath;
std::vector<std::string> Environment::olAutoloadFiles;

#if defined(FEELPP_HAS_HARTS)
hwloc_topology_t Environment::S_hwlocTopology = NULL;
#endif

std::unique_ptr<TimerTable> Environment::S_timers;
std::unique_ptr<Sys::HwSysBase> Environment::S_hwSysInstance;

std::unique_ptr<JournalWatcher> Environment::S_informationObject;
}
