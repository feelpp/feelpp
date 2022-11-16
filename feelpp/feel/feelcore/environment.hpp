//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//! 
//! This file is part of the Feel library
//! 
//! Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! Date: 2010-04-14
//! 
//! Copyright (C) 2010-2012 Université Joseph Fourier (Grenoble I)
//! 
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//! 
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//! 
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//! 
//! 
//! \file environment.hpp
//! \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! \date 2010-04-14
//! 
#ifndef FEELPP_ENVIRONMENT_HPP
#define FEELPP_ENVIRONMENT_HPP 1

#include <cstdlib>
#include <memory>

#include <fmt/core.h>
#include <fmt/format.h>

#include <boost/noncopyable.hpp>
#include <boost/signals2/signal.hpp>

#include <boost/stacktrace.hpp>
#include <boost/exception/all.hpp>

#include <boost/format.hpp>
#include <boost/property_tree/ptree_fwd.hpp>

#include <feel/feelcore/feel.hpp>


#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/uuid/uuid_io.hpp>         // streaming operators etc.
#include <boost/uuid/uuid_serialize.hpp>         // streaming operators etc.

#include <feel/feelcore/parameter.hpp>
#include <feel/feelcore/worldcomm.hpp>
#include <feel/feelcore/worldscomm.hpp>

#include <feel/feelcore/rank.hpp>
#include <feel/feelcore/about.hpp>
#include <feel/feelcore/termcolor.hpp>
#include <feel/options.hpp>
#include <feel/feelcore/repository.hpp>

#if defined ( FEELPP_HAS_PETSC_H )
#include <petscsys.h>
#endif

#if defined(FEELPP_HAS_HARTS)
#include <hwloc.h>
#endif
#include <feel/feelcore/repository.hpp>
#include <feel/feelcore/journalmanager.hpp>
#include <feel/feelhwsys/hwsys.hpp>

namespace pybind11
{
class list;
}
namespace Feel
{
namespace tc = termcolor;
namespace pt =  boost::property_tree;
namespace uuids =  boost::uuids;

// boost::error_info typedef that holds the stacktrace:
using traced = boost::error_info<struct tag_stacktrace, boost::stacktrace::stacktrace>;


/**
 * @brief helper class for throwing any exception with stacktrace:
 *
 * @tparam E exception type
 * @param e exception thown
 */
template <class E>
void throw_with_trace( const E& e )
{
    throw boost::enable_error_info( e )
        << traced( boost::stacktrace::stacktrace() );
}

// forward declarationx@
class TimerTable;

//!
//! @class MemoryUsage
//! @ingroup Core
//! @brief class to query for memory usage
struct FEELPP_EXPORT MemoryUsage
{
    MemoryUsage()
        :
        memory_usage(0)
#if defined ( FEELPP_HAS_PETSC_H )
        , petsc_malloc_usage(0)
        , petsc_malloc_maximum_usage(0)
#endif
        {}
    MemoryUsage(MemoryUsage const& m )
        :
        memory_usage(m.memory_usage)
#if defined ( FEELPP_HAS_PETSC_H )
        , petsc_malloc_usage(m.petsc_malloc_usage)
        , petsc_malloc_maximum_usage(m.petsc_malloc_maximum_usage)
#endif
        {}
    MemoryUsage& operator=(MemoryUsage const& m )
        {
            if ( this != &m )
            {
                memory_usage = m.memory_usage;
#if defined ( FEELPP_HAS_PETSC_H )
                petsc_malloc_usage = m.petsc_malloc_usage;
                petsc_malloc_maximum_usage = m.petsc_malloc_maximum_usage;
#endif
            }
            return *this;
        }
#if defined ( FEELPP_HAS_PETSC_H )
    PetscLogDouble memory_usage;
    PetscLogDouble petsc_malloc_usage;
    PetscLogDouble petsc_malloc_maximum_usage;
#else
    double memory_usage;
#endif

};
//! @ingroup Core
//! default \c makeAbout function to define the \c AboutData structure of the Feel++
//! application
//! @param name name or short name of the application
//! 
FEELPP_EXPORT AboutData makeAboutDefault( std::string name );

//! 
//! @class Environment "Environment"
//! @ingroup Core
//! @brief Initialize, finalize, and query the Feel++ environment.
//! @ingroup Core
//! 
//! The @c Environment class is used to initialize, finalize, and
//! query the Feel++ environment. It will typically be used in the @c
//! main() function of a program, which will create a single instance
//! of @c Environment initialized with the arguments passed to the
//! program:
//! 
//! @code
//! int main(int argc, char* argv[])
//! {
//!    using namespace Feel;
//!    Environment env(argc, argv);
//! }
//! @endcode
//! or more commonly
//! @code
//! int main(int argc, char* argv[])
//! {
//!    using namespace Feel;
//!    po::options_description myoptions( "Laplacian options" );
//!    myoptions.add_options()
//!          ( "mu", po::value<double>()->default_value( 1.0 ), "coeff" );
//!    Environment env( _argc=argc, _argv=argv,
//!                     _desc=myoptions,
//!                    _about=about(_name="myapp",
//!                                 _author="Feel++ Consortium",
//!                                 _email="feelpp-devel@feelpp.org"));
//! }
//! @endcode
//!
//! The instance of @c Environment will initialize Feel++ (by calling @c MPI, @c
//! PETSc, @c SLEPc or Logging initialization routines) in its constructor
//! and finalize in its destructor.
//! 
//! @author Christophe Prud'homme
//! @see Application
//! 
class FEELPP_EXPORT Environment
:   boost::noncopyable,
    public JournalManager
{
public:
    //!
    //! @name Constants
    //!
    //! @{


    //! @}

    /** @name Typedefs
     */
    //@{
    typedef WorldComm worldcomm_type;
    typedef std::shared_ptr<WorldComm> worldcomm_ptrtype;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    /** Initialize the Feel environment.
     *
     *  If the Feel environment has not already been initialized,
     *  initializes Feel
     */
    Environment();

    /** Initialize the Feel environment.
     *
     *  If the Feel environment has not already been initialized,
     *  initializes Feel
     *
     *  @param argc The number of arguments provided in @p argv, as
     *  passed into the program's @c main function.
     *
     *  @param argv The array of argument strings passed to the program
     *  via @c main.
     *
     */
    Environment( int& argc, char** &argv );

    /**
     * @brief Construct a new Environment object
     * 
     * @param argc number of command line aguments
     * @param argv command line aguments
     * @param lvl level of threading in MPI 
     * @param desc command line options for application
     * @param desc_lib command line options for library
     * @param about about data structure
     * @param config configuration of the repository of the resut&lts
     */
    Environment( int argc, char** argv,
                 mpi::threading::level lvl,
                 po::options_description const& desc,
                 po::options_description const& desc_lib,
                 AboutData const& about,
                 Repository::Config const& config );

#if defined(FEELPP_ENABLE_PYTHON_WRAPPING)
    /**
     * @brief Construct a new Environment object
     * 
     * @param arg sys arg
     * @param desc description of the options
     * @param directory directory to save the results
     * @param chdir change directory to FEELPP_REPOSITORY or stay in current directory
     */
    Environment( pybind11::list arg, po::options_description const& desc, Repository::Config const& config );
    Environment( pybind11::list arg );
#endif

    Environment( NA::arguments<
                 typename na::argc::template required_as_t<int>,
                 typename na::argv::template required_as_t<char**>,
                 typename na::threading::template required_as_t<mpi::threading::level>,
                 typename na::desc::template required_as_t<po::options_description const&>,
                 typename na::desc_lib::template required_as_t<po::options_description const&>,
                 typename na::about::template required_as_t<AboutData const&>,
                 typename na::config::template required_as_t<Repository::Config const&>
                 > && args )
        :
        Environment( args.get(_argc), args.get(_argv), args.get(_threading),
                     args.get(_desc), args.get(_desc_lib),
                     args.get(_about), args.get(_config)
                     ) {}
  private :
    template <typename ... Ts>
        decltype(auto) make_args_env_constructor( Ts && ... v )
    {
        auto args0 = NA::make_arguments( std::forward<Ts>(v)... )
            .add_default_arguments( NA::make_default_argument( _threading, mpi::threading::funneled ),
                                    NA::make_default_argument_invocable( _desc, [](){ return feel_nooptions(); } ),
                                    NA::make_default_argument_invocable( _desc_lib, [](){ return feel_options(); } )
                                    );
        char** argv = args0.get(_argv);
        return std::move(args0).add_default_arguments( NA::make_default_argument_invocable( _about, [&argv](){ return makeAboutDefault( argv[0] ); } ),
                                                       NA::make_default_argument_invocable( _config, [&argv](){ return globalRepository(makeAboutDefault( argv[0] ).appName()); } )
                                                       );
    }
  public:

    template <typename ... Ts,typename  = typename std::enable_if_t< sizeof...(Ts) != 0 && ( NA::is_named_argument_v<Ts> && ...) > >
        Environment( Ts && ... v )
            :
            Environment( make_args_env_constructor( std::forward<Ts>(v)... ) )
        {}

    /** Shuts down the Feel environment.
     *
     *  If this @c Environment object was used to initialize the Feel
     *  environment, and the Feel environment has not already been shut
     *  down (finalized), this destructor will shut down the Feel
     *  environment.
     */
    ~Environment() override;

    //@}

    /** @name Operator overloads
     */
    //@{

    //@}

    /** @name Accessors
     */
    //@{

    /** Determine if the MPI environment has already been initialized.
     *
     *  This routine is equivalent to a call to @c MPI_Initialized.
     *
     *  @returns @c true if the MPI environment has been initialized.
     */
    static bool initialized();

    /** Determine if the MPI environment has already been finalized.
     *
     *  The routine is equivalent to a call to @c MPI_Finalized.
     *
     *  @returns @c true if the MPI environment has been finalized.
     */
    static bool finalized();

    /** Determine if the MPI environment has already been aborted.
     *
     *  this occurs if MPI_Abort has been called
     *
     *  @returns @c true if the Feel++ nvironment has been aborted.
     */
    static bool aborted();

    /**
     * @return the shared_ptr WorldComm
     */
    static std::shared_ptr<WorldComm> const& worldCommPtr()
        {
            return S_worldcomm;
        }

    /**
     * return the worldcomm (static)
     */
    static WorldComm& worldComm()
        {
            return *S_worldcomm;
        }
    static WorldComm& worldCommSeq()
    {
        return *S_worldcommSeq;
    }
    static worldcomm_ptr_t& worldCommSeqPtr()
    {
        return S_worldcommSeq;
    }

    /**
     * return n sub world communicators
     */
    static worldscomm_ptr_t &  worldsComm( int n );
    static worldscomm_ptr_t &  worldsCommSeq( int n );

    static worldscomm_ptr_t &  worldsCommGroupBySubspace( int n );

    /**
     * return master world comm associated with a color map of size n
     */
    static worldcomm_t & masterWorldComm( int n );

    /**
     * return number of processors
     */
    static int numberOfProcessors()
    {
        return S_worldcomm->godSize();
    }

    /**
     * return the rank in global mpi communicator
     */
    static rank_type rank()
    {
        return S_worldcomm->globalRank();
    }

    /**
     * return master rank in mpi communicator
     */
    static rank_type masterRank()
    {
        return S_worldcomm->masterRank();
    }

    /** get the thread level that could be set
     * The thread level resquested to mpi may not be supported.
     * This function returns the maximum level of thread support that could be setup
     * @return the maximum thread level supported with respect to the initialization request
     */
    static mpi::threading::level threadLevel();

    /** Are we in the main thread?
     * this function may be useful e.g. in funneled level
     */
    static bool isMainThread();

    /** Abort all MPI processes.
     *  Aborts all MPI processes and returns to the environment. The
     *  precise behavior will be defined by the underlying MPI
     *  implementation. This is equivalent to a call to @c MPI_Abort
     *  with @c MPI_COMM_WORLD.
     *
     *  @param errcode The error code to return to the environment.
     *  @returns Will not return.
     */
    static void abort(int errcode);

    /**
     * @return true if number of process is 1, hence the environment is
     * sequential
     */
    static bool isSequential()
    {
        return numberOfProcessors() == 1;
    }

    /**
     * @return true if the environment is not sequential, that is if the number
     * of process is greater than 1
     */
    static bool isParallel()
    {
        return !isSequential();
    }
    /**
     * rank 0 process is considered the master process
     *
     * the master process can then for example print information in the console
     * or in some files
     */
    static bool isMasterRank()
    {
        return rank() == masterRank();
    }

    static po::command_line_parser const& commandLineParser()
    {
        return *S_commandLineParser;
    }

    static std::vector<std::tuple<std::string,std::istringstream> > & configFiles()
    {
        for ( auto & configFile : S_configFiles )
        {
            std::istringstream & iss = std::get<1>( configFile );
            iss.clear();
            iss.seekg(0, std::ios::beg);
        }
        return S_configFiles;
    }

    /**
     * @brief Set the Configuration from a File 
     * 
     * @param filename filename of the config file
     */
    static void setConfigFile( std::string const& filename );

    /**
     * return variables_map
     */
    static po::variables_map const& vm()
    {
        return S_vm;
    }

    template<typename T>
    static void setOptionValue(std::string s,T val)
    {
        auto it = S_vm.find( s );
        CHECK( it != S_vm.end() ) << "Invalid option " << s << "\n";
        S_vm.at(s).value() = val;
    }

    static AboutData const& about()
    {
        return S_about;
    }

    /**
     * Adds a file to automatically load in Gmsh with Onelab
     */
    static void olLoadInGmsh( std::string filename )
    {
        olAutoloadFiles.push_back( filename );
    }

    /**
     * return options description data structure
     */
    static po::options_description const& optionsDescription()
    {
        return *S_desc;
    }

    /**
     * return options description data structure
     */
    static po::options_description const& optionsDescriptionApplication()
    {
        return *S_desc_app;
    }

    /**
     * return the options description for the Feel++ library
     */
    static po::options_description const& optionsDescriptionLibrary()
    {
        return *S_desc_lib;
    }

    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * @brief get the repository info
     * 
     * @return Repository& 
     */
    static Repository& repository() { return S_repository; }

    /**
     * set the static worldcomm
     */
    static void setWorldComm( WorldComm& worldcomm )
    {
        S_worldcomm = worldcomm.shared_from_this();
    }

#if defined(FEELPP_HAS_HARTS)

    /**
     * Init Hwloc topology structure
     */
    static void initHwlocTopology();

    /**
     * Destroy Hwloc topology structure
     */
    static void destroyHwlocTopology();

    static hwloc_topology_t getHwlocTopology()
    {
        return Environment::S_hwlocTopology;
    }

    /**
     * Binds the current process/thread to the specified core.
     * (Do it early in the application launch, otherwise you might end up accessing data "far away"
     * from the new bound core, thus degrading performance)
     */
    static void bindToCore( unsigned int id );

    /**
     * Counts the number of cores on the current server
     * Calls countCoresInSubtree done on the whole topology
     *
     *  @param logical boolean indicating if we want to include logical cores, i.e. hyperthreading
     */
    static int getNumberOfCores( bool logical = false );

    /**
     * Counts the number of cores under the current hwloc object, using a recursive strategy
     *
     *  @param logical boolean indicating if we want to include logical cores, i.e. hyperthreading
     */
    static int countCoresInSubtree( hwloc_obj_t node, bool logical = false );

    /**
     * Binds the MPI processes in Round Robin on the NUMA nodes
     */
    static void bindNumaRoundRobin( int lazy = false );

    /**
     * Get information about the last CPU bound. You must use --bind-to core with MPI for this feature to work.
     */
    static void getLastBoundCPU( std::vector<int> * cpuAffinity, std::vector<int> * lastCPU );

    /**
     * Writes data about processor affinity and last location of the different processes/threads
     * (last location is not guaranteed to be right, unles you bind the process to a core)
     */
    static void writeCPUData( std::string fname = "CPUData.dat" );

#endif

    //@}

    /** @name  Methods
     */
    //@{

    //! @name directories
    //! @{
    template <typename ... Ts>
    static void changeRepository( Ts && ... v )
    {
        auto args = NA::make_arguments( std::forward<Ts>(v)... );
        boost::format directory = args.get(_directory);
        std::string location_str = args.get_else( _location, "global" );
        std::string const& filename = args.get_else(_filename, "logfile" );
        bool subdir = args.get_else_invocable(_subdir, [](){ return S_vm["npdir"].as<bool>(); } );
        WorldComm const& worldcomm = args.get_else(_worldcomm, Environment::worldComm() );
        bool remove = args.get_else(_remove, false );
        changeRepositoryImpl( directory, filename, location(location_str), subdir, worldcomm, remove );
    }

    //! \return the root repository (default: \c $HOME/feel)
    static std::string const& rootRepository();

    /**
     * Find a file. The lookup is as follows:
     *  - look into current path
     *  - look into paths that went through changeRepository(), it means that we look
     *    for example into the path from which the executable was run
     * If the file has an extension .geo or .msh, try also to
     *  - look into \c localGeoRepository() which is usually $HOME/feel/geo
     *  - look into \c systemGeoRepository() which is usually $FEELPP_DIR/share/feel/geo
     * If \p filename is not found, then the empty string is returned.
     * \return the string containing the filename path
     */
    static std::string findFile( std::string const& filename, std::vector<std::string> paths = {} );

    /**
     * \return the list of paths where Feel++ looks into to find a Gmsh Geo file
     */
    static std::vector<std::string> geoPathList();

    //! \return the local geo files repository (default: \c $HOME/feel/geo)
    static std::string localGeoRepository();

    /**
     * \return a tuple : the system geo files repository (default: \c
     * /usr/share/feel/geo or /usr/local/share/feel/geo) and true or false
     * whether the directory exists or not
     */
    static boost::tuple<std::string,bool> systemGeoRepository();


    //! \return the local config files repository (default: \c $HOME/feel/config)
    static std::string localConfigRepository();

    /**
     * \return a tuple : the system config files repository (default: \c
     * /usr/share/feel/config or /usr/local/share/feel/config) and true or false
     * whether the directory exists or not
     */
    static boost::tuple<std::string,bool> systemConfigRepository();

    //! the application directory, application results are stored there
    //! the directory is controlled by changeRepository
    static std::string appRepository();

    //! the application directory without np, application results are stored there
    //! the directory is controlled by changeRepository
    static std::string appRepositoryWithoutNumProc();

    //! the expressions repository is typically a sub-directory of the \c
    //! appRepository() that contains the expressions generated using Ginac
    static std::string exprRepository();

    //! the logfiles repository is a subdirectory of the \c appRepository containing the logfiles
    static std::string logsRepository();

    //! the exports repository is a subdirectory of the \c appRepository
    //! containing the results exported during the application execution
    static std::string exportsRepository();

    //! the downloads repository is a subdirectory of the \c appRepository
    //! containing the files downloaded during the application execution
    static std::string downloadsRepository();

    //!
    //! Generate a random UUID
    //!
    //! the UUID is very very very very likely to be unique as it is encoded into 128 bits
    //! @param parallel if true generate the same uuid for all MPI process, if false it will be different
    //!
    static uuids::uuid randomUUID( bool parallel=true );

    //! Generate a UUID from a namespace UUID and a name.
    static uuids::uuid nameUUID( uuids::uuid const& dns_namespace_uuid, std::string const& name );

    //! @}
    template <typename ... Ts>
    static std::pair<std::string,po::variable_value> option( Ts && ... v )
    {
        auto args = NA::make_arguments( std::forward<Ts>(v)... );
        std::string const& name  = args.get(_name);
        std::string const& sub = args.get_else(_sub, "" );
        std::string const& prefix = args.get_else(_prefix, "" );
        po::variables_map const& vm = args.get_else(_vm, Environment::vm() );

        auto opt = fmt::memory_buffer();

        if ( !prefix.empty() )
            fmt::format_to( opt, "{}.",prefix);
        if ( !sub.empty() )
            fmt::format_to( std::back_inserter( opt ), "{}-",sub);
        fmt::format_to( std::back_inserter( opt ), "{}",name);
        std::string optname = fmt::to_string(opt);
        auto it = vm.find( optname );
        if ( it == vm.end() )
            throw_with_trace( std::invalid_argument( fmt::format( "{}:{} invalid or missing option {}", __FILE__, __LINE__, optname ) ) );
        ///CHECK( it != vm.end() ) << "Invalid option " << os.str() << "\n";
        return *it;
    }

    template <typename ... Ts>
    static po::variable_value vm( Ts && ... v )
    {
        return option( std::forward<Ts>(v)... ).second;
    }

    /**
     * print resident memory usage as well as PETSc malloc usage in log file
     * \param message message to print to identity the associated memory operation
     */
    static MemoryUsage logMemoryUsage( std::string const& message );

    //! Return the timerstable pointer.
//    static std::unique_ptr<TimerTable> timers();

    /**
     * add timer to a map of timers that can be shown using \c displayTimers()
     */
    static void addTimer( std::string const& msg,
                          std::pair<double,int> const& t,
                          std::string const& uiname );

    /**
     * display and save timers
     */
    static void saveTimers( bool save );
    static void saveTimersMD( std::ostream & os );

    //! get  \c variables_map from \c options_description \p desc
    //static po::variables_map vm( po::options_description const& desc );

    //! get the log verbosity level
    static int logVerbosityLevel() { return FLAGS_v; }

    //! set the verbosity level of the VLOG macro
    static void setLogVerbosityLevel( int logVerbosity );

    //! 
    //!  set log files
    //!  \param prefix prefix for log filenames
    FEELPP_DEPRECATED static void setLogs( std::string const& prefix );

    //!
    //! start logging.
    //! \code
    //! Environment::startLogging();
    //! LOG(INFO) << "Feel++ uses logging";
    //! Environment::stopLogging();
    //! \endcode
    //!
    static void startLogging( std::string decorate );

    //!
    //! stop logging in Feel++.
    //! \param remove deletes the log directory and all its content
    //! \code
    //! Environment::startLogging();
    //! LOG(INFO) << "Feel++ uses logging";
    //! Environment::stopLogging();
    //! \endcode
    //! 
    //! \code
    //! Environment::startLogging();
    //! LOG(INFO) << "Feel++ uses logging";
    //! Environment::stopLogging( true ); // delete the `logs` subdirectory of the current path
    //! \endcode
    //!
    static void stopLogging( bool remove = false );

    template<typename Observer>
    static void
    addDeleteObserver( Observer const& obs )
    {
        S_deleteObservers.connect( obs );
    }
    template<typename Observer>
    static void
    addDeleteObserver( std::shared_ptr<Observer> const& obs )
    {
        S_deleteObservers.connect( boost::bind( &Observer::operator(), obs ) );
    }

    static void clearSomeMemory();

    //!
    //!  \return the scratch directory
    //!
    FEELPP_DEPRECATED static const fs::path& scratchDirectory()
    {
        return S_scratchdir;
    }

    /**
     * @brief expand feel++ pathes in a string
     * @details Feel++ defines some paths
     *  - top_srcdir : top source directory from which feel++ was compiled
     *  - top_builddir : top build directory in which feel++ was compiled
     *  - repository : repository for the results (default: $HOME/feel)
     *  - datadir : repository for data like meshes ($top_srcdir/data)
     *  - exprdbdir : top level directory of ginac expression ($repository/exprDB)
     *  - home : $HOME directory
     * The paths can be used in filename string and expanded using this function by prefixing them by $
     * @code
     * std::string topsrcdir = Environment::expand( "$top_srcdir");
     * std::string topbuilddir = Environment::expand( "$top_builddir");
     * @endcode
     * They can be used in config files or in the command line
     * @code
     * feelpp_qs_laplacian --gmsh.filename=\$home/mesh.geo
     * @endcode
     *
     * @param expr string where paths might be expanded
     * @return the expansion of the feel++ paths defined in string expr
     */
    static std::string expand( std::string const& expr );

    /**
     * try find remotely the file \p fname 
     * \param fname filename 
     * \param subdir ubdirectory to store the file that may be downloaded
     * @return the filename
     */ 
    static std::string findFileRemotely( std::string const& fname, std::string const& subdir = "" ); 
    //@}

private:

    //! Private Methods
    //! @{

    //! change the directory where the results are stored
    static void changeRepositoryImpl( boost::format fmt, std::string const& logfile, Location location, bool add_subdir_np, WorldComm const& worldcomm, bool remove );

#if defined ( FEELPP_HAS_PETSC_H )
    FEELPP_NO_EXPORT void initPetsc( int * argc = 0, char *** argv = NULL );
#endif

    //! process command-line/config-file options
    static FEELPP_NO_EXPORT void doOptions( int argc, char** argv,
                                            po::options_description const& desc,
                                            po::options_description const& desc_lib,
                                            std::string const& appName );

    /**
     * \fn void generateOLFiles( int argc, char ** argv, std::string const& appName )
     * \brief Generate configuration files for interaction with Gmsh through OneLab.
     * \author Carolina Diaz, Jérôme Boeglin, Sébastien Landré
     *
     * @param argc Number of application arguments.
     * @param argv Application arguments.
     * @param appName Name of the application.
     */
    static FEELPP_NO_EXPORT void generateOLFiles( int argc, char ** argv, std::string const& appName );
    static FEELPP_NO_EXPORT void processGenericOptions();
    static FEELPP_NO_EXPORT void parseAndStoreOptions( po::command_line_parser parser, bool extra_parser = false );

    //! update information into ptree
    void updateInformationObject( nl::json & p ) const;

    //! @}

private:
    /// Whether this environment object called MPI_Init
    std::unique_ptr<mpi::environment> M_env;

    //! number of arguments in command line
    static int S_argc;
    //! arguments in command line
    static char** S_argv;

    static std::vector<fs::path> S_paths;

    inline static Repository S_repository = unknownRepository();
    static fs::path S_rootdir;
    static fs::path S_appdir;
    static fs::path S_appdirWithoutNumProc;
    static fs::path S_scratchdir;
    static fs::path S_cfgdir;
    static AboutData S_about;
    static inline bool S_initialized = false;
    static inline bool S_aborted = false;
    static inline bool S_init_python = true;
    static std::shared_ptr<po::command_line_parser> S_commandLineParser;
    static std::vector<std::tuple<std::string,std::istringstream> > S_configFiles;
    static po::variables_map S_vm;
    static std::shared_ptr<po::options_description> S_desc;
    static std::shared_ptr<po::options_description> S_desc_app;
    static std::shared_ptr<po::options_description> S_desc_lib;
    static std::vector<std::string> S_to_pass_further;

    static uuids::random_generator S_generator;

    /**
     * Stores the absolute path and executable name
     */
    static std::string olAppPath;

    /**
     * Stores names of output files for automatic loading in Gmsh with Onelab
     */
    static std::vector<std::string> olAutoloadFiles;

    static boost::signals2::signal<void()> S_deleteObservers;

    static std::shared_ptr<WorldComm> S_worldcomm;
    static std::shared_ptr<WorldComm> S_worldcommSeq;

#if defined(FEELPP_HAS_HARTS)
    static hwloc_topology_t S_hwlocTopology;
#endif
    static std::unique_ptr<TimerTable> S_timers;

    //! Hardware System information instance.
    static std::unique_ptr<Sys::HwSysBase> S_hwSysInstance;

    static std::unique_ptr<JournalWatcher> S_informationObject;
};


template <typename ... Ts>
po::variable_value option( Ts && ... v )
{
    return Environment::vm( std::forward<Ts>(v)... );
}

template <typename ... Ts>
int countoption( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    po::variables_map const& vm = args.get_else(_vm, Environment::vm() );
    try
    {
        return vm.count( Environment::option( std::forward<Ts>(v)... ).first );
    }
    catch( std::invalid_argument const& e )
    {
        return 0;
    }
}


template <typename T,typename ... Ts,typename = typename std::enable_if_t< sizeof...(Ts) != 0 && ( NA::is_named_argument_v<Ts> && ...) > >
T optionT( Ts && ... v )
{
    try
    {
        return Environment::vm( std::forward<Ts>(v)... ).template as<T>();
    }
    catch ( boost::bad_any_cast const& bac )
    {
        auto args = NA::make_arguments( std::forward<Ts>(v)... );
        std::string const& name  = args.get(_name);
        //CHECK( false ) <<"Option "<< name << "  either does not exist or is not a double" <<std::endl;
        CHECK( false ) <<"problem in conversion type of argument "<< name << " : check the option type"<<std::endl;
        return {};
    }
}

template <typename T>
T optionT( std::string const& name )
{
    return optionT<T>(_name=name);
}

template <typename ... Ts>
double doption( Ts && ... v )
{
    return optionT<double>( std::forward<Ts>(v)... );
}
template <typename ... Ts>
bool boption( Ts && ... v )
{
    return optionT<bool>( std::forward<Ts>(v)... );
}
template <typename ... Ts>
int ioption( Ts && ... v )
{
    return optionT<int>( std::forward<Ts>(v)... );
}
template <typename ... Ts>
std::string soption( Ts && ... v )
{
    return optionT<std::string>( std::forward<Ts>(v)... );
}
template <typename ... Ts>
std::vector<std::string> vsoption( Ts && ... v )
{
    return optionT<std::vector<std::string>>( std::forward<Ts>(v)... );
}
template <typename ... Ts>
std::vector<double> vdoption( Ts && ... v )
{
    return optionT<std::vector<double>>( std::forward<Ts>(v)... );
}

/**
 * @brief handle exceptions using a Lippincott function
 * this funnction allows to handle the exceptions thrown in Feel++
 * \code
 * try {
 *   Environment env(),
 *   // Feel++ code here
 * }
 * catch( ... )
 * {
 *   handleExceptions()
 * }
 * \endcode
 */
void handleExceptions();


} // Feel

#include <feel/feelcore/feelio.hpp>
#endif /* FEELPP_ENVIRONMENT_HPP */
