/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-04-14

  Copyright (C) 2010-2012 Université Joseph Fourier (Grenoble I)

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
   \file environment.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-04-14
 */
#ifndef FEELPP_ENVIRONMENT_HPP
#define FEELPP_ENVIRONMENT_HPP 1

#include <cstdlib>
#include <memory>

#include <boost/noncopyable.hpp>
#include <boost/signals2.hpp>
#include <boost/format.hpp>

#include <feel/feelcore/feel.hpp>

#if defined(FEELPP_HAS_BOOST_PYTHON) && defined(FEELPP_ENABLE_PYTHON_WRAPPING)
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>

//#include <mpi4py/mpi4py.h>
#endif

#include <feel/feelcore/parameter.hpp>
#include <feel/feelcore/worldcomm.hpp>
#include <feel/feelcore/worldscomm.hpp>
#include <feel/feelcore/about.hpp>
#include <feel/options.hpp>
#if defined ( FEELPP_HAS_PETSC_H )
#include <petscsys.h>
#endif

#if defined(FEELPP_HAS_HARTS)
#include <hwloc.h>
#endif

namespace Feel
{
struct MemoryUsage
{
    MemoryUsage()
        :
#if defined ( FEELPP_HAS_PETSC_H )
        memory_usage(0),
        petsc_malloc_usage(0),
        petsc_malloc_maximum_usage(0)
#endif
        {}
    MemoryUsage(MemoryUsage const& m )
        :
#if defined ( FEELPP_HAS_PETSC_H )
        memory_usage(m.memory_usage),
        petsc_malloc_usage(m.petsc_malloc_usage),
        petsc_malloc_maximum_usage(m.petsc_malloc_maximum_usage)
#endif
        {}
    MemoryUsage& operator=(MemoryUsage const& m )
        {
            if ( this != &m )
            {
#if defined ( FEELPP_HAS_PETSC_H )
                memory_usage = m.memory_usage;
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
#endif
    
};
/**
 * default \c makeAbout function to define the \c AboutData structure of the Feel++
 * application
 * @param name name or short name of the application
 */
AboutData makeAboutDefault( std::string name );

/**
 *  @class Environment "Environment"
 *  @brief Initialize, finalize, and query the Feel++ environment.
 *
 *  The @c Environment class is used to initialize, finalize, and
 *  query the Feel++ environment. It will typically be used in the @c
 *  main() function of a program, which will create a single instance
 *  of @c Environment initialized with the arguments passed to the
 *  program:
 *
 *  @code
 *  int main(int argc, char* argv[])
 *  {
 *    Feel::Environment env(argc, argv);
 *  }
 *  @endcode
 *
 *  The instance of @c Environment will initialize Feel++ (by calling @c MPI, @c
 *  PETSc, @c SLEPc and @c MAdLib initialization routines) in its constructor
 *  and finalize in its destructor.
 *
 * @author Christophe Prud'homme
 * @see Application
 */
class Environment : boost::noncopyable
{
public:


    /** @name Constants
     */
    //@{


    //@}

    /** @name Typedefs
     */
    //@{
    typedef WorldComm worldcomm_type;
    typedef boost::shared_ptr<WorldComm> worldcomm_ptrtype;

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

    Environment( int argc, char** argv,
#if BOOST_VERSION >= 105500
                 mpi::threading::level lvl,
#endif
                 po::options_description const& desc,
                 po::options_description const& desc_lib,
                 AboutData const& about,
                 std::string directory );
    
#if defined(FEELPP_HAS_BOOST_PYTHON) && defined(FEELPP_ENABLE_PYTHON_WRAPPING)
    Environment( boost::python::list arg );
#endif


    template <class ArgumentPack>
    Environment( ArgumentPack const& args )
        :
        Environment( args[_argc],
                     args[_argv],
#if BOOST_VERSION >= 105500                     
                     args[_threading|mpi::threading::single],
#endif
                     args[_desc|feel_nooptions()],
                     args[_desc_lib | feel_options()],
                     args[_about| makeAboutDefault( args[_argv][0] )],
                     args[_directory|args[_about| makeAboutDefault( args[_argv][0] )].appName()] )
        {}
#if BOOST_VERSION >= 105500                     
    BOOST_PARAMETER_CONSTRUCTOR(
        Environment, ( Environment ), tag,
        ( required
          ( argc,* )
          ( argv,* ) )
        ( optional
          ( desc,* )
          ( desc_lib,* )
          ( about,* )
          ( threading,(mpi::threading::level) )
          ( directory,( std::string ) )
          ) ) // no semicolon
#else
    BOOST_PARAMETER_CONSTRUCTOR(
        Environment, ( Environment ), tag,
        ( required
          ( argc,* )
          ( argv,* ) )
        ( optional
          ( desc,* )
          ( desc_lib,* )
          ( about,* )
          ( directory,( std::string ) )
          ) ) // no semicolon
#endif
    
    /** Shuts down the Feel environment.
     *
     *  If this @c Environment object was used to initialize the Feel
     *  environment, and the Feel environment has not already been shut
     *  down (finalized), this destructor will shut down the Feel
     *  environment.
     */
    ~Environment();

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

    /**
     * return n sub world communicators
     */
    static std::vector<WorldComm> const&  worldsComm( int n );
    static std::vector<WorldComm> const&  worldsCommSeq( int n );

    static std::vector<WorldComm> const&  worldsCommGroupBySubspace( int n );

    /**
     * return master world comm associated with a color map of size n
     */
    static WorldComm const& masterWorldComm( int n );

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
     * rank 0 process is considered the master process
     *
     * the master process can then for example print information in the console
     * or in some files
     */
    static bool isMasterRank()
    {
        return rank() == 0;
    }

    static po::command_line_parser const& commandLineParser()
    {
        return *S_commandLineParser;
    }
    static std::set<std::string> configFileNames()
    {
        return S_configFileNames;
    }

    /**
     * return variables_map
     */
    static po::variables_map const& vm()
    {
        return S_vm;
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
     * Counts the number of cores under the current hwloc object, using a recursive strategy
     */
    static int countCoresInSubtree( hwloc_obj_t node );

    /**
     * Binds the MPI processes in Round Robin on the NUMA nodes
     */
    static void bindNumaRoundRobin( int lazy = false );

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

    BOOST_PARAMETER_MEMBER_FUNCTION(
        ( void ), static changeRepository, tag,
        ( required
          ( directory,( boost::format ) ) )
        ( optional
          ( filename,*( boost::is_convertible<mpl::_,std::string> ),"logfile" )
          ( subdir,*( boost::is_convertible<mpl::_,bool> ),S_vm["npdir"].as<bool>() )
          ( worldcomm, ( WorldComm ), Environment::worldComm() )
          ) )
        {
            changeRepositoryImpl( directory, filename, subdir, worldcomm );
        }

    //! \return the root repository (default: \c $HOME/feel)
    static std::string rootRepository();

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
    static std::string findFile( std::string const& filename );

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

    BOOST_PARAMETER_MEMBER_FUNCTION(
        ( po::variable_value ), static vm, tag,
        ( required
          ( name,( std::string ) ) )
        ( optional
          ( worldcomm, ( WorldComm ), Environment::worldComm() )
          ( sub,( std::string ),"" )
          ( prefix,( std::string ),"" )
          ( vm, ( po::variables_map const& ), Environment::vm() )
        ) )
    {
        std::ostringstream os;

        if ( !prefix.empty() )
            os << prefix << ".";

        if ( !sub.empty() )
            os << sub << "-";

        os << name;
        auto it = vm.find( os.str() );
        CHECK( it != vm.end() ) << "Invalid option " << os.str() << "\n";
        return it->second;
    }

    /**
     * print resident memory usage as well as PETSc malloc usage in log file
     * \param message message to print to identity the associated memory operation
     */
    static MemoryUsage logMemoryUsage( std::string const& message );

    //! get  \c variables_map from \c options_description \p desc
    //static po::variables_map vm( po::options_description const& desc );

    /**
     * set log files
     * \param prefix prefix for log filenames
     */
    static void setLogs( std::string const& prefix );

    template<typename Observer>
    static void
    addDeleteObserver( Observer const& obs )
    {
        S_deleteObservers.connect( obs );
    }
    template<typename Observer>
    static void
    addDeleteObserver( boost::shared_ptr<Observer> const& obs )
    {
        S_deleteObservers.connect( boost::bind( &Observer::operator(), obs ) );
    }

    static void clearSomeMemory();

    /**
     * \return the scratch directory
     */
    static const fs::path& scratchDirectory()
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

    //@}


private:

    //! change the directory where the results are stored
    static void changeRepositoryImpl( boost::format fmt, std::string const& logfile, bool add_subdir_np, WorldComm const& worldcomm );

#if defined ( FEELPP_HAS_PETSC_H )
    void initPetsc( int * argc = 0, char *** argv = NULL );
#endif

    

    //! process command-line/config-file options
    static void doOptions( int argc, char** argv,
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
    static void generateOLFiles( int argc, char ** argv, std::string const& appName );
    static void processGenericOptions();
    static void parseAndStoreOptions( po::command_line_parser parser, bool extra_parser = false );

private:
    /// Whether this environment object called MPI_Init
    bool i_initialized;
    std::unique_ptr<mpi::environment> M_env;

    static std::vector<fs::path> S_paths;

    static  fs::path S_scratchdir;

    static AboutData S_about;
    static boost::shared_ptr<po::command_line_parser> S_commandLineParser;
    static std::set<std::string> S_configFileNames;
    static po::variables_map S_vm;
    static boost::shared_ptr<po::options_description> S_desc;
    static boost::shared_ptr<po::options_description> S_desc_app;
    static boost::shared_ptr<po::options_description> S_desc_lib;
    static std::vector<std::string> S_to_pass_further;

    /**
     * Stores the absolute path and executable name
     */
    static std::string olAppPath;

    /**
     * Stores names of output files for automatic loading in Gmsh with Onelab
     */
    static std::vector<std::string> olAutoloadFiles;

    static boost::signals2::signal<void()> S_deleteObservers;

    static boost::shared_ptr<WorldComm> S_worldcomm;
    static boost::shared_ptr<WorldComm> S_worldcommSeq;

#if defined(FEELPP_HAS_HARTS)
    static hwloc_topology_t S_hwlocTopology;
#endif
};

BOOST_PARAMETER_FUNCTION(
    ( po::variable_value ), option, tag,
    ( required
      ( name,( std::string ) ) )
    ( optional
      ( worldcomm, ( WorldComm ), Environment::worldComm() )
      ( sub,( std::string ),"" )
      ( prefix,( std::string ),"" )
      ( vm, ( po::variables_map const& ), Environment::vm() )
    ) )
{
    return Environment::vm( _name=name,_worldcomm=worldcomm,_sub=sub,_prefix=prefix, _vm=vm );
}

BOOST_PARAMETER_FUNCTION(
    ( double ),
    doption, tag,
    ( required
      ( name,( std::string ) ) )
    ( optional
      ( worldcomm, ( WorldComm ), Environment::worldComm() )
      ( sub,( std::string ),"" )
      ( prefix,( std::string ),"" )
    ) )
{
    double opt;

    try
    {
        opt = Environment::vm( _name=name,_worldcomm=worldcomm,_sub=sub,_prefix=prefix ).template as<double>();
    }

    catch ( boost::bad_any_cast bac )
    {
        CHECK( false ) <<"Option "<< name << "  either does not exist or is not a double" <<std::endl;
    }

    return opt;
}

BOOST_PARAMETER_FUNCTION(
    ( bool ),
    boption, tag,
    ( required
      ( name,( std::string ) ) )
    ( optional
      ( worldcomm, ( WorldComm ), Environment::worldComm() )
      ( sub,( std::string ),"" )
      ( prefix,( std::string ),"" )
    ) )
{
    bool opt;

    try
    {
        opt = Environment::vm( _name=name,_worldcomm=worldcomm,_sub=sub,_prefix=prefix ).template as<bool>();
    }

    catch ( boost::bad_any_cast bac )
    {
        CHECK( false ) <<"Option "<< name << "  either does not exist or is not a boolean" <<std::endl;
    }

    return opt;
}

BOOST_PARAMETER_FUNCTION(
    ( int ),
    ioption, tag,
    ( required
      ( name,( std::string ) ) )
    ( optional
      ( worldcomm, ( WorldComm ), Environment::worldComm() )
      ( sub,( std::string ),"" )
      ( prefix,( std::string ),"" )
    ) )
{
    int opt;

    try
    {
        opt = Environment::vm( _name=name,_worldcomm=worldcomm,_sub=sub,_prefix=prefix ).template as<int>();
    }

    catch ( boost::bad_any_cast bac )
    {
        CHECK( false ) <<"Option "<< name << "  either does not exist or is not an integer" <<std::endl;
    }

    return opt;
}


BOOST_PARAMETER_FUNCTION(
    ( std::string ),
    soption, tag,
    ( required
      ( name,( std::string ) ) )
    ( optional
      ( worldcomm, ( WorldComm ), Environment::worldComm() )
      ( sub,( std::string ),"" )
      ( prefix,( std::string ),"" )
    ) )
{
    std::string opt;

    try
    {
        opt = Environment::vm( _name=name,_worldcomm=worldcomm,_sub=sub,_prefix=prefix ).template as<std::string>();
    }

    catch ( boost::bad_any_cast bac )
    {
        CHECK( false ) <<"Option "<< name << "  either does not exist or is not a string" <<std::endl;
    }

    return opt;
}

BOOST_PARAMETER_FUNCTION(
    ( std::vector<std::string> ),
    vsoption, tag,
    ( required
      ( name,( std::string ) ) )
    ( optional
      ( worldcomm, ( WorldComm ), Environment::worldComm() )
      ( sub,( std::string ),"" )
      ( prefix,( std::string ),"" )
    ) )
{
	std::vector<std::string> opt;

    try
    {
        opt = Environment::vm( _name=name,_worldcomm=worldcomm,_sub=sub,_prefix=prefix ).template as<std::vector<std::string>>();
    }

    catch ( boost::bad_any_cast bac )
    {
        CHECK( false ) <<"Option "<< name << "  either does not exist or is not a string" <<std::endl;
    }

    return opt;
}

namespace detail
{
template<typename Args, typename Tag=tag::opt>
struct option
{
    typedef typename boost::remove_pointer<
    typename boost::remove_const<
    typename boost::remove_reference<
    typename parameter::binding<Args, Tag>::type
    >::type
    >::type
    >::type type;
};

}

BOOST_PARAMETER_FUNCTION(
    ( typename Feel::detail::option<Args>::type ),
    optionT, tag,
    ( required
      ( name,( std::string ) )
      ( in_out( opt ),* ) )
    ( optional
      ( worldcomm, ( WorldComm ), Environment::worldComm() )
      ( sub,( std::string ),"" )
      ( prefix,( std::string ),"" )
    ) )
{
    try
    {
        opt = Environment::vm( _name=name,_worldcomm=worldcomm,_sub=sub,_prefix=prefix ).template as<typename Feel::detail::option<Args>::type>();
    }

    catch ( boost::bad_any_cast bac )
    {
        CHECK( false ) <<"problem in conversion type of argument "<< name << " : check the option type"<<std::endl;
    }

    return opt;
}

}
#endif /* FEELPP_ENVIRONMENT_HPP */
