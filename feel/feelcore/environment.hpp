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
#ifndef __Environment_H
#define __Environment_H 1

#include <cstdlib>

#include <boost/noncopyable.hpp>
#include <boost/signals2.hpp>
#include <boost/format.hpp>



#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/parameter.hpp>
#include <feel/feelcore/worldcomm.hpp>
#include <feel/feelcore/about.hpp>
#include <feel/options.hpp>

namespace Feel
{
namespace detail
{
inline
AboutData
makeAbout( char* name )
{
    AboutData about( name,
                     name,
                     "0.1",
                     name,
                     AboutData::License_GPL,
                     "Copyright (c) 2012 Feel++ Consortium" );

    about.addAuthor( "Feel++ Consortium",
                     "",
                     "feelpp-devel@feelpp.org", "" );
    return about;
}

/** @brief Initialize, finalize, and query the Feel++ environment.
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

    template <class ArgumentPack>
    Environment(ArgumentPack const& args)
        {
            char** argv = args[_argv];
            int argc = args[_argc];
            S_desc = boost::shared_ptr<po::options_description>( new po::options_description( args[_desc | Feel::feel_options()] ) );
            AboutData about = args[_about| makeAbout(argv[0])];
            S_desc->add( file_options( about.appName() ) );

            init( argc, argv, *S_desc, about );
        }

    void init( int argc, char** argv, po::options_description const& desc, AboutData const& about );

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
    static WorldComm& worldComm() { return *S_worldcomm; }

    /**
     * return n sub world communicators
     */
    static std::vector<WorldComm> const&  worldsComm( int n );

    static std::vector<WorldComm> const&  worldsCommGroupBySubspace( int n );

    /**
     * return master world comm associated with a color map of size n
     */
    static WorldComm const& masterWorldComm( int n );

    /**
     * return number of processors
     */
    static int numberOfProcessors()  { return S_worldcomm->godSize(); }

    /**
     * return variables_map
     */
    static po::variables_map const& vm() { return S_vm; }

    static AboutData const& about() { return S_about; }

    /**
     * return options description data structure
     */
    static po::options_description const& optionsDescription() { return *S_desc; }

    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * set the static worldcomm
     */
    static void setWorldComm( WorldComm& worldcomm ) { S_worldcomm = worldcomm.shared_from_this(); }



    //@}

    /** @name  Methods
     */
    //@{

    //! \return the root repository (default: \c $HOME/feel)
    static std::string rootRepository();

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
        (void), static changeRepository, tag,
        (required
         (directory,(boost::format)))
        (optional
         (filename,*( boost::is_convertible<mpl::_,std::string> ),"logfile")
         (subdir,*( boost::is_convertible<mpl::_,bool> ),true)
            ))
        {
            changeRepositoryImpl( directory, filename, subdir );
        }

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
            S_deleteObservers.connect(boost::bind(&Observer::operator(), obs));
        }

    /**
     * \return the scratch directory
     */
    static const fs::path& scratchDirectory() { return S_scratchdir; }

    //@}


private:

    //! change the directory where the results are stored
    static void changeRepositoryImpl( boost::format fmt, std::string const& logfile, bool add_subdir_np );

    //! process command-line/config-file options
    static void doOptions( int argc, char** argv, po::options_description const& desc, std::string const& appName );
    static void processGenericOptions();
    static void parseAndStoreOptions( po::command_line_parser parser, bool extra_parser = false );

private:
    /// Whether this environment object called MPI_Init
    bool i_initialized;
    mpi::environment M_env;

    static  fs::path S_scratchdir;

    static AboutData S_about;
    static po::variables_map S_vm;
    static boost::shared_ptr<po::options_description> S_desc;
    static std::vector<std::string> S_to_pass_further;

    static boost::signals2::signal<void ()> S_deleteObservers;

    static boost::shared_ptr<WorldComm> S_worldcomm;
};
} // detail



class Environment : public detail::Environment
{
public:
    BOOST_PARAMETER_CONSTRUCTOR(
        Environment, (detail::Environment), tag,
        (required
         (argc,*)
         (argv,*))
        (optional
         (desc,*)
         (about,*) )) // no semicolon
};


}
#endif /* __Environment_H */
