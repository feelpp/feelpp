/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-03-17

  Copyright (C) 2005,2006,2009 EPFL
  Copyright (C) 2007,2008 Université Joseph Fourier (Grenoble I)

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
   \file application.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-03-17
 */
#ifndef __application_H
#define __application_H 1

#include <boost/optional.hpp>
#include <boost/format.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifecore/about.hpp>

#include <iostream>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/string.hpp> // Needed to send/receive strings!

#include <boost/mpi.hpp>
#if defined(HAVE_MPI_H)
#include <mpi.h>
#endif /* HAVE_MPI_H */

#if defined(HAVE_TAU)
#include <Profile/Profiler.h>
#endif /* HAVE_TAU */


namespace Life
{
namespace mpi = boost::mpi;
/**
 * \class Application
 *\ingroup Core
 *\brief provides information about the Application
 *
 * @author Christophe Prud'homme
 */
class Application
{
public:


    /** @name Typedefs
     */
    //@{


    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * Construct an MPI Application
     *
     * @param argc number of arguments on the command line
     * @param argv arguments in the command line
     * @param ad \p AboutData structure for this \p Application
     * @param Comm MPI communicator
     */
#if defined( HAVE_MPI )
    Application( int argc, char** argv, AboutData const& ad, MPI_Comm Comm = MPI_COMM_WORLD );
#else
    Application( int argc, char** argv, AboutData const& ad );
#endif

    /**
     * Construct an MPI Application
     *
     * @param argc number of arguments on the command line
     * @param argv arguments in the command line
     * @param ad \p AboutData structure for this \p Application
     * @param od \p po::options_description structure for this \p Application
     * @param Comm MPI communicator
     */
#if defined( HAVE_MPI )
    Application( int argc, char** argv, AboutData const& ad, po::options_description const& od, MPI_Comm Comm = MPI_COMM_WORLD );
#else
    Application( int argc, char** argv, AboutData const& ad, po::options_description const& od );
#endif

    /**
     * copy constructor
     * @param app \p Application to be copy constructed
     */
    Application( Application const & app );

    virtual ~Application();

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /**
     * get the options description
     *
     *
     * @return the options description
     */
    po::options_description const& optionsDescription() const { return _M_desc; }

    /**
     * get the variable map
     *
     *
     * @return the variable map
     */
    po::variables_map const& vm() const { return _M_vm; }

    /**
     * get the about data of the application
     *
     *
     * @return the about data ofthe application
     * @see AboutData
     */
    AboutData const& about() const { return _M_about; }

    /**
     * \return the number of options unrecognized by \p boost::program_options
     */
    int unknownArgc() const { return _M_to_pass_further.size()+1; }

    /**
     * \return the \c char** of unrecognized options
     */
    char** unknownArgv() const;

    /**
     * \return the number of processes
     */
    static uint16_type nProcess() { return uint16_type(_S_n_process); }

    /**
     * \return the id of the current process
     */
    static uint16_type processId() { return uint16_type(_S_process_id); }

     /**
     * \return \p true if MPI is initialized, \p false otherwise
     */
    bool isMPIInitialized() const { return _S_is_mpi_initialized; }

    /** Determine if the MPI environment has already been initialized.
     *
     *  This routine is equivalent to a call to @c MPI_Initialized.
     *
     *  @returns @c true if the MPI environment has been initialized.
     */
    static bool initialized() { return mpi::environment::initialized(); }

    /** Determine if the MPI environment has already been finalized.
     *
     *  The routine is equivalent to a call to @c MPI_Finalized.
     *
     *  @returns @c true if the MPI environment has been finalized.
     */
    static bool finalized() { return mpi::environment::finalized(); }

    /** Retrieve the name of this processor.
     *
     *  This routine returns the name of this processor. The actual form
     *  of the name is unspecified, but may be documented by the
     *  underlying MPI implementation. This routine is implemented as a
     *  call to @c MPI_Get_processor_name.
     *
     *  @returns the name of this processor.
     */
    static std::string processorName()
    {
#if defined( HAVE_MPI )
        return mpi::environment::processor_name();
#else
        // fallback
        return std::string( "localhost" );
#endif
    }

    //! \return the root of life applications (typically $HOME/life)
    std::string rootRepository() const;


    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * name1 represents the first level name
     */
    void setName1( std::string const& name1 );

    /**
     * name2 represents the second level name
     */
    void setName2( std::string const& name2 );

    /**
     * h is the mesh size
     */
    void setH( double h, int precision = 4 );

    /**
     * set the dimension of the problem
     */
    void setDimension( int dim );


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * change to Simulation Repository
     */
    Application& changeRepository( boost::format );

#if defined( HAVE_MPI )
    static MPI_Comm COMM_WORLD;
#endif // HAVE_MPI

    /**
     * @return the communicator
     */
    static mpi::communicator const& comm() { return S_world; }

    /**
     * @return the barrier
     */
    static void barrier() { S_world.barrier(); }

    //@}



protected:

    /**
     * parse and store application options from cmdline
     * @param argc number of arguments
     * @param argv arguments
     */
    void doOptions( int argc, char** argv );

    /**
     * \internal
     * process the generic options passed to the command line
     *
     */
    void processGenericOptions();


    /**
     * \internal
     * parse and store option in \c po::variables_map
     * \param parser the type of parse to be used
     * \param extra_parser \c true if use extra parser for response file, \c false otherwise
     */
    void parseAndStoreOptions( po::command_line_parser parser, bool extra_parser = false );

protected:
    void setLogs();
protected:
    static int _S_n_process;

    static int _S_process_id;

private:

    void initMPI( int, char**, MPI_Comm );
    void initPETSc();
    void initTrilinos();

private:

    AboutData _M_about;

    po::options_description _M_desc;
    po::variables_map _M_vm;

    boost::optional<std::string> _M_name1;
    boost::optional<std::string> _M_name2;
    boost::optional<std::pair<double, int> > _M_h;
    boost::optional<int> _M_dim;

    std::vector<std::string> _M_to_pass_further;


    static bool _S_is_mpi_initialized;
    boost::shared_ptr<mpi::environment> M_env;
    static mpi::communicator S_world;
};

}
#endif /* __Application_H */
