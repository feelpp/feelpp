/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-03-17
 */
#ifndef __application_H
#define __application_H 1

#include <boost/optional.hpp>
#include <boost/format.hpp>
#include <boost/ptr_container/ptr_list.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/about.hpp>
#include <feel/feelcore/simget.hpp>

#include <iostream>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/serialization/string.hpp> // Needed to send/receive strings!

#include <boost/mpi.hpp>
#if defined(FEELPP_HAS_MPI_H)
#include <mpi.h>
#endif /* FEELPP_HAS_MPI_H */

#if defined(FEELPP_HAS_TAU)
#include <Profile/Profiler.h>
#endif /* FEELPP_HAS_TAU */


namespace Feel
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

    //! Simget collection type
    typedef boost::ptr_list<Simget> simgets_type;

    //! Simget iterator over the collection
    typedef simgets_type::iterator  simget_iterator;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    Application();

    /**
     * Construct an MPI Application
     *
     * @param ad \p AboutData structure for this \p Application
     * @param Comm MPI communicator
     */
#if defined( FEELPP_HAS_MPI )
    Application( AboutData const& ad, MPI_Comm Comm = MPI_COMM_WORLD );
#else
    Application( AboutData const& ad );
#endif

    /**
     * Construct an MPI Application
     *
     * @param ad \p AboutData structure for this \p Application
     * @param od \p po::options_description structure for this \p Application
     * @param Comm MPI communicator
     */
#if defined( FEELPP_HAS_MPI )
    Application( AboutData const& ad, po::options_description const& od, MPI_Comm Comm = MPI_COMM_WORLD );
#else
    Application( AboutData const& ad, po::options_description const& od );
#endif

    /**
     * Construct an MPI Application
     *
     * @param argc number of arguments on the command line
     * @param argv arguments in the command line
     * @param ad \p AboutData structure for this \p Application
     * @param Comm MPI communicator
     *
     * \warning this function is marked deprecated and should not be used
     * anymore, use the Environment class instead to initialize the Feel++
     */
#if defined( FEELPP_HAS_MPI )
    Application( int argc, char** argv, AboutData const& ad, MPI_Comm Comm = MPI_COMM_WORLD ) FEELPP_DEPRECATED;
#else
    Application( int argc, char** argv, AboutData const& ad ) FEELPP_DEPRECATED;
#endif

    /**
     * Construct an MPI Application
     *
     * @param argc number of arguments on the command line
     * @param argv arguments in the command line
     * @param ad \p AboutData structure for this \p Application
     * @param od \p po::options_description structure for this \p Application
     * @param Comm MPI communicator
     *
     * \warning this function is marked deprecated and should not be used
     * anymore, use the Environment class instead to initialize the Feel++
     */
#if defined( FEELPP_HAS_MPI )
    Application( int argc, char** argv, AboutData const& ad, po::options_description const& od, MPI_Comm Comm = MPI_COMM_WORLD ) FEELPP_DEPRECATED;
#else
    Application( int argc, char** argv, AboutData const& ad, po::options_description const& od ) FEELPP_DEPRECATED;
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
    po::options_description const& optionsDescription() const
    {
        return M_desc;
    }

    /**
     * get the variable map
     *
     *
     * @return the variable map
     */
    po::variables_map const& vm() const
    {
        return M_vm;
    }

    /**
     * get the about data of the application
     *
     *
     * @return the about data ofthe application
     * @see AboutData
     */
    AboutData const& about() const
    {
        return M_about;
    }

    /**
     * \return the number of options unrecognized by \p boost::program_options
     */
    int unknownArgc() const
    {
        return M_to_pass_further.size()+1;
    }

    /**
     * \return the \c char** of unrecognized options
     */
    char** unknownArgv() const;

    /**
     * \return the number of processes
     */
    uint16_type nProcess()
    {
        return uint16_type( M_comm.size() );
    }

    /**
     * \return the id of the current process
     */
    uint16_type processId()
    {
        return uint16_type( M_comm.rank() );
    }

    /**
    * \return \p true if MPI is initialized, \p false otherwise
    */
    bool isMPIInitialized() const
    {
        return mpi::environment::initialized();
    }

    /** Determine if the MPI environment has already been initialized.
     *
     *  This routine is equivalent to a call to @c MPI_Initialized.
     *
     *  @returns @c true if the MPI environment has been initialized.
     */
    static bool initialized()
    {
        return mpi::environment::initialized();
    }

    /** Determine if the MPI environment has already been finalized.
     *
     *  The routine is equivalent to a call to @c MPI_Finalized.
     *
     *  @returns @c true if the MPI environment has been finalized.
     */
    static bool finalized()
    {
        return mpi::environment::finalized();
    }

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
#if defined( FEELPP_HAS_MPI )
        return mpi::environment::processor_name();
#else
        // fallback
        return std::string( "localhost" );
#endif
    }

    //! \return the root of feel applications (typically $HOME/feel)
    std::string rootRepository() const;

    //! \return the \c begin() iterator
    simget_iterator begin()
    {
        return M_simgets.begin();
    }

    //! \return the \c end() iterator
    simget_iterator end()
    {
        return M_simgets.end();
    }

    //! \return the number of simgets
    size_type nSimgets() const
    {
        return M_simgets.size();
    }

    /**
     * \return true if the verbose command line/config option is used, false
     * otherwise
     */
    bool verbose() const
    {
        return M_vm.count( "verbose" );
    }

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

#if defined( FEELPP_HAS_MPI )
    static MPI_Comm COMM_WORLD;
#endif // FEELPP_HAS_MPI

    /**
     * @return the communicator
     */
    WorldComm& comm()
    {
        return M_comm;
    }

    /**
     * @return the communicator
     */
    WorldComm const& comm() const
    {
        return M_comm;
    }

    /**
     * @return the barrier
     */
    void barrier()
    {
        M_comm.barrier();
    }

    /**
     * add a new simget to the application
     */
    void add( Simget* simget );

    /**
     * execute the set of Simget stored in the Application
     */
    virtual void run();

    /**
     * execute the set of Simget stored in the Application following the
     * input/output model \f$ Y=F(X) \f$. \f$ P\f$ is the number of inputs and
     * \f$ N\f$ the number of outputs. Denote \f$ S \f$ (\c nSimgets()) the
     * number of simgets stored in the Application. \f$ X \f$ and \f$ Y\f$ must
     * be of size \f$ S P\f$ and \f$ S N \f$ respectively.
     */
    virtual void run( const double* X, unsigned long P, double* Y, unsigned long N );

    enum Stats
    {
        FLAT    = 1<<1,
        HEADER  = 1<<2,
        ERRORS  = 1<<3,
        TIME    = 1<<4,
        DATA    = 1<<5,
        NUMBERS = 1<<6,
        ALL     = ERRORS | TIME | DATA | NUMBERS
    };

    /**
     * set statistics to be printed
     */
    void setStats( std::vector<std::string> const& keys );

    /**
     * store statistics \p s in statistics map entry with key \p n
     */
    void storeStats( std::string const&  n, ptree::ptree const& s );

    /**
     * print statistics from applications
     */
    void printStats( std::ostream& out, size_type stats = ALL ) const;


    /**
     * print statistics from applications
     */
    void printStats( std::ostream& out, std::vector<std::string> const& keys, size_type stats = ALL ) const;

    /**
     * result file name
     */
    std::string resultFileName() const;

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

    /**
     * set log files
     * \param prefix prefix for log filenames
     */
    void setLogs();

private:

    void initMPI( int, char**, MPI_Comm );
    void initPETSc();
    void initTrilinos();

private:

    AboutData M_about;

    po::options_description M_desc;
    po::variables_map M_vm;

    boost::optional<std::string> M_name1;
    boost::optional<std::string> M_name2;
    boost::optional<std::pair<double, int> > M_h;
    boost::optional<int> M_dim;

    std::vector<std::string> M_to_pass_further;


    boost::shared_ptr<mpi::environment> M_env;
    WorldComm M_comm;

    simgets_type M_simgets;
    std::map<std::string,std::vector<ptree::ptree> > M_stats;
    std::vector<std::string> M_keys;

};
}
#endif /* __Application_H */
