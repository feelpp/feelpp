/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-10-18

  Copyright (C) 2005,2006,2009 EPFL
  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \date 2005-10-18
 */
#ifndef __Application_H
#define __Application_H 1

#include <iostream>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/string.hpp> // Needed to send/receive strings!

#include <boost/mpi.hpp>
#include <feel/feelconfig.h>

#if defined(FEELPP_HAS_MPI_H)
#include <mpi.h>
#endif /* FEELPP_HAS_MPI_H */



#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>

namespace Feel
{
namespace mpi = boost::mpi;

/**
 * \class Application
 *\ingroup Core
 *\brief MPI Application
 *
 * @author Christophe Prud'homme
 * @see
 */
class Application : public Application
{
    typedef Application super;
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
#if defined( FEELPP_HAS_MPI )
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
#if defined( FEELPP_HAS_MPI )
    Application( int argc, char** argv, AboutData const& ad, po::options_description const& od, MPI_Comm Comm = MPI_COMM_WORLD );
#else
    Application( int argc, char** argv, AboutData const& ad, po::options_description const& od );
#endif


    ~Application();

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return \p true if MPI is initialized, \p false otherwise
     */
    bool isMPIInitialized() const
    {
        return _S_is_mpi_initialized;
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


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{




#if defined( FEELPP_HAS_MPI )
    template<class T>
    static void Send( const T& obj, int proc, int tag, const MPI_Comm& comm = COMM_WORLD )
    {
        std::ostringstream oss;
        boost::archive::binary_oarchive oa( oss );
        oa << obj;

        std::string str( oss.str() );
        int size = str.size();
        MPI_Ssend( &size, 1, MPI_INT, proc, tag, comm );
        MPI_Ssend( &str[0], size, MPI_CHAR, proc, tag, comm );
    }
#else
    template<class T>
    static void Send( const T& obj, int proc, int tag )
    {
        // dummy function, don't have to do anything
    }
#endif // FEELPP_HAS_MPI


#if defined( FEELPP_HAS_MPI )
    template<class T>
    static void Broadcast( T& obj, int root = 0, const MPI_Comm& comm = COMM_WORLD )
    {
        std::ostringstream oss;
        boost::archive::binary_oarchive oa( oss );
        oa << obj;

        std::string str( oss.str() );
        int size = str.size();
        MPI_Bcast( &size, 1, MPI_INT, root, comm );

        if ( Application::processId() != root )
            str.resize( size );

        DVLOG(2) << "[Application::Broadcast] str.size = " << str.size() << "\n";
        DVLOG(2) << "[Application::Broadcast] before str = " << str << "\n";

        MPI_Bcast( &str[0], size, MPI_CHAR, root, comm );
        DVLOG(2) << "[Application::Broadcast] after str = " << str << "\n";

        // deserialize for processId() != root
        if ( Application::processId() != root )
        {
            std::istringstream iss( str );
            boost::archive::binary_iarchive ia( iss );
            ia >> obj;
        }
    }
#else
    template<class T>
    static void Broadcast( T& obj )
    {
        // dummy function, don't have to do anything
    }
#endif // FEELPP_HAS_MP




#if defined( FEELPP_HAS_MPI )
    template<class T>
    static void Recv( T& obj, int proc, int tag, const MPI_Comm& comm = COMM_WORLD )
    {

        int size;
        MPI_Status status;
        MPI_Recv( &size, 1, MPI_INT, proc, tag, comm, &status );
        std::string buf( size, ' ' );
        MPI_Recv( &buf[0], size, MPI_CHAR, proc, tag, comm, &status );

        std::istringstream iss( buf );
        boost::archive::binary_iarchive ia( iss );
        ia >> obj;
    }
#else
    template<class T>
    static void Recv( T& obj, int proc, int tag )
    {
        // dummy function, don't have to do anything
    }
#endif // FEELPP_HAS_MPI
    //@}

#if defined( FEELPP_HAS_MPI )
    static MPI_Comm COMM_WORLD;
#endif // FEELPP_HAS_MPI

    /**
     * @return the communicator
     */
    static mpi::communicator const& comm()
    {
        return S_world;
    }

    /**
     * @return the barrier
     */
    static void barrier()
    {
        S_world.barrier();
    }


protected:


private:

    Application( Application const & );

private:



    static bool _S_is_mpi_initialized;
    boost::shared_ptr<mpi::environment> M_env;
    static mpi::communicator S_world;
};


}
#endif /* __Application_H */
