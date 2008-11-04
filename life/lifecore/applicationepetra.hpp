/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2006-04-25

  Copyright (C) 2006 EPFL

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
   \file application.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-10-18
 */
#ifndef __Application_H
#define __Application_H 1



#if defined( HAVE_TRILINOS_EPETRA )
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
//#include <Epetra_Vector.h>
#if defined(HAVE_MPI)
#include <life/lifecore/application.hpp>
#include <Epetra_MpiComm.h>
#else
#include <life/lifecore/application.hpp>
#include <Epetra_SerialComm.h>
#endif /* HAVE_MPI */
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
namespace Life
{
/**
 * \class Application
 *\ingroup Core
 *\brief Epetra application
 *
 * @author Christophe Prud'homme
 * @see Application
 */
class Application
#if defined(HAVE_MPI)
    : public Application
#else
    : public Application
#endif
{
#if defined(HAVE_MPI)
    typedef Application super;
#else
    typedef Application super;
#endif

public:


    /** @name Typedefs
     */
    //@{

#if defined(HAVE_MPI)
    typedef Epetra_MpiComm comm_type;
#else
    typedef Epetra_SerialComm comm_type;
#endif /* HAVE_MPI */

    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * Initialize the epetra application
     */
#if defined( HAVE_MPI )
    Application( int argc,
                       char** argv,
                       AboutData const& ad,
                       MPI_Comm Comm = MPI_COMM_WORLD );
#else
    Application( int argc,
                       char** argv,
                       AboutData const& ad );
#endif
    /**
     * Initialize the epetra application and pass options to super classes
     */
#if defined( HAVE_MPI )
    Application( int argc,
                       char** argv,
                       AboutData const& ad,
                       po::options_description const& od,
                       MPI_Comm Comm = MPI_COMM_WORLD );
#else
    Application( int argc,
                       char** argv,
                       AboutData const& ad,
                       po::options_description const& od );
#endif
    /**
     * Finalize the epetra application
     */
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
     * \return the Epetra comm type
     */
    static comm_type const& comm() { return *_S_comm; }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{


    //@}



protected:

private:
    Application( Application const & );

private:
    static void init( MPI_Comm& comm );

    static bool _S_is_Initialized;

    static boost::shared_ptr<comm_type> _S_comm;
};
} // Life
#endif /* HAVE_TRILINOS_EPETRA */

#endif /* __Application_H */
