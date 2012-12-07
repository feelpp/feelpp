/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-04-25

  Copyright (C) 2006, 2009 EPFL

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



#if defined( FEELPP_HAS_TRILINOS_EPETRA )
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
//#include <Epetra_Vector.h>
#if defined(FEELPP_HAS_MPI)
#include <feel/feelcore/application.hpp>
#include <Epetra_MpiComm.h>
#else
#include <feel/feelcore/application.hpp>
#include <Epetra_SerialComm.h>
#endif /* FEELPP_HAS_MPI */
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
namespace Feel
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
#if defined(FEELPP_HAS_MPI)
    : public Application
#else
    : public Application
#endif
{
#if defined(FEELPP_HAS_MPI)
    typedef Application super;
#else
    typedef Application super;
#endif

public:


    /** @name Typedefs
     */
    //@{

#if defined(FEELPP_HAS_MPI)
    typedef Epetra_MpiComm comm_type;
#else
    typedef Epetra_SerialComm comm_type;
#endif /* FEELPP_HAS_MPI */

    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * Initialize the epetra application
     */
#if defined( FEELPP_HAS_MPI )
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
#if defined( FEELPP_HAS_MPI )
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
    static comm_type const& comm()
    {
        return *_S_comm;
    }

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
} // Feel
#endif /* FEELPP_HAS_TRILINOS_EPETRA */

#endif /* __Application_H */
