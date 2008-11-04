/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-10-18

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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

#include <life/lifecore/application.hpp>

#if defined( HAVE_PETSC_H )
extern "C"
{
#include <petsc.h>
#include <petscerror.h>
}
#endif /* HAVE_PETSC_H */

namespace Life
{
#if defined( HAVE_PETSC_H )
/**
 * \class Application
 *\ingroup Core
 *\brief petsc application
 *
 * @author Christophe Prud'homme
 * @see Application
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
     * Initialize the petsc application
     */
    Application( int argc, char** argv, AboutData const& ad, MPI_Comm Comm = MPI_COMM_WORLD );

    /**
     * Initialize the petsc application and pass options to super classes
     */
    Application( int argc, char** argv, AboutData const& ad, po::options_description const& od, MPI_Comm Comm = MPI_COMM_WORLD );

    /**
     * Finalize the petsc application
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

};
#endif /* HAVE_PETSC_H */
} // Life
#endif /* __Application_H */
