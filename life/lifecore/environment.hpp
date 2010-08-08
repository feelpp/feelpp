/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2010-04-14

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2010-04-14
 */
#ifndef __Environment_H
#define __Environment_H 1

#include <cstdlib>

#include <boost/format.hpp>

#include <life/lifecore/life.hpp>

namespace Life
{
/**
 * \class Environment
 * \brief Runtime Environment
 *
 * @author Christophe Prud'homme
 * @see Application
 */
class Environment
{
public:


    /** @name Constants
     */
    //@{


    //@}

    /** @name Typedefs
     */
    //@{


    //@}

    /** @name Constructors, destructor
     */
    //@{

    /** Initialize the Life environment.
     *
     *  If the Life environment has not already been initialized,
     *  initializes Life
     */
    Environment();

    /** Initialize the Life environment.
     *
     *  If the Life environment has not already been initialized,
     *  initializes Life
     *
     *  @param argc The number of arguments provided in @p argv, as
     *  passed into the program's @c main function.
     *
     *  @param argv The array of argument strings passed to the program
     *  via @c main.
     *
     */
    Environment(int& argc, char** &argv);

    /** Shuts down the Life environment.
     *
     *  If this @c Environment object was used to initialize the Life
     *  environment, and the Life environment has not already been shut
     *  down (finalized), this destructor will shut down the Life
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

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    //! \return the root repository (default: \c $HOME/life)
    static std::string rootRepository();

    //! \return the local geo files repository (default: \c $HOME/life/geo)
    static std::string localGeoRepository();

    /**
     * \return a tuple : the system geo files repository (default: \c /usr/share/life/geo or /usr/local/share/life/geo) and true or false whether the directory exists or not
     */
    static boost::tuple<std::string,bool> systemGeoRepository();

    //! change the directory where the results are stored
    static void changeRepository( boost::format fmt );

    //@}



protected:

private:
    /// Whether this environment object called MPI_Init
    bool i_initialized;

    mpi::environment M_env;
    mpi::communicator M_comm;
};
}
#endif /* __Environment_H */
