/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2010-04-14
 */
#ifndef __Environment_H
#define __Environment_H 1

#include <cstdlib>

#include <boost/noncopyable.hpp>
#include <boost/format.hpp>
#include <boost/signals2.hpp>

#include <feel/feelcore/feel.hpp>

namespace Feel
{
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

    //@}

    /** @name  Mutators
     */
    //@{


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

    //! change the directory where the results are stored
    static void changeRepository( boost::format fmt, std::string const& = "logfile" );

    //! get  \c variables_map from \c options_description \p desc
    static po::variables_map vm( po::options_description const& desc );

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
    //@}



private:
    /// Whether this environment object called MPI_Init
    bool i_initialized;
    mpi::environment M_env;

    static boost::signals2::signal<void ()> S_deleteObservers;
};
}
#endif /* __Environment_H */
