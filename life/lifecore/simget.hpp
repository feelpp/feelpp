/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2010-07-15

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
   \file simget.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2010-07-15
 */
#ifndef __Simget_H
#define __Simget_H 1

#include <life/lifecore/life.hpp>

namespace Life
{
/**
 * \class Simget
 * \brief Simulation Object
 *
 * A Simget is an object that provides two flavors of \c run() member function
 * - \c run() without any argument which simulates a blackbox \f$ F \f$ whitout
 *   any outputs or inputs
 * - <tt> run( double* X, int P, double* Y, int N ) </tt> which simulates a
 *    blackbox with input/output relationship \f$ Y = F(X) \f$ with \f$ Y \in
 *    \mathbb{R}^N\f$ and \f$ X \in \mathbb{R}^P\f$.
 *
 * @author Christophe Prud'homme
 * @see
 */
class Simget
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

    /**
     * constructor with a \c variables_map
     */
    Simget( po::variables_map const& vm ) : M_vm( vm ), M_about( "", "", "" ) {}

    /**
     * constructor with an \c AboutData that describes the simget
     */
    Simget( AboutData const& about ) : M_vm(), M_about( about ) {}

    /**
     * constructor with a \c variables_map and an \c AboutData that describes
     * the top application
     */
    Simget( po::variables_map const& vm, AboutData const& about ) : M_vm( vm ), M_about( about ) {}

    //! destructor
    virtual ~Simget() {}

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    Simget& operator=( Simget const & o)
        {
            if (this != &o )
            {
            }
            return *this;
        }
    //@}

    /** @name Accessors
     */
    //@{

    //! \return the mpi communicator
    mpi::communicator comm() const { return M_comm; }

    //! \return the \c variables_map
    po::variables_map const& vm() const { return M_vm; }

    //! \return the \c AboutData object
    AboutData const& about() const { return M_about; }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * simply execute the simget
     */
    virtual void run() = 0;

    /**
     * models the input/output relation \f$ Y=F(X) \f$
     */
    virtual void run( const double* X, unsigned long P, double* Y, unsigned long N ) = 0;

    //@}

protected:

private:
    mpi::communicator M_comm;
    po::variables_map M_vm;
    AboutData M_about;
};
}
#endif /* __Simget_H */
