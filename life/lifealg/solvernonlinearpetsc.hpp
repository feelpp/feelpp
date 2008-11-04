/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-07-02

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
   \file solvernonlinearpetsc.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-07-02
 */
#ifndef __SolverNonLinearPetsc_H
#define __SolverNonLinearPetsc_H 1

#include <life/lifecore/life.hpp>
#include <life/lifealg/solvernonlinear.hpp>
#include <life/lifealg/matrixpetsc.hpp>
#include <life/lifealg/vectorpetsc.hpp>

// Petsc include files.
#if defined( HAVE_PETSC_H )

#ifndef USE_COMPLEX_NUMBERS
extern "C" {
# include <petscversion.h>
# include <petsc.h>
# include <petscsnes.h>
}
#else
# include <petscversion.h>
# include <petsc.h>
# include <petscsnes.h>
#endif



/**

 *
 * @author Benjamin Kirk, 2002-2005
 */

namespace Life
{
/**
 * \class SolverNonLinearPetsc
 * \brief Petsc non linear solvers interface
 *
 * This class provides an interface to PETSc iterative solvers that is
 * compatible with the \p SolverNonLinear<> base class
 *
 * @author Christophe Prud'homme
 */
template<typename T>
class SolverNonLinearPetsc
    :
        public SolverNonLinear<T>
{
    typedef SolverNonLinear<T> super;
public:


    /** @name Typedefs
     */
    //@{

    typedef SolverNonLinearPetsc<T> self_type;

    typedef typename super::value_type value_type;
    typedef typename super::real_type real_type;
    typedef typename super::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename super::vector_ptrtype vector_ptrtype;

    typedef typename super::dense_matrix_type dense_matrix_type;
    typedef typename super::dense_vector_type dense_vector_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

      /**
       *  Constructor. Initializes Petsc data structures
       */
    SolverNonLinearPetsc();
    SolverNonLinearPetsc( SolverNonLinearPetsc const & );

    /**
     * Destructor.
     */
    ~SolverNonLinearPetsc();

    /**
     * Initialize data structures if not done so already.
     */
    virtual void init ();

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

    /**
     * Release all memory and clear data structures.
     */
    virtual void clear ();

    /**
     * Call the Petsc solver.  It calls the method below, using the
     * same matrix for the system and preconditioner matrices.
     */
    virtual std::pair<unsigned int, real_type> solve ( sparse_matrix_ptrtype&,    // System Jacobian Matrix
                                                       vector_ptrtype&,          // Solution vector
                                                       vector_ptrtype&,          // Residual vector
                                                       const double,        // Stopping tolerance
                                                       const unsigned int); // N. Iterations

    virtual std::pair<unsigned int, real_type> solve ( dense_matrix_type&,    // System Jacobian Matrix
                                                       dense_vector_type&,          // Solution vector
                                                       dense_vector_type&,          // Residual vector
                                                       const double,        // Stopping tolerance
                                                       const unsigned int); // N. Iterations



    //@}

private:
    /**
     * Nonlinear solver context
     */
    SNES M_snes;

    uint16_type M_prec_mat_structure;
};

template <typename T>
inline
SolverNonLinearPetsc<T>::SolverNonLinearPetsc ()
{
}



template <typename T>
inline
SolverNonLinearPetsc<T>::~SolverNonLinearPetsc ()
{
  this->clear ();
}

} // Life

#endif /* HAVE_PETSC */
#endif /* __SolverNonLinearPetsc_H */
