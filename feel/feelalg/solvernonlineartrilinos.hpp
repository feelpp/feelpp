/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Florent Vielfaure <florent.vielfaure@gmail.com>
       Date: 2009-05-25

  Copyright (C) 2009-2012 Universite Joseph Fourier (Grenoble I)

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
   \file solvernonlineartrilinos.hpp
   \author Florent Vielfaure <florent.vielfaure@gmail.com>
   \date 2009-05-25
 */

#ifndef __SolverNonLinearTrilinos_H
#define __SolverNonLinearTrilinos_H 1

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/solvernonlinear.hpp>
#include <feel/feelalg/matrixepetra.hpp>
#include <feel/feelalg/operatortrilinos.hpp>
#include <feel/feelalg/vectorepetra.hpp>
// #include <feel/feelalg/operatortrilinos.hpp>

#if defined( FEELPP_HAS_TRILINOS )
#include "NOX.H"
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_LinearSystem_AztecOO.H"
#include "NOX_Epetra_Group.H"

namespace Feel
{
/**
 * \class SolverNonLinearTrilinos
 * \brief Trilinos non linear solvers interface
 *
 * This class provides an interface to Trilinos/NOX iterative solvers that is
 * compatible with the \p SolverNonLinear<> base class
 *
 * @author Florent Vielfaure
 */
template<typename T>
class SolverNonLinearTrilinos
    :
public SolverNonLinear<T>
{
    typedef SolverNonLinear<T> super;
public:


    /** @name Typedefs
     */
    //@{

    typedef SolverNonLinearTrilinos<T> self_type;

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
     *  Constructor. Initializes Trilinos data structures
     */
    SolverNonLinearTrilinos( std::string const& prefix = "", WorldComm const& worldComm=Environment::worldComm() );
    SolverNonLinearTrilinos( SolverNonLinearTrilinos const & );

    /**
     * Destructor.
     */
    ~SolverNonLinearTrilinos();

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
     * Call the Trilinos solver.  It calls the method below, using the
     * same matrix for the system and preconditioner matrices.
     */
    virtual std::pair<int, real_type> solve ( sparse_matrix_ptrtype&,    // System Jacobian Matrix
            vector_ptrtype&,          // Solution vector
            vector_ptrtype&,          // Residual vector
            const double,        // Stopping tolerance
            const unsigned int ); // N. Iterations

    virtual std::pair<unsigned int, real_type> solve ( dense_matrix_type&,    // System Jacobian Matrix
            dense_vector_type&,          // Solution vector
            dense_vector_type&,          // Residual vector
            const double,        // Stopping tolerance
            const unsigned int ); // N. Iterations



    //@}

    /**
     * Methods for NOX::Epetra::Interface::Required and NOX::Epetra::Interface::Jacobian
     */
    //bool computeF( const Epetra_Vector & x, Epetra_Vector & f, NOX::Epetra::Interface::Required::FillType F );
    //bool computeJacobian( const Epetra_Vector & x, Epetra_Operator & Jac );
    //bool computePrecMatrix( const Epetra_Vector & x, Epetra_RowMatrix & M );
    //bool computePreconditioner( const Epetra_Vector & x, Epetra_Operator & O );

private:

};

template <typename T>
inline
SolverNonLinearTrilinos<T>::SolverNonLinearTrilinos ( std::string const& prefix, WorldComm const& worldComm )
{}



template <typename T>
inline
SolverNonLinearTrilinos<T>::~SolverNonLinearTrilinos ()
{
    this->clear ();
}

} // Feel
#endif// FEELPP_HAS_TRILINOS_NOX
#endif /* __SolverNonLinearTrilinos_H */
