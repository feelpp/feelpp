/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-11-27

  Copyright (C) 2005,2006 EPFL

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
// $Id: petsc_linear_solver.h,v 1.2 2005/02/22 22:17:34 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
/**
   \file solverlinearpetsc.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-11-27
 */

#ifndef __petsc_linear_solver_h__
#define __petsc_linear_solver_h__


#include <feel/feelcore/feel.hpp>

#include <feel/feelcore/feelpetsc.hpp>
#include <feel/feelalg/solverlinear.hpp>

#include <feel/feelalg/matrixshell.hpp>
#include <feel/feelalg/matrixpetsc.hpp>
#include <feel/feelalg/vectorpetsc.hpp>


/**
 * Petsc include files.
 */
#ifdef FEELPP_HAS_PETSC_H
#include <feel/feelcore/feelpetsc.hpp>

#ifndef USE_COMPLEX_NUMBERS
extern "C" {
# include <petscversion.h>
# if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1)
#   include <petscsles.h>
# else
#   include <petscksp.h>
# endif
}
#else
# include <petscversion.h>
# if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1)
#   include <petscsles.h>
# else
#   include <petscksp.h>
# endif
#endif

//--------------------------------------------------------------------
// Functions with C linkage to pass to PETSc.  PETSc will call these
// methods as needed for preconditioning
//
// Since they must have C linkage they have no knowledge of a namespace.
// Give them an obscure name to avoid namespace pollution.
extern "C"
{
#if PETSC_VERSION_LESS_THAN(3,0,1) && PETSC_VERSION_RELEASE
    /**
     * This function is called by PETSc to initialize the preconditioner.
     * ctx will hold the Preconditioner.
     */
    PetscErrorCode __feel_petsc_preconditioner_setup ( void * ctx );

    /**
     * This function is called by PETSc to acctually apply the preconditioner.
     * ctx will hold the Preconditioner.
     */
    PetscErrorCode __feel_petsc_preconditioner_apply( void *ctx, Vec x, Vec y );

    /**
     * This function is called by PETSc to view the preconditioner.
     * ctx will hold the Preconditioner.
     */
    PetscErrorCode __feel_petsc_preconditioner_view( void *ctx, PetscViewer viewer)
#else
    PetscErrorCode __feel_petsc_preconditioner_setup ( PC );
    PetscErrorCode __feel_petsc_preconditioner_apply( PC, Vec x, Vec y );
    PetscErrorCode __feel_petsc_preconditioner_view( PC, PetscViewer viewer);

#endif
} // end extern "C"



namespace Feel
{


/**
 * This class provides an interface to PETSc
 * iterative solvers that is compatible with the \p libMesh
 * \p LinearSolver<>
 *
 * @author Benjamin Kirk, 2002-2005
 * @author Christophe Prud'homme, 2005
 */

template <typename T>
class SolverLinearPetsc : public SolverLinear<T>
{
    typedef SolverLinear<T> super;

public:

    typedef typename super::value_type value_type;
    typedef typename super::real_type real_type;
    typedef typename super::solve_return_type solve_return_type;

    /**
     *  Constructor. Initializes Petsc data structures
     */
    SolverLinearPetsc ( WorldComm const& worldComm=Environment::worldComm() );

    /**
     *  Constructor. Initializes Petsc data structures
     */
    SolverLinearPetsc ( po::variables_map const& vm, WorldComm const& worldComm=Environment::worldComm() );

    /**
     * Destructor.
     */
    ~SolverLinearPetsc ();

    /**
     * Release all memory and clear data structures.
     */
    void clear ();

    /**
     * Initialize data structures if not done so already.
     */
    void init ();

    /**
     * if \p cns is true, set the null space to be the constant values
     */
    void setConstantNullSpace( bool cns )
    {
        M_constant_null_space = cns;
    }

    /**
     * Call the Petsc solver.  It calls the method below, using the
     * same matrix for the system and preconditioner matrices.
     *
     * \param mat System Matrix
     * \param prec Preconditioning Matrix
     * \param x Solution vector
     * \param b RHS vector
     * \param tolerance Stopping tolerance
     * \param maxit maximum Number of Iterations
     * \param transpose true to solve the transpose system, false otherwise
     */
    solve_return_type
    solve ( MatrixSparse<T>  const &mat,
            Vector<T> & x,
            Vector<T> const& b,
            const double tolerance,
            const unsigned int maxit,
            bool transpose )
    {
        return this->solve( mat, mat, x, b, tolerance, maxit, transpose );
    }

    solve_return_type
    solve ( MatrixShell<T>  const &mat,
            Vector<T> & x,
            Vector<T> const& b,
            const double tolerance,
            const unsigned int maxit,
            bool transpose );

    /**
     * This method allows you to call a linear solver while specifying
     * the matrix to use as the (left) preconditioning matrix.  Note
     * that the linear solver will not compute a preconditioner in
     * this case, and will instead premultiply by the matrix you
     * provide.
     *
     * In PETSc, this is accomplished by calling
     *
     * PCSetType(_pc, PCMAT);
     *
     * before invoking KSPSolve().  Note: this functionality is not
     * implemented in the SolverLinear class since there is not a
     * built-in analog to this method for LasPack -- You could
     * probably implement it by hand if you wanted.
     *
     *
     * \param mat System Matrix
     * \param prec Preconditioning Matrix
     * \param x Solution vector
     * \param b RHS vector
     * \param tolerance Stopping tolerance
     * \param maxit maximum Number of Iterations
     * \param transpose true to solve the transpose system, false otherwise
     */
    solve_return_type
    solve ( MatrixSparse<T>  const& mat,
            MatrixSparse<T>  const& prec,
            Vector<T> & x,
            Vector<T> const& b,
            const double tolerance,
            const unsigned int maxit,
            bool transpose );

    /**
     * @retun the Krylov SubsPace  data structure
     */
    KSP ksp()
        {
            this->init();
            return M_ksp;
        }
    /**
     * Returns the raw PETSc preconditioner context pointer.  This allows
     * you to specify the PCShellSetApply() and PCShellSetSetUp() functions
     * if you desire.  Just don't do anything crazy like calling PCDestroy()!
     */
    PC pc()
    {
        this->init();
        return M_pc;
    }

    /**
     * Fills the input vector with the sequence of residual norms
     * from the latest iterative solve.
     */
    void getResidualHistory( std::vector<double>& hist );

    /**
     * Returns just the initial residual for the solve just
     * completed with this interface.  Use this method instead
     * of the one above if you just want the starting residual
     * and not the entire history.
     */
    real_type getInitialResidual();

protected:
    void check( int err ) { CHKERRABORT( this->worldComm().globalComm(), err ); }
private:

    /**
     * Tells PETSC to use the user-specified solver stored in
     * \p _solver_type
     */
    void setPetscSolverType ();

    /**
     * Tells PETSC to use the user-specified preconditioner stored in
     * \p _preconditioner_type
     */
    void setPetscPreconditionerType ();

    void setPetscConstantNullSpace ();
    void updateNullSpace( Mat A, Vec rhs );
    void updateNearNullSpace( Mat A );

    // SLES removed from >= PETSc 2.2.0
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1)

    /**
     * Linear solver context
     */
    SLES M_sles;

#endif

    /**
     * Preconditioner context
     */
    PC M_pc;

    /**
     * Krylov subspace context
     */
    KSP M_ksp;

    bool M_constant_null_space;
};


/*----------------------- functions ----------------------------------*/
template <typename T>
inline
SolverLinearPetsc<T>::SolverLinearPetsc ( WorldComm const& worldComm )
    :
    super( worldComm ),
    M_constant_null_space( false )
{
    if ( this->worldComm().globalSize() == 1 )
        this->setPreconditionerType(  LU_PRECOND );

    else
        this->setPreconditionerType( BLOCK_JACOBI_PRECOND );
}

template <typename T>
inline
SolverLinearPetsc<T>::SolverLinearPetsc ( po::variables_map const& vm, WorldComm const& worldComm )
    :
    super( vm, worldComm ),
    M_constant_null_space( false )
{
}



template <typename T>
inline
SolverLinearPetsc<T>::~SolverLinearPetsc ()
{
    this->clear ();
}
} // Feel


#endif // #ifdef FEELPP_HAS_PETSC_H
#endif // #ifdef __petsc_linear_solver_h__
