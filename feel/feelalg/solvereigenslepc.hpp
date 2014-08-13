/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-07-04

  Copyright (C) 2007-2011 Universite Joseph Fourier (Grenoble I)

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
   \file solvereigenslepc.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-07-04
 */
#ifndef __SolverEigenSlepc_H
#define __SolverEigenSlepc_H 1

#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelalg/vectorpetsc.hpp>
#include <feel/feelalg/matrixpetsc.hpp>

/**
 * SLEPc include files. SLEPs can only be used
 * together with PETSc.
 */
#if defined(FEELPP_HAS_SLEPC) && defined(FEELPP_HAS_PETSC)
#ifndef USE_COMPLEX_NUMBERS
extern "C"
{
# include <slepceps.h>
#if (SLEPC_VERSION_MAJOR == 3) && (SLEPC_VERSION_MINOR >= 5)
# include <slepcbv.h>
#else
# include <slepcip.h>
#endif
}
#else
# include <slepceps.h>
#if (SLEPC_VERSION_MAJOR == 3) && (SLEPC_VERSION_MINOR >= 5)
# include <slepcbv.h>
#else
# include <slepcip.h>
#endif
#endif


namespace Feel
{
/**
 * \class SolverEigenSlepc
 * \brief Interface to Slepc
 *
 * This class provides an interface to the SLEPc eigenvalue solver
 * library \p www.grycap.upv.es/slepc/.
 *
 *  @author Christophe Prud'homme
 *  @see
 */
template<typename T>
class SolverEigenSlepc : public SolverEigen<T>
{
    typedef SolverEigen<T> super;
public:


    /** @name Typedefs
     */
    //@{

    typedef typename super::value_type value_type;
    typedef typename super::real_type real_type;

    typedef typename super::solvereigen_type solvereigen_type;
    typedef typename super::solvereigen_ptrtype solvereigen_ptrtype;

    typedef typename super::vector_ptrtype vector_ptrtype;
    typedef typename super::sparse_matrix_type sparse_matrix_type;

    typedef typename super::solve_return_type solve_return_type;
    typedef typename super::eigenpair_type eigenpair_type;

    typedef typename super::eigenmodes_type eigenmodes_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     *  Constructor. Initializes Petsc data structures
     */
    SolverEigenSlepc()
    {
        this->M_eigen_solver_type = KRYLOVSCHUR;
    }

    SolverEigenSlepc( po::variables_map const& vm, std::string const& prefix = "" )
        :
        super( vm, prefix )
    {
    }

    SolverEigenSlepc( SolverEigenSlepc const & );

    /**
     * Destructor.
     */
    ~SolverEigenSlepc()
    {
        this->clear ();
    }




    /**
     * Initialize data structures if not done so already.
     */
    void init();

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{


    /**
     * This function returns the real and imaginary part of the ith
     * eigenvalue and copies the respective eigenvector to the solution
     * vector. Note that also in case of purely real matrix entries the
     * eigenpair may be complex values.
     */
    eigenpair_type eigenPair ( unsigned int i );

    /**
     * Returns the eigen modes in a map
     */
    virtual eigenmodes_type eigenModes () ;

    /**
     * @computes and returns the relative error
     * ||A*x-lambda*x||/|lambda*x| of the ith eigenpair.
     */
    real_type relativeError ( unsigned int i );

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
    void clear();


    /**
     * This function calls the SLEPc solver to compute the eigenpairs
     * of matrix matrix_A. \p nev is the number of eigenpairs to be
     * computed and \p ncv is the number of basis vectors to be used
     * in the solution procedure. Return values are the number of
     * converged eigen values and the number of the iterations carried
     * out by the eigen solver.
     */
    solve_return_type  solve ( MatrixSparse<T> &matrix_A,
                               int nev,
                               int ncv,
                               const double tol,
                               const unsigned int m_its );

    /**
     * This function calls the SLEPc solver to compute the eigenpairs
     * of the generalized eigenvalue problem \f$ Ax = \lambda Bx\f$.
     * \p nev is the number of eigenpairs to be computed and \p
     * ncv is the number of basis vectors to be used in the solution
     * procedure. Return values are the number of converged eigen
     * values and the number of the iterations carried out by the
     * eigen solver.
     */
    solve_return_type  solve ( MatrixSparse<T> &matrix_A,
                               MatrixSparse<T> &matrix_B,
                               int nev,
                               int ncv,
                               const double tol,
                               const unsigned int m_its );



    //@}



private:

    /**
     * Tells Slepc to use the user-specified solver stored in
     * \p _eigen_solver_type
     */
    void setSlepcSolverType ();

    /**
     * Tells Slepc to deal with the type of problem stored in
     * \p _eigen_problem_type
     */
    void setSlepcProblemType ();

    /**
     * Tells Slepc to compute the spectrum at the position
     * stored in \p _position_of_spectrum
     */
    void setSlepcPositionOfSpectrum();

    /**
     * set the spectral transforms
     */
    void setSlepcSpectralTransform();

    /**
     * Eigenproblem solver context
     */
    EPS M_eps;

    /**
     * Eigenproblem inner products
     */
#if (SLEPC_VERSION_MAJOR == 3) && (SLEPC_VERSION_MINOR >= 5)
    BV M_ip;
#else
    IP M_ip;
#endif
    /**
     * eigenmode
     */
    Vec M_mode;
};
po::options_description solvereigenslepc_options( std::string const& prefix = "" );

} // Feel

#endif // SLEPC PETSC
#endif /* __SolverEigenSlepc_H */
