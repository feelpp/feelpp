// $Id: linear_solver.h,v 1.2 2005/02/22 22:17:34 jwpeterson Exp $

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

#ifndef __linear_solver_h__
#define __linear_solver_h__

#include <life/lifealg/enums.hpp>
#include <life/lifecore/traits.hpp>

namespace Life
{

template<typename T> class Vector;
template<typename T> class MatrixSparse;

/**
 * This class provides a uniform interface for linear solvers.  This base
 * class is overloaded to provide linear solvers from different packages
 * like LIFE, GMM or PETSC
 *
 * @author Benjamin Kirk, 2003
 * @author Christophe Prud'homme, 2005
 */
template <typename T>
class SolverLinear
{

public:

    typedef SolverLinear<T> self_type;
    typedef boost::shared_ptr<SolverLinear<T> >  self_ptrtype;

    typedef T value_type;
    typedef typename type_traits<T>::real_type real_type;

    /**
     *  Constructor. Initializes Solver data structures
     */
    SolverLinear ();

    /**
     * Destructor.
     */
    virtual ~SolverLinear ();

    /**
     * @returns true if the data structures are
     * initialized, false otherwise.
     */
    bool initialized () const { return _M_is_initialized; }


    /**
     * Release all memory and clear data structures.
     */
    virtual void clear () {}

    /**
     * Initialize data structures if not done so already.
     */
    virtual void init () = 0;

    /**
     * Returns the type of solver to use.
     */
    SolverType solverType () const { return _M_solver_type; }

    /**
     * Sets the type of solver to use.
     */
    void setSolverType (const SolverType st)
    { _M_solver_type = st; }

    /**
     * Returns the type of preconditioner to use.
     */
    PreconditionerType preconditionerType () const
    { return _M_preconditioner_type; }

    /**
     * Sets the type of preconditioner to use.
     */
    void setPreconditionerType (const PreconditionerType pct)
    { _M_preconditioner_type = pct; }

    /**
     * \return the preconditioner matrix structure
     * it may not be relevant to all non linear solvers
     */
    virtual MatrixStructure precMatrixStructure() const { return M_prec_matrix_structure; }

    /**
     * \return the preconditioner matrix structure
     * it may not be relevant to all non linear solvers
     */
    virtual void setPrecMatrixStructure( MatrixStructure mstruct  ) { M_prec_matrix_structure = mstruct; }

    /**
     * This function calls the solver "_M_solver_type" preconditioned
     * with the "_M_preconditioner_type" preconditioner.  Note that
     * this method will compute the preconditioner from the system
     * matrix.
     *
     * \param mat System Matrix
     * \param prec Preconditioning Matrix
     * \param x Solution vector
     * \param b RHS vector
     * \param tolerance Stopping tolerance
     * \param maxit maximum Number of Iterations
     * \param transpose true to solve the transpose system, false otherwise
     */
    virtual std::pair<unsigned int, real_type> solve (MatrixSparse<T> const& mat,
                                                      Vector<T>& x,
                                                      Vector<T> const& b,
                                                      const double tolerance,
                                                      const unsigned int maxit,
                                                      bool transpose
                                                      ) = 0;



    /**
     * This function calls the solver
     * "_M_solver_type" preconditioned with the
     * "_M_preconditioner_type" preconditioner.  Note that this method
     * will compute the preconditioner from the system matrix.
     *
     * \param mat System Matrix
     * \param prec Preconditioning Matrix
     * \param x Solution vector
     * \param b RHS vector
     * \param tolerance Stopping tolerance
     * \param maxit maximum Number of Iterations
     * \param transpose true to solve the transpose system, false otherwise
     */
    virtual std::pair<unsigned int, real_type> solve (MatrixSparse<T> const& mat,
                                                      MatrixSparse<T> const& prec,
                                                      Vector<T>& x,
                                                      Vector<T> const& b,
                                                      const double tolerance,
                                                      const unsigned int maxit,
                                                      bool transpose
                                                      ) = 0;


protected:

    /**
     * set initialized only for subclasses
     */
    void setInitialized( bool init )
    {
        _M_is_initialized = init;
    }

protected:


    /**
     * Enum stating which type of iterative solver to use.
     */
    SolverType _M_solver_type;

    /**
     * Enum statitng with type of preconditioner to use.
     */
    PreconditionerType _M_preconditioner_type;

    /**
     * Flag indicating if the data structures have been initialized.
     */
    bool _M_is_initialized;

    MatrixStructure M_prec_matrix_structure;
};




/*----------------------- inline functions ----------------------------------*/
template <typename T>
inline
SolverLinear<T>::SolverLinear () :

    _M_solver_type         (GMRES),
    _M_preconditioner_type (ILU_PRECOND),
    _M_is_initialized      (false),
    M_prec_matrix_structure( SAME_NONZERO_PATTERN )
{
}



template <typename T>
inline
SolverLinear<T>::~SolverLinear ()
{
    this->clear ();
}


} // Life

#endif // #ifdef __solver_h__
