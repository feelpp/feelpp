/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-07-02

  Copyright (C) 2007, 2009 Université Joseph Fourier (Grenoble I)

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
   \file solvernonlinear.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-07-02
 */
#ifndef __SolverNonLinear_H
#define __SolverNonLinear_H 1

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <feel/feelalg/enums.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelcore/traits.hpp>

namespace Feel
{

template<typename T> class Vector;
template<typename T> class MatrixSparse;

/**
 * \class SolverNonLinear
 * \brief Non linear solver base interface
 *
 * This class provides a uniform interface for nonlinear solvers.
 * This base class is overloaded to provide nonlinear solvers from
 * different packages like PETSC.
 *
 *  @author Christophe Prud'homme
 */
template <typename T>
class SolverNonLinear
{
public:


    /** @name Typedefs
     */
    //@{

    typedef SolverNonLinear<T> self_type;
    typedef boost::shared_ptr<SolverNonLinear<T> > self_ptrtype;
    typedef self_type solvernonlinear_type;
    typedef boost::shared_ptr<self_type> solvernonlinear_ptrtype;

    typedef T value_type;
    typedef typename type_traits<T>::real_type real_type;


    typedef boost::shared_ptr<Vector<value_type> > vector_ptrtype;
    typedef boost::shared_ptr<MatrixSparse<value_type> > sparse_matrix_ptrtype;

    typedef ublas::matrix<value_type> dense_matrix_type;
    typedef ublas::vector<value_type> dense_vector_type;

    typedef boost::function<void (const vector_ptrtype& X,
                                  vector_ptrtype& R)> residual_function_type;
    typedef boost::function<void (const vector_ptrtype& X,
                                  sparse_matrix_ptrtype& J)> jacobian_function_type;
    typedef boost::function<void (const vector_ptrtype& X,
                                  vector_ptrtype& R,
                                  sparse_matrix_ptrtype& J)> matvec_function_type;

    typedef boost::function<void ( dense_vector_type const& X,
                                   dense_vector_type & R)> dense_residual_function_type;
    typedef boost::function<void ( dense_vector_type const& X,
                                   dense_matrix_type& J)> dense_jacobian_function_type;
    typedef boost::function<void ( dense_vector_type const& X,
                                   dense_vector_type& R,
                                   dense_matrix_type& J)> dense_matvec_function_type;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     *  Constructor. Initializes Solver data structures
     */
    SolverNonLinear();

    /**
     * copy constructor
     */
    SolverNonLinear( SolverNonLinear const & );

    /**
     * Destructor.
     */
    virtual ~SolverNonLinear();


    /**
     * Builds a \p NonlinearSolver using the nonlinear solver package specified by
     * the \p variables_map \p vm and  \p prefix
     */
    static solvernonlinear_ptrtype build( po::variables_map const& vm, std::string const& prefix = "" );

    /**
     * Builds a \p NonlinearSolver using the nonlinear solver package specified by
     * \p solver_package
     */
    static solvernonlinear_ptrtype build( SolverPackage solver_package );

    /**
     * Initialize data structures if not done so already.
     */
    virtual void init () = 0;


    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /**
     * @returns true if the data structures are
     * initialized, false otherwise.
     */
    bool initialized () const { return M_is_initialized; }

    /**
     * Release all memory and clear data structures.
     */
    virtual void clear () {}

    /**
     * \return the preconditioner matrix structure
     * it may not be relevant to all non linear solvers
     */
    virtual MatrixStructure precMatrixStructure() const { return M_prec_matrix_structure; }

    SolverNonLinearType getType() const { return M_snl_type; }

    double getAbsoluteResidualTol() const { return M_absoluteResidualTol; }
    double getRelativeResidualTol() const { return M_relativeResidualTol; }
    double getAbsoluteSolutionTol() const { return M_absoluteSolutionTol; }

    uint getNbItMax() const { return M_nbItMax; }
    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * \return the preconditioner matrix structure
     * it may not be relevant to all non linear solvers
     */
    virtual void setPrecMatrixStructure( MatrixStructure mstruct  ) { M_prec_matrix_structure = mstruct; }


    /**
     * Select type of non linear solver : LINEAR_SEARCH, TRUST_REGION, ...
     */
    void setType(SolverNonLinearType snl_type) { M_snl_type=snl_type; }
    /**
     * Returns the type of solver to use.
     */
    SolverNonLinearType nlSolverType () const { return M_snl_type; }

    /**
     * Sets the type of solver to use.
     */
    void setKspSolverType (const SolverType st) { M_kspSolver_type = st; }

    /**
     * Returns the type of solver to use.
     */
    SolverType kspSolverType () const { return M_kspSolver_type; }

    /**
     * Sets the type of preconditioner to use.
     */
    void setPreconditionerType (const PreconditionerType pct) { M_preconditioner_type = pct; }
    /**
     * Returns the type of preconditioner to use.
     */
    PreconditionerType preconditionerType () const { return M_preconditioner_type; }

    /**
     * Define values of tolerance for the non linear solver
     */
    void setRelativeResidualTol(double tol) { M_relativeResidualTol = tol; }
    void setAbsoluteResidualTol(double tol) { M_absoluteResidualTol = tol; }
    void setAbsoluteSolutionTol(double tol) { M_absoluteSolutionTol = tol; }

    /**
     * Define the number max of iteration
     */
    void setNbItMax(uint n) { M_nbItMax=n; }
    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Solves a sparse nonlinear system.
     */
    virtual std::pair<unsigned int, real_type> solve ( sparse_matrix_ptrtype&,  // System Jacobian Matrix
                                                       vector_ptrtype&, // Solution vector
                                                       vector_ptrtype&, // Residual vector
                                                       const double,      // Stopping tolerance
                                                       const unsigned int) = 0; // N. Iterations

    /**
     * Solves a sparse nonlinear system.
     */
    virtual std::pair<unsigned int, real_type> solve ( dense_matrix_type&,  // System Jacobian Matrix
                                                       dense_vector_type&, // Solution vector
                                                       dense_vector_type&, // Residual vector
                                                       const double,      // Stopping tolerance
                                                       const unsigned int) = 0; // N. Iterations

    /**
     * Function that computes the residual \p R(X) of the nonlinear system
     * at the input iterate \p X.
     */
    residual_function_type residual;

    /**
     * Function that computes the Jacobian \p J(X) of the nonlinear system
     * at the input iterate \p X.
     */
    jacobian_function_type jacobian;

    /**
     * Function that computes either the residual \f$ R(X) \f$ or the
     * Jacobian \f$ J(X) \f$ of the nonlinear system at the input
     * iterate \f$ X \f$.  Note that either \p R or \p J could be
     * \p XSNULL.
     */
    matvec_function_type matvec;

    /**
     * Function that computes the residual \p R(X) of the nonlinear system
     * at the input iterate \p X.
     */
    dense_residual_function_type dense_residual;

    /**
     * Function that computes the Jacobian \p J(X) of the nonlinear system
     * at the input iterate \p X.
     */
    dense_jacobian_function_type dense_jacobian;

    /**
     * Function that computes either the residual \f$ R(X) \f$ or the
     * Jacobian \f$ J(X) \f$ of the nonlinear system at the input
     * iterate \f$ X \f$.  Note that either \p R or \p J could be
     * \p XSNULL.
     */
    dense_matvec_function_type dense_matvec;

    //@}



protected:

    /**
     * Flag indicating if the data structures have been initialized.
     */
    bool M_is_initialized;

    MatrixStructure M_prec_matrix_structure;

    /**
     * Define the type of non linear solver
     */
    SolverNonLinearType M_snl_type;

    /**
     * Enum stating which type of iterative linear solver to use.
     */
    SolverType M_kspSolver_type;

    /**
     * Enum statitng with type of preconditioner to use.
     */
    PreconditionerType M_preconditioner_type;

    /**
     * Two differents tolerances on the residual for the resolution of non linear system
     */
    double M_relativeResidualTol;
    double M_absoluteResidualTol;

    /**
     * Absolute tolerances between successive iteration
     */
    double M_absoluteSolutionTol;

    /**
     * The maximum numbers of allowable nonlinear iterations
     */
    uint M_nbItMax;
};


/**
 * command line options
 */
po::options_description nlsolver_options();
} // Feel
#endif /* __SolverNonLinear_H */
