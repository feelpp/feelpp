/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-12-23

  Copyright (C) 2007-2008 Université Joseph Fourier (Grenoble I)

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
   \file backend.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-12-23
 */
#ifndef __Backend_H
#define __Backend_H 1

#include <boost/shared_ptr.hpp>
#include <boost/timer.hpp>
#include <boost/tuple/tuple.hpp>

#include <life/lifealg/enums.hpp>
#include <life/lifealg/vector.hpp>
#include <life/lifealg/matrixsparse.hpp>
#include <life/lifealg/datamap.hpp>

#include <life/lifealg/solvernonlinear.hpp>

namespace Life
{
/**
 * \class Backend
 * \brief base class for all linear algebra backends
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename T>
class Backend
{
public:


    /** @name Typedefs
     */
    //@{
    typedef T value_type;
    typedef typename type_traits<T>::real_type real_type;

    typedef Vector<value_type> vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;
    typedef MatrixSparse<value_type> sparse_matrix_type;
    typedef boost::shared_ptr<sparse_matrix_type> sparse_matrix_ptrtype;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    typedef SolverNonLinear<value_type> solvernonlinear_type;
    typedef boost::shared_ptr<solvernonlinear_type> solvernonlinear_ptrtype;

    typedef boost::tuple<bool, size_type, value_type> solve_return_type;
    typedef boost::tuple<bool, size_type, value_type> nl_solve_return_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    Backend();
    Backend( po::variables_map const& vm, std::string const& prefix = "" );
    Backend( Backend const & );
    virtual ~Backend();

    /**
     * Builds a \p Backend, if Petsc is available, use Petsc by
     * default, otherwise use GMM which is distributed with life
     */
    static backend_ptrtype build(
#if defined( HAVE_PETSC_H )
                                 BackendType = BACKEND_PETSC
#else
                                 BackendType = BACKEND_GMM
#endif
                                 );

    /**
     * Builds a \p Backend
     */
    static backend_ptrtype build( po::variables_map const& vm, std::string const& prefix = "" );

    /**
     * instantiate a new sparse vector
     */
    virtual sparse_matrix_ptrtype newMatrix( DataMap const& dm1, DataMap const& dm2 ) = 0;

    /**
     * helper function
     */
    template<typename DomainSpace, typename ImageSpace>
    sparse_matrix_ptrtype newMatrix( DomainSpace const& dm, ImageSpace const& im  )
    {
        return this->newMatrix( im->map(), dm->map() );
    }

    /**
     * instantiate a new vector
     */
    virtual vector_ptrtype newVector( DataMap const& dm ) = 0;

    /**
     * helper function
     */
    template<typename DomainSpace>
    vector_ptrtype newVector( DomainSpace const& dm  )
    {
        return this->newVector( dm->map() );
    }

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return the type of preconditioner associated to the matrix
     */
    MatrixStructure precMatrixStructure() const { return M_prec_matrix_structure; }

    value_type tolerance() const { return M_tolerance;}

    int maxIterations() const { return M_maxiter; }

    bool converged() const { return M_converged; }

    int nIterations() const { return M_iteration; }

    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * set the type of preconditioner associated to the matrix
     */
    void setPrecMatrixStructure( MatrixStructure mstruct ) {  M_prec_matrix_structure = mstruct; }

    /**
     * \return the non linear solver
     */
    solvernonlinear_ptrtype nlSolver() { return M_nlsolver; }

    void setTolerance( value_type tol ) { M_tolerance = tol; }

    void setMaxIterations( int maxiter ) { M_maxiter = maxiter; }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * clean up
     */
    //virtual void clear() = 0;

    /**
     * \return \f$ r = x^T * y \f$
     */
    virtual real_type dot( vector_type const& x, vector_type const& y ) const;


    /**
     * \return \f$ r = x^T * y \f$
     */
    real_type dot( vector_ptrtype const& x, vector_ptrtype const& y ) const
    {
        return this->dot( *x, *y );
    }
    /**
     * \return \f$ y = A * x \f$
     */
    virtual void prod( sparse_matrix_type const& A, vector_type const& x, vector_type& y ) const = 0;

    /**
     * \return \f$ y = A * x \f$
     */
    void prod( sparse_matrix_ptrtype const& A, vector_ptrtype const& x, vector_ptrtype& y ) const
    {
        this->prod( *A, *x, *y );
    }



    /**
     * solve for \f$P A x = P b\f$ where \f$P\f$ is an approximation of
     * the inverse of \f$A\f$.
     */
    virtual solve_return_type solve( sparse_matrix_ptrtype const& A,
                                     sparse_matrix_ptrtype const& P,
                                     vector_ptrtype& x,
                                     vector_ptrtype const& b ) = 0;

    /**
     * solve for \f$P A x = P b\f$ where \f$P\f$ is an approximation of
     * the inverse of \f$A\f$ with an adaptive strategy to reuse the
     * preconditioner
     */
    virtual solve_return_type solve( sparse_matrix_ptrtype const& A,
                                     sparse_matrix_ptrtype const& P,
                                     vector_ptrtype& x,
                                     vector_ptrtype const& b,
                                     bool reusePC );

    /**
     * solve for the nonlinear problem \f$F( u ) = 0\f$
     */
    virtual nl_solve_return_type nlSolve( sparse_matrix_ptrtype& A,
                                          vector_ptrtype& x,
                                          vector_ptrtype& b,
                                          const double, const int );

    /**
     * solve for the nonlinear problem \f$F( u ) = 0\f$ with an
     * adaptive strategy to reuse the preconditioner
     */
    virtual nl_solve_return_type nlSolve( sparse_matrix_ptrtype& A,
                                          vector_ptrtype& x,
                                          vector_ptrtype& b,
                                          const double, const int,
                                          bool reusePC );
    //@}



protected:

private:

    void start();

    void stop();

    void reset();

private:
    BackendType M_backend;

    solvernonlinear_ptrtype M_nlsolver;

    MatrixStructure M_prec_matrix_structure;

    double M_totalSolveIter;
    double M_lastSolveIter;
    double M_firstSolveTime;
    double M_residual;
    double M_tolerance;
    size_t M_nUsePC;
    bool   M_converged;
    bool   M_reusePC;
    bool   M_reusedPC;
    bool   M_reuseFailed;
    boost::timer M_timer;
    int    M_maxiter;
    int    M_iteration;

};

/**
 * command line options
 */
po::options_description backend_options();
}
#endif /* __Backend_H */
