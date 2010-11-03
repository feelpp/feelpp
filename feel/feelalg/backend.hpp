/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-12-23

  Copyright (C) 2007,2008,2009,2010 Université Joseph Fourier (Grenoble I)

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
   \file backend.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-12-23
 */
#ifndef __Backend_H
#define __Backend_H 1

#include <boost/timer.hpp>
#include <boost/tuple/tuple.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/parameter.hpp>
#include <feel/feelalg/enums.hpp>
#include <feel/feelalg/vector.hpp>
#include <feel/feelalg/matrixsparse.hpp>
#include <feel/feelalg/datamap.hpp>

#include <feel/feelalg/solvernonlinear.hpp>

namespace Feel
{
///! \cond detail
namespace detail
{
template<typename T>
DataMap datamap( T const& t, mpl::true_ )
{
    return t->map();
}
template<typename T>
DataMap datamap( T const& t, mpl::false_ )
{
    return t.map();
}
template<typename T>
DataMap datamap( T const& t )
{
    return datamap( t, detail::is_shared_ptr<T>() );
}

template<typename T>
typename T::reference ref( T t, mpl::true_ )
{
    return *t;
}
template<typename T>
T& ref( T& t, mpl::false_ )
{
    return t;
}
template<typename T>
auto ref( T& t ) -> decltype( ref( t, detail::is_shared_ptr<T>() ) )
{
    return ref( t, detail::is_shared_ptr<T>() );
}

}
///! \endcond detail
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
     * default, otherwise use GMM which is distributed with feel
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
        return this->newMatrix( dm->map(), im->map() );
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

    /**
     * \return the relative tolerance
     */
    value_type rTolerance() const { return M_rtolerance;}

    /**
     * \return the divergence tolerance
     */
    value_type dTolerance() const { return M_dtolerance;}

    /**
     * \return the absolute tolerance
     */
    value_type aTolerance() const { return M_atolerance;}

    /**
     * \return the maximum number of iterations
     */
    size_type maxIterations() const { return M_maxit; }

    bool converged() const { return M_converged; }

    size_type nIterations() const { return M_iteration; }

    bool transpose() const { return M_transpose; }

    /**
     * \return the communicator
     */
    mpi::communicator const& comm() const { return M_comm; }

    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * set tolerances: relative tolerance \p rtol, divergence tolerance \p dtol
     * and absolute tolerance \p atol
     */
    BOOST_PARAMETER_MEMBER_FUNCTION((void),
                                    setTolerances,
                                    tag,
                                    (required
                                     (rtolerance, (double))
                                        )
                                    (optional
                                     (maxit,      (size_type), 1000 )
                                     (atolerance, (double),    1e-50)
                                     (dtolerance, (double),    1e5)
                                        ) )
        {
            M_rtolerance = rtolerance;
            M_dtolerance = dtolerance;
            M_atolerance = atolerance;
            M_maxit = maxit;
        }

    /**
     * set the type of preconditioner associated to the matrix
     */
    void setPrecMatrixStructure( MatrixStructure mstruct ) {  M_prec_matrix_structure = mstruct; }

    /**
     * \return the non linear solver
     */
    solvernonlinear_ptrtype nlSolver() { return M_nlsolver; }

    void setTranspose( bool transpose ) { M_transpose = transpose; }

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
     * solve for \f$P A x = P b\f$ where \f$P\f$ is an approximation
     * of the inverse of \f$A\f$. this interface uses the
     * boost.parameter library to ease the function usage
     *
     * \param A matrix to inverse
     * \param rhs right hand side vector
     * \param solution solution of the system
     * \param P preconditioner
     * \param maxit maximum number of iterations
     * \param tolerance tolerance on the residual
     * \param reuse_prec if true use adaptive preconditioning strategy
     * \param transpose if true solve the transpose problem
     *
     * \warning some parameter may not be meaningful for all backends
     */
    BOOST_PARAMETER_MEMBER_FUNCTION((solve_return_type),
                                    solve,
                                    tag,
                                    (required
                                     (matrix,(sparse_matrix_ptrtype))
                                     (in_out(solution),*(mpl::or_<boost::is_convertible<mpl::_,vector_type&>,
                                                         boost::is_convertible<mpl::_,vector_ptrtype> >))
                                     (rhs,(vector_ptrtype)))
                                    (optional
                                     (prec,(sparse_matrix_ptrtype), matrix )
                                     (maxit,(size_type), 1000 )
                                     (rtolerance,(double), 1e-13)
                                     (atolerance,(double), 1e-50)
                                     (dtolerance,(double), 1e5)
                                     (reuse_prec,(bool), false )
                                     (transpose,(bool), false )
                                     )
                                    )
    {
        this->setTolerances( _dtolerance=dtolerance,
                             _rtolerance=rtolerance,
                             _atolerance=atolerance,
                             _maxit=maxit );
        vector_ptrtype _sol( this->newVector( detail::datamap(solution) ) );
        // initialize
        *_sol = detail::ref(solution);
        this->setTranspose( transpose );
        solve_return_type ret;
        if ( reuse_prec == false )
            ret = solve( matrix, prec, _sol, rhs );
        else
            ret = solve( matrix, prec, _sol, rhs, reuse_prec );
        detail::ref(solution) = *_sol;
        return ret;
    }

    /**
     * solve for \f$P A x = P b\f$ where \f$P\f$ is an approximation of
     * the inverse of \f$A\f$.
     *
     * \param A matrix to inverse
     * \param rhs right hand side vector
     * \param solution solution of the system
     * \param P preconditioner
     * \param maxit maximum number of iterations
     * \param tolerance tolerance on the residual
     * \param transpose if true solve the transpose problem
     *
     * \warning some parameter may not be meaningful for all backends
     */
    virtual solve_return_type solve( sparse_matrix_ptrtype const& A,
                                     sparse_matrix_ptrtype const& P,
                                     vector_ptrtype& x,
                                     vector_ptrtype const& b
                                     ) = 0;

    /**
     * solve for \f$P A x = P b\f$ where \f$P\f$ is an approximation
     * of the inverse of \f$A\f$ using an adaptive preconditioning
     * strategy.
     *
     * \param A matrix to inverse
     * \param rhs right hand side vector
     * \param solution solution of the system
     * \param P preconditioner
     * \param maxit maximum number of iterations
     * \param tolerance tolerance on the residual
     * \param transpose if true solve the transpose problem
     *
     * \warning some parameter may not be meaningful for all backends
     */
    solve_return_type solve( sparse_matrix_ptrtype const& A,
                             sparse_matrix_ptrtype const& P,
                             vector_ptrtype& x,
                             vector_ptrtype const& b,
                             bool reuse_prec
                             );
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

    mpi::communicator M_comm;

    BackendType M_backend;

    solvernonlinear_ptrtype M_nlsolver;

    MatrixStructure M_prec_matrix_structure;

    double M_totalSolveIter;
    double M_lastSolveIter;
    double M_firstSolveTime;
    double M_residual;
    double M_rtolerance;
    double M_dtolerance;
    double M_atolerance;
    size_t M_nUsePC;
    bool   M_converged;
    bool   M_reusePC;
    bool   M_reusedPC;
    bool   M_reuseFailed;
    boost::timer M_timer;
    bool   M_transpose;
    size_type    M_maxit;
    size_type    M_iteration;

};

/**
 * command line options
 */
po::options_description backend_options();
}
#endif /* __Backend_H */
