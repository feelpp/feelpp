/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
             Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2006-11-16

  Copyright (C) 2006 EPFL
  Copyright (C) 2007-2011 Université Joseph Fourier (Grenoble I)

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
   \file backendgmm.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2006-11-16
 */

#ifndef _BACKENDGMM_HPP_
#define _BACKENDGMM_HPP_

#include <boost/program_options/variables_map.hpp>

#include <feel/feelcore/application.hpp>
#include <feel/feelalg/vectorublas.hpp>
#include <feel/feelalg/matrixgmm.hpp>
#include <feel/feelalg/matrixtriplet.hpp>

#include <gmm_dense_Householder.h>
#include <gmm_iter_solvers.h>
#include <gmm_superlu_interface.h>
#include <feel/feelalg/solverumfpack.hpp>
#include <feel/feelalg/backend.hpp>

namespace Feel
{
namespace po = boost::program_options;

// struct to hold default values for Backendgmm
struct BackendGmmDefaults
{
    // Default constructor with default default values
    BackendGmmDefaults()
        :
        solver_type( "umfpack" ),
        pc_type( "ilut" ),
        threshold( 1e-3 ),
        fillin( 2 ),
        restart( 20 ),
        verbose( 0 ),
        maxiter( 1000 ),
        tolerance( 2e-10 )
    {}

    std::string solver_type;
    std::string pc_type;
    double threshold;
    int fillin;
    int restart;
    int verbose;
    int maxiter;
    double tolerance;
};

/**
 * \class BackendGmm
 *
 * this class provides an interface to the GMM linear algebra library
 */
template<typename T>
class BackendGmm : public Backend<T>
{
    typedef Backend<T> super;
public:

    // -- TYPEDEFS --
    typedef typename super::value_type value_type;

    /* defaults */
    typedef BackendGmmDefaults defaults_type;

    /* matrix */
    typedef typename super::sparse_matrix_type sparse_matrix_type;
    typedef typename super::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef MatrixGmm<value_type, gmm::row_major> gmm_sparse_matrix_type;
    typedef boost::shared_ptr<gmm_sparse_matrix_type> gmm_sparse_matrix_ptrtype;

    typedef typename sparse_matrix_type::graph_type graph_type;
    typedef typename sparse_matrix_type::graph_ptrtype graph_ptrtype;

    /* vector */
    typedef typename super::vector_type vector_type;
    typedef typename super::vector_ptrtype vector_ptrtype;
    typedef VectorUblas<value_type> gmm_vector_type;
    typedef boost::shared_ptr<vector_type> gmm_vector_ptrtype;

    typedef typename super::solve_return_type solve_return_type;
    typedef typename super::nl_solve_return_type nl_solve_return_type;

    // -- CONSTRUCTOR --
    BackendGmm();

    BackendGmm( po::variables_map const& vm, std::string const& prefix = "" );

    // -- FACTORY METHODS --
    template<typename DomainSpace, typename DualImageSpace>
    static sparse_matrix_ptrtype newMatrix( boost::shared_ptr<DomainSpace> const& space1,
                                            boost::shared_ptr<DualImageSpace> const& space2,
                                            size_type matrix_properties = NON_HERMITIAN )
    {
        auto A= sparse_matrix_ptrtype( new gmm_sparse_matrix_type( space1->nDof(), space2->nDof() ) );
        A->setMatrixProperties( matrix_properties );
        return A;
    }

    sparse_matrix_ptrtype
    newMatrix(const size_type m,
              const size_type n,
              const size_type m_l,
              const size_type n_l,
              const size_type nnz=30,
              const size_type noz=10,
              size_type matrix_properties = NON_HERMITIAN )
    {
        sparse_matrix_ptrtype mat( new gmm_sparse_matrix_type(m,n) );
        mat->setMatrixProperties( matrix_properties );
        return mat;
    }

    sparse_matrix_ptrtype
    newMatrix(const size_type m,
              const size_type n,
              const size_type m_l,
              const size_type n_l,
              graph_ptrtype const & graph,
              size_type matrix_properties = NON_HERMITIAN)
    {
        sparse_matrix_ptrtype mat( new gmm_sparse_matrix_type(m,n) );
        mat->setMatrixProperties( matrix_properties );
        return mat;
    }

    sparse_matrix_ptrtype
    newMatrix( DataMap const& d1, DataMap const& d2, size_type matrix_properties = NON_HERMITIAN, bool init = true)
    {
        auto A = sparse_matrix_ptrtype( new gmm_sparse_matrix_type( d1.nGlobalElements(), d2.nGlobalElements() ) );
        A->setMatrixProperties( matrix_properties );
        return A;
    }

    sparse_matrix_ptrtype
    newZeroMatrix( const size_type m,
                   const size_type n,
                   const size_type m_l,
                   const size_type n_l )
    {
        auto A = sparse_matrix_ptrtype( new gmm_sparse_matrix_type( m, n ) );
        //A->setMatrixProperties( matrix_properties );
        return A;
    }


    sparse_matrix_ptrtype
    newZeroMatrix( DataMap const& d1, DataMap const& d2 )
    {
        auto A = sparse_matrix_ptrtype( new gmm_sparse_matrix_type( d1.nGlobalElements(), d2.nGlobalElements() ) );
        //A->setMatrixProperties( matrix_properties );
        return A;
    }

    template<typename SpaceT>
    static vector_ptrtype newVector( boost::shared_ptr<SpaceT> const& space )
    {
        return vector_ptrtype( new gmm_vector_type( space->nDof() ) );
    }

    template<typename SpaceT>
    static vector_ptrtype newVector( SpaceT const& space )
    {
        return vector_ptrtype( new gmm_vector_type( space.nDof() ) );
    }

    vector_ptrtype
    newVector( DataMap const& d )
    {
        return vector_ptrtype( new gmm_vector_type( d.nGlobalElements() ) );
    }

    vector_ptrtype newVector( const size_type n, const size_type n_local )
    {
        return vector_ptrtype( new gmm_vector_type( n ) );
    }


    // -- SETTING OF OPTIONS --
    void set_noisy    ( int noisy )        { M_iter.set_noisy( noisy );     }
    void set_maxiter  ( int maxiter )      { M_iter.set_maxiter( maxiter ); }
    void set_fillin   ( int fillin )       { M_fillin    = fillin;          }
    void set_threshold( double threshold ) { M_threshold = threshold;       }
    void set_tol      ( double tol )       { M_iter.set_resmax( tol );      }

    void set_symmetric( bool isSymmetric )
    {
        M_isSymmetric = isSymmetric;
        if ( isSymmetric )
        {
            M_ilut.reset( 0 );
            M_ilutp.reset( 0 );
        }
    }

    void set_direct( bool isDirect )
    {
#if defined(GMM_USES_SUPERLU) || defined(FEELPP_HAS_UMFPACK)
        M_isDirect = isDirect;
        if ( isDirect )
        {
            M_solver_type = "umfpack";
            M_ilut.reset( 0 );
            M_ilutp.reset( 0 );
        }
#else
        if ( isDirect )
            std::cerr << "[BackendGmm] direct solver is not available\n";
#endif
    }

    void set_solver_type( std::string const& solver )
    {
        M_solver_type = solver;
    }

    void set_preconditioner_type( std::string const& prec )
    {
        M_precond_type = prec;
    }

    void set_restart( int restart )
    {
        M_restart = restart;
    }

    // -- LINEAR ALGEBRA INTERFACE --
    void prod( sparse_matrix_type const& A,
               vector_type const& x,
               vector_type& b ) const
    {
        gmm_sparse_matrix_type const& _A = dynamic_cast<gmm_sparse_matrix_type const&>( A );
        gmm_vector_type const& _x = dynamic_cast<gmm_vector_type const&>( x );
        gmm_vector_type& _b = dynamic_cast<gmm_vector_type&>( b );
        gmm::mult( _A.mat(), _x.vec(), _b.vec());
    }

    solve_return_type solve( sparse_matrix_type const& A,
                             vector_type& x,
                             const vector_type& b );

    solve_return_type solve( sparse_matrix_ptrtype const& A,
                             vector_ptrtype& x,
                             const vector_ptrtype& b )
    {
        return this->solve( *A, *x, *b );
    }

    solve_return_type solve( sparse_matrix_ptrtype const& A,
                             sparse_matrix_ptrtype const& P,
                             vector_ptrtype& x,
                             const vector_ptrtype& b )
    {
        return this->solve( *A, *x, *b );
    }


    value_type dot( const gmm_vector_type& f,
                    const gmm_vector_type& x ) const
    {
        value_type result(0);
        typename gmm_vector_type::const_iterator fi;
        typename gmm_vector_type::const_iterator xi;
        for( fi=f.begin(), xi=x.begin();
             ( fi!=f.end() ) && ( xi!=x.end() );
             ++fi, ++xi )
            {
                result += (*xi) * (*fi);
            }
        return result;
    }

    // -- GETTING DETAILS ABOUT SOLVING --
    bool converged() { return M_iter.converged(); }

    size_type get_iteration() { return M_iter.get_iteration(); }

    boost::shared_ptr<MatrixTriplet<T> > toTriplet( sparse_matrix_type const& m );


private:

    /* typedefs preconditioners */
    typedef gmm::ilut_precond<typename gmm_sparse_matrix_type::matrix_type>
    ilut_type;
    typedef gmm::ilutp_precond<typename gmm_sparse_matrix_type::matrix_type>
    ilutp_type;
    typedef std::auto_ptr<ilut_type> ilut_ptrtype;
    typedef std::auto_ptr<ilutp_type> ilutp_ptrtype;
    typedef gmm::diagonal_precond<typename gmm_sparse_matrix_type::matrix_type>
    diagprec_type;

    /* private member data */
    gmm::iteration M_iter;
    int            M_fillin;
    double         M_threshold;
    bool           M_isSymmetric;
    bool           M_isDirect;
    ilut_ptrtype   M_ilut;
    ilutp_ptrtype  M_ilutp;
    std::string    M_solver_type;
    std::string    M_precond_type;
    int            M_restart;

#if defined( FEELPP_HAS_UMFPACK )
    SolverUMFPACK  M_umfpack;
#endif // FEELPP_HAS_UMFPACK

}; // class BackendGmm



po::options_description backendgmm_options( std::string const& prefix = "",
                                            BackendGmmDefaults defaults = BackendGmmDefaults() );

} // Feel

#endif /* _BACKENDGMM_HPP_ */
