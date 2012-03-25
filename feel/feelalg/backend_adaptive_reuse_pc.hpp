/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2006-11-16

  Copyright (C) 2006 EPFL

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
   \file backend_adaptive_reuse_pc.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2006-11-16
 */

#ifndef _BACKEND_ADAPTIVE_REUSE_PC_HPP_
#define _BACKEND_ADAPTIVE_REUSE_PC_HPP_

#include <boost/program_options/variables_map.hpp>
#include <boost/timer.hpp>

namespace Feel
{

// struct to hold default values for BackendAdaptiveReusePC
struct BackendAdaptiveReusePCdefaults
{
    // Default constructor with default default values
    BackendAdaptiveReusePCdefaults()
        :
        maxiter( 1000 )
    {}

    int maxiter;
};

template<class BackendBase>
class BackendAdaptiveReusePC : public BackendBase
{

public:

    // -- TYPEDEFS --
    typedef BackendBase backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename backend_type::value_type value_type;

    /* matrix */
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;

    /* vector */
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    // -- CONSTRUCTOR --
    BackendAdaptiveReusePC() :
        M_backend(),
        M_totalSolveIter(0.),
        M_lastSolveIter(0.),
        M_firstSolveTime(0.),
        M_nUsePC(0),
        M_reusePC(false),
        M_reusedPC(false),
        M_reuseFailed(false),
        M_maxiter(1000),
        M_iteration(0)
    {
    }

    BackendAdaptiveReusePC( po::variables_map const& vm, std::string const& prefix = "" ) :
        M_backend( backend_type::build( vm, prefix ) ),
        M_totalSolveIter(0.),
        M_lastSolveIter(0.),
        M_firstSolveTime(0.),
        M_nUsePC(0),
        M_reusePC(false),
        M_reusedPC(false),
        M_reuseFailed(false),
        M_maxiter(1000),
        M_iteration(0)
    {
        std::string _prefix = prefix;
        if ( !_prefix.empty() )
            _prefix += "-";
        M_maxiter = vm[_prefix+"arpc-maxiter"].template as<int>();
    }

    // -- FACTORY METHODS --
    template<typename DomainSpace, typename DualImageSpace>
    sparse_matrix_ptrtype newMatrix( DomainSpace const& domainSpace,
                                     DualImageSpace const& dualImageSpace)
    {
        return M_backend->newMatrix( domainSpace, dualImageSpace );
    }

    template<typename SpaceT>
    vector_ptrtype newVector( SpaceT const& space )
    {
        return M_backend->newVector( space );
    }

    // -- SETTING OF OPTIONS --
    void set_noisy( int noisy )
    {
        //M_backend->set_noisy( noisy );
    }
    void set_maxiter( int maxiter ) { M_maxiter = maxiter; }
    void set_fillin( int fillin )
    {
        //M_backend->set_fillin( fillin );
    }
    void set_threshold( double threshold )
    {
        //M_backend->set_threshold( threshold );
    }
    void set_tol( double tol )
    {
        //M_backend->set_tol( tol );
    }
    void set_symmetric( bool isSymmetric )
    {
        M_backend->set_symmetric( isSymmetric );
    }
    void set_direct( bool isDirect )
    {
        M_backend->set_direct( isDirect );
    }
    void set_solver_type( std::string const& solver )
    {
        //M_backend->set_solver_type( solver );
    }
    void set_preconditioner_type( std::string const& prec )
    {
        M_backend->set_preconditioner_type( prec );
    }
    void set_restart( int restart )
    {
        M_backend->set_restart( restart );
    }

    // -- LINEAR ALGEBRA INTERFACE --
    template <class Vector>
    inline void applyMatrix( sparse_matrix_type const& A,
                             const Vector& x,
                             vector_type& b )
    {
        M_backend->applyMatrix( A, x, b );
    }

    // if reusePC specified, bypass adaptive behaviour
    template <class Vector>
    inline void solve( sparse_matrix_type const& A,
                       Vector& x,
                       const vector_type& b,
                       bool reusePC )
    {
        if ( !reusePC ) {
            reset();
        }
        start();
        M_backend->solve( A, x, b, reusePC );
        stop();
        M_iteration = M_backend->get_iteration();
        M_reuseFailed = M_reusedPC && !M_backend->converged();
    }

    // adaptive version
    template <class Vector>
    void solve( sparse_matrix_type const& A,
                Vector& x,
                const vector_type& b )
    {
        if ( !M_reusePC ) {
            reset();
        }
        start();
        M_backend->solve( A, x, b, M_reusePC );
        stop();
        M_iteration = M_backend->get_iteration();
        M_reuseFailed = M_reusedPC && !M_backend->converged();
        if ( M_reuseFailed ) {
            reset();
            start();
            M_backend->solve( A, x, b, false );
            stop();
            M_iteration += M_backend->get_iteration();
        }
    }

    template <class Vector>
    inline value_type dot( const vector_type& f,
                                  const Vector& x )
    {
        return M_backend->dot( f, x );
    }

    // -- GETTING DETAILS ABOUT SOLVING --
    inline bool converged()
    {
        return M_backend->converged();
    }

    inline size_type get_iteration()
    {
        return M_iteration;
    }

    inline bool reusePC() const
    {
        return M_reusePC;
    }

    inline bool reusedPC() const
    {
        return M_reusedPC;
    }

    inline bool reuseFailed() const
    {
        return M_reuseFailed;
    }

private:

    void start() {
        M_timer.restart();
    }

    void stop()  {
        double solveTime = M_timer.elapsed();
        double solveIter = M_backend->get_iteration() + 0.01;
        M_reusedPC = M_reusePC;
        ++M_nUsePC;
        if ( M_nUsePC == 1 ) {
            M_reusePC = true;
            M_firstSolveTime = solveTime;
            if ( !M_reuseFailed )
                M_backend->set_maxiter( std::min( M_maxiter,
                                                 (int)(1.5*solveIter + 10.5)) );
        } else {
            double nextSolveIter;
            if ( M_nUsePC == 2 ) {
                M_totalSolveIter = solveIter*(1.0+M_firstSolveTime/solveTime);
                nextSolveIter = solveIter;
            } else {
                M_totalSolveIter += solveIter;
//                 if ( solveIter > M_lastSolveIter )
//                     nextSolveIter = 2*solveIter - M_lastSolveIter;
//                 else
//                     nextSolveIter = solveIter * solveIter / M_lastSolveIter;
                nextSolveIter = solveIter;
            }
            M_reusePC = ( M_totalSolveIter > M_nUsePC * nextSolveIter );
            M_lastSolveIter = solveIter;
            if ( M_reusePC )
            {
                M_backend->set_maxiter( std::min
                                       ( M_maxiter,
                                         (int)( M_totalSolveIter/M_nUsePC
                                                + 0.5 ) ) );
            }
            else
            {
                M_backend->set_maxiter( M_maxiter );
            }
        }
    }

    void reset() {
        M_reusePC = false;
        M_totalSolveIter = 0.0;
        M_nUsePC = 0;
        M_backend->set_maxiter( M_maxiter );
    }

    backend_ptrtype M_backend;

    double M_totalSolveIter;
    double M_lastSolveIter;
    double M_firstSolveTime;
    size_t M_nUsePC;
    bool   M_reusePC;
    bool   M_reusedPC;
    bool   M_reuseFailed;
    boost::timer M_timer;
    int    M_maxiter;
    int    M_iteration;

}; // class BackendAdaptiveReusePC

po::options_description backend_adaptive_reuse_pc_options( std::string const& prefix = "",
                                                           BackendAdaptiveReusePCdefaults defaults = BackendAdaptiveReusePCdefaults() );
} // Feel

#endif /* _BACKEND_ADAPTIVE_REUSE_PC_HPP_ */
