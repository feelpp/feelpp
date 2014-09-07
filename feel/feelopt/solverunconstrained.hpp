/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-14

  Copyright (C) 2008 Universite Joseph Fourier (Grenoble I)

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
   \file solverunconstrained.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
#ifndef SOLVERUNCONSTRAINED_H
#define SOLVERUNCONSTRAINED_H

#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/if.hpp>

#include <clapack.h>

#include <feel/feelcore/feel.hpp>
#include <feel/feelopt/problem.hpp>
#include <feel/feelopt/dirscalingmatrix.hpp>

namespace Feel
{
/**
 * \class SolverUnconstrained
 *
 *  \author Ivan Oliveira and Christophe Prud'homme
 */
template<
typename Data,
         template<class> class Problem = problem
         >
class SolverUnconstrained
{
public:

    //! this solver type
    typedef SolverUnconstrained<Data,Problem> solver_type;

    //! problem data type
    typedef Problem<Data> problem_type;



    enum
    {
        _E_n   = problem_type::_E_n,  //!< number of control variables
        _E_g   = 0,                   //!< number of inequality constraints
        _E_h   = 0,                   //!< number of equality constraints
        _E_nA  = problem_type::_E_nA,  //!< size of the matrix
        _E_nL  = problem_type::_E_nL, //!< number of multipliers
        _E_nAL = problem_type::_E_nAL //!< size of the multipliers matrix
    };

    //! automatic differentiation type of order one
    typedef typename problem_type::ad_0_type ad_0_type;

    //! automatic differentiation type of order one
    typedef typename problem_type::ad_1_type ad_1_type;

    //! automatic differentiation type of order two
    typedef typename problem_type::ad_2_type ad_2_type;

    //! numerical type
    typedef typename ad_0_type::value_type value_type;

    //! functional objective type
    typedef typename problem_type::f_type f_type;

    //! numerical type
    //typedef typename problem_type::value_type value_type;

    //! vector type
    typedef typename problem_type::vector_type vector_type;

    //! matrix type
    typedef typename problem_type::matrix_type matrix_type;
    typedef typename problem_type::symmetric_matrix_type symmetric_matrix_type;
    typedef ublas::banded_matrix<double> banded_matrix_type;

    typedef directionalScalingMatrix<value_type> dsm_type;
    /*!
      \class COptions
      \brief parameters for the solver
    */
    struct COptions
    {
        COptions()
            :
            Delta_init( 1.0 ),
            zeta_min( 0.9 ),
            CGtol( 1e-12 ),
            deps( 1e-9 ),
            max_TR_iter( 100 ),
            TR_tol( 1e-7 ),
            allow_Trust_radius_calculations( true ),
            rho_decrease( 0.3 ),
            rho_big( 0.9 ),
            rho_increase_big( 1.5 ),
            rho_small( 0.3 ),
            rho_increase_small( 1.1 )
        {
            // do nothing here
        }
        value_type Delta_init;
        value_type zeta_min;
        value_type CGtol;
        value_type deps;
        value_type max_TR_iter;
        value_type TR_tol;
        bool allow_Trust_radius_calculations;

        value_type rho_decrease;
        value_type rho_big;
        value_type rho_increase_big;
        value_type rho_small;
        value_type rho_increase_small;
    };

    //! option data type
    typedef COptions options_type;

    /*!
      \class Stats
      \brief Statistics holder for the solver
    */
    class Stats
    {
    public:

        /**
         * \brief default constructor
         * \param \c __c true if collect statistics, false otherwise
         */
        Stats( solver_type* __s, bool __c = false )
            :
            M_s( __s ),
            M_collect( __c ),
            M_x ( _E_n ),
            M_l ( _E_n ),
            M_u ( _E_n )
        {
            // nothing to do here
        }

        //! tells if we collect statistics regarding the solver convergence
        bool collectStats() const
        {
            return M_collect;
        }

        //! reset all statistics
        void clear();

        void push( const vector_type & __x, vector_type const&, vector_type const& );

        void push( value_type norm_Tgrad_fx,
                   int n_CGiter, int n_restarts, int n_indef,
                   int n_crosses_def, int n_crosses_indef,
                   int n_truss_exit_def, int n_truss_exit_indef,
                   value_type Delta, value_type ared_til, value_type phi_til, value_type rho );
        void show() const;

    private:

        solver_type *M_s;
        bool M_collect;
        mutable int iter;

        vector_type M_x;
        vector_type M_l;
        vector_type M_u;

        std::vector<value_type> norm_Tgrad_fx_hstr, norm_err;
        std::vector<int> n_CGiter_hstr, n_restarts_hstr, n_indef_hstr,
            n_crosses_def_hstr, n_crosses_indef_hstr,
            n_truss_exit_def_hstr, n_truss_exit_indef_hstr;
        std::vector<value_type> Delta_hstr, ared_til_hstr, phi_til_hstr, rho_hstr;



    };

    //! statistics data type
    typedef Stats statistics_type;

    //! constructor that takes the control variables bounds
    SolverUnconstrained( value_type x_definitions[_E_n][3] );

    //! constructor that takes a data instance
    SolverUnconstrained( Data const& __data );

    //! destructor
    ~SolverUnconstrained();

    //! redefined the bounds of the control variables
    void redefine_problem( value_type x_definitions[_E_n][3] )
    {
        M_prob.define_problem( x_definitions );
    }

    //! optimization algorithm
    bool optimize( vector_type& __x );

    statistics_type const& stats() const
    {
        return M_solver_stats;
    }

    problem_type& problem()
    {
        return M_prob;
    }

private:

    // make it private to avoid having it called
    SolverUnconstrained();


    //! read the options from a file
    void read_options();

    //!
    void makeCauchyStep( vector_type & _x, value_type _Delta,
                         f_type&,
                         vector_type & _Tgrad_fx,
                         banded_matrix_type & _Hg,
                         symmetric_matrix_type & _Thess_fxT, symmetric_matrix_type & _Htil,
                         vector_type & _neg_grad_fx );

    //!
    void makeStep( vector_type & _x, vector_type & _s, value_type _Delta,
                   f_type&,
                   vector_type & _Tgrad_fx,
                   banded_matrix_type & _Hg,
                   symmetric_matrix_type & _Thess_fxT, symmetric_matrix_type & _Htil );

    //!
    value_type norm_Theta_x_grad_fx( vector_type & _x, value_type _Delta );


    //!
    void CGstep( vector_type & _x, value_type _Delta, vector_type & _sCG, value_type &norm_s_til,
                 int &_CGiter,
                 int &_n_restarts, int &_n_indef,
                 int &_n_crosses_def, int &_n_crosses_indef,
                 int &_n_truss_exit_def, int &_n_truss_exit_indef,
                 value_type &_s_til_x_G_til_x_s_til, value_type &phi_til );

    /*!
      Find \f$\tau \geq 0\f$ such that:

      \f$||\tilde{s} + \tau d|| = \Delta\f$

      in fact it is a second order polynomial in \f$\tau\f$
      and we must find its roots
    */
    value_type tau( vector_type const& _s, vector_type const& _d, value_type const& _Delta );

    /*!
      \brief computes \f$\xi\f$ needed for STEPS 1 and 2.

      \f$\xi\f$ is defined as follows:
      \f$ \xi = \mathrm{min}\{ \xi_1, -(\frac{\tilde{s}}/{d})_i : \frac{\tilde{s}}/{d})_i < 0 \}\f$

      \param s vector at numerator of the fraction
      \param d vector at the denominator of the fraction
      \param xi1 value to compare to for \f$\xi\f$
    */
    value_type xi( vector_type const& __s, vector_type const& __d, value_type const& __xi1 );


    //!
    void lambda_LS( vector_type & _x, vector_type & _lambda_l, vector_type & _lambda_u );

private:

    // problem specification and data
    problem_type M_prob;

    options_type M_options;

    statistics_type M_solver_stats;

    dsm_type M_theta;

};

// SolverUnconstrained

//#define DEBUG_OUT_CG
//#define DEBUG_OUT_TR

template<typename Data,template<class> class Problem>
SolverUnconstrained<Data,Problem>::SolverUnconstrained( value_type x_definitions[_E_n][3] )
    :
    M_prob( x_definitions ),
    M_solver_stats ( this, true ),
    M_theta( _E_n )
{
    read_options();
}
template<typename Data,template<class> class Problem>
SolverUnconstrained<Data,Problem>::SolverUnconstrained( Data const& __data )
    :
    M_prob( __data ),
    M_solver_stats ( this, true ),
    M_theta( M_prob.lowerBounds(), M_prob.upperBounds() )
{
    read_options();
}

template<typename Data,template<class> class Problem>
SolverUnconstrained<Data,Problem>::~SolverUnconstrained()
{
}

template<typename Data,template<class> class Problem>
bool
SolverUnconstrained<Data,Problem>::optimize( vector_type& x )
{
    M_solver_stats.clear();

    // Controlling parameters
    // trust region radius
    value_type Delta = M_options.Delta_init;

    vector_type x_new ( x.size() );
    vector_type stot ( x.size() );
    value_type norm_Tgrad_fx = this->norm_Theta_x_grad_fx( x, Delta );
    value_type norm_s_til = 0;

    M_prob.copy_x0_to_x( x );


    // FIXME: if ( M_options.verbose() )
    //M_prob.print_complete( x );

    if ( M_solver_stats.collectStats() )
        M_solver_stats.push( norm_Tgrad_fx, 0, 0, 0, 0, 0, 0, 0, Delta, 0, 0, 0 );


    try
    {
        int iter = 0;

        //
        // Apply bound constrained Trust region algorithm
        //
        while ( iter == 0 ||
                ( iter < M_options.max_TR_iter && norm_Tgrad_fx > M_options.TR_tol ) )
        {
            iter++;
            int n_CGiter, n_restarts, n_indef,
                n_crosses_def, n_crosses_indef, n_truss_exit_def, n_truss_exit_indef;
            value_type _s_til_x_G_til_x_s_til, phi_til, rho, Delta_used = Delta;

            DVLOG(2) << "\n===================== iter = " << iter << " ===========================";
            //DVLOG(2) << "\nx = " << x << "\n";
            DVLOG(2) << "\n -> norm_Tgrad_fx = " << norm_Tgrad_fx;



            /** find an approximate stot for the step to make
             * solve :
             * Find \f$\tilde{s}^k   = \mathrm{arg}\mathrm{min}_{s \in R^n}{\tilde{\phi}^k : ||s|| < \Delta^k}\f$
             */
            this->CGstep( x, Delta, stot, norm_s_til,
                          n_CGiter,
                          n_restarts, n_indef,
                          n_crosses_def, n_crosses_indef,
                          n_truss_exit_def, n_truss_exit_indef,
                          _s_til_x_G_til_x_s_til, phi_til );

            //
            x_new = x + stot;

            f_type __fx_new;
            M_prob.evaluate( x_new, __fx_new, diff_order<0>() );
            f_type __fx;
            M_prob.evaluate( x, __fx, diff_order<0>() );

            // compute actual merit function reduction
            value_type ared_til = __fx_new.value( 0 ) - __fx.value( 0 ) + 0.5 * _s_til_x_G_til_x_s_til;

            rho = ared_til / phi_til;

            /////////////////////////////////////////////////////////////////////////
            // Trust region radius calculation:
            if ( M_options.allow_Trust_radius_calculations )
            {
                if ( rho <= 1e-8 )
                    Delta = M_options.rho_decrease * Delta;

                else
                {
                    x = x + stot;

                    if     ( rho > M_options.rho_big )
                        Delta = std::max( M_options.rho_increase_big*norm_s_til, Delta );

                    else if ( rho > M_options.rho_small )
                        Delta = std::max( M_options.rho_increase_small*norm_s_til, Delta );
                }

                norm_Tgrad_fx = this->norm_Theta_x_grad_fx( x, Delta );
            }

            else
            {
                x = x + stot;
                norm_Tgrad_fx = this->norm_Theta_x_grad_fx( x, Delta );
            }

            /////////////////////////////////////////////////////////////////////////

            if ( M_solver_stats.collectStats() )
            {
                M_solver_stats.push( norm_Tgrad_fx,
                                      n_CGiter, n_restarts, n_indef,
                                      n_crosses_def, n_crosses_indef,
                                      n_truss_exit_def, n_truss_exit_indef,
                                      Delta_used, ared_til, phi_til, rho );
            }

            if (  norm_Tgrad_fx > 1e-5 )
                M_prob.setAccuracy( std::min( 1e-1, norm_Tgrad_fx ) );

            DVLOG(2) << "norm_Tgrad_fx = " << norm_Tgrad_fx  << "\n";
        }
    }

    catch ( std::exception const& __ex )
    {
        f_type __fx_new;
        M_prob.evaluate ( x, __fx_new, diff_order<2>() );

        if ( norm_inf( __fx_new.gradient( 0 ) ) > 1e-10 )
            throw __ex;
    }

    vector_type __l( _E_n ), __u( _E_n );
    lambda_LS( x, __l, __u );

    if ( M_solver_stats.collectStats() )
        M_solver_stats.push( x, __l, __u );

    return true;
}
#if 0
template<typename Data,template<class> class Problem>    void
SolverUnconstrained<Data,Problem>::history()
{
    M_solver_stats.show();

    vector_type lambda_l ( _E_n ), lambda_u ( _E_n );
    this->lambda_LS( x, lambda_l, lambda_u );
    M_prob.print_stationary_x( x, lambda_l, lambda_u );
}
#endif
template<typename Data,template<class> class Problem>
void
SolverUnconstrained<Data,Problem>::read_options()
{
    value_type num;
    char str[20];
    std::ifstream fin( "t_options.inp" );

    if ( fin.fail() )
    {
        // no option file : give up
        return;
    }

    while ( !fin.eof() )
    {
        fin >> str >> num;

        if ( !strcmp( str, "Delta_init" ) )                      M_options.Delta_init = num;

        if ( !strcmp( str, "zeta_min" ) )                        M_options.zeta_min = num;

        if ( !strcmp( str, "CGtol" ) )                           M_options.CGtol = num;

        if ( !strcmp( str, "deps" ) )                            M_options.deps = num;

        if ( !strcmp( str, "max_TR_iter" ) )                     M_options.max_TR_iter = num;

        if ( !strcmp( str, "TR_tol" ) )                          M_options.TR_tol = num;

        if ( !strcmp( str, "allow_Trust_radius_calculations" ) ) M_options.allow_Trust_radius_calculations = num;

        if ( !strcmp( str, "rho_decrease" ) )                    M_options.rho_decrease = num;

        if ( !strcmp( str, "rho_big" ) )                         M_options.rho_big = num;

        if ( !strcmp( str, "rho_increase_big" ) )                M_options.rho_increase_big = num;

        if ( !strcmp( str, "rho_small" ) )                       M_options.rho_small = num;

        if ( !strcmp( str, "rho_increase_small" ) )              M_options.rho_increase_small = num;
    }

    fin.close();

}


template<typename Data,template<class> class Problem>
void
SolverUnconstrained<Data,Problem>::makeCauchyStep( vector_type & _x, value_type _Delta,
        f_type& __fx,
        vector_type & _Tgrad_fx,
        banded_matrix_type & _Hg,
        symmetric_matrix_type & _Thess_fxT, symmetric_matrix_type & _Htil,
        vector_type & _neg_grad_fx )
{
    M_prob.evaluate ( _x, __fx, diff_order<2>() );


    _neg_grad_fx = -  __fx.gradient( 0 );

    makeStep( _x, _neg_grad_fx, _Delta, __fx, _Tgrad_fx,
              _Hg, _Thess_fxT, _Htil );
}

template<typename Data,template<class> class Problem>
void
SolverUnconstrained<Data,Problem>::makeStep( vector_type & _x, vector_type & _s, value_type _Delta,
        f_type& __fx,
        vector_type & _Tgrad_fx,
        banded_matrix_type & _Hg,
        symmetric_matrix_type & _Thess_fxT, symmetric_matrix_type & _Htil )
{
    M_theta.update( _Delta, _x, _s );

    _Tgrad_fx = prod( M_theta(),  __fx.gradient( 0 ) );

    banded_matrix_type __diag_grad_fx( _E_n, _E_n, 0, 0 );

    matrix<value_type> __m( outer_prod( _Tgrad_fx,
                                        scalar_vector<value_type>( _Tgrad_fx.size(), 1 ) ) );
    __diag_grad_fx = banded_adaptor<matrix<value_type> >( __m, 0, 0 );

    _Hg = prod( M_theta.jacobian(), __diag_grad_fx );

    _Thess_fxT = prod( M_theta(), prod( __fx.hessian( 0 ), M_theta() ) );
    _Htil = _Thess_fxT + _Hg;
}

template<typename Data,template<class> class Problem>
typename SolverUnconstrained<Data,Problem>::value_type
SolverUnconstrained<Data,Problem>::norm_Theta_x_grad_fx( vector_type & _x, value_type _Delta )
{
    f_type __fx;
    M_prob.evaluate ( _x, __fx, diff_order<1>() );

    vector_type __Tgrad_fx( _E_n );
    vector_type __neg_grad_fx( _E_n );
    __neg_grad_fx = - __fx.gradient( 0 );

    M_theta.update( _Delta, _x, __neg_grad_fx, dsm_type::NO_JACOBIAN );

    __Tgrad_fx = prod ( M_theta(), __fx.gradient( 0 ) );

    return norm_2( __Tgrad_fx );
}

template<typename Data,template<class> class Problem>
typename SolverUnconstrained<Data,Problem>::value_type
SolverUnconstrained<Data,Problem>::tau( vector_type const& _s,
                                        vector_type const& _d,
                                        value_type const& _Delta )
{
    value_type a =  inner_prod( _d, _d );
    GST_SMART_ASSERT( a != 0 ).error( "a is 0 and should not be. It will cause a divide by zero." );
    value_type b = 2.0 * inner_prod( _s, _d );
    value_type c = inner_prod( _s, _s )  -  _Delta*_Delta;
    return ( -b + sqrt( b*b - 4.0*a*c ) ) / ( 2.0*a );
}

template<typename Data,template<class> class Problem>
typename SolverUnconstrained<Data,Problem>::value_type
SolverUnconstrained<Data,Problem>::xi( vector_type const& __s,
                                       vector_type const& __d,
                                       value_type const& _xi1 )
{

    vector_type __frac = - element_div ( __s, __d );

    namespace lambda = boost::lambda;

    value_type __xi =  _xi1;
    std::for_each( __frac.begin(), __frac.end(),
                   lambda::if_then( lambda::_1 > 0 && lambda::_1 < lambda::var( __xi ),
                                    lambda::var( __xi ) = lambda::_1 ) );
    return __xi;
}

template<typename Data,template<class> class Problem>
void
SolverUnconstrained<Data,Problem>::CGstep( vector_type & _x, value_type _Delta, vector_type & _sCG,
        value_type &norm_s_til,
        int &_CGiter,
        int &_n_restarts, int &_n_indef,
        int &_n_crosses_def, int &_n_crosses_indef,
        int &_n_truss_exit_def,
        int &_n_truss_exit_indef,
        value_type &_s_til_x_G_til_x_s_til,
        value_type &phi_til )
{
    bool _restart = false;
    bool _done = false;
    _n_restarts = 0;
    _CGiter = 0;
    _n_indef = 0;
    _n_crosses_def = 0;
    _n_crosses_indef = 0;
    _n_truss_exit_def = 0;
    _n_truss_exit_indef = 0;

    f_type __fx;

    vector_type _neg_grad_fx ( _x.size() );
    vector_type _Tgrad_fx ( _x.size() );
    vector_type _Htil_d ( _x.size() ); // Htil * dCG
    vector_type _Htil_s_til ( _x.size() ); // Htil * s_til

    vector_type _s_til ( _x.size() );
    vector_type _s_til_old ( _x.size() );
    vector_type _s_til_eps ( _x.size() );
    vector_type _rCG ( _x.size() );
    vector_type _rCG_old ( _x.size() );
    vector_type _dCG ( _x.size() );

    value_type _alpha = 0;
    value_type _gamma; // d' Htil d
    value_type _beta;
    value_type _xi = 0;
    value_type _tau;

    // diagonal matrices
    banded_matrix_type _Hg ( _x.size(), _x.size(), 0, 0 ); // Gtil

    symmetric_matrix_type _Thess_fxT ( _E_nA, _E_nA ); // Btil
    symmetric_matrix_type _Htil ( _E_nA, _E_nA );  // Htil

    DVLOG(2) << "\n\n[value_type SolverUnconstrained<Data,Problem>::CGstep]...\n";

    // INITIALIZE:
    this->makeCauchyStep( _x, _Delta, __fx, _Tgrad_fx, _Hg, _Thess_fxT, _Htil, _neg_grad_fx );

    DVLOG(2) << "Trust region active (C) : " << M_theta.isTrustRegionActive() << "\n";

    value_type _norm_gtil = norm_2( _Tgrad_fx );
    _s_til = zero_vector<value_type>( _s_til.size() );
    _s_til_old = _s_til;
    _rCG = -_Tgrad_fx;
    _rCG_old = _rCG;
    _dCG = _rCG;

    size_t _inner_iter = 0;

    _CGiter++;

    while ( !_done )
    {
        if ( _restart )
        {
            _n_restarts++;
            _CGiter++;
            // RE-INITIALIZE:
            _inner_iter = 0;

            _s_til_eps = _s_til + M_options.deps * _dCG;

            this->makeStep( _x, _s_til_eps, _Delta,
                            __fx,
                            _Tgrad_fx, _Hg, _Thess_fxT, _Htil );

            _norm_gtil = norm_2( _Tgrad_fx );

            _Htil_s_til = prod( _Htil, _s_til );

            _s_til_old = _s_til;
            _rCG = -_Tgrad_fx - _Htil_s_til;
            _rCG_old  = _rCG ;
            _dCG = _rCG;

            _restart = false;
        }

        //
        // STEP 1:
        //
        _gamma = inner_prod( _dCG, prod( _Htil,_dCG ) );

        if ( _gamma <= 0 )
        {
            _n_indef++;
            _tau = this->tau( _s_til, _dCG, _Delta );

            if ( ( !M_theta.isTrustRegionActive() ) || ( _CGiter == 1 ) )
            {
                _xi = _tau;
            }

            else if ( _inner_iter == 0 )
            {
                _xi = this->xi( _s_til_eps, _dCG, _tau );
            }

            else
            {
                _xi = this->xi( _s_til, _dCG, _tau );
            }

            _s_til_old = _s_til;
            _s_til += _xi * _dCG;

            if ( _xi < _tau )
            {
                _n_crosses_indef++;
                _restart = true;
            }

            else if ( M_theta.isTrustRegionActive() )
            {
                _sCG = _s_til;

                _n_truss_exit_indef++;
                _done = true;
            }
        }

        //
        // STEP 2:
        //
        if ( !_restart && !_done )
        {
            _alpha = inner_prod( _rCG, _rCG ) / _gamma;

            if ( ( !M_theta.isTrustRegionActive() ) || ( _CGiter == 1 ) )
                _xi = _alpha;

            else if ( _inner_iter == 0 )
                _xi = this->xi( _s_til_eps, _dCG, _alpha );

            else
                _xi = this->xi( _s_til, _dCG, _alpha );

            _s_til_old = _s_til;
            _s_til += _xi * _dCG;

            if ( norm_2( _s_til ) > _Delta )
            {
                _tau = this->tau( _s_til_old, _dCG, _Delta );

                _sCG = _s_til_old + _tau * _dCG;

                _n_truss_exit_def++;
                _done = true;
            }
        }

        //
        // STEP 3:
        //
        if ( !_restart && !_done )
        {
            if ( _xi >= _alpha )
            {
                _rCG_old = _rCG;
                _rCG -= _alpha * prod( _Htil,  _dCG );
            }

            else if ( M_theta.isTrustRegionActive() )
            {
                _n_crosses_def++;
                _restart = true;
            }
        }

        if ( norm_2( _rCG ) / _norm_gtil < M_options.CGtol )
        {
            _sCG = _s_til;

            DVLOG(2) << "\n\nNormal CG exit 1: ||rCG||/||g_til|| = " << norm_2( _rCG ) / _norm_gtil << "\n";

            _done = true;
        }

        if ( _CGiter >= _x.size() )
        {
            _sCG = _s_til;

            DVLOG(2) << "\n\nNormal CG exit 2: _CGiter = " << _CGiter << "\n";

            _done = true;
        }

        //
        // STEP 4:
        //
        if ( !_restart && !_done )
        {
            _beta = inner_prod( _rCG, _rCG ) / inner_prod( _rCG_old, _rCG_old );

            DVLOG(2) << "\nbeta = " << _beta << "\n";

            _dCG = _rCG + _beta * _dCG;

            _CGiter++;
            _inner_iter++;
        }
    }

    _s_til_x_G_til_x_s_til = inner_prod( _sCG, prod( _Hg, _sCG ) );

    phi_til = inner_prod( _Tgrad_fx, _sCG ) + 0.5 * inner_prod( _sCG, prod( _Htil,_sCG ) );

    norm_s_til = norm_2( _sCG );

    // restoring to original space
    _sCG = prod( M_theta(), _sCG );
}


template<typename Data, template<class> class Problem>
void
SolverUnconstrained<Data,Problem>::lambda_LS( vector_type & _x, vector_type & _lambda_l, vector_type & _lambda_u )
{
    //symmetric_matrix_type __A ( _E_nL, _E_nL );
    matrix_type __A ( _E_nL, _E_nL );
    vector_type __b ( _E_nL );

    // At = [ -I  (l-x)  0;  I  0  (x-u) ]  -->  AtA = [ I+(l-x)^2  -I;  -I  I+(x-u)^2 ]
    __A = zero_matrix<double>( _E_nL, _E_nL );

    for ( int __i = 0; __i < _E_n; ++__i )
    {
        __A( __i, __i ) = 1.0 + ( M_prob.x_l( __i )-_x ( __i ) )*( M_prob.x_l( __i )-_x ( __i ) );

        __A( __i, _E_n+__i ) = -1.0;
        __A( _E_n+__i,__i ) = -1.0;

        __A( ( _E_n+__i ),  _E_n+__i ) =  1.0+( _x ( __i )-M_prob.x_u( __i ) )*( _x ( __i )-M_prob.x_u( __i ) );
    }

    // b = [ \nabla f;  0;  0 ]  -->  Atb = [ -\nabla f; \nabla f ]
    f_type __f_x;
    M_prob.evaluate( _x, __f_x, diff_order<1>() );

    for ( int i = 0; i < _E_n; i++ )
    {
        value_type val = __f_x.gradient( 0, i );
        __b ( i ) =  val;
        __b ( _E_n+i ) = -val;
    }

    // after solve __b = [ _lambda_l;  _lambda_] = A\b
    //char __uplo = 'U';
    int __info = 0;
    int __nrhs = 1;
    int __N = _E_nL;
    //dppsv_ ( &__uplo, &__N, &__nrhs, __A.data(), __b.data(), &__ldb, &__info );
    //dppsv_ ( &__uplo, &__N, &__nrhs, __A.data().begin(), __b.data().begin(), &__N, &__info );

    ublas::vector<int> __itype( __N );
    //dspsv_( &__uplo, &__N, &__nrhs,  __A.data().begin(), __itype.data().begin(), __b.data().begin(), &__N, &__info );
    dgesv_( &__N,&__nrhs,__A.data().begin(),&__N,__itype.data().begin(),__b.data().begin(),&__N,&__info );

    if ( __info!= 0 )
    {
        std::ostringstream __err;

        __err << "[" << __PRETTY_FUNCTION__ << "] dppsv failed: " << __info << "\n"
              << "A = " << __A << "\n"
              << "B = " << __b << "\n";

        throw std::out_of_range( __err.str() );
    }

    for ( int i = 0; i < _E_n; i++ )
    {
        _lambda_l ( i ) = __b ( i );
        _lambda_u ( i ) = __b ( _E_n+i );
    }
}


//
// Statistics class implementation
//

template<typename Data, template<class> class Problem>
void
SolverUnconstrained<Data,Problem>::Stats::clear()
{
    norm_Tgrad_fx_hstr.clear();
    norm_err.clear();
    n_CGiter_hstr.clear();
    n_restarts_hstr.clear();
    n_indef_hstr.clear();
    n_crosses_def_hstr.clear();
    n_crosses_indef_hstr.clear();
    n_truss_exit_def_hstr.clear();
    n_truss_exit_indef_hstr.clear();
    Delta_hstr.clear();
    ared_til_hstr.clear();
    phi_til_hstr.clear();
    rho_hstr.clear();
}
template<typename Data, template<class> class Problem>
void
SolverUnconstrained<Data,Problem>::Stats::push( const vector_type & __x, vector_type const& __l, vector_type const& __u )
{
    M_x = __x;
    M_l = __l;
    M_u = __u;
}

template<typename Data, template<class> class Problem>
void
SolverUnconstrained<Data,Problem>::Stats::push( value_type norm_Tgrad_fx,
        int n_CGiter, int n_restarts, int n_indef,
        int n_crosses_def, int n_crosses_indef,
        int n_truss_exit_def, int n_truss_exit_indef,
        value_type Delta,
        value_type ared_til,
        value_type phi_til,
        value_type rho )
{
    //push_x_in_hstr( _x );
    norm_Tgrad_fx_hstr.push_back( norm_Tgrad_fx );
    n_CGiter_hstr.push_back( n_CGiter );
    n_restarts_hstr.push_back( n_restarts );
    n_indef_hstr.push_back( n_indef );
    n_crosses_def_hstr.push_back( n_crosses_def );
    n_crosses_indef_hstr.push_back( n_crosses_indef );
    n_truss_exit_def_hstr.push_back( n_truss_exit_def );
    n_truss_exit_indef_hstr.push_back( n_truss_exit_indef );
    Delta_hstr.push_back( Delta );
    ared_til_hstr.push_back( ared_til );
    phi_til_hstr.push_back( phi_til );
    rho_hstr.push_back( rho );
}

template<typename Data, template<class> class Problem>
void
SolverUnconstrained<Data,Problem>::Stats::show() const
{
    M_s->problem().print_stationary_x ( M_x, M_l, M_u );

    iter = norm_Tgrad_fx_hstr.size();
    std::cerr << "\n\nScaled Trust-Region Method Statistics:\n";
    std::cerr << "\nNumber of outter iterations: " << iter-1 << "\n";

    std::cerr << "\n";
    std::cerr << "  k  CGiters restarts indefs crosses trustEx   ||g_til||   Delta^k   ared_til^k  phi_til^k     rho^k \n";
    std::cerr << "--------------------------------------------------------------------------------------------------------\n";

    for ( int k = 0; k < iter; k++ )
    {
        fprintf( stderr, "%3d ", k );

        if ( k == 0 )
            std::cerr << "   N/A     N/A     N/A     N/A     N/A";

        else
        {
            std::cerr << "    " << n_CGiter_hstr[k]
                      << "       "  << n_restarts_hstr[k]
                      << "       "  << n_indef_hstr[k]
                      << "       "  << n_crosses_def_hstr[k] + n_crosses_indef_hstr[k]
                      << "       " << n_truss_exit_def_hstr[k] + n_truss_exit_indef_hstr[k] << " ";
        }

        fprintf ( stderr, " %13.3e %9f", norm_Tgrad_fx_hstr[ k ], Delta_hstr[k] );

        if ( k == 0 )
            std::cerr << "      N/A         N/A          N/A";

        else
        {
            fprintf ( stderr, " %11.3e", ared_til_hstr[ k ] );
            fprintf ( stderr, " %11.3e", phi_til_hstr[ k ] );
            fprintf ( stderr, " %11f", rho_hstr[ k ] );
        }

        std::cerr << std::endl;
    }

    //M_s->lambda_LS(
}

} // Feel
#endif
