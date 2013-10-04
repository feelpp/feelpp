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
   \file solverconstrained.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
#ifndef CONSTRAINED_SOLVER_HPP
#define CONSTRAINED_SOLVER_HPP

#include <feel/feelopt/problem.hpp>

namespace Feel
{
/*!
  \class SolverConstrained

  \author Ivan Oliveira and Christophe Prud'homme
*/
template<
typename Data,
         template<class> class Problem = problem
         >
class SolverConstrained
{
public:

    //! this solver type
    typedef SolverConstrained<Data,Problem> solver_type;

    //! problem data type
    typedef Problem<Data> problem_type;

    enum
    {
        _E_n   = problem_type::_E_n,  //!< number of control variables
        _E_f   = problem_type::_E_f,  //!< number of objective functions
        _E_g   = problem_type::_E_g,  //!< number of inequality constraints
        _E_h   = problem_type::_E_h,  //!< number of equality constraints
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

    typedef typename problem_type::f_type f_type;
    typedef typename problem_type::g_type g_type;
    typedef typename problem_type::h_type h_type;

    /*!
      \class COptions
      \brief holds the options for the SolverConstrained

    */
    struct COptions
    {
        COptions()
            :
            eta ( 1e-8 ),
            tau ( 0.995 ),
            theta_N ( 0.2 ),
            zeta_N ( 0.8 ),
            eTOL ( 1e-6 ),
            rho_N ( 0.3 ),

            Delta_init ( 1 ),
            s_theta_init ( 10 ),
            eMU_init ( 1e-1 ),
            MU_init ( 1e-1 ),
            nu0_init ( 10 ),

            rho_decrease ( 0.1 ),
            rho_big ( 0.9 ),
            rho_increase_big ( 5.0 ),
            rho_small ( 0.3 ),
            rho_increase_small ( 2.0 ),

            Primal_Dual ( true ),

            max_TR_iter ( 4 ),
            max_hom_iter ( 8 )

        {
            // nothing to do here
        }
        double eta;
        double tau;
        double theta_N;
        double zeta_N;
        double eTOL;
        double rho_N;

        double Delta_init;
        double s_theta_init;
        double eMU_init;
        double MU_init;
        double nu0_init;



        double rho_decrease;
        double rho_big;
        double rho_increase_big;
        double rho_small;
        double rho_increase_small;

        bool Primal_Dual;

        int max_TR_iter;
        int max_hom_iter;

    };

    //! option data type
    typedef COptions options_type;

    /*!
      \class Stats
      \brief holds the statistics for SolverConstrained
    */
    struct Stats
    {
        /**
         * \brief default constructor
         * \param \c __c true if collect statistics, false otherwise
         */
        Stats( solver_type* __s, bool __c = false )
            :
            M_s( __s ),
            M_collect( __c )
        {
            // nothing to do here
        }

        //! tells if we collect statistics regarding the solver convergence
        bool collectStats() const
        {
            return M_collect;
        }

        //! clear all statistics
        void clear();

        //! insert x in statistics
        void push_x( const double *__x );

        //! insert all statistics
        void push_all( double mu, double E, double E0,
                       int n_CGiter,
                       bool truss_exit, bool SOC,
                       double Delta, double ared, double pred, double gamma );

        /**
           \brief show statistics
           \param verbose if \c true \c show is more verbose
        */
        void show( bool verbose = false ) const;

    private:
        //! make it private to avoid being used
        Stats();

    private:


        solver_type* M_s;

        bool M_collect;

        mutable int iter;
        std::vector<double> x_hstr, mu_hstr, E_hstr, E0_hstr;
        std::vector<int> n_CGiter_hstr;
        std::vector<bool> truss_exit_hstr, SOC_hstr;
        std::vector<double> Delta_hstr, ared_hstr, pred_hstr, gamma_hstr;
    };

    //! statistics data type
    typedef Stats statistics_type;

    //! default constructor
    /*
      \brief default constructor
      takes the bounds for the control variables
      \todo change this in the future by adding a class for control variables
    */
    SolverConstrained( double x_definitions[_E_n][3] );

    //! destructor
    ~SolverConstrained();

    //! redefined the control variables bounds and if they are avtive
    void redefine_problem( double x_definitions[_E_n][3] )
    {
        M_prob.define_problem( x_definitions );
    }

    /*!
      \brief optimizes the problem defined by \c problem_type
      \return \c true if optimization succeeded, \c false otherwise
    */
    bool optimize( double __x[_E_n] );

    //! \return the problem specifications
    problem_type& problem()
    {
        return M_prob;
    }

    //! \return the statistics of the solver
    statistics_type const& stats() const
    {
        return M_solver_stats;
    }

private:

    //! read options from file if it exists
    void read_options();

    void initialize_solver_with_x( double *_x, double *_s, double *_lambda_h, double *_lambda_g, double MU );

    void update_opt_objs( double *_x, double *_s, double *_lambda_h,
                          double *_lambda_g, double MU, bool update_Hess );

    void solve_AtA_system( double _r[  _E_n + _E_g+2*_E_n ], double _b[ _E_h + _E_g+2*_E_n ],
                           double _g[ _E_n + _E_g+2*_E_n ], double _l[ _E_h + _E_g+2*_E_n ] );

    void dogleg_solve( double *_vCP_til, double *_vN_til, double *_s, double Delta,
                       double *_v_til );

    double model_m( double *_s, double *_v_til );

    double get_vpred( double *_s, double *_v_til );

    //! returns vpred( v_til )
    double get_normal_steps( double *_s, double Delta,
                             double *_dx_normal, double *_ds_normal, double *_v_til );

    void get_Pr( double *_r, double *_Pr );

    bool slack_feasible( double *_w_til, double *_v_til );

    double get_q( double *_v_til, double *_w_til, double MU );

    //! returns q( v_til + w_til )
    double get_tangen_steps( double *_s, double Delta, double *_v_til, double MU,
                             double *_dx_tangen, double *_ds_tangen,
                             bool &truss_exit, int &n_CGiter );

    double get_steps( double *_s, double _Delta, double MU,
                      double *_dx, double *_ds, double &_vpred, double &_q,
                      bool &truss_exit, int &n_CGiter );

    double get_phi( double *_x, double *_s, double MU, double nu, double *rhs );
    double get_phi( double *_x, double *_s, double MU, double nu );
    double get_phi( double *_s, double MU, double nu );

    double get_SOC( double *_s, double *_rhs_xs_p_d, double *_xy, double *_sy );

    //
    // Evaluations
    //

    //!
    void evaluate_Ag( double *x );
    //!
    void evaluate_Ah( double *x );
    //!
    void evaluate_A_hat( double *_x, double *_s );
    //!
    void evaluate_Hessian_L( double *_x, double *_lambda_h, double *_lambda_g );
    //!
    void evaluate_G( double *_x, double *_s, double *_lambda_h, double *_lambda_g,
                     double MU, bool evaluate__Hessians_flag );
    double evaluate_E_current( double *_s, double *_lambda_h, double *_lambda_g, double MU );
    //!
    void evaluate_gradf( double *_x );
    //!
    void evaluate_fxgxhx( double *_x );

    //
    // Print
    //

    //!
    void print_Ag();
    //!
    void print_Ah();
    //!
    void print_lambda_LS( double *_lambda_g, double *_lambda_h );
    //!
    void print_G();
    //!
    void print_A_hat();
    //!
    void print_Hessian_L();

    //! make it private to disable it (force developer to use other constructor)
    SolverConstrained();

private:

    problem_type M_prob;

    //! options data
    options_type options;

    //! statistics data
    statistics_type M_solver_stats;

    double fx, gradf[ _E_n ];
    double gx[ _E_g+2*_E_n ];
    double hx[ _E_h ];
    double Ag[ _E_n ][ _E_g+2*_E_n ];
    double Ah[ _E_n ][ _E_h ];
    double A_hat[ _E_n + _E_g+2*_E_n ][ _E_h + _E_g+2*_E_n ];

    double Hessian_L[ _E_n ][ _E_n ];
    double G[ _E_n+_E_g+2*_E_n ][ _E_n+_E_g+2*_E_n ];

    static const double Max_Allowed_Delta = 1e8;
};


// SolverConstrained

//#define DEBUG_OUT_CG
//#define DEBUG_OUT_TR
//#define DEBUG_OUT_HM
//#define DEBUG_x_OUT

template<typename Data,template<class> class Problem>
SolverConstrained<Data,Problem>::SolverConstrained( double x_definitions[_E_n][3] )
    :
    M_prob( x_definitions ),
    M_solver_stats( this, true ),
    fx( 0 ),
    gradf(),
    gx(),
    hx(),
    Ag(),
    Ah(),
    A_hat(),
    Hessian_L(),
    G()
{
    std::cerr << "SolverConstrained::SolverConstrained()\n";
    read_options();
}

template<typename Data,template<class> class Problem>
SolverConstrained<Data,Problem>::~SolverConstrained()
{
}

template<typename Data,template<class> class Problem>
void SolverConstrained<Data,Problem>::read_options()
{
    double num;
    char str[20];
    std::ifstream fin( "N_options.inp" );

    if ( fin.fail() )
        return;

    options.Primal_Dual = true;

    while ( !fin.eof() )
    {
        fin >> str >> num;

        if ( !strcmp( str, "Delta_init" ) )                      options.Delta_init = num;

        if ( !strcmp( str, "eta" ) )                             options.eta = num;

        if ( !strcmp( str, "tau" ) )                             options.tau = num;

        if ( !strcmp( str, "theta_N" ) )                         options.theta_N = num;

        if ( !strcmp( str, "zeta_N" ) )                          options.zeta_N = num;

        if ( !strcmp( str, "eTOL" ) )                            options.eTOL = num;

        if ( !strcmp( str, "rho_N" ) )                           options.rho_N = num;

        if ( !strcmp( str, "s_theta_init" ) )                    options.s_theta_init = num;

        if ( !strcmp( str, "eMU_init" ) )                        options.eMU_init = num;

        if ( !strcmp( str, "MU_init" ) )                         options.MU_init = num;

        if ( !strcmp( str, "nu0_init" ) )                        options.nu0_init = num;

        if ( !strcmp( str, "Primal_Dual" ) )                     options.Primal_Dual = ( bool )num;

        if ( !strcmp( str, "max_TR_iter" ) )                     options.max_TR_iter = ( int )num;

        if ( !strcmp( str, "max_hom_iter" ) )                    options.max_hom_iter = ( int )num;

        if ( !strcmp( str, "rho_decrease" ) )                    options.rho_decrease = num;

        if ( !strcmp( str, "rho_big" ) )                         options.rho_big = num;

        if ( !strcmp( str, "rho_increase_big" ) )                options.rho_increase_big = num;

        if ( !strcmp( str, "rho_small" ) )                       options.rho_small = num;

        if ( !strcmp( str, "rho_increase_small" ) )              options.rho_increase_small = num;
    }

    fin.close();

}

template<typename Data,template<class> class Problem>
void SolverConstrained<Data,Problem>::print_Hessian_L()
{
    int n = M_prob.n();

    for ( int i = 0; i < n; i++ )
    {
        printf( "\n" );

        for ( int j = 0; j < n; j++ )
            fprintf( stderr, "%10.3f", Hessian_L[i][j] );
    }

    std::cerr << "\n";
}

template<typename Data,template<class> class Problem>
void SolverConstrained<Data,Problem>::evaluate_Hessian_L( double *_x, double *_lambda_h, double *_lambda_g )
{
    f_type __fx;
    g_type __gx;
    h_type __hx;
    M_prob.evaluate ( _x, __fx, diff_order<2>() );
    M_prob.evaluate ( _x, __gx, diff_order<2>() );
    M_prob.evaluate ( _x, __hx, diff_order<2>() );

    for ( int i = 0; i < _E_n; i++ )
    {
        // symetry
        for ( int j = i; j < _E_n; j++ )
        {
            double val = 0;

            for ( int __k = 0; __k < _E_f; ++__k )
            {
                val += __fx.hessian( __k, i, j );
            }

            for ( int m = 0; m < _E_h; m++ )
            {
                val += _lambda_h[ m ] * __hx.hessian( m,i,j );
            }

            for ( int m = 0; m < _E_g; m++ )
            {
                val += _lambda_g[ m ] * __gx.hessian( m,i,j );
            }

            Hessian_L[ i ][ j ] = val;
            Hessian_L[ j ][ i ] = val;
        }
    }
}

template<typename Data,template<class> class Problem>
void SolverConstrained<Data,Problem>::print_G()
{
    int n = M_prob.n(), nv = n+_E_g+2*n;

    for ( int i = 0; i < nv; i++ )
    {
        printf( "\n" );

        for ( int j = 0; j < nv; j++ )
            fprintf( stderr, "%12.5f", G[i][j] );
    }

    std::cerr << "\n";
}

template<typename Data,template<class> class Problem>
void SolverConstrained<Data,Problem>::evaluate_G( double *_x, double *_s, double *_lambda_h,
        double *_lambda_g, double MU, bool evaluate_Hessians_flag )
{
    int n = M_prob.n();
    int nv = n+_E_g+2*n;

    if ( evaluate_Hessians_flag )
        this->evaluate_Hessian_L( _x, _lambda_h, _lambda_g );

    for ( int i = 0; i < nv; i++ )
        for ( int j = 0; j < nv; j++ )
            G[ i ][ j ] = 0;

    for ( int i = 0; i < n; i++ )
        for ( int j = 0; j < n; j++ )
            G[ i ][ j ] = Hessian_L[ i ][ j ];

    for ( int i = 0; i < _E_g+2*n; i++ )
        if ( options.Primal_Dual && _lambda_g[ i ] > 0 )
            G[ n+i ][ n+i ] = _lambda_g[ i ] * _s[ i ];

        else
            G[ n+i ][ n+i ] = MU;

    //this->print_G(); exit(0);

    /*
      std::cerr << "\nHL = ";
      this->print_Hessian_L();
      std::cerr << "\nlambda_g = ";
      __print_vec( _E_g+2*n, _lambda_g );
      std::cerr << "\ns = ";
      __print_vec( _E_g+2*n, _s );
      std::cerr << "\nG = ";
      this->print_G();
      exit(0);
    */
}


template<typename Data,template<class> class Problem>
void SolverConstrained<Data,Problem>::evaluate_Ag( double *_x )
{
    g_type __gx;
    M_prob.evaluate ( _x, __gx, diff_order<1>() );

    int n = M_prob.n();

    assert ( n == _E_n );

    // Ag = [  g_0 ... g_n  -I  I  ]
    for ( int i = 0; i < _E_n; i++ )
        for ( int m = 0; m < _E_g+2*_E_n; m++ )
            Ag[ i ][ m ] = 0;

    for ( int i = 0; i < _E_n; i++ )
    {
        for ( int m = 0; m < _E_g; m++ )
        {
            Ag[ i ][ m ] = __gx.gradient( m,i );
        }

        int __index = _E_g+i;
        Ag[ i ][ __index ]  = -1.0;
        Ag[ i ][ __index + n] = +1.0;
    }
}

template<typename Data,template<class> class Problem>
void SolverConstrained<Data,Problem>::evaluate_Ah( double *_x )
{
    h_type __hx;
    M_prob.evaluate ( _x, __hx, diff_order<1>() );

    int n = M_prob.n();

    // Ah = [  h_0 ... h_n  ]
    for ( int i = 0; i < n; i++ )
        for ( int m = 0; m < _E_h; m++ )
            Ah[ i ][ m ] = 0;

    for ( int i = 0; i < n; i++ )
        for ( int m = 0; m < _E_h; m++ )
            Ah[ i ][ m ] = __hx.gradient( m,i );
}


template<typename Data,template<class> class Problem>
void SolverConstrained<Data,Problem>::print_Ag()
{
    int n = M_prob.n();

    for ( int i = 0; i < n; i++ )
    {
        printf( "\n" );

        for ( int j = 0; j < _E_g+2*n; j++ )
            fprintf( stderr, "%10.1f", Ag[i][j] );
    }

    std::cerr << "\n";
}

template<typename Data,template<class> class Problem>
void SolverConstrained<Data,Problem>::print_Ah()
{
    int n = M_prob.n();

    for ( int i = 0; i < n; i++ )
    {
        printf( "\n" );

        for ( int j = 0; j < _E_h; j++ )
            fprintf( stderr, "%10.1f", Ah[i][j] );
    }

    std::cerr << "\n";
}

template<typename Data,template<class> class Problem>
void SolverConstrained<Data,Problem>::evaluate_fxgxhx( double *_x )
{
    f_type __fx;
    g_type __gx;
    h_type __hx;

    // get the order 0 information only
    M_prob.evaluate ( _x, __fx, diff_order<0>() );
    M_prob.evaluate ( _x, __gx, diff_order<0>() );
    M_prob.evaluate ( _x, __hx, diff_order<0>() );

    fx = __fx.value( 0 );

    for ( uint __i = 0; __i < __gx.value().size(); ++__i )
        gx[__i] = __gx.value( __i );

    for ( uint __i = 0; __i < _E_h; ++__i )
        hx[__i] = __hx.value( __i );
}


template<typename Data,template<class> class Problem>
void
SolverConstrained<Data,Problem>::initialize_solver_with_x( double *_x,
        double *_s,
        double *_lambda_h, double *_lambda_g,
        double MU )
{
    g_type __gx;
    M_prob.evaluate ( _x, __gx, diff_order<0>() );

    for ( uint __i = 0; __i < __gx.value().size(); ++__i )
    {
        gx[__i] = __gx.value( __i );
        _s[ __i ] = std::max( -gx[ __i ], options.s_theta_init );
    }

    /*
      for( int m = _E_g; m < ns; m++ )
      _s[ m ] = options.s_theta_init;
    */

    this->update_opt_objs( _x, _s, _lambda_h, _lambda_g, MU, true );

}

template<typename Data,template<class> class Problem>
void
SolverConstrained<Data,Problem>::evaluate_A_hat( double *_x, double *_s )
{
    int n = M_prob.n();

    this->evaluate_Ag( _x );
    this->evaluate_Ah( _x );

    // A_hat = [ Ah Ag ; 0 S ]
    for ( int i = 0; i < n+_E_g+2*n; i++ )
        for ( int m = 0; m < _E_h+_E_g+2*n; m++ )
            A_hat[ i ][ m ] = 0;

    for ( int i = 0; i < n; i++ )
    {
        for ( int m = 0; m < _E_h; m++ )
            A_hat[ i ][ m ] = Ah[ i ][ m ];

        for ( int m = 0; m < _E_g+2*n; m++ )
            A_hat[ i ][ _E_h+m ] = Ag[ i ][ m ];
    }

    for ( int m = 0; m < _E_g+2*n; m++ )
        A_hat[ n+m ][ _E_h+m ] = _s[ m ];

}

template<typename Data,template<class> class Problem>
void SolverConstrained<Data,Problem>::print_A_hat()
{
    int n = M_prob.n();

    for ( int i = 0; i < n+_E_g+2*n; i++ )
    {
        printf( "\n" );

        for ( int j = 0; j < _E_h+_E_g+2*n; j++ )
            fprintf( stderr, "%14.9f", A_hat[i][j] );
    }

    std::cerr << "\n";
}

template<typename Data,template<class> class Problem>
void SolverConstrained<Data,Problem>::solve_AtA_system( double _r[  _E_n + _E_g+2*_E_n ], double _b[ _E_h + _E_g+2*_E_n ],
        double _g[ _E_n + _E_g+2*_E_n], double _l[ _E_h + _E_g+2*_E_n ] )
{

    int n = M_prob.n();
    int n_A_hat = n + _E_g+2*n;
    int m_A_hat = _E_h + _E_g+2*n;
    int n_big = n_A_hat + m_A_hat;
    int big_N = _E_n + _E_g+2*_E_n + _E_h + _E_g+2*_E_n;

    double bigA[ big_N*big_N ], rhs[ big_N ];

    for ( int i = 0; i < n_big*n_big; i++ )
        bigA[ i ] = 0;

    for ( int i = 0; i < n_A_hat; i++ )
    {
        bigA[ i + i*n_big ] = 1;
        rhs[ i ] = _r[ i ];
    }

    double val;

    for ( int i = 0; i < n_A_hat; i++ )
    {
        for ( int j = 0; j < m_A_hat; j++ )
        {
            val = A_hat[ i ][ j ];
            bigA[ ( i ) + ( n_A_hat+j )*n_big ] = val;
            bigA[ ( n_A_hat+j ) + ( i )*n_big ] = val;
        }
    }

    for ( int j = 0; j < m_A_hat; j++ )
        rhs[ n_A_hat + j ] = -_b[ j ];

    /*
      this->print_A_hat();
      cerr << "\nA = [\n";
      for( int i = 0; i < n_big; i++ )
      {
      for( int j = 0; j < n_big; j++ )
      fprintf( stderr, "%10.1f", bigA[ i + j*n_big ]);
      cerr << "\n";
      }
      cerr << "]\n\nb = [";
      __print_vec( n_big, rhs );
      cerr << "]\n";
    */

    double sol[ big_N ];

#warning disabled Lusolve
    //__LUSolve( n_big, bigA, rhs, sol, STORAGE_GENERAL );

    for ( int i = 0; i < n_A_hat; i++ )
        _g[ i ] = sol[ i ];

    for ( int i = 0; i < m_A_hat; i++ )
        _l[ i ] = sol[ n_A_hat+i ];

}

template<typename Data,template<class> class Problem>
void SolverConstrained<Data,Problem>::evaluate_gradf( double *_x )
{
    f_type __fx;
    M_prob.evaluate ( _x, __fx, diff_order<1>() );

    for ( uint __i = 0; __i < _E_n; ++__i )
    {
        gradf[__i] = __fx.gradient( 0, __i );
    }
}

template<typename Data,template<class> class Problem>
void SolverConstrained<Data,Problem>::update_opt_objs( double *_x, double *_s,
        double *_lambda_h, double *_lambda_g,
        double MU, bool update_Hess )
{
    int n = M_prob.n();
    double _r[ _E_n + _E_g+2*_E_n ], _b[ _E_h + _E_g+2*_E_n ],
           _dummy[ _E_n + _E_g+2*_E_n ];
    double _lambda_LS[ _E_h + _E_g+2*_E_n ];

    this->evaluate_fxgxhx( _x );
    this->evaluate_gradf( _x );

    for ( int i = 0; i < n; i++ )
        _r[ i ] = -gradf[ i ];

    for ( int i = n; i < n+_E_g+2*n; i++ )
        _r[ i ] = MU;

    for ( int i = 0; i < _E_h + _E_g+2*n; i++ )
        _b[ i ] = 0;

    this->evaluate_A_hat( _x, _s );

    this->solve_AtA_system( _r, _b, _dummy, _lambda_LS );

    for ( int i = 0; i < _E_h; i++ )
        _lambda_h[ i ] = _lambda_LS[ i ];

    for ( int i = 0; i < _E_g+2*n; i++ )
        _lambda_g[ i ] = _lambda_LS[ _E_h+i ];

    this->evaluate_G( _x, _s, _lambda_h, _lambda_g, MU, update_Hess );

    /*
      std::cerr << "\nlambda_h = ";
      __print_vec( _E_h, _lambda_h );
      std::cerr << "\nlambda_g = ";
      __print_vec( _E_g+2*n, _lambda_g );
    */

}

template<typename Data,template<class> class Problem>
void SolverConstrained<Data,Problem>::print_lambda_LS( double *_lambda_h, double *_lambda_g )
{
    int n = M_prob.n();
    std::cerr << "\nlambda_h = ";
    __print_vec( _E_h, _lambda_h );
    std::cerr << "\nlambda_g = ";
    __print_vec( _E_g+2*n, _lambda_g );
    std::cerr << "\n";
}

template<typename Data,template<class> class Problem>
double SolverConstrained<Data,Problem>::evaluate_E_current( double *_s, double *_lambda_h, double *_lambda_g, double MU )
{
    int n = M_prob.n();
    double Ahlh[ _E_n ], Aglg[ _E_n ], dfAhlhAglg[ _E_n ], Slg_mue[ _E_g+2*_E_n], gps[ _E_g+2*_E_n];
    double E = 0;

    for ( int i = 0; i < n; i++ )
    {
        Ahlh[ i ] = 0;

        for ( int m = 0; m < _E_h; m++ )
            Ahlh[ i ] += Ah[ i ][ m ] * _lambda_h[ m ];

        Aglg[ i ] = 0;

        for ( int m = 0; m < _E_g+2*n; m++ )
            Aglg[ i ] += Ag[ i ][ m ] * _lambda_g[ m ];

        dfAhlhAglg[ i ] = gradf[ i ] + Ahlh[ i ] + Aglg[ i ];
    }

    E = std::max( E, __Norm_infty( n, dfAhlhAglg ) );

    for ( int m = 0; m < _E_g+2*n; m++ )
        Slg_mue[ m ] = _s[ m ] * _lambda_g[ m ] - MU;

    E = std::max( E, __Norm_infty( _E_g+2*n, Slg_mue ) );

    E = std::max( E, __Norm_infty( _E_h, hx ) );

    for ( int m = 0; m < _E_g+2*n; m++ )
        gps[ m ] = gx[ m ] + _s[ m ];

    E = std::max( E, __Norm_infty( _E_g+2*n, gps ) );

    /*
      this->print_Ah();
      this->print_Ag();
      this->print_A_hat();
      this->print_G();
      __print_vec( _E_g+2*n, _s );
      __print_vec( _E_g+2*n, _lambda_g );
      std::cerr << "\n|| df + Ahlh + Aglg ||_infty = " << __Norm_infty( n, dfAhlhAglg );
      std::cerr << "\nmu = " << MU;
      std::cerr << "\n|| Slg - mue ||_infty        = " << __Norm_infty( _E_g+2*n, Slg_mue );
      std::cerr << "\n|| hx ||_infty               = " << __Norm_infty( _E_h, hx );
      std::cerr << "\n|| gx + s ||_infty           = " << __Norm_infty( _E_g+2*n, gps );
      std::cerr << "\n";
    */


    return E;
}

template<typename Data,template<class> class Problem>
double SolverConstrained<Data,Problem>::model_m( double *_s, double *_v_til )
{
    int n = M_prob.n(), nv = n+_E_g+2*n;
    double A_hatTxv_til[ _E_h + _E_g+2*_E_n ], rhs[ _E_h + _E_g+2*_E_n ];

    for ( int i = 0; i < _E_h+_E_g+2*n; i++ )
    {
        A_hatTxv_til[ i ] = 0;

        for ( int j = 0; j < nv; j++ )
            A_hatTxv_til[ i ] += A_hat[ j ][ i ] * _v_til[ j ];
    }

    for ( int m = 0; m < _E_h; m++ )
        rhs[ m ] = hx[ m ];

    for ( int m = 0; m < _E_g+2*n; m++ )
        rhs[ _E_h + m ] = gx[ m ] + _s[ m ];

    double model = 2. * __Dot( _E_h+_E_g+2*n, rhs, A_hatTxv_til )
                   + __Dot( _E_h+_E_g+2*n, A_hatTxv_til, A_hatTxv_til );

    /*
      this->print_A_hat();
      __print_vec( nv, _v_til );
      __print_vec( _E_h+_E_g+2*n, A_hatTxv_til );
      std::cerr << "\n m = " << model;
      exit(0);
    */


    return model;
}

template<typename Data,template<class> class Problem>
void SolverConstrained<Data,Problem>::dogleg_solve( double *_vCP_til, double *_vN_til, double *_s, double _Delta,
        double *_v_til )
{
    int n = M_prob.n(), nv = n+_E_g+2*n;
    double v_diff[ _E_n+_E_g+2*_E_n ];
    double zeta_N = options.zeta_N, tau = options.tau;
    double theta1 = 1, theta2 = 1, theta3 = 1;

    theta1 = std::min( theta1, zeta_N * _Delta / __Norm_2( nv, _vN_til ) );

    for ( int i = n; i < nv; i++ )
        if ( _vN_til[ i ] < 0 )
            theta1 = std::min( theta1, -tau/2./_vN_til[ i ] );

    assert( theta1 > 0 );

    if ( theta1 == 1 )
        for ( int i = 0; i < nv; i++ )
            _v_til[ i ] = theta1*_vN_til[ i ];

    else
    {
        for ( int i = 0; i < nv; i++ )
            v_diff[ i ] = _vN_til[ i ] - _vCP_til[ i ];

        double a = __Dot( nv, v_diff, v_diff );
        double b = 2. * __Dot( nv, v_diff, _vCP_til );
        double c = __Dot( nv, _vCP_til, _vCP_til ) - ( zeta_N*_Delta )*( zeta_N*_Delta );
        double root = b*b - 4.*a*c;
        bool exist = false;

        if ( root >= 0 )
        {
            double r2 = ( ( -b + std::sqrt( root ) ) / ( 2.*a ) );

            if ( r2 > 1 )
            {
                exist = true;
                theta2 = 1;
            }

            else if ( r2 > 0 )
            {
                exist = true;
                theta2 = r2;
            }

            if ( exist )
                for ( int i = n; i < nv; i++ )
                    if ( v_diff[ i ] < 0 )
                        theta2 = std::min( theta2, ( -tau/2. - _vCP_til[ i ] ) / ( v_diff[ i ] ) );

            if ( theta2 < 0 )
                exist = false;
        }

        if ( !exist )
        {
            theta3 = std::min( theta3, zeta_N * _Delta / __Norm_2( nv, _vCP_til ) );

            for ( int i = n; i < nv; i++ )
                if ( _vCP_til[ i ] < 0 )
                    theta3 = std::min( theta3, -tau/2./_vCP_til[ i ] );

            assert( theta3 > 0 );

            for ( int i = 0; i < nv; i++ )
                _v_til[ i ] = theta3 * _vCP_til[ i ];
        }

        else
            for ( int i = 0; i < nv; i++ )
                _v_til[ i ] = ( 1-theta2 ) * _vCP_til[ i ] + theta2 * _vN_til[ i ];

        double theta1xvN_til[ _E_n+_E_g+2*_E_n ];

        for ( int i = 0; i < nv; i++ )
            theta1xvN_til[ i ] = theta1 * _vN_til[ i ];

        if ( this->model_m( _s, _v_til ) >= this->model_m( _s, theta1xvN_til ) )
            for ( int i = 0; i < nv; i++ )
                _v_til[ i ] = theta1xvN_til[ i ];
    }

    /*
      std::cerr << "\ntheta1 = " << theta1;
      std::cerr << "\ntheta2 = " << theta2;
      std::cerr << "\ntheta3 = " << theta3;
      std::cerr << "\nNorm(_v_til) = " << __Norm_2( nv, _v_til );
      std::cerr << "\n_v_til = ";
      __print_vec( nv, _v_til );
      std::cerr << "\ndx_n = ";
      exit(0);
    */

}

template<typename Data,template<class> class Problem>
double SolverConstrained<Data,Problem>::get_vpred( double *_rhs, double *_v_til )
{
    int n = M_prob.n(), nv = n + _E_g+2*n;
    double A_hatTxv_til[ _E_h+_E_g+2*_E_n ], rhspAv[ _E_h + _E_g+2*_E_n ];

    for ( int i = 0; i < _E_h+_E_g+2*n; i++ )
    {
        A_hatTxv_til[ i ] = 0;

        for ( int j = 0; j < nv; j++ )
            A_hatTxv_til[ i ] += A_hat[ j ][ i ] * _v_til[ j ];

        rhspAv[ i ] = _rhs[ i ] + A_hatTxv_til[ i ];
    }

    return __Norm_2( _E_h+_E_g+2*n, _rhs ) - __Norm_2( _E_h+_E_g+2*n, rhspAv );
}

template<typename Data,template<class> class Problem>
double SolverConstrained<Data,Problem>::get_normal_steps( double *_s, double _Delta,
        double *_dx_normal, double *_ds_normal,
        double *_v_til )
{
    int n = M_prob.n(), nv = n + _E_g+2*n;
    double vCP_til[ _E_n + _E_g+2*_E_n ], vN_til[ _E_n + _E_g+2*_E_n ],
           rhs[ _E_h + _E_g+2*_E_n ], A_hatxrhs[ _E_n + _E_g+2*_E_n ],
           A_hatTxA_hat[ _E_h + _E_g+2*_E_n ][ _E_h + _E_g+2*_E_n ],
           A_hatTxA_hat2[ _E_h + _E_g+2*_E_n ][ _E_h + _E_g+2*_E_n ],
           A_hat2xrhs [ _E_h+_E_g+2*_E_n ];
    double zeros[ _E_n + _E_g+2*_E_n ], dummy[ _E_n + _E_g+2*_E_n ], x_til[ _E_h + _E_g+2*_E_n ];

    for ( int m = 0; m < _E_h; m++ )
        rhs[ m ] = hx[ m ];

    for ( int m = 0; m < _E_g+2*n; m++ )
        rhs[ _E_h + m ] = gx[ m ] + _s[ m ];

    // A^T * A
    for ( int mi = 0; mi < _E_h+_E_g+2*n; mi++ )
    {
        for ( int mj = 0; mj < _E_h+_E_g+2*n; mj++ )
        {
            A_hatTxA_hat[ mi ][ mj ] = 0;

            for ( int i = 0; i < nv; i++ )
                A_hatTxA_hat[ mi ][ mj ] += A_hat[ i ][ mi ] * A_hat[ i ][ mj ];
        }
    }

    // (A^T * A )^2
    for ( int mi = 0; mi < _E_h+_E_g+2*n; mi++ )
    {
        for ( int mj = 0; mj < _E_h+_E_g+2*n; mj++ )
        {
            A_hatTxA_hat2[ mi ][ mj ] = 0;

            for ( int i = 0; i < _E_h+_E_g+2*n; i++ )
                A_hatTxA_hat2[ mi ][ mj ] += A_hatTxA_hat[ i ][ mi ] * A_hatTxA_hat[ i ][ mj ];
        }
    }

    for ( int i = 0; i < nv; i++ )
    {
        zeros[ i ] = 0;
        A_hatxrhs[ i ] = 0;

        for ( int m = 0; m < _E_h+_E_g+2*n; m++ )
            A_hatxrhs[ i ] += A_hat[ i ][ m ] * rhs[ m ];
    }

    for ( int i = 0; i < _E_h+_E_g+2*n; i++ )
    {
        A_hat2xrhs[ i ] = 0;

        for ( int m = 0; m < _E_h+_E_g+2*n; m++ )
            A_hat2xrhs[ i ] += A_hatTxA_hat2[ i ][ m ] * rhs[ m ];
    }

    double norm = __Norm_2( nv, A_hatxrhs );
    double alpha = norm*norm / __Dot( _E_h+_E_g+2*n, rhs, A_hat2xrhs );

    for ( int i = 0; i < nv; i++ )
        vCP_til[ i ] = -alpha * A_hatxrhs[ i ];

    this->solve_AtA_system( zeros, rhs, dummy, x_til );

    for ( int i = 0; i < n+_E_g+2*n; i++ )
    {
        vN_til[ i ] = 0;

        for ( int m = 0; m < _E_h+_E_g+2*n; m++ )
            vN_til[ i ] += -A_hat[ i ][ m ] * x_til[ m ];
    }

    this->dogleg_solve( vCP_til, vN_til, _s, _Delta, _v_til );

    for ( int i = 0; i < n; i++ )
        _dx_normal[ i ] = _v_til[ i ];

    for ( int i = 0; i < _E_g+2*n; i++ )
        _ds_normal[ i ] = _s[ i ] * _v_til[ n+i ];

    /*
      std::cerr << "\nalpha = " << alpha;
      std::cerr << "\nATA = \n";
      for( int i = 0; i < _E_h+_E_g+2*n; i++ )
      {
      printf("\n");
      for( int j = 0; j < _E_h+_E_g+2*n; j++ )
      fprintf( stderr, "%10.1f", A_hatTxA_hat2[i][j]);
      }
      std::cerr << "\n";
      std::cerr << "\nhx = [";
      __print_vec( _E_h, hx );
      std::cerr << "]\n";
      std::cerr << "\ngx = [";
      __print_vec( _E_g+2*n, gx );
      std::cerr << "]\n";
      std::cerr << "\ns = [";
      __print_vec( _E_g+2*n, _s );
      std::cerr << "]\n";
      std::cerr << "\nA_hat = [";
      this->print_A_hat();
      std::cerr << "]\n";
      std::cerr << "\nrhs = [";
      __print_vec( _E_h+_E_g+2*n, rhs );
      std::cerr << "]\n";
      std::cerr << "\nA_hatxrhs = [";
      __print_vec( nv, A_hatxrhs );
      std::cerr << "]\n";
      __print_vec( n, _dx_normal );
      std::cerr << "\nds_n = ";
      __print_vec( _E_g+2*n, _ds_normal );
      exit(0);
    */

    return this->get_vpred( rhs, _v_til );

}

template<typename Data,template<class> class Problem>
double SolverConstrained<Data,Problem>::get_q( double *_v_til, double *_w_til, double MU )
{
    int n = M_prob.n(), nd = n + _E_g+2*n;
    double d_til[ _E_n+_E_g+2*_E_n ], gr[ _E_n+_E_g+2*_E_n ], Gd_til[ _E_n+_E_g+2*_E_n ];

    __SumVecs( nd, _v_til, _w_til, d_til );

    for ( int i = 0; i < n; i++ )
        gr[ i ] = gradf[ i ];

    for ( int i = n; i < nd; i++ )
        gr[ i ] = -MU;

    for ( int i = 0; i < nd; i++ )
    {
        Gd_til[ i ] = 0;

        for ( int j = 0; j < nd; j++ )
            Gd_til[ i ] += G[ i ][ j ] * d_til[ j ];
    }

    return __Dot( nd, gr, d_til ) + 0.5 * __Dot( nd, d_til, Gd_til );
}

template<typename Data,template<class> class Problem>
void SolverConstrained<Data,Problem>::get_Pr( double *_r, double *_Pr )
{
    int n = M_prob.n(), nw = n + _E_g+2*n;
    double zeros[ _E_n + _E_g+2*_E_n ], dummy[ _E_n + _E_g+2*_E_n ],
           Atr[ _E_h + _E_g+2*_E_n ], x_Atr[ _E_h + _E_g+2*_E_n ];

    for ( int i = 0; i < nw; i++ )
    {
        zeros[ i ] = 0;
    }

    for ( int i = 0; i < _E_h+_E_g+2*n; i++ )
    {
        Atr[ i ] = 0;

        for ( int j = 0; j < nw; j++ )
            Atr[ i ] += A_hat[ j ][ i ] * _r[ j ];
    }

    this->solve_AtA_system( zeros, Atr, dummy, x_Atr );

    for ( int i = 0; i < nw; i++ )
    {
        _Pr[ i ] = _r[i];

        for ( int m = 0; m < _E_h+_E_g+2*n; m++ )
            _Pr[ i ] += -A_hat[ i ][ m ] * x_Atr[ m ];
    }
}

template<typename Data,template<class> class Problem>
bool SolverConstrained<Data,Problem>::slack_feasible( double *_w_til, double *_v_til )
{
    int n = M_prob.n();
    double tau = options.tau;
    bool feasible = true;

    for ( int i = n; i < n+_E_g+2*n; i++ )
    {
        if ( _w_til[ i ] + _v_til[ i ] + tau < 1e-13 )
            feasible = false;

        //std::cerr << "\n" << i << ": " << _w_til[ i ] + _v_til[ i ] + tau;
    }

    return feasible;
}

template<typename Data,template<class> class Problem>
double SolverConstrained<Data,Problem>::get_tangen_steps( double *_s, double _Delta, double *_v_til, double MU,
        double *_dx_tangen, double *_ds_tangen,
        bool &truss_exit, int &n_CGiter  )
{
    int n = M_prob.n();
    int nw = n + _E_g+2*n;
    double last_slack_feasible_rg;
    double w_til[ _E_n+_E_g+2*_E_n ], wP_til[ _E_n+_E_g+2*_E_n ],
           last_slack_feasible_w_til[ _E_n+_E_g+2*_E_n ],
           r[ _E_n+_E_g+2*_E_n ], rP[ _E_n+_E_g+2*_E_n ], Gp[ _E_n+_E_g+2*_E_n ],
           g[ _E_n+_E_g+2*_E_n ], gP[ _E_n+_E_g+2*_E_n ],
           p[ _E_n+_E_g+2*_E_n ], pP[ _E_n+_E_g+2*_E_n ],
           last_slack_feasible_p[ _E_n+_E_g+2*_E_n ];

    truss_exit = false;
    n_CGiter = 0;

    for ( int i = 0; i < nw; i++ )
    {
        Gp[ i ] = 0;

        for ( int j = 0; j < nw; j++ )
            Gp[ i ] += G[ i ][ j ] * _v_til[ j ];
    }

    for ( int i = 0; i < n; i++ )
        r[ i ] = gradf[ i ] + Gp[ i ];

    for ( int i = n; i < nw; i++ )
        r[ i ] = -MU + Gp[ i ];

    get_Pr( r, g );

    /*
      this->print_A_hat();
      __print_vec( nw, r );
      __print_vec( nw, g );
      this->print_G();
    */

    for ( int i = 0; i < nw; i++ )
    {
        w_til[ i ] = 0;
        wP_til[ i ] = 0;
        last_slack_feasible_w_til[ i ] = 0;
        p[ i ] = -g[ i ];
        pP[ i ] = -g[ i ];
        last_slack_feasible_p[ i ] = -g[ i ];
    }

    last_slack_feasible_rg = __Dot( nw, g, r );

    double tol = 0.01 * std::sqrt( __Dot( nw, g, r ) );
    int iter = 0;
    bool done = false;
    double gamma, alpha, beta;

    while ( iter < 2*( n-_E_h ) && !done )
    {
#ifdef DEBUG_OUT_CG
        std::cerr << "\nCGiter = " << iter << ",\t  r'g = " << __Dot( nw, g, r );
#endif

        for ( int i = 0; i < nw; i++ )
        {
            Gp[ i ] = 0;

            for ( int j = 0; j < nw; j++ )
                Gp[ i ] += G[ i ][ j ] * p[ j ];
        }

        gamma = __Dot( nw, p, Gp );

        if ( gamma <= 0 )
        {
            double a = __Dot( nw, p, p );
            double b = 2. * __Dot( nw, w_til, p );
            double c = __Dot( nw, w_til, w_til ) + __Dot( nw, _v_til, _v_til ) - _Delta*_Delta;
            double thetaa = ( -b + std::sqrt( b*b - 4.*a*c ) ) / ( 2.*a );
            assert( thetaa > 0 );

            for ( int i = 0; i < nw; i++ )
                wP_til[ i ] = w_til[ i ] + thetaa * p[ i ];

            done = true;
            truss_exit = true;
#ifdef DEBUG_OUT_CG
            std::cerr << "\nCG indefinitness EXIT:  ||wP_til||^2 + ||v_til||^2 = "
                      << __Dot( nw, wP_til, wP_til ) + __Dot( nw, _v_til, _v_til )
                      << " ,\t Delta^2 = "
                      << _Delta*_Delta << "\n";
#endif
        }

        if ( !done )
        {
            alpha = __Dot( nw, r, g ) / gamma;

            for ( int i = 0; i < nw; i++ )
                wP_til[ i ] = w_til[ i ] + alpha * p[ i ];

            if ( __Dot( nw, wP_til, wP_til ) + __Dot( nw, _v_til, _v_til ) > _Delta*_Delta )
            {
                double a = __Dot( nw, p, p );
                double b = 2. * __Dot( nw, w_til, p );
                double c = __Dot( nw, w_til, w_til ) + __Dot( nw, _v_til, _v_til ) - _Delta*_Delta;
                double thetab = ( -b + std::sqrt( b*b - 4.*a*c ) ) / ( 2.*a );
                assert( thetab > 0 );

                for ( int i = 0; i < nw; i++ )
                    wP_til[ i ] = w_til[ i ] + thetab * p[ i ];

                done = true;
                truss_exit = true;
#ifdef DEBUG_OUT_CG
                std::cerr << "\nCG truss EXIT:  ||wP_til||^2 + ||v_til||^2 = "
                          << __Dot( nw, wP_til, wP_til ) + __Dot( nw, _v_til, _v_til )
                          << " ,\t Delta^2 = "
                          << _Delta*_Delta << "\n";
#endif
            }
        }

        if ( !done )
        {
            for ( int i = 0; i < nw; i++ )
                rP[ i ] = r[ i ] + alpha * Gp[ i ];

            get_Pr( rP, gP );

            /*
              this->print_A_hat();
              __print_vec( nw, rP );
              __print_vec( nw, gP );
            */
            if ( __Dot( nw, gP, rP ) < tol )
            {
                done = true;
#ifdef DEBUG_OUT_CG
                std::cerr << "\nCG tol EXIT: \trP'gP = "
                          << __Dot( nw, gP, rP )
                          << " \t ( tol = "
                          << tol << " )\n";
#endif
            }
        }

        if ( !done )
        {
            beta = __Dot( nw, gP, rP ) / __Dot( nw, g, r );

            for ( int i = 0; i < nw; i++ )
            {
                pP[ i ] = -gP[ i ] + beta * p[ i ];
                w_til[ i ] = wP_til[ i ];
                r[ i ] = rP[ i ];
                g[ i ] = gP[ i ];
                p[ i ] = pP[ i ];
            }

            if ( this->slack_feasible( wP_til, _v_til ) )
            {
                for ( int i = 0; i < nw; i++ )
                {
                    last_slack_feasible_w_til[ i ] = wP_til[ i ];
                    last_slack_feasible_p[ i ] = p[ i ];
                }

                last_slack_feasible_rg = __Dot( nw, g, r );
            }
        }

        iter++;
    }

    if ( !( this->slack_feasible( wP_til, _v_til ) ) )
    {
        double a = __Dot( nw, last_slack_feasible_p, last_slack_feasible_p );
        double b = 2. * __Dot( nw, last_slack_feasible_w_til, last_slack_feasible_p );
        double c = __Dot( nw, last_slack_feasible_w_til, last_slack_feasible_w_til )
                   + __Dot( nw, _v_til, _v_til ) - _Delta*_Delta;
        double thetac = ( -b + std::sqrt( b*b - 4.*a*c ) ) / ( 2.*a );
        assert( thetac > 0 );
        double tau = options.tau;
        truss_exit = true;

        for ( int i = n; i < nw; i++ )
            if ( last_slack_feasible_p[ i ] < 0 )
            {
                thetac = std::min( thetac,
                                   ( -tau - _v_til[ i ] - last_slack_feasible_w_til[ i ] )
                                   / last_slack_feasible_p[ i ] );
            }

        assert( thetac > 0 );

        for ( int i = 0; i < nw; i++ )
            wP_til[ i ] = last_slack_feasible_w_til[ i ] + thetac * last_slack_feasible_p[ i ];

#ifdef DEBUG_OUT_CG
        std::cerr << "\nCorrected slack feasibility violation of wP_til: tau = " << tau << ", wP_til_s + v_til_s = " ;

        for ( int i = n; i < nw; i++ )
            std::cerr << "\n " << wP_til[ i ] + _v_til[ i ];

        std::cerr<< "\nr'g = " << last_slack_feasible_rg << "\n";
#endif
    }

    n_CGiter = iter;

    for ( int i = 0; i < n; i++ )
        _dx_tangen[ i ] = wP_til[ i ];

    for ( int i = 0; i < _E_g+2*n; i++ )
        _ds_tangen[ i ] = _s[ i ] * wP_til[ n+i ];

    /*
      __print_vec( nw, _v_til );
      std::cerr << "\nwPtil = ";
      __print_vec( nw, wP_til );

      this->print_G();
      __print_vec( nw, wP_til );
      for( int i = 0; i < nw; i++ )
      {
      Gp[ i ] = 0;
      for( int j = 0; j < nw; j++ )
      Gp[ i ] += G[ i ][ j ] * wP_til[ j ];
      }
      __print_vec( nw, Gp );

      __print_vec( nw, wP_til );
      __print_vec( nw, _v_til );
      std::cerr << "\ntol = "<< tol;
      this->print_G();
      __print_vec( nw, Gp );
      __print_vec( n, gradf );
      this->print_A_hat();
      __print_vec( nw, r );
      __print_vec( nw, g );
      //exit(0);
      */


    return this->get_q( _v_til, wP_til, MU );

}

template<typename Data,template<class> class Problem>
double SolverConstrained<Data,Problem>::get_steps( double *_s, double _Delta, double MU,
        double *_dx, double *_ds, double &_vpred, double &_q,
        bool &truss_exit, int &n_CGiter )
{
    int n = M_prob.n();
    double dx_normal[ _E_n ], ds_normal[ _E_g+2*_E_n ], v_til[ _E_n+_E_g+2*_E_n ];
    double dx_tangen[ _E_n ], ds_tangen[ _E_g+2*_E_n ], d[ _E_n+_E_g+2*_E_n ];

    _vpred = this->get_normal_steps( _s, _Delta, dx_normal, ds_normal, v_til );
    _q = this->get_tangen_steps( _s, _Delta, v_til, MU, dx_tangen, ds_tangen, truss_exit, n_CGiter );

    __SumVecs( n, dx_normal, dx_tangen, _dx );
    __SumVecs( _E_g+2*n, ds_normal, ds_tangen, _ds );

    /*
      std::cerr << "\n dx_normal = ";
      __print_vec( n, dx_normal );
      std::cerr << "\n ds_normal = ";
      __print_vec( _E_g+2*n, ds_normal );
      std::cerr << "\n dx_tangen = ";
      __print_vec( n, dx_tangen );
      std::cerr << "\n ds_tangen = ";
      __print_vec( _E_g+2*n, ds_tangen );
      exit(0);
    */

    for ( int i = 0; i < n; i++ )
        d[ i ] = _dx[ i ];

    for ( int i = 0; i < _E_g+2*n; i++ )
        d[ n+i ] = _ds[ i ];

    return __Norm_2( n+_E_g+2*n, d );
}

template<typename Data,template<class> class Problem>
double SolverConstrained<Data,Problem>::get_phi( double *_x, double *_s, double MU,
        double nu, double *rhs )
{
    double phi_val;

#if 0

    int n = M_prob.n();
    double f_at_x, g_at_x[ _E_g+2*_E_n ], h_at_x[ _E_h ];


    M_prob.fn_double( _x, f_at_x );
    M_prob.gn_double( _x, g_at_x );

    for ( int i = 0; i < n; i++ )
    {
        g_at_x[ _E_g+i ]   = M_prob.x_l( i ) - _x[i];
        g_at_x[ _E_g+n+i ] = _x[i] - M_prob.x_u( i );
    }

    M_prob.hn_double( _x, h_at_x );

    for ( int m = 0; m < _E_h; m++ )
        rhs[ m ] = h_at_x[ m ];

    for ( int m = 0; m < _E_g+2*n; m++ )
        rhs[ _E_h + m ] = g_at_x[ m ] + _s[ m ];

    phi_val = f_at_x;

    for ( int i = 0; i < _E_g+2*n; i++ )
        phi_val += -MU * std::log( _s[ i ] );

    phi_val += nu * __Norm_2( _E_g+2*n, rhs );
#else
    f_type __fx;
    g_type __gx;
    h_type __hx;

    // get the order 0 information only
    M_prob.evaluate ( _x, __fx, diff_order<0>() );
    M_prob.evaluate ( _x, __gx, diff_order<0>() );
    M_prob.evaluate ( _x, __hx, diff_order<0>() );

    double __f_at_x = __fx.value( 0 );
    double __g_at_x[ _E_g+2*_E_n ];

    for ( uint __i = 0; __i < __gx.value().size(); ++__i )
        __g_at_x[__i] = __gx.value( __i );

    double __h_at_x[ _E_h ];

    for ( uint __i = 0; __i < _E_h; ++__i )
        __h_at_x[__i] = __hx.value( __i );

    for ( int m = 0; m < _E_h; m++ )
        rhs[ m ] = __h_at_x[ m ];

    for ( int m = 0; m < _E_g+2*_E_n; m++ )
        rhs[ _E_h + m ] = __g_at_x[ m ] + _s[ m ];

    phi_val = __f_at_x;

    for ( int i = 0; i < _E_g+2*_E_n; i++ )
        phi_val += -MU * std::log( _s[ i ] );

    phi_val += nu * __Norm_2( _E_g+2*_E_n, rhs );
#endif
    return phi_val;
}

template<typename Data,template<class> class Problem>
double SolverConstrained<Data,Problem>::get_phi( double *_x, double *_s, double MU, double nu )
{
    double rhs[ _E_h + _E_g+2*_E_n ];
    double phi_val = this->get_phi( _x, _s, MU, nu, rhs );

    return phi_val;
}

template<typename Data,template<class> class Problem>
double SolverConstrained<Data,Problem>::get_SOC( double *_s, double *_rhs_xs_pd,
        double *_xy, double *_sy )
{
    int n = M_prob.n();

    // Second order corrections y:
    double zeros[ _E_n + _E_g+2*_E_n ], dummy[ _E_n + _E_g+2*_E_n ],
           x_til[ _E_h + _E_g+2*_E_n ], y[ _E_n + _E_g+2*_E_n ];

    for ( int i = 0; i < n+_E_g+2*_E_n; i++ )
        zeros[ i ] = 0;

    this->solve_AtA_system( zeros, _rhs_xs_pd, dummy, x_til );

    for ( int i = 0; i < n+_E_g+2*n; i++ )
    {
        y[ i ] = 0;

        for ( int m = 0; m < _E_h+_E_g+2*n; m++ )
            y[ i ] += A_hat[ i ][ m ] * x_til[ m ];
    }

    for ( int i = 0; i < n; i++ )
        _xy[ i ] = y[ i ];

    for ( int i = 0; i < _E_g+2*n; i++ )
        _sy[ i ] = y[ n+i ];

    return __Norm_2( n+_E_g+2*n, y );
}

template<typename Data,template<class> class Problem>
double SolverConstrained<Data,Problem>::get_phi( double *_s, double MU, double nu )
{
    int n = M_prob.n();
    double rhs[ _E_h + _E_g+2*_E_n ];
    double phi_val;

    for ( int m = 0; m < _E_h; m++ )
        rhs[ m ] = hx[ m ];

    for ( int m = 0; m < _E_g+2*n; m++ )
        rhs[ _E_h + m ] = gx[ m ] + _s[ m ];

    phi_val = fx;

    for ( int i = 0; i < _E_g+2*n; i++ )
        phi_val += -MU * std::log( _s[ i ] );

    phi_val += nu * __Norm_2( _E_g+2*n, rhs );

    /*
      std::cerr << "\nfx = " << fx << "\n";
      std::cerr << "\ns = \n";
      __print_vec( _E_g+2*n, _s );
      std::cerr << "\ngx = \n";
      __print_vec( _E_g+2*n, gx );
      std::cerr << "\nrhs( x, s ) = \n";
      __print_vec( _E_h+_E_g+2*n, rhs );
      std::cerr << "\nhn = \n";
      __print_vec( _E_h, hx );
    */

    return phi_val;
}


template<typename Data,template<class> class Problem>
bool
SolverConstrained<Data,Problem>::optimize( double __x[_E_n] )
{
    int n = M_prob.n();

    M_solver_stats.clear();


    std::cout << "\n[SolverConstrained<Data,Problem>::optimize()]: Optimizing...\n";

    // Controlling parameters
    double Delta = std::min( options.Delta_init, ( double )Max_Allowed_Delta );
    double eMU = options.eMU_init;
    double MU = options.MU_init;
    double MU_last = MU;
    double nu = options.nu0_init;
    double rho_N = options.rho_N;

    double rhs_xs_p_d[ _E_h + _E_g+2*_E_n ];

    double dx[ _E_n ];
    double x_p_dx[ _E_n ];
    double s[ _E_g+2*_E_n ];
    double ds[ _E_g+2*_E_n ];
    double s_p_ds[ _E_g+2*_E_n ];
    double lambda_h[ _E_h ], lambda_g[ _E_g+2*_E_n ];
    double vpred, q, pred, ared, gamma, norm_d;

    M_prob.copy_x0_to_x( __x );
    this->initialize_solver_with_x( __x, s, lambda_h, lambda_g, MU );

    double E = this->evaluate_E_current( s, lambda_h, lambda_g, MU );

    int err_iter = 0, hom_iter = 0;

    while ( err_iter < 500 && hom_iter < options.max_hom_iter && E > options.eTOL )
    {

        double max_TR_iter = options.max_TR_iter;

        if ( hom_iter+1 >= options.max_hom_iter )
            max_TR_iter = 2 * max_TR_iter;

        err_iter++;

#ifdef DEBUG_OUT_HM
        std::cerr << "\n================================================================================"
                  << "\n         homotopy iter = " << hom_iter
                  << ",  MU = " << MU
                  << ",  E(x,y,MU) = " << this->evaluate_E_current( s, lambda_h, lambda_g, MU )
                  << "\n================================================================================\n";
#endif

        int iter = 0;

        while ( iter < max_TR_iter && E > eMU )
        {

            E = this->evaluate_E_current( s, lambda_h, lambda_g, MU );

            double Delta_used = Delta;
            bool truss_exit = false, SOC = false;
            int n_CGiter = 0;

            bool evaluate_Hess = false;

            if ( iter%2 == 0 )
                evaluate_Hess = true;

            norm_d = this->get_steps( s, Delta, MU, dx, ds, vpred, q, truss_exit, n_CGiter );

            if ( q > 0 )
                nu = std::max( 1.5*nu, q / ( ( 1 - rho_N ) * vpred ) );

            pred = -q + nu * vpred;
            __SumVecs( n, __x, dx, x_p_dx );
            __SumVecs( _E_g+2*n, s, ds, s_p_ds );
            ared = get_phi( s, MU, nu ) - get_phi( x_p_dx, s_p_ds, MU, nu, rhs_xs_p_d );
            gamma = ared / pred;

#ifdef DEBUG_OUT_TR
            std::cerr << "\n-------------------------\n"
                      << "\nTRiter = " << iter
                      << "\n  E(x,s,MU) = " << E
                      << "\n  gamma     = " << gamma

                      << "\n  vpred     = " << vpred
                      << "\n  pred      = " << pred
                      << "\n  nu        = " << nu
                      << "\n  ared      = " << ared
                      << "\n  norm_d    = " << norm_d
                      << "\n  Delta     = " << Delta

                      ;
            __print_vec( n, dx );
            __print_vec( _E_g+2*n, ds );
#endif

            if ( gamma >= options.eta )
            {
                __CopyVecs( n, x_p_dx, __x );
                __CopyVecs( _E_g+2*n, s_p_ds, s );

                if ( gamma >= options.rho_big )
                    Delta = std::max( options.rho_increase_big * norm_d , Delta );

                else if ( gamma >= options.rho_small )
                    Delta = std::max( options.rho_increase_small * norm_d , Delta );

                Delta = std::min( Delta, ( double )Max_Allowed_Delta );
            }

            else
            {
                SOC = false;
                double xy[ _E_n ], sy[ _E_g+2*_E_n ];
                double norm_xypsy = this->get_SOC( s_p_ds, rhs_xs_p_d, xy, sy );
#ifdef DEBUG_OUT_TR
                std::cerr << "\n  SECOND ORDER CORRECTION:"
                          << "\n   ||xy||+||sy||   = " << norm_xypsy;
#endif

                if ( norm_xypsy > 0 )
                {
                    bool sy_feasible = true;
                    double tau = options.tau;

                    for ( int i = 0; i < _E_g+2*n; i++ )
                        if ( ds[ i ] + sy[ i ] + tau*s[ i ] < 0 )
                            sy_feasible = false;

                    if ( sy_feasible )
                    {
                        __SumVecs( n, __x, dx, x_p_dx );
                        __SumVecs( _E_g+2*n, s, ds, s_p_ds );

                        __SumVecs( n, x_p_dx, xy, x_p_dx );
                        __SumVecs( _E_g+2*n, s_p_ds, sy, s_p_ds );

                        ared = get_phi( s, MU, nu ) - get_phi( x_p_dx, s_p_ds, MU, nu );
                        gamma = ared / pred;
#ifdef DEBUG_OUT_TR
                        std::cerr << "\n   gamma     = " << gamma;
#endif

                        if ( gamma >= options.eta )
                        {
                            SOC = true;
                            __CopyVecs( n, x_p_dx, __x );
                            __CopyVecs( _E_g+2*n, s_p_ds, s );
                            // (Delta = Delta)
                        }
                    }
                }

                if ( !SOC )
                    Delta = options.rho_decrease * Delta;
            }

            if ( gamma >= options.eta ) // (successful step)
            {
                if ( !options.Primal_Dual )
                    MU_last = MU;

                this->update_opt_objs( __x, s, lambda_h, lambda_g, MU_last, evaluate_Hess );

                M_solver_stats.push_all( MU, E, this->evaluate_E_current( s, lambda_h, lambda_g, 0 ),
                                          n_CGiter, truss_exit, SOC,
                                          Delta_used, ared, pred, gamma );
                iter++;
                MU_last = MU;
            }

#ifdef DEBUG_OUT_TR
            std::cerr << "\n  E(x+Dx,s+ds,MU) = " << E << "  ( MU = " << MU << ")\n";
#endif
        }

        MU = MU * options.theta_N;

        // be careful with reducing eMU because we reach Floating point exception
        // so cut out to 0 when reaches 1e-10
        if ( eMU < 1e-10 )
            eMU = 0;

        else
            eMU = eMU * options.theta_N * eMU;

        // reinitialize Delta
        Delta = options.Delta_init;

        // reinitialize Delta
        nu = options.nu0_init;

        // reevaluate E to start over the homotopy loop
        E = this->evaluate_E_current( s, lambda_h, lambda_g, 0 );

#ifdef DEBUG_OUT_HM
        std::cerr << "\n  E(x+Dx,s+ds,0) = " << E << "\n";
#endif
        // go to next homotopy loop
        hom_iter++;
    }

    M_solver_stats.show();
    //std::cerr << "Number of homotopy levels: " << hom_iter;

    /*
      std::cerr << "\n\nx = \n";
      __print_vec( n, __x );
      std::cerr << "\n\ns = \n";
      __print_vec( _E_g+2*n, s );
      std::cerr << "\n";
    */

    //M_prob.print_stationaryN_x( __x, s, lambda_h, lambda_g );

    //std::cerr << "\n";

    return true;
}

//
// Stats
//
template<typename Data, template<class> class Problem>
void
SolverConstrained<Data,Problem>::Stats::clear()
{
    x_hstr.clear();
    mu_hstr.clear();
    E_hstr.clear();
    E0_hstr.clear();

    n_CGiter_hstr.clear();
    truss_exit_hstr.clear();
    SOC_hstr.clear();

    Delta_hstr.clear();
    ared_hstr.clear();
    pred_hstr.clear();
    gamma_hstr.clear();
}

template<typename Data, template<class> class Problem>
void
SolverConstrained<Data,Problem>::Stats::push_x( const double *__x )
{
    for ( int i = 0; i < solver_type::_E_n; i++ )
        x_hstr.push_back( __x[i] );
}

template<typename Data, template<class> class Problem>
void
SolverConstrained<Data,Problem>::Stats::push_all( double mu, double E, double E0,
        int n_CGiter,
        bool truss_exit, bool SOC,
        double Delta, double ared, double pred, double gamma )
{
    //push_x_in_hstr( _x );
    mu_hstr.push_back( mu );
    E_hstr.push_back( E );
    E0_hstr.push_back( E0 );
    n_CGiter_hstr.push_back( n_CGiter );
    truss_exit_hstr.push_back( truss_exit );
    SOC_hstr.push_back( SOC );
    Delta_hstr.push_back( Delta );
    ared_hstr.push_back( ared );
    pred_hstr.push_back( pred );
    gamma_hstr.push_back( gamma );
}


template<typename Data, template<class> class Problem>
void
SolverConstrained<Data,Problem>::Stats::show( bool verbose ) const
{
    iter = mu_hstr.size();
    std::cerr << "\nScaled Trust-Region Method Statistics:\n\n";

    std::cerr << "  k     mu    E(x,s;0)   E(x,s;mu)  CGiters trustEx SOC  Delta^k   ared^k    pred^k    gamma^k \n";
    std::cerr << "------------------------------------------------------------------------------------------------\n";

    double mu_last = 0;

    for ( int k = 0; k < iter; k++ )
    {
        if ( mu_hstr[ k ] != mu_last )
        {
            fprintf( stderr, "%3d ", k+1 );
            fprintf ( stderr, " %3.1e", mu_hstr[ k ] );
            fprintf ( stderr, " %10.3e", E0_hstr[ k ] );
            fprintf ( stderr, " %10.3e", E_hstr[ k ] );
            fprintf ( stderr, " %5d", n_CGiter_hstr[ k ] );

            if ( truss_exit_hstr[ k ] == true )
                fprintf ( stderr, "      -Y-   " );

            else
                fprintf ( stderr, "       N    " );

            if ( SOC_hstr[ k ] == true )
                fprintf ( stderr, "-Y- " );

            else
                fprintf ( stderr, " N  " );

            fprintf ( stderr, " %3.1e", Delta_hstr[ k ] );
            fprintf ( stderr, " %9.2e", ared_hstr[ k ] );
            fprintf ( stderr, " %9.2e", pred_hstr[ k ] );
            fprintf ( stderr, " %7.1f", gamma_hstr[ k ] );
            std::cerr << std::endl;
            mu_last = mu_hstr[ k ];
        }

        else if ( !verbose )
        {
            fprintf( stderr, "%3d ", k+1 );
            fprintf ( stderr, "        " );
            fprintf ( stderr, " %10.3e", E0_hstr[ k ] );
            fprintf ( stderr, " %10.3e", E_hstr[ k ] );
            fprintf ( stderr, " %5d", n_CGiter_hstr[ k ] );

            if ( truss_exit_hstr[ k ] == true )
                fprintf ( stderr, "      -Y-   " );

            else
                fprintf ( stderr, "       N    " );

            if ( SOC_hstr[ k ] == true )
                fprintf ( stderr, "-Y- " );

            else
                fprintf ( stderr, " N  " );

            fprintf ( stderr, " %3.1e", Delta_hstr[ k ] );
            fprintf ( stderr, " %9.2e", ared_hstr[ k ] );
            fprintf ( stderr, " %9.2e", pred_hstr[ k ] );
            fprintf ( stderr, " %7.1f", gamma_hstr[ k ] );
            std::cerr << std::endl;
        }
    }

    std::cerr << "\nNumber of accepted steps:  " << iter << "\n";
}
} // Feel
#endif
