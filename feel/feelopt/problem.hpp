/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-14

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file problem.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
#ifndef PROBLEM_HPP
#define PROBLEM_HPP

namespace Feel
{
namespace ublas = boost::numeric::ublas;

template<typename T, int>
struct noop
{
    typedef T type;
};

enum
{
    FUNCTIONAL = 0, //!< functional identifier
    EQUALITIES,     //!< equalities identifier
    INEQUALITIES    //!< inequalities identifier
};

template<int O>
struct diff_order
{
    enum { _E_order = O };
};

/*!
  \class core_data

  \author Christophe Prud'homme
*/
template<typename Data,
         int Id = FUNCTIONAL >
class core_data
{
public:
    enum
    {
        _E_n  = Data::_E_n, //!< number of control variables
        _E_f  = 1, //!< number of objective functions
        _E_g  = Data::_E_g, //!< number of inequalities
        _E_ng = Data::_E_g + 2*_E_n, //!< size of informations stored for inequalities
        _E_h  = Data::_E_h, //!< number of equalities
        _E_nA = Data::_E_n + ( _E_n )*( _E_n-1 )/2 //!< number of entries in the matrix
    };
    core_data()
        :
        M_val (),
        M_grad(),
        M_hess()
    {
        for ( uint __i = 0; __i < M_grad.size(); ++__i )
        {
            M_grad[ __i ].resize( ( size_t )_E_n );
            M_hess[ __i ].resize( ( size_t )_E_n, ( size_t )_E_n );
        }
    }
    typedef typename Data::value_type numerical_type;
    typedef ublas::vector<double> numerical_gradient_type;
    //! matrix type
    typedef ublas::symmetric_matrix<double, ublas::upper> numerical_hessian_type;


    //! value type for the functional
    typedef typename IF<Id == FUNCTIONAL,
            typename array_fixed<numerical_type,_E_f>::type,
            typename IF<Id == EQUALITIES,
            typename array_fixed<numerical_type,_E_h>::type,
            typename array_fixed<numerical_type,_E_ng>::type
            >::Result
            >::Result value_type;

    //! value type for the functional gradient
    typedef typename IF<Id == FUNCTIONAL,
            typename array_fixed<numerical_gradient_type,_E_f>::type,
            typename IF<Id == EQUALITIES,
            typename array_fixed<numerical_gradient_type,_E_h>::type,
            typename array_fixed<numerical_gradient_type,_E_g>::type
            >::Result
            >::Result gradient_type;

    //! value type for the functional hessian: in packed format
    typedef typename IF<Id == FUNCTIONAL,
            typename array_fixed<numerical_hessian_type,_E_f>::type,
            typename IF<Id == EQUALITIES,
            typename array_fixed<numerical_hessian_type,_E_h>::type,
            typename array_fixed<numerical_hessian_type,_E_g>::type
            >::Result
            >::Result hessian_type;


    value_type   value() const
    {
        return M_val;
    }
    value_type & value()
    {
        return M_val;
    }

    numerical_type   value( int __i ) const
    {
        return M_val[ __i ];
    }
    numerical_type & value( int __i )
    {
        return M_val[ __i ];
    }

    gradient_type const& gradient() const
    {
        return M_grad;
    }
    gradient_type &      gradient()
    {
        return M_grad;
    }

    numerical_gradient_type & gradient( int __i )
    {
        return M_grad[ __i ];
    }
    numerical_gradient_type   gradient( int __i ) const
    {
        return M_grad[ __i ];
    }

    numerical_type & gradient( int __i, int __j )
    {
        return M_grad[ __i]( __j );
    }
    numerical_type   gradient( int __i, int __j ) const
    {
        return M_grad[ __i]( __j );
    }

    hessian_type const& hessian() const
    {
        return M_hess;
    }
    hessian_type &      hessian()
    {
        return M_hess;
    }

    numerical_hessian_type const& hessian( int __i ) const
    {
        return M_hess[ __i ];
    }
    numerical_hessian_type &      hessian( int __i )
    {
        return M_hess[ __i ];
    }

    numerical_type & hessian( int __i, int __j, int __k )
    {
        return M_hess[__i]( __j, __k );
    }
    numerical_type   hessian( int __i, int __j, int __k ) const
    {
        return M_hess[__i]( __j, __k );
    }

protected:
    value_type M_val;
    gradient_type M_grad;
    hessian_type M_hess;
};
/**
   \class functional
   \brief defines the functional type

   \author Christophe Prud'homme
*/
template<typename Data>
class functional
    :
public core_data<Data, FUNCTIONAL>
{
public:

    typedef core_data<Data, FUNCTIONAL> super;

    template<int Diff = 0>
    struct diff
    {
        enum { diff_order = Diff };
        typedef diff<diff_order> type;
    };
};

/**
   \class inequalities type
   \brief define the inequalities type

   \author Christophe Prud'homme
*/
template<typename Data>
class inequalities
    :
public core_data<Data, INEQUALITIES>
{
public:

    template<int Diff = 0>
    struct diff
    {
        enum { diff_order = Diff };

        typedef diff<diff_order> type;
    };
};



/**
   \class equalities type
   \brief define the equalities type

   \author Christophe Prud'homme
*/
template<typename Data>
class equalities
    :
public core_data<Data, EQUALITIES>
{
    template<int Diff = 0>
    struct diff
    {
        enum { diff_order = Diff };
        typedef diff<diff_order> type;
    };
};

/*!
  \class dummy_data
  \brief dummy data type


  \author Christophe Prud'homme
*/
template<int N>
struct dummy_data
{
    enum { _E_dummy = N };
};

/**
   \class problem
   \brief Optimization problem specifications

   \author Ivan Oliveira and Christophe Prud'homme
*/
template<typename Data>
class problem
    :
public Data
{
public:


    enum
    {
        _E_n   = Data::_E_n, //!< number of control variables
        _E_f   = 1,          //!< number of objective functions
        _E_g   = Data::_E_g, //!< number of inequality constraints
        _E_h   = Data::_E_h, //!< number of equality constraints
        _E_nA  = _E_n + ( _E_n )*( _E_n-1 )/2, //!< size of the matrix
        _E_nL  = 2*_E_n,                      //!< number of multiplers
        _E_nAL = _E_nL + ( _E_nL )*( _E_nL-1 )/2 //!< size of the multipliers matrix
    };

    /*!
      \typedef Data super
      \brief super class for the problem.

      This type defines the cost functional and eventually
      the constraints(inequalities and equalities)
    */
    typedef Data super;
    typedef Data data_type;
    typedef problem<Data> problem_type;

    //! automatic differentiation type of order one
    typedef typename super::ad_0_type ad_0_type;

    //! automatic differentiation type of order one
    typedef typename super::ad_1_type ad_1_type;

    //! automatic differentiation type of order two
    typedef typename super::ad_2_type ad_2_type;

    //! numerical type
    typedef typename super::ad_0_type::value_type value_type;

    //! vector type
    typedef ublas::vector<double> vector_type;

    //! matrix type
    typedef ublas::symmetric_matrix<double, ublas::upper> symmetric_matrix_type;
    typedef ublas::matrix<double> matrix_type;


    typedef functional<Data> f_type;
    typedef typename IF<_E_g!=0,
            inequalities<Data>,
            dummy_data<0> >::Result g_type;
    typedef typename IF<_E_h!=0,
            equalities<Data>,
            dummy_data<1> >::Result h_type;

    //! default constructor: initialize with infinite bounds
    problem();

    //! define specific bounds for each control variable
    problem( value_type x_definitions[ _E_n ][ 3 ] );

    problem( problem<Data> const& p )
        :
        super( p )
    {}

    problem( Data const& p )
        :
        super( p ),
        M_xlower ( _E_n ),
        M_xupper ( _E_n ),
        M_x0( _E_n )
    {
        this->define_problem( p.defs() );
        this->initialize_x0();
    }

    //! destructor
    ~problem();

    //! define the bounds of the problem
    //void define_problem( value_type x_definitions[ _E_n ][3] );
    template<typename Defs>
    void define_problem( Defs );

    //! \return number of control variables which are on
    int n() const
    {
        return _E_n;
    }

    //! \return size of the matrix
    int nA() const
    {
        return _E_nA;
    }

    //! \return lower bound of control variable \c i
    value_type x_l( unsigned int i ) const
    {
        return M_xlower(  i  );
    }

    //! \return upper bound  of control variable \c i
    value_type x_u( unsigned int i ) const
    {
        return M_xupper(  i  );
    }

    vector_type const& lowerBounds() const
    {
        return M_xlower;
    }
    vector_type const& upperBounds() const
    {
        return M_xupper;
    }

    vector_type distanceToLowerBounds( vector_type const& __x ) const
    {
        return __x - M_xlower;
    }
    vector_type distanceToUpperBounds( vector_type const& __x ) const
    {
        return M_xupper - __x;
    }

    /**
       test if __x_definitions is on
       \return \c true if greater in absolute value to \c _S_on_flag, \c false otherwise
    */
    bool isOn ( int __i ) const
    {
        return  __x_definitions[__i][0] >= value_type( .5 );
    }

    /*!
      \class value
      \brief compute the value of the functionals, equalities and inequalities associated with the problem

      \author Christophe Prud'homme

      \c Order is the maximum automatic differentiation order to perform while evaluating the
      functional or the constraints(inequalities and equalities).
    */
    template<int Order>
    struct value
    {
        enum
        {
            _E_n = problem_type::_E_n
        };

        // defines the automatic differentiation depending on \c Order
        typedef typename mpl::if_<( Order == 0 ),
                ad_0_type,
                typename mpl::if_<Order == 1,
                ad_1_type,
                ad_2_type>::type >::type ad_type;

        typedef typename mpl::if_<( _E_g!=0 ),
                typename data_type::template inequalities<ad_type>::type,
                         dummy_data<0> >::Result inequalities_array_type;

        typedef typename mpl::if_<_E_h!=0,
                typename data_type::template equalities<ad_type>::type,
                         dummy_data<0> >::Result equalities_array_type;

        //
        // F : functional
        // \todo : not yet completely fixed for multiple objective functions
        //
        value( problem_type* __pt, vector_type const& __x,  f_type& __fx )
        {
            ad_type __adt;
            // call data::f() : it is supposed to be static
            __pt->f( __x, __adt );
            // FIXME
            to ( __adt, __fx, mpl::int_<Order>() );
        }
        void to ( ad_type const& __adt, f_type& __fx, mpl::int_<0> )
        {
            __fx.value( 0 ) = __adt.value();
        }
        void to ( ad_type const& __adt, f_type& __fx, mpl::int_<1> )
        {
            __fx.value( 0 ) = __adt.value();
            __fx.gradient( 0 ) = __adt.grad();
        }
        void to ( ad_type const& __adt, f_type& __fx, mpl::int_<2> )
        {
            __fx.value( 0 ) = __adt.value();
            __fx.gradient( 0 ) = __adt.grad();
            __fx.hessian( 0 ) = __adt.hessian();
        }
#if 1
        //
        // G : inequalities
        //
        value( problem_type* __pt, vector_type const& __x,  g_type& __gx )
        {
            inequalities_array_type __adt;
            // call data::f() : it is supporsed to be static
            __pt->g( __x, __adt );
            to ( __adt, __gx, mpl::int_<Order>() );

            // in all cases add distances wrt the bounds
            for ( int __i = 0; __i < _E_n; ++__i )
            {
                __gx.value( _E_g + __i )   = __pt->x_l( __i ) - __x(  __i  );
                __gx.value( _E_g + _E_n + __i ) = __x(  __i  ) - __pt->x_u( __i );
            }
        }
        void to ( inequalities_array_type const& __adt, g_type& __gx, mpl::int_<0> )
        {
            for ( int __i = 0; __i < _E_g; ++__i )
            {
                __gx.value( __i ) = __adt[ __i ].value();


            }
        }
        void to ( inequalities_array_type const& __adt, g_type& __gx, mpl::int_<1> )
        {
            for ( int __i = 0; __i < _E_g; ++__i )
            {
                __gx.value( __i ) = __adt[ __i ].value();

                for ( int __j = 0; __j < _E_n; ++__j )
                {

                    __gx.gradient( __i,__j ) = __adt[ __i ].grad( __j );

                }
            }
        }
        void to ( inequalities_array_type const& __adt, g_type& __gx, mpl::int_<2> )
        {
            for ( int __i = 0; __i < _E_g; ++__i )
            {
                __gx.value( __i ) = __adt[ __i ].value();

                for ( int __j = 0; __j < _E_n; ++__j )
                {
                    __gx.gradient( __i, __j ) = __adt[ __i ].grad( __j );

                    for ( int __k = 0; __k <= __j; ++__k )
                    {
                        __gx.hessian( __i, __k, __j ) =__adt[ __i ].hessian( __k, __j );
                    }
                }
            }
        }
        //
        // H: equalities
        //
        value( problem_type* __pt, vector_type const& __x,  h_type& __hx )
        {
            equalities_array_type __adt;
            // call data::f() : it is supposed to be static
            __pt->h( __x, __adt );
            to ( __adt, __hx, mpl::int_<Order>() );
        }
        void to ( equalities_array_type const& __adt, h_type& __hx, mpl::int_<0> )
        {
            for ( int __i = 0; __i < _E_h; ++__i )
            {
                __hx.value( __i ) = __adt[ __i ].value();
            }
        }
        void to ( equalities_array_type const& __adt, h_type& __hx, mpl::int_<1> )
        {
            for ( int __i = 0; __i < _E_h; ++__i )
            {
                __hx.value( __i ) = __adt[ __i ].value();

                for ( int __j = 0; __j < _E_n; ++__j )
                {
                    __hx.gradient( __i,__j ) = __adt[ __i ].grad( __j );
                }
            }
        }
        void to ( equalities_array_type const& __adt, h_type& __hx, mpl::int_<2> )
        {
            for ( int __i = 0; __i < _E_h; ++__i )
            {
                __hx.value( __i ) = __adt[ __i ].value();

                for ( int __j = 0; __j < _E_n; ++__j )
                {
                    __hx.gradient( __i, __j ) = __adt[ __i ].grad( __j );

                    for ( int __k = 0; __k <= __j; ++__k )
                    {
                        __hx.hessian( __i, __k, __j ) =__adt[ __i ].hessian( __k, __j );
                    }
                }
            }
        }
#endif
    };

    /*!
      \brief evaluate components of the problem
      \c Type can be
      #) the objective functions
      #) the inequality contraints
      #) the equality contraints
    */
    template<int Order, typename Type>
    void evaluate( vector_type const& __x, Type& __fx, diff_order<Order> ) const
    {
        value<Order> ( const_cast<problem*>( this ), __x, __fx );
    }

    //! copy initialization x0 into x
    void copy_x0_to_x( vector_type& _x );

    //! copy x to initialization x0
    void copy_x_to_x0( vector_type const& _x );

    //! initialize the starting point \c x0
    void initialize_x0();

    //! use custom starting point
    void custom_initialize_x0( vector_type& __custom );

    //! prints current value of x
    void print( vector_type const& _x ) const;
    void print_complete( std::string const&, vector_type const& _x ) const;
    void print_complete( vector_type const& _x, vector_type const& _s ) const;
    void print_bound_constraints() const;
    void print_stationary_x( vector_type const& _x,
                             vector_type const& _lambda_l,
                             vector_type const& _lambda_u ) const;
    void print_stationaryN_x( vector_type const& _x,
                              vector_type const& _s,
                              vector_type const& _lambda_h,
                              vector_type const& _lambda_g ) const;

private:

    vector_type M_xlower;
    vector_type M_xupper;
    value_type __x_definitions[_E_n][3];

    vector_type M_x0;

    size_t M_n;
    size_t M_nA;

};

template<typename Data>
problem<Data>::problem()
    :
    Data(),
    M_xlower ( _E_n ),
    M_xupper ( _E_n ),
    M_x0( _E_n )
{
    M_xlower = ublas::scalar_vector<value_type>( M_xlower.size(), -inf );
    M_xupper = ublas::scalar_vector<value_type>( M_xupper.size(), +inf );
}


template<typename Data>
problem<Data>::problem( value_type x_definitions[ _E_n][3 ] )
    :
    Data(),
    M_xlower ( _E_n ),
    M_xupper ( _E_n ),
    M_x0( _E_n )
{
    this->define_problem( x_definitions );
    this->initialize_x0();
}


template<typename Data>
problem<Data>::~problem()
{
}

template<typename Data>
template<typename Defs>
void
problem<Data>::define_problem( Defs x_definitions )
{
    /*    if ( !x_definitions )
    	{
    	    for( int i = 0; i < _E_n; i++ )
    	    {
                __x_definitions[ i][0 ] = 1;
                __x_definitions[ i][1 ] = -gopt::inf;
                __x_definitions[ i][2 ] = +gopt::inf;
                M_xlower(  i  ) = -gopt::inf;
                M_xupper(  i  ) = +gopt::inf;
    	    }
        }
        else*/
    {
        //_E_n = 0;
        for ( int i = 0; i < _E_n; i++ )
        {
            __x_definitions[ i][0 ] = x_definitions[ i][0 ];
            __x_definitions[ i][1 ] = x_definitions[ i][1 ];
            __x_definitions[ i][2 ] = x_definitions[ i][2 ];

            // if( isOn ( i ) )
            {
                M_xlower(  i  ) = __x_definitions[ i][1 ];
                M_xupper(  i  ) = __x_definitions[ i][2 ];
                //_E_n++;
            }
        }

        //_E_nA = _E_n + _E_n * (_E_n - 1) / 2;
    }
}

template<typename Data>
void
problem<Data>::copy_x0_to_x( vector_type& __x )
{
    int offset = 0;

    for ( int i = 0; i < _E_n; i++ )
    {
        if ( isOn( i ) )
            __x(  i-offset  ) = M_x0(  i  );

        else
            offset++;
    }
}


template<typename Data>
void
problem<Data>::copy_x_to_x0( vector_type const& __x )
{
    int offset = 0;

    for ( int i = 0; i < _E_n; i++ )
    {
        if ( isOn ( i ) )
            M_x0(  i  ) = __x(  i-offset  );

        else
            offset++;
    }
}


template<typename Data>
void
problem<Data>::initialize_x0()
{
    for ( int i = 0; i < _E_n; i++ )
    {
        if ( isOn ( i ) )
        {
            if ( __x_definitions[ i][1 ] <= -inf && __x_definitions[ i][2 ] >= inf )
                M_x0(  i  ) = 0;

            if ( __x_definitions[ i][1 ] <= -inf && __x_definitions[ i][2 ] < inf )
                M_x0(  i  ) = __x_definitions[ i][2 ] - value_type( 10 );

            if ( __x_definitions[ i][1 ] > -inf && __x_definitions[ i][2 ] >= inf )
                M_x0(  i  ) = __x_definitions[ i][1 ] + value_type( 10 );

            if ( __x_definitions[ i][1 ] > -inf && __x_definitions[ i][2 ] < inf )
                M_x0(  i  ) = 0.5 * ( __x_definitions[ i][1 ] + __x_definitions[ i][2 ] );
        }

        else
            M_x0(  i  ) = __x_definitions[ i][1 ];
    }
}


template<typename Data>
void
problem<Data>::custom_initialize_x0( vector_type& __custom )
{
    for ( int i = 0; i < _E_n; i++ )
    {
        if ( __custom(  i  ) < __x_definitions[ i][1 ] )
        {
            std::ostringstream __ex;
            __ex << "Initialization Error: initial x^0[ " << i << " ] < xlower[ " << i << " ]\n"
                 << " * x^0[ " << i << " ] = " << __custom( i ) << "\n"
                 << " * xlower[ " << i << " ] = " << __x_definitions[ i][1 ] << "\n";

            throw std::invalid_argument( __ex.str() );
        }

        if ( __custom(  i  ) > __x_definitions[ i][2 ] )
        {
            std::ostringstream __ex;
            __ex << "Initialization Error: initial x^0[ " << i << " ] > xupper[ " << i << " ]\n"
                 << " * x^0[ " << i << " ] = " << __custom( i ) << "\n"
                 << " * xupper[ " << i << " ] = " << __x_definitions[ i][2 ] << "\n";

            throw std::invalid_argument( __ex.str() );
        }

        M_x0(  i  ) = __custom(  i  );
    }
}



//
// Printing routines
//

template<typename Data>
void
problem<Data>::print( vector_type const& _x ) const
{
    int offset = 0;
    std::cout << "x:\n";

    for ( int i = 0; i < _E_n; i++ )
    {
        if ( isOn ( i ) )
            std::cout <<  " x [ " << i << " ] = " <<  _x(  i-offset  );

        else
        {
            std::cout <<  " x [ " << i << " ] = " << __x_definitions[ i][1 ]  <<  "   ( variable \"off\" )";
            offset++;
        }

        std::cout << "\n";
    }
}

template<typename Data>
void
problem<Data>::print_bound_constraints() const
{
    //std::cerr << "\nBound constraints with N = " << N() << "\n";
    for ( int i = 0; i < _E_n; i++ )
    {
        if ( isOn ( i ) )
        {
            if ( __x_definitions[ i][1 ] == -inf )
                std::cerr << std::setw( 8 ) << std::right << " -Inf" <<  " <= ";

            else if ( abs( __x_definitions[ i][1 ] ) < 1e4 )
                std::cerr << std::setw( 8 ) << std::right  << __x_definitions[ i][1 ] << " <= ";

            else
                std::cerr << std::setw( 8 ) << std::right << __x_definitions[ i][1 ] << " <= ";

            std::cerr << "x[" << i << "] <= ";

            if ( __x_definitions[ i][2 ] == inf )
                std::cerr << std::left <<  "+Inf";

            else if ( fabs( __x_definitions[ i][2 ] ) < 1e4 )
                std::cerr << std::left << __x_definitions[ i][2 ];

            else
                std::cerr << std::left << __x_definitions[ i][2 ];
        }

        else
        {
            std::cerr << "x[" << i << " ]  := " << __x_definitions[ i][1 ];
        }

        std::cerr << "\n";
    }
}

template<typename Data>
void
problem<Data>::print_complete( std::string const& __mesg, vector_type const& _x ) const
{
    std::cout << __mesg << "\n";
    print( _x );

    print_bound_constraints();
}

template<typename Data>
void
problem<Data>::print_complete( vector_type const& _x,
                               vector_type const& _s ) const
{
    print_bound_constraints();
    std::cerr << "\n";
    print( _x );

    f_type __fx;
    this->evaluate ( _x, __fx, diff_order<2>() );

    std::cout << "fx.gradient:\n"
              << __fx.gradient( 0 ) << "\n";

    std::cout << "fx.hessian:\n"
              << __fx.hessian( 0 ) << "\n";


    g_type __gx;
    this->evaluate ( _x, __gx, diff_order<2>() );

    for ( int m = 0; m < _E_g; m++ )
    {
        std::cout << "gx(" << m << "):\n"
                  << __gx.value( m ) << "\n";

        std::cout << "gx(" << m << ").gradient:\n"
                  << __gx.gradient( m ) << "\n";

        std::cout << "gx(" << m << ").hessian:\n"
                  << __gx.hessian( m ) << "\n";
    }

    h_type __hx;
    this->evaluate ( _x, __hx, diff_order<2>() );

    for ( int m = 0; m < _E_h; m++ )
    {
        std::cout << "hx(" << m << "):\n"
                  << __hx.value( m ) << "\n";

        std::cout << "hx(" << m << ").gradient:\n"
                  << __hx.gradient( m ) << "\n";

        std::cout << "hx(" << m << ").hessian:\n"
                  << __hx.hessian( m ) << "\n";
    }
}


template<typename Data>
void
problem<Data>::print_stationary_x( vector_type const& _x,
                                   vector_type const& _lambda_l,
                                   vector_type const& _lambda_u ) const
{

    vector_type __sum( _E_n );
    value_type __prod;

    this->print_bound_constraints();

    std::cerr.precision ( 5 );

    std::cerr << "\n\nSolution:\n";
    int offset = 0;

    for ( int i = 0; i < _E_n; i++ )
    {
        if ( isOn ( i ) )
        {
            std::cerr << " x*[ " << i << " ] = " << std::setw( 10 ) << _x(  i-offset  );

            if ( std::abs( _lambda_l(  i-offset  ) ) > 1e-6 )
                std::cerr << "   ( lower bound active )";

            if ( std::abs( _lambda_u(  i-offset  ) ) > 1e-6 )
                std::cerr << "   ( upper bound active )";
        }

        else
        {
            std::cerr << " x*[ " << i << " ] = " << std::setw( 10 ) <<  __x_definitions[ i][1 ];
            std::cerr << "   ( variable \"off\" ) ";
            offset++;
        }

        std::cerr << "\n";
    }

    f_type __f_x;
    this->evaluate( _x, __f_x, diff_order<2>() );
    std::cerr << "f(x*)   = " << __f_x.value( 0 ) << "\n"
              << "f'(x*)  = " << __f_x.gradient( 0 ) << "\n"
              << "f''(x*) = " << __f_x.hessian( 0 ) << "\n";

    std::cerr << "\nStationarity conditions (LS):\n"
              << "  \\nabla f   -   \\lambda_l  +   \\lambda_u  =    residual\n"
              << "----------------------------------------------------------\n";

    offset = 0;

    for ( int i = 0; i < _E_n; i++ )
        if ( isOn ( i ) )
        {
            __sum( i-offset ) = __f_x.gradient( 0, i-offset ) - _lambda_l( i-offset ) + _lambda_u( i-offset );

            std::cerr << " " << std::setw( 10 ) << __f_x.gradient( 0, i-offset ) << " - "
                      << " " << std::setw( 10 ) << _lambda_l( i-offset ) << " + "
                      << " " << std::setw( 10 ) << _lambda_u( i-offset ) << " = "
                      << " " << std::setw( 10 ) <<  __sum( i-offset ) << "\n";

        }

        else
        {
            std::cerr << "   x x x x x x x x x x x x x x x x x x x x x x x x x x x\n";
            offset++;
        }

    std::cerr << "\n =========================================================";
    std::cerr << "\n # Norm of residuals ------------------------>  " ;
    std::cerr << " " << std::setprecision( 2 ) << std::setw( 8 ) <<  norm_2( __sum ) << " #";
    std::cerr << "\n =========================================================\n";

    std::cerr << "\nComplimentarity conditions (LS):\n"
              << " \\lambda_l *     ( l  -   x )    =   0.000000 | \\lambda_u *     ( x  -   u )    =   0.00000\n"
              << "----------------------------------------------|----------------------------------------------\n";

    offset = 0;

    for ( int i = 0; i < _E_n; i++ )
        if ( isOn ( i ) )
        {
            __prod = _lambda_l( i-offset ) * ( __x_definitions[ i][1 ] - _x( i-offset ) );

            std::cerr << " " << std::setw( 9 ) <<  _lambda_l( i-offset ) << " * "
                      << " " << std::setw( 15 ) << __x_definitions[ i][1 ] - _x( i-offset ) << " = "
                      << " " << std::setw( 9 ) << __prod << " |";

            __prod = _lambda_u( i-offset ) * ( _x( i-offset ) - __x_definitions[ i][2 ] );

            std::cerr << " " << std::setw( 9 ) <<  _lambda_u( i-offset ) << " * "
                      << " " << std::setw( 15 ) << _x( i-offset )-__x_definitions[ i][2 ] << " = "
                      << " " << std::setw( 9 ) << __prod << "\n";

        }

        else
        {
            std::cerr << "  x x x x x x x x x x x x x x x x x x x x x x |"
                      << "  x x x x x x x x x x x x x x x x x x x x x x\n";
            offset++;
        }
}


template<typename Data>
void
problem<Data>::print_stationaryN_x( vector_type const& _x,
                                    vector_type const& _s,
                                    vector_type const& _lambda_h,
                                    vector_type const& _lambda_g ) const
{
    std::cerr << "\n\nVariables and (variable) constraints:\n";

    this->print_bound_constraints();

    std::cerr << "\n\nSolution:\n\n";

    std::cerr.precision ( 5 );

    int offset = 0;

    for ( int i = 0; i < _E_n; i++ )
    {
        if ( isOn ( i ) )
        {
            std::cerr << " x*[ " << i << " ] = " << std::setw( 10 ) << _x(  i-offset  );

            if ( _x(  i-offset  ) < M_xlower(  i-offset  ) ||
                    _x(  i-offset  ) > M_xupper(  i-offset  ) )
            {
                std::cerr << " <- ( violated ! )";
            }

            else
            {
                value_type lam_l = _lambda_g(  _E_g+i-offset  );
                value_type lam_u = _lambda_g(  _E_g+_E_n+i-offset  );

                if ( std::abs( lam_l ) > 1e-4 )
                    std::cerr << "    ( lower bound active, lambda_l " << std::setw( 10 ) << lam_l;

                if ( std::abs( lam_u ) > 1e-4 )
                    std::cerr << "    ( upper bound active, lambda_u " << std::setw( 10 ) << lam_u;
            }
        }

        else
        {
            std::cerr << " x*[ " << i << " ] = " << std::setw( 10 ) <<  __x_definitions[ i][1 ];
            std::cerr << "    ( variable \"off\" ) ";
            offset++;
        }

        std::cerr << "\n";
    }

    // Functional
    f_type __fx;
    this->evaluate ( _x, __fx, diff_order<1>() );

    std::cerr << " f(x*) = " << std::setw( 10 ) << __fx.value( 0 ) << "\n";

    // inequalities
    g_type __gx;
    this->evaluate ( _x, __gx, diff_order<0>() );

    for ( int m = 0; m < _E_g; m++ )
    {
        if ( __gx.value( m ) < 0 )
        {
            if ( std::abs( std::fabs( _lambda_g(  m  ) ) ) > 1e-3 )
            {
                std::cerr << " g[" << m << "](x*) = " << std::setw( 10 ) << __gx.value( m )
                          << "    ( satisfied - active, lambda = " << std::setw( 10 ) << _lambda_g(  m  );
            }

            else
            {
                std::cerr << " g[" << m << "](x*) = " << std::setw( 10 ) << __gx.value( m )
                          << "    ( satisfied - inactive )";

            }
        }

        else
        {
            std::cerr << " g[" << m << "](x*) = " << std::setw( 10 ) << __gx.value( m )
                      << " <- ( violated ! )";
        }

        std::cerr << "\n";
    }

    // equalities
    h_type __hx;
    this->evaluate ( _x, __hx, diff_order<0>() );

    for ( int m = 0; m < _E_h; m++ )
    {
        std::cerr << " h[" << m << "](x*) = " << std::setw( 10 ) << __hx.value( m );

        if ( std::abs( __hx.value( m ) < 1e-3 ) )
            std::cerr << "    ( satisfied )";

        else
            std::cerr << " <- ( violated ! )";

        std::cerr << "\n";
    }
}
}


#endif

