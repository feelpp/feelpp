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
   \file adtypeorder1.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
#ifndef __ADTypeOrder1_H
#define __ADTypeOrder1_H 1


#include <boost/numeric/ublas/vector.hpp>

#include <boost/numeric/ublas/io.hpp>

namespace Feel
{
/*!
  \class ADType
  \brief Automatically Diffentiated Numerical Type

  The idea is based on the work of Nicolas Di Cesare
  who wrote the FAD<Expr> class. However here there
  are some differences: we can into account the depedencies
  and we compute 1st order or 1st and 2nd order derivatives.

  * T is the numerical Type
  * Nvar is the number of variables
  * Order the order of derivation

  @author Christophe Prud'homme
  @see
*/
template<typename T, int Nvar, int Var>
class ADType<T, Nvar, 1, Var>
{
public:

    enum { nvar = Nvar };
    enum { order = 1 };

    typedef ADVariable<Var> variable_type;

    typedef T value_type;

    typedef boost::numeric::ublas::vector<value_type> gradient_type;
    //typedef typename STinyVector<value_type, Nvar>::type gradient_type;
    //typedef typename STinyVector<bool, Nvar>::type deps_type;
    typedef std::vector<int> deps_type;
    typedef ADType<T,Nvar, 1, Var> This;
    typedef std::set<int> dependency_type;

    template<typename NumT, int NumVar, int Order, int VarNum> friend class ADType;

    /** @name Constructors/Destructor
     */
    //@ {
    ADType()
        :
        M_val( 0. ),
        M_grad( nvar ),
        __dep( nvar, false ),
        __deps()
    {}

    ADType( value_type val )
        :
        M_val( val ),
        M_grad( nvar ),
        __dep( nvar, false ),
        __deps()
    {
        M_grad = boost::numeric::ublas::zero_vector<value_type>( nvar );

        if ( Var != -1 )
        {
            __dep[Var]=true;
            M_grad( Var ) = 1.;
        }
    }

    template<int VarNum>
    ADType( ADType<T,Nvar,1,VarNum> const& sad )
        :
        M_val( sad.M_val ),
        M_grad( sad.M_grad ),
        __dep( sad.__dep ),
        __deps()
    {
    }

    template<typename ExprT>
    ADType ( const ADExpr<ExprT>& expr )
        :
        M_val( 0 ),
        M_grad( nvar ),
        __dep( nvar, false ),
        __deps()
    {
        *this = expr;
    }
    //@}

    /** @name Accessors
     */
    //@{
    value_type value() const
    {
        return M_val;
    }
    deps_type deps() const
    {
        return __dep;
    }

    //! tells if the gradient or the hessian depends on the ith variable
    bool deps( int __i ) const
    {
        ENSURE( __i >= 0 && __i < nvar );
        return __dep[ __i ];
    }

    value_type grad( int __i ) const
    {
        return M_grad( __i );
    }

    gradient_type grad() const
    {
        return M_grad;
    }

    //@}

    /** @name Mutators
     */
    //@{
    value_type& value()
    {
        return M_val;
    }
    //! tells if the gradient or the hessian depends on the ith variable
    bool& deps( int __i )
    {
        ENSURE( __i >= 0 && __i < nvar );
        return __dep[ __i ];
    }
    value_type& grad( int __i )
    {
        return M_grad( __i );
    }

    gradient_type& grad()
    {
        return M_grad;
    }

    //@}

    /** @name Operators
     */
    //@{

    This& operator=( T const& );
    This& operator=( This const& );
    template <class ExprT> This& operator=( const ADExpr<ExprT>& expr );

    //! implements unary + \f$ +x \f$
    ADExpr< ADUnaryPlus< This > > operator+ () const
    {
        typedef ADUnaryPlus<This> expr_t;
        return ADExpr<expr_t> ( expr_t ( *this ) );
    }

    //! implements unary - \f$ -x \f$
    ADExpr< ADUnaryMinus< This > > operator- () const
    {
        typedef ADUnaryMinus<This> expr_t;
        return ADExpr<expr_t> ( expr_t ( *this ) );
    }

#define AD_UNARY_OP( op )                         \
      This& operator op ( value_type val )         \
      {                                            \
        M_val op val;                              \
        return *this;                              \
      }

    AD_UNARY_OP( += );
    AD_UNARY_OP( -= );
    AD_UNARY_OP( *= );
    AD_UNARY_OP( /= );

#undef AD_UNARY_OP

    This& operator += ( This const& sad )
    {
        M_val += sad.M_val;
        M_grad += sad.M_grad;
        return *this;
    }
    This& operator -= ( This const& sad )
    {
        M_val -= sad.M_val;
        M_grad -= sad.M_grad;
        return *this;
    }


    This& operator *= ( This const& sad )
    {
        M_val *= sad.M_val;
        M_grad = M_grad * sad.M_val + M_val * sad.M_grad;
        return *this;
    }

    This& operator /= ( This const& sad )
    {
        M_val /= sad.M_val;
        M_grad = ( M_grad * sad.M_val - M_val * sad.M_grad ) / ( sad.M_val * sad.M_val );
        return *this;
    }

    template<typename Expr>
    This& operator += ( ADExpr<Expr> const& sad )
    {
        *this = *this + sad;
        return *this;
    }

    template<typename Expr>
    This& operator -= ( ADExpr<Expr> const& sad )
    {
        *this = *this - sad;
        return *this;
    }

    template<typename Expr>
    This& operator *= ( ADExpr<Expr> const& sad )
    {
        *this = *this * sad;
        return *this;
    }

    template<typename Expr>
    This& operator /= ( ADExpr<Expr> const& sad )
    {
        *this = *this / sad;
        return *this;
    }
    //@}
protected:

private:
    value_type M_val;
    gradient_type M_grad;

    deps_type __dep;
    dependency_type __deps;
};

template<typename T,int Nvar, int Var>
ADType<T, Nvar, 1, Var>&
ADType<T, Nvar, 1, Var>::operator=( value_type const& val )
{
    M_val = val;
    M_grad = boost::numeric::ublas::zero_vector<value_type>( nvar );
    __dep = false;
    return *this;
}

template<typename T,int Nvar, int Var>
ADType<T, Nvar, 1, Var>&
ADType<T, Nvar, 1, Var>::operator=( This const& sad )
{
    M_val = sad.M_val;
    M_grad = sad.M_grad;
    __dep = sad.__dep;
    return *this;
}
template<typename T,int Nvar, int Var>
template <class ExprT>
ADType<T,Nvar, 1, Var> &
ADType<T,Nvar, 1, Var>::operator=( const ADExpr<ExprT>& expr )
{
    M_val = expr.value();
#if 0
    __deps.clear();

    for ( int __i = 0; __i < nvar; ++__i )
    {

        __dep( __i ) = expr.deps( __i );

        if ( __dep( __i ) )
        {
            //std::cerr << "depends on " << __i << "\n";
            __deps.insert( __i );
        }

    }

    if ( ! __deps.empty() && __deps.size() < nvar )
    {
        //M_grad = 0;
        //__hessian = M_grad;

        dependency_type::iterator __begin = __deps.begin();
        dependency_type::iterator __end = __deps.end();

        dependency_type::iterator __it = __begin;

        while ( __it != __end )
        {
            M_grad( *__it )= expr.grad( *__it );
            ++__it;
        }
    }

    else
    {
#endif

        for ( int __i = 0; __i < nvar; ++__i )
        {
            M_grad( __i )= expr.grad( __i );
        }

#if 0
    }

#endif
    return *this;
}



}


//------------------------------- AD ostream operator ------------------------------------------
template <class T, int Nvar, int Var>
std::ostream&
operator << ( std::ostream& os, const Feel::ADType<T, Nvar, 1, Var>& a )
{
    os.setf( std::ios::fixed,std::ios::floatfield );
    os.width( 12 );
    os << "value    = " << a.value() << "  \n";
    os << "gradient = " << a.grad() << "\n";
    return os;
}


#endif /* __ADTypeOrder1_H */
