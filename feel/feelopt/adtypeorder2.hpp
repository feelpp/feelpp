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
   \file adtypeorder2.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
#ifndef __ADTypeOrder2_H
#define __ADTypeOrder2_H 1

#include <boost/numeric/ublas/symmetric.hpp>
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

  So far only Order = 2 has been implemented
  here is an example
  <pre>
  #include <ADType.hpp>

  ADType<double,3,2, Var> x(1.,0 );
  ADType<double,3,2, Var> y(2.,1 );
  ADType<double,3,2, Var> z(3.,2 );
  Feel:::ADType<double,3,2, Var> __g = x/(y*z);
  std::cout << "g=" << __g << "\n";
  <pre>

  @author Christophe Prud'homme
  @see
  @version $Id: ADType.hpp,v 1.7 2002/06/06 02:22:52 prudhomm Exp $
*/
template<typename T, int Nvar, int Var>
class ADType<T, Nvar, 2, Var>
{
public:

    enum { nvar = Nvar };
    enum { order = 2 };
    enum { depvar = Var };

    typedef ADVariable<Var> variable_type;
    //typedef typename Feel::STL::SFlatten<typename SListGenerator<Nvar,ADVariable>::list_type>::Result variable_list_type;

    typedef T value_type;
    typedef boost::numeric::ublas::vector<value_type> gradient_type;
    //typedef typename STinyVector<value_type, Nvar>::type gradient_type;
    //typedef typename STinyVector<bool, Nvar>::type deps_type;
    typedef std::vector<bool> deps_type;

    //typedef typename STinyVector<typename STinyVector<value_type,Nvar>::type, Nvar>::type hessian_type;
    typedef boost::numeric::ublas::symmetric_matrix<value_type, boost::numeric::ublas::upper> hessian_type;

    typedef ADType<T,Nvar,2, Var> This;

    typedef std::set<int> dependency_type;

    template<typename NumT, int NumVar, int Order, int VarNum> friend class ADType;

    /** @name Constructors/Destructor
     */
    //@ {

    ADType( value_type val  = value_type( 0 ) )
        :
        M_val( val ),
        M_grad( nvar ),
        M_hessian( nvar, nvar ),
        __dep( false ),
        M_deps()
    {
        M_grad = boost::numeric::ublas::zero_vector<double>( nvar );
        M_hessian = boost::numeric::ublas::zero_matrix<double>( nvar, nvar );

        if ( depvar != -1 )
        {
            M_grad[ depvar ] = 1.;
            __dep[ depvar ] = true;
        }
    }

    template<int VarNum>
    ADType( ADType<T,Nvar,2,VarNum> const& sad )
        :
        M_val( sad.M_val ),
        M_grad( sad.M_grad ),
        M_hessian( sad.M_hessian ),
        __dep( sad.__dep ),
        M_deps()
    {

    }
    template<typename ExprT>
    ADType ( const ADExpr<ExprT>& expr )
        :
        M_val( 0 ),
        M_grad( nvar ),
        M_hessian( nvar, nvar ),
        __dep( false ),
        M_deps()
    {
        *this = expr;
    }
    //@}

    /** @name Accessors
     */
    //@{

    //! returns the value
    value_type value() const
    {
        return M_val;
    }

    //! tells if the gradient or the hessian depends on the ith variable
    bool deps( int __i ) const
    {
        //ENSURE( __i >= 0 && __i < nvar );
        return __dep[ __i ];
    }

    //! get the ith corrdinate of the gradient
    value_type grad( int __i ) const
    {
        //ENSURE( __i >= 0 && __i < nvar );
        return M_grad( __i );
    }

    gradient_type grad() const
    {
        return M_grad;
    }

    hessian_type hessian() const
    {
        return M_hessian;
    }

    //! returns the hessian value at position (i,j)
    value_type hessian( int __i, int __j ) const
    {
        //ENSURE( __i >= 0 && __i < nvar && __j >= 0 && __j < nvar );
        return M_hessian( __i, __j );
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
        //ENSURE( __i >= 0 && __i < nvar );
        return __dep[ __i ];
    }
    //! get the ith corrdinate of the gradient
    value_type& grad( int __i )
    {
        //ENSURE( __i >= 0 && __i < nvar );
        return M_grad( __i );
    }

    gradient_type& grad()
    {
        return M_grad;
    }

    hessian_type& hessian()
    {
        return M_hessian;
    }

    //! returns the hessian value at position (i,j)
    value_type& hessian( int __i, int __j )
    {
        //ENSURE( __i >= 0 && __i < nvar && __j >= 0 && __j < nvar );
        return M_hessian( __i, __j );
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

    This& operator += ( This const& sad )
    {
        M_val += sad.M_val;
        M_grad += sad.M_grad;
        M_hessian += sad.M_hessian;
        return *this;
    }
    This& operator -= ( This const& sad )
    {
        M_val -= sad.M_val;
        M_grad -= sad.M_grad;
        M_hessian -= sad.M_hessian;
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


    value_type M_val;
    gradient_type M_grad;
    hessian_type M_hessian;

    deps_type __dep;
    dependency_type M_deps;
};

template<typename T,int Nvar, int Var>
ADType<T, Nvar, 2, Var>&
ADType<T, Nvar, 2, Var>::operator=( value_type const& val )
{
    __dep = false;
    M_val = val;
    M_grad = boost::numeric::ublas::zero_vector<double>( nvar );
    M_hessian = boost::numeric::ublas::zero_matrix<double>( nvar, nvar );;
    return *this;
}

template<typename T,int Nvar, int Var>
ADType<T, Nvar, 2, Var>&
ADType<T, Nvar, 2, Var>::operator=( This const& sad )
{
    __dep = sad.__dep;
    M_val = sad.M_val;
    M_grad = sad.M_grad;
    M_hessian = sad.M_hessian;
    return *this;
}
template<typename T,int Nvar, int Var>
template <class ExprT>
ADType<T,Nvar, 2, Var> &
ADType<T,Nvar, 2, Var>::operator=( const ADExpr<ExprT>& expr )
{

#if 0
    M_val = expr.value();

    typedef typename Feel::STL::SNoDuplicates<typename Feel::STL::SFlatten<typename ADExpr<ExprT>::Var>::Result>::Result Expr;
    typedef typename Feel::Meta::IF<Feel::STL::SIndexOf<Expr,Variable<-1> >::value == -1,
            FOR_EACH<Expr>,
            depManager
            >::Result ForEachExpr;

    // should generate the complement to VarList of Expr
    M_grad = boost::numeric::ublas::zero_vector<double>( nvar );
    ForEachExpr( *this, expr );

    //FOR_EACH<FOR_EACH<Expr> >( M_grad, expr );
    return *this;
#endif
    M_val = expr.value();
#if 0
    M_deps.clear();

    for ( int __i = 0; __i < nvar; ++__i )
    {
        __dep( __i ) = expr.deps( __i );

        if ( __dep( __i ) )
        {
            //std::cerr << "depends on " << __i << "\n";
            M_deps.insert( __i );
        }
    }

    if ( ! M_deps.empty() && M_deps.size() < nvar )
    {
        //M_grad = 0;
        //M_hessian = M_grad;

        dependency_type::iterator __begin = M_deps.begin();
        dependency_type::iterator __end = M_deps.end();

        dependency_type::iterator __it = __begin;
        dependency_type::iterator __it2;

        while ( __it != __end )
        {
            M_grad( *__it )= expr.grad( *__it );
            __it2 = __begin;

            while ( __it2 != __it )
            {
                M_hessian( *__it, *__it2 ) = expr.hessian( *__it, *__it2 );
                //M_hessian( *__it2, *__it ) = M_hessian( *__it][ *__it2 ];

                ++__it2;
            }

            M_hessian( *__it, *__it ) = expr.hessian( *__it, *__it );
            ++__it;
        }
    }

    else
    {
#endif

        for ( int __i = 0; __i < nvar; ++__i )
        {
            M_grad( __i )= expr.grad( __i );
            M_hessian( __i, __i ) = expr.hessian( __i, __i );

            for ( int __j = 0; __j < __i; ++__j )
            {
                M_hessian( __i, __j ) = expr.hessian( __i, __j );
                //M_hessian( __j, __i ) = M_hessian[ __i ][ __j ];
            }
        }

#if 0
    }

#endif
    return *this;

}

} // Feel

//------------------------------- AD ostream operator ------------------------------------------

template <class T, int Nvar, int Var>
std::ostream&
operator << ( std::ostream& os, const Feel::ADType<T, Nvar, 2, Var>& a )
{
    os.setf( std::ios::fixed,std::ios::floatfield );
    //os.width(12);
    os << "value    = " << a.value() << "  \n"
       << "gradient = " << a.grad() << "\n"
       << "hessian  = " << a.hessian() << "\n";
#if 0

    for ( int __i = 0 ; __i < Nvar; ++__i )
    {
        for ( int __j = 0 ; __j < Nvar; ++__j )
            os << a.hessian( __i, __j ) << " ";

        os << "\n";
    }

    os << "]\n";
#endif
    return os;
}

#endif
