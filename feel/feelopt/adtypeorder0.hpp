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
   \file adtypeorder0.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
#ifndef __ADTypeOrder0_H
#define __ADTypeOrder0_H 1

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
  @version $Id: ADType.hpp,v 1.7 2002/06/06 02:22:52 prudhomm Exp $
*/
template<typename T, int Nvar, int Var>
class ADType<T, Nvar, 0, Var>
{
public:

    enum { nvar = Nvar };
    enum { order = 1 };

    typedef ADVariable<Var> variable_type;

    typedef ADType<T, Nvar, 0, Var> This;

    typedef T value_type;

    template<typename NumT, int NumVar, int Order, int VarNum> friend class ADType;

    operator T()
    {
        return M_val;
    }

    ADType (  T __v = 0 )
        :
        M_val (  __v )
    {
    }

    template<int VarNum>
    ADType( ADType<T,Nvar,0,VarNum> const& sad )
        :
        M_val( sad.M_val )
    {
    }

    template<typename ExprT>
    ADType ( const ADExpr<ExprT>& expr )
        :
        M_val( 0 )
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


    //@}

    /** @name Mutators
     */
    //@{
    value_type& value()
    {
        return M_val;
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

#define AD_UNARY_OP( op )                       \
    This& operator op ( value_type val )        \
    {                                           \
        M_val op val;                          \
        return *this;                           \
    }
    AD_UNARY_OP( += );
    AD_UNARY_OP( -= );
    AD_UNARY_OP( *= );
    AD_UNARY_OP( /= );

#undef AD_UNARY_OP

    This& operator += ( This const& sad )
    {
        M_val += sad.M_val;
        return *this;
    }
    This& operator -= ( This const& sad )
    {
        M_val -= sad.M_val;
        return *this;
    }


    This& operator *= ( This const& sad )
    {
        M_val *= sad.M_val;
        return *this;
    }

    This& operator /= ( This const& sad )
    {
        M_val /= sad.M_val;
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
private:

    value_type M_val;

};
template<typename T,int Nvar, int Var>
ADType<T, Nvar, 0, Var>&
ADType<T, Nvar, 0, Var>::operator=( value_type const& val )
{
    M_val = val;
    return *this;
}

template<typename T,int Nvar, int Var>
ADType<T, Nvar, 0, Var>&
ADType<T, Nvar, 0, Var>::operator=( This const& sad )
{
    M_val = sad.M_val;
    return *this;
}
template<typename T,int Nvar, int Var>
template <class ExprT>
ADType<T,Nvar, 0, Var> &
ADType<T,Nvar, 0, Var>::operator=( const ADExpr<ExprT>& expr )
{
    M_val = expr.value();
    return *this;
}


}
//------------------------------- AD ostream operator ------------------------------------------
template <class T, int Nvar, int Var>
std::ostream&
operator << ( std::ostream& os, const Feel::ADType<T, Nvar, 0, Var>& a )
{
    os.setf( std::ios::fixed,std::ios::floatfield );
    os.width( 12 );
    os << "value    = " << a.value() << "  \n";
    return os;
}

#endif
