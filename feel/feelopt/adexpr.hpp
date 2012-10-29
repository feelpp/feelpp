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
   \file adexpr.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
#ifndef __ADExpr_H
#define __ADExpr_H 1

#include <boost/type_traits/is_fundamental.hpp>

namespace Feel
{
/*!
  \class ADExpr
  \brief brief description

  @author Christophe Prud'homme
  @see
*/
template< typename Expr>
class ADExpr
{
public:


    /** @name Typedefs
     */
    //@{
    enum { nvar = Expr::nvar };

    typedef Expr expression_type;
    typedef typename expression_type::value_type value_type;

    typedef ADExpr<Expr> This;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit ADExpr( expression_type const & expr )
        :
        __expression( expr )
    {

    }
    ~ADExpr() {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    value_type value() const
    {
        return __expression.value();
    }

    value_type grad( int __ith ) const
    {
        return __expression.grad( __ith );
    }

    value_type hessian( int __i, int __j ) const
    {
        return __expression.hessian( __i, __j );
    }

    bool deps( int i ) const
    {
        return __expression.deps( i );
    }

    operator bool() const
    {
        return __expression.value() != value_type( 0 );
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{


    //@}



protected:

    ADExpr()
    {
        ;
    }

    expression_type __expression;

private:

};

//------------------------------- AD constant ------------------------------------------
template < class T >
class ADCst
{
public:

    typedef T value_type;

protected:
    ADCst() {}

    const T constant_;

public:
    explicit ADCst( const T& value ) : constant_( value )
    {
        ;
    }

    value_type value()     const
    {
        return constant_;
    }
    value_type grad( int __i ) const
    {
        return 0;
    }
    value_type hessian( int __i, int __j ) const
    {
        return 0;
    }

    bool deps( int ) const
    {
        return false;
    }

};

//------------------------------- AD unary + ------------------------------------------
template < class T >
class ADUnaryPlus
{
public:

    enum { nvar = T::nvar };
    typedef typename T::value_type value_type;

protected:
    ADUnaryPlus() {}

    const T& expr_;

public:
    ADUnaryPlus( const T& value ) : expr_( value )
    {
        ;
    }

    value_type value()     const
    {
        return expr_.value();
    }
    value_type grad( int __i ) const
    {
        return expr_.grad( __i );
    }
    value_type hessian( int __i, int __j ) const
    {
        return expr_.hessian( __i, __j );
    }

    bool deps( int i ) const
    {
        return expr_.deps( i );
    }
};
template <class T> inline
ADExpr< ADUnaryPlus< ADExpr<T> > >
operator + ( const ADExpr<T>& expr )
{
    typedef ADUnaryPlus< ADExpr<T> > expr_t;

    return ADExpr< expr_t >( expr_t( expr ) );
}


//------------------------------- AD unary - ------------------------------------------
template < class T >
class ADUnaryMinus
{
public:
    enum { nvar = T::nvar };
    typedef typename T::value_type value_type;
protected:
    ADUnaryMinus() {}

    const T& expr_;

public:
    ADUnaryMinus( const T& value ) : expr_( value )
    {
        ;
    }

    value_type value()     const
    {
        return - expr_.value();
    }
    value_type grad( int __i ) const
    {
        return - expr_.grad( __i );
    }
    value_type hessian( int __i, int __j ) const
    {
        return - expr_.hessian( __i, __j );
    }

    bool deps( int i ) const
    {
        return expr_.deps( i );
    }
};



template <class T> inline
ADExpr< ADUnaryMinus< ADExpr<T> > >
operator - ( const ADExpr<T>& expr )
{
    typedef ADUnaryMinus< ADExpr<T> > expr_t;

    return ADExpr< expr_t >( expr_t( expr ) );
}

}



#endif /* __ADExpr_H */

