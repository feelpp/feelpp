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
   \file adfunctions.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
template <class Expr> class @NAME@
{
public:
    enum { nvar = Expr::nvar };
    typedef typename Expr::value_type value_type;
protected:
    @NAME@ () {}

    Expr expr_;
public:

    @NAME@ ( const Expr & expr ) : expr_( expr )
    {
        ;
    }


    template<int isFundamental, typename Expr_>
    struct Value
    {

        typedef typename Expr_::value_type value_type;
        static value_type signum( value_type a )
        {
            if ( abs( a ) < std::numeric_limits<value_type>::epsilon() )
                return value_type( 0 );

            return a/abs( a );
        }
        static value_type value( Expr_ const& expr_ )
        {
            return @FCT@( expr_.value() );
        }
        static value_type grad( Expr_ const& expr_, int __i )
        {
            return @GRDI@;
        }
        static value_type hessian( Expr_ const& expr_, int __i, int __j )
        {
            return @HESSIJ@;
        }
    };
    template<typename Expr_>
    struct Value<true, Expr_>
    {
        typedef typename Expr_::value_type value_type;
        static value_type signum( value_type a )
        {
            if ( std::abs( a ) < std::numeric_limits<value_type>::epsilon() )
                return value_type( 0 );

            return a/std::abs( a );
        }
        static value_type value( Expr_ const& expr_ )
        {
            return std::@FCT@( expr_.value() );
        }
        static value_type grad( Expr_ const& expr_, int __i )
        {
            return @GRDI_STD@;
        }
        static value_type hessian( Expr_ const& expr_, int __i, int __j )
        {
            return @HESSIJ_STD@;
        }
    };

    inline value_type value() const
    {
        return Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::value( expr_ );
    }
    inline value_type grad( int __i ) const
    {
        return  Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::grad( expr_, __i );
    }
    inline value_type hessian( int __i, int __j ) const
    {
        return  Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::hessian( expr_, __i, __j );
    }
    inline bool deps( int __i ) const
    {
        return expr_.deps( __i ) ;
    }
};

template <class Expr> inline ADExpr< @NAME@< ADExpr<Expr> > >
@FCT@ ( const ADExpr<Expr>& expr )
{
    typedef @NAME@< ADExpr<Expr> > expr_t;
    return ADExpr< expr_t >(  expr_t( expr ) );
}

template <class T, int Nvar, int Order, int Var> inline ADExpr< @NAME@< ADType<T, Nvar, Order, Var> > >
@FCT@ ( const ADType<T, Nvar, Order, Var>& x )
{
    typedef @NAME@< ADType<T, Nvar, Order, Var> > expr_t;
    return ADExpr< expr_t >(  expr_t( x ) );
}
