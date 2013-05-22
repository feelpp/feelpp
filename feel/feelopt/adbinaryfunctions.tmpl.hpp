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
   \file adbinaryfunctions.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
#ifndef __AdBinaryFunctions_H
#define __AdBinaryFunctions_H 1


# include <boost/preprocessor/comparison/less.hpp>
# include <boost/preprocessor/logical/and.hpp>
# include <boost/preprocessor/control/if.hpp>
# include <boost/preprocessor/list/at.hpp>
# include <boost/preprocessor/list/cat.hpp>
# include <boost/preprocessor/list/for_each_product.hpp>
# include <boost/preprocessor/logical/or.hpp>
# include <boost/preprocessor/tuple/to_list.hpp>
# include <boost/preprocessor/tuple/eat.hpp>
# include <boost/preprocessor/facilities/empty.hpp>
# include <boost/preprocessor/punctuation/comma.hpp>
# include <boost/preprocessor/facilities/identity.hpp>
# include <boost/preprocessor/stringize.hpp>
namespace Feel
{
namespace detail
{
# /* Accessors for the operator datatype. */
# define AD_BINARY_FUNCTION_NAME(O)           BOOST_PP_TUPLE_ELEM(5, 0, O)
# define AD_BINARY_FUNCTION_SYMBOL(O)         BOOST_PP_TUPLE_ELEM(5, 1, O)
# define AD_BINARY_FUNCTION_COMMENT(O)        BOOST_PP_TUPLE_ELEM(5, 2, O)
# define AD_BINARY_FUNCTION_GRAD(O)           BOOST_PP_TUPLE_ELEM(5, 3, O)
# define AD_BINARY_FUNCTION_HESS(O)           BOOST_PP_TUPLE_ELEM(5, 4, O)

# /* List of applicative operators. */
# define AD_BINARY_FUNCTIONS                                            \
    BOOST_PP_TUPLE_TO_LIST(                                             \
                           1,                                           \
                           (                                            \
                            ( ADBinFunPow   , pow, "Power function", (expr2_.value() * expr1_.grad(__i) * math::pow(expr1_.value(),expr2_.value()-1)), (expr2_.value() * ( (expr2_.value()-1) * expr1_.grad(__i) * expr1_.grad(__j)  * math::pow(expr1_.value(),expr2_.value()-2) +  expr1_.hessian(__i,__j) * math::pow(expr1_.value(),expr2_.value()-1))) ) \
                                                                        ) \
                                                                        ) \
    /**/
#

/**
   
   \brief brief description

   @author Christophe Prud'homme
   @see
*/
template <class Expr1, class Expr2> class @NAME@
{
public:

    enum { nvar = Expr1::nvar };
    //enum { nvar2 = Expr2::nvar };

    typedef typename Expr1::value_type value_type;

protected:
    @NAME@ () {}

    Expr1 expr1_;
    Expr2 expr2_;

public:

    @NAME@ ( const Expr1 & expr1, const Expr2& expr2 ) : expr1_( expr1 ), expr2_( expr2 )
    {
        ;
    }


    template<int isFundamental, typename Expr1_, typename Expr2_>
    struct Value
    {

        typedef typename Expr1_::value_type value_type;

        static value_type value( Expr1_ const& expr1_, Expr2_ const& expr2_ )
        {
            return @FCT@( expr1_.value(), expr2_.value() );
        }
        static value_type grad( Expr1_ const& expr1_, Expr2_ const& expr2_, int __i )
        {
            return @GRDI@;
        }
        static value_type hessian( Expr1_ const& expr1_, Expr2_ const& expr2_, int __i, int __j )
        {
            return @HESSIJ@;
        }
    };
    template<typename Expr1_, typename Expr2_>
    struct Value<true, Expr1_, Expr2_>
    {
        typedef typename Expr1_::value_type value_type;

        static value_type value( Expr1_ const& expr1_, Expr2_ const& expr2_ )
        {
            return std::@FCT@( expr1_.value(), expr2_.value() );
        }
        static value_type grad( Expr1_ const& expr1_, Expr2_ const& expr2_, int __i )
        {
            return @GRDI_STD@;
        }
        static value_type hessian( Expr1_ const& expr1_, Expr2_ const& expr2_, int __i, int __j )
        {
            return @HESSIJ_STD@;
        }
    };

    inline value_type value() const
    {
        return Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr1,Expr2>::value( expr1_, expr2_ );
    }
    inline value_type grad( int __i ) const
    {
        return  Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr1,Expr2>::grad( expr1_, expr2_, __i );
    }
    inline value_type hessian( int __i, int __j ) const
    {
        return  Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr1, Expr2>::hessian( expr1_, expr2_, __i, __j );
    }
    inline bool deps( int __i ) const
    {
        return expr1_.deps( __i ) || expr2_.deps( __i );
    }
};

template <class Expr1, class Expr2>
inline
ADExpr< @NAME@< ADExpr<Expr1>, ADExpr<Expr2> > >
@FCT@ ( const ADExpr<Expr1>& expr1, const ADExpr<Expr2>& expr2 )
{
    typedef @NAME@< ADExpr<Expr1>, ADExpr<Expr2> > expr_t;
    return ADExpr< expr_t >(  expr_t( expr1, expr2 ) );
}

template <class T, int Nvar, int Order, int Var>
inline
ADExpr< @NAME@< ADType<T, Nvar, Order, Var>, ADType<T, Nvar, Order, Var> > >
@FCT@ ( const ADType<T, Nvar, Order, Var>& x, const ADType<T, Nvar, Order, Var>& y )
{
    typedef @NAME@< ADType<T, Nvar, Order, Var>, ADType<T, Nvar, Order, Var> > expr_t;
    return ADExpr< expr_t >(  expr_t( x, y ) );
}

template <class W, class T, int Nvar, int Order, int Var>
inline
ADExpr< AdFuncPow< ADType<T, Nvar, Order, Var>, ADCst<W> > >
pow ( const ADType<T, Nvar, Order, Var>& x, W y )
{
    typedef AdFuncPow< ADType<T, Nvar, Order, Var>, ADCst<W> > expr_t;
    ADCst<W> y1 ( y );
    return ADExpr< expr_t >(  expr_t( x, y1 ) );
}
}
}
#endif /* __AdBinaryFunctions_H */
