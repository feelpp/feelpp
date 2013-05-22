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
   
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4


//------------------------------- AD operators ------------------------------------------

template<class E1, class E2>
inline
ADExpr< @STYPE@< ADExpr<E1>, ADExpr<E2> > >
@OP@  ( const ADExpr<E1> &v, const ADExpr<E2> &w )
{
    typedef @STYPE@<ADExpr<E1>, ADExpr<E2> > expr_t;
    return ADExpr<expr_t> ( expr_t ( v , w ) );
}
template<class E, int Nvar, int order, int Var>
inline
ADExpr<@STYPE@<ADExpr<E>,ADType<typename E::value_type, Nvar, order, Var> > >
@OP@ ( const ADExpr<E> &e,const ADType<typename E::value_type, Nvar, order, Var>& v )
{
    typedef typename E::value_type A;
    typedef @STYPE@<ADExpr<E>,ADType<typename E::value_type, Nvar, order, Var> > expr_t;
    return ADExpr<expr_t>( expr_t ( e, v ) );
}

template<typename A, int Nvar, int order, int Var1, int Var2>
inline
ADExpr<@STYPE@<ADType<A, Nvar, order, Var1>,ADType<A, Nvar, order, Var2> > >
@OP@ ( const ADType<A, Nvar, order, Var1> &e1,const ADType<A, Nvar, order, Var2>& e2 )
{
    typedef @STYPE@<ADType<A, Nvar, order, Var1>,ADType<A, Nvar, order, Var2> > expr_t;
    return ADExpr<expr_t>( expr_t ( e1 , e2 ) );
}

template<class E, int Nvar, int order, int Var>
inline
ADExpr<@STYPE@<ADType<typename E::value_type, Nvar, order, Var>,ADExpr<E> > >
@OP@ ( const ADType<typename E::value_type, Nvar, order, Var> &v, const ADExpr<E> &e )
{
    typedef typename E::value_type A;
    typedef @STYPE@<ADType<typename E::value_type, Nvar, order, Var>,ADExpr<E> > expr_t;
    return ADExpr<expr_t> ( expr_t ( v , e ) );
}

template<typename A, typename C, int Nvar, int order, int Var>
inline
ADExpr<@STYPE@<ADCst<C>,ADType<A, Nvar, order, Var> > >
@OP@ ( C a, const ADType<A, Nvar, order, Var> &e )
{
    typedef @STYPE@<ADCst<C>,ADType<A, Nvar, order, Var> > expr_t;
    return ADExpr<expr_t> ( expr_t ( ADCst<C>( a ), e  ) );
}

template<typename A, typename C, int Nvar, int order, int Var>
inline
ADExpr<@STYPE@<ADType<A, Nvar, order, Var>,ADCst<C> > >
@OP@ ( const ADType<A, Nvar, order, Var> &e, C a )
{
    typedef @STYPE@<ADType<A, Nvar, order, Var>,ADCst<C> > expr_t;
    return ADExpr<expr_t>( expr_t ( e ,ADCst<C>( a ) ) );
}

template<class E, class C>
inline
ADExpr<@STYPE@<ADCst<C>,ADExpr<E> > >
@OP@ ( C t, const ADExpr<E> &e )
{
    typedef @STYPE@<ADCst<C>,ADExpr<E> > expr_t;
    return ADExpr<expr_t> ( expr_t ( ADCst<C> ( t ),e ) );
}

template<class E, class C>
inline
ADExpr<@STYPE@<ADExpr<E>,ADCst<C> > >
@OP@ ( const ADExpr<E> &e, C t )
{
    typedef @STYPE@<ADExpr<E>,ADCst<C> > expr_t;
    return ADExpr<expr_t>( expr_t ( e, ADCst<C> ( t ) ) );
}

