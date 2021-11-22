/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-07-19

  Copyright (C) 2005,2006 EPFL

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
   \file context.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-07-19
 */
#ifndef __FEM_Context_HPP
#define __FEM_Context_HPP 1


#include <boost/preprocessor/cat.hpp>
#include <feel/feelcore/context.hpp>
#include <boost/mpl/vector_c.hpp>

namespace Feel
{
namespace vm
{

const size_type JACOBIAN                 = ( 1<<0 );
const size_type KB                       = ( 1<<1 );
const size_type KB2                      = ( 1<<2 );
const size_type FIRST_DERIVATIVE         = ( 1<<3 );
const size_type GRAD                     = ( 1<<4 );
const size_type SECOND_DERIVATIVE        = ( 1<<5 );
const size_type HESSIAN                  = ( 1<<6 );
const size_type LAPLACIAN                = ( 1<<7 );
const size_type MEASURE                  = ( 1<<8 );
const size_type NORMAL                   = ( 1<<9 );
const size_type TANGENT                  = ( 1<<10 );
const size_type FIRST_DERIVATIVE_NORMAL  = ( 1<<11 );
const size_type POINT                    = ( 1<<12 );
const size_type SYMM                     = ( 1<<13 );
const size_type UNSYMM                   = ( 1<<14 );
const size_type DIV                      = ( 1<<15 );
const size_type CURL                     = ( 1<<16 );
const size_type INTERPOLANT              = ( 1<<17 );
const size_type BASIS_FUNCTION           = ( 1<<18 );
const size_type MASS                     = ( 1<<20 );
const size_type STIFFNESS                = ( 1<<21 );
const size_type NORMAL_COMPONENT         = ( 1<<22 );
const size_type LOCAL_BASIS              = ( 1<<23 );
const size_type TRACE                    = ( 1<<24 );
const size_type DYNAMIC                  = ( 1<<25 );
const size_type DYNAMIC_BASIS_FUNCTION   = ( 1<<26 );

#define FEELPP_DEFINE_CONTEXT(ctx_v,ctx)                                \
    template<size_type Context>                                         \
    using BOOST_PP_CAT(has_,ctx) = has_value<Context, ctx_v>;           \
    template<size_type Context>                                         \
    using BOOST_PP_CAT(BOOST_PP_CAT(has_,ctx),_t) = has_value<Context, ctx_v>; \
    template<size_type Context>                                         \
    constexpr bool BOOST_PP_CAT(has_,BOOST_PP_CAT(ctx,_v)) = has_value<Context, ctx_v>::value; \
    inline bool BOOST_PP_CAT(has,ctx_v)( size_type c ) { return hasValue<ctx_v>( c ); }

FEELPP_DEFINE_CONTEXT(JACOBIAN,jacobian)
FEELPP_DEFINE_CONTEXT(KB,kb)
FEELPP_DEFINE_CONTEXT(KB2,kb2)
FEELPP_DEFINE_CONTEXT(FIRST_DERIVATIVE,first_derivative)
FEELPP_DEFINE_CONTEXT(GRAD,grad)
FEELPP_DEFINE_CONTEXT(SECOND_DERIVATIVE,second_derivative)
FEELPP_DEFINE_CONTEXT(HESSIAN,hessian)
FEELPP_DEFINE_CONTEXT(LAPLACIAN,laplacian)
FEELPP_DEFINE_CONTEXT(MEASURE,measure)
FEELPP_DEFINE_CONTEXT(NORMAL,normal)
FEELPP_DEFINE_CONTEXT(TANGENT,tangent)
FEELPP_DEFINE_CONTEXT(FIRST_DERIVATIVE_NORMAL,first_derivative_normal)
FEELPP_DEFINE_CONTEXT(POINT,point)
FEELPP_DEFINE_CONTEXT(SYMM,symm)
FEELPP_DEFINE_CONTEXT(UNSYMM,unsymm)
FEELPP_DEFINE_CONTEXT(DIV,div)
FEELPP_DEFINE_CONTEXT(CURL,curl)
FEELPP_DEFINE_CONTEXT(INTERPOLANT,interpolant)
FEELPP_DEFINE_CONTEXT(BASIS_FUNCTION,basis_function)
FEELPP_DEFINE_CONTEXT(MASS,mass)
FEELPP_DEFINE_CONTEXT(STIFFNESS,stiffness)
FEELPP_DEFINE_CONTEXT(NORMAL_COMPONENT,normal_component)
FEELPP_DEFINE_CONTEXT(LOCAL_BASIS,local_basis)
FEELPP_DEFINE_CONTEXT(TRACE,trace)
FEELPP_DEFINE_CONTEXT(DYNAMIC,dynamic)
FEELPP_DEFINE_CONTEXT(DYNAMIC_BASIS_FUNCTION,dynamic_basis_function)

} // vm
using namespace vm;

} // Feel
#endif /* __FEM_Context_HPP */
