/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2005-2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
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

inline constexpr size_type JACOBIAN                 = ( 1<<0 );
inline constexpr size_type KB                       = ( 1<<1 );
inline constexpr size_type KB2                      = ( 1<<2 );
inline constexpr size_type FIRST_DERIVATIVE         = ( 1<<3 );
inline constexpr size_type GRAD                     = ( 1<<4 );
inline constexpr size_type SECOND_DERIVATIVE        = ( 1<<5 );
inline constexpr size_type HESSIAN                  = ( 1<<6 );
inline constexpr size_type LAPLACIAN                = ( 1<<7 );
inline constexpr size_type MEASURE                  = ( 1<<8 );
inline constexpr size_type NORMAL                   = ( 1<<9 );
inline constexpr size_type TANGENT                  = ( 1<<10 );
inline constexpr size_type FIRST_DERIVATIVE_NORMAL  = ( 1<<11 );
inline constexpr size_type POINT                    = ( 1<<12 );
inline constexpr size_type SYMM                     = ( 1<<13 );
inline constexpr size_type UNSYMM                   = ( 1<<14 );
inline constexpr size_type DIV                      = ( 1<<15 );
inline constexpr size_type CURL                     = ( 1<<16 );
inline constexpr size_type INTERPOLANT              = ( 1<<17 );
inline constexpr size_type BASIS_FUNCTION           = ( 1<<18 );
inline constexpr size_type MASS                     = ( 1<<20 );
inline constexpr size_type STIFFNESS                = ( 1<<21 );
inline constexpr size_type NORMAL_COMPONENT         = ( 1<<22 );
inline constexpr size_type LOCAL_BASIS              = ( 1<<23 );
inline constexpr size_type TRACE                    = ( 1<<24 );
inline constexpr size_type DYNAMIC                  = ( 1<<25 );
inline constexpr size_type DYNAMIC_BASIS_FUNCTION   = ( 1<<26 );
inline constexpr size_type GRADDIV                  = ( 1<<27 );
inline constexpr size_type CURLCURL                 = ( 1<<28 );

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
FEELPP_DEFINE_CONTEXT(GRADDIV,graddiv)
FEELPP_DEFINE_CONTEXT(CURLCURL,curlcurl)

/** Mapping derivatives needed by both Hessians and natural second operators. */
template<size_type Context>
inline constexpr bool has_mapping_hessian_v =
    ( Context & ( SECOND_DERIVATIVE | HESSIAN | GRADDIV | CURLCURL ) ) != 0;

inline bool hasMAPPING_HESSIAN( size_type context )
{
    return ( context & ( SECOND_DERIVATIVE | HESSIAN | GRADDIV | CURLCURL ) ) != 0;
}

} // vm
using namespace vm;

} // Feel
#endif /* __FEM_Context_HPP */
