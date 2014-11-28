/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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
const size_type MASS                     = ( 1<<20 );
const size_type STIFFNESS                = ( 1<<21 );



typedef mpl::vector_c<size_type,
                      JACOBIAN, KB, KB2, FIRST_DERIVATIVE, GRAD, SECOND_DERIVATIVE, HESSIAN, LAPLACIAN,
                      MEASURE, NORMAL, TANGENT, FIRST_DERIVATIVE_NORMAL, POINT,
                      SYMM, UNSYMM,
                      DIV,CURL,
                      MASS, STIFFNESS> contexts;

template<size_type Context>
struct has_jacobian
{
    static const bool value = has_value<Context, JACOBIAN>::value;
};
template<size_type Context>
struct has_kb
{
    static const bool value = has_value<Context, KB>::value;
};
template<size_type Context>
struct has_kb2
{
    static const bool value = has_value<Context, KB2>::value;
};
template<size_type Context>
struct has_first_derivative
{
    static const bool value = has_value<Context, FIRST_DERIVATIVE>::value;
};
template<size_type Context>
struct has_grad
{
    static const bool value = has_value<Context, GRAD>::value;
};
template<size_type Context>
struct has_second_derivative
{
    static const bool value = has_value<Context, SECOND_DERIVATIVE>::value;
};
template<size_type Context>
struct has_hessian
{
    static const bool value = has_value<Context, HESSIAN>::value;
};
template<size_type Context>
struct has_laplacian
{
    static const bool value = has_value<Context, LAPLACIAN>::value;
};

template<size_type Context>
struct has_normal
{
    static const bool value = has_value<Context, NORMAL>::value;
};
template<size_type Context>
struct has_tangent
{
    static const bool value = has_value<Context, TANGENT>::value;
};

template<size_type Context>
struct has_first_derivative_normal
{
    static const bool value = has_value<Context, FIRST_DERIVATIVE_NORMAL>::value;
};

template<size_type Context>
struct has_point
{
    static const bool value = has_value<Context, POINT>::value;
};
template<size_type Context>
struct has_symm
{
    static const bool value = has_value<Context, SYMM>::value;
};
template<size_type Context>
struct has_unsymm
{
    static const bool value = has_value<Context, UNSYMM>::value;
};
template<size_type Context>
struct has_div
{
    static const bool value = has_value<Context, DIV>::value;
};
template<size_type Context>
struct has_curl
{
    static const bool value = has_value<Context, CURL>::value;
};

template<size_type Context>
struct has_mass
{
    static const bool value = has_value<Context, MASS>::value;
};
template<size_type Context>
struct has_stifness
{
    static const bool value = has_value<Context, STIFFNESS>::value;
};

} // vm
} // Feel
#endif /* __FEM_Context_HPP */
