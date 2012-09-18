/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-14

  Copyright (C) 2008-2011 Université Joseph Fourier (Grenoble I)

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
   \file adtype.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
#ifndef __ADType_H
#define __ADType_H 1

#include <set>
#include <boost/type_traits.hpp>
#include <boost/type_traits/is_fundamental.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelopt/adexpr.hpp>
#include <feel/feelopt/adtraits.hpp>

namespace Feel
{
template< typename Expr > class ADExpr;

template<int var_>
struct ADVariable
{
    enum { var = var_ };


};
#if 0
template<int N, template<int> class T>
struct SListGenerator
{
    typedef STYPELIST_1( T<N> ) VL;
    typedef typename ::St::STL::SAppend<VL,typename SListGenerator<N-1,T>::list_type>::Result list_type;
};
template<template<int> class T>
struct SListGenerator<0, T>
{
    typedef T<0> list_type;
};
#endif
/**
  \class ADType
  \brief Automatically Diffentiated Numerical Type

  The idea is based on the work of Nicolas Di Cesare
  who wrote the FAD<Expr> class. However here there
  are some differences: we can into account the depedencies
  and we compute 1st order or 1st and 2nd order derivatives.

  * T is the numerical Type
  * Nvar is the number of variables
  * Order the order of derivation

  So  only up to Order == 2 has been implemented
  here is an example
  <pre>
  #include <ADType.hpp>

  ADType<double,3,2, 0> x(1. );
  ADType<double,3,2, 1> y(2. );
  ADType<double,3,2, 2> z(3. );
  Feel::ADType<double,3,2> __g = x/(y*z);
  std::cout << "g=" << __g << "\n";
  <pre>

  @author Christophe Prud'homme
  @see
*/
template<typename T, int Nvar, int Order = 1, int Var = -1>
class ADType
{
public:

protected:

private:

};

} // Feel
#include <feel/feelopt/adtypeorder0.hpp>
#include <feel/feelopt/adtypeorder1.hpp>
#include <feel/feelopt/adtypeorder2.hpp>

#include <feel/feelopt/adoperators.hpp>
#include <feel/feelopt/adfunctions.hpp>
#include <feel/feelopt/adbinaryfunctions.hpp>

#endif /* __ADType_H */

