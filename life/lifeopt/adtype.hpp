/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-02-14

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-02-14
 */
#ifndef __ADType_H
#define __ADType_H 1

#include <set>

#include <life/lifeopt/adexpr.hpp>
#include <life/lifeopt/adtraits.hpp>

namespace Life
{
template< typename Expr > class ADExpr;

template<int var_>
struct ADVariable
{
    enum { var = var_ };


};

template<int N, template<int> class T>
struct SListGenerator
{
    typedef STYPELIST_1(T<N>) VL;
    typedef typename ::St::STL::SAppend<VL,typename SListGenerator<N-1,T>::list_type>::Result list_type;
};
template<template<int> class T>
struct SListGenerator<0, T>
{
    typedef T<0> list_type;
};

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

  So far only Order = 2 has been implemented
  here is an example
  <pre>
  #include <ADType.hpp>

  ADType<double,3,2, Var> x(1.,0 );
  ADType<double,3,2, Var> y(2.,1 );
  ADType<double,3,2, Var> z(3.,2 );
  St::Array::ADType<double,3,2, Var> __g = x/(y*z);
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

}
} // Life
#include <life/lifeopt/adtypeorder0.hpp>
#include <life/lifeopt/adtypeorder1.hpp>
#include <life/lifeopt/adtypeorder2.hpp>

#include <life/lifeopt/adoperators.hpp>
#include <life/lifeopt/adfunctions.hpp>
#include <life/lifeopt/adbinaryfunctions.hpp>

#endif /* __ADType_H */

