/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-05-30

  Copyright (C) 2007 Universite Joseph Fourier (Grenoble I)

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
   \file form.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-05-30
 */
#ifndef FEELPP_VF_FORM_H
#define FEELPP_VF_FORM_H

#include <feel/feelcore/parameter.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/vector.hpp>
#include <feel/feelalg/matrixsparse.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feelvf/bilinearform.hpp>
#include <feel/feelvf/linearform.hpp>
#include <feel/feelvf/integrate.hpp>

namespace Feel
{
using namespace std::string_literals;
template<typename T, typename SizeT> class Vector;

//
// free form functions
//
template<typename X1, typename X2>
inline
vf::detail::BilinearForm<X1, X2>
form( std::string name,
      std::shared_ptr<X1> const& __X1,
      std::shared_ptr<X2> const& __X2,
      std::shared_ptr<MatrixSparse<double> > __M,
      size_type rowstart = 0,
      size_type colstart = 0,
      bool init = false,
      bool do_threshold = false,
      typename X1::value_type threshold = type_traits<double>::epsilon(),
      size_type pattern  = Pattern::COUPLED )
{
    return vf::detail::BilinearForm<X1, X2>( name, __X1, __X2, __M, rowstart, colstart, init, do_threshold, threshold, pattern );
}

template<typename X1, typename RepType>
inline
vf::detail::LinearForm<X1, RepType, RepType>
form( std::string name,
      std::shared_ptr<X1> const& __X1,
      std::shared_ptr<RepType> __M,
      size_type rowstart = 0,
      bool init = false,
      bool do_threshold = false,
      typename X1::value_type threshold = type_traits<typename RepType::value_type>::epsilon() )
{
    return vf::detail::LinearForm<X1, RepType, RepType>( name, __X1, __M, rowstart, init, do_threshold, threshold );
}


/**
 * @addtogroup FreeFunction
 * @{
 */

template <typename ... Ts>
auto form1( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && test = args.get(_test );
    auto && backend = args.get_else_invocable(_backend, [&test](){ return Feel::backend(_worldcomm=test->worldCommPtr()); } );
    auto && vector = args.get_else_invocable(_vector, [&test,&backend](){ return backend->newVector( _test=test ); } );
    bool init = args.get_else(_init, false );
    bool do_threshold = args.get_else(_do_threshold, false );
    double threshold = args.get_else(_threshold, type_traits<double>::epsilon() );
    size_type rowstart = args.get_else(_rowstart, 0 );
    std::string const& name = args.get_else(_name, "linearform.f" );

    return form( name, test, vector, rowstart, init, do_threshold, threshold );
}


/**
 * @}
 */

template <typename ... Ts>
auto form2( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && test = args.get(_test );
    auto && trial = args.get(_trial );
    bool init = args.get_else(_init, false );
    size_type properties = args.get_else(_properties, NON_HERMITIAN );
    size_type pattern = args.get_else(_pattern,Pattern::COUPLED );
    auto && backend = args.get_else_invocable(_backend, [&test](){ return Feel::backend(_worldcomm=test->worldCommPtr() ); } );
    std::shared_ptr<MatrixSparse<double>> matrix = args.get_else_invocable(_matrix, [&backend,&test,&trial,&pattern,&properties](){ return backend->newMatrix( _test=test, _trial=trial, _pattern=pattern, _properties=properties ); } );
    size_type rowstart = args.get_else(_rowstart, 0 );
    size_type colstart = args.get_else(_colstart, 0 );
    std::string const& name = args.get_else(_name, "bilinearform.a" );

    bool do_threshold = false;
    double threshold = 1e-16;
    return form( name, test, trial, matrix, rowstart, colstart, init, do_threshold, threshold, pattern );
}



} // Feel
#endif /* __Form_H */
