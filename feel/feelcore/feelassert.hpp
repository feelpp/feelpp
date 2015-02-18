/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-02-19

  Copyright (C) 2009 Université de Grenoble 1
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
   \file feelassert.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-02-19
 */
#ifndef FEELASSERT_HPP
#define FEELASSERT_HPP 1

/*!  \page Macros Feel Macros
  \section assert_macros Assertion Macros

  \subsection assertions Assertions

  Feel defines a few macros to test pre and post conditions. To
  enable them all you define the preprocessing variable \c
  FEELPP_CHECK_ALL . This is done easily by configuring Feel with the
  flag \c --enable-debug which will include the option
  \c -DFEELPP_CHECK_ALL

  -# #FEELPP_ASSERT

  \subsection hints Feel C++ Compiler Hints

  -# #INLINE
  -# #FEELPP_RESTRICT

  \subsection attribute_macro Feel Attribute Macros

  -# #FEELPP_EXPORT and #FEELPP_NO_EXPORT
  -# #FEELPP_PACKED
  -# #FEELPP_DEPRECATED
  -# #FEELPP_ISLIKELY and #FEELPP_ISUNLIKELY
*/

// access to smart assertion from feel.hpp
#include <feel/feelcore/smartassert.hpp>

#endif
