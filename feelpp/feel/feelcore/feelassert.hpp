/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-02-19

  Copyright (C) 2010-2022 Feel++ Consortium
  Copyright (C) 2010-2022 University of Strasbourg
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

/*!  \page Macros Macros
  \ingroup Core
  \section assert_macros Assertion Macros

  \subsection assertions Assertions
  It is a good practice to check expected conditions in your program frequently to detect errors as early as possible. 
  Feel++ uses [google-glog](https://github.com/google/glog#check-macros) and its CHECK macros for this purpose.

  The CHECK macro provides the ability to abort the application when a condition is not met, similar to the assert 
  macro defined in the standard C library.

  \par CHECK
  CHECK aborts the application if a condition is not true. Unlike assert, it is *not* controlled by NDEBUG, 
  so the check will be executed regardless of compilation mode. Therefore, fp->Write(x) in the following example is always executed:
  \code
  CHECK(fp->Write(x) == 4) << "Write failed!";
  \endcode 

  \par CHECK_EQ, CHECK_NE, CHECK_LE, CHECK_LT, CHECK_GE, and CHECK_GT
  There are various helper macros for equality/inequality checks - 
  CHECK_EQ, CHECK_NE, CHECK_LE, CHECK_LT, CHECK_GE, and CHECK_GT. 
  They compare two values, and log a FATAL message including the two values when the result is not as expected. 
  The values must have operator<<(ostream, ...) defined.
  You may append to the error message like so:
  \code
  CHECK_NE(1, 2) << ": The world must be ending!";
  \endcode
  We are very careful to ensure that each argument is evaluated exactly once, and that anything which is legal to pass as a function argument is legal here. 
  In particular, the arguments may be temporary expressions which will end up being destroyed at the end of the apparent statement, for example:
  \code {.cpp}
  CHECK_EQ(string("abc")[1], ’b’);
  \endcode
  The compiler reports an error if one of the arguments is a pointer and the other is NULL. 
  To work around this, simply static_cast NULL to the type of the desired pointer.
  \code
  CHECK_EQ(some_ptr, static_cast<SomeType*>(NULL));
  \endcode
  Better yet, use the CHECK_NOTNULL macro:
  \code
  CHECK_NOTNULL(some_ptr);
  some_ptr->DoSomething();
  \endcode


  \par DCHECK
  Some check may be enabled only in Debug mode. To do this use DCHECK which is enabled is NDEBUG is not:
  \code
  DCHECK(fp->Write(x) == 4) << "Write failed!";
  \endcode

*/

// access to smart assertion from feel.hpp
#include <feel/feelcore/smartassert.hpp>

#endif
