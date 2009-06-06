/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-02-19

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
   \file lifeassert.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-02-19
 */
#ifndef LIFEASSERT_HPP
#define LIFEASSERT_HPP 1

/*!  \page macros_page Life Macros
  \section assert_macros Assertion Macros

  \subsection assertions Assertions

  Life defines a few macros to test pre and post conditions. To
  enable them all you define the preprocessing variable \c
  LIFE_CHECK_ALL . This is done easily by configuring Life with the
  flag \c --enable-debug which will include the option
  \c -DLIFE_CHECK_ALL

  -# #LIFE_ASSERT

  \subsection hints Life C++ Compiler Hints

  -# #INLINE
  -# #LIFE_RESTRICT

  \subsection attribute_macro Life Attribute Macros

  -# #LIFE_EXPORT and #LIFE_NO_EXPORT
  -# #LIFE_PACKED
  -# #LIFE_DEPRECATED
  -# #LIFE_ISLIKELY and #LIFE_ISUNLIKELY
*/

// access to smart assertion from life.hpp
#include <life/lifecore/smartassert.hpp>

#if LIFE_IS_VERSION(0,9,0)

#define ERROR_MSG(A) LIFE_ASSERT( 0 ).error( A );
#define ASSERT0(X,A) LIFE_ASSERT( X ).error( A );
#define ASSERT_PRE0(X,A) LIFE_ASSERT( X ).error( "Precondition Error"  );
#define ASSERT_POS0(X,A) LIFE_ASSERT( X ).error( "Postcondition Error"  );
#define ASSERT_INV0(X,A) LIFE_ASSERT( X ).error( "Invariant Error : "  );
#define ASSERT_BD0(X)    LIFE_ASSERT( X ).error( "Array bounds error" );

#else

# define ERROR_MSG(A)  \
   do { std::cerr << std::endl << std::endl << A << std::endl << std::endl ; ABORT() ; } while (0)



# define ASSERT0(X,A) if ( !(X) ) \
ERROR_MSG(A << std::endl << "Error in file" << __FILE__ << " line " << __LINE__) ;


# define ASSERT_PRE0(X,A) if ( !(X) ) \
ERROR_MSG(A << std::endl << "Precondition Error " << "in file " << __FILE__ \
     << " line " << __LINE__) ;


# define ASSERT_POS0(X,A) if ( !(X) ) \
ERROR_MSG(A << std::endl <<"Postcondition Error " << "in file " << __FILE__ \
     << " line " << __LINE__) ;


# define ASSERT_INV0(X,A)  if ( !(X) ) \
ERROR_MSG(A <<std::endl <<  "Invariant Error " << "in file " << __FILE__  \
   << " line " << __LINE__) ;

# define ASSERT_BD0(X)  if ( !(X) ) \
ERROR_MSG("Array bound error " << "in file " << __FILE__  \
   << " line " << __LINE__) ;

#endif /* 0 */

#if 0
#ifdef  LIFE_CHECK_ALL
#define CHECK_KN
#define TEST_PRE
#define TEST_POS
#define TEST_INV
#define TEST_BOUNDS
#define NOINLINE
#undef  NDEBUG
#endif /* LIFE_CHECK_ALL */

#ifdef NDEBUG
#define ASSERT(X,A)
#else
#define ASSERT(X,A) ASSERT0(X,A)
#endif

#ifdef TEST_PRE
#define ASSERT_PRE(X,A) ASSERT_PRE0(X,A)
#else
#define ASSERT_PRE(X,A)
#endif

#ifdef TEST_POS
#define ASSERT_POS(X,A) ASSERT_POS0(X,A)
#else
#define ASSERT_POS(X,A)
#endif

#ifdef TEST_INV
#define ASSERT_INV(X,A) ASSERT_INV0(X,A)
#else
#define ASSERT_INV(X,A)
#endif

#ifdef TEST_BOUNDS
#define ASSERT_BD(X) ASSERT_BD0(X)
#else
#define ASSERT_BD(X)
#endif

#endif

#endif /* LIFEASSERT_HPP */

