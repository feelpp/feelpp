/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-01-24

  Copyright (C) 2005,2006 EPFL

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
   \file lifemacros.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-01-24
 */
#ifndef LIFEMACROS_HPP
#define LIFEMACROS_HPP 1

/*!  \page macros_page Life Macros
  \section macros Macros

   \subsection hints Life C++ Compiler Hints

  -# #INLINE
  -# #LIFE_RESTRICT

  \subsection attribute_macro Life Attribute Macros

  -# #LIFE_EXPORT and #LIFE_NO_EXPORT
  -# #LIFE_PACKED
  -# #LIFE_DEPRECATED
  -# #LIFE_ISLIKELY and #LIFE_ISUNLIKELY
*/
/**
   \def LIFE_CONSTRUCTOR_BEGIN(Area,x)
   Inform that the constructor of the class x has started
 */
#define LIFE_CONSTRUCTOR_BEGIN(Area, A) Debug( Area ) << "Constructor of " << A << " begins\n";
#define LIFE_CONSTRUCTOR(Area,A) LIFE_CONSTRUCTOR_BEGIN(Area,A)
#define CONSTRUCTOR(A) LIFE_CONSTRUCTOR_BEGIN(20000,A)

/**
   \def LIFE_CONSTRUCTOR_END(Area,x)
   Inform that the constructor of the class x has ended
 */
#define LIFE_CONSTRUCTOR_END(Area,A) Debug( Area ) << "Constructor of " << A << " ends\n";

/**
   \def LIFE_DESTRUCTOR_BEGIN(Area,x)
   Inform that the destructor of the class x has started
 */
#define LIFE_DESTRUCTOR_BEGIN(Area,A) Debug( Area ) << "Destructor of " << A << " begins\n";
#define LIFE_DESTRUCTOR(Area,A) LIFE_DESTRUCTOR_END(Area,A)
#define DESTRUCTOR(A) LIFE_DESTRUCTOR_BEGIN(20000,A)

/**
   \def LIFE_DESTRUCTOR_END(Area,x)
   Inform that the destructor of the class x has started
 */
#define LIFE_DESTRUCTOR_END(Area,A) Debug( Area ) << "Destructor of " << A << " ends\n";


/**
   \def INLINE

   Alias to the C/C++ keyword \c inline
 */
#define INLINE inline

/**
   \def LIFE_RESTRICT
   \brief C99 feature of restricted(not aliased) pointers and references

   As with gcc, g++ understands the C99 feature of restricted
   pointers, specified with the \c __restrict__, or __restrict type
   qualifier. Because you cannot compile C++ by specifying the
   -std=c99 language flag, restrict is not a keyword in C++.


   In addition to allowing restricted pointers, you can specify
   restricted references, which indicate that the reference is not
   aliased in the local context.

   \code
   void fn (int *__restrict__ rptr, int &__restrict__ rref)
   {
   ...
   }
   \endcode

   In the body of \c fn, \c rptr points to an unaliased integer and \c rref
   refers to a (different) unaliased integer.


   You may also specify whether a member function's this pointer is
   unaliased by using \c __restrict__ as a member function qualifier.

   \code
   void T::fn () __restrict__
   {
   ...
   }
   \endcode

   Within the body of T::fn, this will have the effective definition
   <tt>T* __restrict__</tt> const this. Notice that the interpretation of a
   \c __restrict__ member function qualifier is different to that of
   const or volatile qualifier, in that it is applied to the pointer
   rather than the object. This is consistent with other compilers
   which implement restricted pointers.

   As with all outermost parameter qualifiers, \c __restrict__ is ignored
   in function definition matching. This means you only need to
   specify \c __restrict__ in a function definition, rather than in a
   function prototype as well.

   In order to ensure that the code is portable to other compiler than
   gcc/g++ a macro has been defined LIFE_RESTRICT that is equal to
   __restrict__ if the compiler supports it.
 */
#define LIFE_RESTRICT __restrict__





/**
   \def LIFE_EXPORT
   \brief Load time improvements for DSO libraries

   Here are a few explanations why this is useful.  For more info
   checkout http://www.nedprod.com/programs/gccvisibility.html

   -# It very substantially improves load times of your DSO (Dynamic
   Shared Object) For example, the TnFOX Boost.Python bindings library
   now loads in eight seconds rather than over six minutes!

   -# It lets the optimiser produce better code PLT indirections (when a
   function call or variable access must be looked up via the Global
   Offset Table such as in PIC code) can be completely avoided, thus
   substantially avoiding pipeline stalls on modern processors and thus
   much faster code. Furthermore when most of the symbols are bound
   locally, they can be safely elided (removed) completely through the
   entire DSO. This gives greater latitude especially to the inliner
   which no longer needs to keep an entry point around "just in case".

   -# It reduces the size of your DSO by 5-20% ELF's exported symbol
   table format is quite a space hog, giving the complete mangled symbol
   name which with heavy template usage can average around 1000
   bytes. C++ templates spew out a huge amount of symbols and a typical
   C++ library can easily surpass 30,000 symbols which is around 5-6Mb!
   Therefore if you cut out the 60-80% of unnecessary symbols, your DSO
   can be megabytes smaller!

   -# Much lower chance of symbol collision The old woe of two libraries
   internally using the same symbol for different things is finally
   behind us with this patch. Hallelujah!

   here is an example on how to use them
   \code
   int LIFE_NO_EXPORT foo;
   int LIFE_EXPORT bar;

   extern "C" LIFE_EXPORT void function(int a);

   class LIFE_EXPORT SomeClass
   {
     int c;

     // Only for use within this DSO
     LIFE_NO_EXPORT void privateMethod();

    public:

     Person(int _c) : c(_c) { }
     static void foo(int a);
    };
   \endcode
*/
/**
   \def LIFE_NO_EXPORT

   Counterpart to #LIFE_EXPORT.
 */
#if __GNUC__ - 0 > 3 || (__GNUC__ - 0 == 3 && __GNUC_MINOR__ - 0 > 2)
#define LIFE_EXPORT __attribute__ ((visibility("default")))

#define LIFE_NO_EXPORT __attribute__ ((visibility("hidden")))
#else
#define LIFE_EXPORT
#define LIFE_NO_EXPORT
#endif

/**
   \def LIFE_PACKED
   The LIFE_PACKED can be used to hint the compiler that a particular
   structure or class should not contain unnecessary paddings.

   Here is an explanation from http://sig9.com/articles/gcc-packed-structures

   GCC allows you to specify attributes of variables and structures
   using the keyword \c __attribute__, the syntax of which is
   \c __attribute__((attribute list)). One such attribute is \c __packed__
   which specifies that

   a variable or structure field should have the smallest possible
   alignment--one byte for a variable, and one bit for a field, unless
   you specify a larger value with the aligned attribute.


   which means that GCC will not add any of the zero's for padding (for
   memory alignement) and make variables or fields immediately next to
   each other. For example, here are some things I tried out -- I created
   a C source file - \c test.c

   \code
   struct test_t {
   int  a;
   char b;
   int  c;
   } ;

   struct test_t test = { 10, 20, 30};
   \endcode

   And compiled it with the -S option (ie to generate the assembly
   equivalent of the code generated).

   \code
      .file "t.cpp"
      .globl test
        .data
        .align 4
        .type test, @object
        .size test, 12
      <b>test:
      .long 10
      .byte 20
      .zero 3
      .long 30</b>
      .section .note.GNU-stack,"",@progbits
      .ident   "GCC: (GNU) 3.3.5 (Debian 1:3.3.5-6)"
   \endcode

   Notice the emphasized code. You can see that the structure "test"
   is being declared. First the field "a" (int) as .long 10 followed
   by "b" (char) as .byte 20. To keep the fields' word alignment,
   notice that GCC has added 3 zero bytes (.zero 3) before field "c"
   (int) which is declared as .long 30. This makes the effective
   sizeof struct test_t as 12 instead of the expected 9. Then I tried
   with the __packed__ attribute -

   \code
   struct test_t {
   int  a;
   char b;
   int  c;
   } LIFE_PACKED

   struct test_t test = { 10, 20, 30};
   \endcode

   and the "-S" output I got after compiling was

   \code
   .file "t.cpp"
   .globl test
     .data
     .type test, @object
     .size test, 9
   test:
     .long 10
     .byte 20
     .long 30
   .section .note.GNU-stack,"",@progbits
   .ident   "GCC: (GNU) 3.3.5 (Debian 1:3.3.5-6)"
   \endcode

   in which the zeros are missing making the sizeof structure test_t =
   9. Always remember that memory alignment is *good* even if it
   compromises space, so think twice before using this attribute. It
   is generally useful when you want to assign a structure to a block
   of memory and manipulate it through the fields of a structure.
 */
#ifdef __GNUC__
#define LIFE_PACKED __attribute__((__packed__))
#else
#define LIFE_PACKED
#endif

/**
   The LIFE_DEPRECATED macro can be used to trigger compile-time warnings
   with gcc >= 3.2 when deprecated functions are used.

   For non-inline functions, the macro gets inserted at the very end of the
   function declaration, right before the semicolon:

   \code
   DeprecatedConstructor() LIFE_DEPRECATED;
   void deprecatedFunctionA() LIFE_DEPRECATED;
   int deprecatedFunctionB() const LIFE_DEPRECATED;
   \endcode

   Functions which are implemented inline are handled differently: for them,
   the LIFE_DEPRECATED macro is inserted at the front, right before the return
   type, but after "static" or "virtual":

   \code
   LIFE_DEPRECATED void deprecatedInlineFunctionA() { .. }
   virtual LIFE_DEPRECATED int deprecatedInlineFunctionB() { .. }
   static LIFE_DEPRECATED bool deprecatedInlineFunctionC() { .. }
   \end

   You can also mark whole structs or classes as deprecated, by inserting the
   LIFE_DEPRECATED macro after the struct/class keyword, but before the
   name of the struct/class:

   \code
   class LIFE_DEPRECATED DeprecatedClass { };
   struct LIFE_DEPRECATED DeprecatedStruct { };
   \endcode
*/
#if __GNUC__ - 0 > 3 || (__GNUC__ - 0 == 3 && __GNUC_MINOR__ - 0 >= 2)
# define LIFE_DEPRECATED __attribute__ ((deprecated))
#else
# define LIFE_DEPRECATED
#endif

/**
   \def LIFE_ISLIKELY(x)
   The LIFE_ISLIKELY macro tags a boolean expression as likely to evaluate to
   'true'. When used in an if ( ) statement, it gives a hint to the compiler
   that the following codeblock is likely to get executed. Providing this
   information helps the compiler to optimize the code for better performance.
   Using the macro has an insignificant code size or runtime memory footprint impact.
   The code semantics is not affected.

   \note
   Providing wrong information ( like marking a condition that almost never
   passes as 'likely' ) will cause a significant runtime slowdown. Therefore only
   use it for cases where you can be sure about the odds of the expression to pass
   in all cases ( independent from e.g. user configuration ).

   \par
   The LIFE_ISUNLIKELY macro tags an expression as unlikely evaluating to 'true'.

   \note
   Do NOT use ( !LIFE_ISLIKELY(foo) ) as an replacement for LIFE_ISUNLIKELY !

   \code
   if ( LIFE_ISUNLIKELY( testsomething() ) )
       abort();     // assume its unlikely that the application aborts
   \endcode
*/
/**
   \def LIFE_ISUNLIKELY(x)
   Counterpart to #LIFE_ISLIKELY
   The LIFE_ISUNLIKELY macro tags an expression as unlikely evaluating to 'true'.
 */
#if __GNUC__ - 0 >= 3
# define LIFE_ISLIKELY( x )    __builtin_expect(!!(x),1)
# define LIFE_ISUNLIKELY( x )  __builtin_expect(!!(x),0)
#else
# define LIFE_ISLIKELY( x )   ( x )
# define LIFE_ISUNLIKELY( x )  ( x )
#endif

#endif /* LIFEMACROS_HPP */
