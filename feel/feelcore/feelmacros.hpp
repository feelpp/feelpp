/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-01-24

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
   \file feelmacros.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-01-24
 */
#ifndef FEELMACROS_HPP
#define FEELMACROS_HPP 1

/*!  \page Macros Feel Macros
  \section macros Macros

   \subsection hints Feel C++ Compiler Hints

  -# #INLINE
  -# #FEELPP_RESTRICT

  \subsection attribute_macro Feel Attribute Macros

  -# #FEELPP_EXPORT and #FEELPP_NO_EXPORT
  -# #FEELPP_PACKED
  -# #FEELPP_DEPRECATED
  -# #FEELPP_ISLIKELY and #FEELPP_ISUNLIKELY
*/
#ifdef __GNUC__
#define FEELPP_GNUC_AT_LEAST(x,y) ((__GNUC__>=x && __GNUC_MINOR__>=y) || __GNUC__>x)
#else
#define FEELPP_GNUC_AT_LEAST(x,y) 0
#endif

#ifdef __clang__
#define FEELPP_CLANG_AT_LEAST(x,y) ((__clang_major__>=x && __clang_minor__>=y) || __clang_major__>x)
#else
#define FEELPP_CLANG_AT_LEAST(x,y) 0
#endif

/**
   \def FEELPP_CONSTRUCTOR_BEGIN(x)
   Inform that the constructor of the class x has started
 */
#define FEELPP_CONSTRUCTOR_BEGIN(A) DVLOG(3) << "Constructor of " << A << " begins\n";
#define FEELPP_CONSTRUCTOR(A) FEELPP_CONSTRUCTOR_BEGIN(A)
#define CONSTRUCTOR(A) FEELPP_CONSTRUCTOR_BEGIN(A)

/**
   \def FEELPP_CONSTRUCTOR_END(x)
   Inform that the constructor of the class x has ended
 */
#define FEELPP_CONSTRUCTOR_END(A) DVLOG(3) << "Constructor of " << A << " ends\n";

/**
   \def FEELPP_DESTRUCTOR_BEGIN(x)
   Inform that the destructor of the class x has started
 */
#define FEELPP_DESTRUCTOR_BEGIN(A) DVLOG(3) << "Destructor of " << A << " begins\n";
#define FEELPP_DESTRUCTOR(A) FEELPP_DESTRUCTOR_END(A)
#define DESTRUCTOR(A) FEELPP_DESTRUCTOR_BEGIN(A)

/**
   \def FEELPP_DESTRUCTOR_END(x)
   Inform that the destructor of the class x has started
 */
#define FEELPP_DESTRUCTOR_END(A) DVLOG(3) << "Destructor of " << A << " ends\n";


/**
   \def INLINE

   Alias to the C/C++ keyword \c inline
 */
#define INLINE inline

/**
   \def FEELPP_RESTRICT
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
   gcc/g++ a macro has been defined FEELPP_RESTRICT that is equal to
   __restrict__ if the compiler supports it.
 */
#define FEELPP_RESTRICT __restrict__





/**
   \def FEELPP_EXPORT
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
   int FEELPP_NO_EXPORT foo;
   int FEELPP_EXPORT bar;

   extern "C" FEELPP_EXPORT void function(int a);

   class FEELPP_EXPORT SomeClass
   {
     int c;

     // Only for use within this DSO
     FEELPP_NO_EXPORT void privateMethod();

    public:

     Person(int _c) : c(_c) { }
     static void foo(int a);
    };
   \endcode
*/
/**
   \def FEELPP_NO_EXPORT

   Counterpart to #FEELPP_EXPORT.
 */
#if __GNUC__ - 0 > 3 || (__GNUC__ - 0 == 3 && __GNUC_MINOR__ - 0 > 2)
#define FEELPP_EXPORT __attribute__ ((visibility("default")))

#define FEELPP_NO_EXPORT __attribute__ ((visibility("hidden")))
#else
#define FEELPP_EXPORT
#define FEELPP_NO_EXPORT
#endif

/**
   \def FEELPP_PACKED
   The FEELPP_PACKED can be used to hint the compiler that a particular
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
   } FEELPP_PACKED

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
#define FEELPP_PACKED __attribute__((__packed__))
#else
#define FEELPP_PACKED
#endif

/**
   The FEELPP_DEPRECATED macro can be used to trigger compile-time warnings
   with gcc >= 3.2 when deprecated functions are used.

   For non-inline functions, the macro gets inserted at the very end of the
   function declaration, right before the semicolon:

   \code
   DeprecatedConstructor() FEELPP_DEPRECATED;
   void deprecatedFunctionA() FEELPP_DEPRECATED;
   int deprecatedFunctionB() const FEELPP_DEPRECATED;
   \endcode

   Functions which are implemented inline are handled differently: for them,
   the FEELPP_DEPRECATED macro is inserted at the front, right before the return
   type, but after "static" or "virtual":

   \code
   FEELPP_DEPRECATED void deprecatedInlineFunctionA() { .. }
   virtual FEELPP_DEPRECATED int deprecatedInlineFunctionB() { .. }
   static FEELPP_DEPRECATED bool deprecatedInlineFunctionC() { .. }
   \end

   You can also mark whole structs or classes as deprecated, by inserting the
   FEELPP_DEPRECATED macro after the struct/class keyword, but before the
   name of the struct/class:

   \code
   class FEELPP_DEPRECATED DeprecatedClass { };
   struct FEELPP_DEPRECATED DeprecatedStruct { };
   \endcode
*/
#if __GNUC__ - 0 > 3 || (__GNUC__ - 0 == 3 && __GNUC_MINOR__ - 0 >= 2)
# define FEELPP_DEPRECATED __attribute__ ((deprecated))
#else
# define FEELPP_DEPRECATED
#endif

/**
   \def FEELPP_ISLIKELY(x)
   The FEELPP_ISLIKELY macro tags a boolean expression as likely to evaluate to
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
   The FEELPP_ISUNLIKELY macro tags an expression as unlikely evaluating to 'true'.

   \note
   Do NOT use ( !FEELPP_ISLIKELY(foo) ) as an replacement for FEELPP_ISUNLIKELY !

   \code
   if ( FEELPP_ISUNLIKELY( testsomething() ) )
       abort();     // assume its unlikely that the application aborts
   \endcode
*/
/**
   \def FEELPP_ISUNLIKELY(x)
   Counterpart to #FEELPP_ISLIKELY
   The FEELPP_ISUNLIKELY macro tags an expression as unlikely evaluating to 'true'.
 */
#if __GNUC__ - 0 >= 3
# define FEELPP_ISLIKELY( x )    __builtin_expect(!!(x),1)
# define FEELPP_ISUNLIKELY( x )  __builtin_expect(!!(x),0)
#else
# define FEELPP_ISLIKELY( x )   ( x )
# define FEELPP_ISUNLIKELY( x )  ( x )
#endif

/**
 * \def FEELPP_PREFETCH(x)
 * \brief Prefetching

 Another important method of improving performance is through caching
 of necessary data close to the processor. Caching minimizes the
 amount of time it takes to access the data. Most modern processors
 have three classes of memory:

 -# Level 1 cache commonly supports single-cycle access

 -# Level 2 cache supports two-cycle access

 -# System memory supports longer access times

 To to minimize access latency, and thus improve performance, it's
 best to have your data in the closest memory. Performing this task
 manually is called prefetching. GCC supports manual prefetching of
 data through a built-in function called __builtin_prefetch. You use
 this function to pull data into the cache shortly before it's
 needed. As shown below, the __builtin_prefetch function takes three
 arguments:

 -# The address of the data

 -# The rw parameter, which you use to indicate whether the data is
 being pulled in for Read or preparing for a Write operation

 -# The locality parameter, which you use to define whether the data
 should be left in cache or purged after use

 \code
 void  __builtin_prefetch( const void *addr, int rw, int locality );
 \endcode


 Prefetching is used extensively by the Linux kernel. Most often it is
 used through macros and wrapper functions. Listing below is an example of
 a helper function that uses a wrapper over the built-in function
 (from ./linux/include/linux/prefetch.h). The function implements a
 preemptive look-ahead mechanism for streamed operations. Using this
 function can generally result in better performance by minimizing
 cache misses and stalls.

 \code
 #ifndef ARCH_HAS_PREFETCH
 #define prefetch(x) __builtin_prefetch(x)
 #endif

 static inline void prefetch_range(void *addr, size_t len)
 {
 #ifdef ARCH_HAS_PREFETCH
 char *cp;
 char *end = addr + len;

 for (cp = addr; cp < end; cp += PREFETCH_STRIDE)
   prefetch(cp);
 #endif
 }
 \endcode
 */
#if __GNUC__ - 0 >= 3
# define FEELPP_PREFETCH( x, rw, locality )   __builtin_prefetch( (x), rw, locality )
#else
# define FEELPP_PREFETCH( x, rw, locality )
#endif // __GNUC__

/**
 * \def FEELPP_IS_CONSTANT(x)
 * \brief detect at compile if it is a constant

 GCC provides a built-in function that you can use to determine
 whether a value is a constant at compile-time. This is valuable
 information because you can construct expressions that can be
 optimized through constant folding. The __builtin_constant_p function
 is used to test for constants.

 The prototype for __builtin_constant_p is shown below. Note that
 __builtin_constant_p cannot verify all constants, because some are not
 easily proven by GCC.

 \code
 int __builtin_constant_p( exp )
 \endcode


 Linux uses constant detection quite frequently. In the example shown
 in Listing 3 (from ./linux/include/linux/log2.h), constant detection
 is used to optimize the roundup_pow_of_two macro. If the expression
 can be verified as a constant, then a constant expression (which is
 available for optimization) is used. Otherwise, if the expression is
 not a constant, another macro function is called to round up the
 value to a power of two.


 Listing. Constant detection to optimize a macro function
\code
#define roundup_pow_of_two(n)			\
(						\
	__builtin_constant_p(n) ? (		\
		(n == 1) ? 1 :			\
		(1UL << (ilog2((n) - 1) + 1))	\
				   ) :		\
	__roundup_pow_of_two(n)			\
)
\endcode
 */
#if __GNUC__ - 0 >= 3
# define FEELPP_IS_CONSTANT( n ) __builtin_constant_p( n )
#else
# define FEELPP_IS_CONSTANT( n )
#endif // __GNUC__

#define FEELPP_DEBUG_VAR(x) std::cerr << #x << " = " << x << std::endl;

#ifdef NDEBUG
# ifndef FEELPP_NO_DEBUG
#  define FEELPP_NO_DEBUG
# endif
#endif

// FEELPP_ALWAYS_INLINE_ATTRIB should be use in the declaration of function
// which should be inlined even in debug mode.
#if FEELPP_GNUC_AT_LEAST(4,0)
#define FEELPP_ALWAYS_INLINE_ATTRIB __attribute__((always_inline))
#else
#define FEELPP_ALWAYS_INLINE_ATTRIB
#endif

// FEELPP_FORCE_INLINE means "inline as much as possible"
#if (defined _MSC_VER) || (defined __intel_compiler)
#define FEELPP_STRONG_INLINE __forceinline
#else
#define FEELPP_STRONG_INLINE inline
#endif

#if (defined __GNUC__)
#define FEELPP_DONT_INLINE __attribute__((noinline))
#elif (defined _MSC_VER)
#define FEELPP_DONT_INLINE __declspec(noinline)
#else
#define FEELPP_DONT_INLINE
#endif

#endif /* FEELMACROS_HPP */
