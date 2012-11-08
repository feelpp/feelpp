// CLN internal inlining hints

#ifndef _CL_MAYBE_INLINE_H
#define _CL_MAYBE_INLINE_H

#include "cl_config.h"

/*
 * Selectively inline a function in *some* translation units.
 *
 * The need to inline a function in some places and not in others came from
 * three situations:
 *
 * 1) Some functions, like cl_SF_zerop or cl_FF_zerop, are just a single
 * machine instruction when inlined. Putting their definitions into a public
 * header would expose too many internals of the library, leading violation
 * of abstraction and increased compilation times. Still, it would be nice
 * to use the inline version of these functions in the library itself.
 *
 * 2) Some functions, like cl_{SF,FF,DF,LF}_idecode, are usually only
 * invoked through a dispatcher cl_F_idecode that does nothing but dispatch
 * the call to the right function. Here inlining is used, regardless of
 * the size of the inlined functions, because it removes one function call
 * from the chain of function calls. A compiler cannot know that this
 * caller is the main caller for the 4 inlined functions.
 *
 * 3) Similarly, cl_I_from_NDS would be a bottleneck if not inlined: every
 * creation of a new cl_I goes through this function. A compiler cannot
 * know a priori the bottlenecks.
 *
 * Hence, there is a set of macros which help to create inline and
 * non-inline versions of a function without duplicating the code.
 *
 * Usage:
 *
 * 1. In the public header, declare function as usual:
 *
 * extern cl_bar cl_foo(const cl_baz&);
 *
 * 2. Put the definition into a separate file, say, cl_foo.cc, in the
 * following way:
 *
 * // cl_foo.cc
 *
 * #include "cl_macros.h"
 * #include "whatever/you/need.h"
 *
 * CL_INLINE cl_bar CL_INLINE_DECL(cl_foo)(const cl_baz& x)
 * {
 *   // the actual code goes here
 * }
 *
 * This provides normal (non-inline) version of a function cl_foo.
 *
 * 3. In order to use the inline version, do
 *
 * // cl_blah.cc
 *
 * #include "cl_inline.h"
 * #include "path/to/cl_foo.cc"
 *
 * This will declare and define function cl_foo_inline, which is an inline
 * version of cl_foo.
 *
 * XXX:
 * The name of the inline version *really* has to be different, since ISO C++
 * demands (in 7.1.2.4)
 *
 * "If a function with external linkage is declared inline in one translation
 * unit, it shall be declared inline in all translation units in which it
 * appears; no diagnostic is required."
 *
 * Feel free to implement this functionality in a better *standard-compliant*
 * way. Or submit a DR (defect report) to the standard committee.
 */
#define CL_INLINE
#define CL_INLINE_DECL(fcn) fcn

/*
 * Use these macros to provide inline and non-inline versions of a function
 * which uses an inline version of other function(s). 
 */
#define CL_INLINE2
#define CL_INLINE2_DECL(fcn) fcn

/*
 * Some functions (zerop, signum, etc) just dispatch the call to the
 * appropriate type-specific functions. It would be nice to have these
 * type-specific functions inlined. However, the compiler can not know that,
 * unless one gives it a hint.
 *
 * Usage:
 *
 * const cl_foo CL_FLATTEN cl_bar(const cl_R& x)
 * {
 *   // the actual code
 * }
 *
 * Please note: 
 *
 * 1. This is only a *hint*, it's always up to the compiler to NOT inline
 *    a function.
 * 2. It's ignored if the optimization is switched off.
 */
#if defined(CL_HAVE_ATTRIBUTE_FLATTEN)
#define CL_FLATTEN __attribute__((flatten))
#else
#define CL_FLATTEN
#endif

/*
 * Tell the compiler to inline a function more aggressively, i.e. even if the
 * optimization is switched off.
 *
 * Usage:
 *
 * cl_blah CL_INLINE_HINT cl_foo(const cl_baz& x)
 * {
 *   // the actual code
 * }
 *
 * Notes:
 * 1. This is only a hint, it does NOT guarantee the function will be
 *    actually always inlined.
 * 2. CL_INLINE and CL_INLINE2 macros set this attribute automagically.
 */
#ifdef __GNUC__
#define CL_INLINE_HINT __attribute__((always_inline))
#else
#define CL_INLINE_HINT
#endif

#endif /* _CL_MAYBE_INLINE_H */
