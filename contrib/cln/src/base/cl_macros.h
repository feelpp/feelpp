// CLN internal macros

#ifndef _CL_MACROS_H
#define _CL_MACROS_H

#include "cln/types.h"
#include "cln/exception.h"

// Concatenation of macroexpanded tokens.
// Example:
//   #undef x
//   #define y 16
//   CONCAT(x,y)        ==>  'x16' (not 'xy' !)
  #define CONCAT_(xxx,yyy)  xxx##yyy
  #define CONCAT3_(aaa,bbb,ccc)  aaa##bbb##ccc
  #define CONCAT4_(aaa,bbb,ccc,ddd)  aaa##bbb##ccc##ddd
  #define CONCAT5_(aaa,bbb,ccc,ddd,eee)  aaa##bbb##ccc##ddd##eee
  #define CONCAT6_(aaa,bbb,ccc,ddd,eee,fff)  aaa##bbb##ccc##ddd##eee##fff
  #define CONCAT7_(aaa,bbb,ccc,ddd,eee,fff,ggg)  aaa##bbb##ccc##ddd##eee##fff##ggg
  #define CONCAT(xxx,yyy)  CONCAT_(xxx,yyy)
  #define CONCAT3(aaa,bbb,ccc)  CONCAT3_(aaa,bbb,ccc)
  #define CONCAT4(aaa,bbb,ccc,ddd)  CONCAT4_(aaa,bbb,ccc,ddd)
  #define CONCAT5(aaa,bbb,ccc,ddd,eee)  CONCAT5_(aaa,bbb,ccc,ddd,eee)
  #define CONCAT6(aaa,bbb,ccc,ddd,eee,fff)  CONCAT6_(aaa,bbb,ccc,ddd,eee,fff)
  #define CONCAT7(aaa,bbb,ccc,ddd,eee,fff,ggg)  CONCAT7_(aaa,bbb,ccc,ddd,eee,fff,ggg)

// Convert tokens to strings.
// STRING(token)  ==>  "token"
  #define STRING(token) #token
  #define STRINGIFY(token) STRING(token)

// Declare functions that don't return.
// nonreturning_function(extern,exit,(void)); == extern void exit (void);
  #ifdef __GNUC__
    #if (__GNUC__ >= 3) || ((__GNUC__ == 2) && (__GNUC_MINOR__ >= 9))
      #define nonreturning_function(storclass,funname,arguments)  \
        storclass void funname arguments __attribute__((__noreturn__))
    #else
      #define nonreturning_function(storclass,funname,arguments)  \
        typedef void CONCAT3(funname,_function_,__LINE__) arguments; \
        storclass __volatile__ CONCAT3(funname,_function_,__LINE__) funname
    #endif
  #else
    #define nonreturning_function(storclass,funname,arguments)  \
      storclass void funname arguments
  #endif

// Declaration of variables.
  #define var

// `if' with more than one clause:
// if (cond1) ... {elif (condi) ...} [else ...]
  #define elif  else if

// Endless loop, leave with  break;  or return...;
  #define loop  while (1)

// Reversed end condition.
// Allows   until (expression) statement
// and      do statement until (expression);
  #define until(expression)  while(!(expression))

// Boolean values.
  #define FALSE  0
  #define TRUE   1

// Ignore a value (instead of assigning it to a variable).
// unused ...
  #if defined(__GNUC__) || defined(__KCC) // avoid a gcc warning "statement with no effect"
    #define unused  (void)
  #else
    #define unused
  #endif

// Denotes a point where control flow can never arrive.
// NOTREACHED
  #define NOTREACHED  throw notreached_exception(__FILE__,__LINE__);

// Check an arithmetic expression.
// ASSERT(expr)
  #define ASSERT(expr)  { if (!(expr)) { NOTREACHED } }

// alloca()
  #if defined(__GNUC__) && !defined(__riscos) && !defined(__convex__)
    #undef alloca
    #define alloca  __builtin_alloca
  #elif defined(_MSC_VER)
    #include <malloc.h>
    #define alloca  _alloca
  #elif defined(HAVE_ALLOCA_H) || defined(__riscos)
    #include <alloca.h>
    #ifndef alloca // Sometimes `alloca' is defined as a macro...
      #if defined(__osf__)
        extern "C" char* alloca (int size);
      #else
        extern "C" void* alloca (size_t size);
      #endif
    #endif
  #elif defined(_AIX)
    #pragma alloca // AIX requires this to be the first thing in the file.
  #elif defined(WATCOM)
    #include <malloc.h> // defines `alloca' as a macro
  #elif !defined(NO_ALLOCA)
    extern "C" void* alloca (size_t size);
  #endif

// NULL pointer.
  #undef NULL
  #define NULL  0

// Bit number n (0<=n<32 or 0<=n<64)
  #ifdef HAVE_FAST_LONGLONG
    #define bit(n)  (1LL<<(n))
  #else
    #define bit(n)  (1L<<(n))
  #endif
// Bit number n (0<n<=32) mod 2^32
  #ifdef HAVE_FAST_LONGLONG
    #define bitm(n)  (2LL<<((n)-1))
  #else
    #define bitm(n)  (2L<<((n)-1))
  #endif
// Test bit n in x, n constant, x a cl_uint:
  #if !(defined(__sparc__) || defined(__sparc64__))
    #define bit_test(x,n)  ((x) & bit(n))
  #else
    // On Sparcs long constants are slower than shifts.
    #if !defined(__GNUC__)
      #define bit_test(x,n)  \
        ((n)<12 ? ((x) & bit(n)) : ((sint32)((uint32)(x) << (31-(n))) < 0))
    #else // gcc optimizes boolean expressions better this way:
      #define bit_test(x,n)  \
        (   ( ((n)<12) && ((x) & bit(n)) )                           \
         || ( ((n)>=12) && ((sint32)((uint32)(x) << (31-(n))) < 0) ) \
        )
    #endif
  #endif
// minus bit number n (0<=n<32 or 0<=n<64)
  #ifdef HAVE_FAST_LONGLONG
    #define minus_bit(n)  (-1LL<<(n))
  #else
    #define minus_bit(n)  (-1L<<(n))
  #endif
// minus bit number n (0<n<=32) mod 2^32
  #ifdef HAVE_FAST_LONGLONG
    #define minus_bitm(n)  (-2LL<<((n)-1))
  #else
    #define minus_bitm(n)  (-2L<<((n)-1))
  #endif

// Return 2^n, n a constant expression.
// Same as bit(n), but undefined if n<0 or n>={long_}long_bitsize.
  #if defined(HAVE_FAST_LONGLONG) || defined(intQsize)
    #define bitc(n)  (1ULL << (((n) >= 0 && (n) < long_long_bitsize) ? (n) : 0))
  #else
    #define bitc(n)  (1UL << (((n) >= 0 && (n) < long_bitsize) ? (n) : 0))
  #endif

// floor(a,b) for a>=0, b>0 returns floor(a/b).
// b should be a constant expression.
  #define floor(a_from_floor,b_from_floor)  ((a_from_floor) / (b_from_floor))
// Save the macro in case we need to include <cmath>.
  #define cln_floor(a_from_floor,b_from_floor)  ((a_from_floor) / (b_from_floor))

// ceiling(a,b) for a>=0, b>0 returns ceiling(a/b) = floor((a+b-1)/b).
// b should be a constant expression.
  #define ceiling(a_from_ceiling,b_from_ceiling)  \
    (((a_from_ceiling) + (b_from_ceiling) - 1) / (b_from_ceiling))

// round_down(a,b) decreases a>=0 such that it becomes divisible by b>0.
// b should be a constant expression.
  #define round_down(a_from_round,b_from_round)  \
    (floor(a_from_round,b_from_round)*(b_from_round))

// round_up(a,b) increases a>=0 such that it becomes divisible by b>0.
// b should be a constant expression.
  #define round_up(a_from_round,b_from_round)  \
    (ceiling(a_from_round,b_from_round)*(b_from_round))

// We never call malloc(0), so no need to handle it.
  #define __MALLOC_0_RETURNS_NULL

// Loop which executes a statement a given number of times.
// dotimesC(countvar,count,statement);
// countvar must be of type `uintC'. It is modified!
  #define dotimesC(countvar_from_dotimesC,count_from_dotimesC,statement_from_dotimesC)  \
    { countvar_from_dotimesC = (count_from_dotimesC);         \
      until (countvar_from_dotimesC==0)                       \
        {statement_from_dotimesC; countvar_from_dotimesC--; } \
    }
  #define dotimespC(countvar_from_dotimespC,count_from_dotimespC,statement_from_dotimespC)  \
    { countvar_from_dotimespC = (count_from_dotimespC);                   \
      do {statement_from_dotimespC} until (--countvar_from_dotimespC==0); \
    }

// doconsttimes(count,statement);
// führt statement count mal aus (count mal der Code!),
// wobei count eine constant-expression >=0, <=8 ist.
  #define doconsttimes(count_from_doconsttimes,statement_from_doconsttimes)  \
    { if (0 < (count_from_doconsttimes)) { statement_from_doconsttimes; } \
      if (1 < (count_from_doconsttimes)) { statement_from_doconsttimes; } \
      if (2 < (count_from_doconsttimes)) { statement_from_doconsttimes; } \
      if (3 < (count_from_doconsttimes)) { statement_from_doconsttimes; } \
      if (4 < (count_from_doconsttimes)) { statement_from_doconsttimes; } \
      if (5 < (count_from_doconsttimes)) { statement_from_doconsttimes; } \
      if (6 < (count_from_doconsttimes)) { statement_from_doconsttimes; } \
      if (7 < (count_from_doconsttimes)) { statement_from_doconsttimes; } \
    }

// DOCONSTTIMES(count,macroname);
// ruft count mal den Macro macroname auf (count mal der Code!),
// wobei count eine constant-expression >=0, <=8 ist.
// Dabei bekommt macroname der Reihe nach die Werte 0,...,count-1 übergeben.
  #define DOCONSTTIMES(count_from_DOCONSTTIMES,macroname_from_DOCONSTTIMES)  \
    { if (0 < (count_from_DOCONSTTIMES)) { macroname_from_DOCONSTTIMES((0 < (count_from_DOCONSTTIMES) ? 0 : 0)); } \
      if (1 < (count_from_DOCONSTTIMES)) { macroname_from_DOCONSTTIMES((1 < (count_from_DOCONSTTIMES) ? 1 : 0)); } \
      if (2 < (count_from_DOCONSTTIMES)) { macroname_from_DOCONSTTIMES((2 < (count_from_DOCONSTTIMES) ? 2 : 0)); } \
      if (3 < (count_from_DOCONSTTIMES)) { macroname_from_DOCONSTTIMES((3 < (count_from_DOCONSTTIMES) ? 3 : 0)); } \
      if (4 < (count_from_DOCONSTTIMES)) { macroname_from_DOCONSTTIMES((4 < (count_from_DOCONSTTIMES) ? 4 : 0)); } \
      if (5 < (count_from_DOCONSTTIMES)) { macroname_from_DOCONSTTIMES((5 < (count_from_DOCONSTTIMES) ? 5 : 0)); } \
      if (6 < (count_from_DOCONSTTIMES)) { macroname_from_DOCONSTTIMES((6 < (count_from_DOCONSTTIMES) ? 6 : 0)); } \
      if (7 < (count_from_DOCONSTTIMES)) { macroname_from_DOCONSTTIMES((7 < (count_from_DOCONSTTIMES) ? 7 : 0)); } \
    }

// AT_INITIALIZATION(id) { ... }
// executes the given code at initialization time of the file.
// The id is something unique.
  #define AT_INITIALIZATION(id)  \
    class CONCAT3(INIT_CLASS_,id,__LINE__) {				\
      public: CONCAT3(INIT_CLASS_,id,__LINE__) (void);			\
    } CONCAT4(INIT_CLASS_,id,__LINE__,_DUMMY);				\
    inline CONCAT3(INIT_CLASS_,id,__LINE__)::CONCAT3(INIT_CLASS_,id,__LINE__) (void)

// AT_DESTRUCTION(id) { ... }
// executes the given code at destruction time of the file.
// The id is something unique.
  #define AT_DESTRUCTION(id)  \
    class CONCAT3(DESTR_CLASS_,id,__LINE__) {				\
      public: ~CONCAT3(DESTR_CLASS_,id,__LINE__) (void);		\
    } CONCAT4(DESTR_CLASS_,id,__LINE__,_DUMMY);				\
    CONCAT3(DESTR_CLASS_,id,__LINE__)::~CONCAT3(DESTR_CLASS_,id,__LINE__) (void)

// Inside a class definition:
// Overload `new' so that a class object can be allocated anywhere.
#if !((defined(__rs6000__) || defined(__alpha__)) && !defined(__GNUC__))
#define ALLOCATE_ANYWHERE(classname)  \
    /* Ability to place an object at a given address. */		\
public:									\
    void* operator new (size_t size) { return malloc_hook(size); }	\
    void* operator new (size_t size, classname* ptr) { unused size; return ptr; } \
    void operator delete (void* ptr) { free_hook(ptr); }
#else
// For some compilers, work around template problem with "classname".
#define ALLOCATE_ANYWHERE(classname)  \
    /* Ability to place an object at a given address. */		\
public:									\
    void* operator new (size_t size) { return malloc_hook(size); }	\
    void* operator new (size_t size, void* ptr) { unused size; return ptr; } \
    void operator delete (void* ptr) { free_hook(ptr); }
#endif

// init1(type, object) (value);
// initializes `object' with `value', by calling `type''s constructor.
// (The identifiers `init' and `Init' are already in use by <streambuf.h>,
// it's a shame!)
#define init1(type,lvalue)  (void) new (&(lvalue)) type

#include "base/cl_maybe_inline.h"

#endif /* _CL_MACROS_H */
