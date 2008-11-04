/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Dynamic Array Library (dal)                                  */
/* File    :  dal_std.h : Compatibility Header.                            */
/*                                                                         */
/* Date : June 01, 1995.                                                   */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr,                     */
/*          Julien Pommier, Julien.Pommier@gmm.insa-tlse.fr.               */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 1995-2005  Yves Renard.                                   */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */

#ifndef DAL_STD_H__
#define DAL_STD_H__

#ifndef NOGETFEM_VERIFY
#  define GETFEM_VERIFY
#endif

#ifndef __USE_STD_IOSTREAM
  #define __USE_STD_IOSTREAM
#endif

#ifndef __USE_BSD
  #define __USE_BSD
#endif

/* ********************************************************************** */
/*	Compilers detection.						  */
/* ********************************************************************** */

/* for sun CC 5.0 ...
#if defined(__SUNPRO_CC) && __SUNPRO_CC >= 0x500
# include <stdcomp.h>
# undef _RWSTD_NO_CLASS_PARTIAL_SPEC
# undef _RWSTD_NO_NAMESPACE
#endif 
*/
/* for VISUAL C++ ...
   #if defined(_MSC_VER) //  && !defined(__MWERKS__)
   #define _GETFEM_MSVCPP_ _MSC_VER
   #endif
*/

#if defined(__GNUC__)
#  if (__GNUC__ < 3)
#    error : PLEASE UPDATE g++ TO AT LEAST 3.0 VERSION
#  endif
#endif

/* ********************************************************************** */
/*	C++ Standard Headers.						  */
/* ********************************************************************** */
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <iostream>
//#include <ios> essai
#include <fstream>
#include <ctime>

/* ********************************************************************** */
/*	S.T.L. Headers.						          */
/* ********************************************************************** */

// #include <cstdlib>  CC de SGI ne reconnait pas ce header.
#include <exception>
#include <typeinfo>
#include <stdexcept>
#include <iterator>
#include <algorithm>
#include <vector>
#include <deque>
#include <string>
#include <complex>
#include <limits>
#include <sstream>

using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;

/* ********************************************************************** */
/*	Math functions.                     			          */
/* ********************************************************************** */

namespace dal
{
  template <typename T> inline T sqr(T a) { return a * a; }
  template <typename T> inline T abs(T a) { return (a < T(0)) ? T(-a) : a; }
  template <typename T> inline T abs(std::complex<T> a) { return std::abs(a); }
  template <typename T> inline T abs_sqr(T a) { return a*a; }
  template <typename T> inline T abs_sqr(std::complex<T> a)
  { return dal::sqr(a.real()) + dal::sqr(a.imag()); }
  template <typename T> inline T pos(T a) { return (a < T(0)) ? T(0) : a; }
  template <typename T> inline T neg(T a) { return (a < T(0)) ? T(-a) : T(0); }
  template <typename T> inline T sgn(T a) { return (a < T(0)) ? T(-1) : T(1); }
  inline double random() { return double(rand())/(RAND_MAX+0.5); }
  template <typename T> inline T random(T)
  { return T(rand()*2.0)/(T(RAND_MAX)+T(1)/T(2)) - T(1); }
  template <typename T> inline std::complex<T> random(std::complex<T>)
  { return std::complex<T>(dal::random(T()), dal::random(T())); }
  template <typename T> inline T irandom(T max)
  { return T(dal::random() * max); }
  template <typename T> inline T conj(T a) { return a; }
  template <typename T> inline std::complex<T> conj(std::complex<T> a)
  { return std::conj(a); }
  template <typename T> inline T real(T a) { return a; }
  template <typename T> inline T real(std::complex<T> a) { return a.real(); }
  template <typename T> inline T imag(T ) { return T(0); }
  template <typename T> inline T imag(std::complex<T> a) { return a.imag(); }  
  template <typename T> inline T sqrt(T a) { return ::sqrt(a); }
  template <typename T> inline std::complex<T> sqrt(std::complex<T> a) {
    T x = a.real(), y = a.imag();
    if (x == T(0)) {
      T t = ::sqrt(dal::abs(y) / T(2));
      return std::complex<T>(t, y < T(0) ? -t : t);
    }
    T t = ::sqrt(T(2) * (dal::abs(a) + dal::abs(x))), u = t / T(2);
    return x > T(0) ? std::complex<T>(u, y / t)
      : std::complex<T>(dal::abs(y) / t, y < T(0) ? -u : u);
  }


}


#ifndef M_PI
# define	M_E		2.7182818284590452354       /* e          */
# define	M_LOG2E		1.4426950408889634074       /* 1/ln(2)    */
# define	M_LOG10E	0.43429448190325182765      /* 1/ln(10)   */
# define	M_LN2		0.69314718055994530942      /* ln(2)      */
# define	M_LN10		2.30258509299404568402      /* ln(10)     */
# define	M_PI		3.14159265358979323846      /* pi         */
# define	M_PI_2		1.57079632679489661923      /* pi/2       */
# define	M_PI_4		0.78539816339744830962      /* pi/4       */
# define	M_1_PI		0.31830988618379067154      /* 1/pi       */
# define	M_2_PI		0.63661977236758134308      /* 2/pi       */
# define	M_2_SQRTPI	1.12837916709551257390      /* 2/sqrt(pi) */
# define	M_SQRT2		1.41421356237309504880      /* sqrt(2)    */
# define	M_SQRT1_2	0.70710678118654752440      /* sqrt(2)/2  */
#endif 

#ifndef M_PIl
# define M_PIl       3.1415926535897932384626433832795029L  /* pi         */
# define M_PI_2l     1.5707963267948966192313216916397514L  /* pi/2       */
# define M_PI_4l     0.7853981633974483096156608458198757L  /* pi/4       */
# define M_1_PIl     0.3183098861837906715377675267450287L  /* 1/pi       */
# define M_2_PIl     0.6366197723675813430755350534900574L  /* 2/pi       */
# define M_2_SQRTPIl 1.1283791670955125738961589031215452L  /* 2/sqrt(pi) */
#endif

namespace dal {

  /* ******************************************************************* */
  /*       Clock functions.                                              */
  /* ******************************************************************* */
  
# if  defined(HAVE_SYS_TIMES)
  inline double uclock_sec(void) {
    static double ttclk = 0.;
    if (ttclk == 0.) ttclk = sysconf(_SC_CLK_TCK);
    tms t; times(&t); return double(t.tms_utime) / ttclk;
  }
# else
  inline double uclock_sec(void)
  { return double(clock())/double(CLOCKS_PER_SEC); }
# endif
  
  /* ******************************************************************** */
  /*	Fixed size integer types.                     			  */
  /* ******************************************************************** */


typedef signed char    int8_type;
typedef unsigned char uint8_type;

#if INT_MAX == 32767
  typedef signed int    int16_type;
  typedef unsigned int uint16_type;
#elif  SHRT_MAX == 32767
  typedef signed short int    int16_type;
  typedef unsigned short int uint16_type;
#else
# error "impossible to build a 16bits integer"
#endif

#if INT_MAX == 2147483647
  typedef signed int    int32_type;
  typedef unsigned int uint32_type;
#elif  SHRT_MAX == 2147483647
  typedef signed short int    int32_type;
  typedef unsigned short int uint32_type;
#elif LONG_MAX == 2147483647
  typedef signed long int    int32_type;
  typedef unsigned long int uint32_type;
#else
# error "impossible to build a 32bits integer"
#endif

#if INT_MAX == 9223372036854775807L || INT_MAX == 9223372036854775807
  typedef signed int    int64_type;
  typedef unsigned int uint64_type;
#elif LONG_MAX == 9223372036854775807L || LONG_MAX == 9223372036854775807
  typedef signed long int    int64_type;
  typedef unsigned long int uint64_type;
  // utiliser long long dans le futur.
  // #elif LLONG_MAX == 9223372036854775807LL || LLONG_MAX == 9223372036854775807L || LLONG_MAX == 9223372036854775807
  // typedef signed long long int int64_type;
  // typedef unsigned long long int uint64_type;
#endif

#ifdef __GNUC__
/* 
   g++ can issue a warning at each usage of a function declared with this special attribute 
   (also works with typedefs and variable declarations)
*/
# define IS_DEPRECATED __attribute__ ((__deprecated__))
/*
   the specified function is inlined at any optimization level 
*/
# define ALWAYS_INLINE __attribute__((always_inline))
#else
# define IS_DEPRECATED
# define ALWAYS_INLINE
#endif

  /* Pour forcer l'instanciation dans libgetfem de tous les
     type qu'on est suceptible d' afficher. Si il en manque un (c'est
     le cas pour unsigned et float) dans libgetfem.so qui est utilise
     dans libgetfemint.so alors on retrouve un warning
     comme quoi le symbole T_Q13std8ios_base__ et T_Q13std11logic_error__
     sont dupliqués
     visiblement il reinstancie un peu trop de choses .. 

     pour s'assurer que le probleme n'est pas de retour, faire
     nm libgetfem.so | grep 'td::basic_ostream<char,
                              std::char_traits<char> >::operator <<'
     et faire de même dans tous les .o de matlabint, et s'assurer que
     tous ceux de matlabint sont bien inclus dans libgetfem.so

     enfin c'est pas tout a fait ça puisque tous ces trucs sont deja
     instancies dans libcxxstd.a, et c'est le noeud du probleme j'imagine:
     nm /usr/lib/cmplrs/cxx/V6.3-008/libcxxstd.a
       | grep 'td::basic_ostream<char, std::char_traits<char> >::operator <<' 
     ---> std::basic_ostream<char, std::char_traits<char> >::operator
            <<(long long) | 0000000000006368 | T | 0000000000000008
     std::basic_ostream<char, std::char_traits<char> >::operator
         <<(const void*) | 0000000000006368 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator
         <<(std::basic_ostream<char, std::char_traits<char> >&
	 (*)(std::basic_ostream<char, std::char_traits<char> >&))
	 | 0000000000000000 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator
         <<(unsigned long long) | 0000000000006368 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator
         <<(unsigned int) | 0000000000006368 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator
         <<(unsigned long) | 0000000000006368 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator
         <<(unsigned short) | 0000000000006368 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator 
         <<(double) | 0000000000006368 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator
         <<(float) | 0000000000006368 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator
         <<(int) | 0000000000006368 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator
         <<(bool) | 0000000000006368 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator
         <<(long double128) | 0000000000006368 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator
         <<(long) | 0000000000006368 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator
         <<(short) | 0000000000006368 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator
         <<(long double64) | 0000000000006368 | T | 0000000000000008
  */
  struct just_for_the_fine_cxx {
    static void f() {
      long double z(1.0);
      std::stringstream s;
      s << int(1) << double(2.0) 
	<< "hello" << std::string("hello") << unsigned(1) << float(2) 
	<< char('a') << (unsigned char)('b') << short(1)
	<< (unsigned short)(2)  << long(1) << (unsigned long)(2)
	<< (const void*)NULL << bool(1) << z;
    }
  };
  
} /* end of namespace dal.                                                */

#endif /* DAL_STD_H__ */

