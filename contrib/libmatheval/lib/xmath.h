/* 
 * Copyright (C) 1999, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2011
 * Free Software Foundation, Inc.
 * 
 * This file is part of GNU libmatheval
 * 
 * GNU libmatheval is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * GNU libmatheval is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU libmatheval.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef XMATH_H
#define XMATH_H 1

#if HAVE_CONFIG_H
#include "config.h"
#endif

#if HAVE_MATH_H
#include <math.h>
#else
#error no <math.h> avaialable
#endif

#ifdef INFINITY
#define MATH_INFINITY INFINITY
#elif defined HUGE_VAL
#define MATH_INFINITY HUGE_VAL
#else
#define MATH_INFINITY (1.0/0.0)
#endif

#ifdef NAN
#define MATH_NAN NAN
#elif defined INFINITY
#define MATH_NAN (INFINITY/INFINITY)
#else
#define MATH_NAN (0.0/0.0)
#endif

/* Calculate cotangent of value x.  */
double          math_cot(double x);

/* Calculate secant of value x.  */
double          math_sec(double x);

/* Calculate cosecant of value x.  */
double          math_csc(double x);

/* Calculate inverse cotangent of value x.  */
double          math_acot(double x);

/* Calculate inverse secant of value x.  */
double          math_asec(double x);

/* Calculate inverse cosecant of value x.  */
double          math_acsc(double x);

/* Calculate hyperbolical cotangent of value x.  */
double          math_coth(double x);

/* Calculate hyperbolical secant of value x.  */
double          math_sech(double x);

/* Calculate hyperbolical cosecant of value x.  */
double          math_csch(double x);

/* Calculate inverse hyperbolical sine of value x.  */
double          math_asinh(double x);

/* Calculate inverse hyperbolical cosine of value x.  */
double          math_acosh(double x);

/* Calculate inverse hyperbolical tangent of value x.  */
double          math_atanh(double x);

/* Calculate inverse hyperbolical cotangent of value x.  */
double          math_acoth(double x);

/* Calculate inverse hyperbolical secant of value x.  */
double          math_asech(double x);

/* Calculate inverse hyperbolical cosecant of value x.  */
double          math_acsch(double x);

/* Calculate Heaviside step function value for given value x.  */
double          math_step(double x);

/* Calculate Dirac delta function value for given value x.  */
double          math_delta(double x);

/* Calculate variation of Dirac delta function (with not-a-number instead
 * of infinity value for x= 0) value for given value x. */
double          math_nandelta(double x);

#endif
