/** @file ginac-extern-math.h
 *
 *  Declaration of some functions from math.h required with using ginac-excompiler
 *  on systems where it is not possible to compile on running nodes.
 *
 */

//extern const constant Pi;
//extern const constant Catalan;
//extern const constant Euler;

#if defined( USE_STANDARD_HEADERS_IN_GINAC_EXCOMPILER )
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#define Pi M_PI
#define pi M_PI
#else

extern double sin(double x);
extern double cos(double x);
extern double tan(double x);
extern double acos(double x);
extern double asin(double x);
extern double atan2(double y, double x);
extern double hypot(double x, double y);

extern double exp(double x);
extern double log(double x);

extern double sqrt(double x);
extern double pow(double x, double y);

#endif
