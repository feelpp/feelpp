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
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <random>
#include <iostream>
#define Pi M_PI
#define pi M_PI

inline std::pair<std::uniform_real_distribution<>, std::mt19937&> uniformDistribution( double a, double b )
{
    static std::random_device rd;
    static std::mt19937 gen( rd() );
    std::uniform_real_distribution<> dis( a, b );
    return std::pair<std::uniform_real_distribution<>&, std::mt19937&>{ dis, gen };
}
inline std::pair<std::normal_distribution<>, std::mt19937&> normalDistribution( double m, double s )
{
    static std::random_device rd;
    static std::mt19937 gen( rd() );
    std::normal_distribution<> dis{ m, s };
    return std::pair<std::normal_distribution<>&, std::mt19937&>{ dis, gen };
}
inline std::pair<std::lognormal_distribution<>, std::mt19937&> lognormalDistribution( double m, double s )
{
    static std::random_device rd;
    static std::mt19937 gen( rd() );
    std::lognormal_distribution<> dis{ m, s };
    return std::pair<std::lognormal_distribution<>&, std::mt19937&>{ dis, gen };
}

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
