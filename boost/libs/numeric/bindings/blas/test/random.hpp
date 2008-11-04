//  Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.
//  Copyright Toon Knapen, Karl Meerbergen

#include <cstdlib>
#include <complex>
#include <assert.h>


template <typename T>
T random_value() {
   assert( false );
   return 0;
}

template <>
float random_value<float>() {return float(std::rand()) / float(RAND_MAX) - 0.5;}

template <>
double random_value<double>() {return double(std::rand()) / double(RAND_MAX) - 0.5;}

template <>
std::complex<float> random_value< std::complex<float> >() {return std::complex<float>(random_value<float>(), random_value<float>() );}

template <>
std::complex<double> random_value< std::complex<double> >() {return std::complex<double>(random_value<double>(), random_value<double>() );}

