//  Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.
//  Copyright Toon Knapen, Karl Meerbergen
//

#include "../../blas/test/random.hpp"

#include <boost/numeric/bindings/lapack/geqrf.hpp>
#include <boost/numeric/bindings/lapack/ormqr.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <iostream>
#include <limits>


namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;


// Randomize a matrix
template <typename M>
void randomize(M& m) {
   typedef typename M::size_type  size_type ;
   typedef typename M::value_type value_type ;

   size_type size1 = m.size1() ;
   size_type size2 = m.size2() ;

   for (size_type i=0; i<size2; ++i) {
      for (size_type j=0; j<size1; ++j) {
         m(j,i) = random_value< value_type >() ;
      }
   }
} // randomize()


template <typename T>
struct transpose {
   static const char value ;
};

template <typename T>
const char transpose<T>::value = 'T';



template <typename T>
struct transpose< std::complex<T> > {
   static const char value ;
};

template <typename T>
const char transpose< std::complex<T> >::value = 'C';

template <typename M>
ublas::triangular_adaptor<const M, ublas::upper> upper_part(const M& m) {
   return ublas::triangular_adaptor<const M, ublas::upper>( m );
}

template <typename T, typename W>
int do_memory_type(int n, W workspace) {
   typedef typename boost::numeric::bindings::traits::type_traits<T>::real_type real_type ;
   typedef std::complex< real_type >                                            complex_type ;

   typedef ublas::matrix<T, ublas::column_major> matrix_type ;
   typedef ublas::vector<T>                      vector_type ;

   // Set matrix
   matrix_type a( n, n );
   vector_type tau( n );

   randomize( a );
   matrix_type a2( a );

   // Compute QR factorization.
   lapack::geqrf( a, tau, workspace ) ;

   // Apply the orthogonal transformations to a2
   lapack::ormqr( 'L', transpose<T>::value, a, tau, a2, workspace );

   // The upper triangular parts of a and a2 must be equal.
   if (norm_frobenius( upper_part( a - a2 ) )
            > std::numeric_limits<real_type>::epsilon() * 10.0 * norm_frobenius( upper_part( a ) ) ) return 255 ;

   return 0 ;
} // do_value_type()



template <typename T>
int do_value_type() {
   const int n = 8 ;
   
   if (do_memory_type<T,lapack::optimal_workspace>( n, lapack::optimal_workspace() ) ) return 255 ;
   if (do_memory_type<T,lapack::minimal_workspace>( n, lapack::minimal_workspace() ) ) return 255 ;

   ublas::vector<T> work( n );
   do_memory_type<T, lapack::detail::workspace1<ublas::vector<T> > >( n, lapack::workspace(work) );
   return 0;
} // do_value_type()


int main() {
   // Run tests for different value_types
   if (do_value_type<float>()) return 255;
   if (do_value_type<double>()) return 255;
   if (do_value_type< std::complex<float> >()) return 255;
   if (do_value_type< std::complex<double> >()) return 255;

   std::cout << "Regression test succeeded\n" ;
   return 0;
}

