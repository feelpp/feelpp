//  Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.
//  Copyright Toon Knapen, Karl Meerbergen
//

#include "../../blas/test/random.hpp"

#include <boost/numeric/bindings/lapack/hbev.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_banded.hpp>
#include <boost/numeric/bindings/traits/ublas_hermitian.hpp>
#include <boost/numeric/ublas/hermitian.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <iostream>
#include <limits>


namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;


// Fill a matrix
template <typename M>
void fill(M& m) {
   typedef typename M::size_type  size_type ;
   typedef typename M::value_type value_type ;

   typedef typename boost::numeric::bindings::traits::type_traits<value_type>::real_type real_type ;

   int size = m.size2() ;
   int band = m.upper() ;

   for (int i=0; i<size; ++i) {
      for (int j=std::max(0,i-band); j<i; ++j) m(j,i) = random_value<value_type>();
      m(i,i) = random_value<real_type>();
   }
} // randomize()


template <typename U>
int lower() {
   return 0;
}

template <>
int lower<ublas::lower>() { return 2; }


template <typename U>
int upper() {
   return 0;
}

template <>
int upper<ublas::upper>() { return 2; }

template <typename H, typename E, typename Z>
int check_residual(H const& h, E const& e, Z const& z) {
   typedef typename H::value_type value_type ;
   typedef typename E::value_type real_type ;

   // Check eigen decomposition
   int n = h.size1();
   ublas::matrix<value_type> error( n, n ); error.clear();

   // Copy band matrix in error
   error.assign( h );
   assert( norm_frobenius( error - herm( error ) ) == 0.0 ) ;

   for (int i=0; i<n; ++i) {
      error .minus_assign( outer_prod( z.column(i), e(i) * conj( z.column(i) ) ) ) ;
   }
   return (norm_frobenius( error )
           >= n* norm_2( e ) * std::numeric_limits< real_type >::epsilon() ) ;
} // check_residual()


template <typename T, typename W, typename UPLO, typename Orientation>
int do_memory_uplo(int n, W& workspace ) {
   typedef typename boost::numeric::bindings::traits::type_traits<T>::real_type real_type ;

   typedef ublas::banded_matrix<T, Orientation>      banded_type ;
   typedef ublas::matrix<T, ublas::column_major>     matrix_type ;
   typedef ublas::vector<real_type>                  vector_type ;

   // Set matrix
   banded_type a( n, n, lower<UPLO>(), upper<UPLO>() ); a.clear(); // Make upper triangular band matrix
   matrix_type z( n, n );
   vector_type e1( n );
   vector_type e2( n );

   fill( a );
   banded_type a2( a );

   ublas::hermitian_adaptor<banded_type, UPLO> h( a ), h2( a2 );

   // Compute Schur decomposition.
   lapack::hbev( h, e1, z, workspace ) ;

   if (check_residual( h2, e1, z )) return 255 ;

   lapack::hbev( h2, e2, workspace ) ;
   if (norm_2( e1 - e2 ) > n * norm_2( e1 ) * std::numeric_limits< real_type >::epsilon()) return 255 ;

   // Test for a matrix range
   fill( a ); a2.assign( a );

   typedef ublas::matrix_range< banded_type > banded_range ;

   ublas::range r(1,n-1) ;
   banded_range a_r( a, r, r );
   ublas::hermitian_adaptor< banded_range, UPLO> h_r( a_r );
   ublas::vector_range< vector_type> e_r( e1, r );
   ublas::matrix_range< matrix_type> z_r( z, r, r );

   lapack::hbev( h_r, e_r, z_r, workspace );

   if (check_residual( h2(r,r), e_r, z_r )) return 255 ;

   return 0 ;
} // do_memory_uplo()


template <typename T, typename W>
int do_memory_type(int n, W workspace) {
   std::cout << "  upper\n" ;
   if (do_memory_uplo<T,W,ublas::upper, ublas::row_major>(n, workspace)) return 255 ;
   std::cout << "  lower\n" ;
   if (do_memory_uplo<T,W,ublas::lower, ublas::row_major>(n, workspace)) return 255 ;
   return 0 ;
}


template <typename T>
struct Workspace {
   typedef ublas::vector<T>                         array_type ;
   typedef lapack::detail::workspace1< array_type > type ;

   Workspace(size_t n)
   : work_( 3*n-2 )
   {}

   type operator() () {
      return type( work_ );
   }

   array_type work_ ;
};


template <typename T>
struct Workspace< std::complex<T> > {
   typedef ublas::vector<T>                                                 real_array_type ;
   typedef ublas::vector< std::complex<T> >                                 complex_array_type ;
   typedef lapack::detail::workspace2< complex_array_type,real_array_type > type ;

   Workspace(size_t n)
   : work_( n )
   , rwork_( 3*n-2 )
   {}

   type operator() () {
      return type( work_, rwork_ );
   }

   complex_array_type work_ ;
   real_array_type    rwork_ ;
};


template <typename T>
int do_value_type() {
   const int n = 8 ;

   std::cout << " optimal workspace\n";
   if (do_memory_type<T,lapack::optimal_workspace>( n, lapack::optimal_workspace() ) ) return 255 ;

   std::cout << " minimal workspace\n";
   if (do_memory_type<T,lapack::minimal_workspace>( n, lapack::minimal_workspace() ) ) return 255 ;

   std::cout << " workspace array\n";
   Workspace<T> work( n );
   do_memory_type<T,typename Workspace<T>::type >( n, work() );
   return 0;
} // do_value_type()


int main() {
   // Run tests for different value_types
   std::cout << "float\n" ;
   if (do_value_type<float>()) return 255;

   std::cout << "double\n" ;
   if (do_value_type<double>()) return 255;

   std::cout << "complex<float>\n" ;
   if (do_value_type< std::complex<float> >()) return 255;

   std::cout << "complex<double>\n" ;
   if (do_value_type< std::complex<double> >()) return 255;

   std::cout << "Regression test succeeded\n" ;
   return 0;
}

