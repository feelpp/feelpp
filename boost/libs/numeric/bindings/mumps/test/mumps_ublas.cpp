#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/traits/ublas_vector2.hpp>
#include <boost/numeric/bindings/traits/ublas_sparse.hpp>
#include <boost/numeric/bindings/mumps/mumps_driver.hpp>
#include <iostream>
#include <fstream>

int main() {
  namespace ublas = ::boost::numeric::ublas ;
  namespace mumps = ::boost::numeric::bindings::mumps ;


  int const n = 10 ;

  typedef ublas::coordinate_matrix<double, ublas::column_major, 1, ublas::unbounded_array<int> > coo_type ;

  coo_type coo( n, n, n + 6 ) ;

  for (int i=0; i<n; ++i) coo(i,i) = i+1.0 ;
  coo(2,3) = 1.0 ;
  coo(2,4) = 1.0 ;
  coo(5,6) = -1.0 ;
  coo(2,6) = 1.0 ;
  coo(9,0) = 1.0 ;
  coo(2,7) = -1.0 ;

  coo.sort() ;
  std::cout << "matrix " << coo << std::endl ;

  ublas::vector<double> v( 10 ) ;
  ublas::vector<double> w( 10 ) ;

  std::fill( w.begin(), w.end(), 1.0 ) ;

  for (int i=1; i<n; ++i) {
    w(i) += w(i-1) ;
  }

  for (int i=0; i<n; ++i) {
    v(i) = coo(i,i) * w(i) ;
  }
  v(2) += coo(2,3) * w(3) ;
  v(2) += coo(2,4) * w(4) ;
  v(5) += coo(5,6) * w(6) ;
  v(2) += coo(2,6) * w(6) ;
  v(9) += coo(9,0) * w(0) ;
  v(2) += coo(2,7) * w(7) ;

  mumps::mumps< coo_type > mumps_coo ;

  mumps_coo.icntl[2]=mumps_coo.icntl[3] = 0 ;

  // Analysis
  mumps_coo.job = 1 ;
  matrix_integer_data( mumps_coo, coo ) ;
  driver( mumps_coo ) ;

  // Factorization
  mumps_coo.job = 2 ;
  matrix_value_data( mumps_coo, coo ) ;
  driver( mumps_coo ) ;

  // Solve
  mumps_coo.job = 3 ;
  rhs_sol_value_data( mumps_coo, v ) ;
  driver( mumps_coo ) ;

  std::cout << "w : " << w << std::endl ;
  std::cout << "v : " << v << std::endl ;

}
