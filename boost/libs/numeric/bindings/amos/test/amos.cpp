#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/amos/amos.hpp>

int main()
{
  boost::numeric::ublas::vector< std::complex< double > > z(10) ; 
  boost::numeric::ublas::vector< std::complex< double > > cy(10) ; 
  int nz;
  boost::numeric::bindings::amos::besi( z, 1.0, 1, cy, nz ) ;

  return 0 ;
}
