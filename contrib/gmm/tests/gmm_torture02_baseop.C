// SQUARED_MATRIX_PARAM;
// VECTOR_PARAM;
// ENDPARAM;

#include <gmm_kernel.h>

using gmm::size_type;
bool print_debug = false;

template <typename MAT1, typename VECT1>
bool test_procedure(const MAT1 &m1_, const VECT1 &v1_) {
  VECT1 &v1 = const_cast<VECT1 &>(v1_);
  MAT1  &m1 = const_cast<MAT1  &>(m1_);
  typedef typename gmm::linalg_traits<MAT1>::value_type T;
  typedef typename gmm::number_traits<T>::magnitude_type R;
  R prec = gmm::default_tol(R());

  size_type m = gmm::mat_nrows(m1);

  R norm = gmm::vect_norm2_sqr(v1);

  R normtest(0);

  for (size_type i = 0; i < m; ++i) {
    T x(1), y = v1[i];;
    x *= v1[i];
    x += v1[i];
    x += v1[i];
    x -= v1[i];
    x -= y;
    x *= v1[i];
    x /= v1[i];
    if (y != v1[i])
      DAL_THROW(gmm::failure_error, "Error in basic operations");
    if (!(y == v1[i]))
      DAL_THROW(gmm::failure_error, "Error in basic operations");
    normtest += gmm::abs_sqr(x);
  }
  
  if (gmm::abs(norm - normtest) > prec * R(100))
    DAL_THROW(gmm::failure_error, "Error in basic operations");
  
  return true;
}
