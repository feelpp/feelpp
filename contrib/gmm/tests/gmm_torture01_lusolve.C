// SQUARED_MATRIX_PARAM;
// DENSE_VECTOR_PARAM;
// VECTOR_PARAM;
// ENDPARAM;


#include <gmm_kernel.h>
#include <gmm_dense_lu.h>
#include <gmm_condition_number.h>

using gmm::size_type;
bool print_debug = false;

template <typename MAT1, typename VECT1, typename VECT2>
bool test_procedure(const MAT1 &m1_, const VECT1 &v1_, const VECT2 &v2_) {
  VECT1 &v1 = const_cast<VECT1 &>(v1_);
  VECT2 &v2 = const_cast<VECT2 &>(v2_);
  MAT1  &m1 = const_cast<MAT1  &>(m1_);
  typedef typename gmm::linalg_traits<MAT1>::value_type T;
  typedef typename gmm::number_traits<T>::magnitude_type R;
  R prec = gmm::default_tol(R());
  static size_type nb_iter = 0;

  size_type m = gmm::mat_nrows(m1);
  std::vector<T> v3(m);

  R det = gmm::abs(gmm::lu_det(m1)), error;
  R cond = gmm::condition_number(m1);

  if (print_debug) cout << "cond = " << cond << " det = " << det << endl;
  if (det == R(0) && cond < R(0.01) / prec && cond != R(0))
    DAL_THROW(gmm::failure_error, "Inconsistent condition number: " << cond);

  if (prec * cond < R(1)/R(10000) && det != R(0)) {
    ++nb_iter;

    gmm::lu_solve(m1, v1, v2);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);

    error = gmm::vect_norm2(v3);
    if (!(error <= prec * cond * R(20000)))
      DAL_THROW(gmm::failure_error, "Error too large: " << error);

    gmm::lu_inverse(m1);
    gmm::mult(m1, v2, v1);
    gmm::lu_inverse(m1);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    
    error = gmm::vect_norm2(v3);
    if (!(error <= prec * cond * R(20000)))
      DAL_THROW(gmm::failure_error, "Error too large: "<< error);

    if (nb_iter == 100) return true;
  }
  return false;
}
