// SQUARED_MATRIX_PARAM
// VECTOR_PARAM;
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
  R error, det, cond;
  static  size_type nb_iter(0);
  ++nb_iter;

  size_type m = gmm::vect_size(v1), n = m/2;
  std::vector<T> v3(n);

  det = gmm::abs(gmm::lu_det(gmm::sub_matrix(m1, gmm::sub_interval(0,n))));
  cond = gmm::condition_number(gmm::sub_matrix(m1, gmm::sub_interval(0,n)));
  if (prec * cond < R(1)/R(10000) && det != R(0)) {
    gmm::lu_solve(gmm::sub_matrix(m1, gmm::sub_interval(0,n)), v3,
		  gmm::sub_vector(v2, gmm::sub_interval(0,n)));
    gmm::mult(gmm::sub_matrix(m1, gmm::sub_interval(0,n)), v3,
	      gmm::sub_vector(v1, gmm::sub_interval(0,n)));
    gmm::add(gmm::scaled(gmm::sub_vector(v1, gmm::sub_interval(0,n)), T(-1)),
	     gmm::sub_vector(v2, gmm::sub_interval(0,n)), v3);
    if (!((error = gmm::vect_norm2(v3)) <= prec * R(20000) * cond))
      DAL_THROW(gmm::failure_error, "Error too large: "<< error);
  }

  det = gmm::abs(gmm::lu_det(gmm::sub_matrix(m1, gmm::sub_slice(0,n,1))));
  cond = gmm::condition_number(gmm::sub_matrix(m1, gmm::sub_slice(0,n,1)));
  if (prec * cond < R(1)/R(10000) && det != R(0)) {
    gmm::lu_solve(gmm::sub_matrix(m1, gmm::sub_slice(0,n,1)), v3,
		  gmm::sub_vector(v2, gmm::sub_slice(0,n,1)));
    gmm::mult(gmm::sub_matrix(m1, gmm::sub_slice(0,n,1)), v3,
	      gmm::sub_vector(v1, gmm::sub_slice(0,n,1)));
    gmm::add(gmm::scaled(gmm::sub_vector(v1, gmm::sub_slice(0,n,1)), T(-1)),
	     gmm::sub_vector(v2, gmm::sub_slice(0,n,1)), v3);
    if (!((error = gmm::vect_norm2(v3)) <= prec * R(20000) * cond))
      DAL_THROW(gmm::failure_error, "Error too large: "<< error);
  }
  
  gmm::copy(gmm::identity_matrix(), gmm::sub_matrix(gmm::transposed(m1),
		        gmm::sub_interval(0,n), gmm::sub_interval(0,m)));
  gmm::clear(gmm::sub_vector(v2, gmm::sub_interval(n, m-n)));
  if (print_debug) {
    cout << "sub matrix of m1 : "
	 << gmm::sub_matrix(gmm::transposed(m1), gmm::sub_interval(0,n))
	 << endl;
  }

  gmm::mult(gmm::sub_matrix(m1, gmm::sub_interval(0,m), 
			    gmm::sub_interval(0,n)),
	    gmm::sub_vector(v2, gmm::sub_interval(0,n)),
	    gmm::scaled(v2, T(-1)), v1);
  if (!((error = gmm::vect_norm2(v1)) <= prec * R(2000)))
    DAL_THROW(gmm::failure_error, "Error too large: " << error);
  
  if (nb_iter == 100) return true;
  return false;
}
