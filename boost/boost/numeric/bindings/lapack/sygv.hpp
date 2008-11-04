//
// Copyright Fabien Dekeyser, Quoc-Cuong Pham 2008
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_SYGV_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_SYGV_HPP

#include <boost/numeric/bindings/traits/traits.hpp>
#include <boost/numeric/bindings/lapack/lapack.h>
#include <boost/numeric/bindings/lapack/workspace.hpp>
#include <boost/numeric/bindings/traits/detail/array.hpp>
// #include <boost/numeric/bindings/traits/std_vector.hpp>

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK 
#  include <boost/static_assert.hpp>
#  include <boost/type_traits.hpp>
#endif 

#include <cassert>


namespace boost { namespace numeric { namespace bindings { 

  namespace lapack {

    ///////////////////////////////////////////////////////////////////
    //
    // sygv
    // 
    ///////////////////////////////////////////////////////////////////

	  /* 
	  *  sygv() computes all the eigenvalues, and optionally, the eigenvectors
	  *  of a real generalized symmetric-definite eigenproblem, of the form
	  *  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
	  *  Here A and B are assumed to be symmetric and B is also
	  *  positive definite.
	  * TYPE   (input) INTEGER
	  *          Specifies the problem type to be solved:
	  *          = 1:  A*x = (lambda)*B*x
	  *          = 2:  A*B*x = (lambda)*x
	  *          = 3:  B*A*x = (lambda)*x
	  *
	  *  JOBZ    (input) CHARACTER*1
	  *          = 'N':  Compute eigenvalues only;
	  *          = 'V':  Compute eigenvalues and eigenvectors.
	  *
	  *  UPLO    (input) CHARACTER*1
	  *          = 'U':  Upper triangles of A and B are stored;
	  *          = 'L':  Lower triangles of A and B are stored.
	  *
	  *  N       (input) INTEGER
	  *          The order of the matrices A and B.  N >= 0.
	  *
	  *  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
	  *          On entry, the symmetric matrix A.  If UPLO = 'U', the
	  *          leading N-by-N upper triangular part of A contains the
	  *          upper triangular part of the matrix A.  If UPLO = 'L',
	  *          the leading N-by-N lower triangular part of A contains
	  *          the lower triangular part of the matrix A.
	  *
	  *          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
	  *          matrix Z of eigenvectors.  The eigenvectors are normalized
	  *          as follows:
	  *          if ITYPE = 1 or 2, Z**T*B*Z = I;
	  *          if ITYPE = 3, Z**T*inv(B)*Z = I.
	  *          If JOBZ = 'N', then on exit the upper triangle (if UPLO='U')
	  *          or the lower triangle (if UPLO='L') of A, including the
	  *          diagonal, is destroyed.
	  *
	  *  LDA     (input) INTEGER
	  *          The leading dimension of the array A.  LDA >= max(1,N).
	  *
	  *  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
	  *          On entry, the symmetric positive definite matrix B.
	  *          If UPLO = 'U', the leading N-by-N upper triangular part of B
	  *          contains the upper triangular part of the matrix B.
	  *          If UPLO = 'L', the leading N-by-N lower triangular part of B
	  *          contains the lower triangular part of the matrix B.
	  *
	  *          On exit, if INFO <= N, the part of B containing the matrix is
	  *          overwritten by the triangular factor U or L from the Cholesky
	  *          factorization B = U**T*U or B = L*L**T.
	  *
	  *  LDB     (input) INTEGER
	  *          The leading dimension of the array B.  LDB >= max(1,N).
	  *
	  *  W       (output) DOUBLE PRECISION array, dimension (N)
	  *          If INFO = 0, the eigenvalues in ascending order.
	  *
	  *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
	  *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
	  *
	  *  LWORK   (input) INTEGER
	  *          The length of the array WORK.  LWORK >= max(1,3*N-1).
	  *          For optimal efficiency, LWORK >= (NB+2)*N,
	  *          where NB is the blocksize for DSYTRD returned by ILAENV.
	  *
	  *          If LWORK = -1, then a workspace query is assumed; the routine
	  *          only calculates the optimal size of the WORK array, returns
	  *          this value as the first entry of the WORK array, and no error
	  *          message related to LWORK is issued by XERBLA.
	  *
	  *  INFO    (output) INTEGER
	  *          = 0:  successful exit
	  *          < 0:  if INFO = -i, the i-th argument had an illegal value
	  *          > 0:  DPOTRF or DSYEV returned an error code:
	  *             <= N:  if INFO = i, DSYEV failed to converge;
	  *                    i off-diagonal elements of an intermediate
	  *                    tridiagonal form did not converge to zero;
	  *             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
	  *                    minor of order i of B is not positive definite.
	  *                    The factorization of B could not be completed and
	  *                    no eigenvalues or eigenvectors were computed.
	  *
	  */ 

    namespace detail {

      inline 
      void sygv (int const itype, char const jobz, char const uplo, int const n, 
					float *a, int const lda, float *b, int const ldb, 
					float *w, float *work, int const lwork, int& info)
	  {
	 
        LAPACK_SSYGV (&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, &info);
      }

      inline 
      void sygv (int const itype, char const jobz, char const uplo, int const n, 
					double *a, int const lda, double *b, int const ldb, 
					double *w, double *work, int const lwork, int& info)
      {
        LAPACK_DSYGV (&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, &info);
      }


      template <typename A, typename B, typename W, typename Work>
      inline
      int sygv (int itype, char jobz, char uplo, A& a, B& b, W& w, Work& work) {

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK 
        BOOST_STATIC_ASSERT((boost::is_same<
          typename traits::matrix_traits<A>::matrix_structure, 
          traits::general_t
        >::value)); 
#endif 

        int const n = traits::matrix_size1 (a);
        assert ( n>0 );
        assert (traits::matrix_size2 (a)==n); 
        assert (traits::leading_dimension (a)>=n); 
        assert (traits::vector_size (w)==n); 
        assert (3*n-1 <= traits::vector_size (work)); 

		int const nb = traits::matrix_size1 (b);
        assert ( nb>0 );
		assert (traits::matrix_size2 (b)==nb); 
        assert (traits::leading_dimension (b)>=nb); 
		assert ( n== nb);

        assert ( uplo=='U' || uplo=='L' );
        assert ( jobz=='N' || jobz=='V' );

		assert( itype==1 || itype==2 || itype==3);


        int info; 
        detail::sygv (itype, jobz, uplo, n,
                     traits::matrix_storage (a), 
                     traits::leading_dimension (a),
					 traits::matrix_storage (b), 
                     traits::leading_dimension (b),
                     traits::vector_storage (w),  
                     traits::vector_storage (work),
                     traits::vector_size (work),
                     info);
        return info; 
      }
    }  // namespace detail


    // Function that allocates work arrays
    template <typename A, typename B, typename W>
    inline
    int sygv (int itype, char jobz, char uplo, A& a, B& b, W& w, optimal_workspace ) {
       typedef typename A::value_type value_type ;

       int const n = traits::matrix_size1 (a);
       traits::detail::array<value_type> work( std::max<int>(1,34*n) );
       return detail::sygv(itype, jobz, uplo, a, b, w, work);
    } // sygv()


    // Function that allocates work arrays
    template <typename A, typename B, typename W>
    inline
    int sygv (int itype, char jobz, char uplo, A& a, B& b, W& w, minimal_workspace ) {
       typedef typename A::value_type value_type ;

       int const n = traits::matrix_size1 (a);

       traits::detail::array<value_type> work( std::max<int>(1,3*n-1) );
       return detail::sygv(itype, jobz, uplo, a, b, w, work);
    } // sygv()


    // Function that allocates work arrays
    template <typename A, typename B, typename W, typename Work>
    inline
    int sygv (int itype, char jobz, char uplo, A& a, B& b, W& w, detail::workspace1<Work> workspace ) {
       typedef typename A::value_type value_type ;

       return detail::sygv(itype, jobz, uplo, a, b, w, workspace.w_);
    } // sygv()

	template <typename A, typename B, typename W>
    inline
    int sygv (int itype, char jobz, char uplo, A& a, B& b, W& w) {
       return sygv(itype, jobz, uplo, a, b, w, optimal_workspace());
    } // sygv()

  }

}}}

#endif 
