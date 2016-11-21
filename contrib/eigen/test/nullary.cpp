// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2010-2011 Jitse Niesen <jitse@maths.leeds.ac.uk>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "main.h"

template<typename MatrixType>
bool equalsIdentity(const MatrixType& A)
{
  typedef typename MatrixType::Scalar Scalar;
  Scalar zero = static_cast<Scalar>(0);

  bool offDiagOK = true;
  for (Index i = 0; i < A.rows(); ++i) {
    for (Index j = i+1; j < A.cols(); ++j) {
      offDiagOK = offDiagOK && (A(i,j) == zero);
    }
  }
  for (Index i = 0; i < A.rows(); ++i) {
    for (Index j = 0; j < (std::min)(i, A.cols()); ++j) {
      offDiagOK = offDiagOK && (A(i,j) == zero);
    }
  }

  bool diagOK = (A.diagonal().array() == 1).all();
  return offDiagOK && diagOK;
}

template<typename VectorType>
void testVectorType(const VectorType& base)
{
  typedef typename VectorType::Scalar Scalar;

  const Index size = base.size();
  
  Scalar high = internal::random<Scalar>(-500,500);
  Scalar low = (size == 1 ? high : internal::random<Scalar>(-500,500));
  if (low>high) std::swap(low,high);

  const Scalar step = ((size == 1) ? 1 : (high-low)/(size-1));

  // check whether the result yields what we expect it to do
  VectorType m(base);
  m.setLinSpaced(size,low,high);

  if(!NumTraits<Scalar>::IsInteger)
  {
    VectorType n(size);
    for (int i=0; i<size; ++i)
      n(i) = low+i*step;
    VERIFY_IS_APPROX(m,n);
  }

  VectorType n(size);
  for (int i=0; i<size; ++i)
    n(i) = size==1 ? low : (low + ((high-low)*Scalar(i))/(size-1));
  VERIFY_IS_APPROX(m,n);

  // random access version
  m = VectorType::LinSpaced(size,low,high);
  VERIFY_IS_APPROX(m,n);

  VERIFY( internal::isApprox(m(m.size()-1),high) );
  VERIFY( size==1 || internal::isApprox(m(0),low) );

  // sequential access version
  m = VectorType::LinSpaced(Sequential,size,low,high);
  VERIFY_IS_APPROX(m,n);

  VERIFY( internal::isApprox(m(m.size()-1),high) );
  VERIFY( size==1 || internal::isApprox(m(0),low) );

  // check whether everything works with row and col major vectors
  Matrix<Scalar,Dynamic,1> row_vector(size);
  Matrix<Scalar,1,Dynamic> col_vector(size);
  row_vector.setLinSpaced(size,low,high);
  col_vector.setLinSpaced(size,low,high);
  // when using the extended precision (e.g., FPU) the relative error might exceed 1 bit
  // when computing the squared sum in isApprox, thus the 2x factor.
  VERIFY( row_vector.isApprox(col_vector.transpose(), Scalar(2)*NumTraits<Scalar>::epsilon()));

  Matrix<Scalar,Dynamic,1> size_changer(size+50);
  size_changer.setLinSpaced(size,low,high);
  VERIFY( size_changer.size() == size );

  typedef Matrix<Scalar,1,1> ScalarMatrix;
  ScalarMatrix scalar;
  scalar.setLinSpaced(1,low,high);
  VERIFY_IS_APPROX( scalar, ScalarMatrix::Constant(high) );
  VERIFY_IS_APPROX( ScalarMatrix::LinSpaced(1,low,high), ScalarMatrix::Constant(high) );

  // regression test for bug 526 (linear vectorized transversal)
  if (size > 1) {
    m.tail(size-1).setLinSpaced(low, high);
    VERIFY_IS_APPROX(m(size-1), high);
  }
}

template<typename MatrixType>
void testMatrixType(const MatrixType& m)
{
  using std::abs;
  const Index rows = m.rows();
  const Index cols = m.cols();
  typedef typename MatrixType::Scalar Scalar;
  typedef typename MatrixType::RealScalar RealScalar;

  Scalar s1;
  do {
    s1 = internal::random<Scalar>();
  } while(abs(s1)<RealScalar(1e-5) && (!NumTraits<Scalar>::IsInteger));

  MatrixType A;
  A.setIdentity(rows, cols);
  VERIFY(equalsIdentity(A));
  VERIFY(equalsIdentity(MatrixType::Identity(rows, cols)));


  A = MatrixType::Constant(rows,cols,s1);
  Index i = internal::random<Index>(0,rows-1);
  Index j = internal::random<Index>(0,cols-1);
  VERIFY_IS_APPROX( MatrixType::Constant(rows,cols,s1)(i,j), s1 );
  VERIFY_IS_APPROX( MatrixType::Constant(rows,cols,s1).coeff(i,j), s1 );
  VERIFY_IS_APPROX( A(i,j), s1 );
}

void test_nullary()
{
  CALL_SUBTEST_1( testMatrixType(Matrix2d()) );
  CALL_SUBTEST_2( testMatrixType(MatrixXcf(internal::random<int>(1,300),internal::random<int>(1,300))) );
  CALL_SUBTEST_3( testMatrixType(MatrixXf(internal::random<int>(1,300),internal::random<int>(1,300))) );
  
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_4( testVectorType(VectorXd(internal::random<int>(1,300))) );
    CALL_SUBTEST_5( testVectorType(Vector4d()) );  // regression test for bug 232
    CALL_SUBTEST_6( testVectorType(Vector3d()) );
    CALL_SUBTEST_7( testVectorType(VectorXf(internal::random<int>(1,300))) );
    CALL_SUBTEST_8( testVectorType(Vector3f()) );
    CALL_SUBTEST_8( testVectorType(Vector4f()) );
    CALL_SUBTEST_8( testVectorType(Matrix<float,8,1>()) );
    CALL_SUBTEST_8( testVectorType(Matrix<float,1,1>()) );

    CALL_SUBTEST_9( testVectorType(VectorXi(internal::random<int>(1,300))) );
    CALL_SUBTEST_9( testVectorType(Matrix<int,1,1>()) );
  }

#ifdef EIGEN_TEST_PART_6
  // Assignment of a RowVectorXd to a MatrixXd (regression test for bug #79).
  VERIFY( (MatrixXd(RowVectorXd::LinSpaced(3, 0, 1)) - RowVector3d(0, 0.5, 1)).norm() < std::numeric_limits<double>::epsilon() );
#endif

#ifdef EIGEN_TEST_PART_10
  // check some internal logic
  VERIFY((  internal::has_nullary_operator<internal::scalar_constant_op<double> >::value ));
  VERIFY(( !internal::has_unary_operator<internal::scalar_constant_op<double> >::value ));
  VERIFY(( !internal::has_binary_operator<internal::scalar_constant_op<double> >::value ));
  VERIFY((  internal::functor_has_linear_access<internal::scalar_constant_op<double> >::ret ));

  VERIFY(( !internal::has_nullary_operator<internal::scalar_identity_op<double> >::value ));
  VERIFY(( !internal::has_unary_operator<internal::scalar_identity_op<double> >::value ));
  VERIFY((  internal::has_binary_operator<internal::scalar_identity_op<double> >::value ));
  VERIFY(( !internal::functor_has_linear_access<internal::scalar_identity_op<double> >::ret ));

  VERIFY(( !internal::has_nullary_operator<internal::linspaced_op<float,float,false> >::value ));
  VERIFY((  internal::has_unary_operator<internal::linspaced_op<float,float,false> >::value ));
  VERIFY(( !internal::has_binary_operator<internal::linspaced_op<float,float,false> >::value ));
  VERIFY((  internal::functor_has_linear_access<internal::linspaced_op<float,float,false> >::ret ));

  // Regression unit test for a weird MSVC bug.
  // Search "nullary_wrapper_workaround_msvc" in CoreEvaluators.h for the details.
  // See also traits<Ref>::match.
  {
    MatrixXf A = MatrixXf::Random(3,3);
    Ref<const MatrixXf> R = 2.0*A;
    VERIFY_IS_APPROX(R, A+A);

    Ref<const MatrixXf> R1 = MatrixXf::Random(3,3)+A;

    VectorXi V = VectorXi::Random(3);
    Ref<const VectorXi> R2 = VectorXi::LinSpaced(3,1,3)+V;
    VERIFY_IS_APPROX(R2, V+Vector3i(1,2,3));

    VERIFY((  internal::has_nullary_operator<internal::scalar_constant_op<float> >::value ));
    VERIFY(( !internal::has_unary_operator<internal::scalar_constant_op<float> >::value ));
    VERIFY(( !internal::has_binary_operator<internal::scalar_constant_op<float> >::value ));
    VERIFY((  internal::functor_has_linear_access<internal::scalar_constant_op<float> >::ret ));

    VERIFY(( !internal::has_nullary_operator<internal::linspaced_op<int,int,false> >::value ));
    VERIFY((  internal::has_unary_operator<internal::linspaced_op<int,int,false> >::value ));
    VERIFY(( !internal::has_binary_operator<internal::linspaced_op<int,int,false> >::value ));
    VERIFY((  internal::functor_has_linear_access<internal::linspaced_op<int,int,false> >::ret ));
  }
#endif
}
