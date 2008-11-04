/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Gilles Steiner <gilles.steiner@epfl.ch>
       Date: 2005-12-22

  Copyright (C) 2005,2006 EPFL

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file test_qm.cpp
   \author Gilles Steiner <gilles.steiner@epfl.ch>
   \date 2005-12-22
 */

/**
    The purpose of this file is to test the quadrature methods
    on tensorised geometries.
**/

// Boost.Test
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;

#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <fstream>
#include <life/lifepoly/im.hpp>


namespace mpl = boost::mpl;


/** Streams **/
//@{
std::ofstream QK_log("QK_error.log");
std::ofstream ost("QK_results.txt");
//@}


double pow(double x, int i)
{
  double res = 1.0;
  for(int l=1; l <= i; ++l)
    res *= x;
  return res;

}

template<Life::uint16_type N, typename T>
void QK_find_N_opt()
{
  using namespace Life;
  typedef T value_type;

  Gauss<Quad, N ,value_type > im;

  ublas::vector<value_type> x_i(ublas::row(im.points(),0));
  ublas::vector<value_type> y_i(ublas::row(im.points(),1));
  const value_type tol = value_type(7.0)*type_traits<value_type>::epsilon();

  /* 1D number of Nodes */
  //@{
  uint16_type Q = (uint16_type)sqrt(im.nPoints());
  //@}

  value_type error = 0.0;
  value_type res;
  uint16_type i=1;
  value_type sum=0.0;

  ost << "Nbre of Points on the Quadrangle : " << Q << "^2 = "<< im.nPoints() <<  std::endl;

  do{
    if (i%2 == 0)
      res = value_type(4.0)/value_type(i+1.0)/value_type(i+1.0);
    else
      res = 0.0;

    sum = 0.0;

    for(uint16_type l=0;l< x_i.size();++l)
      {
        sum += pow(x_i(l),i)*pow(y_i(l),i)*im.weight(l);
      }
    error = math::abs(sum - res);
    ost << "i = " << i <<" Error = "<< error << "\n";
    ++i;
  }while(error <= tol);

  if ( i-2 < 2*Q-1  )
    QK_log << "Q = " << Q << " ; i = " << i << " ; Error" << error << std::endl;

  BOOST_CHECK( i-2 >= 2*Q-1  );
}

template<Life::int16_type P>
void add_tests( test_suite* test )
{
  using namespace Life;
  using namespace std;

  test->add( BOOST_TEST_CASE( ( QK_find_N_opt<P,double> ) ) );

#if defined(LIFE_HAVE_QD_REAL)
  //test->add( BOOST_TEST_CASE( ( QK_find_N_opt<P,dd_real> ) ) );
  //test->add( BOOST_TEST_CASE( ( QK_find_N_opt<P,qd_real> ) ) );
#endif

  add_tests<P+2>(test);
}

template<>
void add_tests<80>( test_suite* test )
{
  using namespace Life;
  using namespace std;

  test->add( BOOST_TEST_CASE( ( QK_find_N_opt<80,double> ) ) );

#if defined(LIFE_HAVE_QD_REAL)
  //test->add( BOOST_TEST_CASE( ( QK_find_N_opt<80,dd_real> ) ) );
  //test->add( BOOST_TEST_CASE( ( QK_find_N_opt<80,qd_real> ) ) );
#endif
}

test_suite*
init_unit_test_suite( int /*argc*/, char** /*argv*/ )
{

    test_suite* test = BOOST_TEST_SUITE( "Quadrature methods on tensorised geometries test suite" );

    QK_log << "This file attempt to save the errors encountered while numerically integrating $x^iy^i$ on the square !" << "\n\n";

    ost << "This file save the details of the execution of 'test_qm_qk' !" << "\n\n";

#if defined(LIFE_HAVE_QD_REAL)
    unsigned int old_cw;
    fpu_fix_start(&old_cw);
#endif
    add_tests<2>( test );

    return test;
}
