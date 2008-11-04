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
   \file test_qm_pk_3D.cpp
   \author Gilles Steiner <gilles.steiner@epfl.ch>
   \date 2005-12-22
 */

/** The purpose of this file is to test the quadrature methods
    on the tetraedra.
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
std::ofstream PK_x_log("PK_x_log_3D.log");

std::ofstream ost_x("Monom_x_3D.txt");

//@}


double pow(double x, int i)
{
  double res = 1.0;
  for(int l=1; l <= i; ++l)
    res *= x;
  return res;

}

template<Life::uint16_type N, typename T>
void PK_Monom_N_opt()
{
  using namespace Life;

  typedef T value_type;

  GaussLobatto<Simplex<3,1>, N, value_type> im;

  ublas::vector<value_type> x_i(ublas::row(im.points(),0));

  const value_type tol = value_type(2.5) * type_traits<value_type>::epsilon();

  int Q = (N+3)/2+1;

  value_type error_x = 0.0;
  value_type res_x;
  int i_x=1;
  value_type sum_x=0.0;

  ost_x << "Number of Points on the tetraedra : " << Q << "^3 = "<< im.nPoints() <<  std::endl;

  do{
    if (i_x%2 == 0)
      res_x = value_type(1.0)/value_type(i_x+3.0) + value_type(1.0)/value_type(i_x+1.0);
    else
      res_x = -value_type(2.0)/value_type(i_x+2.0);

    sum_x = 0.0;

    for(uint16_type l=0;l< x_i.size();++l)
      sum_x += pow(x_i(l),i_x)*im.weight(l);

    error_x = math::abs(sum_x - res_x);

    //  std::cout << "i = " << i_x <<"\nError = "<< error_x << "\nRes = "<< res_x << " - sum_x = "<< sum_x <<"\n";
    ost_x << "i = " << i_x <<" Error = "<< error_x << "\n";
    ++i_x;
  }while(error_x <= tol);

  if ( i_x-2 < 2*Q-3  )
    PK_x_log << "Q = " << Q << " ; i = " << i_x << " ; Error = " << error_x << std::endl;

   BOOST_CHECK( (i_x-2 >= 2*Q-3));
}

template<Life::int16_type P>
void add_tests( test_suite* test )
{
  using namespace Life;
  using namespace std;

  test->add( BOOST_TEST_CASE( ( PK_Monom_N_opt<P,double> ) ) );
#if defined( LIFE_HAVE_QD_REAL)
  //test->add( BOOST_TEST_CASE( ( PK_Monom_N_opt<P,dd_real> ) ) );
  //test->add( BOOST_TEST_CASE( ( PK_Monom_N_opt<P,qd_real> ) ) );
#endif

  add_tests<P+2>(test);
}

template<>
void add_tests<40>( test_suite* test )
{
  using namespace Life;
  using namespace std;

  test->add( BOOST_TEST_CASE( ( PK_Monom_N_opt<40,double> ) ) );

#if defined( LIFE_HAVE_QD_REAL)
  //test->add( BOOST_TEST_CASE( ( PK_Monom_N_opt<80,dd_real> ) ) );
  //test->add( BOOST_TEST_CASE( ( PK_Monom_N_opt<80,qd_real> ) ) );
#endif
}

test_suite*
init_unit_test_suite( int /*argc*/, char** /*argv*/ )
{

    test_suite* test = BOOST_TEST_SUITE( "Integration methods on simplicies test suite" );

    PK_x_log << "This file attempt to save the errors encountered while numerically integrating $x^i$ on the tetraedra !" << "\n\n";
#if defined( LIFE_HAVE_QD_REAL)
    unsigned int old_cw;
    fpu_fix_start(&old_cw);
#endif
    add_tests<2>( test );

    return test;
}
