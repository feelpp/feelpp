/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil;
   c-basic-offset: 4; show-trailing-whitespace: t -*-
vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

This file is part of the Feel library

Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
Date: 2014-01-14

Copyright (C) 2014-2016 Feel++ Consortium

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#define USE_BOOST_TEST 1
#if defined(USE_BOOST_TEST)
#define BOOST_TEST_MODULE test_fit
#include <testsuite/testsuite.hpp>
#endif

#include <feel/feel.hpp>
#include <feel/feelfit/fit.hpp>
#include <feel/feelfit/fitdiff.hpp>

/** use Feel namespace */
using namespace Feel;

inline
  AboutData
makeAbout()
{
  AboutData about( "test_fit" ,
      "test_fit" ,
      "0.1",
      "nD(n=2,3) test testfit",
      Feel::AboutData::License_GPL,
      "Copyright (c) 2014 Feel++ Consortium" );

  about.addAuthor( "Vincent Huber", "developer", "vincent.huber@cemosis.fr", "" );
  return about;

}

double f(double t = 0.) { return 5. * t + sin(t); }

#if defined(USE_BOOST_TEST)
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), feel_options() );
BOOST_AUTO_TEST_SUITE( testfit_suite )

BOOST_AUTO_TEST_CASE( test_main_fit )
{
  Feel::Environment::changeRepository( boost::format( "testsuite/feelfit/%1%/" )
      % Feel::Environment::about().appName()
      );
  auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
  auto Xh = Pch<1>(mesh);
  auto T = Xh->element(); // the temperature (say)
  auto K = Xh->element(); // K(T) - the dependant of the temperature conductivity
  auto Kd= Xh->element(); // K'(T)
  auto T_K = Xh->element(); // T_ = true 
  auto T_Kd= Xh->element(); // 
  auto D_K = Xh->element(); // D_ = difference
  auto D_Kd= Xh->element(); // 
  T.on(_range=elements(mesh), _expr=Px());
  T_K.on(_range=elements(mesh),_expr=(5*Px()+sin(Px())));
  T_Kd.on(_range=elements(mesh),_expr=(5+cos(Px())));

#if 0
  auto e = exporter(_mesh = mesh );
#endif
  // Generates the datafile
  // we assume an unitsquare as mesh
  std::ofstream datafile(Environment::expand(soption("fit.datafile")));
  for(double t = -1; t < 2; t+=0.32)
    datafile << t << "\t" << f(t) << "\n";
  datafile.close();


  for(int i = 0; i < 4; i++)  
  {
    BOOST_TEST_MESSAGE( boost::format("test %1%")% i );
    // evaluate K(T) with the interpolation from the datafile
    K.on(_range=elements(mesh), _expr=fit(idv(T),Environment::expand(soption("fit.datafile")),i));
    Kd.on(_range=elements(mesh), _expr=fitDiff(idv(T),Environment::expand(soption("fit.datafile")), i) );

    D_K.on(_range=elements(mesh),_expr=vf::abs(idv(K)-idv(T_K)));
    D_Kd.on(_range=elements(mesh),_expr=vf::abs(idv(Kd)-idv(T_Kd)));

    auto max_K = D_K.max();
    auto max_Kd= D_Kd.max();
#if 0
    e->step(i)->add("T",T);
    e->step(i)->add("K",K);
    e->step(i)->add("Kd",Kd);
    e->step(i)->add("T_K",T_K);
    e->step(i)->add("T_Kd",T_Kd);
    e->step(i)->add("D_K",D_K);
    e->step(i)->add("D_Kd",D_Kd);
    e->save();
#endif
    /// Note : the error has nothing to do with the mesh size but the step on the datafile
    switch(i)
    {
        case 0: //P0 interpolation
            {
                BOOST_CHECK_SMALL(max_K, 0.95);
                BOOST_CHECK_SMALL(max_Kd, 6.0001); // the derivative is null
                break;
            }
        case 1: // P1 interpolation
            {
                BOOST_CHECK_SMALL(max_K, 0.01);
                BOOST_CHECK_SMALL(max_Kd, 0.15); 
                break;
            }
        case 2: // CSpline interpolation
            {
                BOOST_CHECK_SMALL(max_K, 0.01);
                BOOST_CHECK_SMALL(max_Kd, 0.15); 
                break;
            }
        case 3: // Akima interpolation
            {
                BOOST_CHECK_SMALL(max_K, 0.016);
                BOOST_CHECK_SMALL(max_Kd, 0.03); 
                break;
            }
    }

    BOOST_TEST_MESSAGE( boost::format("test %1% done")% i );
  }
}

BOOST_AUTO_TEST_SUITE_END()
#endif

