/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Stephane Veys <stephane.veys@imag.fr>
             Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-12-28

  Copyright (C) 2011-2014 Feel++ Consortium

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
//#define USE_BOOST_TEST 1

#define BOOST_TEST_MODULE test_bdf2
#include <fmt/core.h>
#include <fmt/ostream.h>
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

/** use Feel namespace */
using namespace Feel;
//using Feel::project;

template<int Dim>
class Test:
    public Simget
{
public :

    template<typename T>
    void checkBDFTimeLoop( T& ts )
    {
        for( int o = 1; o < 5; ++o)
        {
            ts->setOrder( o );
            ts->setTimeInitial( 0 );
            ts->setTimeStep( 0.1 );
            fmt::print("times: {}", ts->times() );
            ts->setTimeFinal( 1 );
            for( double t = ts->start(); !ts->isFinished(); t=ts->next() )
            {
                fmt::print("===========================================================\n");
                fmt::print("t={} iteration:{} times={}\n",t, ts->iteration(), ts->times());
                //fmt::print("bdf: {}\n", ts->backwardDifferences(ts->times().tail(3),2) );
            }
        }
    }
    template<typename T>
    void checkBackwardDifferences( T const& mybdf )
    {
        eigen_vector_x_type<double> Ts(5);
        Ts << .1,.2,.3,.4,.5;
        double dt=0.1;
        auto [alpha, d, eta]  = mybdf->backwardDifferences(Ts,1);
        BOOST_TEST_MESSAGE( fmt::format("backward differences 1: {}\n", d  ) );
        BOOST_TEST_MESSAGE( fmt::format("               alpha 1: {}\n", alpha  ) );
        BOOST_TEST_MESSAGE( fmt::format("                 eta 1: {}\n", eta  ) );
        BOOST_CHECK_CLOSE(eta,dt,1e-12);
        BOOST_CHECK_CLOSE(alpha(4),1/dt,1e-12);
        BOOST_CHECK_CLOSE(alpha(3),-1/dt,1e-12);
        BOOST_CHECK_SMALL(alpha(2),1e-12);
        BOOST_CHECK_SMALL(alpha(1),1e-12);
        BOOST_CHECK_SMALL(alpha(0),1e-12);

        auto [alpha2, d2, eta2]  = mybdf->backwardDifferences(Ts,2);
        //BOOST_TEST_MESSAGE( fmt::print("backward differences: {}\n", d  ) );
        //BOOST_TEST_MESSAGE( fmt::print("alpha: {}\n", alpha  ) );
        //BOOST_TEST_MESSAGE( fmt::print("eta: {}\n", eta  ) );
        BOOST_CHECK_CLOSE(eta2,2./(3.*100),1e-12);
        BOOST_CHECK_CLOSE(alpha2(4),3./(2*dt),1e-12);
        BOOST_CHECK_CLOSE(alpha2(3),-2./dt,1e-12);
        BOOST_CHECK_CLOSE(alpha2(2),1./(2*dt),1e-12);
        BOOST_CHECK_SMALL(alpha2(1),1e-12);
        BOOST_CHECK_SMALL(alpha2(0),1e-12);


        auto [alpha3, d3, eta3]  = mybdf->backwardDifferences(Ts,3);
        //BOOST_TEST_MESSAGE( fmt::print("backward differences: {}\n", d  ) );
        //BOOST_TEST_MESSAGE( fmt::print("alpha: {}\n", alpha  ) );
        //BOOST_TEST_MESSAGE( fmt::print("eta: {}\n", eta  ) );
        BOOST_CHECK_CLOSE(eta3,0.001090909090909091,1e-12);
        BOOST_CHECK_CLOSE(alpha3(4),11./(6*dt),1e-12);
        BOOST_CHECK_CLOSE(alpha3(3),-3./dt,1e-12);
        BOOST_CHECK_CLOSE(alpha3(2),3./(2*dt),1e-12);
        BOOST_CHECK_CLOSE(alpha3(1),-1./(3*dt),1e-12);
        BOOST_CHECK_SMALL(alpha3(0),1e-12);

        auto [alpha4, d4, eta4]  = mybdf->backwardDifferences(Ts,4);
        //BOOST_TEST_MESSAGE( fmt::print("backward differences: {}\n", d  ) );
        //BOOST_TEST_MESSAGE( fmt::print("alpha: {}\n", alpha  ) );
        //BOOST_TEST_MESSAGE( fmt::print("eta: {}\n", eta  ) );
        BOOST_CHECK_CLOSE(eta4,0.00028799999999999995,1e-12);
        BOOST_CHECK_CLOSE(alpha4(4),25./(12*dt),1e-12);
        BOOST_CHECK_CLOSE(alpha4(3),-4./dt,1e-12);
        BOOST_CHECK_CLOSE(alpha4(2),3./(dt),1e-12);
        BOOST_CHECK_CLOSE(alpha4(1),-4./(3*dt),1e-12);
        BOOST_CHECK_CLOSE(alpha4(0), 1./(4.*dt),1e-12);

    }

    void run() override
    {
        auto mesh = createGMSHMesh( _mesh=new Mesh<Simplex<Dim,1>>,
                                    _desc=domain( _name=( boost::format( "%1%-%2%" ) % soption(_name="gmsh.domain.shape") % Dim ).str() ,
                                                  _dim=Dim ) );

        auto Xh = Pch<2>( mesh );
        auto u = Xh->element();
        auto ue = Xh->element();
        auto v = Xh->element();
        auto solution = Xh->element();
        auto mybdf = bdf( _space=Xh, _name="mybdf" );

        checkBackwardDifferences( mybdf );

        checkBDFTimeLoop( mybdf );
    }
};




FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( bdf2 )

BOOST_AUTO_TEST_CASE( test_1 )
{
    Test<2> test;
    test.run();
}

BOOST_AUTO_TEST_SUITE_END()
