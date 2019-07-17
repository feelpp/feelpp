/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

   Author(s): Thomas Lantz
   Date: 2015-04-27

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file test_multiscaleimage.cpp
   \author Thomas Lantz
   \date 2015-04-27
 */

#define USE_BOOST_TEST 1
#if defined(USE_BOOST_TEST)
#define BOOST_TEST_MODULE test_multiscaleimage
#include <feel/feelcore/testsuite.hpp>
#endif

#include <feel/feeldiscr/multiscaleimage.hpp>

namespace Feel
{

inline
AboutData
makeAbout()
{
    AboutData about( "test_multiscaleimage" ,
                     "test_multiscaleimage" ,
                     "8.9e-3",
                     "Test multiscale image acces",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2014 Feel++ Consortium" );

    about.addAuthor( "Thomas Lantz", "student", "", "" );
    return about;
}


template < typename T  >
using holo3_image = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> ;

class Test_MultiScale
{
    public :


    Test_MultiScale(int nx,int ny)
    {
        ima = holo3_image<float>(ny,nx);
        for (int i=0;i<nx;i++)
        {
            for (int j=0;j<ny;j++)
            {
                ima(j,i)=nx*j+i;
            }
        }
    }
    
    
    void run (double x, double y, double l) 
    {
      
        if ( x*l >= ima.cols()*8.9e-3 || y*l >= ima.rows()*8.9e-3 )
        {
            std::cout << " Coord out of the image " << std::endl;
        }
            else
            {
                ublas::vector<double> v (2);
                v(0)=x;
                v(1)=y;
                MultiScaleImage<float> m(ima,l);
                int tmp=m(v,v);
                std::cout << "coord x : " << x << " / coord y : " << y << " / Res :" << tmp << std::endl;
            }
     
    }

/*
    void runLevel (double x, double y, double level) 
    {
      
        if ( x*level >= ima.cols()*8.9e-3 || y*level >= ima.rows()*8.9e-3 )
        {
            std::cout << " Coord out of the image " << std::endl;
        }
            else
            {
                ublas::vector<double> v (2);
                v(0)=x;
                v(1)=y;
                MultiScaleImage m(ima);
                int tmp=m(v,level);
                std::cout << "coord x : " << level*x << " / coord y : " << level*y << " / Res :" << tmp << std::endl;
            }
     
    }
    */
    private :

    holo3_image<float> ima;
    
};

#if defined(USE_BOOST_TEST)
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), feel_options() )
BOOST_AUTO_TEST_SUITE( multiscale_suite )

double ps=8.9e-3;

BOOST_AUTO_TEST_CASE( test_run0 )
{
    Test_MultiScale tms0(12,12);
    tms0.run(2*ps,2*ps,4);
}

BOOST_AUTO_TEST_CASE( test_run1 )
{
    Test_MultiScale tms1(12,12);
    tms1.run(11*ps,4*ps,1);
}

BOOST_AUTO_TEST_CASE( test_run2 )
{
    Test_MultiScale tms2(32,7);
    tms2.run(11*ps,6*ps,2);
}

BOOST_AUTO_TEST_CASE( test_run3 )
{
    Test_MultiScale tms3(6,19);
    tms3.run(2*ps,10*ps,2);
}

/*
this size is really too big for a test
BOOST_AUTO_TEST_CASE( test_run4 )
{
    Test_MultiScale tms4(12345,67890);
    tms4.run(90*ps,1234*ps,4);
}
 */

BOOST_AUTO_TEST_CASE( test_run5 )
{
    Test_MultiScale tms5(5,7);
    tms5.run(0*ps,2*ps,1);
}
/*
BOOST_AUTO_TEST_CASE( test_runL0 )
{
    Test_MultiScale tms0(12,12);
    tms0.runLevel(2*ps,2*ps,4);
}

BOOST_AUTO_TEST_CASE( test_runL1 )
{
    Test_MultiScale tms1(12,12);
    tms1.runLevel(1*ps,2*ps,8);
}

BOOST_AUTO_TEST_CASE( test_runL2 )
{
    Test_MultiScale tms2(32,7);
    tms2.runLevel(11*ps,2*ps,2);
}

BOOST_AUTO_TEST_CASE( test_runL3 )
{
    Test_MultiScale tms3(6,19);
    tms3.runLevel(2*ps,10*ps,1);
}

BOOST_AUTO_TEST_CASE( test_runL4 )
{
    Test_MultiScale tms4(12345,67890);
    tms4.runLevel(90*ps,1234*ps,4);
}

BOOST_AUTO_TEST_CASE( test_runL5 )
{
    Test_MultiScale tms5(5,7);
    tms5.runLevel(1*ps,3*ps,2);
}
*/
BOOST_AUTO_TEST_SUITE_END()
#else
std::cout << "USE_BOOST_TEST non define" << std::endl;
#endif


}
