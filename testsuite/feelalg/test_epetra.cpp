/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s):
  Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  Klaus Sapelza <klaus.sapelza@epfl.ch>
       Date: 2006-09-14

  Copyright (C) 2006 EPFL

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
   \file test_epetra.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Klaus.Sapelza <klaus.sapelza@epfl.ch>
   \date 2006-09-14
 */
#define FEELPP_HAS_BOOST_TEST 1
#include <cstdlib>
#include <cassert>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/traits.hpp>


#include <feel/feelconfig.h>

#if defined(FEELPP_HAS_BOOST_TEST)
// Boost.Test
#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/bind.hpp>

using boost::unit_test_framework::test_suite;


#if defined(FEELPP_HAS_TRILINOS_EPETRA) && defined(FEELPP_HAS_PARMETIS_H)
#include <feel/feelcore/application.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/feelalg/vectorepetra.hpp>
#include <feel/feelalg/matrixepetra.hpp>


namespace Feel
{
AboutData
makeAbout()
{
    AboutData about( "test_epetra" ,
                     "test_epetra" ,
                     "0.1",
                     "test_epetra",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2006 EPFL" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Klaus Sapelza", "developer", "klaus.sapelza@epfl.ch", "" );
    return about;

}

class MyEpetraApp : public Application
{
    typedef Application super;
public:

    MyEpetraApp( int argc, char**argv, AboutData const& about ) : super( argc, argv, about ) {}

    void run()
    {

        std::cout << " \n\nTesting Epetra Vector Interface\n\n" << std::endl;

        Epetra_MpiComm Comm( Application::COMM_WORLD );
        Epetra_Map Map4Vector( 4,0,Comm );



        VectorEpetra<double> x( Map4Vector );
        VectorEpetra<double> y( Map4Vector );
        VectorEpetra<double> z( Map4Vector );
        VectorEpetra<double> w( Map4Vector );

        int MyLength = x.localSize();
        //int MyMaxLID = Map.MaxLID();

        //cout << Comm << endl;
        //cout << map << endl;

        //    std::cout << "\nMyLength = " << MyLength <<  std::endl;


        x.set( 0,1.0 );
        x.set( 1,2.0 );
        x.set( 2,3.0 );
        x.set( 3,4.0 );
        //x.set(1.0);
        y.set( 2.0 );

        z=x;
        z+=y;

        x.close();
        y.close();
        z.close();
        w.close();

        x.printMatlab( "test_x.m" );
        w.printMatlab( "test_w.m" );

        for ( int i=0 ; i<MyLength ; ++i )
        {
            std::cout << "x(" << i << ") = " <<  x( i ) <<  std::endl;
        }

        std::cout <<  std::endl;

        for ( int i=0 ; i<MyLength ; ++i )
        {
            std::cout << "y(" << i << ") = " <<  y( i ) <<  std::endl;
        }

        std::cout << std::endl;

        for ( int i=0 ; i<MyLength ; ++i )
        {
            std::cout << "z(" << i << ") = " <<  z( i ) <<  std::endl;
        }

        std::cout <<  std::endl;



        std::cout << " \n\nTesting Epetra Matrix Interface\n\n" << std::endl;

        int NumRowElements = 4;
        int NumColElements = 6;

        //create a map to create a distributed Matrix

        Epetra_Map RowMap( NumRowElements,0,Comm );
        Epetra_Map ColMap( NumColElements,NumColElements,0,Comm );

        Epetra_Map RowMap1( 4,2,0,Comm );
        Epetra_Map ColMap1( 4,4,0,Comm );

        //std::cout << "Map: " <<  Map4Matrix <<  std::endl;

        Feel::MatrixEpetra A( RowMap, 3 );
        Feel::MatrixEpetra B( 6,4,3,2,2,0 );
        Feel::MatrixEpetra C( RowMap,ColMap );
        Feel::MatrixEpetra D( RowMap1,ColMap1 );



        switch ( Application::processId() )
        {
        case 0:
            B.set( 0,0,10.0 );
            B.set( 0,1,12.0 );
            B.set( 0,2,13.0 );
            B.set( 1,0,14.0 );
            B.set( 1,1,15.0 );
            B.set( 1,2,16.0 );

            C.set( 0,0,20.0 );
            C.set( 0,1,22.0 );
            C.set( 0,2,23.0 );
            C.set( 1,0,24.0 );
            C.set( 1,1,25.0 );
            C.set( 1,2,26.0 );

            D.set( 0,0,1.0 );
            D.set( 1,1,2.0 );
            D.set( 1,2,3.0 );
            D.set( 1,3,4.0 );
            break;

        case 1:
            B.set( 1,1,17.0 );
            B.set( 2,0,17.0 );
            B.set( 2,1,18.0 );
            B.set( 5,3,19.0 );

            C.set( 1,2,-20.0 );
            C.set( 2,0,27.0 );
            C.set( 2,1,28.0 );
            C.set( 3,5,29.0 );

            D.set( 0,0,-2.0 );
            D.set( 2,0,8.0 );
            D.set( 2,2,5.0 );
            D.set( 3,2,7.0 );
            break;
        }


        A.set( 0,0,0.0 );
        A.set( 0,1,2.0 );
        A.set( 0,2,3.0 );
        A.set( 1,0,4.0 );
        A.set( 1,1,5.0 );
        A.set( 1,2,6.0 );
        A.set( 2,0,7.0 );
        A.set( 2,1,8.0 );
        A.set( 2,2,9.0 );

        A.add( 1,1,31.0 );





        if ( A.EpetraIndicesAreLocal() == 0 )
        {
            std::cout << "Indices HAVE NOT been transformed to Local Indexspace, ie. A.close HAS NOT been called! " << std::endl;
        }

        else std::cout << "Indices HAVE been transformed to Local Indexspace, ie. A.close HAS been called!  " << std::endl;


        int spalten = A.size1();
        int zeilen = A.size2();
        int rowstart = A.rowStart();
        int rowstop = A.rowStop();

        std::cout << "Nb of rows: " <<  spalten <<  std::endl;
        std::cout << "Nb of columns: " <<  zeilen <<  std::endl;

        std::cout << "First local row index: " <<  rowstart <<  std::endl;
        std::cout << "First local row index: " <<  rowstop <<  std::endl;

        if ( A.closed() == 0 )
        {
            std::cout << "Matrix NOT closed! " << std::endl;
        }

        else std::cout << "Matrix IS closed! " << std::endl;


        //     for( int i=0 ; i<NumGlobalElements ; ++i )
        //     {
        //       for( int j=0 ; j<NumGlobalElements ; ++j )
        //       {
        //         std::cout << "A(" << i << "," << j << ") = " << A(i,j) << " ";
        //       }
        //       std::cout << std::endl;
        //     }

        //

        A.close();
        B.close();
        C.close();
        D.close();

        // Print out Matrices after closing them: see components assembled together!
        //A.printKonsole();
        //B.printKonsole();
        //C.printKonsole();
        //D.printKonsole();

        // Works only if matrices have been closed!
        A.printMatlab( "test_A.m" );
        B.printMatlab( "test_B.m" );
        C.printMatlab( "test_C.m" );
        D.printMatlab( "test_D.m" );



        if ( A.closed() == 0 )
        {
            std::cout << "Matrix NOT closed, will close NOW! " << std::endl;
            A.close();
        }

        else
        {
            std::cout << "Matrix IS closed! " << std::endl;
        }

        std::cout << "l1 Norm = " <<  A.l1Norm() <<  std::endl;
        std::cout << "linf Norm = " <<  A.linftyNorm() <<  std::endl;

    } //end run()



};
} // Feel
void epetra_manager()
{
}

test_suite*
init_unit_test_suite( int argc, char** argv )
{
    test_suite* test= BOOST_TEST_SUITE( "EPETRA Unit Test" );

    Feel::MyEpetraApp myapp( argc, argv, Feel::makeAbout() );

    myapp.run();

    // this example will pass cause we know ahead of time number of expected failures
    //test->add( BOOST_TEST_CASE( &epetra_manager ), 0 );

    return test;
}
#else
test_suite*
init_unit_test_suite( int /*argc*/, char** /*argv*/ )
{
    test_suite* test= BOOST_TEST_SUITE( "EPETRA Unit Test" );
    return test;
}
#endif /* FEELPP_HAS_TRILINOS_EPETRA */
#else
int main()
{
    return EXIT_SUCCESS;
}
#endif
