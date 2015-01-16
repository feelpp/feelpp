/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
Date: 2007-12-19

Copyright (C) 2007-2012 Université Joseph Fourier (Grenoble I)

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
  \file test_interpolation.cpp
  \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  \date 2007-12-19
  */
#define USE_BOOST_TEST 1

// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE interpolation testsuite
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN

#include <testsuite/testsuite.hpp>

#include <feel/options.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feelmesh/geoentity.hpp>
#include <feel/feelmesh/refentity.hpp>
#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/interpolate.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>

#include <feel/feelfilters/exporterquick.hpp>
#include <feel/feelvf/vf.hpp>

const double DEFAULT_MESH_SIZE=0.5;

using namespace Feel;

template<int Dim, int Order=1, int RDim=Dim>
struct imesh
{
  typedef Mesh<Simplex<Dim,Order,RDim>, double > type;
  typedef boost::shared_ptr<type> ptrtype;
};

template<int Dim, int Order, int RDim>
typename imesh<Dim, Order, RDim>::ptrtype
createMesh( double hsize )
{
    return createGMSHMesh( _mesh=new typename imesh<Dim, Order, RDim>::type,
                           _desc=domain( _name=( boost::format( "%1%-%2%-%3%" )  % "hypercube" % Dim % Order ).str() ,
                                         _addmidpoint=false,
                                         _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES,
                                         _shape="hypercube",
                                         _dim=Dim,
                                         _order=Order,
                                         _h=doption(_name="gmsh.hsize") ) );
    //typedef typename imesh<Dim, Order, RDim>::type meshType;
    //return loadMesh( _mesh=new meshType,
    //                       _desc=domain(_addmidpoint=false,
    //                                     _h=hsize ) );

}

template<int Dim, int Order, int GeoOrder=1>
struct test_interpolation
{
    typedef typename imesh<Dim,GeoOrder>::type mesh_type;
    typedef double value_type;

    test_interpolation( double meshSize_=DEFAULT_MESH_SIZE )
        :
        meshSize( meshSize_ ),
        mesh( createMesh<Dim,GeoOrder,Dim>( meshSize ) )
        {}
    void operator()()
        {
            using namespace Feel;
            using namespace Feel::vf;

            BOOST_MESSAGE( "= geometric order = " << GeoOrder );
            const value_type eps = (GeoOrder == 1 ? meshSize*meshSize : 1000*Feel::type_traits<value_type>::epsilon() );

            auto Xh = Pch<Order>(mesh);
            auto u = Xh->element();

            u = vf::project( _space=Xh, _range=elements( mesh ), _expr=constant( 1.0 ) );

            node_type pt( Dim );
            pt[0] = 0.11;

            if ( Dim >= 2 )
                pt[1] = 0.11;

            if ( Dim >= 3 )
                pt[2] = 0.11;

            double v0 = u( pt )( 0,0,0 );
#if defined(USE_BOOST_TEST)
            BOOST_CHECK_SMALL( v0-1.0, eps );
#else
            FEELPP_ASSERT( math::abs( v0-1.0 ) < eps )( v0 )( math::abs( v0-1.0 ) )( eps ).warn ( "v0 != 1" );
#endif /* USE_BOOST_TEST */

            u = vf::project( _space=Xh, _range=elements( mesh ), _expr=constant( 2.0 ) - Px()*Px()-Py()*Py()-Pz()*Pz() );
            pt[0] = 0.5;

            if ( Dim >= 2 )
                pt[1] = 0.5;

            if ( Dim >= 3 )
                pt[2] = 0.5;

            double v1 = u( pt )( 0,0,0 );
            double v1_ex = 2-Dim*0.5*0.5;
#if defined(USE_BOOST_TEST)
            BOOST_MESSAGE(  "[test_interpolation] v1    = " << v1 );
            BOOST_MESSAGE(  "[test_interpolation] v1_ex = " << v1_ex );
            BOOST_CHECK_SMALL( v1-v1_ex, eps );
#else
            FEELPP_ASSERT( math::abs( v1-v1_ex ) < eps )( v1 )( math::abs( v1-v1_ex ) )( eps ).warn ( "v1 != v0_ex" );
#endif /* USE_BOOST_TEST */
            auto gradient( u.grad( pt ) );
            double g_v1_x = gradient( 0,0,0 );
            double g_v1_ex_x = -2*0.5;
#if defined(USE_BOOST_TEST)
            BOOST_MESSAGE( "[test_interpolation] g_v1_x    = " << g_v1_x );
            BOOST_MESSAGE( "[test_interpolation] g_v1_ex_x = " << g_v1_ex_x );
            if(GeoOrder == 1)
                BOOST_CHECK_SMALL( g_v1_x-g_v1_ex_x, meshSize );
            else
                BOOST_CHECK_SMALL( g_v1_x-g_v1_ex_x, eps );
#else
            FEELPP_ASSERT( math::abs( g_v1_x-g_v1_ex_x ) < eps )( g_v1_x )( math::abs( g_v1_x-g_v1_ex_x ) )( eps ).warn ( "g_v1 != g_v1_ex" );
#endif /* USE_BOOST_TEST */

            if ( Dim >= 2 )
            {
                double g_v1_y = gradient( 0,1,0 );
#if defined(USE_BOOST_TEST)
                if(GeoOrder == 1)
                    BOOST_CHECK_SMALL( g_v1_y-g_v1_ex_x, 2*meshSize );
                else
                    BOOST_CHECK_SMALL( g_v1_y-g_v1_ex_x, eps );
#endif /* USE_BOOST_TEST */
            }

            if ( Dim >= 3 )
            {
                double g_v1_z = gradient( 0,2,0 );
#if defined(USE_BOOST_TEST)
                if(GeoOrder == 1)
                    BOOST_CHECK_SMALL( g_v1_z-g_v1_ex_x, meshSize );
                else
                    BOOST_CHECK_SMALL( g_v1_z-g_v1_ex_x, eps );
#endif /* USE_BOOST_TEST */
            }


        }
    double meshSize;



    typename imesh<Dim,GeoOrder>::ptrtype mesh;

};


inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description integrationoptions( "Test Integration options" );
    integrationoptions.add_options()
        ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "h value" )
        ;
    return integrationoptions.add( Feel::feel_options() );
}

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "test_interpolation" ,
                           "test_interpolation" ,
                           "0.1",
                           "interpolation tests",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2007, 2008 Universite Joseph Fourier (Grenoble I)\n"
                           "Copyright (C) 2013 Universite de Strasbourg" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

#if defined(USE_BOOST_TEST)
#if 0
    boost::shared_ptr<Feel::Application> mpiapp;
    test_suite*
        init_unit_test_suite( int argc, char** argv )
    {
        //boost::mpi::environment( argc, argv );
        mpiapp = boost::shared_ptr<Feel::Application>( new Feel::Application( argc, argv, makeAbout(), makeOptions() ) );
        Feel::Assert::setLog( "test_integration.assert" );
        test_suite* test = BOOST_TEST_SUITE( "Interpolation test suite" );

#if 1
        test->add( BOOST_TEST_CASE( ( test_interpolation<1,2>( mpiapp->vm()["hsize"].as<double>() ) ) ) );
        test->add( BOOST_TEST_CASE( ( test_interpolation<2,2>( mpiapp->vm()["hsize"].as<double>() ) ) ) );
        test->add( BOOST_TEST_CASE( ( test_interpolation<3,2>( mpiapp->vm()["hsize"].as<double>() ) ) ) );

        test->add( BOOST_TEST_CASE( ( test_interpolation<2,2,2>( mpiapp->vm()["hsize"].as<double>() ) ) ) );

        test->add( BOOST_TEST_CASE( ( test_interpolation_op<1,1,1>( mpiapp->vm()["hsize"].as<double>() ) ) ) );
        test->add( BOOST_TEST_CASE( ( test_interpolation_op<1,1,2>( mpiapp->vm()["hsize"].as<double>() ) ) ) );

        test->add( BOOST_TEST_CASE( ( test_lagrange_p1_op<2,1>( mpiapp->vm()["hsize"].as<double>() ) ) ) );
        test->add( BOOST_TEST_CASE( ( test_lagrange_p1_op<2,2>( mpiapp->vm()["hsize"].as<double>() ) ) ) );
        test->add( BOOST_TEST_CASE( ( test_lagrange_p1_op<2,5>( mpiapp->vm()["hsize"].as<double>() ) ) ) );
        test->add( BOOST_TEST_CASE( ( test_lagrange_p1_op<2,3,2>( mpiapp->vm()["hsize"].as<double>() ) ) ) );
        //test->add( BOOST_TEST_CASE( ( test_lagrange_p1_op<3,1>( mpiapp->vm()["hsize"].as<double>() ) ) ) );
        //test->add( BOOST_TEST_CASE( ( test_lagrange_p1_op<3,4>( mpiapp->vm()["hsize"].as<double>() ) ) ) );
#else
        test->add( BOOST_TEST_CASE( ( test_lagrange_p1_op<2,1>( mpiapp->vm()["hsize"].as<double>() ) ) ) );
#endif
        return test;
    }
#else
    FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );

    BOOST_AUTO_TEST_SUITE( interpolation_suite )

        typedef Feel::Application Application_type;
    typedef boost::shared_ptr<Application_type> Application_ptrtype;

    BOOST_AUTO_TEST_CASE( test_interpolation12 )
    {
        BOOST_MESSAGE( "================================================================================" );
        BOOST_MESSAGE( "== test_interpolation<1,2>" );
        test_interpolation<1,2> t;
        t();
    }
    BOOST_AUTO_TEST_CASE( test_interpolation22 )
    {
        BOOST_MESSAGE( "================================================================================" );
        BOOST_MESSAGE( "== test_interpolation<2,2>" );
        test_interpolation<2,2> t;
        t();
    }
    BOOST_AUTO_TEST_CASE( test_interpolation32 )
    {
        BOOST_MESSAGE( "================================================================================" );
        BOOST_MESSAGE( "== test_interpolation<3,2>" );
        test_interpolation<3,2> t;
        t();
    }
    BOOST_AUTO_TEST_CASE( test_interpolation222 )
    {
        BOOST_MESSAGE( "================================================================================" );
        BOOST_MESSAGE( "== test_interpolation<2,2,2>" );
        test_interpolation<2,2,2> t;
        t();
    }
    BOOST_AUTO_TEST_SUITE_END()

#endif // 0
#else
        int
        main( int argc, char** argv )
    {
        Feel::Environment env( _argc=argc, _argv=argv,
                               _desc=makeOptions(),
                               _about=makeAbout() );
        Feel::Application mpiapp;
        Feel::Assert::setLog( "test_interpolation.assert" );

        //test_interpolation<2,1,2> t212( mpiapp.vm()["hsize"].as<double>() );
        //t212();


    }
#endif /* USE_BOOST_TEST */
#if 0
    int BOOST_TEST_CALL_DECL
        main( int argc, char* argv[] )
    {
        Feel::Environment env( argc, argv );
        int ret = ::boost::unit_test::unit_test_main( &init_unit_test, argc, argv );

        return ret;
    }
#endif
