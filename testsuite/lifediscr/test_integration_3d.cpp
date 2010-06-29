/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2006-08-25

  Copyright (C) 2006 EPFL
  Copyright (C) 2007,2008 Université Joseph Fourier (Grenoble I)

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
   \file test_integration.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2006-08-25
 */
#define USE_BOOST_TEST 1
// Boost.Test

// make sure that the init_unit_test function is defined by UTF
#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE 3D integration testsuite
// disable the main function creation, use our own
#define BOOST_TEST_NO_MAIN


#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;
#include <boost/test/floating_point_comparison.hpp>

#include <life/options.hpp>
#include <life/lifecore/environment.hpp>
#include <life/lifemesh/geoentity.hpp>
#include <life/lifemesh/refentity.hpp>
#include <life/lifediscr/functionspace.hpp>
#include <life/lifediscr/mesh.hpp>
#include <life/lifemesh/filters.hpp>
#include <life/lifepoly/im.hpp>
#include <life/lifealg/matrixgmm.hpp>
#include <life/lifefilters/gmshtensorizeddomain.hpp>
#include <life/lifefilters/gmsh.hpp>
#include <life/lifefilters/gmsh.hpp>

#include <life/lifevf/vf.hpp>

const double DEFAULT_MESH_SIZE=0.1;

namespace Life
{
template<typename T, int Dim, int Order = 1>
struct imesh
{
    typedef Mesh<Simplex<Dim, Order>, T > type;
    typedef boost::shared_ptr<type> ptrtype;
};

template<typename T, int Dim>
typename imesh<T, Dim>::ptrtype
createMesh( double hsize )
{
    double meshSize = hsize;
    typedef typename imesh<T,Dim>::type mesh_type;
    typename imesh<T,Dim>::ptrtype mesh( new typename imesh<T,Dim>::type );

    GmshTensorizedDomain<Dim,1,Dim,Simplex> td;
    td.setCharacteristicLength( meshSize );
    //td.setX( std::make_pair( -1, 1 ) );
    ImporterGmsh<typename imesh<T,Dim>::type> import( td.generate( ( boost::format( "hypercube-%1%D" ) % Dim ).str() ) );
    mesh->accept( import );
    mesh->components().set( MESH_RENUMBER | MESH_UPDATE_FACES | MESH_UPDATE_EDGES | MESH_PARTITION );
    mesh->updateForUse();

    return mesh;
}
}
template<typename value_type = double, int Dim=2>
struct test_integration_internal_faces_v
{
    test_integration_internal_faces_v( double meshSize_=DEFAULT_MESH_SIZE ): meshSize(meshSize_), mesh( Life::createMesh<value_type,Dim>( meshSize ) )
    {}

    void operator()()
    {
        using namespace Life;
        using namespace Life::vf;


        typedef typename imesh<value_type,Dim>::type mesh_type;
        typedef typename imesh<value_type,Dim>::ptrtype mesh_ptrtype;
        typename imesh<value_type,Dim>::ptrtype mesh( createMesh<value_type,Dim>( meshSize ) );

        const value_type eps = 1000*Life::type_traits<value_type>::epsilon();

        // int ([-1,1],[-1,x]) 1 dx
        value_type meas = integrate( elements(mesh), cst(1.) ).evaluate()( 0, 0 );
        value_type v0 = integrate( elements(mesh), vf::min(constant(1.0),constant(2.0)) ).evaluate()( 0, 0 );

#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v0, meas, eps );
#else
        LIFE_ASSERT( math::abs( v0-1.0) < eps )( v0 )( math::abs( v0-1.0) )( eps ).warn ( "v0 != 1" );
#endif /* USE_BOOST_TEST */

        value_type v1 = integrate( internalfaces(mesh), jumpv(trans(2*P())) ).evaluate()( 0, 0 );

#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v1, eps );
#else
        LIFE_ASSERT( math::abs( v1-0.0) < eps )( v1 )( math::abs( v1-0.0) )( eps ).warn ( "v1 != 0" );
#endif /* USE_BOOST_TEST */

        value_type v2 = integrate( internalfaces(mesh),
                                   leftfacev(vf::sqrt(trans(P())*P()))-rightfacev(vf::sqrt(trans(P())*P())) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v2, eps );
#else
        LIFE_ASSERT( math::abs( v2-0.0) < eps )( v2 )( math::abs( v2-0.0) )( eps ).warn ( "v2 != 0" );
#endif /* USE_BOOST_TEST */

        value_type v3 = integrate( internalfaces(mesh),
                                   leftfacev(vf::sqrt(trans(P()*P())*(P()*P())))-rightfacev(vf::sqrt(trans(P()*P())*(P()*P()))) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v3, eps );
#else
        LIFE_ASSERT( math::abs( v3-0.0) < eps )( v3 )( math::abs( v3-0.0) )( eps ).warn ( "v3 != 0" );
#endif /* USE_BOOST_TEST */

        typedef FunctionSpace<mesh_type, fusion::vector<Lagrange<3, Scalar> >, double> space_type;
        typedef boost::shared_ptr<space_type> space_ptrtype;
        space_ptrtype Xh = space_type::New( mesh );
        typedef typename space_type::element_type element_type;
        element_type u( Xh, "u" );
        //auto u_exact = Px()+Py()+Pz();
        //auto u_exact = Px()*Px()+Py()*Py()+Pz()*Pz();
        auto u_exact = Px()*Px()*Pz()+Py()*Py()*Px()+Pz()*Pz()*Py();
        //auto u_exact = Px();
        u = vf::project( Xh, elements( mesh ), u_exact );

        double avgv = integrate( internalfaces(mesh), averagev(idv(u)-(u_exact)) ).evaluate()( 0, 0 );
        BOOST_TEST_MESSAGE( "avg(v-uexact)=" << avgv << "\n" );
        double avgv1 = integrate( internalfaces(mesh), .5*(leftfacev(idv(u)-(u_exact))+rightfacev(idv(u)-(u_exact))) ).evaluate()( 0, 0 );
        BOOST_TEST_MESSAGE( ".5*(left(v-uexact)+right(v-uexact))=" << avgv1 << "\n" );
        double avgv2 = integrate( internalfaces(mesh), .5*(leftfacev(idv(u))-leftfacev(u_exact)+rightfacev(idv(u))-rightfacev(u_exact))).evaluate()( 0, 0 );
        BOOST_CHECK_SMALL( avgv, eps );
        BOOST_CHECK_SMALL( avgv1, eps );
        BOOST_CHECK_SMALL( avgv2, eps );


        double leftv = integrate( internalfaces(mesh), leftfacev(idv(u)-(u_exact)) ).evaluate()( 0, 0);
        double rightv = integrate( internalfaces(mesh), rightfacev(idv(u)-(u_exact))).evaluate()( 0, 0);
        BOOST_TEST_MESSAGE( "leftv=" << leftv << "\n" );
        BOOST_CHECK_SMALL( leftv, eps );
        BOOST_TEST_MESSAGE( "rightv=" << rightv << "\n" );
        BOOST_CHECK_SMALL( rightv, eps );


        double n_jumpun = integrate( internalfaces(mesh), trans(leftfacev(N()))*jumpv(idv(u)) ).evaluate()( 0, 0 );
        BOOST_CHECK_SMALL( n_jumpun, eps );
        double n_leftun_rightun = integrate( internalfaces(mesh), trans(leftfacev(N()))*(leftfacev(idv(u)*N())+rightfacev(idv(u)*N())) ).evaluate()( 0, 0 );
        BOOST_CHECK_SMALL( n_leftun_rightun, eps );

        BOOST_CHECK_CLOSE( n_leftun_rightun, n_jumpun, eps );

        double left_n = integrate( internalfaces(mesh), leftfacev(N())).evaluate()( 0, 0 );
        double right_n = integrate( internalfaces(mesh), rightfacev(N())).evaluate()( 0, 0);
        BOOST_CHECK_CLOSE( left_n, -right_n, eps );

        u = vf::project( Xh, elements( mesh ), cst(1.));
        double leftv_1 = integrate( internalfaces(mesh), leftfacev(idv(u))).evaluate()( 0, 0);
        double rightv_1 = integrate( internalfaces(mesh), rightfacev(idv(u))).evaluate()( 0, 0);
        double sumv_1 = integrate( internalfaces(mesh),  rightfacev(idv(u))+leftfacev(idv(u))).evaluate()( 0, 0);
        double avgv_1 = integrate( internalfaces(mesh),  averagev(idv(u))).evaluate()( 0, 0);
        BOOST_CHECK_CLOSE( leftv_1, rightv_1, eps );
        BOOST_CHECK_CLOSE( leftv_1+rightv_1, sumv_1, eps );
        BOOST_CHECK_CLOSE( 2*avgv_1, sumv_1, eps );

    }
    double meshSize;
    typename Life::imesh<value_type,Dim>::ptrtype mesh;
};

template<typename value_type = double, int Dim=2>
struct test_integration_internal_faces_lf
{
    test_integration_internal_faces_lf( double meshSize_=DEFAULT_MESH_SIZE ): meshSize(meshSize_), mesh( Life::createMesh<value_type,Dim>( meshSize ) )
    {}

    void operator()()
    {
        using namespace Life;
        using namespace Life::vf;


        typedef typename imesh<value_type,Dim>::type mesh_type;
        typedef typename imesh<value_type,Dim>::ptrtype mesh_ptrtype;
        typename imesh<value_type,Dim>::ptrtype mesh( createMesh<value_type,Dim>( meshSize ) );

        const value_type eps = 1000*Life::type_traits<value_type>::epsilon();


        typedef FunctionSpace<mesh_type, fusion::vector<Lagrange<3, Scalar> >, double> space_type;
        typedef boost::shared_ptr<space_type> space_ptrtype;
        space_ptrtype Xh = space_type::New( mesh );
        typedef typename space_type::element_type element_type;
        element_type u( Xh, "u" );
        //auto u_exact = Px()+Py()+Pz();
        //auto u_exact = Px()*Px()+Py()*Py()+Pz()*Pz();
        auto u_exact = Px()*Px()*Pz()+Py()*Py()*Px()+Pz()*Pz()*Py();
        //auto u_exact = Px();
        u = vf::project( Xh, elements( mesh ), u_exact );


        boost::shared_ptr<VectorUblas<double> > F( new VectorUblas<double>( u.size() ) );
        std::fill( F->begin(), F->end(), (double)0 );
        form1( Xh, F ) = integrate( internalfaces(mesh),
                                    trans(leftface(id(u)*N())+(rightface(id(u)*N())))*(leftfacev(N()))
            );
        F->close();
        double jumpu_F = inner_product( u, *F );
        BOOST_CHECK_SMALL( jumpu_F, eps );


        form1( Xh, F, _init=true ) = integrate( internalfaces(mesh), leftface(id(u))-(rightface(id(u))) );
        double u_left_right_F = inner_product( u, *F );
        BOOST_CHECK_SMALL( u_left_right_F, eps );

        form1( Xh, F, _init=true ) = integrate( internalfaces(mesh),
                                    (jump(grad(u)))
                                    );
        double jump_gradu_F = inner_product( u, *F );
        BOOST_TEST_MESSAGE ( "jump(grad(u) u^T F = " << jump_gradu_F << "\n" );
        BOOST_CHECK_SMALL( jump_gradu_F, eps );

        form1( Xh, F, _init=true ) = integrate( internalfaces(mesh), leftface(grad(u)*N()));
        double left_gradu_n = inner_product( u, *F );
        form1( Xh, F, _init=true ) = integrate( internalfaces(mesh), rightface(grad(u)*N()));
        double right_gradu_n = inner_product( u, *F );
        BOOST_TEST_MESSAGE(  "jump(left(grad(u)*N)) u^T F = " << left_gradu_n << "\n" );
        BOOST_TEST_MESSAGE(  "jump(right(grad(u)*N)) u^T F = " << right_gradu_n << "\n" );
        BOOST_CHECK_CLOSE( left_gradu_n, -right_gradu_n, eps*10 );

        u = vf::project( Xh, elements( mesh ), cst(1.));
        form1( Xh, F, _init=true ) = integrate( internalfaces(mesh), leftface(id(u)));
        double left_1 = inner_product( u, *F );
        BOOST_TEST_MESSAGE(  "left(id(u)) u^T F = " << left_1 << "\n" );
        form1( Xh, F, _init=true ) = integrate( internalfaces(mesh), rightface(id(u)));
        double right_1 = inner_product( u, *F );
        BOOST_TEST_MESSAGE(  "right(id(u)) u^T F = " << right_1 << "\n" );
        BOOST_CHECK_CLOSE( left_1, right_1, eps );

        form1( Xh, F, _init=true ) = integrate( internalfaces(mesh), trans(N())*jump(id(u)));
        BOOST_CHECK_SMALL( inner_product( u, *F ), eps );

        form1( Xh, F, _init=true ) = integrate( internalfaces(mesh), jump(grad(u)));
        BOOST_CHECK_SMALL( inner_product( u, *F ), eps );

    }
    double meshSize;
    typename Life::imesh<value_type,Dim>::ptrtype mesh;
};

inline
Life::po::options_description
makeOptions()
{
    Life::po::options_description integrationoptions("Test Integration 2D/3D options");
    integrationoptions.add_options()
        ("hsize", Life::po::value<double>()->default_value( 3 ), "h value")
        ;
    return integrationoptions.add( Life::life_options() );
}

inline
Life::AboutData
makeAbout()
{
    Life::AboutData about( "test_integration" ,
                           "test_integration" ,
                            "0.1",
                           "2D/3D integration tests",
                           Life::AboutData::License_GPL,
                           "Copyright (C) 2006-2010 Université Joseph Fourier (Grenoble I)");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}

BOOST_AUTO_TEST_CASE( test_integration_internal_faces_v_double_2 ) { test_integration_internal_faces_v<double,2> t2( 1 ); t2(); }
BOOST_AUTO_TEST_CASE( test_integration_internal_faces_v_double_3 ) { test_integration_internal_faces_v<double,3> t2( 1 ); t2(); }
BOOST_AUTO_TEST_CASE( test_integration_internal_faces_lf_double_2 ) { test_integration_internal_faces_lf<double,2> t2( 1 ); t2(); }
BOOST_AUTO_TEST_CASE( test_integration_internal_faces_lf_double_3 ) { test_integration_internal_faces_lf<double,3> t2( 1 ); t2(); }

#if 0
#if defined(USE_BOOST_TEST)
boost::shared_ptr<Life::Application> mpi;
test_suite*
init_unit_test_suite( int argc, char** argv )
{
    mpi = boost::shared_ptr<Life::Application>( new Life::Application( argc, argv, makeAbout(), makeOptions() ) );
    Life::Assert::setLog( "test_integration.assert");
    test_suite* test = BOOST_TEST_SUITE( "2D Generic finite element solver test suite" );

    test->add( BOOST_TEST_CASE( ( test_integration_internal_faces<double>( mpi->vm()["hsize"].as<double>() ) ) ) );

    return test;
}
#else
int
main( int argc, char** argv )
{
    Life::Application mpi( argc, argv, makeAbout(), makeOptions() );
    Life::Assert::setLog( "test_integration_3d.assert");

    test_integration_internal_faces<double> c ( mpi.vm()["hsize"].as<double>() );
    c();
}
#endif /* USE_BOOST_TEST */
#endif


int BOOST_TEST_CALL_DECL
main( int argc, char* argv[] )
{
    Life::Environment env( argc, argv );
    int ret = ::boost::unit_test::unit_test_main( &init_unit_test, argc, argv );

    return ret;
}
