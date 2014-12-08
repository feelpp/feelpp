/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-08-25

  Copyright (C) 2006 EPFL
  Copyright (C) 2007-2010 Université Joseph Fourier (Grenoble I)

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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-08-25
 */
#define USE_BOOST_TEST 1
// Boost.Test

// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE 3D integration testsuite
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN

#include <testsuite/testsuite.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/options.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>

#include <feel/feelvf/vf.hpp>

const double DEFAULT_MESH_SIZE=0.1;

namespace Feel
{
template<typename T, int Dim, int Order = 1>
struct imesh
{
    typedef Simplex<Dim, Order> convex_type;
    typedef Mesh<convex_type, T > type;
    typedef boost::shared_ptr<type> ptrtype;
};

template<typename value_type = double, int Dim=2>
struct test_integration_internal_faces_v: public Application
{
    typedef typename imesh<value_type,Dim>::convex_type convex_type;
    typedef typename imesh<value_type,Dim>::type mesh_type;
    typedef typename imesh<value_type,Dim>::ptrtype mesh_ptrtype;
    typedef FunctionSpace<mesh_type, bases<Lagrange<3, Scalar> >, double> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    test_integration_internal_faces_v()
        :
        Application(),
        backend( Backend<double>::build( soption( _name="backend" ) ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( this->vm()["shape"].template as<std::string>() ),
        mesh()
    {
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % Dim ).str() ,
                                             _usenames=true,
                                             _convex=( convex_type::is_hypercube )?"Hypercube":"Simplex",
                                             _shape=shape,
                                             _dim=Dim,
                                             _h=meshSize ) );
    }
    void operator()()
    {
        using namespace Feel::vf;

        const value_type eps = 1000*Feel::type_traits<value_type>::epsilon();

        space_ptrtype Xh = space_type::New( _mesh=mesh,_extended_doftable=std::vector<bool>(1,true) );

        // int ([-1,1],[-1,x]) 1 dx
        value_type meas = integrate( elements( mesh ), cst( 1. ) ).evaluate()( 0, 0 );
        value_type v0 = integrate( elements( mesh ), vf::min( constant( 1.0 ),constant( 2.0 ) ) ).evaluate()( 0, 0 );

#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v0, meas, eps );
#else
        FEELPP_ASSERT( math::abs( v0-1.0 ) < eps )( v0 )( math::abs( v0-1.0 ) )( eps ).warn ( "v0 != 1" );
#endif /* USE_BOOST_TEST */

        value_type v1 = integrate( internalfaces( mesh ), jumpv( trans( 2*P() ) ) ).evaluate()( 0, 0 );
        value_type v1l = integrate( internalfaces( mesh ), leftfacev( trans( 2*P() )*N() ) ).evaluate()( 0, 0 );
        value_type v1r = integrate( internalfaces( mesh ), rightfacev( trans( 2*P() )*N() ) ).evaluate()( 0, 0 );

#if defined(USE_BOOST_TEST)
        BOOST_TEST_MESSAGE( "int jump(2*X^t) =" << v1 << "\n" );
        BOOST_CHECK_SMALL( v1, eps );
        BOOST_TEST_MESSAGE( "int left(2*X^t) =" << v1l << "\n" );
        BOOST_TEST_MESSAGE( "int right(2*X^t) =" << v1r << "\n" );
        BOOST_CHECK_CLOSE( v1l, -v1r, eps );
#else
        FEELPP_ASSERT( math::abs( v1-0.0 ) < eps )( v1 )( math::abs( v1-0.0 ) )( eps ).warn ( "v1 != 0" );
#endif /* USE_BOOST_TEST */

        auto normp = vf::sqrt( trans( P() )*P() );
        auto vnormp = normp*unitX()+normp*unitY()+normp*unitZ();
        value_type v2 = integrate( internalfaces( mesh ), leftfacev( trans( vnormp )*N() )+rightfacev( trans( vnormp )*N() ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_TEST_MESSAGE( "int jump(vnormp) =" << v2 << "\n" );
        BOOST_CHECK_SMALL( v2, eps );
#else
        FEELPP_ASSERT( math::abs( v2-0.0 ) < eps )( v2 )( math::abs( v2-0.0 ) )( eps ).warn ( "v2 != 0" );
#endif /* USE_BOOST_TEST */

        value_type v3 = integrate( internalfaces( mesh ),
                                   leftfacev( vf::sqrt( trans( P()*P() )*( P()*P() ) ) )-rightfacev( vf::sqrt( trans( P()*P() )*( P()*P() ) ) ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v3, eps );
#else
        FEELPP_ASSERT( math::abs( v3-0.0 ) < eps )( v3 )( math::abs( v3-0.0 ) )( eps ).warn ( "v3 != 0" );
#endif /* USE_BOOST_TEST */

        element_type u( Xh, "u" );
        //auto u_exact = Px()+Py()+Pz();
        //auto u_exact = Px()*Px()+Py()*Py()+Pz()*Pz();
        auto u_exact = Px()*Px()*Pz()+Py()*Py()*Px()+Pz()*Pz()*Py();
        //auto u_exact = Px();
        u = vf::project( Xh, elements( mesh ), u_exact );

        double avgv = integrate( internalfaces( mesh ), averagev( idv( u )-( u_exact ) ) ).evaluate()( 0, 0 );
        BOOST_TEST_MESSAGE( "avg(v-uexact)=" << avgv << "\n" );
        double avgv1 = integrate( internalfaces( mesh ), .5*( leftfacev( idv( u )-( u_exact ) )+rightfacev( idv( u )-( u_exact ) ) ) ).evaluate()( 0, 0 );
        BOOST_TEST_MESSAGE( ".5*(left(v-uexact)+right(v-uexact))=" << avgv1 << "\n" );
        double avgv2 = integrate( internalfaces( mesh ), .5*( leftfacev( idv( u ) )-leftfacev( u_exact )+rightfacev( idv( u ) )-rightfacev( u_exact ) ) ).evaluate()( 0, 0 );
        BOOST_CHECK_SMALL( avgv, eps );
        BOOST_CHECK_SMALL( avgv1, eps );
        BOOST_CHECK_SMALL( avgv2, eps );


        double leftv = integrate( internalfaces( mesh ), leftfacev( idv( u )-( u_exact ) ) ).evaluate()( 0, 0 );
        double rightv = integrate( internalfaces( mesh ), rightfacev( idv( u )-( u_exact ) ) ).evaluate()( 0, 0 );
        BOOST_TEST_MESSAGE( "leftv=" << leftv << "\n" );
        BOOST_CHECK_SMALL( leftv, eps );
        BOOST_TEST_MESSAGE( "rightv=" << rightv << "\n" );
        BOOST_CHECK_SMALL( rightv, eps );


        double n_jumpun = integrate( internalfaces( mesh ), trans( leftfacev( N() ) )*jumpv( idv( u ) ) ).evaluate()( 0, 0 );
        BOOST_CHECK_SMALL( n_jumpun, eps );
        double n_leftun_rightun = integrate( internalfaces( mesh ), trans( leftfacev( N() ) )*( leftfacev( idv( u )*N() )+rightfacev( idv( u )*N() ) ) ).evaluate()( 0, 0 );
        BOOST_CHECK_SMALL( n_leftun_rightun, eps );

        BOOST_CHECK_CLOSE( n_leftun_rightun, n_jumpun, eps );

        double left_n = integrate( internalfaces( mesh ), leftfacev( N() ) ).evaluate()( 0, 0 );
        double right_n = integrate( internalfaces( mesh ), rightfacev( N() ) ).evaluate()( 0, 0 );
        BOOST_CHECK_CLOSE( left_n, -right_n, eps );

        u = vf::project( Xh, elements( mesh ), cst( 1. ) );
        double leftv_1 = integrate( internalfaces( mesh ), leftfacev( idv( u ) ) ).evaluate()( 0, 0 );
        double rightv_1 = integrate( internalfaces( mesh ), rightfacev( idv( u ) ) ).evaluate()( 0, 0 );
        double sumv_1 = integrate( internalfaces( mesh ),  rightfacev( idv( u ) )+leftfacev( idv( u ) ) ).evaluate()( 0, 0 );
        double avgv_1 = integrate( internalfaces( mesh ),  averagev( idv( u ) ) ).evaluate()( 0, 0 );
        BOOST_CHECK_CLOSE( leftv_1, rightv_1, eps );
        BOOST_CHECK_CLOSE( leftv_1+rightv_1, sumv_1, eps );
        BOOST_CHECK_CLOSE( 2*avgv_1, sumv_1, eps );

    }
    boost::shared_ptr<Feel::Backend<double> > backend;
    double meshSize;
    std::string shape;
    mesh_ptrtype mesh;

};

template<typename value_type = double, int Dim=2>
struct test_integration_internal_faces_lf : public Application
{
    typedef typename imesh<value_type,Dim>::convex_type convex_type;
    typedef typename imesh<value_type,Dim>::type mesh_type;
    typedef typename imesh<value_type,Dim>::ptrtype mesh_ptrtype;
    typedef FunctionSpace<mesh_type, bases<Lagrange<3, Scalar> >, double> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    test_integration_internal_faces_lf()
        :
        Application(),
        backend( Backend<double>::build( soption( _name="backend" ) ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( this->vm()["shape"].template as<std::string>() ),
        mesh()
    {
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % Dim ).str() ,
                                             _usenames=true,
                                             _convex=( convex_type::is_hypercube )?"Hypercube":"Simplex",
                                             _shape=shape,
                                             _dim=Dim,
                                             _h=meshSize ) );
    }

    void operator()()
    {
        using namespace Feel::vf;

        const value_type eps = 1e-10;


        space_ptrtype Xh = space_type::New( _mesh=mesh,_extended_doftable=std::vector<bool>(1,true) );

        auto u = Xh->element( "u" );
        //auto u_exact = Px()+Py()+Pz();
        //auto u_exact = Px()*Px()+Py()*Py()+Pz()*Pz();
        auto u_exact = Px()*Px()*Pz()+Py()*Py()*Px()+Pz()*Pz()*Py();
        auto v_exact = u_exact *unitX() + u_exact*unitY()+ u_exact*unitZ();
        u = vf::project( Xh, elements( mesh ), u_exact );


        auto F = backend->newVector( Xh );
        auto _F_ = integrate( internalfaces( mesh ),
                trans( leftface( id( u )*N() )+rightface( id( u )*N() ) )*leftfacev( N() ) );
        form1( _test=Xh, _vector=F, _init=true ) = _F_;
        //  integrate( internalfaces( mesh ),
        //        //print(trans(print(leftface(id(u)*print(N(),"leftN:")),"leftuN=")+print(rightface(id(u)*print(N(),"rightN:")),"rightuN=")),"leftuN+rightuN=" )*print(leftfacev(N()),"leftN=")
        //        trans( leftface( id( u )*N() )+rightface( id( u )*N() ) )*leftfacev( N() )
        //                                                    );

        F->close();
        F->printMatlab( "F.m" );
        u.printMatlab( "u.m" );
        double jumpu_F = inner_product( u, *F );
        BOOST_TEST_MESSAGE ( "jump(u) = " << jumpu_F << "\n" );
        BOOST_CHECK_SMALL( jumpu_F, eps );

#if 1
        form1( _test=Xh, _vector=F, _init=true ) = integrate( internalfaces( mesh ),
                ( jump( grad( u ) ) ) );
        double jump_gradu_F = inner_product( u, *F );
        BOOST_TEST_MESSAGE ( "jump(grad(u) u^T F = " << jump_gradu_F << "\n" );
        BOOST_CHECK_SMALL( jump_gradu_F, eps );

        form1(_test=Xh, _vector=F, _init=true ) = integrate( internalfaces( mesh ), leftface( grad( u )*N() ) );
        double left_gradu_n = inner_product( u, *F );
        form1( _test=Xh, _vector=F, _init=true ) = integrate( internalfaces( mesh ), rightface( grad( u )*N() ) );
        double right_gradu_n = inner_product( u, *F );
        BOOST_TEST_MESSAGE(  "jump(left(grad(u)*N)) u^T F = " << left_gradu_n << "\n" );
        BOOST_TEST_MESSAGE(  "jump(right(grad(u)*N)) u^T F = " << right_gradu_n << "\n" );
        BOOST_CHECK_CLOSE( left_gradu_n, -right_gradu_n, eps*100 );

        u = vf::project( Xh, elements( mesh ), cst( 1. ) );
        form1( _test=Xh, _vector=F, _init=true ) = integrate( internalfaces( mesh ), leftface( id( u ) ) );
        double left_1 = inner_product( u, *F );
        BOOST_TEST_MESSAGE(  "left(id(u)) u^T F = " << left_1 << "\n" );
        form1( _test=Xh, _vector=F, _init=true ) = integrate( internalfaces( mesh ), rightface( id( u ) ) );
        double right_1 = inner_product( u, *F );
        BOOST_TEST_MESSAGE(  "right(id(u)) u^T F = " << right_1 << "\n" );
        BOOST_CHECK_CLOSE( left_1, right_1, eps );

        form1( _test=Xh, _vector=F, _init=true ) = integrate( internalfaces( mesh ), trans( N() )*jump( id( u ) ) );
        BOOST_CHECK_SMALL( inner_product( u, *F ), eps );

        form1( _test=Xh, _vector=F, _init=true ) = integrate( internalfaces( mesh ), jump( grad( u ) ) );
        BOOST_CHECK_SMALL( inner_product( u, *F ), eps );
#endif

    }
    boost::shared_ptr<Feel::Backend<double> > backend;
    double meshSize;
    std::string shape;
    mesh_ptrtype mesh;
};
} // Feel
inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description integrationoptions( "Test Integration 2D/3D options" );
    integrationoptions.add_options()
    ( "hsize", Feel::po::value<double>()->default_value( 3 ), "h value" )
    ( "shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (hypercube, simplex, ellipsoid)" )
    ;
    return integrationoptions.add( Feel::feel_options() );
}

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "test_integration_ifaces" ,
                           "test_integration_ifaces" ,
                           "0.2",
                           "1D/2D/3D internal faces integration tests",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2006-2010 Université Joseph Fourier (Grenoble I)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( integration )

//typedef boost::mpl::list<boost::mpl::int_<1>,boost::mpl::int_<2>,boost::mpl::int_<3> > dim_types;
typedef boost::mpl::list<boost::mpl::int_<1>,boost::mpl::int_<2> > dim_types;
//typedef boost::mpl::list<boost::mpl::int_<3> > dim_types;
//typedef boost::mpl::list<boost::mpl::int_<2>,boost::mpl::int_<3>,boost::mpl::int_<1> > dim_types;
#if 1
BOOST_AUTO_TEST_CASE_TEMPLATE( test_integration_ifaces_v, T, dim_types )
{
    BOOST_TEST_MESSAGE( "Test integration on internal faces v (" << T::value << "D)" );
    Feel::test_integration_internal_faces_v<double,T::value> t;
    t();
    BOOST_TEST_MESSAGE( "Test integration on internal faces v (" << T::value << "D) done." );
}
#endif // 0
BOOST_AUTO_TEST_CASE_TEMPLATE( test_integration_ifaces_lf, T, dim_types )
{


    BOOST_TEST_MESSAGE( "Test integration on internal faces in linear forms (" << T::value << "D)" );
    Feel::test_integration_internal_faces_lf<double,T::value> t;
    t();
    BOOST_TEST_MESSAGE( "Test integration on internal faces in linear forms (" << T::value << "D) done" );
}
BOOST_AUTO_TEST_SUITE_END()

#if 0
int BOOST_TEST_CALL_DECL
main( int argc, char* argv[] )
{
    Feel::Environment env( argc, argv );
    Feel::Assert::setLog( "test_integration_ifaces.assert" );
    int ret = ::boost::unit_test::unit_test_main( &init_unit_test, argc, argv );

    return ret;
}
#endif

#if 0
#if defined(USE_BOOST_TEST)
boost::shared_ptr<Feel::Application> mpi;
test_suite*
init_unit_test_suite( int argc, char** argv )
{
    mpi = boost::shared_ptr<Feel::Application>( new Feel::Application( argc, argv, makeAbout(), makeOptions() ) );
    Feel::Assert::setLog( "test_integration.assert" );
    test_suite* test = BOOST_TEST_SUITE( "2D Generic finite element solver test suite" );

    test->add( BOOST_TEST_CASE( ( test_integration_internal_faces<double>( mpi->vm()["hsize"].as<double>() ) ) ) );

    return test;
}
#else
int
main( int argc, char** argv )
{
    Feel::Application mpi( argc, argv, makeAbout(), makeOptions() );
    Feel::Assert::setLog( "test_integration_ifaces.assert" );

    test_integration_internal_faces<double> c ( mpi.vm()["hsize"].as<double>() );
    c();
}
#endif /* USE_BOOST_TEST */
#endif
