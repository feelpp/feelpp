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
#define BOOST_TEST_MODULE interpolation operator testsuite
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
        _h=hsize ) );
  //typedef typename imesh<Dim, Order, RDim>::type meshType;
  //return loadMesh( _mesh=new meshType,
  //                       _desc=domain(_addmidpoint=false,
  //                                     _h=hsize ) );

}


template<int Dim, int Order, int RDim, int GeoOrder=1, typename value_type = double>
struct test_interpolation_op
{
  typedef typename imesh<Dim, GeoOrder>::type mesh_type;
  test_interpolation_op( double meshSize_=DEFAULT_MESH_SIZE )
    :
      meshSize( meshSize_ ),
      mesh( createMesh<Dim,GeoOrder,Dim>( meshSize ) )
  {}
  void operator()()
  {
    using namespace Feel;
    using namespace Feel::vf;

    BOOST_MESSAGE(  "[test_interpolation_op]   nDim : " << Dim << "\n" );
    BOOST_MESSAGE(  "[test_interpolation_op] nOrder : " << Order << "\n" );

    const value_type eps = 1000*Feel::type_traits<value_type>::epsilon();

    typedef fusion::vector<Lagrange<Order+2+GeoOrder-1, Scalar> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    boost::shared_ptr<space_type> Xh( new space_type( mesh ) );

    typename space_type::element_type u( Xh, "u" );

    auto u_exact = Px()*Px();
    u = vf::project( Xh, elements( *mesh ), u_exact );

    typename imesh<Dim,GeoOrder,RDim>::ptrtype mesh1( createMesh<Dim,GeoOrder,RDim>( meshSize/2 ) );
    typedef typename imesh<Dim,GeoOrder,RDim>::type mesh1_type;

    typedef fusion::vector<Lagrange<Order+GeoOrder-1, Scalar> > imagebasis_type;
    typedef FunctionSpace<mesh1_type, imagebasis_type> imagespace_type;
    boost::shared_ptr<imagespace_type> Yh( new imagespace_type( mesh1 ) );
    typename imagespace_type::element_type v( Yh, "v" );
    BOOST_MESSAGE(  "[test_interpolation_op] functionspace allocated\n" );

    typedef Backend<double> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
#if 0
#if defined( FEELPP_HAS_PETSC_H )
    backend_ptrtype backend( backend_type::build( BACKEND_PETSC ) );
#else
    backend_ptrtype backend( backend_type::build( BACKEND_GMM ) );
#endif
#else
    backend_ptrtype backend( backend_type::build( soption( _name="backend" ) ) );
#endif
    BOOST_MESSAGE(  "[test_interpolation_op] backend allocated\n" );

    auto I = opInterpolation( _domainSpace=Xh, _imageSpace=Yh, _backend=backend );


    BOOST_MESSAGE(  "[test_interpolation_op] OI allocated\n" );
    FsFunctionalLinear<imagespace_type> fsv( Yh );
    BOOST_MESSAGE(  "[test_interpolation_op] FSV allocated\n" );
    I->apply( u, fsv );
    BOOST_MESSAGE(  "[test_interpolation_op] applied OI \n" );

    auto y = Yh->element();
    auto Ih = opInterpolation( _domainSpace=Xh, _imageSpace=Yh, _range=elements( Yh->mesh() ) );
    Ih->apply( u, y );

    interpolate( Yh, u, v );
    BOOST_MESSAGE(  "[test_interpolation_op] interpolate\n" );


    typename imagespace_type::element_type w( Yh, "w" );
    w = fsv.container();

    double err = math::sqrt( integrate( elements( Yh->mesh() ), ( idv( w )-idv( y ) )*( idv( w )-idv( y ) ) ).evaluate()( 0, 0 ) );
    BOOST_MESSAGE(  "[err] ||w-y||_2 = " << err << "\n" );
    BOOST_CHECK_SMALL( err, 1e-12 );

    //std::cout << "w=" << w << "\n" );
    value_type xw = math::sqrt( integrate( elements( mesh1 ), ( u_exact-idv( w ) )*( u_exact-idv( w ) ) ).evaluate()( 0, 0 ) );
    value_type vw = math::sqrt( integrate( elements( mesh1 ), ( idv( v )-idv( w ) )*( idv( v )-idv( w ) ) ).evaluate()( 0, 0 ) );
    BOOST_MESSAGE(  "[test_interpolation_op] ||x-w||_2 = " << xw << "\n" );
    BOOST_MESSAGE(  "[test_interpolation_op] ||v-w||_2 = " << vw << "\n" );

    std::ostringstream ostr;
    ostr << "ointerpu-" << Dim << "." << Order;
    ExporterQuick<mesh_type> exp( ostr.str(), "ensight" );
    exp.save( 0, u );

    std::ostringstream ostr2;
    ostr2 << "ointerpv-" << Dim << "." << Order;
    ExporterQuick<mesh1_type> exp2( ostr2.str(), "ensight" );
    exp2.save( 0, v, w );
  }
  double meshSize;
  typename imesh<Dim>::ptrtype mesh;
};

template<int DimDomain, int OrderDomain, int RealDimDomain,
  int DimImage, int OrderImage, int RealDimImage>
  struct test_interpolation_op_2
{
  typedef typename imesh<DimDomain, 1, RealDimDomain>::type domain_mesh_type;
  typedef double value_type;
  test_interpolation_op_2( double meshSize_=DEFAULT_MESH_SIZE )
    :
      meshSize( meshSize_ ) ,
      mesh( createMesh<DimDomain,1,RealDimDomain>( meshSize ) )
  {}
  void operator()()
  {
    using namespace Feel;
    using namespace Feel::vf;

    BOOST_MESSAGE(  "[test_interpolation_op_2] domain   nDim : " << DimDomain << "\n" );
    BOOST_MESSAGE(  "[test_interpolation_op_2] domain nOrder : " << OrderDomain << "\n" );
    BOOST_MESSAGE(  "[test_interpolation_op_2]  image   nDim : " << DimImage << "\n" );
    BOOST_MESSAGE(  "[test_interpolation_op_2]  image nOrder : " << OrderImage << "\n" );




#if 0
    typedef fusion::vector<Lagrange<OrderDomain, Scalar> > domain_basis_type;
    typedef FunctionSpace<mesh_type, domain_basis_type> domain_space_type;
    boost::shared_ptr<domain_space_type> Xh( new domain_space_type( mesh ) );

    typename domain_space_type::element_type u( Xh, "u" );

    u = vf::project( Xh, elements( mesh ), Px() );

    typename imesh<DimImage,1,RealDimImage>::ptrtype image_mesh( createMesh<DimImage,1,RealDimImage>( meshSize/2 ) );

    typedef fusion::vector<Lagrange<Order, Scalar> > imagebasis_type;
    typedef FunctionSpace<mesh_type, imagebasis_type> imagespace_type;
    boost::shared_ptr<imagespace_type> Yh( new imagespace_type( mesh1 ) );
    typename imagespace_type::element_type v( Yh, "v" );

    typedef Backend<double> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    backend_ptrtype backend( Backend<double>::build( BACKEND_GMM ) );
    OperatorInterpolation<space_type, imagespace_type> I( Xh, Yh, backend );


    FsFunctionalLinear<imagespace_type> fsv( Yh );
    I.apply( u, fsv );

    interpolate( Yh, u, v );


    typename imagespace_type::element_type w( Yh, "w" );
    w = fsv.container();
    //std::cout << "w=" << w << "\n" );
    value_type xw = math::sqrt( integrate( elements( mesh1 ), IM<Dim,Order,value_type,Simplex>(), ( Px()-idv( w ) )*( Px()-idv( w ) ) ).evaluate()( 0, 0 ) );
    value_type vw = math::sqrt( integrate( elements( mesh1 ), IM<Dim,Order,value_type,Simplex>(), ( idv( v )-idv( w ) )*( idv( v )-idv( w ) ) ).evaluate()( 0, 0 ) );
    BOOST_MESSAGE(  "[test_interpolation_op] ||x-w||_2 = " << xw << "\n" );
    BOOST_MESSAGE(  "[test_interpolation_op] ||v-w||_2 = " << vw << "\n" );

    std::ostringstream ostr;
    ostr << "ointerpu-" << Dim << "." << Order;
    ExporterQuick<mesh_type> exp( ostr.str(), "ensight" );
    exp.save( 0, u );

    std::ostringstream ostr2;
    ostr2 << "ointerpv-" << Dim << "." << Order;
    ExporterQuick<mesh_type> exp2( ostr2.str(), "ensight" );
    exp2.save( 0, v, w );
#endif
  }
  double meshSize;
  typename imesh<DimDomain,1,RealDimDomain>::ptrtype mesh;
};

template<int Dim, int Order, int GeoOrder = 1, typename value_type = double>
struct test_lagrange_p1_op
{
  typedef typename imesh<Dim,GeoOrder,Dim>::type mesh_type;
  test_lagrange_p1_op( double meshSize_=DEFAULT_MESH_SIZE )
    :
      meshSize( meshSize_ ),
      mesh( createMesh<Dim, GeoOrder, Dim>( meshSize ) )

  {}
  void operator()()
  {
    using namespace Feel;
    using namespace Feel::vf;
    BOOST_MESSAGE(  "[test_lagrange_p1_op]   nDim : " << Dim << "\n" );
    BOOST_MESSAGE(  "[test_lagrange_p1_op] nOrder : " << Order << "\n" );




    //const value_type eps = 1000*Feel::type_traits<value_type>::epsilon();

    typedef fusion::vector<Lagrange<Order, Scalar> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    boost::shared_ptr<space_type> Xh( new space_type( mesh ) );

    typename space_type::element_type u( Xh, ( boost::format( "u_%1%.%2%.%3%" ) % Dim % Order % GeoOrder ).str() );

    u = vf::project( Xh, elements( *mesh ), Px() );

#if 0
    std::ostringstream ostr1;
    ostr1 << "olagp1-init-" << Dim << "." << Order << "." << GeoOrder;
    ExporterQuick<mesh_type> exp1( ostr1.str(), "ensight" );
    exp1.save( 0, u );
#endif

    typedef Backend<value_type> backend_type;
#if defined( FEELPP_HAS_PETSC_H )
    boost::shared_ptr<backend_type> backend( backend_type::build( BACKEND_PETSC ) );
#else
    boost::shared_ptr<backend_type> backend( backend_type::build( BACKEND_GMM ) );
#endif
    OperatorLagrangeP1<space_type> I( Xh, backend );
    typedef typename OperatorLagrangeP1<space_type>::dual_image_space_type::mesh_type image_mesh_type;
    typename OperatorLagrangeP1<space_type>::dual_image_space_ptrtype Yh( I.dualImageSpace() );
    typename OperatorLagrangeP1<space_type>::dual_image_space_type::element_type w( Yh, ( boost::format( "w_%1%.%2%.%3%" ) % Dim % Order % GeoOrder ).str() );
    typename OperatorLagrangeP1<space_type>::dual_image_space_type::element_type e( Yh, ( boost::format( "w_%1%.%2%.%3%" ) % Dim % Order % GeoOrder ).str() );
    typename OperatorLagrangeP1<space_type>::dual_image_space_type::element_type yy( Yh, ( boost::format( "yy_%1%.%2%.%3%" ) % Dim % Order % GeoOrder ).str() );

    FsFunctionalLinear<typename OperatorLagrangeP1<space_type>::dual_image_space_type> fsv( Yh );
    I.apply( u, fsv );

    w = fsv.container();
    //std::cout << "u=" << u << "\n" );
    //std::cout << "w=" << w << "\n" );

    value_type xw = math::sqrt( integrate( elements( Yh->mesh() ),
          IM<Dim,2*Order,value_type,Simplex>(),
          ( Px()-idv( w ) )*( Px()-idv( w ) ) ).evaluate()( 0, 0 ) );
    BOOST_MESSAGE(  "[test_lagrange_p1_op] ||x-w||_2 = " << xw << "\n" );
    e=w;
    e.setName( ( boost::format( "e_%1%.%2%.%3% " ) % Dim % Order % GeoOrder ).str() );
    e-=u;
    //std::cout << "e=" << e << "\n" );
    BOOST_MESSAGE(  "[test_lagrange_p1_op] ||x-w||_infty = " << e.linftyNorm() << "\n" );

    yy = vf::project( Yh, elements( Yh->mesh() ), Px() );
    std::ostringstream ostr;
    ostr << "olagp1-" << Dim << "." << Order << "." << GeoOrder;
    ExporterQuick<image_mesh_type> exp( ostr.str(), "ensight" );
    exp.save( 0, w, e, yy );
  }
  double meshSize;
  typename imesh<Dim,GeoOrder,Dim>::ptrtype mesh;
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
  Feel::AboutData about( "test_interpolation_op" ,
                         "test_interpolation_op" ,
                         "0.1",
                         "interpolation operator tests",
                         Feel::AboutData::License_GPL,
                         "Copyright (C) 2007-2013 Université Joseph Fourier (Grenoble I)\n"
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

BOOST_AUTO_TEST_SUITE( interpolation_op_suite )

  typedef Feel::Application Application_type;
  typedef boost::shared_ptr<Application_type> Application_ptrtype;

BOOST_AUTO_TEST_CASE( test_interpolation_op_121 )
{
  auto myApp = Application_ptrtype(new Feel::Application);
  BOOST_MESSAGE( "test_interpolation_op<1,2,1>" );
  test_interpolation_op<1,2,1> t( myApp->vm()["hsize"].as<double>() );
  t();
}
BOOST_AUTO_TEST_CASE( test_interpolation_op_212 )
{
  auto myApp = Application_ptrtype(new Feel::Application);
  BOOST_MESSAGE( "test_interpolation_op<2,1,2>" );
  test_interpolation_op<2,1,2> t( myApp->vm()["hsize"].as<double>() );
  t();
}
BOOST_AUTO_TEST_CASE( test_interpolation_op_313 )
{
  auto myApp = Application_ptrtype(new Feel::Application);
  BOOST_MESSAGE( "test_interpolation_op<3,1,3>" );
  test_interpolation_op<3,1,3> t( myApp->vm()["hsize"].as<double>() );
  t();
}
//BOOST_AUTO_TEST_CASE( test_interpolation_op112 ) { BOOST_MESSAGE( "test_interpolation_op<1,1,2>"); test_interpolation_op<1,1,2> t( 0.1 ); t(); }
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

#if 1
  test_interpolation_op<1,1,1,1> t( mpiapp.vm()["hsize"].as<double>() );
  t();
  test_interpolation_op<2,1,2,1> t21( mpiapp.vm()["hsize"].as<double>() );
  t21();


  test_interpolation_op<3,1,3,1> t31( mpiapp.vm()["hsize"].as<double>() );
  t31();
#else
  //test_interpolation<2,1,2> t212( mpiapp.vm()["hsize"].as<double>() );
  //t212();
#endif
#if 0
  test_lagrange_p1_op<1,1> tlp1_11( mpiapp->vm()["hsize"].as<double>() );
  tlp1_11();
  test_lagrange_p1_op<1,2> tlp1_12( mpiapp->vm()["hsize"].as<double>() );
  tlp1_12();
  test_lagrange_p1_op<1,3> tlp1_13( mpiapp->vm()["hsize"].as<double>() );
  tlp1_13();
  test_lagrange_p1_op<1,4> tlp1_14( mpiapp->vm()["hsize"].as<double>() );
  tlp1_14();
#endif
  //test_lagrange_p1_op<2,6> t23( mpiapp->vm()["hsize"].as<double>() );
  //t23();
  //test_lagrange_p1_op<3,1> t31( mpiapp->vm()["hsize"].as<double>() );
  //t31();
  //test_lagrange_p1_op<2,2,2> t21( mpiapp.vm()["hsize"].as<double>() );t21();
  //test_interpolation_op<1,1,2> t21( mpiapp.vm()["hsize"].as<double>() );t21();
  //test_lagrange_p1_op<2,2> t22( mpiapp.vm()["hsize"].as<double>() );t22();

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
