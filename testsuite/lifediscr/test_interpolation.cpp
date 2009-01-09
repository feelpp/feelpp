/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-12-19

  Copyright (C) 2007-2008 Université Joseph Fourier (Grenoble I)

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
   \file test_interpolation.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-12-19
 */
//#define USE_BOOST_TEST 1

// Boost.Test
#define BOOST_TEST_MAIN
#if defined(USE_BOOST_TEST)
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;
#include <boost/test/floating_point_comparison.hpp>
#endif

#include <life/options.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifecore/application.hpp>
#include <life/lifealg/backend.hpp>
#include <life/lifemesh/geoentity.hpp>
#include <life/lifemesh/refentity.hpp>
#include <life/lifediscr/functionspace.hpp>

#include <life/lifediscr/mesh.hpp>
#include <life/lifediscr/operatorlagrangep1.hpp>
#include <life/lifediscr/interpolate.hpp>
#include <life/lifemesh/filters.hpp>
#include <life/lifepoly/im.hpp>
#include <life/lifefilters/importergmsh.hpp>
#include <life/lifefilters/gmsh.hpp>
#include <life/lifefilters/gmshtensorizeddomain.hpp>
#include <life/lifefilters/exporterquick.hpp>
#include <life/lifevf/vf.hpp>

const double DEFAULT_MESH_SIZE=0.5;

using namespace Life;

template<int Dim, int Order=1, int RDim=Dim>
struct imesh
{
    typedef Mesh<GeoEntity<Simplex<Dim,Order,RDim> >, double > type;
    typedef boost::shared_ptr<type> ptrtype;
};

template<int Dim, int Order, int RDim>
typename imesh<Dim, Order, RDim>::ptrtype
createMesh( double hsize )
{
    Gmsh __gmsh;
    std::string fname;
    std::ostringstream ostr;
    std::ostringstream nameStr;


    typename imesh<Dim, Order, RDim>::ptrtype mesh( new typename imesh<Dim, Order, RDim>::type );

    GmshTensorizedDomain<Dim,Order,RDim,Simplex> td;
    td.setCharacteristicLength( hsize );
    td.setX( std::make_pair( -1, 1 ) );
    //td.setX( std::make_pair( -1, 1 ) );
    td.setOrder( Order );
    fname = td.generate( Simplex<Dim,Order,RDim>::name() );

    ImporterGmsh<typename imesh<Dim,Order,RDim>::type> import( fname );
    mesh->accept( import );

    mesh->components().set( MESH_RENUMBER | MESH_UPDATE_FACES | MESH_UPDATE_EDGES );
    mesh->updateForUse();
    return mesh;
}

template<int Dim, int Order, int GeoOrder=1>
struct test_interpolation
{
    typedef typename imesh<Dim,GeoOrder>::type mesh_type;
    typedef double value_type;

    test_interpolation( double meshSize_=DEFAULT_MESH_SIZE )
        :
        meshSize(meshSize_),
        mesh( createMesh<Dim,GeoOrder,Dim>( meshSize ) )
    {}
    void operator()()
    {
        using namespace Life;
        using namespace Life::vf;

        const value_type eps = 1000*Life::type_traits<value_type>::epsilon();

        typedef fusion::vector<fem::Lagrange<Dim, Order, Scalar, Continuous, double> > basis_type;
        typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
        boost::shared_ptr<space_type> Xh( new space_type(mesh) );
        typename space_type::element_type u( Xh );

        u = project( Xh, elements(*mesh), constant(1.0) );

        node_type pt(Dim);
        pt[0] = 0.11;
        if ( Dim >= 2 )
            pt[1] = 0.11;
        if ( Dim >= 3 )
            pt[2] = 0.11;
        double v0 = u( pt )(0,0,0);
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v0-1.0, eps );
#else
        LIFE_ASSERT( math::abs( v0-1.0) < eps )( v0 )( math::abs( v0-1.0) )( eps ).warn ( "v0 != 1" );
#endif /* USE_BOOST_TEST */

        u = project( Xh, elements(*mesh), constant(2.0) - Px()*Px()-Py()*Py()-Pz()*Pz() );
        pt[0] = 0.5;
        if ( Dim >= 2 )
            pt[1] = 0.5;
        if ( Dim >= 3 )
            pt[2] = 0.5;

        double v1 = u( pt )(0,0,0);
        double v1_ex = 2-Dim*0.5*0.5;
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v1-v1_ex, eps );
#else
        LIFE_ASSERT( math::abs( v1-v1_ex) < eps )( v1 )( math::abs( v1-v1_ex) )( eps ).warn ( "v1 != v0_ex" );
#endif /* USE_BOOST_TEST */
        typename space_type::element_type::grad_type gradient( u.grad( pt ) );
        double g_v1_x = gradient(0,0,0);
        double g_v1_ex_x = -2*0.5;
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( g_v1_x-g_v1_ex_x, eps );
#else
        LIFE_ASSERT( math::abs( g_v1_x-g_v1_ex_x) < eps )( g_v1_x )( math::abs( g_v1_x-g_v1_ex_x) )( eps ).warn ( "g_v1 != g_v1_ex" );
#endif /* USE_BOOST_TEST */

        if ( Dim >= 2 )
            {
                double g_v1_y = gradient(0,1,0);
#if defined(USE_BOOST_TEST)
                BOOST_CHECK_SMALL( g_v1_y-g_v1_ex_x, eps );
#endif /* USE_BOOST_TEST */
            }
        if ( Dim >= 3 )
            {
                double g_v1_z = gradient(0,2,0);
#if defined(USE_BOOST_TEST)
                BOOST_CHECK_SMALL( g_v1_z-g_v1_ex_x, eps );
#endif /* USE_BOOST_TEST */
            }


    }
    double meshSize;



    typename imesh<Dim,GeoOrder>::ptrtype mesh;

};

template<int Dim, int Order, int RDim, typename value_type = double>
struct test_interpolation_op
{
    typedef typename imesh<Dim>::type mesh_type;
    test_interpolation_op( double meshSize_=DEFAULT_MESH_SIZE )
        :
        meshSize(meshSize_),
        mesh( createMesh<Dim,1,Dim>( meshSize ) )
    {}
    void operator()()
    {
        using namespace Life;
        using namespace Life::vf;

        Debug() << "[test_interpolation_op]   nDim : " << Dim << "\n";
        Debug() << "[test_interpolation_op] nOrder : " << Order << "\n";

        const value_type eps = 1000*Life::type_traits<value_type>::epsilon();

        typedef fusion::vector<fem::Lagrange<Dim, Order+2, Scalar, Continuous, double> > basis_type;
        typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
        boost::shared_ptr<space_type> Xh( new space_type(mesh) );

        typename space_type::element_type u( Xh, "u" );

        u = project( Xh, elements(*mesh), Px() );

        typename imesh<Dim,1,RDim>::ptrtype mesh1( createMesh<Dim,1,RDim>( meshSize/2 ) );

        typedef fusion::vector<fem::Lagrange<Dim, Order, Scalar, Continuous, double> > imagebasis_type;
        typedef FunctionSpace<mesh_type, imagebasis_type, value_type> imagespace_type;
        boost::shared_ptr<imagespace_type> Yh( new imagespace_type(mesh1) );
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
        //std::cout << "w=" << w << "\n";
        value_type xw = math::sqrt( integrate( elements( mesh1 ), IM<Dim,Order,value_type,Simplex>(), (Px()-idv(w))*(Px()-idv(w)) ).evaluate()( 0, 0 ) );
        value_type vw = math::sqrt( integrate( elements( mesh1 ), IM<Dim,Order,value_type,Simplex>(), (idv(v)-idv(w))*(idv(v)-idv(w)) ).evaluate()( 0, 0) );
        Debug() << "[test_interpolation_op] ||x-w||_2 = " << xw << "\n";
        Debug() << "[test_interpolation_op] ||v-w||_2 = " << vw << "\n";

        std::ostringstream ostr;
        ostr << "ointerpu-" << Dim << "." << Order;
        ExporterQuick<mesh_type> exp( ostr.str(), "ensight" );
        exp.save( 0, u );

        std::ostringstream ostr2;
        ostr2 << "ointerpv-" << Dim << "." << Order;
        ExporterQuick<mesh_type> exp2( ostr2.str(), "ensight" );
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
        meshSize(meshSize_) ,
        mesh( createMesh<DimDomain,1,RealDimDomain>( meshSize ) )
    {}
    void operator()()
    {
        using namespace Life;
        using namespace Life::vf;

        Debug() << "[test_interpolation_op_2] domain   nDim : " << DimDomain << "\n";
        Debug() << "[test_interpolation_op_2] domain nOrder : " << OrderDomain << "\n";
        Debug() << "[test_interpolation_op_2]  image   nDim : " << DimImage << "\n";
        Debug() << "[test_interpolation_op_2]  image nOrder : " << OrderImage << "\n";




#if 0
        typedef fusion::vector<fem::Lagrange<DimDomain, OrderDomain, Scalar, Continuous, double> > domain_basis_type;
        typedef FunctionSpace<mesh_type, domain_basis_type, value_type> domain_space_type;
        boost::shared_ptr<domain_space_type> Xh( new domain_space_type(mesh) );

        typename domain_space_type::element_type u( Xh, "u" );

        u = project( Xh, elements(mesh), Px() );

        typename imesh<DimImage,1,RealDimImage>::ptrtype image_mesh( createMesh<DimImage,1,RealDimImage>( meshSize/2 ) );

        typedef fusion::vector<fem::Lagrange<DimImage, Order, Scalar, Continuous, double> > imagebasis_type;
        typedef FunctionSpace<mesh_type, imagebasis_type, value_type> imagespace_type;
        boost::shared_ptr<imagespace_type> Yh( new imagespace_type(mesh1) );
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
        //std::cout << "w=" << w << "\n";
        value_type xw = math::sqrt( integrate( elements( mesh1 ), IM<Dim,Order,value_type,Simplex>(), (Px()-idv(w))*(Px()-idv(w)) ).evaluate()( 0, 0 ) );
        value_type vw = math::sqrt( integrate( elements( mesh1 ), IM<Dim,Order,value_type,Simplex>(), (idv(v)-idv(w))*(idv(v)-idv(w)) ).evaluate()( 0, 0) );
        Debug() << "[test_interpolation_op] ||x-w||_2 = " << xw << "\n";
        Debug() << "[test_interpolation_op] ||v-w||_2 = " << vw << "\n";

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
        meshSize(meshSize_),
        mesh( createMesh<Dim, GeoOrder, Dim>( meshSize ) )

    {}
    void operator()()
    {
        using namespace Life;
        using namespace Life::vf;
        Debug() << "[test_lagrange_p1_op]   nDim : " << Dim << "\n";
        Debug() << "[test_lagrange_p1_op] nOrder : " << Order << "\n";




        //const value_type eps = 1000*Life::type_traits<value_type>::epsilon();

        typedef fusion::vector<fem::Lagrange<Dim, Order, Scalar, Continuous, double> > basis_type;
        typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
        boost::shared_ptr<space_type> Xh( new space_type(mesh) );

        typename space_type::element_type u( Xh, (boost::format( "u_%1%.%2%.%3%" ) % Dim % Order % GeoOrder ).str() );

        u = project( Xh, elements(*mesh), Px() );

#if 0
        std::ostringstream ostr1;
        ostr1 << "olagp1-init-" << Dim << "." << Order << "." << GeoOrder;
        ExporterQuick<mesh_type> exp1( ostr1.str(), "ensight" );
        exp1.save( 0, u );
#endif

        typedef Backend<value_type> backend_type;
        boost::shared_ptr<backend_type> backend( backend_type::build( BACKEND_GMM ) );
        OperatorLagrangeP1<space_type> I( Xh, backend );
        typedef typename OperatorLagrangeP1<space_type>::dual_image_space_type::mesh_type image_mesh_type;
        typename OperatorLagrangeP1<space_type>::dual_image_space_ptrtype Yh( I.dualImageSpace() );
        typename OperatorLagrangeP1<space_type>::dual_image_space_type::element_type w( Yh, (boost::format( "w_%1%.%2%.%3%" ) % Dim % Order % GeoOrder ).str() );
        typename OperatorLagrangeP1<space_type>::dual_image_space_type::element_type e( Yh, (boost::format( "w_%1%.%2%.%3%" ) % Dim % Order % GeoOrder ).str() );
        typename OperatorLagrangeP1<space_type>::dual_image_space_type::element_type yy( Yh, (boost::format( "yy_%1%.%2%.%3%" ) % Dim % Order % GeoOrder ).str() );

        FsFunctionalLinear<typename OperatorLagrangeP1<space_type>::dual_image_space_type> fsv( Yh );
        I.apply( u, fsv );

        w = fsv.container();
        //std::cout << "u=" << u << "\n";
        //std::cout << "w=" << w << "\n";

        value_type xw = math::sqrt( integrate( elements( Yh->mesh() ),
                                               IM<Dim,2*Order,value_type,Simplex>(),
                                               (Px()-idv(w))*(Px()-idv(w)) ).evaluate()( 0, 0 ) );
        Debug() << "[test_lagrange_p1_op] ||x-w||_2 = " << xw << "\n";
        e=w;
        e.setName( (boost::format( "e_%1%.%2%.%3% " ) % Dim % Order % GeoOrder ).str() );
        e-=u;
        //std::cout << "e=" << e << "\n";
        Debug() << "[test_lagrange_p1_op] ||x-w||_infty = " << e.linftyNorm() << "\n";

        yy = project( Yh, elements(Yh->mesh()), Px() );
        std::ostringstream ostr;
        ostr << "olagp1-" << Dim << "." << Order << "." << GeoOrder;
        ExporterQuick<image_mesh_type> exp( ostr.str(), "ensight" );
        exp.save( 0, w, e, yy );
    }
    double meshSize;
    typename imesh<Dim,GeoOrder,Dim>::ptrtype mesh;
};

inline
Life::po::options_description
makeOptions()
{
    Life::po::options_description integrationoptions("Test Integration options");
    integrationoptions.add_options()
        ("hsize", Life::po::value<double>()->default_value( 0.3 ), "h value")
        ;
    return integrationoptions.add( Life::life_options() );
}

inline
Life::AboutData
makeAbout()
{
    Life::AboutData about( "test_interpolation" ,
                           "test_interpolation" ,
                            "0.1",
                           "interpolation tests",
                           Life::AboutData::License_GPL,
                           "Copyright (C) 2007, 2008 Université Joseph Fourier (Grenoble I)");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}

#if defined(USE_BOOST_TEST)
boost::shared_ptr<Life::Application> mpiapp;
test_suite*
init_unit_test_suite( int argc, char** argv )
{
    //boost::mpi::environment( argc, argv );
    mpiapp = boost::shared_ptr<Life::Application>( new Life::Application( argc, argv, makeAbout(), makeOptions() ) );
    Life::Assert::setLog( "test_integration.assert");
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
int
main( int argc, char** argv )
{
    Life::Application mpiapp( argc, argv, makeAbout(), makeOptions() );
    Life::Assert::setLog( "test_interpolation.assert");

#if 0
    test_interpolation_op<1,1> t( mpiapp->vm()["hsize"].as<double>() );
    t();
    test_interpolation_op<2,1> t21( mpiapp->vm()["hsize"].as<double>() );
    t21();


    test_interpolation_op<3,1> t31( mpiapp->vm()["hsize"].as<double>() );
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

