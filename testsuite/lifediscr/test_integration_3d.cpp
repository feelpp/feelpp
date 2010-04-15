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
    ImporterGmsh<typename imesh<T,Dim>::type> import( td.generate( "Simplex" ) );
    mesh->accept( import );


    return mesh;
}
}
template<typename value_type = double>
struct test_integration_internal_faces
{
    test_integration_internal_faces( double meshSize_=DEFAULT_MESH_SIZE ): meshSize(meshSize_), mesh( Life::createMesh<value_type,3>( meshSize ) )
    {}

    void operator()()
    {
        using namespace Life;
        using namespace Life::vf;


        typedef typename imesh<value_type,3>::type mesh_type;
        typedef typename imesh<value_type,3>::ptrtype mesh_ptrtype;
        typename imesh<value_type,3>::ptrtype mesh( createMesh<value_type,3>( meshSize ) );

        const value_type eps = 1000*Life::type_traits<value_type>::epsilon();
#if 0
        // int ([-1,1],[-1,x]) 1 dx
        value_type v0 = 1;//integrate( elements(*mesh), IM<3,1,value_type,Simplex>(),
        //vf::min(constant(1.0),constant(2.0)) ).evaluate()( 0, 0 );
        std::cout << "v0 = " << v0 << "\n";
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v0-8.0, eps );
#else
        LIFE_ASSERT( math::abs( v0-1.0) < eps )( v0 )( math::abs( v0-1.0) )( eps ).warn ( "v0 != 1" );
#endif /* USE_BOOST_TEST */

        value_type v1 = integrate( internalfaces(*mesh), IM<3,1,value_type,Simplex>(),
                                   jumpv(trans(2*P())) ).evaluate()( 0, 0 );
        std::cout << "v1 = " << v1 << "\n";
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v1-0.0, eps );
#else
        LIFE_ASSERT( math::abs( v1-0.0) < eps )( v1 )( math::abs( v1-0.0) )( eps ).warn ( "v1 != 0" );
#endif /* USE_BOOST_TEST */

        value_type v2 = integrate( internalfaces(*mesh), IM<3,3,value_type,Simplex>(),
                                   leftfacev(vf::sqrt(trans(P())*P()))-rightfacev(vf::sqrt(trans(P())*P())) ).evaluate()( 0, 0 );
        std::cout << "v2 = " << v2 << "\n";
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v2-0.0, eps );
#else
        LIFE_ASSERT( math::abs( v2-0.0) < eps )( v2 )( math::abs( v2-0.0) )( eps ).warn ( "v2 != 0" );
#endif /* USE_BOOST_TEST */

        value_type v3 = integrate( internalfaces(*mesh), IM<3,3,value_type,Simplex>(),
                                   leftfacev(vf::sqrt(trans(P()*P())*(P()*P())))-rightfacev(vf::sqrt(trans(P()*P())*(P()*P()))) ).evaluate()( 0, 0 );
        std::cout << "v3 = " << v3 << "\n";
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v3-0.0, eps );
#else
        LIFE_ASSERT( math::abs( v3-0.0) < eps )( v3 )( math::abs( v3-0.0) )( eps ).warn ( "v3 != 0" );
#endif /* USE_BOOST_TEST */
#endif

        typedef FunctionSpace<mesh_type, fusion::vector<Lagrange<3, Scalar> >, double> space_type;
        typedef boost::shared_ptr<space_type> space_ptrtype;
        space_ptrtype Xh = space_type::New( mesh );
        typedef typename space_type::element_type element_type;
        element_type u( Xh, "u" );
        u = project( Xh, elements( *mesh ), Px()+Py()+Pz());

        boost::shared_ptr<VectorUblas<double> > F( new VectorUblas<double>( u.size() ) );
        std::fill( F->begin(), F->end(), (double)0 );
        form1( Xh, F ) = integrate( internalfaces(*mesh), IM<3,6,value_type,Simplex>(),
                                    //print( trans( jump( id( u ) ) )*jump( id( u ) ), "[u]^T [u]" ) );
                                    //(print( leftfacev( trans(P())*P() ), "left P(): " ) - print( rightfacev( trans(P())*P() ), "right P(): " ))*id(u) );
                                    //trans(leftface( id( u )*N() )+ rightface( id( u )*N() ))*
                                    //(leftface( id( u )*N() )+ rightface( id( u )*N() ))
                                    (2*leftface( grad( u )*N() )+ 2*rightface( grad( u )*N() ))
                                    //trans(one())*(jump(id(u)))

                                    );
        F->printMatlab( "F.m" );
        std::cout << "u^T F = " << inner_product( u, *F ) << "\n";

#if 0
        MatrixGmm<double,gmm::row_major> M;
        form( Xh, Xh, M ) = integrate( internalfaces(*mesh), IM<3,2,value_type,Simplex>(),
                                       trans(leftfacet( idt( u )*N() )+ rightfacet( idt( u )*N() ))*
                                       (jump(id(u)))

                                       );
        M.close();
        M.printMatlab( "M.m" );
        //std::cout << "u^T F = " << ublas::inner_prod( u, F ) << "\n";

#endif
    }
    double meshSize;
    typename Life::imesh<value_type,3>::ptrtype mesh;
};

inline
Life::po::options_description
makeOptions()
{
    Life::po::options_description integrationoptions("Test Integration 3D options");
    integrationoptions.add_options()
        ("hsize", Life::po::value<double>()->default_value( 0.3 ), "h value")
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
                           "3D integration tests",
                           Life::AboutData::License_GPL,
                           "Copyright (C) 2006,2007,2008 Université Joseph Fourier (Grenoble I)");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}

BOOST_AUTO_TEST_CASE( test_integration_internal_faces_double )
{
    test_integration_internal_faces<double> t( 0.3 );
    t();

}
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
