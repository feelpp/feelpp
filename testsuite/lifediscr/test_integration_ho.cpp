/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2006-08-25

  Copyright (C) 2006 EPFL
  Copyright (C) 2006,2007,2008 Université Joseph Fourier (Grenoble I)

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
   \file test_integration.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2006-08-25
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
#include <life/lifealg/backendgmm.hpp>
#include <life/lifemesh/geoentity.hpp>
#include <life/lifemesh/refentity.hpp>
#include <life/lifediscr/functionspace.hpp>
#include <life/lifediscr/mesh.hpp>
#include <life/lifemesh/filters.hpp>
#include <life/lifepoly/im.hpp>
#include <life/lifefilters/importergmsh.hpp>
#include <life/lifefilters/gmsh.hpp>
#include <life/lifefilters/gmshsimplexdomain.hpp>
#include <life/lifevf/vf.hpp>

const double DEFAULT_MESH_SIZE=0.1;

namespace Life
{
struct f_Px
{
    static const size_type context = vm::JACOBIAN|vm::POINT;
    typedef double value_type;
    static const uint16_type rank = 0;
    double operator()( uint16_type, uint16_type, ublas::vector<double> const& x, ublas::vector<double> const& /*n*/ ) const
    {
        return x[0];
    }
};
struct f_Nx
{
    static const size_type context = vm::JACOBIAN|vm::POINT|vm::NORMAL;
    typedef double value_type;
    static const uint16_type rank = 0;
    double operator()( uint16_type, uint16_type, ublas::vector<double> const& /*x*/, ublas::vector<double> const& n ) const
    {
        return n[0];
    }
};
struct f_Ny
{
    static const size_type context = vm::JACOBIAN|vm::POINT|vm::NORMAL;
    typedef double value_type;
    static const uint16_type rank = 0;
    double operator()( uint16_type, uint16_type, ublas::vector<double> const& /*x*/, ublas::vector<double> const& n ) const
    {
        return n[1];
    }
};
struct f_sinPx
{
    static const size_type context = vm::JACOBIAN|vm::POINT;
    typedef double value_type;
    static const uint16_type rank = 0;
    double operator()( uint16_type, uint16_type, ublas::vector<double> const& x, ublas::vector<double> const& /*n*/ ) const
    {
        return math::sin( x[0] );
    }
};

template<typename T, int Order = 1>
struct imesh
{
    typedef Mesh<GeoEntity<Simplex<2, Order> >, T > type;
    typedef boost::shared_ptr<type> ptrtype;
};

template<typename T>
typename imesh<T>::ptrtype
createMesh( double hsize )
{
    double meshSize = hsize;
    //std::cout << "hsize = " << meshSize << std::endl;

    Gmsh __gmsh;
    std::string fname;
    std::ostringstream ostr;
    std::ostringstream nameStr;

    //std::cout <<"Mesh generation ... ";
#if 0
    ostr << "h=" << meshSize << ";\n"
         << "Point(1) = {-1, -1,0.0,h};\n"
         << "Point(2) = { 1, -1,0.0,h};\n"
         << "Point(3) = {-1,  1,0.0,h};\n"
         << "Line(1) = {2,3};\n"
         << "Line(2) = {3,1};\n"
         << "Line(3) = {1,2};\n"
         << "Line Loop(4) = {1,2,3};\n"
         << "Plane Surface(5) = {4};\n"
         << "Physical Surface(30) = {5};\n"
         << "Physical Line(31) = {1};\n"
         << "Physical Line(32) = {2};\n"
         << "Physical Line(33) = {3};\n";

    nameStr << "triangle." << meshSize;

    fname = __gmsh.generate( nameStr.str(), ostr.str() );
#else
    ostr << "Mesh.MshFileVersion = 1;\n"
         << "h=" << meshSize << ";\n"
         << "Point(1) = {-1,-1,0,h};\n"
         << "Point(2) = {1,-1,0,h};\n"
         << "Point(3) = {1,1,0,h};\n"
         << "Point(4) = {-1,1,0,h};\n"
         << "Line(1) = {1,2};\n"
         << "Line(2) = {2,3};\n"
         << "Line(3) = {3,4};\n"
         << "Line(4) = {4,1};\n"
         << "Line Loop(4) = {1,2,3,4};\n"
         << "Plane Surface(5) = {4};\n"
         << "Physical Line(10) = {1,3};\n"
         << "Physical Line(20) = {2,4};\n"
         << "Physical Surface(6) = {5};\n";
    nameStr << "square." << meshSize;
    fname = __gmsh.generate( nameStr.str(), ostr.str() );
#endif


    /* Mesh */


    typename imesh<T>::ptrtype mesh( new typename imesh<T>::type );

    ImporterGmsh<typename imesh<T>::type> import( fname );
    mesh->accept( import );

    mesh->components().set( MESH_RENUMBER | MESH_UPDATE_FACES | MESH_UPDATE_EDGES );
    mesh->updateForUse();
    return mesh;
}
template<typename T, int Order>
typename imesh<T,Order>::ptrtype
createCircle( double hsize )
{
    double meshSize = hsize;
    typename imesh<T,Order>::ptrtype mesh( new typename imesh<T,Order>::type );
    typedef typename imesh<T,Order>::type mesh_type;
#if 1

    //std::cout << "hsize = " << meshSize << std::endl;

    Gmsh __gmsh;
    std::string fname;
    std::ostringstream ostr;
    std::ostringstream nameStr;

    //std::cout <<"Mesh generation ... ";
    ostr << "h=" << meshSize << ";\n"
         << "Point(1) = {0,0,0,h};\n"
         << "Point(2) = {1,0,0,h};\n"
         << "Point(3) = {0,1,0,h};\n"
         << "Point(4) = {-1,0,0,h};\n"
         << "Point(5) = {0,-1,0,h};\n"
         << "Circle(10) = {2,1,3};\n"
#if 1
         << "Circle(11) = {3,1,4};\n"
         << "Circle(12) = {4,1,5};\n"
         << "Circle(13) = {5,1,2};\n"
         << "Line Loop(14) = {10,11,12,13};\n"
         << "Plane Surface(15) = {14};\n"
         << "Physical Line(20) = {10,11};\n"
         << "Physical Line(21) = {12,13};\n"
         << "Physical Surface(30) = {15};\n";
#else
         << "Line(1) = {1,2};\n"
         << "Line(2) = {3,1};\n"
         << "Line Loop(14) = {1,10,2};\n"
         << "Plane Surface(15) = {14};\n"
         << "Physical Line(20) = {1,10,2};\n"
         << "Physical Surface(30) = {15};\n";
#endif



    nameStr << "circle." << meshSize;
    __gmsh.setOrder( Order );
    fname = __gmsh.generate( nameStr.str(), ostr.str() );



    ImporterGmsh<typename imesh<T, Order>::type> import( fname );
    mesh->accept( import );
#else

    typedef typename mesh_type::point_type point_type;
    typedef typename point_type::node_type node_type;
    typedef typename mesh_type::edge_type edge_type;
    typedef typename mesh_type::face_type face_type;
    typedef typename mesh_type::element_type element_type;

    node_type __nd( 2 );
    __nd[0] = 0;
    __nd[1] = -1;
    point_type __pt0( 0,__nd, true );
    __pt0.marker() = 0;
    mesh->addPoint( __pt0 );


    __nd[0] = 1;
    __nd[1] = 0;
    point_type __pt1( 1,__nd, true );
    __pt1.marker() = 0;
    mesh->addPoint( __pt1 );

    __nd[0] = 0;
    __nd[1] = 1;
    point_type __pt2( 2,__nd, true );
    __pt2.marker() = 0;
    mesh->addPoint( __pt2 );

    __nd[0] = math::sqrt(2.)/2;
    __nd[1] = math::sqrt(2.)/2;
    point_type __pt3( 3,__nd, true );
    __pt3.marker() = 0;
    mesh->addPoint( __pt3 );

    __nd[0] = -1;
    __nd[1] = 0;
    point_type __pt4( 4,__nd, true );
    __pt4.marker() = 0;
    mesh->addPoint( __pt4 );

    __nd[0] = math::sqrt(2.)/2;;
    __nd[1] = -math::sqrt(2.)/2;;
    point_type __pt5( 5,__nd, true );
    __pt5.marker() = 0;
    mesh->addPoint( __pt5 );

    element_type * pf = new element_type;

    pf->setId( 0 );
    pf->setMarker( 0 );

    // Warning : Vtk orientation is not the same as Life orientation !

    pf->setPoint( 0, mesh->point( 0 ) );
    pf->setPoint( 1, mesh->point( 1 ) );
    pf->setPoint( 2, mesh->point( 2 ) );
    pf->setPoint( 3, mesh->point( 3 ) );
    pf->setPoint( 4, mesh->point( 4 ) );
    pf->setPoint( 5, mesh->point( 5 ) );

    mesh->addElement( *pf );
    delete pf;


#endif
    mesh->components().set( MESH_RENUMBER | MESH_UPDATE_FACES | MESH_UPDATE_EDGES );
    mesh->updateForUse();
    return mesh;
}
template<typename T, int Order>
typename imesh<T,Order>::ptrtype
createSimplex( double hsize )
{
    double meshSize = hsize;
    //std::cout << "hsize = " << meshSize << std::endl;

    GmshSimplexDomain<2,Order> ts;
    ts.setCharacteristicLength( meshSize );
    ts.setOrder( Order );
    std::string fname = ts.generate( "simplex" );
    typename imesh<T,Order>::ptrtype mesh( new typename imesh<T,Order>::type );

    ImporterGmsh<typename imesh<T,Order>::type> import( fname );
    mesh->accept( import );

    return mesh;
}
template<typename T, int Order>
typename imesh<T,Order>::ptrtype
createSin( double hsize )
{
    double meshSize = hsize;
    typename imesh<T,Order>::ptrtype mesh( new typename imesh<T,Order>::type );
    typedef typename imesh<T,Order>::type mesh_type;

    //std::cout << "hsize = " << meshSize << std::endl;

    Gmsh __gmsh;
    std::string fname;
    std::ostringstream ostr;
    std::ostringstream nameStr;

    //std::cout <<"Mesh generation ... ";
    ostr << "Mesh.MshFileVersion = 2;\n"
         << "h=" << meshSize << ";\n"
         << "Include \"test_integration_ho_geometry.geo\";\n";



    nameStr << "test_integration_ho";
    __gmsh.setOrder( Order );
    fname = __gmsh.generate( nameStr.str(), ostr.str(), true );

    ImporterGmsh<typename imesh<T, Order>::type> import( fname, "2.0" );
    mesh->accept( import );

    mesh->components().set( MESH_RENUMBER | MESH_UPDATE_FACES | MESH_UPDATE_EDGES );
    mesh->updateForUse();
    return mesh;
}

}

template<int Order, typename value_type = double>
struct test_integration_circle
{
    test_integration_circle( double meshSize_=DEFAULT_MESH_SIZE ): meshSize(meshSize_) {}
    void operator()()
    {
        using namespace Life;
        using namespace Life::vf;


        double t = 0.0;
        AUTO( mycst, cst_ref( t ) );
        typename imesh<value_type,Order>::ptrtype mesh( createCircle<value_type,Order>( meshSize ) );

        t = 1.0;
        value_type v0 = integrate( elements(*mesh), IM<2,Order-1,value_type,Simplex>(), mycst ).evaluate()( 0, 0 );
        std::cout.setf( std::ios::scientific );
        std::cout.precision( 15 );
        std::cout << "v0=" << v0 << "\n";
        value_type v00 = ( integrate( boundaryelements(*mesh), IM<2,Order,value_type,Simplex>(), mycst ).evaluate()( 0, 0 )+
                           integrate( internalelements(*mesh), IM<2,Order,value_type,Simplex>(), mycst ).evaluate()( 0, 0 ) );
        std::cout << "v00=" << v00 << "\n";
        std::cout << "[circle] v0 0 = " << integrate( boundaryfaces(*mesh), IM<2,Order,value_type,Simplex>(), N() ).evaluate() << "\n";
        std::cout << "[circle] v0 0 = " << integrate( boundaryfaces(*mesh), IM<2,Order,value_type,Simplex>(), N() ).evaluate() << "\n";
        std::cout << "[circle] v0 1 = " << integrate( boundaryfaces(*mesh), IM<2,Order,value_type,Simplex>(), vec( idf(f_Nx()), idf(f_Ny() ) ) ).evaluate() << "\n";
        std::cout << "[circle] v0 2 = " << integrate( boundaryfaces(*mesh), IM<2,Order,value_type,Simplex>(),
                                                        trans(vec(constant(1),Ny()))*one() ).evaluate() << "\n";

        std::cout << "[circle] v0 3 = " << integrate( boundaryfaces(*mesh), IM<2,Order,value_type,Simplex>(),
                                                      mat<2,2>(Nx(),constant(1),constant(1),Ny())*vec(constant(1),constant(1)) ).evaluate() << "\n";
        std::cout << "[circle] v0 4 (==v0 3) = " << integrate( boundaryfaces(*mesh), IM<2,Order,value_type,Simplex>(),
                                                               mat<2,2>( idf(f_Nx()),constant(1),
                                                                         constant(1),idf(f_Ny()) )*vec(constant(1),constant(1)) ).evaluate() << "\n";
        value_type pi = 4.0*math::atan(1.0);
        const value_type eps = 1000*Life::type_traits<value_type>::epsilon();
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v0-pi, 1e-2 );
        BOOST_CHECK_SMALL( v0-v00, eps  );
#else
        LIFE_ASSERT( math::abs( v0-pi) < 1e-2 )( v0 )( math::abs( v0-pi) )( 1e-2 ).warn ( "v0 != pi" );
        LIFE_ASSERT( math::abs( v0-v00) < 1e-2 )( v0 )( v00 )( math::abs( v0-v00) )( eps ).warn ( "v0 != pi" );
#endif /* USE_BOOST_TEST */

        typedef typename imesh<value_type,Order>::type mesh_type;
        typedef fusion::vector<fem::Lagrange<2, Order+1, Scalar, Continuous, double> > basis_type;
        typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
        boost::shared_ptr<space_type> Xh( new space_type(mesh) );
        typename space_type::element_type u( Xh );

        // int ([-1,1],[-1,x]) 1 dx
        u = project( Xh, elements(*mesh), constant(1.0) );
        v0 = integrate( elements(*mesh), IM<2,(Order-1)*(Order-1),value_type,Simplex>(), idv( u ) ).evaluate()( 0, 0 );
        std::cout << "int( 1 )=" << v0 << "  (should be pi)\n";
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v0-pi, math::pow(meshSize,2*Order) );
#else
        LIFE_ASSERT( math::abs( v0-pi) < math::pow(meshSize,2*Order) )( v0 )( math::abs( v0-pi) )( math::pow(meshSize,2*Order) ).warn ( "v0 != pi" );
#endif /* USE_BOOST_TEST */

        u = project( Xh, elements(*mesh), Px()+Py()+1 );
        node_type pt(2);
        pt[0] = 0.111;
        pt[1] = 0.111;

        node_type pt2(2);
        pt2[0] = 0.111;
        pt2[1] = math::sqrt(1-pt2[0]*pt2[0]-0.001);

        node_type pt3(2);
        pt3[0] = 0;
        pt3[1] = 1;

        u = project( Xh, elements(*mesh), Px() );
        std::cout << "error x=" << integrate( elements(mesh), IM<2,8*Order,value_type,Simplex>(), (idv(u)-Px())*(idv(u)-Px()) ).evaluate()( 0, 0 ) << "\n";

        u = project( Xh, elements(*mesh), Py() );
        std::cout << "error y=" << integrate( elements(mesh), IM<2,8*Order,value_type,Simplex>(), (idv(u)-Py())*(idv(u)-Py()) ).evaluate()( 0, 0 ) << "\n";

        u = project( Xh, elements(*mesh), Px()*Px() );
        std::cout << "error x^2=" << integrate( elements(mesh), IM<2,8*Order,value_type,Simplex>(), (idv(u)-Px()*Px())*(idv(u)-Px()*Px()) ).evaluate()( 0, 0 ) << "\n";

        u = project( Xh, elements(*mesh), Px()*Py() );
        std::cout << "error xy=" << integrate( elements(mesh), IM<2,8*Order,value_type,Simplex>(), (idv(u)-Px()*Py())*(idv(u)-Px()*Py()) ).evaluate()( 0, 0 ) << "\n";

        u = project( Xh, elements(*mesh), Py()*Py() );
        std::cout << "error y^2=" << integrate( elements(mesh), IM<2,8*Order,value_type,Simplex>(), (idv(u)-Py()*Py())*(idv(u)-Py()*Py()) ).evaluate()( 0, 0 ) << "\n";

#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( (u(pt)(0,0,0)-(1+pt[0]+pt[1])), math::pow(meshSize,2*Order) );
#else
        std::cout << "u(" << pt << ")-(1+pt[0]+pt[1]) = " << math::abs( u(pt)(0,0,0)-(1+pt[0]+pt[1]) ) << "\n";
        std::cout << "u(" << pt2 << ")-(1+pt2[0]+pt2[1]) = " << math::abs( u(pt2)(0,0,0)-(1+pt2[0]+pt2[1]) ) << "\n";
        std::cout << "u(" << pt3 << ")-(1+pt3[0]+pt3[1]) = " << math::abs( u(pt3)(0,0,0)-(1+pt3[0]+pt3[1]) ) << "\n";
#endif
        v0 = integrate( elements(*mesh), IM<2,Order-1,value_type,Simplex>(), gradv( u )*vec(constant(1.0),constant(0.)) ).evaluate()( 0, 0 );
        std::cout << "int( dx(Px()) )=" << v0 << " (should be pi)\n";
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v0-pi, math::pow(meshSize,2*Order) );
#else
        LIFE_ASSERT( math::abs( v0-pi) < math::pow(meshSize,2*Order) )( v0 )( math::abs( v0-pi) )( math::pow(meshSize,2*Order) ).warn ( "v0 != pi" );
#endif /* USE_BOOST_TEST */


    }
    double meshSize;
};

template<int Order, typename value_type = double>
struct test_integration_sin
{
    test_integration_sin( double meshSize_=DEFAULT_MESH_SIZE ): meshSize(meshSize_) {}
    void operator()()
    {
        std::cout.setf( std::ios::scientific );
        std::cout.precision( 15 );

        using namespace Life;
        using namespace Life::vf;

        typedef typename imesh<value_type,Order>::type mesh_type;

        typename imesh<value_type,Order>::ptrtype mesh( createSin<value_type,Order>( meshSize ) );

        typedef fusion::vector<fem::Lagrange<2, Order+1, Scalar, Continuous, double> > basis_type;
        typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
        boost::shared_ptr<space_type> Xh( new space_type(mesh) );
        typename space_type::element_type u( Xh );

        // int ([-1,1],[-1,x]) 1 dx
        u = project( Xh, elements(*mesh), constant(1.0) );
        double v0 = integrate( elements(*mesh), IM<2,(Order-1)*(Order-1),value_type,Simplex>(), idv( u ) ).evaluate()( 0, 0 );
        //double v0 = integrate( elements(*mesh), IM<2,2*Order,value_type,Simplex>(), idv( u ) ).evaluate()( 0, 0 );
        std::cout << "int( 1 )=" << v0 << "  (should be 1+1/2*pi = 1.15915494309190)\n";


    }
    double meshSize;
};

inline
Life::po::options_description
makeOptions()
{
    Life::po::options_description integrationoptions("Test Integration options");
    integrationoptions.add_options()
        ("hsize", Life::po::value<double>()->default_value( 0.3 ), "h value")
        ("order", Life::po::value<int>()->default_value( 1 ), "geometry order ")
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
                           "integration tests",
                           Life::AboutData::License_GPL,
                           "Copyright (C) 2006,2007 Université Joseph Fourier (Grenoble I)");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}

#if defined(USE_BOOST_TEST)
boost::shared_ptr<Life::Application> mpi;
test_suite*
init_unit_test_suite( int argc, char** argv )
{
    //boost::mpi::environment( argc, argv );
    mpi = boost::shared_ptr<Life::Application>( new Life::Application( argc, argv, makeAbout(), makeOptions() ) );
    Life::Assert::setLog( "test_integration.assert");
    test_suite* test = BOOST_TEST_SUITE( "2D Generic finite element solver test suite" );

    test->add( BOOST_TEST_CASE( ( test_integration_circle<double>( mpi->vm()["hsize"].as<double>() ) ) ) );
    return test;
}
#else
int
main( int argc, char** argv )
{
    Life::Application mpi( argc, argv, makeAbout(), makeOptions() );
    Life::Assert::setLog( "test_integration.assert");


    std::cout << "Order = " << mpi.vm()["order"].as<int>()<< "\n";

    if ( mpi.vm()["order"].as<int>() == 1 )
        {
            test_integration_sin<1,double> t1( mpi.vm()["hsize"].as<double>() ); t1();
        }
#if LIFE_MESH_MAX_ORDER >= 2
    else if ( mpi.vm()["order"].as<int>() == 2 )
        {
            test_integration_sin<2,double> t2( mpi.vm()["hsize"].as<double>() ); t2();
        }
#elif LIFE_MESH_MAX_ORDER >= 3
    else if ( mpi.vm()["order"].as<int>() == 3 )
        {
            test_integration_sin<3,double> t2( mpi.vm()["hsize"].as<double>() ); t2();
        }
#elif LIFE_MESH_MAX_ORDER >= 4
    else if ( mpi.vm()["order"].as<int>() == 4 )
        {
            test_integration_sin<4,double> t2( mpi.vm()["hsize"].as<double>() ); t2();
        }
#elif LIFE_MESH_MAX_ORDER >= 5
    else if ( mpi.vm()["order"].as<int>() == 5 )
        {
            test_integration_sin<5,double> t2( mpi.vm()["hsize"].as<double>() ); t2();
        }
#endif
}
#endif // USE_BOOST_TEST
