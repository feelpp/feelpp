/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-08-25

  Copyright (C) 2006 EPFL
  Copyright (C) 2006-2012 Université Joseph Fourier (Grenoble I)

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
   \file test_integration_ho.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-08-25
 */

#define USE_BOOST_TEST 1
#define BOOST_TEST_MAIN
#if 0
// Boost.Test

#if defined(USE_BOOST_TEST)
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;
#include <boost/test/floating_point_comparison.hpp>
#endif
#else
#include <testsuite/testsuite.hpp>
#endif

#include <boost/preprocessor/comparison/greater_equal.hpp>
#include <feel/options.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feelmesh/geoentity.hpp>
#include <feel/feelmesh/refentity.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/importergmsh.hpp>
#include <feel/feelfilters/straightenmesh.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/gmshsimplexdomain.hpp>
#include <feel/feelvf/vf.hpp>

const double DEFAULT_MESH_SIZE=0.1;

namespace Feel
{
struct f_Px
{
    static const size_type context = vm::JACOBIAN|vm::POINT;
    typedef double value_type;
    typedef value_type evaluate_type;
    static const uint16_type rank = 0;
    static const uint16_type imorder = 1;
    static const bool imIsPoly = true;
    double operator()( uint16_type, uint16_type, ublas::vector<double> const& x, ublas::vector<double> const& /*n*/ ) const
    {
        return x[0];
    }
};
struct f_Nx
{
    static const size_type context = vm::JACOBIAN|vm::POINT|vm::NORMAL;
    typedef double value_type;
    typedef value_type evaluate_type;
    static const uint16_type rank = 0;
    static const uint16_type imorder = 1;
    static const bool imIsPoly = true;
    double operator()( uint16_type, uint16_type, ublas::vector<double> const& /*x*/, ublas::vector<double> const& n ) const
    {
        return n[0];
    }
};
struct f_Ny
{
    static const size_type context = vm::JACOBIAN|vm::POINT|vm::NORMAL;
    typedef double value_type;
    typedef value_type evaluate_type;
    static const uint16_type rank = 0;
    static const uint16_type imorder = 1;
    static const bool imIsPoly = true;
    double operator()( uint16_type, uint16_type, ublas::vector<double> const& /*x*/, ublas::vector<double> const& n ) const
    {
        return n[1];
    }
};
struct f_sinPx
{
    static const size_type context = vm::JACOBIAN|vm::POINT;
    typedef double value_type;
    typedef value_type evaluate_type;
    static const uint16_type rank = 0;
    static const uint16_type imorder = 2;
    static const bool imIsPoly = false;
    double operator()( uint16_type, uint16_type, ublas::vector<double> const& x, ublas::vector<double> const& /*n*/ ) const
    {
        return math::sin( x[0] );
    }
};

template<typename T, int Dim = 2,int Order = 1>
struct imesh
{
    typedef Mesh<Simplex<Dim, Order>, T > type;
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
    fname = __gmsh.generate( nameStr.str(), ostr.str() ).get<0>();
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
    typename imesh<T,2,Order>::ptrtype mesh( new typename imesh<T,2,Order>::type );
    typedef typename imesh<T,2,Order>::type mesh_type;
#if 1

    //std::cout << "hsize = " << meshSize << std::endl;

    Gmsh __gmsh( 2,Order );
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
    fname = __gmsh.generate( nameStr.str(), ostr.str() ).get<0>();



    ImporterGmsh<typename imesh<T, 2, Order>::type> import( fname );
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

    __nd[0] = math::sqrt( 2. )/2;
    __nd[1] = math::sqrt( 2. )/2;
    point_type __pt3( 3,__nd, true );
    __pt3.marker() = 0;
    mesh->addPoint( __pt3 );

    __nd[0] = -1;
    __nd[1] = 0;
    point_type __pt4( 4,__nd, true );
    __pt4.marker() = 0;
    mesh->addPoint( __pt4 );

    __nd[0] = math::sqrt( 2. )/2;;
    __nd[1] = -math::sqrt( 2. )/2;;
    point_type __pt5( 5,__nd, true );
    __pt5.marker() = 0;
    mesh->addPoint( __pt5 );

    element_type * pf = new element_type;

    pf->setId( 0 );
    pf->setMarker( 0 );

    // Warning : Vtk orientation is not the same as Feel orientation !

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
typename imesh<T,2,Order>::ptrtype
createSimplex( double hsize )
{
    double meshSize = hsize;
    //std::cout << "hsize = " << meshSize << std::endl;

    GmshSimplexDomain ts( 2,Order );
    ts.setCharacteristicLength( meshSize );
    ts.setOrder( Order );
    std::string fname = ts.generate( "simplex" );
    typename imesh<T,2,Order>::ptrtype mesh( new typename imesh<T,2,Order>::type );

    ImporterGmsh<typename imesh<T,2,Order>::type> import( fname );
    mesh->accept( import );

    return mesh;
}
template<typename T, int Order>
typename imesh<T,2,Order>::ptrtype
createSin( double hsize )
{
    double meshSize = hsize;
    typename imesh<T,2,Order>::ptrtype mesh( new typename imesh<T,2,Order>::type );
    typedef typename imesh<T,2,Order>::type mesh_type;

    //std::cout << "hsize = " << meshSize << std::endl;

    Gmsh __gmsh( 2,Order );
    std::string fname;
    std::ostringstream ostr;
    std::ostringstream nameStr;

    //std::cout <<"Mesh generation ... ";
    ostr << "Mesh.MshFileVersion = 2;\n"
         << "Mesh.ElementOrder=" << Order << ";\n"
         << "Mesh.SecondOrderIncomplete = 0;\n"
         << "h=" << meshSize << ";\n"
         << "Include \"test_integration_ho_geometry.geo\";\n";



    nameStr << "test_integration_ho";
    __gmsh.setOrder( Order );
    fname = __gmsh.generate( nameStr.str(), ostr.str(), true ).get<0>();

    ImporterGmsh<typename imesh<T, 2,Order>::type> import( fname, "2.0" );
    mesh->accept( import );

    mesh->components().set( MESH_RENUMBER | MESH_UPDATE_FACES | MESH_UPDATE_EDGES );
    mesh->updateForUse();
    return mesh;
}

}

template<int Order, typename value_type = double>
struct test_integration_circle
{
    test_integration_circle( double meshSize_=DEFAULT_MESH_SIZE ): meshSize( meshSize_ ) {}
    void operator()()
    {
        using namespace Feel;
        using namespace Feel::vf;


        double t = 0.0;
        AUTO( mycst, cst_ref( t ) );
        typename imesh<value_type,2,Order>::ptrtype mesh( createCircle<value_type,Order>( meshSize ) );

        t = 1.0;
        value_type v0 = integrate( elements( mesh ), mycst, _Q<Order-1>() ).evaluate()( 0, 0 );
        std::cout.setf( std::ios::scientific );
        std::cout.precision( 15 );
        std::cout << "v0=" << v0 << "\n";
        value_type v00 = ( integrate( boundaryelements( mesh ), mycst ).evaluate()( 0, 0 )+
                           integrate( internalelements( mesh ), mycst ).evaluate()( 0, 0 ) );
        std::cout << "v00=" << v00 << "\n";
        std::cout << "[circle] v0 0 = " << integrate( boundaryfaces( mesh ), N() ).evaluate() << "\n";
        std::cout << "[circle] v0 0 = " << integrate( boundaryfaces( mesh ),  N() ).evaluate() << "\n";
        std::cout << "[circle] v0 1 = " << integrate( boundaryfaces( mesh ),  vec( idf( f_Nx() ), idf( f_Ny() ) ) ).evaluate() << "\n";
        std::cout << "[circle] v0 2 = " << integrate( boundaryfaces( mesh ),
                  trans( vec( constant( 1 ),Ny() ) )*one() ).evaluate() << "\n";

        std::cout << "[circle] v0 3 = " << integrate( boundaryfaces( mesh ),
                  mat<2,2>( Nx(),constant( 1 ),constant( 1 ),Ny() )*vec( constant( 1 ),constant( 1 ) ) ).evaluate() << "\n";
        std::cout << "[circle] v0 4 (==v0 3) = " << integrate( boundaryfaces( mesh ),
                  mat<2,2>( idf( f_Nx() ),constant( 1 ),
                            constant( 1 ),idf( f_Ny() ) )*vec( constant( 1 ),constant( 1 ) ) ).evaluate() << "\n";
        value_type pi = 4.0*math::atan( 1.0 );
        const value_type eps = 1000*Feel::type_traits<value_type>::epsilon();
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v0-pi, 1e-2 );
        BOOST_CHECK_SMALL( v0-v00, eps  );
#else
        FEELPP_ASSERT( math::abs( v0-pi ) < 1e-2 )( v0 )( math::abs( v0-pi ) )( 1e-2 ).warn ( "v0 != pi" );
        FEELPP_ASSERT( math::abs( v0-v00 ) < 1e-2 )( v0 )( v00 )( math::abs( v0-v00 ) )( eps ).warn ( "v0 != pi" );
#endif /* USE_BOOST_TEST */

        typedef typename imesh<value_type,2,Order>::type mesh_type;
        typedef fusion::vector<Lagrange<Order+1, Scalar> > basis_type;
        typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
        boost::shared_ptr<space_type> Xh( new space_type( mesh ) );
        typename space_type::element_type u( Xh );

        // int ([-1,1],[-1,x]) 1 dx
        u = vf::project( Xh, elements( mesh ), constant( 1.0 ) );
        v0 = integrate( elements( mesh ), idv( u ) ).evaluate()( 0, 0 );
        std::cout << "int( 1 )=" << v0 << "  (should be pi)\n";
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v0-pi, math::pow( meshSize,2*Order ) );
#else
        FEELPP_ASSERT( math::abs( v0-pi ) < math::pow( meshSize,2*Order ) )( v0 )( math::abs( v0-pi ) )( math::pow( meshSize,2*Order ) ).warn ( "v0 != pi" );
#endif /* USE_BOOST_TEST */

        u = vf::project( Xh, elements( mesh ), Px()+Py()+1 );
        node_type pt( 2 );
        pt[0] = 0.111;
        pt[1] = 0.111;

        node_type pt2( 2 );
        pt2[0] = 0.111;
        pt2[1] = math::sqrt( 1-pt2[0]*pt2[0]-0.001 );

        node_type pt3( 2 );
        pt3[0] = 0;
        pt3[1] = 1;

        u = vf::project( Xh, elements( mesh ), Px() );
        std::cout << "error x=" << integrate( elements( mesh ), ( idv( u )-Px() )*( idv( u )-Px() ) ).evaluate()( 0, 0 ) << "\n";

        u = vf::project( Xh, elements( mesh ), Py() );
        std::cout << "error y=" << integrate( elements( mesh ), ( idv( u )-Py() )*( idv( u )-Py() ) ).evaluate()( 0, 0 ) << "\n";

        u = vf::project( Xh, elements( mesh ), Px()*Px() );
        std::cout << "error x^2=" << integrate( elements( mesh ), ( idv( u )-Px()*Px() )*( idv( u )-Px()*Px() ) ).evaluate()( 0, 0 ) << "\n";

        u = vf::project( Xh, elements( mesh ), Px()*Py() );
        std::cout << "error xy=" << integrate( elements( mesh ), ( idv( u )-Px()*Py() )*( idv( u )-Px()*Py() ) ).evaluate()( 0, 0 ) << "\n";

        u = vf::project( Xh, elements( mesh ), Py()*Py() );
        std::cout << "error y^2=" << integrate( elements( mesh ), ( idv( u )-Py()*Py() )*( idv( u )-Py()*Py() ) ).evaluate()( 0, 0 ) << "\n";

#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( ( u( pt )( 0,0,0 )-( 1+pt[0]+pt[1] ) ), math::pow( meshSize,2*Order ) );
#else
        std::cout << "u(" << pt << ")-(1+pt[0]+pt[1]) = " << math::abs( u( pt )( 0,0,0 )-( 1+pt[0]+pt[1] ) ) << "\n";
        std::cout << "u(" << pt2 << ")-(1+pt2[0]+pt2[1]) = " << math::abs( u( pt2 )( 0,0,0 )-( 1+pt2[0]+pt2[1] ) ) << "\n";
        std::cout << "u(" << pt3 << ")-(1+pt3[0]+pt3[1]) = " << math::abs( u( pt3 )( 0,0,0 )-( 1+pt3[0]+pt3[1] ) ) << "\n";
#endif
        v0 = integrate( elements( mesh ), gradv( u )*vec( constant( 1.0 ),constant( 0. ) ) ).evaluate()( 0, 0 );
        std::cout << "int( dx(Px()) )=" << v0 << " (should be pi)\n";
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v0-pi, math::pow( meshSize,2*Order ) );
#else
        FEELPP_ASSERT( math::abs( v0-pi ) < math::pow( meshSize,2*Order ) )( v0 )( math::abs( v0-pi ) )( math::pow( meshSize,2*Order ) ).warn ( "v0 != pi" );
#endif /* USE_BOOST_TEST */


    }
    double meshSize;
};

template<int Dim, int Order, typename value_type = double>
struct test_integration_sin
{
    test_integration_sin( Feel::po::variables_map const& _vm ): vm( _vm ) {}
    void operator()()
    {
        std::cout.setf( std::ios::scientific );
        std::cout.precision( 15 );

        using namespace Feel;
        using namespace Feel::vf;
        double meshSize = vm["hsize"].template as<double>();
        GeomapStrategyType geomap = ( GeomapStrategyType )vm["geomap"].template as<int>();
        int straighten = vm["straighten"].template as<int>();
        typedef typename imesh<value_type,Dim,Order>::type mesh_type;

        //typename imesh<value_type,Order>::ptrtype mesh( createSin<value_type,Order>( meshSize ) );
        auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                    _desc=domain( _name=( boost::format( "ellipsoid-%1%-%2%-%3%" ) % Dim % Order % straighten ).str() ,
                                            _usenames=true,
                                            _shape="ellipsoid",
                                            _xmin=-1,_xmax=1,
                                            _ymin=-1,_ymax=1,
                                            _zmin=-1,_zmax=1,
                                            _dim=Dim,
                                            _order=Order,
                                            _h=meshSize,
                                            _straighten=straighten  ) );

        typedef fusion::vector<Lagrange<1, Scalar> > basis_type;
        typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
        boost::shared_ptr<space_type> Xh( new space_type( mesh ) );
        typename space_type::element_type u( Xh );

        // int ([-1,1],[-1,x]) 1 dx
        u = vf::project( Xh, elements( mesh ), constant( 1.0 ) );
        double v3 = integrate( _range=elements( mesh ), _expr=cst( 1.0 ), _quad=_Q<15>(), _geomap=geomap  ).evaluate()( 0, 0 );
        double v0 = integrate( _range=elements( mesh ), _expr=idv( u ), _quad=_Q<15>(), _geomap=geomap ).evaluate()( 0, 0 );
        double v2 = integrate( _range=boundaryfaces( mesh ), _expr=idv( u ), _quad=_Q<15>(), _geomap=geomap  ).evaluate()( 0, 0 );
        //double v0 = integrate( elements(mesh), idv( u ) ).evaluate()( 0, 0 );
        std::cout << "int( 1 )=" << v3 << "  (=pi) error= " << math::abs( v3 - M_PI ) << std::endl;
        std::cout << "int( u=1 )=" << v0 << "  (=pi) error= " << math::abs( v0 - M_PI ) << std::endl;
        std::cout << "int(boundary, 1 )=" << v2 << "  (=2*pi) error= " << math::abs( v2 - 2*M_PI ) << std::endl;
        double v1 = integrate( _range=boundaryfaces( mesh ), _expr=trans( vec( cst( 1. ),cst( 1. ) ) )*N(), _quad=_Q<( Order-1 )*( Order-1 )>()  ).evaluate()( 0, 0 );
        std::cout << "int( 1 .N )=" << v1 << "  (=pi) error= " << math::abs( v1 ) << std::endl;
        BOOST_CHECK_SMALL( v1, 1e-12 );
        //BOOST_CHECK_SMALL( math::abs( v0 - M_PI ), 3*exp(-4*Order));
    }
    Feel::po::variables_map vm;
};

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description integrationoptions( "Test Integration options" );
    integrationoptions.add_options()
    ( "hsize", Feel::po::value<double>()->default_value( 0.3 ), "h value" )
    ( "order", Feel::po::value<int>()->default_value( 1 ), "geometry order " )
    ( "straighten", Feel::po::value<int>()->default_value( 0 ), "straighten mesh " )
    ( "geomap", Feel::po::value<int>()->default_value( 2 ), "high order integration (0=opt, 1=p1, 2=ho " )
    ;
    return integrationoptions.add( Feel::feel_options() );
}

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "test_integration_ho" ,
                           "test_integration_ho" ,
                           "0.1",
                           "integration tests",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2006,2007 Université Joseph Fourier (Grenoble I)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

#if defined(USE_BOOST_TEST)
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )
BOOST_AUTO_TEST_SUITE( integration )

//typedef boost::mpl::list<boost::mpl::int_<1>,boost::mpl::int_<2>,boost::mpl::int_<3>,boost::mpl::int_<4>,boost::mpl::int_<5>  > order_types;
typedef boost::mpl::list<boost::mpl::int_<1>,boost::mpl::int_<2>,boost::mpl::int_<3>,boost::mpl::int_<4>> order_types;
//typedef boost::mpl::list<boost::mpl::int_<1>  > order_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( test_integration_ho, T, order_types )
{
    std::cout << "============================================================\n";
    Feel::Application mpi;
    Feel::Assert::setLog( "test_integration_ho.assert" );
    std::cout << "Order=" << T::value << "/5,"
              << " hsize=" << mpi.vm()["hsize"].as<double>() << ","
              << " straighten=" << mpi.vm()["straighten"].as<int>() << ","
              << " geomap=" << mpi.vm()["geomap"].as<int>() << "\n";

    test_integration_sin<2,T::value,double> t(  mpi.vm() );
    t();
}
BOOST_AUTO_TEST_SUITE_END()

#else

int
main( int argc, char** argv )
{
    Feel::Application mpi( argc, argv, makeAbout(), makeOptions() );
    Feel::Assert::setLog( "test_integration.assert" );

    std::cout << "Order = " << mpi.vm()["order"].as<int>() << " / " << FEELPP_MESH_MAX_ORDER << "\n";

    if ( mpi.vm()["order"].as<int>() == 1 )
    {
        test_integration_sin<2,1,double> t1( mpi.vm()["hsize"].as<double>() );
        t1();
    }

    else if ( mpi.vm()["order"].as<int>() == 2 )
    {
#if BOOST_PP_GREATER_EQUAL(FEELPP_MESH_MAX_ORDER, 2)
        test_integration_sin<2,2,double> t2( mpi.vm()["hsize"].as<double>() );
        t2();
#endif
    }

    else if ( mpi.vm()["order"].as<int>() == 3 )
    {
#if BOOST_PP_GREATER_EQUAL(FEELPP_MESH_MAX_ORDER, 3)
        test_integration_sin<2,3,double> t2( mpi.vm()["hsize"].as<double>() );
        t2();
#endif
    }

    else if ( mpi.vm()["order"].as<int>() == 4 )
    {
#if BOOST_PP_GREATER_EQUAL(FEELPP_MESH_MAX_ORDER, 4)
        test_integration_sin<2,4,double> t2( mpi.vm()["hsize"].as<double>() );
        t2();
#endif
    }

    // disable 
#if 0
    else if ( mpi.vm()["order"].as<int>() == 5 )
    {
#if BOOST_PP_GREATER_EQUAL(FEELPP_MESH_MAX_ORDER, 5)
        test_integration_sin<2,5,double> t2( mpi.vm()["hsize"].as<double>() );
        t2();
#endif
    }
#endif

}
#endif // USE_BOOST_TEST
