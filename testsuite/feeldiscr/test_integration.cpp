/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2006-08-25

  Copyright (C) 2006 EPFL
  Copyright (C) 2006-2009 Université Joseph Fourier (Grenoble I)

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

// make sure that the init_unit_test function is defined by UTF
#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE 3D integration testsuite
// disable the main function creation, use our own
#define BOOST_TEST_NO_MAIN

#if defined(USE_BOOST_TEST)
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;
#include <boost/test/floating_point_comparison.hpp>
#endif

#include <feel/options.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/feelalg/backendgmm.hpp>
#include <feel/feelmesh/geoentity.hpp>
#include <feel/feelmesh/refentity.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/gmshsimplexdomain.hpp>
#include <feel/feelvf/vf.hpp>

const double DEFAULT_MESH_SIZE=0.05;

namespace Feel
{
struct f_Px
{
    static const size_type context = vm::JACOBIAN|vm::POINT;
    typedef double value_type;
    typedef Feel::uint16_type uint16_type;
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
    typedef Feel::uint16_type uint16_type;
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
    typedef Feel::uint16_type uint16_type;
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
    typedef Feel::uint16_type uint16_type;
    static const uint16_type rank = 0;
    double operator()( uint16_type, uint16_type, ublas::vector<double> const& x, ublas::vector<double> const& /*n*/ ) const
    {
        return math::sin( x[0] );
    }
};

template<typename T, int Order = 1>
struct imesh
{
    typedef Mesh<Simplex<2, Order>, T > type;
    typedef boost::shared_ptr<type> ptrtype;
};

template<typename T>
typename imesh<T>::ptrtype
createMesh( double hsize )
{
    double meshSize = hsize;
    //BOOST_TEST_MESSAGE( "hsize = " << meshSize << std::endl;

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
    ostr << "Mesh.MshFileVersion = 2;\n"
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
    //BOOST_TEST_MESSAGE( "hsize = " << meshSize << std::endl;

    Gmsh __gmsh;
    std::string fname;
    std::ostringstream ostr;
    std::ostringstream nameStr;

    //std::cout <<"Mesh generation ... ";
    ostr << "Mesh.MshFileVersion = 2;\n"
         << "h=" << meshSize << ";\n"
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
#else
         << "Line(1) = {1,2};\n"
         << "Line(2) = {3,1};\n"
         << "Line Loop(14) = {1,10,2};\n"
#endif

         << "Plane Surface(15) = {14};\n"
         << "Physical Line(20) = {10,11,12,13};\n"
         << "Physical Surface(30) = {15};\n";

    nameStr << "circle." << meshSize;
    if ( Order == 2 )
        __gmsh.setOrder( GMSH_ORDER_TWO );
    fname = __gmsh.generate( nameStr.str(), ostr.str() );

    typename imesh<T,Order>::ptrtype mesh( new typename imesh<T,Order>::type );

    ImporterGmsh<typename imesh<T, Order>::type> import( fname );
    mesh->accept( import );

    mesh->components().set( MESH_RENUMBER | MESH_UPDATE_FACES | MESH_UPDATE_EDGES );
    mesh->updateForUse();
    return mesh;
}
template<typename T, int Order>
typename imesh<T,Order>::ptrtype
createSimplex( double hsize )
{
    double meshSize = hsize;
    //BOOST_TEST_MESSAGE( "hsize = " << meshSize << std::endl;

    GmshSimplexDomain ts(2,Order);
    ts.setCharacteristicLength( meshSize );
    ts.setX( std::make_pair(-1,1) );
    ts.setY( std::make_pair(-1,1) );
    std::string fname = ts.generate( "simplex" );
    typename imesh<T,Order>::ptrtype mesh( new typename imesh<T>::type );

    ImporterGmsh<typename imesh<T>::type> import( fname );
    mesh->accept( import );

    return mesh;
}
}

template<typename value_type = double>
struct test_integration_circle
{
    test_integration_circle( double meshSize_=DEFAULT_MESH_SIZE ): meshSize(meshSize_) {}
    void operator()()
    {
        using namespace Feel;
        using namespace Feel::vf;

        int Order=2;
        double t = 0.0;
        AUTO( mycst, cst_ref( t ) );
        typename imesh<value_type,1>::ptrtype mesh( createCircle<value_type,1>( meshSize ) );

        t = 1.0;
        value_type v0 = integrate( elements(mesh), mycst ).evaluate()( 0, 0 );
        BOOST_TEST_MESSAGE( "v0=" << v0 << "\n" );
        value_type v00 = ( integrate( boundaryelements(mesh), IM<2,2,value_type,Simplex>(), mycst ).evaluate()( 0, 0 )+
                           integrate( internalelements(mesh), IM<2,2,value_type,Simplex>(), mycst ).evaluate()( 0, 0 ) );
        BOOST_TEST_MESSAGE( "v00=" << v00 << "\n" );
        BOOST_TEST_MESSAGE( "[circle] v0 0 = " << integrate( boundaryfaces(mesh), IM<2,2,value_type,Simplex>(), N() ).evaluate() << "\n" );
        BOOST_TEST_MESSAGE( "[circle] v0 0 = " << integrate( boundaryfaces(mesh), IM<2,2,value_type,Simplex>(), N() ).evaluate() << "\n" );
        BOOST_TEST_MESSAGE( "[circle] v0 1 = " << integrate( boundaryfaces(mesh), IM<2,2,value_type,Simplex>(), vec( idf(f_Nx()), idf(f_Ny() ) ) ).evaluate() << "\n" );
        BOOST_TEST_MESSAGE( "[circle] v0 2 = " << integrate( boundaryfaces(mesh), IM<2,2,value_type,Simplex>(),
                                                        trans(vec(constant(1),Ny()))*one() ).evaluate() << "\n" );

        BOOST_TEST_MESSAGE( "[circle] v0 3 = " << integrate( boundaryfaces(mesh), IM<2,2,value_type,Simplex>(),
                                                      mat<2,2>(Nx(),constant(1),constant(1),Ny())*vec(constant(1),constant(1)) ).evaluate() << "\n" );
        BOOST_TEST_MESSAGE( "[circle] v0 4 (==v0 3) = " << integrate( boundaryfaces(mesh), IM<2,2,value_type,Simplex>(),
                                                               mat<2,2>( idf(f_Nx()),constant(1),
                                                                         constant(1),idf(f_Ny()) )*vec(constant(1),constant(1)) ).evaluate() << "\n" );
        value_type pi = 4.0*math::atan(1.0);
        const value_type eps = 1000*Feel::type_traits<value_type>::epsilon();
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v0, pi, 2e-1 );
        BOOST_CHECK_SMALL( v0-v00, eps  );
#else
        FEEL_ASSERT( math::abs( v0-pi) < math::pow( meshSize, 2*Order ) )( v0 )( math::abs( v0-pi) )( math::pow( meshSize, 2*Order ) ).warn ( "v0 != pi" );
        FEEL_ASSERT( math::abs( v0-v00) < eps )( v0 )( v00 )( math::abs( v0-v00) )( eps ).warn ( "v0 != pi" );
#endif /* USE_BOOST_TEST */

        typedef typename imesh<value_type,1>::type mesh_type;
        typedef fusion::vector<Lagrange<2, Scalar> > basis_type;
        typedef FunctionSpace<mesh_type, basis_type> space_type;
        boost::shared_ptr<space_type> Xh( new space_type(mesh) );
        typename space_type::element_type u( Xh );

        u = vf::project( Xh, elements(mesh), constant(1.0) );
        v0 = integrate( elements(mesh), idv( u ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v0, pi, 2e-1);
#else
        FEEL_ASSERT( math::abs( v0-pi) < math::pow( meshSize, 2*Order ) )( v0 )( math::abs( v0-pi) )( math::pow( meshSize, 2*Order ) ).warn ( "v0 != pi" );
#endif /* USE_BOOST_TEST */

        typedef fusion::vector<Lagrange<2, Vectorial> > vector_basis_type;
        typedef FunctionSpace<mesh_type, vector_basis_type, value_type> vector_space_type;
        boost::shared_ptr<vector_space_type> Xvh( new vector_space_type(mesh) );
        typename vector_space_type::element_type U( Xvh );

        U = vf::project( Xvh, elements(mesh), vec(constant(1.0),constant(1.0)) );
        v0 = integrate( boundaryfaces(mesh), trans(idv( U ))*N() ).evaluate()( 0, 0 );
        v00 = integrate( elements(mesh), divv(U) ).evaluate()( 0, 0 );

        BOOST_TEST_MESSAGE( "[circle] v0=" << v0 << " v00=" << v00 << " should be equal thanks to gauss int 1.N() != int div 1\n");
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v0-v00, eps );
#else
        FEEL_ASSERT( math::abs( v0-v00) < eps )( v0 )(v00)( math::abs( v0-v00) ).warn ( "int 1.N() != int div 1" );
#endif /* USE_BOOST_TEST */

        U = vf::project( Xvh, elements(mesh), vec(Px(),Py()) );
        v0 = integrate( boundaryfaces(mesh), trans(idv( U ))*N() ).evaluate()( 0, 0 );
        v00 = integrate( elements(mesh), divv(U) ).evaluate()( 0, 0 );

        BOOST_TEST_MESSAGE( "[circle] v0=" << v0 << " v00=" << v00 << " should be equal thanks to gauss int (x,y).N() != int div (x,y)\n");
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v0, v00, eps );
#else
        FEEL_ASSERT( math::abs( v0-v00) < eps )( v0 )(v00)( math::abs( v0-v00) ).warn ( "int (x,y).N() != int div (x,y)" );
#endif /* USE_BOOST_TEST */

    }
    double meshSize;
};
template<typename value_type = double>
struct test_integration_simplex
{
    test_integration_simplex( double meshSize_=DEFAULT_MESH_SIZE ): meshSize(meshSize_) {}
    void operator()()
    {
        using namespace Feel;
        using namespace Feel::vf;


        const value_type eps = 1e-9;
        typename imesh<value_type,1>::ptrtype mesh( createSimplex<value_type,1>( meshSize ) );
        typedef typename imesh<value_type>::type mesh_type;

        typedef fusion::vector<Lagrange<3, Scalar> > basis_type;
        typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
        boost::shared_ptr<space_type> Xh( new space_type(mesh) );
        typename space_type::element_type u( Xh );

        value_type meas = integrate( elements(mesh), constant(1.0) ).evaluate()( 0, 0 );
        value_type v0 = integrate( elements(mesh), constant(1.0) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( meas, 2, eps );
#else
        const value_type eps = 1000*Feel::type_traits<value_type>::epsilon();

        FEEL_ASSERT( math::abs( v0-pi) < 1e-2 )( v0 )( math::abs( v0-pi) )( eps ).warn ( "v0 != pi" );
#endif /* USE_BOOST_TEST */

        auto in1 = integrate( boundaryfaces(mesh), N() ).evaluate();
        auto in2 = integrate( boundaryfaces(mesh), _Q<2>(), vec( idf(f_Nx()), idf(f_Ny() ) ) ).evaluate();
        BOOST_CHECK_CLOSE( in1(0,0), in2(0,0), eps );
        BOOST_CHECK_CLOSE( in1(1,0), in2(1,0), eps );
        BOOST_TEST_MESSAGE( "[simplex] v0 0 = " << in1  << "\n" );
        BOOST_TEST_MESSAGE( "[simplex] v0 1 = " << in2  << "\n" );
        BOOST_TEST_MESSAGE( "[simplex] v0 2 = " << integrate( boundaryfaces(mesh), trans(vec(constant(1),Ny()))*one() ).evaluate() << "\n" );

        auto in3 = integrate( boundaryfaces(mesh), _Q<2>(), mat<2,2>(Nx(),constant(1),constant(1),Ny())*vec(constant(1),constant(1)) ).evaluate();
        auto in4 = integrate( boundaryfaces(mesh), _Q<2>(), mat<2,2>( idf(f_Nx()),constant(1),
                                                                        constant(1),idf(f_Ny()) )*vec(constant(1),constant(1)) ).evaluate();
        BOOST_CHECK_CLOSE( in3(0,0), in4(0,0), eps );
        BOOST_CHECK_CLOSE( in3(1,0), in4(1,0), eps );

        BOOST_TEST_MESSAGE( "[simplex] v0 3 = " << in3 << "\n" );
        BOOST_TEST_MESSAGE( "[simplex] v0 4 (==v0 3) = " <<  in4 << "\n" );

        // int ([-1,1],[-1,x]) 1 dx
        u = vf::project( Xh, elements(mesh), Px()+Py() );
        double v1 = integrate( elements(mesh), gradv(u)*trans(gradv(u) ) ).evaluate()(0,0) ;
        BOOST_CHECK_CLOSE( v1, 2*meas, eps );

        u = vf::project( Xh, elements(mesh), Px()*Px() );
        double lapu = integrate( elements(mesh),  trace( hessv(u) ) ).evaluate()( 0, 0 );
        BOOST_CHECK_CLOSE( lapu, 2*meas, eps );
        double lapu1 = integrate( elements(mesh),  trace( hessv(u)*trans(hessv(u) ) )).evaluate()( 0, 0 );
        BOOST_CHECK_CLOSE( lapu1, 4*meas, eps );
        auto hessu = integrate( elements(mesh), hessv(u) ).evaluate();
        BOOST_CHECK_CLOSE( hessu(0,0), 2*meas, eps );
        BOOST_CHECK_SMALL( hessu(1,0), eps );
        BOOST_CHECK_SMALL( hessu(0,1), eps );
        BOOST_CHECK_SMALL( hessu(1,1), eps );

        u = vf::project( Xh, elements(mesh), Px()*Px()+Py()*Py() );
        lapu = integrate( elements(mesh),  trace( hessv(u) ) ).evaluate()( 0, 0 );
        BOOST_CHECK_CLOSE( lapu, 4*meas, eps );
        lapu1 = integrate( elements(mesh),  trace( hessv(u)*trans(hessv(u) ) )).evaluate()( 0, 0 );
        BOOST_CHECK_CLOSE( lapu1, 8*meas, eps );
        hessu = integrate( elements(mesh), hessv(u) ).evaluate();
        BOOST_CHECK_CLOSE( hessu(0,0), 2*meas, eps );
        BOOST_CHECK_SMALL( hessu(1,0), eps );
        BOOST_CHECK_SMALL( hessu(0,1), eps );
        BOOST_CHECK_CLOSE( hessu(1,1), 2*meas, eps );

        u = vf::project( Xh, elements(mesh), Px()*Px()+Py()*Py()+3*Px()*Py() );
        lapu = integrate( elements(mesh),  trace( hessv(u) ) ).evaluate()( 0, 0 );
        BOOST_CHECK_CLOSE( lapu, 4*meas, eps );
        lapu1 = integrate( elements(mesh),  trace( hessv(u)*trans(hessv(u) ) )).evaluate()( 0, 0 );
        BOOST_CHECK_CLOSE( lapu1, 26*meas, eps );
        hessu = integrate( elements(mesh), hessv(u) ).evaluate();
        BOOST_CHECK_CLOSE( hessu(0,0), 2*meas, eps );
        BOOST_CHECK_CLOSE( hessu(1,0), 3*meas, eps );
        BOOST_CHECK_CLOSE( hessu(0,1), 3*meas, eps );
        BOOST_CHECK_CLOSE( hessu(1,1), 2*meas, eps );


        typedef fusion::vector<Lagrange<3, Vectorial> > v_basis_type;
        typedef FunctionSpace<mesh_type, v_basis_type, value_type> v_space_type;
        boost::shared_ptr<v_space_type> Yh( new v_space_type(mesh) );
        typename v_space_type::element_type v( Yh );

        auto p2 = Px()*Px()*Py()+Py()*Py()+cos(Px());
        v = vf::project( Yh, elements(mesh), vec(p2,p2) );
        double divp2 = integrate( elements(mesh), divv(v) ).evaluate()( 0, 0 );
        double unp2 = integrate( boundaryfaces(mesh), trans(idv(v))*N() ).evaluate()( 0, 0 );
        BOOST_CHECK_CLOSE( divp2, unp2, eps );
    }
    double meshSize;
};
template<typename value_type = double>
struct test_integration_domain
{
    test_integration_domain( double meshSize_=DEFAULT_MESH_SIZE ): meshSize(meshSize_) {}
    void operator()()
    {
        using namespace Feel;
        using namespace Feel::vf;


        typename imesh<value_type>::ptrtype mesh( createMesh<value_type>( meshSize ) );

        const value_type eps = 1e-9;

        // int ([-1,1],[-1,x]) 1 dx
        value_type v0 = integrate( elements(mesh),
                                   vf::min(constant(1.0),constant(2.0)) ).evaluate()( 0, 0 );
        value_type v00 = ( integrate( boundaryelements(mesh),
                                      vf::min(constant(1.0),constant(2.0)) ).evaluate()( 0, 0 )+
                           integrate( internalelements(mesh),
                                      vf::min(constant(1.0),constant(2.0)) ).evaluate()( 0, 0 ) );

        BOOST_TEST_MESSAGE( "[domain] v0 = " << v0 << "\n" );
        BOOST_TEST_MESSAGE( "[domain] v00 = " << v00 << "\n" );
        BOOST_TEST_MESSAGE( "[domain] int(1*N()) = " << integrate( boundaryfaces(mesh),
                                                   constant(1.0)*N() ).evaluate() << "\n" );
        BOOST_TEST_MESSAGE( "[domain] int(tr(N())**N()) = " << integrate( boundaryfaces(mesh),
                                                          N()*trans(N()) ).evaluate() << "\n" );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v0, 4.0, eps );
        BOOST_CHECK_CLOSE( v0, v00, eps );
#else
        FEEL_ASSERT( math::abs( v0-4.0) < eps )( v0 )( math::abs( v0-4.0) )( eps ).warn ( "v0 != 4" );
        FEEL_ASSERT( math::abs( v0-v00) < eps )( v0 )(v00)( math::abs( v0-v00) )( eps ).warn ( "v0 != v00" );
#endif /* USE_BOOST_TEST */


        // int ([-1,1],[-1,1]) x dx
        value_type v1 = integrate( elements(mesh), Px() ).evaluate()( 0, 0 );
        value_type v11 = integrate( elements(mesh), _Q<2>(), idf(f_Px()) ).evaluate()( 0, 0 );

        BOOST_TEST_MESSAGE( "[domain] int(P()) = " << integrate( elements(mesh),
                                                 P() ).evaluate() << "\n" );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v1- 0.0, eps );
        BOOST_CHECK_SMALL( v11- 0.0, eps );
#else
        FEEL_ASSERT( math::abs( v1-0.0) < eps )( v1 )( math::abs( v1-0.0) )( eps ).warn ( "v1 != 0" );
        FEEL_ASSERT( math::abs( v11-0.0) < eps )( v11 )( math::abs( v11-0.0) )( eps ).warn ( "v11 != 0" );
#endif /* USE_BOOST_TEST */

         // int ([-1,1],[-1,1]) abs(x) dx
        value_type vsin = integrate( elements(mesh), _Q<5>(), sin(Px()) ).evaluate()( 0, 0 );
        value_type vsin1 = integrate( elements(mesh), _Q<5>(), idf(f_sinPx()) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( vsin- 2*(-math::cos(1.0)+math::cos(-1.0)), eps );
        BOOST_CHECK_SMALL( vsin1- 2*(-math::cos(1.0)+math::cos(-1.0)), eps );
#else
        FEEL_ASSERT( math::abs( vsin-2.0*(-math::cos(1.0)+math::cos(-1.0))) < eps )
            ( vsin )
            ( math::abs( vsin-2.0*(-math::cos(1.0)+math::cos(-1.0))) )( eps ).warn ( "vsin != 2*(cos(1)-cos(-1))" );
        FEEL_ASSERT( math::abs( vsin1-2.0*(-math::cos(1.0)+math::cos(-1.0))) < eps )
            ( vsin1 )
            ( math::abs( vsin1-2.0*(-math::cos(1.0)+math::cos(-1.0))) )( eps ).warn ( "vsin1 != 2*(cos(1)-cos(-1))" );
#endif /* USE_BOOST_TEST */

        // int ([-1,1],[-1,1]) abs(x) dx
        value_type vabs = integrate( elements(mesh), _Q<5>(), abs(Px())+abs(Py()) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( vabs, 4.0, 1e-2 );
#else
        FEEL_ASSERT( math::abs( vabs-4.0) < 1e-2 )( vabs )( math::abs( vabs-4.0) )( 1e-2 ).warn ( "vabs != 4" );
#endif /* USE_BOOST_TEST */



        // int (\partial ([-1,1],[-1,1]) 1 dx
        value_type v2 = integrate( boundaryfaces(mesh), constant(1.0) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v2, 8.0, eps );
#else
        FEEL_ASSERT( math::abs( v2-8.0) < eps )( v2 )( math::abs( v2-8.0) )( eps ).warn ( "v2 != 8" );
#endif /* USE_BOOST_TEST */


    }
    double meshSize;
};

template<typename value_type = double>
struct test_integration_boundary
{
    test_integration_boundary( double meshSize_=DEFAULT_MESH_SIZE ): meshSize(meshSize_) {}
    void operator()()
    {
        using namespace Feel;
        using namespace Feel::vf;


        typename imesh<value_type>::ptrtype mesh( createMesh<value_type>( meshSize ) );

        const value_type eps = 1000*Feel::type_traits<value_type>::epsilon();

        // int (\partial ([-1,1],[-1,1]) 1 dx
        value_type v2 = integrate( boundaryfaces(mesh), constant(1.0) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v2, 8.0, eps );
#else
        FEEL_ASSERT( math::abs( v2-8.0) < eps )( v2 )( math::abs( v2-8.0) )( eps ).warn ( "v2 != 8" );
#endif /* USE_BOOST_TEST */

        // int (\partial ([-1,1],[-1,1]) x * y dx dy
        value_type v3 = integrate( boundaryfaces(mesh), Px()*Py() ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v3- 0.0, eps );
#else
        FEEL_ASSERT( math::abs( v3-0.0) < eps )( v3 )( math::abs( v3-0.0) )( eps ).warn ( "v3 != 0" );
#endif /* USE_BOOST_TEST */

        // int (\partial ([-1,1],[-1,1]) -2*y dx dy
        value_type v4 = integrate( boundaryfaces(mesh), -2*Py() ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v4- 0.0, eps );
#else
        FEEL_ASSERT( math::abs( v4-0.0) < eps )( v4 )( math::abs( v4-0.0) )( eps ).warn ( "v4 != 0" );
#endif /* USE_BOOST_TEST */

        // int_10 exp(Py) dx
        value_type v5 = integrate( markedfaces(*mesh,10), exp(Py()) ).evaluate()( 0, 0 );
        value_type v5_ex = 2.0*(math::exp(1.0)+math::exp(-1.0));
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v5, v5_ex, eps );
#else
        FEEL_ASSERT( math::abs( v5-v5_ex) < eps )( v5 )(v5_ex)( math::abs( v5-v5_ex) )( eps ).warn ( "v5 != v5_ex" );
#endif /* USE_BOOST_TEST */

        // int_20 exp(Py) dx
        value_type v6 = integrate( markedfaces(*mesh,20), _Q<5>(), exp(Py()) ).evaluate()( 0, 0 );
        value_type v6_ex = 2.0*(math::exp(1.0)-math::exp(-1.0));
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v6, v6_ex, eps );
#else
        FEEL_ASSERT( math::abs( v6-v6_ex) < eps )( v6 )(v6_ex)( math::abs( v6-v6_ex) )( eps ).warn ( "v6 != v6_ex" );
#endif /* USE_BOOST_TEST */


    }
    double meshSize;
};

template<int Order, typename value_type = double>
struct test_integration_functions
{
    test_integration_functions( double meshSize_=DEFAULT_MESH_SIZE ): meshSize(meshSize_) {}
    void operator()()
    {
        using namespace Feel;
        using namespace Feel::vf;


        typedef typename imesh<value_type>::type mesh_type;
        typename imesh<value_type>::ptrtype mesh( createMesh<value_type>( meshSize ) );

        const value_type eps = 1e-9;

        typedef fusion::vector<Lagrange<Order, Scalar> > basis_type;
        typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
        boost::shared_ptr<space_type> Xh( new space_type(mesh) );
        typename space_type::element_type u( Xh );

        // int ([-1,1],[-1,x]) 1 dx
        u = vf::project( Xh, elements(mesh), constant(1.0) );
        value_type v0 = integrate( elements(mesh), idv( u ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v0, 4.0, eps );
#else
        FEEL_ASSERT( math::abs( v0-4.0) < eps )( v0 )( math::abs( v0-4.0) )( eps ).warn ( "v0 != 4" );
#endif /* USE_BOOST_TEST */

        //
        u = vf::project( Xh, elements(mesh), Px() );
        value_type v1 = integrate( elements(mesh), idv( u ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v1- 0.0, eps );
#else
        FEEL_ASSERT( math::abs( v1-0.0) < eps )( v1 )( math::abs( v1-0.0) )( eps ).warn ( "v1 != 0" );
#endif /* USE_BOOST_TEST */
        value_type v2 = integrate( elements(mesh), dxv( u ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v2, 4.0, eps );
#else
        FEEL_ASSERT( math::abs( v2-4.0) < eps )( v2 )( math::abs( v2-4.0) )( eps ).warn ( "v2 != 4" );
#endif /* USE_BOOST_TEST */

        value_type v3 = integrate( elements(mesh), dyv( u ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v3- 0.0, eps );
#else
        FEEL_ASSERT( math::abs( v3-0.0) < eps )( v3 )( math::abs( v3-0.0) )( eps ).warn ( "v3 != 0" );
#endif /* USE_BOOST_TEST */

        u = vf::project( Xh, elements(mesh), exp(Px())*exp(Py()) );
        value_type v4 = integrate( elements(mesh), idv( u ) ).evaluate()( 0, 0 );
        value_type v4_ex = (math::exp(1.0)-math::exp(-1.0))*(math::exp(1.0)-math::exp(-1.0));
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v4, v4_ex, std::pow(10.0,-2.0*Order) );
#else
        FEEL_ASSERT( math::abs( v4-v4_ex) < std::pow(10.0,-2.0*Order) )
            ( v4 )
            ( v4_ex )
            ( math::abs( v4-v4_ex ) )( std::pow(10.0,-2.0*Order) ).warn ( "v4 != v4_ex" );
#endif /* USE_BOOST_TEST */
        value_type v4_x = integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), dxv( u ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v4_x, v4_ex, std::pow(10.0,-2.0*Order) );
#else
        FEEL_ASSERT( math::abs( v4_x-v4_ex) < std::pow(10.0,-2.0*Order) )
            ( v4_x )
            ( v4_ex )
            ( math::abs( v4_x-v4_ex ) )( std::pow(10.0,-2.0*Order) ).warn ( "v4_x != v4_ex" );
#endif /* USE_BOOST_TEST */

        value_type v4_y = integrate( elements(mesh), dyv( u ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v4_y, v4_ex, std::pow(10.0,-2.0*Order) );
#else
        FEEL_ASSERT( math::abs( v4_y-v4_ex) < std::pow(10.0,-2.0*Order) )
            ( v4_y )
            ( v4_ex )
            ( math::abs( v4_y-v4_ex ) )( std::pow(10.0,-2.0*Order) ).warn ( "v4_y != v4_ex" );
#endif /* USE_BOOST_TEST */

        value_type v5 = integrate( markedfaces(*mesh,10), idv(u) ).evaluate()( 0, 0 );
        value_type v5_ex = (math::exp(1.0)-math::exp(-1.0))*(math::exp(1.0)+math::exp(-1.0));
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v5, v5_ex, std::pow(10.0,-2.0*Order) );
#else
        FEEL_ASSERT( math::abs( v5-v5_ex) < std::pow(10.0,-2.0*Order) )
            ( v5 )
            ( v5_ex )
            ( math::abs( v5-v5_ex ) )( std::pow(10.0,-2.0*Order) ).warn ( "v5 != v5_ex" );
#endif /* USE_BOOST_TEST */

#if 1
        double meas = integrate( elements(mesh), cst(1.) ).evaluate()( 0, 0 );
        u = vf::project( Xh, elements(mesh), Px()*Px() + Py()*Py() +Px()*Py()  );
        auto hess6 = integrate( elements(mesh), hessv(u) ).evaluate();
#if defined(USE_BOOST_TEST)
        BOOST_TEST_MESSAGE( "hess6 =" << hess6 << "\n" );
        BOOST_CHECK_CLOSE( hess6(0,0), 2*meas, eps );
        BOOST_CHECK_CLOSE( hess6(0,1), meas, eps );
        BOOST_CHECK_CLOSE( hess6(1,0), meas, eps );
        BOOST_CHECK_CLOSE( hess6(1,0), hess6(0, 1), eps );
        BOOST_CHECK_CLOSE( hess6(1,1), 2*meas, eps );
#else
        FEEL_ASSERT( math::abs( v6-v6_ex) < std::pow(10.0,-2.0*Order) )
            ( v6 )
            ( v6_ex )
            ( math::abs( v6-v6_ex ) )( std::pow(10.0,-2.0*Order) ).warn ( "v6 != v6_ex" );
#endif /* USE_BOOST_TEST */
#endif


#if 0
        BOOST_TEST_MESSAGE( "idv(u) = " << integrate( markedfaces(*mesh,10), IM<2,Order+1,value_type,Simplex>(), idv(u) ).evaluate() << "\n" );
        BOOST_TEST_MESSAGE( "gradv(u) = " << integrate( markedfaces(*mesh,10), IM<2,Order+1,value_type,Simplex>(), gradv(u) ).evaluate() << "\n" );
        BOOST_TEST_MESSAGE( "trans(gradv(u))*N() = " << integrate( markedfaces(*mesh,10), IM<2,Order+1,value_type,Simplex>(), trans(gradv(u))*N() ).evaluate() << "\n" );
        BOOST_TEST_MESSAGE( "u*Nx()+u*Ny() = " << integrate( markedfaces(*mesh,10), IM<2,Order+1,value_type,Simplex>(), idv(u)*Nx() + idv(u)*Ny() ).evaluate() << "\n" );
        value_type v6 = integrate( markedfaces(*mesh,10), IM<2,Order+1,value_type,Simplex>(), trans(gradv(u))*N() ).evaluate()( 0, 0 );
        value_type v6_ex = (math::exp(1.0)-math::exp(-1.0))*(math::exp(1.0)+math::exp(-1.0));
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v6-v6_ex, std::pow(10.0,-2.0*Order) );
#else
        FEEL_ASSERT( math::abs( v6-v6_ex) < std::pow(10.0,-2.0*Order) )
            ( v6 )
            ( v6_ex )
            ( math::abs( v6-v6_ex ) )( std::pow(10.0,-2.0*Order) ).warn ( "v6 != v6_ex" );
#endif /* USE_BOOST_TEST */
#endif
    }
    double meshSize;
};

template<int Order, typename value_type = double>
struct test_integration_vectorial_functions
{
    test_integration_vectorial_functions( double meshSize_=DEFAULT_MESH_SIZE ): meshSize(meshSize_) {}
    void operator()()
    {
        using namespace Feel;
        using namespace Feel::vf;


        typedef typename imesh<value_type>::type mesh_type;
        typename imesh<value_type>::ptrtype mesh( createMesh<value_type>( meshSize ) );

        const value_type eps = 1000*Feel::type_traits<value_type>::epsilon();

        typedef fusion::vector<Lagrange<Order, Vectorial> > basis_type;
        typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
        boost::shared_ptr<space_type> Xh( new space_type(mesh) );
        typename space_type::element_type u( Xh );

        u = vf::project( Xh, elements(mesh), P() );
        BOOST_TEST_MESSAGE( "int(proj P() = " << integrate( elements(mesh), idv( u ) ).evaluate() << "\n" );
        u = vf::project( Xh, elements(mesh), one() );
        BOOST_TEST_MESSAGE( "int(proj one() = " << integrate( elements(mesh), idv( u ) ).evaluate() << "\n" );

        u = vf::project( Xh, elements(mesh), P() );
        ublas::matrix<value_type> m(  integrate( elements(mesh), gradv( u ) ).evaluate() );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( m(0,0), 4, std::pow(10.0,-2.0*Order) );
        BOOST_CHECK_SMALL( m(0,1)- 0, std::pow(10.0,-2.0*Order) );
        BOOST_CHECK_SMALL( m(1,0)- 0, std::pow(10.0,-2.0*Order) );
        BOOST_CHECK_CLOSE( m(1,1), 4, std::pow(10.0,-2.0*Order) );
#endif
        ublas::matrix<value_type> mx( integrate( elements(mesh), dxv( u ) ).evaluate());
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( mx(0,0), 4, std::pow(10.0,-2.0*Order) );
        BOOST_CHECK_SMALL( mx(1,0)- 0, std::pow(10.0,-2.0*Order) );
#endif

        ublas::matrix<value_type> my( integrate( elements(mesh), dyv( u ) ).evaluate());
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( my(0,0)- 0, std::pow(10.0,-2.0*Order) );
        BOOST_CHECK_CLOSE( my(1,0), 4, std::pow(10.0,-2.0*Order) );
#endif

        u = vf::project( Xh, elements(mesh), Py()*oneX() + Px()*oneY() );
        ublas::matrix<value_type> int_divu( integrate( elements(mesh), divv( u ) ).evaluate());
#if defined(USE_BOOST_TEST)
        value_type norm_int_divu = ublas::norm_frobenius(int_divu);
        BOOST_CHECK_SMALL( norm_int_divu, eps );
#endif
        u = vf::project( Xh, elements(mesh), P() );
        int_divu = integrate( elements(mesh), divv( u ) ).evaluate();
#if defined(USE_BOOST_TEST)
        norm_int_divu = ublas::norm_frobenius(int_divu);
        BOOST_CHECK_CLOSE( norm_int_divu, 8, eps );
#endif

        // check the divergence theorem
        u = vf::project( Xh, elements(mesh), P() );
        int_divu = integrate( elements(mesh), divv( u ) ).evaluate();
        BOOST_TEST_MESSAGE( "int_divu = " << int_divu << "\n" );
        ublas::matrix<value_type> int_un = integrate( boundaryfaces(mesh), trans(idv( u ))*N() ).evaluate();
        BOOST_TEST_MESSAGE( "(1, N) = " << integrate( boundaryfaces(mesh), IM<2,Order+1,value_type,Simplex>(), trans(one())*N() ).evaluate() << "\n" );
        BOOST_TEST_MESSAGE( "(P, N) = " << integrate( boundaryfaces(mesh), IM<2,Order+1,value_type,Simplex>(), trans(P())*N() ).evaluate() << "\n" );
        BOOST_TEST_MESSAGE( "int_un = " << int_un << "\n" );
#if defined(USE_BOOST_TEST)
        value_type norm_divergence = ublas::norm_frobenius(int_divu-int_un);
        BOOST_CHECK_SMALL( norm_divergence, eps );
#endif

        // check the stokes theorem
        ublas::matrix<value_type> int_curlu = integrate( elements(mesh), curlzv( u ) ).evaluate();
        BOOST_TEST_MESSAGE( "int_curlu = " << int_curlu << "\n" );
        ublas::matrix<value_type> int_ut = integrate( boundaryfaces(mesh), trans( idv( u ) )*(Ny()*oneX()-Nx()*oneY() ) ).evaluate();
        BOOST_TEST_MESSAGE( "int_ut = " << int_ut << "\n" );
#if defined(USE_BOOST_TEST)
        value_type norm_stokes = ublas::norm_frobenius(int_curlu-int_ut);
        BOOST_CHECK_SMALL( norm_stokes, eps );
#endif
    }
    double meshSize;
};


template<int Order, typename value_type = double>
struct test_integration_composite_functions
{
    test_integration_composite_functions( double meshSize_=DEFAULT_MESH_SIZE ): meshSize(meshSize_) {}
    void operator()()
    {
        using namespace Feel;
        using namespace Feel::vf;


        typedef typename imesh<value_type>::type mesh_type;
        typename imesh<value_type>::ptrtype mesh( createMesh<value_type>( meshSize ) );

        const value_type eps = 1000*Feel::type_traits<value_type>::epsilon();

        typedef fusion::vector<Lagrange<Order, Vectorial>,
            Lagrange<Order-1, Scalar> > basis_type;
        typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
        boost::shared_ptr<space_type> Xh( new space_type(mesh) );
        typename space_type::element_type u( Xh );

        u.template element<0>() = vf::project( Xh->template functionSpace<0>(), elements(mesh), P() );
        u.template element<1>() = vf::project( Xh->template functionSpace<1>(), elements(mesh), constant(1) );
        BOOST_TEST_MESSAGE( "int(proj P() = " << integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), idv( u.template element<0>() ) ).evaluate() << "\n" );
        BOOST_TEST_MESSAGE( "int(1 = " << integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), idv( u.template element<1>() ) ).evaluate() << "\n" );

        ublas::matrix<value_type> m(  integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), gradv( u.template element<0>() ) ).evaluate() );
        BOOST_TEST_MESSAGE( "int(grad(P()) = " << m << "\n" );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( m(0,0)-4, std::pow(10.0,-2.0*Order) );
        BOOST_CHECK_SMALL( m(0,1), std::pow(10.0,-2.0*Order) );
        BOOST_CHECK_SMALL( m(1,0), std::pow(10.0,-2.0*Order) );
        BOOST_CHECK_SMALL( m(1,1)-4, std::pow(10.0,-2.0*Order) );
#endif
        m = integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), trace( gradv( u.template element<0>() )*trans(gradv( u.template element<0>() ) ) ) ).evaluate();
        BOOST_TEST_MESSAGE( "int(grad(P()*grad^T(P())) = " << m << "\n" );

        m = integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), trace( trans(gradv( u.template element<0>() ))*(gradv( u.template element<0>() ) ) ) ).evaluate();
        BOOST_TEST_MESSAGE( "int(grad(P()^T*grad(P())) = " << m << "\n" );

#if 0
        AUTO( u_exact,(P()^(2))*(Px()+Py()));
        AUTO( grad_exact, (mat<2,2>( 3*Px()*Px()+2*Px()*Py(), (Px()^(2)), (Py()^(2)), 3*Py()*Py()+2*Py()*Px() ) ) );
        AUTO( div_grad_exact, vec( 6*Px()+2*Py(), 6*Py()+2*Px() ) );
#else
        //AUTO( u_exact, vec(sin(Px())*sin(Py()), cos(Px())*cos(Py()) ) );
        AUTO( u_exact, sin(Px())*sin(Py())*oneX() + cos(Px())*cos(Py())*oneY() );
        AUTO( p_exact, val(cos(Px())*sin(Py())));
        AUTO( du_dx, val(cos(Px())*sin(Py())));
        AUTO( du_dy, val(sin(Px())*cos(Py())));

        AUTO( dv_dx, val(-sin(Px())*cos(Py())));
        AUTO( dv_dy, val(-cos(Px())*sin(Py())));

        AUTO( grad_exact, (mat<2,2>(du_dx, du_dy, dv_dx, dv_dy)) );
        AUTO( div_grad_exact, (vec(-sin(Px())*sin(Py())-sin(Px())*sin(Py()),
                                   -cos(Px())*cos(Py())-cos(Px())*cos(Py()) ) ) );

#endif

        u.template element<0>() = vf::project( Xh->template functionSpace<0>(), elements(mesh), u_exact );
        u.template element<1>() = vf::project( Xh->template functionSpace<1>(), elements(mesh), Px()+Py() );
        BOOST_TEST_MESSAGE( "int(u) = " << integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), idv( u.template element<0>() ) ).evaluate() << "\n" );
        BOOST_TEST_MESSAGE( "int(u - u_exact) = " << integrate( elements(mesh), IM<2,Order,value_type,Simplex>(),
                                                     trans(idv( u.template element<0>() )-u_exact)*(idv( u.template element<0>() )-u_exact) ).evaluate() << "\n" );

        BOOST_TEST_MESSAGE( "int(u<1>) = " << integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), idv( u.template element<1>() ) ).evaluate() << "\n" );

        m= integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), gradv( u.template element<0>() ) ).evaluate();
        BOOST_TEST_MESSAGE( "int(grad(u)) = " << m << "\n" );

        m= integrate( boundaryfaces(mesh), IM<2,Order,value_type,Simplex>(), gradv( u.template element<0>() )*N() ).evaluate();
        BOOST_TEST_MESSAGE( "int_bfaces(grad_u*N()) = " << m << "\n" );

        m= integrate( boundaryfaces(mesh), IM<2,Order,value_type,Simplex>(), grad_exact*N() ).evaluate();
        BOOST_TEST_MESSAGE( "int_bfaces(grad_exact*N) = " << m << "\n" );

        m= integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), div_grad_exact ).evaluate();
        BOOST_TEST_MESSAGE( "int(div_grad_exact) = " << m << "\n" );

        m= integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), grad_exact  ).evaluate();
        BOOST_TEST_MESSAGE( "int(grad_exact) = " << m << "\n" );
        m= integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), gradv( u.template element<0>() ) - grad_exact ).evaluate();
        BOOST_TEST_MESSAGE( "int( grad(u)-grad_exact) = " << m << "\n" );


        m= integrate( elements(mesh), IM<2,Order,value_type,Simplex>(),
                      trace( (gradv( u.template element<0>())-grad_exact)*trans(gradv( u.template element<0>())-grad_exact) )
                      ).evaluate();
        BOOST_TEST_MESSAGE( "|grad(u)-grad_exact|_0^2 = " << m << "\n" );




        typename  BackendGmm<value_type>::vector_ptrtype F = BackendGmm<value_type>::newVector( Xh );
        form( Xh, *F ) = integrate( elements( *mesh ),
                                    IM<2,Order,value_type,Simplex>(),
                                    ( trans(vec(constant(1.0), constant(1.0) )) * id(u.template element<0>())+
                                      id(u.template element<1>()) ) );
        F->printMatlab( "composite_F.m" );
        F->close();
        value_type vf = inner_product( *F, u );
        value_type vv1 = integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), idv( u.template element<0>() ) ).evaluate()( 0, 0 );
        value_type vv2 = integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), idv( u.template element<0>() ) ).evaluate()( 1, 0 );
        value_type vv3 = integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), idv( u.template element<1>() ) ).evaluate()( 0, 0 );


#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( vf-(vv1+vv2+vv3), eps );
#endif
    }
    double meshSize;
};


inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description integrationoptions("Test Integration options");
    integrationoptions.add_options()
        ("hsize", Feel::po::value<double>()->default_value( 0.3 ), "h value")
        ;
    return integrationoptions.add( Feel::feel_options() );
}

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "test_integration" ,
                           "test_integration" ,
                            "0.1",
                           "integration tests",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2006,2007 Université Joseph Fourier (Grenoble I)");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}

#if defined(USE_BOOST_TEST)

BOOST_AUTO_TEST_CASE( test_integration_1 ) { test_integration_circle<double> t( 0.02 ); t(); }
BOOST_AUTO_TEST_CASE( test_integration_2 ) { test_integration_domain<double> t( 0.1 ); t(); }
BOOST_AUTO_TEST_CASE( test_integration_3 ){ test_integration_boundary<double> t( 0.1 ); t(); }
BOOST_AUTO_TEST_CASE( test_integration_4 ) { test_integration_functions<2,double> t( 0.1 ); t();}
BOOST_AUTO_TEST_CASE( test_integration_5 ) { test_integration_vectorial_functions<2,double> t( 0.1 ); t(); }
BOOST_AUTO_TEST_CASE( test_integration_6 ) { test_integration_composite_functions<2,double> t(0.1); t(); }
BOOST_AUTO_TEST_CASE( test_integration_7 ) { test_integration_simplex<double> t( 0.1 ); t(); }

int BOOST_TEST_CALL_DECL
main( int argc, char* argv[] )
{
    Feel::Environment env( argc, argv );
    int ret = ::boost::unit_test::unit_test_main( &init_unit_test, argc, argv );

    return ret;
}

#else
int
main( int argc, char** argv )
{
    Feel::Application mpi( argc, argv, makeAbout(), makeOptions() );
    Feel::Assert::setLog( "test_integration.assert");

    test_integration_circle<double> t1( mpi.vm()["hsize"].as<double>() ); t1();
    test_integration_domain<double> t2( mpi.vm()["hsize"].as<double>() ); t2();
    test_integration_boundary<double> t3( mpi.vm()["hsize"].as<double>() ); t3();
    test_integration_functions<2,double> t4( mpi.vm()["hsize"].as<double>() ); t4();
    test_integration_vectorial_functions<2,double> t5(mpi.vm()["hsize"].as<double>() ); t5();
    test_integration_composite_functions<2,double> t6(mpi.vm()["hsize"].as<double>() ); t6();
}
#endif // USE_BOOST_TEST
