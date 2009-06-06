/* -*- mode: c++ -*-

  This file is part of the Life library

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
    typedef Life::uint16_type uint16_type;
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
    typedef Life::uint16_type uint16_type;
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
    typedef Life::uint16_type uint16_type;
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
    typedef Life::uint16_type uint16_type;
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
    //std::cout << "hsize = " << meshSize << std::endl;

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
    //std::cout << "hsize = " << meshSize << std::endl;

    GmshSimplexDomain<2,Order> ts;
    ts.setCharacteristicLength( meshSize );
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
        using namespace Life;
        using namespace Life::vf;

        int Order=2;
        double t = 0.0;
        AUTO( mycst, cst_ref( t ) );
        typename imesh<value_type,1>::ptrtype mesh( createCircle<value_type,1>( meshSize ) );

        t = 1.0;
        value_type v0 = integrate( elements(mesh), IM<2,2,value_type,Simplex>(), mycst ).evaluate()( 0, 0 );
        std::cout << "v0=" << v0 << "\n";
        value_type v00 = ( integrate( boundaryelements(mesh), IM<2,2,value_type,Simplex>(), mycst ).evaluate()( 0, 0 )+
                           integrate( internalelements(mesh), IM<2,2,value_type,Simplex>(), mycst ).evaluate()( 0, 0 ) );
        std::cout << "v00=" << v00 << "\n";
        std::cout << "[circle] v0 0 = " << integrate( boundaryfaces(mesh), IM<2,2,value_type,Simplex>(), N() ).evaluate() << "\n";
        std::cout << "[circle] v0 0 = " << integrate( boundaryfaces(mesh), IM<2,2,value_type,Simplex>(), N() ).evaluate() << "\n";
        std::cout << "[circle] v0 1 = " << integrate( boundaryfaces(mesh), IM<2,2,value_type,Simplex>(), vec( idf(f_Nx()), idf(f_Ny() ) ) ).evaluate() << "\n";
        std::cout << "[circle] v0 2 = " << integrate( boundaryfaces(mesh), IM<2,2,value_type,Simplex>(),
                                                        trans(vec(constant(1),Ny()))*one() ).evaluate() << "\n";

        std::cout << "[circle] v0 3 = " << integrate( boundaryfaces(mesh), IM<2,2,value_type,Simplex>(),
                                                      mat<2,2>(Nx(),constant(1),constant(1),Ny())*vec(constant(1),constant(1)) ).evaluate() << "\n";
        std::cout << "[circle] v0 4 (==v0 3) = " << integrate( boundaryfaces(mesh), IM<2,2,value_type,Simplex>(),
                                                               mat<2,2>( idf(f_Nx()),constant(1),
                                                                         constant(1),idf(f_Ny()) )*vec(constant(1),constant(1)) ).evaluate() << "\n";
        value_type pi = 4.0*math::atan(1.0);
        const value_type eps = 1000*Life::type_traits<value_type>::epsilon();
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v0-pi, math::pow( meshSize, 2*Order ) );
        BOOST_CHECK_SMALL( v0-v00, eps  );
#else
        LIFE_ASSERT( math::abs( v0-pi) < math::pow( meshSize, 2*Order ) )( v0 )( math::abs( v0-pi) )( math::pow( meshSize, 2*Order ) ).warn ( "v0 != pi" );
        LIFE_ASSERT( math::abs( v0-v00) < eps )( v0 )( v00 )( math::abs( v0-v00) )( eps ).warn ( "v0 != pi" );
#endif /* USE_BOOST_TEST */

        typedef typename imesh<value_type,1>::type mesh_type;
        typedef fusion::vector<Lagrange<2, Scalar> > basis_type;
        typedef FunctionSpace<mesh_type, basis_type> space_type;
        boost::shared_ptr<space_type> Xh( new space_type(mesh) );
        typename space_type::element_type u( Xh );

        // int ([-1,1],[-1,x]) 1 dx
        u = project( Xh, elements(mesh), constant(1.0) );
        v0 = integrate( elements(mesh), IM<2,3,value_type,Simplex>(), idv( u ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v0-pi, math::pow( meshSize, 2*Order ) );
#else
        LIFE_ASSERT( math::abs( v0-pi) < math::pow( meshSize, 2*Order ) )( v0 )( math::abs( v0-pi) )( math::pow( meshSize, 2*Order ) ).warn ( "v0 != pi" );
#endif /* USE_BOOST_TEST */

        typedef fusion::vector<Lagrange<2, Vectorial> > vector_basis_type;
        typedef FunctionSpace<mesh_type, vector_basis_type, value_type> vector_space_type;
        boost::shared_ptr<vector_space_type> Xvh( new vector_space_type(mesh) );
        typename vector_space_type::element_type U( Xvh );

        U = project( Xvh, elements(mesh), vec(constant(1.0),constant(1.0)) );
        v0 = integrate( boundaryfaces(mesh), IM<2,3,value_type,Simplex>(), trans(idv( U ))*N() ).evaluate()( 0, 0 );
        v00 = integrate( elements(mesh), IM<2,3,value_type,Simplex>(), divv(U) ).evaluate()( 0, 0 );

        std::cout << "[circle] v0=" << v0 << " v00=" << v00 << " should be equal thanks to gauss int 1.N() != int div 1\n" ;
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v0-v00, eps );
#else
        LIFE_ASSERT( math::abs( v0-v00) < eps )( v0 )(v00)( math::abs( v0-v00) ).warn ( "int 1.N() != int div 1" );
#endif /* USE_BOOST_TEST */

        U = project( Xvh, elements(mesh), vec(Px(),Py()) );
        v0 = integrate( boundaryfaces(mesh), IM<2,3,value_type,Simplex>(), trans(idv( U ))*N() ).evaluate()( 0, 0 );
        v00 = integrate( elements(mesh), IM<2,3,value_type,Simplex>(), divv(U) ).evaluate()( 0, 0 );

        std::cout << "[circle] v0=" << v0 << " v00=" << v00 << " should be equal thanks to gauss int (x,y).N() != int div (x,y)\n" ;
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v0-v00, eps );
#else
        LIFE_ASSERT( math::abs( v0-v00) < eps )( v0 )(v00)( math::abs( v0-v00) ).warn ( "int (x,y).N() != int div (x,y)" );
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
        using namespace Life;
        using namespace Life::vf;



        typename imesh<value_type,1>::ptrtype mesh( createSimplex<value_type,1>( meshSize ) );
        typedef typename imesh<value_type>::type mesh_type;

        typedef fusion::vector<Lagrange<3, Scalar> > basis_type;
        typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
        boost::shared_ptr<space_type> Xh( new space_type(mesh) );
        typename space_type::element_type u( Xh );

        value_type v0 = integrate( elements(mesh), IM<2,2,value_type,Simplex>(), constant(1.0) ).evaluate()( 0, 0 );
        std::cout << "[simplex] v0 0 = " << integrate( boundaryfaces(mesh), IM<2,2,value_type,Simplex>(), N() ).evaluate() << "\n";
        std::cout << "[simplex] v0 1 = " << integrate( boundaryfaces(mesh), IM<2,2,value_type,Simplex>(), vec( idf(f_Nx()), idf(f_Ny() ) ) ).evaluate() << "\n";
        std::cout << "[simplex] v0 2 = " << integrate( boundaryfaces(mesh), IM<2,2,value_type,Simplex>(),
                                                        trans(vec(constant(1),Ny()))*one() ).evaluate() << "\n";

        std::cout << "[simplex] v0 3 = " << integrate( boundaryfaces(mesh), IM<2,2,value_type,Simplex>(),
                                                      mat<2,2>(Nx(),constant(1),constant(1),Ny())*vec(constant(1),constant(1)) ).evaluate() << "\n";
        std::cout << "[simplex] v0 4 (==v0 3) = " << integrate( boundaryfaces(mesh), IM<2,2,value_type,Simplex>(),
                                                               mat<2,2>( idf(f_Nx()),constant(1),
                                                                         constant(1),idf(f_Ny()) )*vec(constant(1),constant(1)) ).evaluate() << "\n";
        value_type pi = 4.0*math::atan(1.0);
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v0-pi, 1e-2 );
#else
        const value_type eps = 1000*Life::type_traits<value_type>::epsilon();

        LIFE_ASSERT( math::abs( v0-pi) < 1e-2 )( v0 )( math::abs( v0-pi) )( eps ).warn ( "v0 != pi" );
#endif /* USE_BOOST_TEST */

        // int ([-1,1],[-1,x]) 1 dx
        u = project( Xh, elements(mesh), Px() );
        std::cout << "[simplex] int grad(X)=" << integrate( elements(mesh), IM<2,1,value_type,Simplex>(), gradv(u) ).evaluate() << "\n";
        std::cout << "[simplex] int hess(X)=" << integrate( elements(mesh), IM<2,1,value_type,Simplex>(), hessv(u) ).evaluate() << "\n";
        u = project( Xh, elements(mesh), Px()*Px()+Px()*Py()+Py()*Py() );
        std::cout << "[simplex] int grad(X^2)=" << integrate( elements(mesh), IM<2,1,value_type,Simplex>(), gradv(u) ).evaluate() << "\n";
        std::cout << "[simplex] int hess(X^2)=" << integrate( elements(mesh), IM<2,1,value_type,Simplex>(), hessv(u) ).evaluate() << "\n";
        u = project( Xh, elements(mesh), Px()*Px()*Px()+Py()*Py()*Py() );
        std::cout << "[simplex] int grad(X^3)=" << integrate( elements(mesh), IM<2,1,value_type,Simplex>(), gradv(u) ).evaluate() << "\n";
        std::cout << "[simplex] int hess(X^3)=" << integrate( elements(mesh), IM<2,1,value_type,Simplex>(), hessv(u) ).evaluate() << "\n";
    }
    double meshSize;
};
template<typename value_type = double>
struct test_integration_domain
{
    test_integration_domain( double meshSize_=DEFAULT_MESH_SIZE ): meshSize(meshSize_) {}
    void operator()()
    {
        using namespace Life;
        using namespace Life::vf;


        typename imesh<value_type>::ptrtype mesh( createMesh<value_type>( meshSize ) );

        const value_type eps = 1000*Life::type_traits<value_type>::epsilon();

        // int ([-1,1],[-1,x]) 1 dx
        value_type v0 = integrate( elements(mesh), IM<2,1,value_type,Simplex>(),
                                   vf::min(constant(1.0),constant(2.0)) ).evaluate()( 0, 0 );
        value_type v00 = ( integrate( boundaryelements(mesh), IM<2,1,value_type,Simplex>(),
                                      vf::min(constant(1.0),constant(2.0)) ).evaluate()( 0, 0 )+
                           integrate( internalelements(mesh), IM<2,1,value_type,Simplex>(),
                                      vf::min(constant(1.0),constant(2.0)) ).evaluate()( 0, 0 ) );

        std::cout << "[domain] v0 = " << v0 << "\n";
        std::cout << "[domain] v00 = " << v00 << "\n";
        std::cout << "[domain] int(1*N()) = " << integrate( boundaryfaces(mesh), IM<2,1,value_type,Simplex>(),
                                                   constant(1.0)*N() ).evaluate() << "\n";
        std::cout << "[domain] int(tr(N())**N()) = " << integrate( boundaryfaces(mesh), IM<2,1,value_type,Simplex>(),
                                                          N()*trans(N()) ).evaluate() << "\n";
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v0-4.0, eps );
        BOOST_CHECK_SMALL( v0-v00, eps );
#else
        LIFE_ASSERT( math::abs( v0-4.0) < eps )( v0 )( math::abs( v0-4.0) )( eps ).warn ( "v0 != 4" );
        LIFE_ASSERT( math::abs( v0-v00) < eps )( v0 )(v00)( math::abs( v0-v00) )( eps ).warn ( "v0 != v00" );
#endif /* USE_BOOST_TEST */


        // int ([-1,1],[-1,1]) x dx
        value_type v1 = integrate( elements(mesh), IM<2,2,value_type,Simplex>(), Px() ).evaluate()( 0, 0 );
        value_type v11 = integrate( elements(mesh), IM<2,2,value_type,Simplex>(), idf(f_Px()) ).evaluate()( 0, 0 );

        std::cout << "[domain] int(P()) = " << integrate( elements(mesh), IM<2,1,value_type,Simplex>(),
                                                 P() ).evaluate() << "\n";
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v1-0.0, eps );
        BOOST_CHECK_SMALL( v11-0.0, eps );
#else
        LIFE_ASSERT( math::abs( v1-0.0) < eps )( v1 )( math::abs( v1-0.0) )( eps ).warn ( "v1 != 0" );
        LIFE_ASSERT( math::abs( v11-0.0) < eps )( v11 )( math::abs( v11-0.0) )( eps ).warn ( "v11 != 0" );
#endif /* USE_BOOST_TEST */

         // int ([-1,1],[-1,1]) abs(x) dx
        value_type vsin = integrate( elements(mesh), IM<2,5,value_type,Simplex>(), sin(Px()) ).evaluate()( 0, 0 );
        value_type vsin1 = integrate( elements(mesh), IM<2,5,value_type,Simplex>(), idf(f_sinPx()) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( vsin-2*(-math::cos(1.0)+math::cos(-1.0)), eps );
        BOOST_CHECK_SMALL( vsin1-2*(-math::cos(1.0)+math::cos(-1.0)), eps );
#else
        LIFE_ASSERT( math::abs( vsin-2.0*(-math::cos(1.0)+math::cos(-1.0))) < eps )
            ( vsin )
            ( math::abs( vsin-2.0*(-math::cos(1.0)+math::cos(-1.0))) )( eps ).warn ( "vsin != 2*(cos(1)-cos(-1))" );
        LIFE_ASSERT( math::abs( vsin1-2.0*(-math::cos(1.0)+math::cos(-1.0))) < eps )
            ( vsin1 )
            ( math::abs( vsin1-2.0*(-math::cos(1.0)+math::cos(-1.0))) )( eps ).warn ( "vsin1 != 2*(cos(1)-cos(-1))" );
#endif /* USE_BOOST_TEST */

        // int ([-1,1],[-1,1]) abs(x) dx
        value_type vabs = integrate( elements(mesh), IM<2,5,value_type,Simplex>(), abs(Px())+abs(Py()) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( vabs-4.0, 1e-2 );
#else
        LIFE_ASSERT( math::abs( vabs-4.0) < 1e-2 )( vabs )( math::abs( vabs-4.0) )( 1e-2 ).warn ( "vabs != 4" );
#endif /* USE_BOOST_TEST */



        // int (\partial ([-1,1],[-1,1]) 1 dx
        value_type v2 = integrate( boundaryfaces(mesh), IM<2,1,value_type,Simplex>(), constant(1.0) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v2-8.0, eps );
#else
        LIFE_ASSERT( math::abs( v2-8.0) < eps )( v2 )( math::abs( v2-8.0) )( eps ).warn ( "v2 != 8" );
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
        using namespace Life;
        using namespace Life::vf;


        typename imesh<value_type>::ptrtype mesh( createMesh<value_type>( meshSize ) );

        const value_type eps = 1000*Life::type_traits<value_type>::epsilon();

        // int (\partial ([-1,1],[-1,1]) 1 dx
        value_type v2 = integrate( boundaryfaces(mesh), IM<2,1,value_type,Simplex>(), constant(1.0) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v2-8.0, eps );
#else
        LIFE_ASSERT( math::abs( v2-8.0) < eps )( v2 )( math::abs( v2-8.0) )( eps ).warn ( "v2 != 8" );
#endif /* USE_BOOST_TEST */

        // int (\partial ([-1,1],[-1,1]) x * y dx dy
        value_type v3 = integrate( boundaryfaces(mesh), IM<2,2,value_type,Simplex>(), Px()*Py() ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v3-0.0, eps );
#else
        LIFE_ASSERT( math::abs( v3-0.0) < eps )( v3 )( math::abs( v3-0.0) )( eps ).warn ( "v3 != 0" );
#endif /* USE_BOOST_TEST */

        // int (\partial ([-1,1],[-1,1]) -2*y dx dy
        value_type v4 = integrate( boundaryfaces(mesh), IM<2,2,value_type,Simplex>(), -2*Py() ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v4-0.0, eps );
#else
        LIFE_ASSERT( math::abs( v4-0.0) < eps )( v4 )( math::abs( v4-0.0) )( eps ).warn ( "v4 != 0" );
#endif /* USE_BOOST_TEST */

        // int_10 exp(Py) dx
        value_type v5 = integrate( markedfaces(*mesh,10), IM<2,2,value_type,Simplex>(), exp(Py()) ).evaluate()( 0, 0 );
        value_type v5_ex = 2.0*(math::exp(1.0)+math::exp(-1.0));
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v5-v5_ex, eps );
#else
        LIFE_ASSERT( math::abs( v5-v5_ex) < eps )( v5 )(v5_ex)( math::abs( v5-v5_ex) )( eps ).warn ( "v5 != v5_ex" );
#endif /* USE_BOOST_TEST */

        // int_20 exp(Py) dx
        value_type v6 = integrate( markedfaces(*mesh,20), IM<2,5,value_type,Simplex>(), exp(Py()) ).evaluate()( 0, 0 );
        value_type v6_ex = 2.0*(math::exp(1.0)-math::exp(-1.0));
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v6-v6_ex, eps );
#else
        LIFE_ASSERT( math::abs( v6-v6_ex) < eps )( v6 )(v6_ex)( math::abs( v6-v6_ex) )( eps ).warn ( "v6 != v6_ex" );
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
        using namespace Life;
        using namespace Life::vf;


        typedef typename imesh<value_type>::type mesh_type;
        typename imesh<value_type>::ptrtype mesh( createMesh<value_type>( meshSize ) );

        const value_type eps = 1000*Life::type_traits<value_type>::epsilon();

        typedef fusion::vector<Lagrange<Order, Scalar> > basis_type;
        typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
        boost::shared_ptr<space_type> Xh( new space_type(mesh) );
        typename space_type::element_type u( Xh );

        // int ([-1,1],[-1,x]) 1 dx
        u = project( Xh, elements(mesh), constant(1.0) );
        value_type v0 = integrate( elements(mesh), IM<2,1,value_type,Simplex>(), idv( u ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v0-4.0, eps );
#else
        LIFE_ASSERT( math::abs( v0-4.0) < eps )( v0 )( math::abs( v0-4.0) )( eps ).warn ( "v0 != 4" );
#endif /* USE_BOOST_TEST */

        //
        u = project( Xh, elements(mesh), Px() );
        value_type v1 = integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), idv( u ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v1-0.0, eps );
#else
        LIFE_ASSERT( math::abs( v1-0.0) < eps )( v1 )( math::abs( v1-0.0) )( eps ).warn ( "v1 != 0" );
#endif /* USE_BOOST_TEST */
        value_type v2 = integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), dxv( u ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v2-4.0, eps );
#else
        LIFE_ASSERT( math::abs( v2-4.0) < eps )( v2 )( math::abs( v2-4.0) )( eps ).warn ( "v2 != 4" );
#endif /* USE_BOOST_TEST */

        value_type v3 = integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), dyv( u ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v3-0.0, eps );
#else
        LIFE_ASSERT( math::abs( v3-0.0) < eps )( v3 )( math::abs( v3-0.0) )( eps ).warn ( "v3 != 0" );
#endif /* USE_BOOST_TEST */

        u = project( Xh, elements(mesh), exp(Px())*exp(Py()) );
        value_type v4 = integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), idv( u ) ).evaluate()( 0, 0 );
        value_type v4_ex = (math::exp(1.0)-math::exp(-1.0))*(math::exp(1.0)-math::exp(-1.0));
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v4-v4_ex, std::pow(10.0,-2.0*Order) );
#else
        LIFE_ASSERT( math::abs( v4-v4_ex) < std::pow(10.0,-2.0*Order) )
            ( v4 )
            ( v4_ex )
            ( math::abs( v4-v4_ex ) )( std::pow(10.0,-2.0*Order) ).warn ( "v4 != v4_ex" );
#endif /* USE_BOOST_TEST */
        value_type v4_x = integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), dxv( u ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v4_x-v4_ex, std::pow(10.0,-2.0*Order) );
#else
        LIFE_ASSERT( math::abs( v4_x-v4_ex) < std::pow(10.0,-2.0*Order) )
            ( v4_x )
            ( v4_ex )
            ( math::abs( v4_x-v4_ex ) )( std::pow(10.0,-2.0*Order) ).warn ( "v4_x != v4_ex" );
#endif /* USE_BOOST_TEST */

        value_type v4_y = integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), dyv( u ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v4_y-v4_ex, std::pow(10.0,-2.0*Order) );
#else
        LIFE_ASSERT( math::abs( v4_y-v4_ex) < std::pow(10.0,-2.0*Order) )
            ( v4_y )
            ( v4_ex )
            ( math::abs( v4_y-v4_ex ) )( std::pow(10.0,-2.0*Order) ).warn ( "v4_y != v4_ex" );
#endif /* USE_BOOST_TEST */

        value_type v5 = integrate( markedfaces(*mesh,10), IM<2,Order+1,value_type,Simplex>(), idv(u) ).evaluate()( 0, 0 );
        value_type v5_ex = (math::exp(1.0)-math::exp(-1.0))*(math::exp(1.0)+math::exp(-1.0));
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v5-v5_ex, std::pow(10.0,-2.0*Order) );
#else
        LIFE_ASSERT( math::abs( v5-v5_ex) < std::pow(10.0,-2.0*Order) )
            ( v5 )
            ( v5_ex )
            ( math::abs( v5-v5_ex ) )( std::pow(10.0,-2.0*Order) ).warn ( "v5 != v5_ex" );
#endif /* USE_BOOST_TEST */

#if 1
        u = project( Xh, elements(mesh), Px()*Px() );
        std::cout << "hess(X^2)="
                  << integrate( elements(mesh), IM<2,Order-2,value_type,Simplex>(), hessv(u) ).evaluate()
                  << "\n";
        value_type v6 = integrate( elements(mesh), IM<2,Order-2,value_type,Simplex>(),
                                   trace( hessv(u)*trans(hessv(u)) ) ).evaluate()( 0, 0 );
        std::cout << "int(hess(x^2):hess(x^2)) = " << v6 << "\n";
        value_type v6_ex = 16.0;
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v6-v6_ex, std::pow(10.0,-2.0*Order) );
#else
        LIFE_ASSERT( math::abs( v6-v6_ex) < std::pow(10.0,-2.0*Order) )
            ( v6 )
            ( v6_ex )
            ( math::abs( v6-v6_ex ) )( std::pow(10.0,-2.0*Order) ).warn ( "v6 != v6_ex" );
#endif /* USE_BOOST_TEST */
#endif


#if 0
        std::cout << "idv(u) = " << integrate( markedfaces(*mesh,10), IM<2,Order+1,value_type,Simplex>(), idv(u) ).evaluate() << "\n";
        std::cout << "gradv(u) = " << integrate( markedfaces(*mesh,10), IM<2,Order+1,value_type,Simplex>(), gradv(u) ).evaluate() << "\n";
        std::cout << "trans(gradv(u))*N() = " << integrate( markedfaces(*mesh,10), IM<2,Order+1,value_type,Simplex>(), trans(gradv(u))*N() ).evaluate() << "\n";
        std::cout << "u*Nx()+u*Ny() = " << integrate( markedfaces(*mesh,10), IM<2,Order+1,value_type,Simplex>(), idv(u)*Nx() + idv(u)*Ny() ).evaluate() << "\n";
        value_type v6 = integrate( markedfaces(*mesh,10), IM<2,Order+1,value_type,Simplex>(), trans(gradv(u))*N() ).evaluate()( 0, 0 );
        value_type v6_ex = (math::exp(1.0)-math::exp(-1.0))*(math::exp(1.0)+math::exp(-1.0));
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v6-v6_ex, std::pow(10.0,-2.0*Order) );
#else
        LIFE_ASSERT( math::abs( v6-v6_ex) < std::pow(10.0,-2.0*Order) )
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
        using namespace Life;
        using namespace Life::vf;


        typedef typename imesh<value_type>::type mesh_type;
        typename imesh<value_type>::ptrtype mesh( createMesh<value_type>( meshSize ) );

        const value_type eps = 1000*Life::type_traits<value_type>::epsilon();

        typedef fusion::vector<Lagrange<Order, Vectorial> > basis_type;
        typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
        boost::shared_ptr<space_type> Xh( new space_type(mesh) );
        typename space_type::element_type u( Xh );

        u = project( Xh, elements(mesh), P() );
        std::cout << "int(proj P() = " << integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), idv( u ) ).evaluate() << "\n";
        u = project( Xh, elements(mesh), one() );
        std::cout << "int(proj one() = " << integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), idv( u ) ).evaluate() << "\n";

        u = project( Xh, elements(mesh), P() );
        ublas::matrix<value_type> m(  integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), gradv( u ) ).evaluate() );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( m(0,0)-4, std::pow(10.0,-2.0*Order) );
        BOOST_CHECK_SMALL( m(0,1), std::pow(10.0,-2.0*Order) );
        BOOST_CHECK_SMALL( m(1,0), std::pow(10.0,-2.0*Order) );
        BOOST_CHECK_SMALL( m(1,1)-4, std::pow(10.0,-2.0*Order) );
#endif
        ublas::matrix<value_type> mx( integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), dxv( u ) ).evaluate());
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( mx(0,0)-4, std::pow(10.0,-2.0*Order) );
        BOOST_CHECK_SMALL( mx(1,0), std::pow(10.0,-2.0*Order) );
#endif

        ublas::matrix<value_type> my( integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), dyv( u ) ).evaluate());
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( my(0,0), std::pow(10.0,-2.0*Order) );
        BOOST_CHECK_SMALL( my(1,0)-4, std::pow(10.0,-2.0*Order) );
#endif

        u = project( Xh, elements(mesh), Py()*oneX() + Px()*oneY() );
        ublas::matrix<value_type> int_divu( integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), divv( u ) ).evaluate());
#if defined(USE_BOOST_TEST)
        value_type norm_int_divu = ublas::norm_frobenius(int_divu);
        BOOST_CHECK_SMALL( norm_int_divu, eps );
#endif
        u = project( Xh, elements(mesh), P() );
        int_divu = integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), divv( u ) ).evaluate();
#if defined(USE_BOOST_TEST)
        norm_int_divu = ublas::norm_frobenius(int_divu);
        BOOST_CHECK_SMALL( norm_int_divu-8, eps );
#endif

        // check the divergence theorem
        u = project( Xh, elements(mesh), P() );
        int_divu = integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), divv( u ) ).evaluate();
        std::cout << "int_divu = " << int_divu << "\n";
        ublas::matrix<value_type> int_un = integrate( boundaryfaces(mesh), IM<2,Order+3,value_type,Simplex>(), trans(idv( u ))*N() ).evaluate();
        std::cout << "(1, N) = " << integrate( boundaryfaces(mesh), IM<2,Order+1,value_type,Simplex>(), trans(one())*N() ).evaluate() << "\n";
        std::cout << "(P, N) = " << integrate( boundaryfaces(mesh), IM<2,Order+1,value_type,Simplex>(), trans(P())*N() ).evaluate() << "\n";
        std::cout << "int_un = " << int_un << "\n";
#if defined(USE_BOOST_TEST)
        value_type norm_divergence = ublas::norm_frobenius(int_divu-int_un);
        BOOST_CHECK_SMALL( norm_divergence, eps );
#endif

        // check the stokes theorem
        ublas::matrix<value_type> int_curlu = integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), curlzv( u ) ).evaluate();
        std::cout << "int_curlu = " << int_curlu << "\n";
        ublas::matrix<value_type> int_ut = integrate( boundaryfaces(mesh), IM<2,Order+1,value_type,Simplex>(), trans( idv( u ) )*(Ny()*oneX()-Nx()*oneY() ) ).evaluate();
        std::cout << "int_ut = " << int_ut << "\n";
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
        using namespace Life;
        using namespace Life::vf;


        typedef typename imesh<value_type>::type mesh_type;
        typename imesh<value_type>::ptrtype mesh( createMesh<value_type>( meshSize ) );

        const value_type eps = 1000*Life::type_traits<value_type>::epsilon();

        typedef fusion::vector<Lagrange<Order, Vectorial>,
            Lagrange<Order-1, Scalar> > basis_type;
        typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
        boost::shared_ptr<space_type> Xh( new space_type(mesh) );
        typename space_type::element_type u( Xh );

        u.template element<0>() = project( Xh->template functionSpace<0>(), elements(mesh), P() );
        u.template element<1>() = project( Xh->template functionSpace<1>(), elements(mesh), constant(1) );
        std::cout << "int(proj P() = " << integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), idv( u.template element<0>() ) ).evaluate() << "\n";
        std::cout << "int(1 = " << integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), idv( u.template element<1>() ) ).evaluate() << "\n";

        ublas::matrix<value_type> m(  integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), gradv( u.template element<0>() ) ).evaluate() );
        std::cout << "int(grad(P()) = " << m << "\n";
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( m(0,0)-4, std::pow(10.0,-2.0*Order) );
        BOOST_CHECK_SMALL( m(0,1), std::pow(10.0,-2.0*Order) );
        BOOST_CHECK_SMALL( m(1,0), std::pow(10.0,-2.0*Order) );
        BOOST_CHECK_SMALL( m(1,1)-4, std::pow(10.0,-2.0*Order) );
#endif
        m = integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), trace( gradv( u.template element<0>() )*trans(gradv( u.template element<0>() ) ) ) ).evaluate();
        std::cout << "int(grad(P()*grad^T(P())) = " << m << "\n";

        m = integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), trace( trans(gradv( u.template element<0>() ))*(gradv( u.template element<0>() ) ) ) ).evaluate();
        std::cout << "int(grad(P()^T*grad(P())) = " << m << "\n";

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

        u.template element<0>() = project( Xh->template functionSpace<0>(), elements(mesh), u_exact );
        u.template element<1>() = project( Xh->template functionSpace<1>(), elements(mesh), Px()+Py() );
        std::cout << "int(u) = " << integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), idv( u.template element<0>() ) ).evaluate() << "\n";
        std::cout << "int(u - u_exact) = " << integrate( elements(mesh), IM<2,Order,value_type,Simplex>(),
                                                     trans(idv( u.template element<0>() )-u_exact)*(idv( u.template element<0>() )-u_exact) ).evaluate() << "\n";

        std::cout << "int(u<1>) = " << integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), idv( u.template element<1>() ) ).evaluate() << "\n";

        m= integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), gradv( u.template element<0>() ) ).evaluate();
        std::cout << "int(grad(u)) = " << m << "\n";

        m= integrate( boundaryfaces(mesh), IM<2,Order,value_type,Simplex>(), gradv( u.template element<0>() )*N() ).evaluate();
        std::cout << "int_bfaces(grad_u*N()) = " << m << "\n";

        m= integrate( boundaryfaces(mesh), IM<2,Order,value_type,Simplex>(), grad_exact*N() ).evaluate();
        std::cout << "int_bfaces(grad_exact*N) = " << m << "\n";

        m= integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), div_grad_exact ).evaluate();
        std::cout << "int(div_grad_exact) = " << m << "\n";

        m= integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), grad_exact  ).evaluate();
        std::cout << "int(grad_exact) = " << m << "\n";
        m= integrate( elements(mesh), IM<2,Order,value_type,Simplex>(), gradv( u.template element<0>() ) - grad_exact ).evaluate();
        std::cout << "int( grad(u)-grad_exact) = " << m << "\n";


        m= integrate( elements(mesh), IM<2,Order,value_type,Simplex>(),
                      trace( (gradv( u.template element<0>())-grad_exact)*trans(gradv( u.template element<0>())-grad_exact) )
                      ).evaluate();
        std::cout << "|grad(u)-grad_exact|_0^2 = " << m << "\n";




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
    test->add( BOOST_TEST_CASE( ( test_integration_domain<double>( mpi->vm()["hsize"].as<double>() ) ) ) );
    test->add( BOOST_TEST_CASE( ( test_integration_boundary<double>( mpi->vm()["hsize"].as<double>() ) ) ) );
    test->add( BOOST_TEST_CASE( ( test_integration_functions<2,double>( mpi->vm()["hsize"].as<double>() ) ) ) );
    test->add( BOOST_TEST_CASE( ( test_integration_vectorial_functions<2,double>(mpi->vm()["hsize"].as<double>() ) ) ) );
    test->add( BOOST_TEST_CASE( ( test_integration_composite_functions<2,double>(mpi->vm()["hsize"].as<double>() ) ) ) );

    return test;
}
#else
int
main( int argc, char** argv )
{
    Life::Application mpi( argc, argv, makeAbout(), makeOptions() );
    Life::Assert::setLog( "test_integration.assert");

    test_integration_circle<double> t1( mpi.vm()["hsize"].as<double>() ); t1();
    test_integration_domain<double> t2( mpi.vm()["hsize"].as<double>() ); t2();
    test_integration_boundary<double> t3( mpi.vm()["hsize"].as<double>() ); t3();
    test_integration_functions<2,double> t4( mpi.vm()["hsize"].as<double>() ); t4();
    test_integration_vectorial_functions<2,double> t5(mpi.vm()["hsize"].as<double>() ); t5();
    test_integration_composite_functions<2,double> t6(mpi.vm()["hsize"].as<double>() ); t6();
}
#endif // USE_BOOST_TEST
