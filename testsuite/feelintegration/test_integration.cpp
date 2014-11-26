/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-08-25

  Copyright (C) 2006 EPFL
  Copyright (C) 2006-2009 Universite Joseph Fourier (Grenoble I)

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


// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE integration testsuite
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN

#include <testsuite/testsuite.hpp>

#include <feel/options.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/savegmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/geotool.hpp>

//#define USE_BOOST_TEST 1
//#undef USE_BOOST_TEST

const double DEFAULT_MESH_SIZE=0.05;

namespace Feel
{
struct f_Px
{
    static const size_type context = vm::JACOBIAN|vm::POINT;
    typedef double value_type;
    typedef value_type evaluate_type;
    typedef Feel::uint16_type uint16_type;
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
    typedef Feel::uint16_type uint16_type;
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
    typedef Feel::uint16_type uint16_type;
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
    typedef Feel::uint16_type uint16_type;
    static const uint16_type rank = 0;
    static const uint16_type imorder = 2;
    static const bool imIsPoly = false;
    double operator()( uint16_type, uint16_type, ublas::vector<double> const& x, ublas::vector<double> const& /*n*/ ) const
    {
        return math::sin( x[0] );
    }
};

#if 0
#include <matheval.h>
struct f_matheval
{
    static const size_type context = vm::JACOBIAN|vm::POINT;
    typedef double value_type;
    typedef Feel::uint16_type uint16_type;
    static const uint16_type rank = 0;
    static const uint16_type imorder = 2;
    static const bool imIsPoly = false;
    double operator()( uint16_type, uint16_type, ublas::vector<double> const& x, ublas::vector<double> const& /*n*/ ) const
    {
        /* introduce matheval function */
        void *f, *f_prim;		/* Evaluators for function and function derivative.  */
        char **names;			/* Function variables names. */
        int count;			/* Number of function variables. */

        f = evaluator_create ("sin(x)");
        evaluator_get_variables (f, &names, &count);
        double val = evaluator_evaluate_x (f, x[0]);
        evaluator_destroy (f);
        return val;
    }
};
#endif
template<typename T, int Dim, int Order = 1>
struct imesh
{
    typedef Simplex<Dim, Order> convex_type;
    typedef Mesh<convex_type, T > type;
    typedef boost::shared_ptr<type> ptrtype;
};

template<typename T>
typename imesh<T,2>::ptrtype
createMesh( double hsize )
{
    double meshSize = hsize;
    //BOOST_TEST_MESSAGE( "hsize = " << meshSize << std::endl;

    Gmsh __gmsh;
    std::string fname;
    std::ostringstream ostr;
    std::ostringstream nameStr;

    //std::cout <<"Mesh generation ... ";

    GeoTool::Node x1(-1, -1);
    GeoTool::Node x2( 1,  1);
    GeoTool::Rectangle R( meshSize,"MyRectangle",x1,x2);
    R.setMarker(_type="line",_name="Gamma1",_marker1=true,_marker3=true);
    R.setMarker(_type="line",_name="Gamma2",_marker2=true,_marker4=true);
    R.setMarker(_type="surface",_name="Omega",_markerAll=true);

    auto mesh = R.createMesh(_mesh = new typename imesh<T,2>::type,
                             _name="square" );

    return mesh;
}

template<typename value_type = double>
struct test_integration_circle: public Application
{
    typedef typename imesh<value_type,2,4>::convex_type convex_type;
    typedef typename imesh<value_type,2,4>::type mesh_type;
    typedef typename imesh<value_type,2,4>::ptrtype mesh_ptrtype;
    typedef FunctionSpace<mesh_type, bases<Lagrange<4, Scalar> >, double> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef bases<Lagrange<4, Vectorial> > vector_basis_type;
    typedef FunctionSpace<mesh_type, vector_basis_type, value_type> vector_space_type;

    test_integration_circle()
        :
        Application(),
        M_backend( Backend<double>::build( soption( _name="backend" ) ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( "ellipsoid" ),
        mesh()
    {
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % 2 ).str() ,
                                             _usenames=true,
                                             _convex=( convex_type::is_hypercube )?"Hypercube":"Simplex",
                                             _shape=shape,
                                             _dim=2,
                                             _order=4,
                                             _xmin=-1.,_ymin=-1.,
                                             _h=meshSize ),
                               _update=MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES );
    }

    void operator()()
    {
        using namespace Feel;
        using namespace Feel::vf;

        saveGMSHMesh(_mesh=mesh,_filename="ellipsoid.msh");

        int Order=2;
        double t = 0.0;
        AUTO( mycst, cst_ref( t ) );

        t = 1.0;
        //value_type v0 = integrate( elements(mesh), mycst, _Q<2>() ).evaluate()( 0, 0 );
        value_type v0 = integrate( _range=elements( mesh ), _expr=mycst,_quad=_Q<10>() ).evaluate()( 0, 0 );
        value_type v0x = integrate( _range=elements( mesh ), _expr=mycst, _geomap=GeomapStrategyType::GEOMAP_O1 ).evaluate()( 0, 0 );
        value_type v0y = integrate( _range=elements( mesh ), _expr=mycst, _geomap=GeomapStrategyType::GEOMAP_OPT).evaluate()( 0, 0 );
        value_type v0z = integrate( _range=elements( mesh ), _expr=mycst, _geomap=GeomapStrategyType::GEOMAP_HO ).evaluate()( 0, 0 );
        BOOST_TEST_MESSAGE( "v0=" << v0 << "\n" );
        value_type v00 = ( integrate( _range=boundaryelements( mesh ), _expr=mycst, _quad=_Q<10>() ).evaluate()( 0, 0 )+
                           integrate( _range=internalelements( mesh ), _expr=mycst, _quad=_Q<10>() ).evaluate()( 0, 0 ) );
        BOOST_TEST_MESSAGE( "v00=" << v00 << "\n" );

        BOOST_TEST_MESSAGE( "[circle] v0 0 = " << integrate( boundaryfaces( mesh ),  N(), _Q<4>() ).evaluate() << "\n" );
        BOOST_TEST_MESSAGE( "[circle] v0 1 = " << integrate( boundaryfaces( mesh ),  vec( idf( f_Nx() ), idf( f_Ny() ) ), _Q<4>() ).evaluate() << "\n" );
        BOOST_TEST_MESSAGE( "[circle] v0 2 = " << integrate( boundaryfaces( mesh ),
                            trans( vec( constant( 1 ),Ny() ) )*one(), _Q<4>() ).evaluate() << "\n" );

        BOOST_TEST_MESSAGE( "[circle] v0 3 = " << integrate( boundaryfaces( mesh ),
                            mat<2,2>( Nx(),constant( 1 ),constant( 1 ),Ny() )*vec( constant( 1 ),constant( 1 ) ), _Q<4>() ).evaluate() << "\n" );
        BOOST_TEST_MESSAGE( "[circle] v0 4 (==v0 3) = " << integrate( boundaryfaces( mesh ),
                            mat<2,2>( idf( f_Nx() ),constant( 1 ),
                                      constant( 1 ),idf( f_Ny() ) )*vec( constant( 1 ),constant( 1 ) ), _Q<4>() ).evaluate() << "\n" );
        value_type pi = M_PI;//4.0*math::atan(1.0);
        const value_type eps = 1000*Feel::type_traits<value_type>::epsilon();
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v0, pi, 2e-1 );
        BOOST_CHECK_CLOSE( v0x, pi, 2e-1 );
        BOOST_CHECK_CLOSE( v0y, pi, 2e-1 );
        BOOST_CHECK_CLOSE( v0z, pi, 2e-1 );
        BOOST_CHECK_CLOSE( v0, v00, eps  );
        BOOST_CHECK_CLOSE( v00, pi, eps  );
#else
        FEELPP_ASSERT( math::abs( v0-pi ) < math::pow( meshSize, 2*Order ) )( v0 )( math::abs( v0-pi ) )( math::pow( meshSize, 2*Order ) ).warn ( "v0 != pi" );
        FEELPP_ASSERT( math::abs( v0-v00 ) < eps )( v0 )( v00 )( math::abs( v0-v00 ) )( eps ).warn ( "v0 != pi" );
#endif /* USE_BOOST_TEST */

        auto v1 = integrate( elements( mesh ), Px()*Px()+Py()*Py(), _Q<2>() ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v1, pi/2, 2e-1 );
        BOOST_CHECK_CLOSE( v1, v0/2, 2e-1 );
#endif

        boost::shared_ptr<space_type> Xh( new space_type( mesh ) );
        auto u = Xh->element();

        u = project( Xh, elements( mesh ), constant( 1.0 ) );
        v0 = integrate( elements( mesh ), idv( u ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_TEST_MESSAGE( "[circle] v0(~pi)=" << v0 << " pi=" << pi << " should be equal \n" );
        BOOST_CHECK_CLOSE( v0, pi, 2e-1 );
#else
        FEELPP_ASSERT( math::abs( v0-pi ) < math::pow( meshSize, 2*Order ) )( v0 )( math::abs( v0-pi ) )( math::pow( meshSize, 2*Order ) ).warn ( "v0 != pi" );
#endif /* USE_BOOST_TEST */

        u = project( Xh, elements( mesh ), Px()*Px()+Py()*Py() );
        v0 = integrate( elements( mesh ), idv( u ), _Q<2>() ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_TEST_MESSAGE( "[circle] v0(~pi/2)=" << v0 << " pi/2=" << pi/2 << " should be equal \n" );
        BOOST_CHECK_CLOSE( v0, pi/2, 2e-1 );
#endif /* USE_BOOST_TEST */


        boost::shared_ptr<vector_space_type> Xvh( new vector_space_type( mesh ) );
        auto U = Xvh->element();

        U = project( Xvh, elements( mesh ), vec( constant( 1.0 ),constant( 1.0 ) ) );
        v0 = integrate( boundaryfaces( mesh ), trans( idv( U ) )*N() ).evaluate()( 0, 0 );
        v00 = integrate( elements( mesh ), divv( U ) ).evaluate()( 0, 0 );

        BOOST_TEST_MESSAGE( "[circle] v0=" << v0 << " v00=" << v00 << " should be equal thanks to gauss int 1.N() != int div 1\n" );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v0, eps );
        BOOST_CHECK_SMALL( v00, eps );
#else
        FEELPP_ASSERT( math::abs( v0-v00 ) < eps )( v0 )( v00 )( math::abs( v0-v00 ) ).warn ( "int 1.N() != int div 1" );
#endif /* USE_BOOST_TEST */

        U = project( Xvh, elements( mesh ), vec( Px(),Py() ) );
        v0 = integrate( boundaryfaces( mesh ), trans( idv( U ) )*N() ).evaluate()( 0, 0 );
        v00 = integrate( elements( mesh ), divv( U ) ).evaluate()( 0, 0 );

        BOOST_TEST_MESSAGE( "[circle] v0=" << v0 << " v00=" << v00 << " should be equal thanks to gauss int (x,y).N() != int div (x,y)\n" );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v0, v00, eps );
#else
        FEELPP_ASSERT( math::abs( v0-v00 ) < eps )( v0 )( v00 )( math::abs( v0-v00 ) ).warn ( "int (x,y).N() != int div (x,y)" );
#endif /* USE_BOOST_TEST */

    }
    boost::shared_ptr<Feel::Backend<double> > M_backend;
    double meshSize;
    std::string shape;
    mesh_ptrtype mesh;
};
template<typename value_type = double>
struct test_integration_simplex: public Application
{
    typedef typename imesh<value_type,2>::convex_type convex_type;
    typedef typename imesh<value_type,2>::type mesh_type;
    typedef typename imesh<value_type,2>::ptrtype mesh_ptrtype;
    typedef FunctionSpace<mesh_type, bases<Lagrange<2, Scalar> >, double> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef bases<Lagrange<2, Vectorial> > vector_basis_type;
    typedef FunctionSpace<mesh_type, vector_basis_type, value_type> vector_space_type;

    test_integration_simplex()
        :
        Application(),
        M_backend( Backend<double>::build( soption( _name="backend" ) ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( "simplex" ),
        mesh()
    {
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % 2 ).str() ,
                                             _usenames=true,
                                             _convex=( convex_type::is_hypercube )?"Hypercube":"Simplex",
                                             _shape=shape,
                                             _dim=2,
                                             _xmin=-1.,_ymin=-1.,
                                             _h=meshSize ),
                               _update=MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES );
    }

    void operator()()
    {
        using namespace Feel;
        using namespace Feel::vf;


        const value_type eps = 1e-9;
        //typename imesh<value_type,2,1>::ptrtype mesh( createSimplex<value_type,1>( meshSize ) );
        //typedef typename imesh<value_type,2>::type mesh_type;

        typedef bases<Lagrange<3, Scalar> > basis_type;
        typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
        boost::shared_ptr<space_type> Xh( new space_type( mesh ) );
        typename space_type::element_type u( Xh );

        value_type meas = integrate( elements( mesh ), constant( 1.0 ) ).evaluate()( 0, 0 );
        value_type v0 = integrate( elements( mesh ), constant( 1.0 ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( meas, 2, eps );
#else
        CHECK( math::abs( v0-pi ) < 1e-2 ) << "v0=" <<  v0 << " error(v0-pi)=" << math::abs( v0-pi ) << " eps =" <<  eps;
#endif /* USE_BOOST_TEST */

        auto in1 = integrate( boundaryfaces( mesh ), N() ).evaluate();
        auto in2 = integrate( boundaryfaces( mesh ), vec( idf( f_Nx() ), idf( f_Ny() ) ), _Q<2>() ).evaluate();
#if defined( USE_BOOST_TEST )
        BOOST_CHECK_CLOSE( in1( 0,0 ), in2( 0,0 ), eps );
        BOOST_CHECK_CLOSE( in1( 1,0 ), in2( 1,0 ), eps );
#endif
        BOOST_TEST_MESSAGE( "[simplex] v0 0 = " << in1  << "\n" );
        BOOST_TEST_MESSAGE( "[simplex] v0 1 = " << in2  << "\n" );
        BOOST_TEST_MESSAGE( "[simplex] v0 2 = " << integrate( boundaryfaces( mesh ), trans( vec( constant( 1 ),Ny() ) )*one() ).evaluate() << "\n" );

        auto in3 = integrate( boundaryfaces( mesh ), mat<2,2>( Nx(),constant( 1 ),constant( 1 ),Ny() )*vec( constant( 1 ),constant( 1 ) ) ).evaluate();
        auto in4 = integrate( boundaryfaces( mesh ), mat<2,2>( idf( f_Nx() ),constant( 1 ),
                              constant( 1 ),idf( f_Ny() ) )*vec( constant( 1 ),constant( 1 ) ) ).evaluate();

#if defined( USE_BOOST_TEST )
        BOOST_CHECK_CLOSE( in3( 0,0 ), in4( 0,0 ), eps );
        BOOST_CHECK_CLOSE( in3( 1,0 ), in4( 1,0 ), eps );
#endif

        BOOST_TEST_MESSAGE( "[simplex] v0 3 = " << in3 << "\n" );
        BOOST_TEST_MESSAGE( "[simplex] v0 4 (==v0 3) = " <<  in4 << "\n" );

        // int ([-1,1],[-1,x]) 1 dx
        u = project( Xh, elements( mesh ), Px()+Py() );
        double v1 = integrate( elements( mesh ), gradv( u )*trans( gradv( u ) ) ).evaluate()( 0,0 ) ;

#if defined(USE_BOOST_TEST )
        BOOST_CHECK_CLOSE( v1, 2*meas, eps );
#endif
        u = project( Xh, elements( mesh ), Px()*Px() );
        double lapu = integrate( elements( mesh ),  trace( hessv( u ) ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST )
        BOOST_CHECK_CLOSE( lapu, 2*meas, eps );
#endif
        double lapu1 = integrate( elements( mesh ),  trace( hessv( u )*trans( hessv( u ) ) ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST )
        BOOST_CHECK_CLOSE( lapu1, 4*meas, eps );
#endif
        auto hessu = integrate( elements( mesh ), hessv( u ) ).evaluate();
#if defined(USE_BOOST_TEST )
        BOOST_CHECK_CLOSE( hessu( 0,0 ), 2*meas, eps );
        BOOST_CHECK_SMALL( hessu( 1,0 ), eps );
        BOOST_CHECK_SMALL( hessu( 0,1 ), eps );
        BOOST_CHECK_SMALL( hessu( 1,1 ), eps );
#endif
        u = project( Xh, elements( mesh ), Px()*Px()+Py()*Py() );
        lapu = integrate( elements( mesh ),  trace( hessv( u ) ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST )
        BOOST_CHECK_CLOSE( lapu, 4*meas, eps );
#endif
        lapu1 = integrate( elements( mesh ),  trace( hessv( u )*trans( hessv( u ) ) ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST )
        BOOST_CHECK_CLOSE( lapu1, 8*meas, eps );
#endif
        hessu = integrate( elements( mesh ), hessv( u ) ).evaluate();
#if defined(USE_BOOST_TEST )
        BOOST_CHECK_CLOSE( hessu( 0,0 ), 2*meas, eps );
        BOOST_CHECK_SMALL( hessu( 1,0 ), eps );
        BOOST_CHECK_SMALL( hessu( 0,1 ), eps );
        BOOST_CHECK_CLOSE( hessu( 1,1 ), 2*meas, eps );
#endif
        u = project( Xh, elements( mesh ), Px()*Px()+Py()*Py()+3*Px()*Py() );
        lapu = integrate( elements( mesh ),  trace( hessv( u ) ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST )
        BOOST_CHECK_CLOSE( lapu, 4*meas, eps );
#endif
        lapu1 = integrate( elements( mesh ),  trace( hessv( u )*trans( hessv( u ) ) ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST )
        BOOST_CHECK_CLOSE( lapu1, 26*meas, eps );
#endif
        hessu = integrate( elements( mesh ), hessv( u ) ).evaluate();
#if defined(USE_BOOST_TEST )
        BOOST_CHECK_CLOSE( hessu( 0,0 ), 2*meas, eps );
        BOOST_CHECK_CLOSE( hessu( 1,0 ), 3*meas, eps );
        BOOST_CHECK_CLOSE( hessu( 0,1 ), 3*meas, eps );
        BOOST_CHECK_CLOSE( hessu( 1,1 ), 2*meas, eps );
#endif

        typedef bases<Lagrange<3, Vectorial> > v_basis_type;
        typedef FunctionSpace<mesh_type, v_basis_type, value_type> v_space_type;
        boost::shared_ptr<v_space_type> Yh( new v_space_type( mesh ) );
        typename v_space_type::element_type v( Yh );

        auto p2 = Px()*Px()*Py()+Py()*Py()+cos( Px() );
        v = project( Yh, elements( mesh ), vec( p2,p2 ) );
        double divp2 = integrate( elements( mesh ), divv( v ), _Q<2>() ).evaluate()( 0, 0 );
        double unp2 = integrate( boundaryfaces( mesh ), trans( idv( v ) )*N(), _Q<2>() ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST )
        BOOST_CHECK_CLOSE( divp2, unp2, eps );
#endif
    }
    boost::shared_ptr<Feel::Backend<double> > M_backend;
    double meshSize;
    std::string shape;
    mesh_ptrtype mesh;
};
template<typename value_type = double>
struct test_integration_domain: public Application
{
    typedef typename imesh<value_type,2>::convex_type convex_type;
    typedef typename imesh<value_type,2>::type mesh_type;
    typedef typename imesh<value_type,2>::ptrtype mesh_ptrtype;
    typedef FunctionSpace<mesh_type, bases<Lagrange<2, Scalar> >, double> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef bases<Lagrange<2, Vectorial> > vector_basis_type;
    typedef FunctionSpace<mesh_type, vector_basis_type, value_type> vector_space_type;

    test_integration_domain()
        :
        Application(),
        M_backend( Backend<double>::build( soption( _name="backend" ) ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( "hypercube" ),
        mesh()
    {
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % 2 ).str() ,
                                             _usenames=true,
                                             _convex=( convex_type::is_hypercube )?"Hypercube":"Simplex",
                                             _shape=shape,
                                             _dim=2,
                                             _xmin=-1.,_ymin=-1.,
                                             _h=meshSize ),
                               _update=MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES );
    }

    void operator()()
    {
        using namespace Feel;
        using namespace Feel::vf;

        typename imesh<value_type,2>::ptrtype mesh( createMesh<value_type>( meshSize ) );

        const value_type eps = 1e-9;

        // int ([-1,1],[-1,x]) 1 dx
        value_type v0 = integrate( elements( mesh ),
                                   vf::min( constant( 1.0 ),constant( 2.0 ) ) ).evaluate()( 0, 0 );
        value_type v00 = ( integrate( boundaryelements( mesh ),
                                      vf::min( constant( 1.0 ),constant( 2.0 ) ) ).evaluate()( 0, 0 )+
                           integrate( internalelements( mesh ),
                                      vf::min( constant( 1.0 ),constant( 2.0 ) ) ).evaluate()( 0, 0 ) );

        BOOST_TEST_MESSAGE( "[domain] v0 = " << v0 << "\n" );
        BOOST_TEST_MESSAGE( "[domain] v00 = " << v00 << "\n" );
        BOOST_TEST_MESSAGE( "[domain] int(1*N()) = " << integrate( boundaryfaces( mesh ),
                            constant( 1.0 )*N() ).evaluate() << "\n" );
        BOOST_TEST_MESSAGE( "[domain] int(tr(N())**N()) = " << integrate( boundaryfaces( mesh ),
                            N()*trans( N() ) ).evaluate() << "\n" );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v0, 4.0, eps );
        BOOST_CHECK_CLOSE( v0, v00, eps );
#else
        FEELPP_ASSERT( math::abs( v0-4.0 ) < eps )( v0 )( math::abs( v0-4.0 ) )( eps ).warn ( "v0 != 4" );
        FEELPP_ASSERT( math::abs( v0-v00 ) < eps )( v0 )( v00 )( math::abs( v0-v00 ) )( eps ).warn ( "v0 != v00" );
#endif /* USE_BOOST_TEST */


        // int ([-1,1],[-1,1]) x dx
        value_type v1 = integrate( elements( mesh ), Px() ).evaluate()( 0, 0 );
        value_type v11 = integrate( elements( mesh ), idf( f_Px() ) ).evaluate()( 0, 0 );
        value_type vx2 = integrate( elements( mesh ), Px() * idf( f_Px() ) ).evaluate()( 0, 0 );
        std::cout << "vx2=" << vx2 << "\n" << std::flush;

        BOOST_TEST_MESSAGE( "[domain] int(P()) = " << integrate( elements( mesh ),
                            P() ).evaluate() << "\n" );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v1- 0.0, eps );
        BOOST_CHECK_SMALL( v11- 0.0, eps );
#else
        FEELPP_ASSERT( math::abs( v1-0.0 ) < eps )( v1 )( math::abs( v1-0.0 ) )( eps ).warn ( "v1 != 0" );
        FEELPP_ASSERT( math::abs( v11-0.0 ) < eps )( v11 )( math::abs( v11-0.0 ) )( eps ).warn ( "v11 != 0" );
#endif /* USE_BOOST_TEST */


        // int ([-1,1],[-1,1]) abs(x) dx
        value_type vsin = integrate( elements( mesh ), sin( Px() ), _Q<5>() ).evaluate()( 0, 0 );
        value_type vsin1 = integrate( elements( mesh ), idf( f_sinPx() ), _Q<5>() ).evaluate()( 0, 0 );
        std::cout << vsin << " " << vsin1 << " [f_sinPx()] " << std::flush;
#if 0
        value_type vmatheval = integrate( elements( mesh ), idf( f_matheval() ), _Q<5>() ).evaluate()( 0, 0 );
        std::cout << vmatheval << " [matheval]\n";
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( vsin- 2*( -math::cos( 1.0 )+math::cos( -1.0 ) ), eps );
        BOOST_CHECK_SMALL( vsin1- 2*( -math::cos( 1.0 )+math::cos( -1.0 ) ), eps );
#else
        FEELPP_ASSERT( math::abs( vsin-2.0*( -math::cos( 1.0 )+math::cos( -1.0 ) ) ) < eps )
        ( vsin )
        ( math::abs( vsin-2.0*( -math::cos( 1.0 )+math::cos( -1.0 ) ) ) )( eps ).warn ( "vsin != 2*(cos(1)-cos(-1))" );
        FEELPP_ASSERT( math::abs( vsin1-2.0*( -math::cos( 1.0 )+math::cos( -1.0 ) ) ) < eps )
        ( vsin1 )
        ( math::abs( vsin1-2.0*( -math::cos( 1.0 )+math::cos( -1.0 ) ) ) )( eps ).warn ( "vsin1 != 2*(cos(1)-cos(-1))" );
#endif /* USE_BOOST_TEST */
#endif
        // int ([-1,1],[-1,1]) abs(x) dx
        value_type vabs = integrate( elements( mesh ), abs( Px() )+abs( Py() ), _Q<5>() ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( vabs, 4.0, 1e-2 );
#else
        FEELPP_ASSERT( math::abs( vabs-4.0 ) < 1e-2 )( vabs )( math::abs( vabs-4.0 ) )( 1e-2 ).warn ( "vabs != 4" );
#endif /* USE_BOOST_TEST */



        // int (\partial ([-1,1],[-1,1]) 1 dx
        value_type v2 = integrate( boundaryfaces( mesh ), constant( 1.0 ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v2, 8.0, eps );
#else
        FEELPP_ASSERT( math::abs( v2-8.0 ) < eps )( v2 )( math::abs( v2-8.0 ) )( eps ).warn ( "v2 != 8" );
#endif /* USE_BOOST_TEST */

    }
    boost::shared_ptr<Feel::Backend<double> > M_backend;
    double meshSize;
    std::string shape;
    mesh_ptrtype mesh;
};

template<typename value_type = double>
struct test_integration_boundary: public Application
{
    typedef typename imesh<value_type,2>::convex_type convex_type;
    typedef typename imesh<value_type,2>::type mesh_type;
    typedef typename imesh<value_type,2>::ptrtype mesh_ptrtype;
    typedef FunctionSpace<mesh_type, bases<Lagrange<2, Scalar> >, double> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef bases<Lagrange<2, Vectorial> > vector_basis_type;
    typedef FunctionSpace<mesh_type, vector_basis_type, value_type> vector_space_type;

    test_integration_boundary()
        :
        Application(),
        M_backend( Backend<double>::build( soption( _name="backend" ) ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( "hypercube" ),
        mesh()
    {
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _desc=domain( _name=( boost::format( "test_integration_boundary_%1%-%2%" ) % shape % 2 ).str() ,
                                             _usenames=true,
                                             _convex=( convex_type::is_hypercube )?"Hypercube":"Simplex",
                                             _shape=shape,
                                             _dim=2,
                                             _xmin=-1.,_ymin=-1.,
                                             _h=meshSize ),
                               _update=MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES );
    }
    void operator()()
    {
        using namespace Feel;
        using namespace Feel::vf;


        const value_type eps = 1000*Feel::type_traits<value_type>::epsilon();

        // int (\partial ([-1,1],[-1,1]) 1 dx
        value_type v2 = integrate( boundaryfaces( mesh ), constant( 1.0 ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v2, 8.0, eps );
#else
        FEELPP_ASSERT( math::abs( v2-8.0 ) < eps )( v2 )( math::abs( v2-8.0 ) )( eps ).warn ( "v2 != 8" );
#endif /* USE_BOOST_TEST */

        // int (\partial ([-1,1],[-1,1]) x * y dx dy
        value_type v3 = integrate( boundaryfaces( mesh ), Px()*Py() ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v3- 0.0, eps );
#else
        FEELPP_ASSERT( math::abs( v3-0.0 ) < eps )( v3 )( math::abs( v3-0.0 ) )( eps ).warn ( "v3 != 0" );
#endif /* USE_BOOST_TEST */

        // int (\partial ([-1,1],[-1,1]) -2*y dx dy
        value_type v4 = integrate( boundaryfaces( mesh ), -2*Py() ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v4- 0.0, eps );
#else
        FEELPP_ASSERT( math::abs( v4-0.0 ) < eps )( v4 )( math::abs( v4-0.0 ) )( eps ).warn ( "v4 != 0" );
#endif /* USE_BOOST_TEST */

        // int_10 exp(Py) dx
        value_type v5 = integrate( markedfaces( mesh,"Neumann" ), exp( Py() ) ).evaluate()( 0, 0 );
        value_type v5_ex = 2.0*( math::exp( 1.0 )+math::exp( -1.0 ) );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v5, v5_ex, eps );
#else
        FEELPP_ASSERT( math::abs( v5-v5_ex ) < eps )( v5 )( v5_ex )( math::abs( v5-v5_ex ) )( eps ).warn ( "v5 != v5_ex" );
#endif /* USE_BOOST_TEST */

        // int_20 exp(Py) dx
        value_type v6 = integrate( markedfaces( mesh,"Dirichlet" ), exp( Py() ), _Q<5>() ).evaluate()( 0, 0 );
        value_type v6_ex = 2.0*( math::exp( 1.0 )-math::exp( -1.0 ) );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v6, v6_ex, eps );
#else
        FEELPP_ASSERT( math::abs( v6-v6_ex ) < eps )( v6 )( v6_ex )( math::abs( v6-v6_ex ) )( eps ).warn ( "v6 != v6_ex" );
#endif /* USE_BOOST_TEST */


    }
    boost::shared_ptr<Feel::Backend<double> > M_backend;
    double meshSize;
    std::string shape;
    mesh_ptrtype mesh;
};

template<int Order, typename value_type = double>
struct test_integration_functions: public Application
{
    typedef typename imesh<value_type,2>::convex_type convex_type;
    typedef typename imesh<value_type,2>::type mesh_type;
    typedef typename imesh<value_type,2>::ptrtype mesh_ptrtype;
    typedef FunctionSpace<mesh_type, bases<Lagrange<Order, Scalar> >, double> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef bases<Lagrange<Order, Vectorial> > vector_basis_type;
    typedef FunctionSpace<mesh_type, vector_basis_type, value_type> vector_space_type;

    test_integration_functions()
        :
        Application(),
        M_backend( Backend<double>::build( soption( _name="backend" ) ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( "hypercube" ),
        mesh()
    {
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % 2 ).str() ,
                                             _usenames=true,
                                             _convex=( convex_type::is_hypercube )?"Hypercube":"Simplex",
                                             _shape=shape,
                                             _dim=2,
                                             _xmin=-1.,_ymin=-1.,
                                             _h=meshSize ),
                               _update=MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES );
    }
    void operator()()
    {
        using namespace Feel;
        using namespace Feel::vf;

        const value_type eps = 1e-9;

        boost::shared_ptr<space_type> Xh( new space_type( mesh ) );
        typename space_type::element_type u( Xh );

        // int ([-1,1],[-1,x]) 1 dx
        u = project( Xh, elements( mesh ), constant( 1.0 ) );
        value_type v0 = integrate( elements( mesh ), idv( u ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v0, 4.0, eps );
#else
        FEELPP_ASSERT( math::abs( v0-4.0 ) < eps )( v0 )( math::abs( v0-4.0 ) )( eps ).warn ( "v0 != 4" );
#endif /* USE_BOOST_TEST */

        //
        u = project( Xh, elements( mesh ), Px() );
        value_type v1 = integrate( elements( mesh ), idv( u ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v1- 0.0, eps );
#else
        FEELPP_ASSERT( math::abs( v1-0.0 ) < eps )( v1 )( math::abs( v1-0.0 ) )( eps ).warn ( "v1 != 0" );
#endif /* USE_BOOST_TEST */
        value_type v2 = integrate( elements( mesh ), dxv( u ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v2, 4.0, eps );
#else
        FEELPP_ASSERT( math::abs( v2-4.0 ) < eps )( v2 )( math::abs( v2-4.0 ) )( eps ).warn ( "v2 != 4" );
#endif /* USE_BOOST_TEST */

        value_type v3 = integrate( elements( mesh ), dyv( u ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v3- 0.0, eps );
#else
        FEELPP_ASSERT( math::abs( v3-0.0 ) < eps )( v3 )( math::abs( v3-0.0 ) )( eps ).warn ( "v3 != 0" );
#endif /* USE_BOOST_TEST */

        u = project( Xh, elements( mesh ), exp( Px() )*exp( Py() ) );
        value_type v4 = integrate( elements( mesh ), idv( u ), _Q<2>() ).evaluate()( 0, 0 );
        value_type v4_ex = ( math::exp( 1.0 )-math::exp( -1.0 ) )*( math::exp( 1.0 )-math::exp( -1.0 ) );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v4, v4_ex, std::pow( 10.0,-2.0*Order ) );
#else
        FEELPP_ASSERT( math::abs( v4-v4_ex ) < std::pow( 10.0,-2.0*Order ) )
        ( v4 )
        ( v4_ex )
        ( math::abs( v4-v4_ex ) )( std::pow( 10.0,-2.0*Order ) ).warn ( "v4 != v4_ex" );
#endif /* USE_BOOST_TEST */
        value_type v4_x = integrate( elements( mesh ), dxv( u ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v4_x, v4_ex, std::pow( 10.0,-2.0*Order ) );
#else
        FEELPP_ASSERT( math::abs( v4_x-v4_ex ) < std::pow( 10.0,-2.0*Order ) )
        ( v4_x )
        ( v4_ex )
        ( math::abs( v4_x-v4_ex ) )( std::pow( 10.0,-2.0*Order ) ).warn ( "v4_x != v4_ex" );
#endif /* USE_BOOST_TEST */

        value_type v4_y = integrate( elements( mesh ), dyv( u ) ).evaluate()( 0, 0 );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v4_y, v4_ex, std::pow( 10.0,-2.0*Order ) );
#else
        FEELPP_ASSERT( math::abs( v4_y-v4_ex ) < std::pow( 10.0,-2.0*Order ) )
        ( v4_y )
        ( v4_ex )
        ( math::abs( v4_y-v4_ex ) )( std::pow( 10.0,-2.0*Order ) ).warn ( "v4_y != v4_ex" );
#endif /* USE_BOOST_TEST */

        value_type v5 = integrate( markedfaces( mesh,"Neumann" ), idv( u ), _Q<2>() ).evaluate()( 0, 0 );
        value_type v5_ex = ( math::exp( 1.0 )-math::exp( -1.0 ) )*( math::exp( 1.0 )+math::exp( -1.0 ) );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( v5, v5_ex, std::pow( 10.0,-2.0*Order ) );
#else
        FEELPP_ASSERT( math::abs( v5-v5_ex ) < std::pow( 10.0,-2.0*Order ) )
        ( v5 )
        ( v5_ex )
        ( math::abs( v5-v5_ex ) )( std::pow( 10.0,-2.0*Order ) ).warn ( "v5 != v5_ex" );
#endif /* USE_BOOST_TEST */

#if 1
        double meas = integrate( elements( mesh ), cst( 1. ) ).evaluate()( 0, 0 );
        u = project( Xh, elements( mesh ), Px()*Px() + Py()*Py() +Px()*Py()  );
        auto hess6 = integrate( elements( mesh ), hessv( u ) ).evaluate();
#if defined(USE_BOOST_TEST)
        BOOST_TEST_MESSAGE( "hess6 =" << hess6 << "\n" );
        BOOST_CHECK_CLOSE( hess6( 0,0 ), 2*meas, eps );
        BOOST_CHECK_CLOSE( hess6( 0,1 ), meas, eps );
        BOOST_CHECK_CLOSE( hess6( 1,0 ), meas, eps );
        BOOST_CHECK_CLOSE( hess6( 1,0 ), hess6( 0, 1 ), eps );
        BOOST_CHECK_CLOSE( hess6( 1,1 ), 2*meas, eps );
#else
        // TODO
#endif /* USE_BOOST_TEST */
#endif


#if 0
        BOOST_TEST_MESSAGE( "idv(u) = " << integrate( markedfaces( mesh,"Neumann" ), idv( u ) ).evaluate() << "\n" );
        BOOST_TEST_MESSAGE( "gradv(u) = " << integrate( markedfaces( mesh,"Neumann" ), gradv( u ) ).evaluate() << "\n" );
        BOOST_TEST_MESSAGE( "trans(gradv(u))*N() = " << integrate( markedfaces( mesh, "Neumann" ), trans( gradv( u ) )*N() ).evaluate() << "\n" );
        BOOST_TEST_MESSAGE( "u*Nx()+u*Ny() = " << integrate( markedfaces( mesh, "Neumann" ), idv( u )*Nx() + idv( u )*Ny() ).evaluate() << "\n" );
        value_type v6 = integrate( markedfaces( mesh, "Neumann" ), trans( gradv( u ) )*N() ).evaluate()( 0, 0 );
        value_type v6_ex = ( math::exp( 1.0 )-math::exp( -1.0 ) )*( math::exp( 1.0 )+math::exp( -1.0 ) );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( v6-v6_ex, std::pow( 10.0,-2.0*Order ) );
#else
        FEELPP_ASSERT( math::abs( v6-v6_ex ) < std::pow( 10.0,-2.0*Order ) )
        ( v6 )
        ( v6_ex )
        ( math::abs( v6-v6_ex ) )( std::pow( 10.0,-2.0*Order ) ).warn ( "v6 != v6_ex" );
#endif /* USE_BOOST_TEST */
#endif
    }
    boost::shared_ptr<Feel::Backend<double> > M_backend;
    double meshSize;
    std::string shape;
    mesh_ptrtype mesh;
};

template<int Order, typename value_type = double>
struct test_integration_vectorial_functions: public Application
{
    typedef typename imesh<value_type,2>::convex_type convex_type;
    typedef typename imesh<value_type,2>::type mesh_type;
    typedef typename imesh<value_type,2>::ptrtype mesh_ptrtype;
    typedef FunctionSpace<mesh_type, bases<Lagrange<Order, Vectorial> >, double> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    test_integration_vectorial_functions()
        :
        Application(),
        M_backend( Backend<double>::build( soption( _name="backend" ) ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( "hypercube" ),
        mesh()
    {
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % 2 ).str() ,
                                             _usenames=true,
                                             _convex=( convex_type::is_hypercube )?"Hypercube":"Simplex",
                                             _shape=shape,
                                             _dim=2,
                                             _xmin=-1.,_ymin=-1.,
                                             _h=meshSize ),
                               _update=MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES );
    }
    void operator()()
    {
        using namespace Feel;
        using namespace Feel::vf;

        const value_type eps = 1000*Feel::type_traits<value_type>::epsilon();

        boost::shared_ptr<space_type> Xh( new space_type( mesh ) );
        typename space_type::element_type u( Xh );

        u = project( Xh, elements( mesh ), P() );
        BOOST_TEST_MESSAGE( "int(proj P() = " << integrate( elements( mesh ), idv( u ) ).evaluate() << "\n" );
        u = project( Xh, elements( mesh ), one() );
        BOOST_TEST_MESSAGE( "int(proj one() = " << integrate( elements( mesh ), idv( u ) ).evaluate() << "\n" );

        u = project( Xh, elements( mesh ), P() );
        auto m=  integrate( elements( mesh ), gradv( u ) ).evaluate();
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( m( 0,0 ), 4, std::pow( 10.0,-2.0*Order ) );
        BOOST_CHECK_SMALL( m( 0,1 )- 0, std::pow( 10.0,-2.0*Order ) );
        BOOST_CHECK_SMALL( m( 1,0 )- 0, std::pow( 10.0,-2.0*Order ) );
        BOOST_CHECK_CLOSE( m( 1,1 ), 4, std::pow( 10.0,-2.0*Order ) );
#endif
        auto mx= integrate( elements( mesh ), dxv( u ) ).evaluate();
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( mx( 0,0 ), 4, std::pow( 10.0,-2.0*Order ) );
        BOOST_CHECK_SMALL( mx( 1,0 )- 0, std::pow( 10.0,-2.0*Order ) );
#endif

        auto my = integrate( elements( mesh ), dyv( u ) ).evaluate();
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( my( 0,0 )- 0, std::pow( 10.0,-2.0*Order ) );
        BOOST_CHECK_CLOSE( my( 1,0 ), 4, std::pow( 10.0,-2.0*Order ) );
#endif

        u = project( Xh, elements( mesh ), Py()*oneX() + Px()*oneY() );
        auto int_divu = integrate( elements( mesh ), divv( u ) ).evaluate();
#if defined(USE_BOOST_TEST)
        value_type norm_int_divu = int_divu.norm();
        BOOST_CHECK_SMALL( norm_int_divu, eps );
#endif
        u = project( Xh, elements( mesh ), P() );
        int_divu = integrate( elements( mesh ), divv( u ) ).evaluate();
#if defined(USE_BOOST_TEST)
        norm_int_divu = int_divu.norm();
        BOOST_CHECK_CLOSE( norm_int_divu, 8, eps );
#endif

        // check the divergence theorem
        u = project( Xh, elements( mesh ), P() );
        int_divu = integrate( elements( mesh ), divv( u ) ).evaluate();
        BOOST_TEST_MESSAGE( "int_divu = " << int_divu << "\n" );
        auto int_un = integrate( boundaryfaces( mesh ), trans( idv( u ) )*N() ).evaluate();
        BOOST_TEST_MESSAGE( "(1, N) = " << integrate( boundaryfaces( mesh ), trans( one() )*N() ).evaluate() << "\n" );
        BOOST_TEST_MESSAGE( "(P, N) = " << integrate( boundaryfaces( mesh ), trans( P() )*N() ).evaluate() << "\n" );
        BOOST_TEST_MESSAGE( "int_un = " << int_un << "\n" );
#if defined(USE_BOOST_TEST)
        value_type norm_divergence = ( int_divu-int_un ).norm();
        BOOST_CHECK_SMALL( norm_divergence, eps );
#endif

        // check the stokes theorem
        auto int_curlu = integrate( elements( mesh ), curlzv( u ) ).evaluate()( 0,0 );
        BOOST_TEST_MESSAGE( "int_curlu = " << int_curlu << "\n" );
        auto int_ut = integrate( boundaryfaces( mesh ), trans( idv( u ) )*( Ny()*oneX()-Nx()*oneY() ) ).evaluate()( 0,0 );
        BOOST_TEST_MESSAGE( "int_ut = " << int_ut << "\n" );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( int_curlu, eps );
        BOOST_CHECK_SMALL( int_ut, eps );
#endif
    }
    boost::shared_ptr<Feel::Backend<double> > M_backend;
    double meshSize;
    std::string shape;
    mesh_ptrtype mesh;
};

template<int Order, typename value_type = double>
struct test_integration_matricial_functions: public Application
{
    typedef typename imesh<value_type,2>::convex_type convex_type;
    typedef typename imesh<value_type,2>::type mesh_type;
    typedef typename imesh<value_type,2>::ptrtype mesh_ptrtype;

    test_integration_matricial_functions()
        :
        Application(),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( "hypercube" ),
        mesh()
    {
        BOOST_TEST_MESSAGE( "[test_integration_matricial_functions::test_integration_matricial_functions]\n" );
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % 2 ).str() ,
                                             _usenames=true,
                                             _convex=( convex_type::is_hypercube )?"Hypercube":"Simplex",
                                             _shape=shape,
                                             _dim=2,
                                             _xmin=-1.,_ymin=-1.,
                                             _h=meshSize ),
                               _update=MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES );
    }
    void operator()()
    {
        BOOST_TEST_MESSAGE( "[test_integration_matricial_functions::operator()]\n" );
        using namespace Feel;
        using namespace Feel::vf;

        const value_type eps = 1000*Feel::type_traits<value_type>::epsilon();
        auto xxT = P()*trans( P() );
        BOOST_TEST_MESSAGE( "[test_integration_matricial_functions::operator()] v0\n" );
        auto v0 = integrate( _range=elements( mesh ), _expr=sym( xxT ) ).evaluate();
        auto v000 = integrate( _range=elements( mesh ), _expr=sym( xxT )( 0,0 ) ).evaluate()( 0,0 );
        auto v011 = integrate( _range=elements( mesh ), _expr=sym( xxT )( 1,1 ) ).evaluate()( 0,0 );
        auto v001 = integrate( _range=elements( mesh ), _expr=sym( xxT )( 0,1 ) ).evaluate()( 0,0 );
        auto v010 = integrate( _range=elements( mesh ), _expr=sym( xxT )( 1,0 ) ).evaluate()( 0,0 );
        decltype( v0 ) mat;
        mat << v000, v001, v010, v011;
        BOOST_TEST_MESSAGE( "[test_integration_matricial_functions::operator()] v1\n" );
        auto v1 = integrate( _range=elements( mesh ), _expr=0.5*( xxT+trans( xxT ) ) ).evaluate();
        auto v101 = integrate( _range=elements( mesh ), _expr=( 0.5*( xxT+trans( xxT ) ) )( 0,1 ) ).evaluate()( 0,0 );
        auto v110 = integrate( _range=elements( mesh ), _expr=( 0.5*( xxT+trans( xxT ) ) )( 1,0 ) ).evaluate()( 0,0 );
        BOOST_CHECK_SMALL( ( v0-v1 ).norm(), 1e-12 );
        BOOST_CHECK_SMALL( ( v0-mat ).norm(), 1e-12 );
        BOOST_CHECK_CLOSE( v001, v010, 1e-12 );
        BOOST_CHECK_CLOSE( v001, v110, 1e-12 );
        BOOST_CHECK_CLOSE( v001, v101, 1e-12 );

        BOOST_TEST_MESSAGE( "[sym check] v0=" << v0 << "\n v1=" << v1 << "\n mat=" << mat << "\n error=" << ( v0-v1 ).norm() << "\n" );
        BOOST_TEST_MESSAGE( "[test_integration_matricial_functions::operator()] v2\n" );
        auto v2 = integrate( _range=elements( mesh ), _expr=antisym( xxT ) ).evaluate();
        BOOST_TEST_MESSAGE( "[test_integration_matricial_functions::operator()] v3\n" );
        auto v3 = integrate( _range=elements( mesh ), _expr=0.5*( xxT-trans( xxT ) ) ).evaluate();
        BOOST_CHECK_SMALL( ( v2-v3 ).norm(), 1e-12 );
        BOOST_TEST_MESSAGE( "[antisym check] v2=" << v2 << " v3=" << v3 << " error=" << ( v2-v3 ).norm() << "\n" );
        BOOST_TEST_MESSAGE( "[test_integration_matricial_functions::operator()] done\n" );
    }
    double meshSize;
    std::string shape;
    mesh_ptrtype mesh;
};


template<int Order, typename value_type = double>
struct test_integration_composite_functions: public Application
{
    typedef typename imesh<value_type,2>::convex_type convex_type;
    typedef typename imesh<value_type,2>::type mesh_type;
    typedef typename imesh<value_type,2>::ptrtype mesh_ptrtype;
    typedef bases<Lagrange<Order, Vectorial>,Lagrange<Order-1, Scalar> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef typename space_type::element_type element_type;
    test_integration_composite_functions()
        :
        Application(),
        M_backend( Backend<double>::build( soption( _name="backend" ) ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( "hypercube" ),
        mesh()
    {
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % 2 ).str() ,
                                             _usenames=true,
                                             _convex=( convex_type::is_hypercube )?"Hypercube":"Simplex",
                                             _shape=shape,
                                             _dim=2,
                                             _xmin=-1.,_ymin=-1.,
                                             _h=meshSize ),
                               _update=MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES );
    }

    void operator()()
    {
        using namespace Feel;

        const value_type eps = 1000*Feel::type_traits<value_type>::epsilon();

        boost::shared_ptr<space_type> Xh( new space_type( mesh ) );
        element_type u( Xh );

        u.template element<0>() = project( _space=Xh->template functionSpace<0>(), _range=elements( mesh ), _expr=P() );
        u.template element<1>() = project( _space=Xh->template functionSpace<1>(), _range=elements( mesh ), _expr=constant( 1. ) );
        BOOST_TEST_MESSAGE( "int(proj P() = " << integrate( elements( mesh ), idv( u.template element<0>() ) ).evaluate() << "\n" );
        BOOST_TEST_MESSAGE( "int(1 = " << integrate( elements( mesh ), idv( u.template element<1>() ) ).evaluate() << "\n" );

        auto m(  integrate( elements( mesh ), gradv( u.template element<0>() ) ).evaluate() );
        BOOST_TEST_MESSAGE( "int(grad(P()) = " << m << "\n" );
#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( m( 0,0 )-4, std::pow( 10.0,-2.0*Order ) );
        BOOST_CHECK_SMALL( m( 0,1 ), std::pow( 10.0,-2.0*Order ) );
        BOOST_CHECK_SMALL( m( 1,0 ), std::pow( 10.0,-2.0*Order ) );
        BOOST_CHECK_SMALL( m( 1,1 )-4, std::pow( 10.0,-2.0*Order ) );
#endif
        auto m1 = integrate( elements( mesh ), trace( gradv( u.template element<0>() )*trans( gradv( u.template element<0>() ) ) ) ).evaluate();
        BOOST_TEST_MESSAGE( "int(grad(P()*grad^T(P())) = " << m1 << "\n" );

        auto m2 = integrate( elements( mesh ), trace( trans( gradv( u.template element<0>() ) )*( gradv( u.template element<0>() ) ) ) ).evaluate();
        BOOST_TEST_MESSAGE( "int(grad(P()^T*grad(P())) = " << m2 << "\n" );

#if 0
        AUTO( u_exact,( P()^( 2 ) )*( Px()+Py() ) );
        AUTO( grad_exact, ( mat<2,2>( 3*Px()*Px()+2*Px()*Py(), ( Px()^( 2 ) ), ( Py()^( 2 ) ), 3*Py()*Py()+2*Py()*Px() ) ) );
        AUTO( div_grad_exact, vec( 6*Px()+2*Py(), 6*Py()+2*Px() ) );
#else
        //AUTO( u_exact, vec(sin(Px())*sin(Py()), cos(Px())*cos(Py()) ) );
        AUTO( u_exact, sin( Px() )*sin( Py() )*oneX() + cos( Px() )*cos( Py() )*oneY() );
        AUTO( p_exact, val( cos( Px() )*sin( Py() ) ) );
        AUTO( du_dx, val( cos( Px() )*sin( Py() ) ) );
        AUTO( du_dy, val( sin( Px() )*cos( Py() ) ) );

        AUTO( dv_dx, val( -sin( Px() )*cos( Py() ) ) );
        AUTO( dv_dy, val( -cos( Px() )*sin( Py() ) ) );

        AUTO( grad_exact, ( mat<2,2>( du_dx, du_dy, dv_dx, dv_dy ) ) );
        AUTO( div_grad_exact, ( vec( -sin( Px() )*sin( Py() )-sin( Px() )*sin( Py() ),
                                     -cos( Px() )*cos( Py() )-cos( Px() )*cos( Py() ) ) ) );

#endif

        u.template element<0>() = project( Xh->template functionSpace<0>(), elements( mesh ), u_exact );
        u.template element<1>() = project( Xh->template functionSpace<1>(), elements( mesh ), Px()+Py() );
        BOOST_TEST_MESSAGE( "int(u) = " << integrate( elements( mesh ), idv( u.template element<0>() ) ).evaluate() << "\n" );
        BOOST_TEST_MESSAGE( "int(u - u_exact) = " << integrate( elements( mesh ),
                            trans( idv( u.template element<0>() )-u_exact )*( idv( u.template element<0>() )-u_exact ) ).evaluate() << "\n" );

        BOOST_TEST_MESSAGE( "int(u<1>) = " << integrate( elements( mesh ), idv( u.template element<1>() ) ).evaluate() << "\n" );

        m= integrate( elements( mesh ), gradv( u.template element<0>() ) ).evaluate();
        BOOST_TEST_MESSAGE( "int(grad(u)) = " << m << "\n" );
        auto du_dx1 = integrate( _range=elements( mesh ), _expr=gradv( u.template element<0>() )( 0,0 ) ).evaluate()( 0,0 );
        auto du_dy1 = integrate( _range=elements( mesh ), _expr=gradv( u.template element<0>() )( 0,1 ) ).evaluate()( 0,0 );
        auto du_dx2 = integrate( elements( mesh ), du_dx ).evaluate()( 0,0 );
        auto du_dy2 = integrate( elements( mesh ), du_dy ).evaluate()( 0,0 );
        BOOST_CHECK_SMALL( du_dx1,  1e-12 );
        BOOST_CHECK_SMALL( du_dy1, 1e-12 );
        BOOST_CHECK_SMALL( du_dx2,  1e-12 );
        BOOST_CHECK_SMALL( du_dy2, 1e-12 );
        auto dv_dx1 = integrate( _range=elements( mesh ), _expr=gradv( u.template element<0>() )( 1,0 ) ).evaluate()( 0,0 );
        auto dv_dy1 = integrate( _range=elements( mesh ), _expr=gradv( u.template element<0>() )( 1,1 ) ).evaluate()( 0,0 );
        auto dv_dx2 = integrate( elements( mesh ), dv_dx ).evaluate()( 0,0 );
        auto dv_dy2 = integrate( elements( mesh ), dv_dy ).evaluate()( 0,0 );
        BOOST_CHECK_SMALL( dv_dx1,  1e-12 );
        BOOST_CHECK_SMALL( dv_dy1, 1e-12 );
        BOOST_CHECK_SMALL( dv_dx2,  1e-12 );
        BOOST_CHECK_SMALL( dv_dy2, 1e-12 );


        auto m11= integrate( boundaryfaces( mesh ), gradv( u.template element<0>() )*N() ).evaluate();
        BOOST_TEST_MESSAGE( "int_bfaces(grad_u*N()) = " << m11 << "\n" );

        m11= integrate( boundaryfaces( mesh ), grad_exact*N() ).evaluate();
        BOOST_TEST_MESSAGE( "int_bfaces(grad_exact*N) = " << m11 << "\n" );

        auto m22= integrate( elements( mesh ), div_grad_exact ).evaluate();
        BOOST_TEST_MESSAGE( "int(div_grad_exact) = " << m22 << "\n" );

        auto m3= integrate( elements( mesh ), grad_exact  ).evaluate();
        BOOST_TEST_MESSAGE( "int(grad_exact) = " << m3 << "\n" );
        auto m4= integrate( elements( mesh ), gradv( u.template element<0>() ) - grad_exact ).evaluate();
        BOOST_TEST_MESSAGE( "int( grad(u)-grad_exact) = " << m4 << "\n" );


        auto m5= integrate( elements( mesh ),
                            trace( ( gradv( u.template element<0>() )-grad_exact )*trans( gradv( u.template element<0>() )-grad_exact ) )
                          ).evaluate();
        BOOST_TEST_MESSAGE( "|grad(u)-grad_exact|_0^2 = " << m5 << "\n" );

        auto F = M_backend->newVector( Xh );
        BOOST_TEST_MESSAGE( "[6] F built" );
        form1( _test=Xh, _vector=F ) = integrate( elements( mesh ),
                                                  ( trans( vec( constant( 1.0 ), constant( 1.0 ) ) ) * id( u.template element<0>() )+
                                                    id( u.template element<1>() ) ) );
        BOOST_TEST_MESSAGE( "[6] integration done" );
        F->printMatlab( "composite_F.m" );
        BOOST_TEST_MESSAGE( "[6] vector saved" );
        F->close();
        BOOST_TEST_MESSAGE( "[6] start tests" );
        value_type vf = inner_product( *F, u );
        value_type vv1 = integrate( elements( mesh ), idv( u.template element<0>() ) ).evaluate()( 0, 0 );
        value_type vv2 = integrate( elements( mesh ), idv( u.template element<0>() ) ).evaluate()( 1, 0 );
        value_type vv3 = integrate( elements( mesh ), idv( u.template element<1>() ) ).evaluate()( 0, 0 );


#if defined(USE_BOOST_TEST)
        BOOST_CHECK_CLOSE( vf,vv1+vv2+vv3, eps );
#endif
    }
    boost::shared_ptr<Feel::Backend<double> > M_backend;
    double meshSize;
    std::string shape;
    mesh_ptrtype mesh;
};

} // namespace Feel

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description integrationoptions( "Test Integration options" );
    integrationoptions.add_options()
    ( "hsize", Feel::po::value<double>()->default_value( 0.05 ), "h value" )
    ( "shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (hypercube, simplex, ellipsoid)" )
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
                           "Copyright (C) 2006,2007,2008,2009,2010 Universite Joseph Fourier (Grenoble I)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

#if defined(USE_BOOST_TEST)

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( integration )

#if 1
BOOST_AUTO_TEST_CASE( test_integration_1 )
{
    BOOST_TEST_MESSAGE( "Test integration Circle" );
    Feel::test_integration_circle<double> t;
#if defined( FEELPP_HAS_TBB )
    //int n = tbb::task_scheduler_init::default_num_threads();
    int n = 1 ;

    for ( int p=1; p<=n; ++p )
    {
        BOOST_TEST_MESSAGE( "[test_integration_1] start tests with nthreads = " << p );

        tbb::task_scheduler_init init( p );
        tbb::tick_count t0 = tbb::tick_count::now();

        t();

        tbb::tick_count t1 = tbb::tick_count::now();
        double t = ( t1-t0 ).seconds();

        BOOST_TEST_MESSAGE( "[test_integration_1] start tests with " << p << " threads, time=" << t << "seconds\n" );
    }

#else
    t();
#endif // FEELPP_HAS_TBB
    BOOST_TEST_MESSAGE( "Test integration Circle Done" );
}

BOOST_AUTO_TEST_CASE( test_integration_2 )
{
    BOOST_TEST_MESSAGE( "test_integration_2" );
    Feel::test_integration_domain<double> t;
    t();
}
BOOST_AUTO_TEST_CASE( test_integration_3 )
{
    BOOST_TEST_MESSAGE( "test_integration_3" );
    Feel::test_integration_boundary<double> t;
    t();
}
BOOST_AUTO_TEST_CASE( test_integration_4 )
{
    BOOST_TEST_MESSAGE( "test_integration_4" );
    Feel::test_integration_functions<2,double> t;
    t();
}
BOOST_AUTO_TEST_CASE( test_integration_5 )
{
    BOOST_TEST_MESSAGE( "test_integration_5" );
    Feel::test_integration_vectorial_functions<2,double> t;
    t();
}

BOOST_AUTO_TEST_CASE( test_integration_6 )
{
    BOOST_TEST_MESSAGE( "test_integration_6" );
    Feel::test_integration_composite_functions<2,double> t;
    t();
}

BOOST_AUTO_TEST_CASE( test_integration_7 )
{
    BOOST_TEST_MESSAGE( "test_integration_7" );
    Feel::test_integration_simplex<double> t;
    t();
}
BOOST_AUTO_TEST_CASE( test_integration_8 )
{
    BOOST_TEST_MESSAGE( "test_integration_8" );
    Feel::test_integration_matricial_functions<2,double> t;
    t();

}
#else
BOOST_AUTO_TEST_CASE( test_integration_3 )
{
    BOOST_TEST_MESSAGE( "test_integration_3" );
    Feel::test_integration_boundary<double> t;
    t();
}

#endif
BOOST_AUTO_TEST_SUITE_END()

#if 0
int BOOST_TEST_CALL_DECL
main( int argc, char* argv[] )
{
    Feel::Environment env( argc, argv );
    Feel::Assert::setLog( "test_integration.assert" );
    int ret = ::boost::unit_test::unit_test_main( &init_unit_test, argc, argv );

    return ret;
}
#endif

#else
int
main( int argc, char** argv )
{
    using namespace Feel;
    Feel::Environment env( _argc=argc, _argv=argv, _about=makeAbout(), _desc=makeOptions() );


    test_integration_circle<double> t1;
    t1();
#if 0
    test_integration_domain<double> t2( mpi.vm()["hsize"].as<double>() );
    t2();
    test_integration_boundary<double> t3( mpi.vm()["hsize"].as<double>() );
    t3();
    test_integration_functions<2,double> t4( mpi.vm()["hsize"].as<double>() );
    t4();
    test_integration_vectorial_functions<2,double> t5( mpi.vm()["hsize"].as<double>() );
    t5();
    test_integration_composite_functions<2,double> t6( mpi.vm()["hsize"].as<double>() );
    t6();
#endif
}
#endif // USE_BOOST_TEST
