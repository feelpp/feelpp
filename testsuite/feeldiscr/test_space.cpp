/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-11-08

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006-2010 Université Joseph Fourier (Grenoble I)

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
   \file test_space.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-11-08
 */
#define USE_BOOST_TEST 1

// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE function space testsuite
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN

#include <testsuite/testsuite.hpp>

#include <boost/timer.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelpoly/equispaced.hpp>
#include <feel/feelpoly/lagrange.hpp>
#include <feel/feelpoly/raviartthomas.hpp>
#include <feel/feelpoly/boundaryadaptedpolynomialset.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelvf/vf.hpp>

namespace Feel
{
template<int Dim, int N, typename T>
class TestSpace1
{
public:
    TestSpace1()
    {}

    void operator()() const
    {
        using namespace Feel;

        typedef Mesh<Simplex<Dim, 1> > mesh_type;
        typedef Simplex<Dim, 1> convex_type;

        typedef FunctionSpace<mesh_type, bases<Lagrange<N, Scalar> >, T> scalar_space_type;
        BOOST_STATIC_ASSERT( scalar_space_type::nDim == Dim );
        BOOST_STATIC_ASSERT( scalar_space_type::N_COMPONENTS == 1 );

        typedef FunctionSpace<mesh_type, bases<Lagrange<N, Vectorial> >, T> vectorial_space_type;
        BOOST_STATIC_ASSERT( vectorial_space_type::nDim == Dim );
        BOOST_STATIC_ASSERT( vectorial_space_type::N_COMPONENTS == Dim );

        typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
        mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                            _desc=domain( _name=( boost::format( "%1%-%2%" ) % "hypercube" % Dim ).str() ,
                                                    _usenames=true,
                                                    _shape="hypercube",
                                                    _dim=Dim,
                                                    _h=0.2 ) );

        boost::shared_ptr<scalar_space_type> Xh( new scalar_space_type( mesh ) );
        auto u = Xh->element( "u" );

        boost::shared_ptr<vectorial_space_type> Xhv( new vectorial_space_type( mesh ) );
        auto U = Xhv->element( "U" );


        const T eps = 1000*Feel::type_traits<T>::epsilon();

#if USE_BOOST_TEST
        BOOST_CHECK( Dim*u.nDof() == U.nDof() );

        for ( int i = 0; i < Dim; ++i )
        {
            BOOST_CHECK( U.comp( ( ComponentType )i ).is_scalar );
            U.comp( ( ComponentType )i ) = ublas::scalar_vector<T>( U.comp( ( ComponentType )i ).nDof(), i+1 );
            BOOST_CHECK( u.nDof() == U.comp( ( ComponentType )i ).nDof() );
            BOOST_CHECK( u.size() == U.comp( ( ComponentType )i ).size() );
            u = U.comp( ( ComponentType )i );
            BOOST_CHECK_SMALL( ublas::norm_inf( u.vec() - ublas::scalar_vector<T>( u.size(), i+1 ) ), eps );
        }

#endif

    }
};

template<int Dim, typename T>
class TestSpace2
{
public:
    typedef Mesh<Simplex<Dim, 1> > mesh_type;
    typedef bases<Lagrange<2, Vectorial>,
            Lagrange<1, Scalar>,
            Lagrange<1, Scalar> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type, T> space_type;

    TestSpace2()
    {}

    void operator()() const
    {
        using namespace Feel;

        typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
        mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                            _desc=domain( _name=( boost::format( "%1%-%2%" ) % "hypercube" % Dim ).str() ,
                                                    _usenames=true,
                                                    _shape="hypercube",
                                                    _dim=Dim,
                                                    _h=0.2 ) );


        boost::shared_ptr<space_type> Xh( space_type::New( mesh ) );
        auto U = Xh->element( "U" );

#if !defined( USE_BOOST_TEST )
#define BOOST_CHECK( e ) FEELPP_ASSERT( e ).warn( "BOOST_CHECK assertion failed" );
#endif

        BOOST_CHECK( Xh->nDof() == ( Xh->template functionSpace<0>()->nDof() +
                                     Xh->template functionSpace<1>()->nDof() +
                                     Xh->template functionSpace<2>()->nDof() ) );


        BOOST_CHECK( Xh->template functionSpace<1>()->nDof() == Xh->template functionSpace<2>()->nDof() );
        BOOST_CHECK( Xh->nDof() == U.nDof() );
        BOOST_CHECK( Xh->template functionSpace<0>()->nDof() == U.template functionSpace<0>()->nDof() );
        BOOST_CHECK( Xh->template functionSpace<1>()->nDof() == U.template functionSpace<1>()->nDof() );
        BOOST_CHECK( Xh->template functionSpace<2>()->nDof() == U.template functionSpace<2>()->nDof() );

        //
        // local(per processor) dof
        //
        BOOST_CHECK( Xh->nLocalDof() == ( Xh->template functionSpace<0>()->nLocalDof() +
                                          Xh->template functionSpace<1>()->nLocalDof() +
                                          Xh->template functionSpace<2>()->nLocalDof() ) );
        BOOST_CHECK( Xh->template functionSpace<1>()->nLocalDof() == Xh->template functionSpace<2>()->nLocalDof() );
        BOOST_CHECK( Xh->nLocalDof() == U.nLocalDof() );
        BOOST_CHECK( Xh->template functionSpace<0>()->nLocalDof() == U.template functionSpace<0>()->nLocalDof() );
        BOOST_CHECK( Xh->template functionSpace<1>()->nLocalDof() == U.template functionSpace<1>()->nLocalDof() );
        BOOST_CHECK( Xh->template functionSpace<2>()->nLocalDof() == U.template functionSpace<2>()->nLocalDof() );


        std::cout << "<0> start = " << U.template element<0>().start() << "\n"
                  << "<0> size = " << U.template element<0>().size() << "\n"
                  << "<0> ndof = " << U.template functionSpace<0>()->nDof() << "\n"
                  << "<1> start = " << U.template element<1>().start() << "\n"
                  << "<1> size = " << U.template element<1>().size() << "\n"
                  << "<1> ndof = " << U.template functionSpace<1>()->nDof() << "\n"
                  << "<2> start = " << U.template element<2>().start() << "\n"
                  << "<2> size = " << U.template element<2>().size() << "\n"
                  << "<2> ndof = " << U.template functionSpace<2>()->nDof() << "\n"
                  << "size = " << U.size() << "\n"
                  << "ndof = " << U.nDof() << "\n";
#if USE_BOOST_TEST
        // check elements
        BOOST_CHECK( U.template element<0>().start() == 0 &&
                     U.template element<0>().size() == U.template functionSpace<0>()->nDof() );
        BOOST_CHECK( U.template element<1>().start() == U.template functionSpace<0>()->nDof() &&
                     U.template element<1>().size() == U.template functionSpace<1>()->nDof() );
        BOOST_CHECK( U.template element<2>().start() == ( U.template functionSpace<0>()->nDof() + U.template functionSpace<1>()->nDof() ) &&
                     U.template element<2>().size() == U.template functionSpace<2>()->nDof() );
#else
        FEELPP_ASSERT( U.template element<0>().start() == 0 &&
                       U.template element<0>().size() == U.template functionSpace<0>()->nDof() )
        ( U.template element<0>().start() )
        ( U.template element<0>().size() )
        ( U.template functionSpace<0>()->nDof() ).warn( "invalid size" );
        FEELPP_ASSERT( U.template element<1>().start() == U.template functionSpace<0>()->nDof() &&
                       U.template element<1>().size() == U.template functionSpace<1>()->nDof() )
        ( U.template element<1>().start() )
        ( U.template functionSpace<0>()->nDof() )
        ( U.template element<1>().size() )
        ( U.template functionSpace<1>()->nDof() ).warn( "invalid size 2" );
        FEELPP_ASSERT( U.template element<2>().start() == ( U.template functionSpace<0>()->nDof() + U.template functionSpace<1>()->nDof() ) &&
                       U.template element<2>().size() == U.template functionSpace<2>()->nDof() )
        ( U.template element<2>().start() )
        ( U.template functionSpace<0>()->nDof() + U.template functionSpace<1>()->nDof() )
        ( U.template element<2>().size() )
        ( U.template functionSpace<2>()->nDof()  ).warn( "invalid size 3" );

#endif

#if 0

        if ( Dim == 2 )
        {
            auto ux = U.template element<0>().comp( X );
            auto uy = U.template element<0>().comp( Y );
            ux = ublas::scalar_vector<double>( ux.size(), -1 );
            //std::cout << "ux=" << ux << "\n";
            uy = ublas::scalar_vector<double>( uy.size(), -2 );
            //std::cout << "uy=" << uy << "\n";
            U.template element<1>() = ublas::scalar_vector<double>( U.template element<1>().size(), -3 );
            //std::cout << "U<1>=" << U.template element<1>() << "\n";
            U.template element<2>() = ublas::scalar_vector<double>( U.template element<2>().size(), -4 );
            //std::cout << "U<2>=" << U.template element<2>() << "\n";
            //std::cout << "U=" << U << "\n";
        }

#endif

    }
};

template<int Dim>
class TestSpaceRT
{
public:
    typedef Mesh<Simplex<Dim> > mesh_type;
    typedef bases<RaviartThomas<0> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;

    TestSpaceRT()
    {}

    void operator()() const
    {
        using namespace Feel;

        typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
        mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                            _desc=domain( _name=( boost::format( "%1%-%2%" ) % "hypercube" % Dim ).str() ,
                                                    _usenames=true,
                                                    _shape="hypercube",
                                                    _dim=Dim,
                                                    _h=0.2 ) );



        boost::shared_ptr<space_type> Xh( space_type::New( mesh ) );
        typename space_type::element_type U( Xh, "U" );

#if !defined( USE_BOOST_TEST )
#define BOOST_CHECK( e ) FEELPP_ASSERT( e ).warn( "BOOST_CHECK assertion failed" );
#define BOOST_CHECK_EQUAL( a, b ) FEELPP_ASSERT( a==b )(a)(b).warn( "BOOST_CHECK assertion failed" );
#endif
        BOOST_CHECK( Xh->is_scalar == false );
        BOOST_CHECK( Xh->is_vectorial == true );

        if ( Dim == 2 )
            BOOST_CHECK_EQUAL( Xh->nDof(), mesh->numEdges() );

        if ( Dim == 3 )
            BOOST_CHECK_EQUAL( Xh->nDof(), mesh->numFaces() );


    }
};


template<int Dim, int N, typename T>
class TestBASpace
{
public:
    TestBASpace()
    {
        // std::cout << "Testing " << Dim << " D Boundary adapted construction of order " << N << "\n";
    }

    void operator()() const
    {
        using namespace Feel;

        typedef Mesh<Simplex<Dim, 1> > mesh_type;
        typedef Simplex<Dim, 1> convex_type;

        typedef FunctionSpace<mesh_type, bases<BoundaryAdaptedPolynomialSet<N, Scalar> >, T> scalar_space_type;
        BOOST_STATIC_ASSERT( scalar_space_type::nDim == Dim );
        BOOST_STATIC_ASSERT( scalar_space_type::N_COMPONENTS == 1 );

        //std::cout << "antes vectorial\n";
        //        typedef FunctionSpace<mesh_type, bases<detail::BoundaryAdaptedPolynomialSet<Dim, N, Vectorial, T> >, T> vectorial_space_type;
        //	BOOST_STATIC_ASSERT( vectorial_space_type::nDim == Dim );
        //        BOOST_STATIC_ASSERT( vectorial_space_type::N_COMPONENTS == Dim );

        typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

        std::string fname;
#if 1
        std::cout << std::setfill( '=' ) << std::setw( 81 ) << "\n";
        std::cout << std::setfill( ' ' ) << "Testing " << Dim << " D Boundary adapted construction of order " << N << "\n";
#endif
        std::cout << " Generate Meshes !\n";

        mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                            _desc=domain( _name=( boost::format( "%1%-%2%" ) % "hypercube" % Dim ).str() ,
                                                    _usenames=true,
                                                    _shape="hypercube",
                                                    _dim=Dim,
                                                    _h=0.2 ) );


        std::cout << "Construct Space ...\n";

        boost::shared_ptr<scalar_space_type> Xh( new scalar_space_type( mesh ) );
        auto u = Xh->element( "u" );

        Xh.get()->dof().get()->showMe();

    }
};

template<int Dim, int N, typename T>
class TestSpaceProd
{
public:
    TestSpaceProd()
    {}

    void operator()() const
    {
        operator()( mpl::int_<Dim>() );
    }
    void operator()( mpl::int_<2> ) const
    {
        using namespace Feel;

        typedef Simplex<Dim, 1> convex1_type;
        typedef Mesh<Simplex<Dim, 1> > mesh1_type;
        typedef Simplex<Dim-1, 1,Dim> convex2_type;
        typedef Mesh<convex2_type> mesh2_type;

        typedef FunctionSpace<meshes<mesh1_type,mesh2_type>, bases<Lagrange<N, Scalar>,Lagrange<N,Scalar> > > space_type;
        auto mesh1 = createGMSHMesh( _mesh=new mesh1_type,
                                     _desc=domain( _name=( boost::format( "%1%-%2%" ) % "hypercube" % Dim ).str() ,
                                                   _usenames=true,
                                                   _shape="hypercube",
                                                   _dim=Dim,
                                                   _h=0.2 ) );
        auto mesh2 = mesh1->trace( boundaryfaces(mesh1) );

        auto m = fusion::make_vector(mesh1, mesh2 );
        BOOST_CHECK_EQUAL( fusion::at_c<0>(m), mesh1 );
        BOOST_CHECK_EQUAL( fusion::at_c<1>(m), mesh2 );
        auto Xh = space_type::New( _mesh=m );
        using namespace vf;
        auto res1 = integrate( elements(Xh->template mesh<0>()), cst(1.)).evaluate();
        BOOST_TEST_MESSAGE( "int 1 = " << res1 );
        BOOST_CHECK_CLOSE( res1(0,0), 1, 1e-12 );
        auto res2 = integrate( elements(Xh->template mesh<1>()), cst(1.)).evaluate();
        BOOST_TEST_MESSAGE( "int_bdy 1 = " << res2 );
        BOOST_CHECK_CLOSE( res2(0,0), 4, 1e-13 );

        auto U=Xh->element();
        if ( N == 1 )
        {
            auto u = U.template element<0>();
            auto l = U.template element<1>();
            BOOST_CHECK_EQUAL( u.size(), mesh1->numVertices() );
            BOOST_CHECK_EQUAL( l.size(), mesh2->numVertices() );
            u = vf::project( Xh->template functionSpace<0>(), elements(Xh->template mesh<0>() ), Px() );
            l = vf::project( Xh->template functionSpace<1>(), elements(Xh->template mesh<1>() ), Py() );
            auto int1_1 = integrate( elements(Xh->template mesh<0>()), Px()).evaluate()(0,0);
            auto int1_2 = integrate( elements(Xh->template mesh<0>()), idv(u)).evaluate()(0,0);
            BOOST_CHECK_CLOSE( int1_1, int1_2, 1e-13 );
            auto int2_1 = integrate( elements(Xh->template mesh<1>()), Py()).evaluate()(0,0);
            auto int2_2 = integrate( elements(Xh->template mesh<1>()), idv(l)).evaluate()(0,0);
            BOOST_CHECK_CLOSE( int1_1, int1_2, 1e-13 );
        }

        // Mh(domain), Lh(trace)
        //auto Xh = Mh*Lh;
        // typedef decltype( Mh*Lh ) space_type
        // auto Xh = FunctionSpace<decltype(meshes(Mh->mesh(),Lh->mesh())),
    }
    void operator()( mpl::int_<3> ) const
    {
        using namespace Feel;

        typedef Simplex<Dim, 1> convex1_type;
        typedef Mesh<Simplex<Dim, 1> > mesh1_type;
        typedef Simplex<Dim-1, 1,Dim> convex2_type;
        typedef Mesh<convex2_type> mesh2_type;
        typedef Simplex<Dim-2, 1,Dim> convex3_type;
        typedef Mesh<convex3_type> mesh3_type;

        typedef FunctionSpace<meshes<mesh1_type,mesh2_type,mesh3_type>, bases<Lagrange<N, Scalar>,Lagrange<N,Scalar>,Lagrange<N,Scalar> > > space_type;
        auto mesh1 = createGMSHMesh( _mesh=new mesh1_type,
                                     _desc=domain( _name=( boost::format( "%1%-%2%" ) % "hypercube" % Dim ).str() ,
                                                   _usenames=false,
                                                   _shape="hypercube",
                                                   _dim=Dim,
                                                   _h=0.2 ) );
        BOOST_TEST_MESSAGE( "N elements : " << mesh1->numElements() );
        //auto mesh2 = mesh1->trace( boundaryfaces(mesh1) );
        auto mesh2 = mesh1->trace( markedfaces(mesh1,6) );
        BOOST_TEST_MESSAGE( "N elements trace : " << mesh2->numElements() );
        auto mesh3 = mesh2->trace( boundaryfaces(mesh2) );
        BOOST_CHECK( !mesh3->elements().empty() );
        BOOST_TEST_MESSAGE( "N elements trace of trace: " << mesh3->numElements() );
        auto m = fusion::make_vector(mesh1, mesh2, mesh3 );
        BOOST_CHECK_EQUAL( fusion::at_c<0>(m), mesh1 );
        BOOST_CHECK_EQUAL( fusion::at_c<1>(m), mesh2 );
        BOOST_CHECK_EQUAL( fusion::at_c<2>(m), mesh3 );
        auto Xh = space_type::New( _mesh=m );
        using namespace vf;
        auto res1 = integrate( elements(Xh->template mesh<0>()), cst(1.)).evaluate();
        BOOST_TEST_MESSAGE( "int 1 = " << res1 );
        BOOST_CHECK_CLOSE( res1(0,0), 1, 5e-13 );
        auto res2 = integrate( elements(Xh->template mesh<1>()), cst(1.)).evaluate();
        BOOST_TEST_MESSAGE( "int_trace 1 = " << res2 );
        BOOST_CHECK_CLOSE( res2(0,0), 1, 5e-13 );
        auto res3 = integrate( elements(Xh->template mesh<2>()), cst(1.)).evaluate();
        BOOST_TEST_MESSAGE( "int_trace_trace 1 = " << res3 );
        BOOST_CHECK_CLOSE( res3(0,0), 4, 5e-13 );

        // Mh(domain), Lh(trace)
        //auto Xh = Mh*Lh;
        // typedef decltype( Mh*Lh ) space_type
        // auto Xh = FunctionSpace<decltype(meshes(Mh->mesh(),Lh->mesh())),
    }
};

}

#if USE_BOOST_TEST

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( space )

BOOST_AUTO_TEST_CASE( test_space1_11 )
{
    BOOST_TEST_MESSAGE( "test_space1_11" );
    Feel::TestSpace1<1, 1, double> t;
    t();
    BOOST_TEST_MESSAGE( "test_space1_11 done" );
}
BOOST_AUTO_TEST_CASE( test_space1_12 )
{
    BOOST_TEST_MESSAGE( "test_space1_12" );
    Feel::TestSpace1<1, 2, double> t;
    t();
    BOOST_TEST_MESSAGE( "test_space1_12 done" );
}
BOOST_AUTO_TEST_CASE( test_space1_21 )
{
    BOOST_TEST_MESSAGE( "test_space1_21" );
    Feel::TestSpace1<2, 1, double> t;
    t();
    BOOST_TEST_MESSAGE( "test_space1_21 done" );
}
BOOST_AUTO_TEST_CASE( test_space1_22 )
{
    BOOST_TEST_MESSAGE( "test_space1_22" );
    Feel::TestSpace1<2, 2, double> t;
    t();
    BOOST_TEST_MESSAGE( "test_space1_22 done" );
}
BOOST_AUTO_TEST_CASE( test_space1_32 )
{
    BOOST_TEST_MESSAGE( "test_space1_32" );
    Feel::TestSpace1<3, 2, double> t;
    t();
    BOOST_TEST_MESSAGE( "test_space1_32 done" );
}
BOOST_AUTO_TEST_CASE( test_space2_2 )
{
    BOOST_TEST_MESSAGE( "test_space2_2" );
    Feel::TestSpace2<2, double> t;
    t();
    BOOST_TEST_MESSAGE( "test_space2_2 done" );
}
BOOST_AUTO_TEST_CASE( test_space2_3 )
{
    BOOST_TEST_MESSAGE( "test_space2_3" );
    Feel::TestSpace2<3, double> t;
    t();
    BOOST_TEST_MESSAGE( "test_space2_3 done" );
}
BOOST_AUTO_TEST_CASE( test_spaceprod2_1 )
{
    BOOST_TEST_MESSAGE( "test_spaceprod2_1" );
    Feel::TestSpaceProd<2, 1, double> t;
    t();
    BOOST_TEST_MESSAGE( "test_spaceprod2_1 done" );
}
BOOST_AUTO_TEST_CASE( test_spaceprod2_3 )
{
    BOOST_TEST_MESSAGE( "test_spaceprod2_3" );
    Feel::TestSpaceProd<2, 3, double> t;
    t();
    BOOST_TEST_MESSAGE( "test_spaceprod2_3 done" );
}

BOOST_AUTO_TEST_CASE( test_spaceprod3_1 )
{
    BOOST_TEST_MESSAGE( "test_spaceprod3_1" );
    Feel::TestSpaceProd<3, 1, double> t;
    t();
    BOOST_TEST_MESSAGE( "test_spaceprod3_1 done" );
}

//BOOST_AUTO_TEST_CASE( test_space_rt_1 ) { BOOST_TEST_MESSAGE( "test_space_rt_1" );   Feel::TestSpaceRT<2> t;t(); BOOST_TEST_MESSAGE( "test_space_rt_1 done" );}
//BOOST_AUTO_TEST_CASE( test_space_rt_2 ) { BOOST_TEST_MESSAGE( "test_space_rt_2" );   Feel::TestSpaceRT<3> t;t(); BOOST_TEST_MESSAGE( "test_space_rt_2 done" );}

BOOST_AUTO_TEST_SUITE_END()
#if 0
int BOOST_TEST_CALL_DECL
main( int argc, char* argv[] )
{
    Feel::Environment env( argc, argv );
    Feel::Assert::setLog( "test_space.assert" );
    int ret = ::boost::unit_test::unit_test_main( &init_unit_test, argc, argv );

    return ret;
}
#endif

#else
int main( int argc, char** argv )
{
    Feel::Environment env( argc, argv );

    Feel::Assert::setLog( "test_space.assert" );
#if 0
    Feel::TestSpace1<1, 2, double> t11;
    t11();
    Feel::TestSpace1<2, 2, double> t12;
    t12();
    Feel::TestSpace1<3, 2, double> t13;
    t13();
    Feel::TestSpace2<2, double> t21;
    t21();
    Feel::TestSpace2<3, double> t22;
    t22();
#endif
    Feel::TestSpaceRT<3> trt;
    trt();

}
#endif

