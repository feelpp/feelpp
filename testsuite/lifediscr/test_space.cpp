/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-11-08

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file test_space.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-11-08
 */
//#define USE_BOOST_TEST 1
#define BOOST_TEST_MAIN
#include <boost/timer.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <life/lifepoly/equispaced.hpp>
#include <life/lifepoly/lagrange.hpp>

#include <life/lifediscr/mesh.hpp>
#include <life/lifediscr/functionspace.hpp>

#include <life/lifefilters/importergmsh.hpp>
#include <life/lifefilters/gmsh.hpp>

#include <testsuite/testsuite.hpp>


namespace Life
{
template<int Dim, int N, typename T>
class TestSpace1
{
public:
    TestSpace1()
    {}

    void operator()() const
    {
        using namespace Life;

        typedef Mesh<Simplex<Dim, 1> > mesh_type;
        typedef Simplex<Dim, 1> convex_type;

        typedef FunctionSpace<mesh_type, fusion::vector<Lagrange<N, Scalar> >, T> scalar_space_type;
        BOOST_STATIC_ASSERT( scalar_space_type::nDim == Dim );
        BOOST_STATIC_ASSERT( scalar_space_type::N_COMPONENTS == 1 );

        typedef FunctionSpace<mesh_type, fusion::vector<Lagrange<N, Vectorial> >, T> vectorial_space_type;
        BOOST_STATIC_ASSERT( vectorial_space_type::nDim == Dim );
        BOOST_STATIC_ASSERT( vectorial_space_type::N_COMPONENTS == Dim );

        boost::shared_ptr<mesh_type> mesh( new mesh_type );

        std::string fname;

        if ( Dim == 1 )
        {
            Gmsh gen;
            fname =  gen.generateLine("line", 0.2);
        }

        if ( Dim == 2 )
        {
            Gmsh gen;
            fname =  gen.generateSquare("square", 0.2);
        }

        if ( Dim == 3 )
        {
            Gmsh gen;
            fname = gen.generateCube("cube", 0.2);
        }

        ImporterGmsh<mesh_type> import( fname );
        mesh->accept( import );
        mesh->components().set( MESH_CHECK | MESH_RENUMBER | MESH_UPDATE_EDGES | MESH_UPDATE_FACES );
        mesh->updateForUse();

        boost::shared_ptr<scalar_space_type> Xh( new scalar_space_type( mesh ) );
        typename scalar_space_type::element_type u( Xh, "u" );

        boost::shared_ptr<vectorial_space_type> Xhv( new vectorial_space_type( mesh ) );
        typename vectorial_space_type::element_type U( Xhv, "U" );

        const T eps = 1000*Life::type_traits<T>::epsilon();

#if USE_BOOST_TEST
        BOOST_CHECK( Dim*u.nDof() == U.nDof() );
        for ( int i = 0; i < Dim; ++i )
            {
                BOOST_CHECK( U.comp((ComponentType)i).is_scalar );
                U.comp((ComponentType)i) = ublas::scalar_vector<T>( U.comp((ComponentType)i).nDof(), i+1 );
                BOOST_CHECK( u.nDof() == U.comp((ComponentType)i).nDof() );
                BOOST_CHECK( u.size() == U.comp((ComponentType)i).size() );
                u = U.comp((ComponentType)i);
                BOOST_CHECK_SMALL( ublas::norm_inf( u - ublas::scalar_vector<T>( u.size(), i+1 ) ), eps );
            }
#endif

    }
};

template<int Dim, typename T>
class TestSpace2
{
public:
    typedef Mesh<Simplex<Dim, 1> > mesh_type;
    typedef fusion::vector<Lagrange<2, Vectorial>,
                           Lagrange<1, Scalar>,
                           Lagrange<1, Scalar> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type, T> space_type;

    TestSpace2()
    {}

    void operator()() const
    {
        using namespace Life;

        boost::shared_ptr<mesh_type> mesh( new mesh_type );

        std::string fname;
        if ( Dim == 2 )
        {
            Gmsh gen;
            fname =  gen.generateSquare("square", 0.2);
        }
        if ( Dim == 3 )
        {
            Gmsh gen;
            fname = gen.generateCube("cube", 0.2);
        }
        ImporterGmsh<mesh_type> import( fname );
        mesh->accept( import );
        mesh->components().set( MESH_CHECK | MESH_RENUMBER | MESH_UPDATE_EDGES | MESH_UPDATE_FACES );
        mesh->updateForUse();

        boost::shared_ptr<space_type> Xh( space_type::New( mesh ) );
        typename space_type::element_type U( Xh, "U" );

#if !defined( USE_BOOST_TEST )
#define BOOST_CHECK( e ) LIFE_ASSERT( e ).warn( "BOOST_CHECK assertion failed" );
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
        BOOST_CHECK( U.template element<2>().start() == (U.template functionSpace<0>()->nDof() + U.template functionSpace<1>()->nDof()) &&
                     U.template element<2>().size() == U.template functionSpace<2>()->nDof() );
#else
        LIFE_ASSERT( U.template element<0>().start() == 0 &&
                     U.template element<0>().size() == U.template functionSpace<0>()->nDof() )
            ( U.template element<0>().start() )
            ( U.template element<0>().size() )
            ( U.template functionSpace<0>()->nDof() ).warn( "invalid size" );
        LIFE_ASSERT( U.template element<1>().start() == U.template functionSpace<0>()->nDof() &&
                     U.template element<1>().size() == U.template functionSpace<1>()->nDof() )
            ( U.template element<1>().start() )
            ( U.template functionSpace<0>()->nDof() )
            ( U.template element<1>().size() )
            ( U.template functionSpace<1>()->nDof() ).warn( "invalid size 2");
        LIFE_ASSERT( U.template element<2>().start() == (U.template functionSpace<0>()->nDof() + U.template functionSpace<1>()->nDof()) &&
                     U.template element<2>().size() == U.template functionSpace<2>()->nDof() )
            ( U.template element<2>().start() )
            ( U.template functionSpace<0>()->nDof() + U.template functionSpace<1>()->nDof() )
            ( U.template element<2>().size() )
            ( U.template functionSpace<2>()->nDof()  ).warn( "invalid size 3" );

#endif

#if 0
        if ( Dim == 2 )
            {
                typename space_type::element_type::template sub_element<0>::type::component_type ux = U.template element<0>().comp(X);
                typename space_type::element_type::template sub_element<0>::type::component_type uy = U.template element<0>().comp(Y);
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
        using namespace Life;

        typedef Mesh<Simplex<Dim, 1> > mesh_type;
        typedef Simplex<Dim, 1> convex_type;

        typedef FunctionSpace<mesh_type, fusion::vector<BoundaryAdaptedPolynomialSet<Dim, N, Scalar, T> >, T> scalar_space_type;
        BOOST_STATIC_ASSERT( scalar_space_type::nDim == Dim );
        BOOST_STATIC_ASSERT( scalar_space_type::N_COMPONENTS == 1 );

        //std::cout << "antes vectorial\n";
	//        typedef FunctionSpace<mesh_type, fusion::vector<BoundaryAdaptedPolynomialSet<Dim, N, Vectorial, T> >, T> vectorial_space_type;
	//	BOOST_STATIC_ASSERT( vectorial_space_type::nDim == Dim );
	//        BOOST_STATIC_ASSERT( vectorial_space_type::N_COMPONENTS == Dim );

        boost::shared_ptr<mesh_type> mesh( new mesh_type );

        std::string fname;
#if 1
	std::cout << std::setfill('=') << std::setw(81) << "\n";
	std::cout << std::setfill(' ') << "Testing " << Dim << " D Boundary adapted construction of order " << N << "\n";
#endif
	std::cout << " Generate Meshes !\n";

        if ( Dim == 1 )
        {
            Gmsh gen;
            fname =  gen.generateLine("line", 0.3);
        }

        if ( Dim == 2 )
        {
            Gmsh gen;
            fname =  gen.generateSquare("square", 0.3);
        }

        if ( Dim == 3 )
        {
            Gmsh gen;
            fname = gen.generateCube("cube", 0.8);
        }

        ImporterGmsh<mesh_type> import( fname );
        mesh->accept( import );
        mesh->components().set( MESH_CHECK | MESH_RENUMBER | MESH_UPDATE_EDGES | MESH_UPDATE_FACES );
        mesh->updateForUse();

	std::cout << "Construct FESpace ...\n";

        boost::shared_ptr<scalar_space_type> Xh( new scalar_space_type( mesh ) );
	typename scalar_space_type::element_type u( Xh, "u" );

	Xh.get()->dof().get()->showMe();

    }
};
}

#if USE_BOOST_TEST
boost::unit_test::test_suite*
init_unit_test_suite( int /*argc*/, char* /*argv*/[] )
// main( int argc, char* argv[] )
{
    boost::unit_test::test_suite* test = BOOST_TEST_SUITE( "FunctionSpace test suite" );

#if 1 // Finite Element

#if 0
    test->add( BOOST_TEST_CASE( ( Life::TestSpace1<1, 1, double>() )  ) );
    test->add( BOOST_TEST_CASE( ( Life::TestSpace1<1, 2, double>() )  ) );
    test->add( BOOST_TEST_CASE( ( Life::TestSpace1<1, 3, double>() )  ) );
    test->add( BOOST_TEST_CASE( ( Life::TestSpace1<1, 4, double>() )  ) );

    test->add( BOOST_TEST_CASE( ( Life::TestSpace1<2, 1, double>() )  ) );
    test->add( BOOST_TEST_CASE( ( Life::TestSpace1<2, 2, double>() )  ) );
    test->add( BOOST_TEST_CASE( ( Life::TestSpace1<2, 3, double>() )  ) );
    test->add( BOOST_TEST_CASE( ( Life::TestSpace1<2, 4, double>() )  ) );

    test->add( BOOST_TEST_CASE( ( Life::TestSpace1<3, 1, double>() )  ) );
    test->add( BOOST_TEST_CASE( ( Life::TestSpace1<3, 2, double>() )  ) );
    test->add( BOOST_TEST_CASE( ( Life::TestSpace1<3, 3, double>() )  ) );
    test->add( BOOST_TEST_CASE( ( Life::TestSpace1<3, 4, double>() )  ) );
#endif
    test->add( BOOST_TEST_CASE( ( Life::TestSpace2<2, double>() )  ) );
    test->add( BOOST_TEST_CASE( ( Life::TestSpace2<3, double>() )  ) );

#else // Boundary Adapted

#if 0 // 1D tests
    test->add( BOOST_TEST_CASE( ( Life::TestBASpace<1, 1, double>() )  ) );
    test->add( BOOST_TEST_CASE( ( Life::TestBASpace<1, 2, double>() )  ) );
    test->add( BOOST_TEST_CASE( ( Life::TestBASpace<1, 3, double>() )  ) );
    test->add( BOOST_TEST_CASE( ( Life::TestBASpace<1, 6, double>() )  ) );
#endif

#if 0 // 2D tests
    test->add( BOOST_TEST_CASE( ( Life::TestBASpace<2, 1, double>() )  ) );
    test->add( BOOST_TEST_CASE( ( Life::TestBASpace<2, 2, double>() )  ) );
    test->add( BOOST_TEST_CASE( ( Life::TestBASpace<2, 3, double>() )  ) );
    test->add( BOOST_TEST_CASE( ( Life::TestBASpace<2, 6, double>() )  ) );
#endif

#if 1 // 3D tests
    test->add( BOOST_TEST_CASE( ( Life::TestBASpace<3, 1, double>() )  ) );
    test->add( BOOST_TEST_CASE( ( Life::TestBASpace<3, 2, double>() )  ) );
    test->add( BOOST_TEST_CASE( ( Life::TestBASpace<3, 3, double>() )  ) );
    test->add( BOOST_TEST_CASE( ( Life::TestBASpace<3, 6, double>() )  ) );
#endif

#endif // Finite Element vs Boundary Adapted

    //boost::unit_test::framework::run( test );
    return test;
}
#else
int main()
{
    Life::Assert::setLog( "test_space.assert");
    Life::TestSpace1<1, 2, double> t11;t11();
    Life::TestSpace1<2, 2, double> t12;t12();
    Life::TestSpace1<3, 2, double> t13;t13();
    Life::TestSpace2<2, double> t21;t21();
    Life::TestSpace2<3, double> t22;t22();

}
#endif

