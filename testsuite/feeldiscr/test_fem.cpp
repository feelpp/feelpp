/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-08-18

  Copyright (C) 2005,2006 EPFL

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
   \file test_fem.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-08-18
 */
// Boost.Test
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/rational.hpp>

#include <feel/feelpoly/lagrange.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feelalg/matrixublas.hpp>
#include <feel/feelalg/vectorublas.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/gmsh.hpp>


using namespace Feel;

const int Dim = 2;
const int Order = 1;

template<typename FE>
struct
        test_lagrange
{
    void operator()() const
    {

        typedef typename FE::value_type value_type;

        Gmsh __gmsh;
        std::string fname;

        if ( Dim == 2 )
        {
            std::ostringstream ostr;
            ostr << "h=" << 0.1 << ";\n"
                 << "Point(1) = {0.0,0.0,0.0,h};\n"
                 << "Point(2) = {0.0,1,0.0,h};\n"
                 << "Point(3) = {1,1,0.0,h};\n"
                 << "Point(4) = {1,0,0.0,h};\n"
                 << "Line(1) = {1,4};\n"
                 << "Line(2) = {4,3};\n"
                 << "Line(3) = {3,2};\n"
                 << "Line(4) = {2,1};\n"
                 << "Line Loop(5) = {3,4,1,2};\n"
                 << "Plane Surface(6) = {5};\n"
                 << "Physical Line(10) = {1,3};\n"
                 << "Physical Line(20) = {2,4};\n"
                 << "Physical Surface(30) = {6};\n";

            fname = __gmsh.generate( "square", ostr.str() );
        }

        else
        {
            std::ostringstream ostr;
            ostr << "h=" << 0.1 << ";\n"
                 << "Point(1) = {0.0,0.0,0.0,h};\n"
                 << "Point(2) = {0.0,1,0.0,h};\n"
                 << "Point(3) = {1,1,0.0,h};\n"
                 << "Point(4) = {1,0,0.0,h};\n"
                 << "Line(1) = {1,4};\n"
                 << "Line(2) = {4,3};\n"
                 << "Line(3) = {3,2};\n"
                 << "Line(4) = {2,1};\n"
                 << "Line Loop(5) = {3,4,1,2};\n"
                 << "Plane Surface(6) = {5};\n"
                 << "Extrude Surface {6, {0.0,0.0,1.0}};\n"
                 << "Physical Surface(10) = {15,23,6,28};\n"
                 << "Physical Surface(20) = {19,27};\n"
                 << "Surface Loop(31) = {28,15,-6,19,23,27};\n"
                 << "\n"
                 << "// volume\n"
                 << "Volume(1) = {31};\n"
                 << "Physical Volume(2) = {1};\n";

            fname = __gmsh.generate( "cube", ostr.str() );
        }

        /*mesh*/
        typedef Mesh<GeoEntity<Simplex<Dim, 1> > > mesh_type;
        typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

        mesh_ptrtype mesh( new mesh_type );

        ImporterGmsh<mesh_type> import( fname );
        mesh->accept( import );

        typedef FunctionSpace<mesh_type, FE, value_type > space_type;
        boost::shared_ptr<space_type> Xh( new space_type( mesh ) );

        typedef MatrixUBlas<value_type,ublas::column_major> csc_matrix_type;
        csc_matrix_type M;
        vf::BilinearForm<space_type, space_type, csc_matrix_type> a( Xh, Xh, M );

        VectorUblas<double> F( Xh->nDof() );
        vf::LinearForm<space_type, VectorUblas<double> > l( Xh, F );
        using namespace Feel::vf;
        typename space_type::element_type u( Xh.get() );
        typename space_type::element_type v( Xh.get() );
        a = integrate( elements( *mesh ), IM_PK<2, 2>(), id( v )*idt( u ) )+
            integrate( markedfaces( *mesh, 20 ), IM_PK<2, 2>(), 2.2*id( v )*idt( u ) );
        l = integrate( elements( *mesh ), IM_PK<2, 2>(), id( v ) )+
            integrate( markedfaces( *mesh, 10 ), IM_PK<2, 2>(), 3.3*id( v ) );

        ublas::vector<value_type> one( ublas::scalar_vector<value_type>( M.size2(), 1.0 ) );
        std::cout << "1^t M 1 = " << ublas::inner_prod( one, ublas::prod( M.mat(), one ) ) << "\n";
        std::cout << "1^t L = " << ublas::inner_prod( one, F ) << "\n";

    }

};

#if 1
test_suite*
init_unit_test_suite( int /*argc*/, char** /*argv*/ )
{
    test_suite* test = BOOST_TEST_SUITE( "Finite element test suite" );


    test_lagrange<fem::Lagrange<2, 2, Scalar, Continuous> > a;
    a();

    // 1D simplex
    //     test->add( BOOST_TEST_CASE( ( test_lagrange<fem::PK<1, 1, real96_type> >() ) ) );
    //     test->add( BOOST_TEST_CASE( ( test_lagrange<fem::PK<1, 2, real96_type> >() ) ) );
    //     test->add( BOOST_TEST_CASE( ( test_lagrange<fem::PK<1, 3, real96_type> >() ) ) );
    //     test->add( BOOST_TEST_CASE( ( test_lagrange<fem::PK<1, 4, real96_type> >() ) ) );
    //     test->add( BOOST_TEST_CASE( ( test_lagrange<fem::PK<1, 5, real96_type> >() ) ) );
    //
    //     // 2D simplex
    //     test->add( BOOST_TEST_CASE( ( test_lagrange<fem::PK<2, 1, real96_type> >() ) ) );
    //     test->add( BOOST_TEST_CASE( ( test_lagrange<fem::PK<2, 2, real96_type> >() ) ) );
    //test->add( BOOST_TEST_CASE( (test_lagrange<fem::Lagrange<2, 1, real96_type> >())  ) );
    //     test->add( BOOST_TEST_CASE( ( test_lagrange<fem::PK<2, 4, real96_type> >() ) ) );
    //     test->add( BOOST_TEST_CASE( ( test_lagrange<fem::PK<2, 5, real96_type> >() ) ) );
    //
    //     // 3D simplex
    //     test->add( BOOST_TEST_CASE( ( test_lagrange<fem::PK<3, 1, real96_type> >() ) ) );
    //     test->add( BOOST_TEST_CASE( ( test_lagrange<fem::PK<3, 2, real96_type> >() ) ) );
    //     test->add( BOOST_TEST_CASE( ( test_lagrange<fem::PK<3, 3, real96_type> >() ) ) );
    //     test->add( BOOST_TEST_CASE( ( test_lagrange<fem::PK<3, 4, real96_type> >() ) ) );
    //     test->add( BOOST_TEST_CASE( ( test_lagrange<fem::PK<3, 5, real96_type> >() ) ) );

    //     test->add( BOOST_TEST_CASE( ( test_lagrange<fem::P3p<2, dd_real,  Simplex> >() ) ) );



    return test;
}
#else
int main( int /*argc*/, char** /*argv*/ )
{
    test_lagrange<fem::Lagrange<2, 1, Scalar, Continuous, double> > a;
    a();

}
#endif

