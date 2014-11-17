/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-10-21

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
   \file test_disc.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-10-21
 */
#define USE_BOOST_TEST 1
// Boost.Test

// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE discontinuity testsuite
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN


#include <testsuite/testsuite.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/options.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelmesh/geoentity.hpp>
#include <feel/feelmesh/refentity.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>

#include <feel/feelvf/vf.hpp>

const double DEFAULT_MESH_SIZE=0.1;

namespace Feel
{
template<typename T, int Dim, int Order = 1>
struct imesh
{
    typedef Simplex<Dim, Order> convex_type;
    typedef Mesh<convex_type, T > type;
    typedef boost::shared_ptr<type> ptrtype;
};


template<typename value_type = double, int Dim=2, int Order =2>
struct test_disc: public Application
{
    typedef typename imesh<value_type,Dim>::convex_type convex_type;
    typedef typename imesh<value_type,Dim>::type mesh_type;
    typedef typename imesh<value_type,Dim>::ptrtype mesh_ptrtype;
    typedef DiscontinuousInterfaces<fusion::vector<mpl::vector<mpl::int_<4>, mpl::int_<6>, mpl::int_<7> > > > discontinuity_type;
    typedef bases<Lagrange<Order, Scalar, discontinuity_type> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    test_disc()
        :
        Application(),
        backend( Backend<double>::build( soption( _name="backend" ) ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( this->vm()["shape"].template as<std::string>() ),
        mesh()
    {
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _desc=createMesh(),
                               _update=MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES );
        Xh = space_type::New( mesh );
    }
    gmsh_ptrtype createMesh()
    {
        std::ostringstream ostr;
        ostr << "Mesh.MshFileVersion = " << 2 << ";\n"
             << "a=" << -1 << ";\n"
             << "b=" << 1 << ";\n"
             << "c=" << -1 << ";\n"
             << "d=" << 1 << ";\n"
             << "h=" << meshSize << ";\n"
             << "Point(1) = {a,c,0.0,h};\n"
             << "Point(2) = {b,c,0.0,h};\n"
             << "Point(3) = {b,d,0.0,h};\n"
             << "Point(4) = {a,d,0.0,h};\n"
             << "Point(5) = {0,c,0.0,h};\n"
             << "Point(6) = {0,d,0.0,h};\n"
             << "Point(7) = {a,0,0.0,h};\n"
             << "Point(8) = {b,0,0.0,h};\n"
             << "Point(9) = {0,0,0.0,h};\n"
             << "Line(1) = {1,5};\n"
             << "Line(2) = {5,2};\n"
             << "Line(3) = {2,8};\n"
             << "Line(4) = {8,3};\n"
             << "Line(5) = {3,6};\n"
             << "Line(6) = {6,4};\n"
             << "Line(7) = {4,7};\n"
             << "Line(8) = {7,1};\n"
             << "/* discontinuity (vertical line) */\n"
             << "Line(9) = {5, 9};\n"
             << "Line(10) = {9, 6};\n"
             << "/* horizontal line through square */\n"
             << "Line(11) = {7, 9};\n"
             << "Line(12) = {9, 8};\n"
             << "\n"
             << "Line Loop(19) = {3, -12, -9, 2};\n"
             << "Plane Surface(20) = {19};\n"
             << "Line Loop(21) = {4, 5, -10, 12};\n"
             << "Plane Surface(22) = {21};\n"
             << "Line Loop(23) = {6, 7, 11, 10};\n"
             << "Plane Surface(24) = {23};\n"
             << "Line Loop(25) = {8, 1, 9, -11};\n"
             << "Plane Surface(26) = {25};\n"
             << "\n"
             << "Physical Line(\"Tflux\") = {3, 4};\n"
             << "Physical Line(\"Tfixed\") = {8, 7};\n"
             << "Physical Line(\"Tinsulated\") = {1, 2, 6, 5};\n"
             << "Physical Line(\"Tdiscontinuity\") = {10, 9};\n"
             << "Physical Line(\"Tline\") = {11, 12};\n"
             << "\n"
             << "Physical Surface(\"k1\") = {20, 22};\n"
             << "Physical Surface(\"k2\") = {24, 26};\n";

        gmsh_ptrtype gmshp( new gmsh_type );
        gmshp->setPrefix( "hypercube-2" );
        gmshp->setDescription( ostr.str() );

        return gmshp;
    }

    void operator()()
    {
        using namespace Feel::vf;
        auto u = Xh->element();

        // this is the way to take into account the discontinuity at x = 0
        // properly
        u = vf::project( Xh, elements( mesh ),
                         ( emarker()==mesh->markerName( "k2" ) )*( 2-Px() )+
                         ( emarker()==mesh->markerName( "k1" ) )*Px() );
        BOOST_TEST_MESSAGE( "here\n" );
        double len1 = integrate( markedfaces( mesh, "Tdiscontinuity" ), cst( 1.0 ) ).evaluate()( 0,0 );
        BOOST_CHECK_CLOSE( len1, 2, 1e-12 );
        BOOST_TEST_MESSAGE( "here 1\n" );
        double len11 = integrate( markedfaces( mesh, "Tdiscontinuity" ), leftfacev( cst( 1.0 ) )+rightfacev( cst( 1. ) ) ).evaluate()( 0,0 );
        BOOST_CHECK_CLOSE( len11, 4, 1e-12 );
        BOOST_TEST_MESSAGE( "here 2\n" );
        double len2 = integrate( markedfaces( mesh, "Tline" ), cst( 1.0 ) ).evaluate()( 0,0 );
        BOOST_CHECK_CLOSE( len2, 2, 1e-12 );
        auto int1 = integrate( markedfaces( mesh, "Tdiscontinuity" ), jumpv( idv( u ) ) ).evaluate();
        BOOST_CHECK_CLOSE( int1( 0,0 ), 4, 1e-12 );
        BOOST_CHECK_SMALL( int1( 1,0 ), 1e-12 );
        auto int2 = integrate( markedfaces( mesh, "Tdiscontinuity" ), leftfacev( idv( u ) )+rightfacev( idv( u ) ) ).evaluate();
        BOOST_CHECK_CLOSE( int2( 0,0 ), 4, 1e-12 );

        u = vf::project( Xh, elements( mesh ),
                         ( emarker()==mesh->markerName( "k2" ) )*( 2-Px() )*Py()-
                         ( emarker()==mesh->markerName( "k1" ) )*Px()*Py() );
        auto int3 = integrate( markedfaces( mesh, "Tdiscontinuity" ), jumpv( idv( u ) ) ).evaluate();
        BOOST_CHECK_SMALL( int3( 0,0 ), 1e-12 );
        BOOST_CHECK_SMALL( int3( 1,0 ), 1e-12 );

        u = vf::project( Xh, elements( mesh ), sin( Px() ) );
        auto int4 = integrate( markedfaces( mesh, "Tdiscontinuity" ), jumpv( idv( u ) ) ).evaluate();
        BOOST_CHECK_SMALL( int3( 0,0 ), 1e-12 );
        BOOST_CHECK_SMALL( int3( 1,0 ), 1e-12 );
    }
    boost::shared_ptr<Feel::Backend<double> > backend;
    double meshSize;
    std::string shape;
    mesh_ptrtype mesh;
    space_ptrtype Xh;
};

} // Feel
inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description integrationoptions( "Test Function Space features options" );
    integrationoptions.add_options()
    ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "h value" )
    ( "shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (hypercube, simplex, ellipsoid)" )
    ;
    return integrationoptions.add( Feel::feel_options() );
}

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "test_disc" ,
                           "test_disc" ,
                           "0.2",
                           "1D/2D/3D functionspace features checks",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2010 Université Joseph Fourier (Grenoble I)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );

BOOST_AUTO_TEST_SUITE( disc )

//typedef boost::mpl::list<boost::mpl::int_<1>,boost::mpl::int_<2>,boost::mpl::int_<3> > dim_types;
//typedef boost::mpl::list<boost::mpl::int_<1>,boost::mpl::int_<2> > dim_types;
typedef boost::mpl::list<boost::mpl::int_<2> > dim_types;
//typedef boost::mpl::list<boost::mpl::int_<2>,boost::mpl::int_<3>,boost::mpl::int_<1> > dim_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( test_disc, T, dim_types )
{
    BOOST_TEST_MESSAGE( "Test disc (" << T::value << "D)" );
    Feel::test_disc<double,T::value> t;
    t();
    BOOST_TEST_MESSAGE( "Test disc (" << T::value << "D) done." );
}

BOOST_AUTO_TEST_SUITE_END()



#if 0
int BOOST_TEST_CALL_DECL
main( int argc, char* argv[] )
{
    Feel::Environment env( argc, argv );
    Feel::Assert::setLog( "test_disc.assert" );
    int ret = ::boost::unit_test::unit_test_main( &init_unit_test, argc, argv );

    return ret;
}


#endif
