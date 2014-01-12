/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-08-13

  Copyright (C) 2008,2009,2010 Université Joseph Fourier (Grenoble I)

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
   \file test_mesh_codim1.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-08-13
 */
#define USE_BOOST_TEST 1
// Boost.Test
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE mesh_codim1 testsuite
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN

#include <testsuite/testsuite.hpp>


#include <feel/feelcore/feel.hpp>
#include <feel/feelmesh/geoentity.hpp>
#include <feel/feelmesh/refentity.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/importergmsh.hpp>

namespace Feel
{
namespace detail
{
typedef Mesh<Simplex<1, 1, 2> > mesh_type;
typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
}
Feel::detail::mesh_ptrtype
createMesh( double hsize )
{
    BOOST_TEST_MESSAGE( "create mesh" );
    double meshSize = hsize;
    //std::cout << "hsize = " << meshSize << std::endl;

    Gmsh __gmsh;
    std::string fname;
    std::ostringstream ostr;
    std::ostringstream nameStr;

    ostr << "Mesh.MshFileVersion = 2;\n"
         << "h=" << meshSize << ";\n"
         << "Point(1) = {-1, 0,0.0,h};\n"
         << "Point(2) = { 1, 0,0.0,h};\n"
         << "Line(1) = {1,2};\n"
         << "Line Loop(2) = {1};\n"
         << "Physical Point(31) = {1};\n"
         << "Physical Point(32) = {2};\n"
         << "Physical Line(3) = {1};\n";


    nameStr << "line." << meshSize;
    std::cout <<"Mesh generation ... \n";
    fname = __gmsh.generate( nameStr.str(), ostr.str() ).get<0>();

    /* Mesh */


    Feel::detail::mesh_ptrtype mesh( new Feel::detail::mesh_type );

    ImporterGmsh<Feel::detail::mesh_type> import( fname );
    import.setVersion( "2.0" );
    mesh->accept( import );
    mesh->components().set( MESH_CHECK | MESH_RENUMBER | MESH_UPDATE_EDGES | MESH_UPDATE_FACES );
    mesh->updateForUse();
    BOOST_TEST_MESSAGE( "create mesh done" );
    return mesh;
}
}
struct test_mesh_filters
{
    test_mesh_filters( double meshSize_=1 ): meshSize( meshSize_ ), mesh( Feel::createMesh( meshSize ) )
    {
    }
    void operator()()
    {
        BOOST_TEST_MESSAGE( "test_mesh_filters starts" );
        Feel::Assert::setLog( "test_mesh_filters.assert" );


        using namespace Feel;

        // elements
        {

            Feel::detail::mesh_type::gm_ptrtype __gm = mesh->gm();
            typedef Feel::detail::mesh_type::gm_type gm_type;
            typedef gm_type::precompute_ptrtype geopc_ptrtype;
            typedef gm_type::precompute_type geopc_type;
            typedef gm_type::Context<vm::POINT, Feel::detail::mesh_type::element_type> gmc_type;
            typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
            BOOST_TEST_MESSAGE( "test_mesh_filters check elements" );
            //
            // Precompute some data in the reference element for
            // geometric mapping and reference finite element
            //1
            Feel::detail::mesh_type::reference_convex_type ref_conv;
            geopc_ptrtype __geopc( new geopc_type( __gm, ref_conv.points() ) );
            Feel::MeshTraits<Feel::detail::mesh_type>::element_const_iterator it = mesh->beginElement();
            Feel::MeshTraits<Feel::detail::mesh_type>::element_const_iterator en = mesh->endElement();

            BOOST_TEST_MESSAGE( "Checking " << std::distance( it, en ) << " elements...\n" );

            for ( ; it != en; ++it )
            {
                // check that the geometric transformation from
                // the current gives back the vertices of the
                // element
                gmc_ptrtype __c( new gmc_type( __gm, *it, __geopc ) );

                FEELPP_ASSERT( ublas::norm_frobenius( __c->xReal() - it->G() ) < 1e-15 )( it->id() )( __c->xReal() )( it->G() ).error( "invalid element" );
                FEELPP_ASSERT( it->marker().value() == 3 )( it->id() )( it->marker().value() ).error( "invalid element marker" );

                //BOOST_CHECK_SMALL( ublas::norm_frobenius( __c->xReal() - it->G() ), 1e-15 );
                BOOST_CHECK_EQUAL( it->marker().value(), 3 );

            }
        }
        // location faces
        {

            Feel::MeshTraits<Feel::detail::mesh_type>::location_face_const_iterator it = mesh->beginFaceOnBoundary();
            Feel::MeshTraits<Feel::detail::mesh_type>::location_face_const_iterator en = mesh->endFaceOnBoundary();

            BOOST_TEST_MESSAGE( "Checking " << std::distance( it, en ) << " boundary faces...\n" );

            for ( ; it != en; ++it )
            {
#if defined(USE_BOOST_TEST)
                BOOST_CHECK( it->isConnectedTo0() &&
                             !it->isConnectedTo1() );
                BOOST_CHECK( it->marker().value() == 31 ||
                             it->marker().value() == 32 ||
                             it->marker().value() == 33 );
#endif
            }


            it = mesh->beginInternalFace();
            en = mesh->endInternalFace();
            std::cout << "Checking " << std::distance( it, en ) << " internal faces...\n";

#if defined(USE_BOOST_TEST)
            //BOOST_CHECK( std::distance( it, en ) == 1 );
#endif

            for ( ; it != en; ++it )
            {
#if defined(USE_BOOST_TEST)
                // the face must be connected with two elements
                BOOST_CHECK( it->isConnectedTo0() &&
                             it->isConnectedTo1() );
#endif
                // check that the points coordinates are the same for the face vertices
                int face_0 = it->pos_first();
                int face_1 = it->pos_second();
                Feel::node<double>::type n00 = it->element( 0 ).point( it->element( 0 ).fToP( face_0, 0 ) ).node();
                Feel::node<double>::type n10 = it->element( 1 ).point( it->element( 1 ).fToP( face_1, 0 ) ).node();
                FEELPP_ASSERT( ublas::norm_2( n00 - n10 ) < 1e-15 )
                ( it->id() )
                ( it->element( 0 ).G() )
                ( face_0 )
                ( it->element( 0 ).fToP( face_0, 0 ) )
                ( it->element( 1 ).G() )
                ( face_1 )
                ( it->element( 1 ).fToP( face_1, 0 ) )
                ( n00 )
                ( n10 )
                ( ublas::norm_2( n00 - n10 ) ).warn( "check failed" );

                if ( !( ublas::norm_2( n00 - n10 ) < 1e-15 ) )
                {
                    std::cout << "n00= " << n00 << "\n"
                              << "n10= " << n10 << "\n";
                }

#if defined(USE_BOOST_TEST)
                BOOST_CHECK_SMALL( ublas::norm_2( n00 - n10 ), 1e-15 );
#endif
                //FEELPP_ASSERT( it->element(0).nGeometricFaces() == 2 )( it->element(0).nGeometricFaces() ).error( "invalid number of faces" );
                Feel::node<double>::type n01 = it->element( 0 ).point( it->element( 0 ).fToP( face_0, 0 ) ).node();
                Feel::node<double>::type n11 = it->element( 1 ).point( it->element( 1 ).fToP( face_1, 0 ) ).node();
                FEELPP_ASSERT( ublas::norm_2( n01 - n11 ) < 1e-15 )
                ( it->id() )
                ( it->element( 0 ).G() )
                ( face_0 )
                ( it->element( 0 ).fToP( face_0, 0 ) )
                ( it->element( 1 ).G() )
                ( face_1 )
                ( it->element( 1 ).fToP( face_1, 0 ) )
                ( face_1 )
                ( n01 )( n11 )( ublas::norm_2( n01 - n11 ) ).warn( "check failed" );
#if defined(USE_BOOST_TEST)
                BOOST_CHECK_SMALL( ublas::norm_2( n01 - n11 ), 1e-15 );

                BOOST_CHECK_EQUAL( it->marker().value(), 3 );
#endif
            }
        }



    }
    double meshSize;
    Feel::detail::mesh_ptrtype mesh;
};

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( mesh_codim1 )

BOOST_AUTO_TEST_CASE( test_simple_mesh1d )
{
    using namespace Feel;


    typedef Mesh<Simplex<2,1> >  mesh_type;
    mesh_type mesh;


#if 0
    std::vector<Feel::detail::mesh_type::point_type*> vector_p;
    boost::sub_range<std::vector<Feel::detail::mesh_type::point_type*> > range_p;
    std::vector<std::vector<Feel::detail::mesh_type::point_type*> > vector_face_p( 4 );

    node_type n1( 2 );
    n1( 0 ) = 0;
    n1( 1 ) = 0;
    vector_p.push_back( mesh.add( Feel::detail::mesh_type::point_type( 0, n1, true, 0 ) ) );

    node_type n2( 2 );
    n2( 0 ) = 1;
    n2( 1 ) = 0;
    vector_p.push_back( mesh.add( Feel::detail::mesh_type::point_type( 1, n2, true, 1 ) ) );

    node_type n3( 2 );
    n3( 0 ) = 0;
    n3( 1 ) = 1;
    vector_p.push_back( mesh.add( Feel::detail::mesh_type::point_type( 2, n3, true, 2 ) ) );

    node_type n4( 2 );
    n4( 0 ) = 1;
    n4( 1 ) = 1;
    vector_p.push_back( mesh.add( Feel::detail::mesh_type::point_type( 3, n4, true, 2 ) ) );

    // face points
    vector_face_p[0].push_back( vector_p[0] );
    vector_face_p[0].push_back( vector_p[1] );
    vector_face_p[1].push_back( vector_p[1] );
    vector_face_p[1].push_back( vector_p[2] );
    vector_face_p[2].push_back( vector_p[2] );
    vector_face_p[2].push_back( vector_p[0] );

    Feel::detail::mesh_type::element_type* elt;
    elt = mesh.add( Feel::detail::mesh_type::element_type( 0, range_p( vector_p.begin(), vector_p.end() ), 0 ) );

    Feel::detail::mesh_type::face_type* f1;
    f1 = mesh.add( Feel::detail::mesh_type::face_type( 0, range_p( vector_face_p[0].begin(), vector_face_p[0].end() ), 0 ) );

    Feel::detail::mesh_type::face_type* f2;
    f2 = mesh.add( Feel::detail::mesh_type::face_type( 1, range_p( vector_face_p[1].begin(), vector_face_p[1].end() ), 1 ) );

    Feel::detail::mesh_type::face_type* f3;
    f3 = mesh.add( Feel::detail::mesh_type::face_type( 1, range_p( vector_face_p[2].begin(), vector_face_p[2].end() ), 2 ) );

#endif
    //mesh.updateElementFaces();

}



BOOST_AUTO_TEST_CASE( test_mesh_filters_ )
{

    test_mesh_filters tmf1;
    //test_mesh_filters tmf1; tmf1();

    //test_mesh_filters tmf2(0.2);

}

BOOST_AUTO_TEST_SUITE_END()
