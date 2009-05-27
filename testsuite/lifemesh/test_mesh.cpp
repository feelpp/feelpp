/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-09-03

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2009 Université de Grenoble 1 (Joseph Fourier)

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
   \file test_mesh.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-09-03
 */
// Boost.Test
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

using boost::unit_test::test_suite;

#include <life/lifecore/life.hpp>
#include <life/lifemesh/geoentity.hpp>
#include <life/lifemesh/refentity.hpp>

#include <life/lifediscr/mesh.hpp>
#include <life/lifemesh/filters.hpp>
#include <life/lifefilters/importergmsh.hpp>
#include <life/lifefilters/gmsh.hpp>

namespace Life
{
typedef Mesh<Simplex<2> > mesh_type;
typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

}
struct test_mesh_filters
{
    test_mesh_filters( double meshSize_=1 ): meshSize(meshSize_), mesh()
    {
        BOOST_TEST_MESSAGE( "setup mesh" );
        BOOST_CHECK( meshSize_ <= 1 );

        mesh = this->createMesh( meshSize );

        BOOST_CHECK( mesh != 0 );
        BOOST_TEST_MESSAGE( "setup mesh done" );
    }
    Life::mesh_ptrtype
    createMesh( double hsize )
    {
        using namespace Life;
        double meshSize = hsize;
        //std::cout << "hsize = " << meshSize << std::endl;

        Gmsh __gmsh;
        std::string fname;
        std::ostringstream ostr;
        std::ostringstream nameStr;

        BOOST_TEST_CHECKPOINT( "Gmsh generator instantiated" );

        ostr << "Mesh.MshFileVersion = 2;\n"
             << "h=" << meshSize << ";\n"
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

        BOOST_TEST_CHECKPOINT( "Described mesh geometry" );

        nameStr << "triangle." << meshSize;
        BOOST_TEST_CHECKPOINT( "Described mesh name" );
        try {
            fname = __gmsh.generate( nameStr.str(), ostr.str() );
        }
        catch( ... )
            {
                std::cout << "Caught exception\n";
            }
        BOOST_TEST_CHECKPOINT( "Generating mesh with h=" << meshSize );
        /* Mesh */


        mesh_ptrtype mesh( new mesh_type );
        BOOST_TEST_CHECKPOINT( "Instantiating mesh" );


        ImporterGmsh<mesh_type> import( fname );
        BOOST_TEST_CHECKPOINT( "Importer instantiated" );
        import.setVersion( "2.0" );
        mesh->accept( import );
        BOOST_TEST_CHECKPOINT( "mesh imported" );
        mesh->components().set( MESH_CHECK | MESH_RENUMBER | MESH_UPDATE_EDGES | MESH_UPDATE_FACES );
        mesh->updateForUse();
        BOOST_TEST_CHECKPOINT( "mesh ready for use" );
        return mesh;
    }

    void operator()()
    {
        BOOST_TEST_MESSAGE( "testing mesh for h=" << meshSize );
        Life::Assert::setLog( "test_mesh_filters.assert" );


        using namespace Life;

        BOOST_TEST_MESSAGE( "testing mesh faces" );
        // location faces
        {
            Life::MeshTraits<mesh_type>::location_face_const_iterator it = mesh->beginInternalFace();
            Life::MeshTraits<mesh_type>::location_face_const_iterator en = mesh->endInternalFace();
            //BOOST_CHECK( std::distance( it, en ) == 1 );
            for( ; it != en; ++it )
                {
                    // the face must be connected with two elements
                    BOOST_CHECK( it->isConnectedTo0() &&
                                 it->isConnectedTo1() );
                    // check that the points coordinates are the same for the face vertices
                    int face_0 = it->pos_first();
                    int face_1 = it->pos_second();
                    Life::node<double>::type n00 = it->element(0).point( it->element(0).fToP( face_0, 0 ) ).node();
                    Life::node<double>::type n10 = it->element(1).point( it->element(1).fToP( face_1, 1 ) ).node();
                    LIFE_ASSERT( ublas::norm_2( n00 - n10 ) < 1e-15 )
                        ( it->id() )
                        ( it->element( 0 ).G() )
                        ( face_0 )
                        ( it->element(0).fToP( face_0, 0 ) )
                        ( it->element( 1 ).G() )
                        ( face_1 )
                        ( it->element(1).fToP( face_1, 1 ) )
                        ( n00 )
                        ( n10 )
                        ( ublas::norm_2( n00 - n10 ) ).warn( "check failed" );
                    BOOST_CHECK( ublas::norm_2( n00 - n10 ) < 1e-15 );
                    Life::node<double>::type n01 = it->element(0).point( it->element(0).fToP( face_0, 1 ) ).node();
                    Life::node<double>::type n11 = it->element(1).point( it->element(1).fToP( face_1, 0 ) ).node();
                    LIFE_ASSERT( ublas::norm_2( n01 - n11 ) < 1e-15 )
                        ( it->id() )
                        ( it->element( 0 ).G() )
                        ( face_0 )
                        ( it->element(0).fToP( face_0, 1 ) )
                        ( it->element( 1 ).G() )
                        ( face_1 )
                        ( it->element(1).fToP( face_1, 0 ) )
                        ( face_1 )
                        ( n01 )( n11 )( ublas::norm_2( n01 - n11 ) ).warn( "check failed" );
                    BOOST_CHECK( ublas::norm_2( n01 - n11 ) < 1e-15 );
                }
            it = mesh->beginFaceOnBoundary();
            en = mesh->endFaceOnBoundary();
            //BOOST_CHECK( std::distance( it, en ) == 4 );
            for( ; it != en; ++it )
                {
                    BOOST_CHECK( it->isConnectedTo0() &&
                                 !it->isConnectedTo1() );
                    BOOST_CHECK( it->marker().value() == 31 ||
                                 it->marker().value() == 32 ||
                                 it->marker().value() == 33 );
                }
        }
        BOOST_TEST_MESSAGE( "testing mesh elements" );
        // elements
        {
            mesh_type::gm_ptrtype __gm = mesh->gm();
            typedef mesh_type::gm_type gm_type;
            typedef gm_type::precompute_ptrtype geopc_ptrtype;
            typedef gm_type::precompute_type geopc_type;
            typedef gm_type::Context<vm::POINT, mesh_type::element_type> gmc_type;
            typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
            //
            // Precompute some data in the reference element for
            // geometric mapping and reference finite element
            //1
            mesh_type::reference_convex_type ref_conv;
            geopc_ptrtype __geopc( new geopc_type( __gm, ref_conv.points() ) );
            Life::MeshTraits<mesh_type>::element_const_iterator it = mesh->beginElement();
            Life::MeshTraits<mesh_type>::element_const_iterator en = mesh->endElement();
            //BOOST_CHECK( std::distance( it, en ) == 1 );
            for( ; it != en; ++it )
                {
                    // check that the geometric transformation from
                    // the current gives back the vertices of the
                    // element
                    gmc_ptrtype __c( new gmc_type( __gm, *it, __geopc ) );

                    BOOST_CHECK( ublas::norm_frobenius( __c->xReal() - it->G() ) < 1e-15 );
                }
        }
        BOOST_TEST_MESSAGE( "testing mesh for h=" << meshSize << " done" );
    }
    double meshSize;
    Life::mesh_ptrtype mesh;
};
BOOST_AUTO_TEST_CASE( test_mesh_filters_ )
{
    test_mesh_filters tmf;
    tmf();
}

BOOST_AUTO_TEST_CASE( test_simple_mesh2d )
{
    using namespace Life;

    typedef Mesh<Simplex<2,1> >  mesh_type;
    mesh_type mesh;


#if 0
    std::vector<mesh_type::point_type*> vector_p;
    boost::sub_range<std::vector<mesh_type::point_type*> > range_p;
    std::vector<std::vector<mesh_type::point_type*> > vector_face_p( 4 );

    node_type n1( 2 );
    n1( 0 ) = 0;
    n1( 1 ) = 0;
    vector_p.push_back( mesh.add( mesh_type::point_type( 0, n1, true, 0 ) ) );

    node_type n2( 2 );
    n2( 0 ) = 1;
    n2( 1 ) = 0;
    vector_p.push_back( mesh.add( mesh_type::point_type( 1, n2, true, 1 ) ) );

    node_type n3( 2 );
    n3( 0 ) = 0;
    n3( 1 ) = 1;
    vector_p.push_back( mesh.add( mesh_type::point_type( 2, n3, true, 2 ) ) );

    node_type n4( 2 );
    n4( 0 ) = 1;
    n4( 1 ) = 1;
    vector_p.push_back( mesh.add( mesh_type::point_type( 3, n4, true, 2 ) ) );

    // face points
    vector_face_p[0].push_back( vector_p[0] );
    vector_face_p[0].push_back( vector_p[1] );
    vector_face_p[1].push_back( vector_p[1] );
    vector_face_p[1].push_back( vector_p[2] );
    vector_face_p[2].push_back( vector_p[2] );
    vector_face_p[2].push_back( vector_p[0] );

    mesh_type::element_type* elt;
    elt = mesh.add( mesh_type::element_type( 0, range_p( vector_p.begin(), vector_p.end() ), 0 ) );

    mesh_type::face_type* f1;
    f1 = mesh.add( mesh_type::face_type( 0, range_p( vector_face_p[0].begin(), vector_face_p[0].end() ), 0 ) );

    mesh_type::face_type* f2;
    f2 = mesh.add( mesh_type::face_type( 1, range_p( vector_face_p[1].begin(), vector_face_p[1].end() ), 1 ) );

    mesh_type::face_type* f3;
    f3 = mesh.add( mesh_type::face_type( 1, range_p( vector_face_p[2].begin(), vector_face_p[2].end() ), 2 ) );

#endif
    //mesh.updateElementFaces();

}
