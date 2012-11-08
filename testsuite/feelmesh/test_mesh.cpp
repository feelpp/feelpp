/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-09-03

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2009 Université de Grenoble 1 (Joseph Fourier)

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
   \file test_mesh.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-09-03
 */

// Boost.Test
// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE mesh testsuite
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN

#include <testsuite/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelmesh/geoentity.hpp>
#include <feel/feelmesh/refentity.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/geotool.hpp>

namespace Feel
{
namespace detail
{
typedef Mesh<Simplex<2> > mesh_type;
typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
}
}
struct test_mesh_filters
{
    test_mesh_filters( double meshSize_=1 ): meshSize( meshSize_ ), mesh()
    {
        BOOST_TEST_MESSAGE( "setup mesh" );
        BOOST_CHECK( meshSize_ <= 1 );

        mesh = this->createMesh( meshSize );

        BOOST_CHECK( mesh != 0 );
        BOOST_TEST_MESSAGE( "setup mesh done" );
    }
    Feel::detail::mesh_ptrtype
    createMesh( double hsize )
    {
        BOOST_TEST_MESSAGE( "create mesh" );
        using namespace Feel;
        double meshSize = hsize;
        //std::cout << "hsize = " << meshSize << std::endl;

        Gmsh __gmsh;
        std::string fname;
        std::ostringstream ostr;
        std::ostringstream nameStr;

        BOOST_TEST_CHECKPOINT( "Gmsh generator instantiated" );

        GeoTool::Node x1(-1, -1);
        GeoTool::Node x2( 1, -1);
        GeoTool::Node x3(-1,  1);
        GeoTool::Triangle T( meshSize,"MyTriangle",x1,x2,x3);
        T.setMarker(_type="line",_name="Gamma1",_marker1=true);
        T.setMarker(_type="line",_name="Gamma2",_marker2=true);
        T.setMarker(_type="line",_name="Gamma3",_marker3=true);
        T.setMarker(_type="surface",_name="Omega",_markerAll=true);

        auto mesh = T.createMesh(_mesh = new Feel::detail::mesh_type,
                                 _name="triangle" );

        BOOST_TEST_CHECKPOINT( "mesh ready for use" );
        BOOST_TEST_MESSAGE( "create mesh done" );
        return mesh;
    }

    void operator()()
    {
        BOOST_TEST_MESSAGE( "testing mesh for h=" << meshSize );
        Feel::Assert::setLog( "test_mesh_filters.assert" );


        using namespace Feel;

        BOOST_TEST_MESSAGE( "testing mesh faces" );
        // location faces
        {
            Feel::MeshTraits<Feel::detail::mesh_type>::location_face_const_iterator it = mesh->beginInternalFace();
            Feel::MeshTraits<Feel::detail::mesh_type>::location_face_const_iterator en = mesh->endInternalFace();

            //BOOST_CHECK( std::distance( it, en ) == 1 );
            for ( ; it != en; ++it )
            {
                // the face must be connected with two elements
                BOOST_CHECK( it->isConnectedTo0() &&
                             it->isConnectedTo1() );
                // check that the points coordinates are the same for the face vertices
                int face_0 = it->pos_first();
                int face_1 = it->pos_second();
                Feel::node<double>::type n00 = it->element( 0 ).point( it->element( 0 ).fToP( face_0, 0 ) ).node();
                Feel::node<double>::type n10 = it->element( 1 ).point( it->element( 1 ).fToP( face_1, 1 ) ).node();
                FEELPP_ASSERT( ublas::norm_2( n00 - n10 ) < 1e-15 )
                ( it->id() )
                ( it->element( 0 ).G() )
                ( face_0 )
                ( it->element( 0 ).fToP( face_0, 0 ) )
                ( it->element( 1 ).G() )
                ( face_1 )
                ( it->element( 1 ).fToP( face_1, 1 ) )
                ( n00 )
                ( n10 )
                ( ublas::norm_2( n00 - n10 ) ).warn( "check failed" );
                BOOST_CHECK( ublas::norm_2( n00 - n10 ) < 1e-15 );
                Feel::node<double>::type n01 = it->element( 0 ).point( it->element( 0 ).fToP( face_0, 1 ) ).node();
                Feel::node<double>::type n11 = it->element( 1 ).point( it->element( 1 ).fToP( face_1, 0 ) ).node();
                FEELPP_ASSERT( ublas::norm_2( n01 - n11 ) < 1e-15 )
                ( it->id() )
                ( it->element( 0 ).G() )
                ( face_0 )
                ( it->element( 0 ).fToP( face_0, 1 ) )
                ( it->element( 1 ).G() )
                ( face_1 )
                ( it->element( 1 ).fToP( face_1, 0 ) )
                ( face_1 )
                ( n01 )( n11 )( ublas::norm_2( n01 - n11 ) ).warn( "check failed" );
                BOOST_CHECK( ublas::norm_2( n01 - n11 ) < 1e-15 );
            }

            it = mesh->beginFaceOnBoundary();
            en = mesh->endFaceOnBoundary();

            //BOOST_CHECK( std::distance( it, en ) == 4 );
            for ( ; it != en; ++it )
            {
                BOOST_CHECK( it->isConnectedTo0() &&
                             !it->isConnectedTo1() );
                BOOST_CHECK( it->marker().value() == mesh->markerName("Gamma1") ||
                             it->marker().value() == mesh->markerName("Gamma2") ||
                             it->marker().value() == mesh->markerName("Gamma3") );
            }
        }
        BOOST_TEST_MESSAGE( "testing mesh elements" );
        // elements
        {
            Feel::detail::mesh_type::gm_ptrtype __gm = mesh->gm();
            typedef Feel::detail::mesh_type::gm_type gm_type;
            typedef gm_type::precompute_ptrtype geopc_ptrtype;
            typedef gm_type::precompute_type geopc_type;
            typedef gm_type::Context<vm::POINT, Feel::detail::mesh_type::element_type> gmc_type;
            typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
            //
            // Precompute some data in the reference element for
            // geometric mapping and reference finite element
            //1
            Feel::detail::mesh_type::reference_convex_type ref_conv;
            geopc_ptrtype __geopc( new geopc_type( __gm, ref_conv.points() ) );
            Feel::MeshTraits<Feel::detail::mesh_type>::element_const_iterator it = mesh->beginElement();
            Feel::MeshTraits<Feel::detail::mesh_type>::element_const_iterator en = mesh->endElement();

            //BOOST_CHECK( std::distance( it, en ) == 1 );
            for ( ; it != en; ++it )
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
    Feel::detail::mesh_ptrtype mesh;
};

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( mesh )


BOOST_AUTO_TEST_CASE( test_mesh_filters_ )
{
    test_mesh_filters tmf;
    tmf();
}
BOOST_AUTO_TEST_CASE( test_mesh_comp )
{
    using namespace Feel;
    typedef Mesh<Simplex<2,1> >  mesh_type;
    mesh_type mesh;

    mesh.components().reset();
    BOOST_CHECK( mesh.components().test( MESH_CHECK ) == false );
    BOOST_CHECK( mesh.components().test( MESH_RENUMBER ) == false );
    BOOST_CHECK( mesh.components().test( MESH_UPDATE_FACES ) == false );
    BOOST_CHECK( mesh.components().test( MESH_UPDATE_EDGES ) == false );

    mesh.components().reset();
    mesh.components().set( MESH_CHECK );
    BOOST_TEST_MESSAGE( "check MESH_CHECK comp: " << mesh.components().context()  << "\n" );
    BOOST_CHECK( mesh.components().test( MESH_CHECK ) == true );
    BOOST_CHECK( mesh.components().test( MESH_RENUMBER ) == false );
    BOOST_CHECK( mesh.components().test( MESH_UPDATE_FACES ) == false );
    BOOST_CHECK( mesh.components().test( MESH_UPDATE_EDGES ) == false );


    mesh.components().reset();
    mesh.components().set( MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES );
    BOOST_TEST_MESSAGE( "check MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES comp: " << mesh.components().context() << "\n" );
    BOOST_CHECK( mesh.components().test( MESH_CHECK ) == true );
    BOOST_CHECK( mesh.components().test( MESH_RENUMBER ) == false );
    BOOST_CHECK( mesh.components().test( MESH_UPDATE_FACES ) == true );
    BOOST_CHECK( mesh.components().test( MESH_UPDATE_EDGES ) == true );


    mesh.components().reset();
    mesh.components().set( MESH_RENUMBER|MESH_UPDATE_FACES );
    BOOST_TEST_MESSAGE( "check MESH_RENUMBER|MESH_UPDATE_FACES comp: " << mesh.components().context() << "\n" );
    BOOST_CHECK( mesh.components().test( MESH_CHECK ) == false );
    BOOST_CHECK( mesh.components().test( MESH_RENUMBER ) == true );
    BOOST_CHECK( mesh.components().test( MESH_UPDATE_FACES ) == true );
    BOOST_CHECK( mesh.components().test( MESH_UPDATE_EDGES ) == false );


    mesh.components().reset();
    mesh.components().set( MESH_RENUMBER );
    BOOST_TEST_MESSAGE( "check MESH_RENUMBER comp: " << mesh.components().context() << "\n" );
    BOOST_CHECK( mesh.components().test( MESH_CHECK ) == false );
    BOOST_CHECK( mesh.components().test( MESH_RENUMBER ) == true );
    BOOST_CHECK( mesh.components().test( MESH_UPDATE_FACES ) == false );
    BOOST_CHECK( mesh.components().test( MESH_UPDATE_EDGES ) == false );
}
BOOST_AUTO_TEST_CASE( test_mesh_lmethod )
{

    using namespace Feel;
    typedef Mesh<Simplex<2,1> >  mesh_type;
    mesh_type mesh;

    auto pit = mesh.beginPoint();
    auto pen = mesh.endPoint();


    while ( pit != pen )
    {
#if 0
        auto eit = pit->beginElement();
        auto een = pit->endElement();

        while ( eit != een )
        {
            auto element = mesh->element( *eit );

            for ( int f = 0; f < eit->nDim; ++f )
            {
                // plocal local id of the vertex in the element
                auto face = mesh->face( element->v2f( plocal, f ) );


            }

            ++eit;
        }

#endif
        ++pit;
    }


}
BOOST_AUTO_TEST_CASE( test_simple_mesh2d )
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

BOOST_AUTO_TEST_SUITE_END()

