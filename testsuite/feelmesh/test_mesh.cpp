/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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

// give a name to the testsuite
#define BOOST_TEST_MODULE mesh testsuite
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN
#include <fmt/core.h>
#include <fmt/chrono.h>
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelmesh/geoentity.hpp>
#include <feel/feelmesh/refentity.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feelfilters/loadmesh.hpp>

namespace Feel
{
namespace detail
{
typedef Mesh<Simplex<2> > mesh_type;
typedef std::shared_ptr<mesh_type> mesh_ptrtype;
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


        using namespace Feel;

        BOOST_TEST_MESSAGE( "testing mesh faces" );
        // location faces
        {
            auto rangeInternalFaces = mesh->internalFaces();
            auto it = std::get<0>( rangeInternalFaces );
            auto en = std::get<1>( rangeInternalFaces );
            //BOOST_CHECK( std::distance( it, en ) == 1 );
            for ( ; it != en; ++it )
            {
                auto const& iface = boost::unwrap_ref( *it );
                // the face must be connected with two elements
                BOOST_CHECK( iface.isConnectedTo0() &&
                             iface.isConnectedTo1() );
                // check that the points coordinates are the same for the face vertices
                int face_0 = iface.pos_first();
                int face_1 = iface.pos_second();
                Feel::node<double>::type n00 = iface.element( 0 ).point( iface.element( 0 ).fToP( face_0, 0 ) ).node();
                Feel::node<double>::type n10 = iface.element( 1 ).point( iface.element( 1 ).fToP( face_1, 1 ) ).node();
                FEELPP_ASSERT( ublas::norm_2( n00 - n10 ) < 1e-15 )
                ( iface.id() )
                ( iface.element( 0 ).G() )
                ( face_0 )
                ( iface.element( 0 ).fToP( face_0, 0 ) )
                ( iface.element( 1 ).G() )
                ( face_1 )
                ( iface.element( 1 ).fToP( face_1, 1 ) )
                ( n00 )
                ( n10 )
                ( ublas::norm_2( n00 - n10 ) ).warn( "check failed" );
                BOOST_CHECK( ublas::norm_2( n00 - n10 ) < 1e-15 );
                Feel::node<double>::type n01 = iface.element( 0 ).point( iface.element( 0 ).fToP( face_0, 1 ) ).node();
                Feel::node<double>::type n11 = iface.element( 1 ).point( iface.element( 1 ).fToP( face_1, 0 ) ).node();
                FEELPP_ASSERT( ublas::norm_2( n01 - n11 ) < 1e-15 )
                ( iface.id() )
                ( iface.element( 0 ).G() )
                ( face_0 )
                ( iface.element( 0 ).fToP( face_0, 1 ) )
                ( iface.element( 1 ).G() )
                ( face_1 )
                ( iface.element( 1 ).fToP( face_1, 0 ) )
                ( face_1 )
                ( n01 )( n11 )( ublas::norm_2( n01 - n11 ) ).warn( "check failed" );
                BOOST_CHECK( ublas::norm_2( n01 - n11 ) < 1e-15 );
            }

            auto rangeBoundaryFaces = mesh->facesOnBoundary();
            it = std::get<0>( rangeBoundaryFaces );
            en = std::get<1>( rangeBoundaryFaces );
            //BOOST_CHECK( std::distance( it, en ) == 4 );
            for ( ; it != en; ++it )
            {
                auto const& bface = boost::unwrap_ref( *it );
                BOOST_CHECK( bface.isConnectedTo0() &&
                             !bface.isConnectedTo1() );
                BOOST_CHECK( bface.marker().value() == mesh->markerName("Gamma1") ||
                             bface.marker().value() == mesh->markerName("Gamma2") ||
                             bface.marker().value() == mesh->markerName("Gamma3") );
            }
        }
        BOOST_TEST_MESSAGE( "testing mesh elements" );
        // elements
        {
            Feel::detail::mesh_type::gm_ptrtype __gm = mesh->gm();
            //
            // Precompute some data in the reference element for
            // geometric mapping and reference finite element
            //1
            Feel::detail::mesh_type::reference_convex_type ref_conv;
            auto __geopc = __gm->preCompute( ref_conv.points() );
            Feel::MeshTraits<Feel::detail::mesh_type>::element_const_iterator it = mesh->beginElement();
            Feel::MeshTraits<Feel::detail::mesh_type>::element_const_iterator en = mesh->endElement();

            //BOOST_CHECK( std::distance( it, en ) == 1 );
            for ( ; it != en; ++it )
            {
                auto const& elt = it->second;
                // check that the geometric transformation from
                // the current gives back the vertices of the
                // element
                auto __c = __gm->template context<vm::POINT>( elt, __geopc );

                BOOST_CHECK( ublas::norm_frobenius( __c->xReal() - elt.G() ) < 1e-15 );
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
typedef boost::mpl::list<
    //    std::pair<boost::mpl::int_<2>, boost::mpl::int_<2>>,
    std::pair<boost::mpl::int_<1>, boost::mpl::int_<3>>,
    std::pair<boost::mpl::int_<2>, boost::mpl::int_<3>>,
    std::pair<boost::mpl::int_<3>, boost::mpl::int_<3>>>
    dim_types;
BOOST_AUTO_TEST_CASE_TEMPLATE( test_mesh_measure, T, dim_types )
{
    using namespace Feel;
    constexpr int topodim = T::first_type::value;
    constexpr int realdim = T::second_type::value;
    using mesh_t = Mesh<Simplex<topodim,1,realdim>>;
    using mesh_ptr_t = std::shared_ptr<mesh_t>;
    std::string desc = fmt::format( R"(
        SetFactory( "OpenCASCADE" );
        h = {};
        Box( 1 ) = {{ 0, 0, 0, 1, 1, 1 }};
        Characteristic Length{{ PointsOf{{ Surface{{ 1:6 }}; }} }} = h;
        Physical Volume( "Box" ) = {{ 1 }};
        Physical Surface( "Walls" ) = {{  1:5 }};
        Physical Curve( "Edges" ) = {{  1:12 }};
        Mesh {};)", doption("gmsh.hsize"), topodim);
    BOOST_TEST_MESSAGE( desc );

    auto [ fname, generated ] = Gmsh{topodim}.generate( fmt::format("hypercube{}{}",topodim,realdim), desc );
    BOOST_TEST_MESSAGE( fmt::format( "mesh filename: {}", fname ) );
    auto mesh = loadMesh(_mesh=new mesh_t,_filename=fname);
    BOOST_CHECK( nelements(elements(mesh)) > 0 );
    auto beg = std::chrono::high_resolution_clock::now();
    double meas = 0, meas_int = 0, meas_surf = 0, meas_surf_int = 0;
    for ( auto const& wf : elements(mesh) )
    {
        auto const& f = boost::unwrap_ref( wf );
        meas += f.measure();
    }
    auto end = std::chrono::high_resolution_clock::now();
    BOOST_TEST_MESSAGE( fmt::format( "meas : {} in {}", meas, std::chrono::duration_cast<std::chrono::milliseconds>( end - beg ) ) );
    beg = std::chrono::high_resolution_clock::now();
    meas_int = integrate( _range = elements( mesh ), _expr = cst( 1. ) ).evaluate()( 0, 0 );
    end = std::chrono::high_resolution_clock::now();
    BOOST_TEST_MESSAGE( fmt::format( "meas int: {} in {}", meas_int, std::chrono::duration_cast<std::chrono::milliseconds>( end - beg ) ) );
    if constexpr ( topodim > 1 )
    {
        beg = std::chrono::high_resolution_clock::now();
        for ( auto const& wf : boundaryfaces( mesh ) )
        {
            auto const& f = boost::unwrap_ref( wf );
            meas_surf += f.measure();
        }
        end = std::chrono::high_resolution_clock::now();
        BOOST_TEST_MESSAGE( fmt::format( "meas surf: {} in {}", meas_surf, std::chrono::duration_cast<std::chrono::milliseconds>( end - beg ) ) );

        beg = std::chrono::high_resolution_clock::now();
        meas_surf_int = integrate( _range = boundaryfaces( mesh ), _expr = cst( 1. ) ).evaluate()( 0, 0 );
        end = std::chrono::high_resolution_clock::now();
        BOOST_TEST_MESSAGE( fmt::format( "meas surf int: {} in {}", meas_surf_int, std::chrono::duration_cast<std::chrono::milliseconds>( end - beg ) ) );
    }
    
    
    if constexpr ( topodim == 1 && realdim == 3 )
    {
        BOOST_CHECK_CLOSE( meas, 12, 1e-10 );
        BOOST_CHECK_SMALL( meas_surf, 1e-15 );
    }
    if constexpr( topodim == 2 && realdim == 3)
    {
        // cube with one face removed
        BOOST_CHECK_CLOSE( meas, 5, 1e-10);
        BOOST_CHECK_CLOSE( meas_surf, 4, 1e-10 );
    }
    if constexpr ( topodim == 3 && realdim == 3 )
    {
        BOOST_CHECK_CLOSE( meas, 1, 1e-10 );
        BOOST_CHECK_CLOSE( meas_surf, 6, 1e-10 );
    }
    BOOST_CHECK_CLOSE( meas, meas_int, 1e-12);
    if ( meas_surf > 1e-10 )
        BOOST_CHECK_CLOSE( meas_surf, meas_surf_int, 1e-12 );
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

