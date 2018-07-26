/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
   -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-07

  Copyright (C) 2008-2012 Universite Joseph Fourier (Grenoble I)

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   \file test_upmarker_faces.cpp
   \author Guillaume Dolle <gdolle at unistra.fr>
   \date 2014-01-21
 */

//#define USE_BOOST_TEST 1
#if defined(USE_BOOST_TEST)
#define BOOST_TEST_MODULE onimpl
#include <testsuite.hpp>
#endif

#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

inline
po::options_description
makeOptions()
{
    po::options_description opts( "test_updatemarker2faces" );
    opts.add_options()
    ( "radius", po::value<double>()->default_value( 0.2 ), "circle or sphere radius" )
    ;
    return opts.add( Feel::feel_options() );
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_updatemarker2faces" ,
                     "test_updatemarker2faces" ,
                     "0.1",
                     "test onimpl",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2015 Feel++ Consortium" );

    about.addAuthor( "Guillaume Dolle", "developer", "gdolle@unistra.fr", "" );
    about.addAuthor( "Vincent Doyeux", "developer", "vincent.doyeux@ujf-grenoble.fr", "" );

    return about;
}

/// Method to create a marker on faces though a P0 discontinuous function.
///     \tparam DIM         Topological dimension.
///     \tparam H_ORDER     Space polynomial order..
///     \tparam G_ORDER     Geometrical polynomial order.
template<int DIM, int H_ORDER, int G_ORDER>
class TestUpMarker2Faces
{
public: 
    using mesh_type = Mesh< Simplex<DIM,G_ORDER> >;
    using mesh_ptrtype = std::shared_ptr< mesh_type >;
    using face_type=typename mesh_type::face_type;
    using cont_range_type=std::vector< typename mesh_type::face_type >;

    // Create a test on default cube geometry.
    TestUpMarker2Faces() :
        M_mesh( loadMesh( _mesh=new mesh_type ) )
    {}

    /// Create a marker on selected boundary faces.
    //
    /// It has no mathematical sense to project a P0 discontinuous function on
    /// boundary faces (no dof). Therefore to create a similar interface to
    /// `updateMarker2` the idea is:
    ///     * Project a P0 discontinuous function on boundary elements such that its
    ///       value denotes the marker id (marker2 or marker3 must be available);
    ///     * Update one of the available marker;
    ///     * Create a range (iterator) on faces from the marker2 (or 3) view.
    void test()
    {
        auto Xh = Pch<H_ORDER>(M_mesh);
        auto Xdh0 = Pdh<0>(M_mesh);
        auto mark = Xdh0->element();
        int fid = 1;
        double R = doption("radius");
        auto test_boundary = Xdh0->element();
        auto test_elements0 = Xdh0->element();
        auto parts = Xdh0->element();
        auto test_elements1 = Xh->element();
        auto test_faces = Xh->element();

        // Select elements in the circle with radius 0.5 at the center of 4 faces.
        //auto e = chi( ((Px()-0.5)*(Px()-0.5) + (Py()-0.5)*(Py()-0.5)) < R*R );
        auto e1 = chi( ((Px()-0.5)*(Px()-0.5) + (Py()-0.5)*(Py()-0.5)) < R*R );
        auto e2 = chi( ((Py()-0.5)*(Py()-0.5) + (Pz()-0.5)*(Pz()-0.5)) < R*R );

        int rank = Environment::worldComm().rank();
        parts.on( elements(M_mesh), cst(rank) );
        mark.on( boundaryelements( M_mesh ), chi( (e1 + e2) > 0 ) );

        M_mesh->updateMarker2( mark );
        M_mesh->addMarkerName( "select_elements", fid, DIM );

        auto range = marked2elements(M_mesh,"select_elements");

        auto it=begin(range);
        auto en=end(range);

        std::shared_ptr<cont_range_type> myfaces( new cont_range_type );

        // Create faces range: Iterate on selected elements, select only
        // boundary faces, add them in a vector.
        for(;it!=en;++it)
        {
            auto const& elt = boost::unwrap_ref( *it );
            for(uint16_type f=0;f<elt.numTopologicalFaces; ++f)
            {
                auto face = elt.face(f);
                if( face.isOnBoundary() )
                    myfaces->push_back(boost::cref(face));
            }
        }

        // Create a range of faces.
        auto myrangefaces = boost::make_tuple( mpl::size_t<MESH_FACES>(),
                                        myfaces->begin(),
                                        myfaces->end() );



        M_mesh->updateMarker2WithRangeFaces( myrangefaces, 666 );
        M_mesh->addMarkerName( "select_faces", 666, DIM-1 );

        test_boundary.on( boundaryelements(M_mesh), cst(42) );
        test_elements0.on( marked2elements(M_mesh,"select_elements"), cst(42) );
        test_elements1.on( marked2elements(M_mesh,"select_elements"), cst(42) );
        test_faces.on( marked2faces(M_mesh,"select_faces"), cst(42) );

#if defined(USE_BOOST_TEST)
        // Max error around mesh size h.
        //BOOST_CHECK_CLOSE( phi.max(), M_radius, h() );
        //BOOST_CHECK_SMALL( err.max(), h() )
#else
        auto fit = myfaces->begin();
        auto fen = myfaces->end();
#if 0
        for( ;fit!=fen;++fit)
        {
            std::cout << "rank: " << Environment::worldComm().rank() << " face id: " << fit->id() << "\n";
        }
#endif
        // ----------------------------------
        // CHECK
        int size=0;
        //for(;it!=en;++it)
        //    size++;
        it = boost::get<1>(range);
        for(;it!=en;++it)
        {
            auto const& elt = boost::unwrap_ref( *it );
            for(uint16_type f=0;f<elt.numTopologicalFaces; ++f)
            {
                auto face = elt.face(f);
                if( face.isOnBoundary() )
                {
                    std::cout
                        << "elt mk2: "<< elt.marker2().value()
                        << " elt mk3: "<< elt.marker3().value()
                        << " | face mk2: " << face.marker2().value()
                        << " | face mk3: " << face.marker3().value()
                        << "\n";
                    size++;
                }
            }
        }

        auto exp = exporter(_mesh=M_mesh, _name="testsuite_upmarker_faces");
        exp->step(0)->setMesh(M_mesh);
        exp->step(0)->add("test_elements_boundary", test_boundary);
        exp->step(0)->add("mark", mark);
        exp->step(0)->add("partitions", parts);
        exp->step(0)->add("test_elements_Xdh0", test_elements0);
        exp->step(0)->add("test_elements_Xh", test_elements1);
        exp->step(0)->add("test_faces", test_faces);
        exp->save();
#endif
    }

private:
    mesh_ptrtype M_mesh;
};

#if defined(USE_BOOST_TEST)
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
BOOST_AUTO_TEST_SUITE( upmarker2faces )

//BOOST_AUTO_TEST_CASE( test_1d )
//{
//    TestUpMarker2Faces<3,1,1> tumf;
//    tumf.test();
//}

BOOST_AUTO_TEST_SUITE_END()

#else
int main(int argc, char** argv )
{
    Feel::Environment env( _argc=argc, _argv=argv,
                           _about=makeAbout(),
                           _desc=makeOptions() );
    TestUpMarker2Faces<3,1,1> tumf;
    tumf.test();

    return 0;
}
#endif
