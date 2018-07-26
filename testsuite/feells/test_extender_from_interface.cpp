/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 07 juil. 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
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



//#define USE_BOOST_TEST 1
#if defined(USE_BOOST_TEST)
#define BOOST_TEST_MODULE levelset
#include <testsuite.hpp>
#endif

#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feells/reinit_fms.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feells/extenderfrominterface.hpp>

using namespace Feel;


inline
po::options_description
makeOptions()
{
    po::options_description opts( "test_extenderFromInterface" );
    opts.add_options()
    ( "radius", po::value<double>()->default_value( 1.0 ), "circle or sphere radius" )
    ;
    return opts.add( Feel::feel_options() );
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_extenderFromInterface" ,
                     "test_extenderFromInterface" ,
                     "0.1",
                     "test extenderFromInterface",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2015 Feel++ Consortium" );

    about.addAuthor( "Vincent Doyeux", "developer", "vincent.doyeux@ujf-grenoble.fr", "" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );

    return about;
}

/// Test the reinitialisation method (Fast Marching, ...)
///     \tparam DIM         Topological dimension.
///     \tparam H_ORDER     Space polynomial order..
///     \tparam G_ORDER     Geometrical polynomial order.
template<int DIM, int H_ORDER, int G_ORDER>
class TestExtenderFromInterface
{
public:
    using mesh_type = Mesh< Simplex<DIM,G_ORDER> >;
    using mesh_ptrtype = std::shared_ptr< mesh_type >;

    /// Init the geometry with a circle/sphere from radius and characteristic length
    ///     \param radius   Circle or sphere radius.
    ///     \param h        Mesh size.
    TestExtenderFromInterface( double radius=doption("radius") ) :
        M_mesh( createGMSHMesh( _mesh=new mesh_type,
                                _desc=domain( _name="ellipsoid_nd",
                                              _shape="ellipsoid",
                                              _dim=DIM,
                                              _order=G_ORDER,
                                              _xmin=-radius,
                                              _ymin=-radius,
                                              _zmin=-radius,
                                              _xmax=radius,
                                              _ymax=radius,
                                              _zmax=radius,
                                              _h=doption("gmsh.hsize") ),
                                _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES
                              ) ),
        M_radius( radius )
    {}


    void operator()()
    {
        auto Xh = Pch<H_ORDER>(M_mesh);
        auto thefms = fms( Xh );
        auto phio = Xh->element();
        phio.on( _range=elements(M_mesh), _expr=h());
        phio.on( _range=boundaryfaces(M_mesh), _expr= -h()/100. );
        auto phi = thefms->march(phio);
        auto dist = vf::project( Xh, elements(M_mesh), M_radius - sqrt( Px()*Px()+Py()*Py()+Pz()*Pz() ) );
        auto err = vf::project( Xh, elements(M_mesh), abs( idv(phi) - idv(dist) ) );

        auto psi = Xh->element();
        psi.on( _range=boundaryelements(M_mesh), _expr=Px() );

        ExtenderFromInterface<DIM> ext( Xh,  phi );
        ext.extendFromInterface( psi );

        psi.printMatlab("psi.m" );
#if defined(USE_BOOST_TEST)
        // Max error around mesh size h.
        //BOOST_CHECK_CLOSE( phi.max(), M_radius, h() );
        BOOST_CHECK_SMALL( err.max(), h() )
#else

        LOG(INFO) << "phi1 max" << phi.max();
        LOG(INFO) << "err1 max: " << err.max();
        auto exp = exporter(_mesh=M_mesh, _name="extend");
        exp->add("phi", phi);
        exp->add("phio", phio);
        exp->add("dist", dist);
        exp->add("psi", psi);
        exp->add("marker", ext.marker());
        exp->save();
#endif
    }

    
private:
    mesh_ptrtype M_mesh;
    double M_radius;

};

#if defined(USE_BOOST_TEST)
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
BOOST_AUTO_TEST_SUITE( extenderFromInterface )

// Test wall distance on :
// Unit 2D circle P1G1.
BOOST_AUTO_TEST_CASE( test_2d_p1g1 )
{
    TestExtenderFromInterface<2,1,1> tls;
    tls();
}

// Unit 3D sphere P1G1.
BOOST_AUTO_TEST_CASE( test_3d_p1g1 )
{
    TestExtenderFromInterface<3,1,1> tls;
    tls();
}




BOOST_AUTO_TEST_SUITE_END()

#else
int main(int argc, char** argv )
{
    Feel::Environment env( _argc=argc, _argv=argv,
                           _about=makeAbout(),
                           _desc=makeOptions() );
#if 0
    TestExtenderFromInterface<2,1,1> tls;
    TestExtenderFromInterface<3,1,1> tls;
    TestExtenderFromInterface<2,1,1> tls( 10, 0.1 );
    TestExtenderFromInterface<3,1,1> tls( 10, 1 );
#endif
    //TestExtenderFromInterface<3,1,1> tls( 1, 0.1 );
    //TestExtenderFromInterface<2,1,2> tls( 1, 0.1 );
    TestExtenderFromInterface<2,1,1> tls;
    tls();

    return 0;
}
#endif


