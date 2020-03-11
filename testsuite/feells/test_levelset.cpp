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
//#define USE_BOOST_TEST 1
#undef USE_BOOST_TEST
#if defined(USE_BOOST_TEST)
#define BOOST_TEST_MODULE levelset
#include <feel/feelcore/testsuite.hpp>
#endif

#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feells/reinit_fms.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/domain.hpp>

using namespace Feel;


inline
po::options_description
makeOptions()
{
    po::options_description opts( "test_levelset" );
    opts.add_options()
    ( "radius", po::value<double>()->default_value( 1.0 ), "circle or sphere radius" )
    ;
    return opts.add( Feel::feel_options() );
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_levelset" ,
                     "test_levelset" ,
                     "0.1",
                     "test levelset",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2015 Feel++ Consortium" );

    about.addAuthor( "Guillaume Dolle", "developer", "gdolle@unistra.fr", "" );
    about.addAuthor( "Vincent Doyeux", "developer", "vincent.doyeux@ujf-grenoble.fr", "" );

    return about;
}

/// Test the reinitialisation method (Fast Marching, ...)
///     \tparam DIM         Topological dimension.
///     \tparam H_ORDER     Space polynomial order..
///     \tparam G_ORDER     Geometrical polynomial order.
template<int DIM, int H_ORDER, int G_ORDER>
class TestLevelSet
{
public:
    using mesh_type = Mesh< Simplex<DIM,G_ORDER> >;
    using mesh_ptrtype = std::shared_ptr< mesh_type >;

    /// Init the geometry with a circle/sphere from radius and characteristic length
    ///     \param radius   Circle or sphere radius.
    ///     \param h        Mesh size.
    TestLevelSet( double radius=doption("radius") ) :
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

    /// First method: let the FMM search for the elements crossed by the interface
    /// and the maximum distance is checked for a circle/sphere.
    void wallDist_1()
    {
        auto Xh = Pch<H_ORDER>(M_mesh);
        auto thefms = fms( Xh );
        auto phio = Xh->element();
        phio = vf::project(_space=Xh, _range=elements(M_mesh), _expr=h() );
        phio.on( _range=boundaryfaces(M_mesh), _expr= -h()/100. );
        auto phi = thefms->march(phio);
        auto dist = vf::project( _space=Xh,  _range=elements(M_mesh), _expr=M_radius - sqrt( Px()*Px()+Py()*Py()+Pz()*Pz() ) );
        auto err = vf::project( _space=Xh,  _range=elements(M_mesh), _expr=abs( idv(phi) - idv(dist) ) );

#if defined(USE_BOOST_TEST)
        // Max error around mesh size h.
        //BOOST_CHECK_CLOSE( phi.max(), M_radius, h() );
        BOOST_CHECK_SMALL( err.max(), h() )
#else

        LOG(INFO) << "phi1 max" << phi.max();
        LOG(INFO) << "err1 max: " << err.max();
        auto exp = exporter(_mesh=M_mesh, _name="testsuite_levelset_distw1");
        exp->step(0)->add("phi1", phi);
        exp->step(0)->add("phio1", phio);
        exp->step(0)->add("dist1", dist);
        exp->step(0)->add("err1", err);
        exp->save();
#endif
    }

    /// Second method: give to the FMM the elements having the good value to start with,
    /// and the maximum distance is checked for a circle/sphere.
    void wallDist_2()
    {
        auto Xh = Pch<H_ORDER>(M_mesh);
        auto Xh0 = Pdh<0>(M_mesh);
        auto thefms = fms( Xh );
        auto phio = Xh->element();
        auto mark = Xh0->element();

        phio.on( _range=boundaryelements(M_mesh), _expr=h() );
        phio.on( _range=boundaryfaces(M_mesh), _expr=cst(0) );

        mark.on( _range=boundaryelements(M_mesh), _expr=cst(1) );
        M_mesh->updateMarker2( mark );

        auto phi = thefms->march(phio, boundaryelements(M_mesh));
        auto dist = vf::project( _space=Xh, _range=elements(M_mesh), _expr=M_radius - sqrt( Px()*Px()+Py()*Py()+Pz()*Pz() ) );
        auto err = vf::project( _space=Xh, _range=elements(M_mesh), _expr=abs( idv(phi) - idv(dist) ) );

#if defined(USE_BOOST_TEST)
        // Max error around mesh size h.
        //BOOST_CHECK_CLOSE( phi.max(), M_radius, h() );
        BOOST_CHECK_SMALL( err.max(), h() )
#else
        LOG(INFO) << "phi2_max" << phi.max();
        LOG(INFO) << "err2 max: " << err.max();
        auto exp = exporter(_mesh=M_mesh, _name="testsuite_levelset_distw2");
        exp->step(0)->add("phio2", phio);
        exp->step(0)->add("phi2", phi);
        exp->step(0)->add("dist2", dist);
        exp->step(0)->add("err2", err);
        exp->step(0)->add("mark2", mark);
        exp->save();
#endif
    }

    /// Third test: distance from a circle inside the domain
    void circle()
    {
        auto Xh = Pch<H_ORDER>(M_mesh);
        auto Xh0 = Pdh<0>(M_mesh);
        auto thefms = fms( Xh );
        auto x_center =  0.;
        auto y_center =  0.;
        auto radius = doption("radius")*0.25;

        auto phio = vf::project(_space=Xh, _range=elements(M_mesh), _expr=(Px()-x_center)*(Px()-x_center) + (Py()-y_center)*(Py()-y_center) - radius);
        auto dist = vf::project(_space=Xh, _range=elements(M_mesh), _expr=sqrt((Px()-x_center)*(Px()-x_center) + (Py()-y_center)*(Py()-y_center)) - radius);
        
        auto phi1 = thefms->march(phio) ;

        auto mark = vf::project(_space=Xh0, _range=elements(M_mesh), _expr=chi(abs(idv(phio)) <= 2.*doption("gmsh.hsize"))) ;
        M_mesh->updateMarker2( mark ) ;

        auto phi2 = thefms->march(phio, marked2elements(M_mesh, 1)) ;
        
        auto err1 = vf::project( _space=Xh, _range=elements(M_mesh), _expr=abs( idv(phi1) - idv(dist) ) );
        auto err2 = vf::project( _space=Xh, _range=elements(M_mesh), _expr=abs( idv(phi2) - idv(dist) ) );
#if defined(USE_BOOST_TEST)
        // Max error around mesh size h.
        //BOOST_CHECK_CLOSE( phi.max(), M_radius, h() );
        BOOST_CHECK_SMALL( err1.max(), h() )
        BOOST_CHECK_SMALL( err2.max(), h() )
#else
        auto exp = exporter(_mesh=M_mesh, _name="testsuite_levelset_circle2");
        for(int i = 0; i < 3; i++){
        exp->step(i)->add("phio", phio);
        exp->step(i)->add("phi1", phi1);
        exp->step(i)->add("phi2", phi2);
        exp->step(i)->add("dist", dist);
        exp->step(i)->add("err1", err1);
        exp->step(i)->add("err2", err2);
        exp->step(i)->add("mark", mark);
        exp->save();
        }
#endif
    }
private:
    mesh_ptrtype M_mesh;
    double M_radius;

};

#if defined(USE_BOOST_TEST)
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
BOOST_AUTO_TEST_SUITE( levelset )

// Test wall distance on :
// Unit 2D circle P1G1.
BOOST_AUTO_TEST_CASE( test_2d_p1g1 )
{
    TestLevelSet<2,1,1> tls;
    tls.wallDist_1();
}

// Unit 3D sphere P1G1.
BOOST_AUTO_TEST_CASE( test_3d_p1g1 )
{
    TestLevelSet<3,1,1> tls;
    tls.wallDist_1();
}

// Unit 2D circle P1G2.
BOOST_AUTO_TEST_CASE( test_2d_p1g2 )
{
    TestLevelSet<2,1,2> tls;
    tls.wallDist_1();
}

// Unite 3D sphere P1G2.
BOOST_AUTO_TEST_CASE( test_3D_p1g2 )
{
    TestLevelSet<3,1,2> tls;
    tls.wallDist_1();
}

// Unit 2D circle P1G1
BOOST_AUTO_TEST_CASE( test_2d_circle )
{
    TestLevelSet<2,1,1> tls;
    tls.circle();
}

BOOST_AUTO_TEST_SUITE_END()

#else
int main(int argc, char** argv )
{
    Feel::Environment env( _argc=argc, _argv=argv,
                           _about=makeAbout(),
                           _desc=makeOptions() );
#if 0
    TestLevelSet<2,1,1> tls;
    TestLevelSet<3,1,1> tls;
    TestLevelSet<2,1,1> tls( 10, 0.1 );
    TestLevelSet<3,1,1> tls( 10, 1 );
#endif
    //TestLevelSet<3,1,1> tls( 1, 0.1 );
    //TestLevelSet<2,1,2> tls( 1, 0.1 );
    TestLevelSet<2,1,1> tls;
    tls.wallDist_1();
    tls.wallDist_2();
    tls.circle();

    return 0;
}
#endif
