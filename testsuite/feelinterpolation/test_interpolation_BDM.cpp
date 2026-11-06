/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2026-02-18

   Copyright (C) 2026 Feel++ Consortium

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
#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE interpolation bdm testsuite

#include <feel/feelcore/testsuite.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/print.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/brezzidouglasmarini.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>

namespace Feel
{
inline po::options_description
makeOptions()
{
    po::options_description testHdivInterpolationOptions( "test h_div BDM options" );
    testHdivInterpolationOptions.add_options()
        ( "meshes-2d", po::value< std::vector<std::string> >(), "vector containing mesh names" )
        ( "meshes-3d", po::value< std::vector<std::string> >(), "vector containing mesh names" );
    return testHdivInterpolationOptions.add( Feel::feel_options() );
}

inline AboutData
makeAbout()
{
    AboutData about( "test_interpolation_BDM",
                     "test_interpolation_BDM",
                     "0.1",
                     "Test for interpolation with BDM h_div space",
                     AboutData::License_GPL,
                     "Copyright (c) 2026 Feel++ Consortium" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}

template<int Dim>
class TestInterpolationBDM
    :
    public Application
{
    typedef Application super;

public:
    typedef double value_type;
    typedef Backend<value_type> backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;

    typedef Simplex<Dim,1> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    typedef bases<BrezziDouglasMarini<0>> basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef std::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    TestInterpolationBDM()
        :
        super(),
        M_backend( backend_type::build( soption( _name="backend" ) ) )
    {
        this->changeRepository( boost::format( "%1%/" ) % this->about().appName() );
    }

    void testInterpolationOneElt( std::string one_element_mesh );
    void testInterpolation();

private:
    backend_ptrtype M_backend;
};

template<int Dim>
void
TestInterpolationBDM<Dim>::testInterpolationOneElt( std::string one_element_mesh )
{
    int is3D = ( Dim == 3 ) ? 1 : 0;
    auto myexpr = unitX() + unitY() + is3D * unitZ();

    auto mesh_name = one_element_mesh + ".msh";
    fs::path mesh_path( mesh_name );

    mesh_ptrtype oneelement_mesh = loadMesh( _mesh=new mesh_type,
                                             _filename=mesh_name );

    auto refine_level = std::floor( 1 - math::log( 0.1 ) );
    mesh_ptrtype mesh = loadMesh( _mesh=new mesh_type,
                                  _filename=mesh_name,
                                  _refine=( int )refine_level );

    space_ptrtype Xh = space_type::New( oneelement_mesh );

    std::vector<std::string> faces;
    if ( Dim == 2 )
        faces = { "hypo", "vert", "hor" };
    else if ( Dim == 3 )
        faces = { "xzFace", "xyFace", "xyzFace", "yzFace" };

    element_type U_h_int = Xh->element();
    element_type U_h_on = Xh->element();

    const uint16_type nFacetDof = ( Dim == 2 ) ? Xh->fe()->nDofPerEdge : Xh->fe()->nDofPerFace;
    for ( int f = 0; f < convex_type::numTopologicalFaces; ++f )
    {
        for ( uint16_type l = 0; l < nFacetDof; ++l )
        {
            const uint16_type ldof = static_cast<uint16_type>( f * nFacetDof + l );
            CHECK( ldof < Xh->nLocalDof() );
            CHECK( mesh->hasMarkers( { faces[f] } ) );
            U_h_int( ldof ) = integrate( _range=markedfaces( oneelement_mesh, faces[f] ),
                                         _expr=trans( N() ) * myexpr ).evaluate()( 0, 0 ) / nFacetDof;
        }
    }

    U_h_on.zero();
    U_h_on.on( _range=elements( oneelement_mesh ), _expr=myexpr );

    auto exporter_proj = exporter( _mesh=mesh,
                                   _name=( boost::format( "%1%-%2%" )
                                           % this->about().appName()
                                           % mesh_path.stem().string() ).str() );
    exporter_proj->step( 0 )->add( "U_interpolation_handly-" + mesh_path.stem().string(), U_h_int );
    exporter_proj->step( 0 )->add( "U_interpolation_on-" + mesh_path.stem().string(), U_h_on );
    exporter_proj->save();

    auto error = vf::project( _space=Xh,
                              _range=elements( oneelement_mesh ),
                              _expr=idv( U_h_int ) - idv( U_h_on ) );
    double L2error = error.l2Norm();
    std::cout << "L2 error [BDM one-elt] = " << L2error << std::endl;
}

template<int Dim>
void
TestInterpolationBDM<Dim>::testInterpolation()
{
    int is3D = ( Dim == 3 ) ? 1 : 0;
    auto myexpr = unitX() + unitY() + is3D * unitZ();

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<Dim>> );
    space_ptrtype Xh = space_type::New( mesh );

    auto u_on = Xh->element();
    auto u_proj = Xh->element();

    u_on.on( _range=elements( mesh ), _expr=myexpr );
    u_proj = vf::project( _space=Xh, _range=elements( mesh ), _expr=myexpr );

    auto error_on = vf::project( _space=Xh, _range=elements( mesh ), _expr=myexpr - idv( u_on ) );
    auto error_proj = vf::project( _space=Xh, _range=elements( mesh ), _expr=myexpr - idv( u_proj ) );
    std::cout << "[BDM on] L2 error  = " << error_on.l2Norm() << std::endl;
    std::cout << "[BDM proj] L2 error = " << error_proj.l2Norm() << std::endl;
}

} // namespace Feel

#if USE_BOOST_TEST

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAbout(), Feel::makeOptions() )

BOOST_AUTO_TEST_SUITE( BDM_INTERPOLANT )

BOOST_AUTO_TEST_CASE( test_bdm_interpolant_1 )
{
    using namespace Feel;
    TestInterpolationBDM<2> t2;
    std::vector<std::string> mygeoms2d = vsoption( _name="meshes-2d" );
    for ( std::string const& geo2d : mygeoms2d )
    {
        BOOST_TEST_MESSAGE( "*** BDM interpolant [one-elt 2D] on " << geo2d << " ***" );
        t2.testInterpolationOneElt( geo2d );
    }
    BOOST_TEST_MESSAGE( "*** BDM interpolant [2D] ***" );
    t2.testInterpolation();

    TestInterpolationBDM<3> t3;
    std::vector<std::string> mygeoms3d = vsoption( _name="meshes-3d" );
    for ( std::string const& geo3d : mygeoms3d )
    {
        BOOST_TEST_MESSAGE( "*** BDM interpolant [one-elt 3D] on " << geo3d << " ***" );
        t3.testInterpolationOneElt( geo3d );
    }
    BOOST_TEST_MESSAGE( "*** BDM interpolant [3D] ***" );
    t3.testInterpolation();
}

BOOST_AUTO_TEST_SUITE_END()
#endif
