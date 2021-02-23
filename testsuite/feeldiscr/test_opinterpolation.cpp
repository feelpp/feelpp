/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Guillaume Dollé <gdolle@unistra.fr>
 Date: 26 Mar 2015

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

#if defined(USE_BOOST_TEST)
#define BOOST_TEST_MODULE test_opinterpolation
#include <feel/feelcore/testsuite.hpp>
#endif

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/unithypercube.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description opts( "Test Environment options" );
    opts.add_options()
        ( "d", Feel::po::value<double>()->default_value( 1 ), "Value" )
    ;
    return opts.add( Feel::feel_options() );
}

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "test_opinterpolation" ,
                           "test_opinterpolation" ,
                           "0.1",
                           "On keyword interpolation test",
                           Feel::AboutData::License_LGPL,
                           "Copyright (c) 2015 Feel++ Consortium" );

    about.addAuthor( "Guillaume Dolle", "developer", "gdolle@unistra.fr", "" );
    return about;
}

template<int DIM, int H_ORDER=1>
class TestOpInterpolation
{
public:
    using mesh_type = Mesh< Simplex<DIM> >;
    using mesh_ptrtype = std::shared_ptr< mesh_type >;

    // Create a test on default cube geometry.
    TestOpInterpolation() :
        M_mesh( loadMesh( _mesh=new mesh_type ) )
    {}

    /// Test the opInterpolator on edges
    /// This test create a 1D submesh, project a constant, then interpolate
    /// function on the 2D/3D domain.
    void testEdges()
    {
        auto mesh1d = createSubmesh( _mesh=M_mesh, _range=markededges( M_mesh, "L" ) );

        auto Xh = Pch<H_ORDER>(M_mesh);
        auto Xh1d = Pch<H_ORDER>(mesh1d);

        auto op1d = opInterpolation( _domainSpace=Xh, _imageSpace=Xh1d );
        auto op1dT = op1d->adjoint( MATRIX_TRANSPOSE_UNASSEMBLED );

        auto l = Xh->element();
        auto l1d = Xh1d->element();

        l.on( _range=markededges( M_mesh, "L" ), _expr=cst(42) );
        l1d.on( _range=elements( mesh1d ), _expr=cst(42) );

        auto l1d_interp = op1dT->operator()( l1d );

#if !defined(USE_BOOST_TEST)
        auto e = exporter( _mesh=M_mesh, _name="edges" );
        e->add( "l", l );
        e->add( "l1d_interp", l1d_interp );
        e->save();
#endif

        auto l1d_l2error = normL2( elements(M_mesh), idv(l)-idv(l1d_interp) );

        VLOG(3) << "L2 error on line: " << l1d_l2error;

#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( l1d_l2error, 1e-14 );
#endif
    }

    /// Test the opInterpolator on surfaces
    /// This test create a 2D submesh, project a constant, then interpolate
    /// function on the 2D/3D domain.
    void testSurfaces()
    {
        auto mesh2d = createSubmesh( _mesh=M_mesh, _range=markedfaces(M_mesh,"S") );

        auto Xh = Pch<H_ORDER>(M_mesh);
        auto Xh2d = Pch<H_ORDER>(mesh2d);

        auto op2d = opInterpolation( _domainSpace=Xh, _imageSpace=Xh2d );
        auto op2dT = op2d->adjoint(MATRIX_TRANSPOSE_UNASSEMBLED);

        auto s = Xh->element();
        auto s2d = Xh2d->element();

        s.on( _range=markefaces( M_mesh,"S"), _expr=cst(42) );
        s2d.on( _range=elements( mesh2d ), _expr=cst(42) );

        auto s2d_interp = op2dT->operator()( s2d );

#if !defined(USE_BOOST_TEST)
        auto e = exporter( _mesh=M_mesh, _name="surfaces" );
        e->add( "s", s);
        e->add( "s2d_interp", s2d_interp );
        e->save();
#endif

        auto s2d_l2error = normL2( elements(M_mesh), idv(s)-idv(s2d_interp) );

        VLOG(3) << "L2 error on surface: " << s2d_l2error;

#if defined(USE_BOOST_TEST)
        BOOST_CHECK_SMALL( s2d_l2error, 1e-14 );
#endif
    }

private:
    mesh_ptrtype M_mesh;
};

#if defined(USE_BOOST_TEST)
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
BOOST_AUTO_TEST_SUITE( opinterpolation )

BOOST_AUTO_TEST_CASE( test_edges )
{
    TestOpInterpolation<3> top;
    top.testEdges();
}

BOOST_AUTO_TEST_CASE( test_surfaces )
{
    TestOpInterpolation<3> top;
    top.testSurfaces();
}

BOOST_AUTO_TEST_SUITE_END()

#else
int main(int argc, char** argv )
{
    Feel::Environment env( _argc=argc, _argv=argv,
                           _about=makeAbout(),
                           _desc=makeOptions() );

    TestOpInterpolation<3> top;
    top.testEdges();

    return 0;
}
#endif
