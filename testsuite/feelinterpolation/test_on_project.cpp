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

#define BOOST_TEST_MODULE test_on
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/vf.hpp>
//#include <feel/feelfilters/exporter.hpp>

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
    Feel::AboutData about( "test_on_project" ,
                           "test_on_project" ,
                           "0.1",
                           "On keyword interpolation test",
                           Feel::AboutData::License_LGPL,
                           "Copyright (c) 2015 Feel++ Consortium" );

    about.addAuthor( "Guillaume Dolle", "developer", "gdolle@unistra.fr", "" );
    return about;
}

template<int DIM, int H_ORDER=1>
class TestOnProject
{
public:
    using mesh_type = Mesh< Simplex<DIM> >;
    using mesh_ptrtype = std::shared_ptr< mesh_type >;

    // Create a test on default cube geometry.
    TestOnProject() :
        M_mesh( loadMesh( _mesh=new mesh_type ) )
    {}

    // Check the on keyword for different topological entities.
    void test()
    {
        auto Xh = Pch<H_ORDER>(M_mesh);

        auto p = Xh->element();
        auto l = Xh->element();
        auto s = Xh->element();
        auto v = Xh->element();

        size_type np = nelements( markedpoints( M_mesh,"P"),true );
        BOOST_CHECK( np > 0 );
        p.on( _range=markedpoints( M_mesh,"P"), _expr=cst(42.) );
        sync( p );
        double sump = p.sum();
        BOOST_CHECK_SMALL( sump - 42, 1e-12 );

        size_type ne = nelements( markededges( M_mesh,"L"),true );
        BOOST_CHECK( ne > 0 );
        l.on( _range=markededges( M_mesh,"L"), _expr=cst(42.) );
        sync( l );
        double maxl = l.max();
        BOOST_CHECK_SMALL( maxl - 42, 1e-12 );

        s.on( _range=markedfaces( M_mesh,"S"), _expr=cst(42.) );
        double ints = integrate(_range=markedfaces( M_mesh,"S"),_expr=idv(s) ).evaluate()(0,0);
        BOOST_CHECK_SMALL( ints - 42, 1e-12 );

        v.on( _range=elements( M_mesh ), _expr=cst(42.) );
        double intv = integrate(_range=elements( M_mesh ),_expr=idv(v) ).evaluate()(0,0);
        BOOST_CHECK_SMALL( intv - 42*27, 1e-10 );

#if 0
        auto e = exporter( _mesh=M_mesh );
        e->add( "p", p);
        e->add( "l", l);
        e->add( "s", s);
        e->add( "v", v);
        e->save();
#endif
    }

private:
    mesh_ptrtype M_mesh;
};


FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )
BOOST_AUTO_TEST_SUITE( on_project )

BOOST_AUTO_TEST_CASE( test_3d )
{
    TestOnProject<3> top;
    top.test();
}

// Not supported yet
//BOOST_AUTO_TEST_CASE( test_2d )
//{
//    TestOnProject<2> top;
//    top.test();
//}

BOOST_AUTO_TEST_SUITE_END()

