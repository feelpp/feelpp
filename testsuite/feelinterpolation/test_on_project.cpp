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
#include <testsuite/testsuite.hpp>

#include <feel/feel.hpp>

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
    using mesh_ptrtype = boost::shared_ptr< mesh_type >;

    // Create a test on default cube geometry.
    TestOnProject() :
        M_mesh( loadMesh( _mesh=new mesh_type ) )
    {}

    // Check the on keyword for different topological entities.
    void test()
    {
        //auto Xh = Pch<H_ORDER>(M_mesh);
        auto Xh = Pch<H_ORDER>(M_mesh);

        auto p = Xh->element();
        auto l = Xh->element();
        auto s = Xh->element();
        auto v = Xh->element();

        size_type np = nelements( markedpoints( M_mesh,"P") );
        LOG_IF(WARNING, np == 0 ) << "no points marked P:"  << np ;
        //p.on( _range=markedpoints( M_mesh,"P"), _expr=cst(42) );
        if ( Environment::isMasterRank() )
            p.printMatlab("P.m");
        size_type ne = nelements( markededges( M_mesh,"L") );
        LOG_IF(WARNING, ne  == 0 ) << "no edges marked L:"  << ne;
        l.on( _range=markededges( M_mesh,"L"), _expr=cst(42) );
        if ( Environment::isMasterRank() )
            l.printMatlab("L.m");
        s.on( _range=markedfaces( M_mesh,"S"), _expr=cst(42) );
        v.on( _range=elements( M_mesh ), _expr=cst(42) );

        auto e = exporter( _mesh=M_mesh );
        e->add( "p", p);
        e->add( "l", l);
        e->add( "s", s);
        e->add( "v", v);
        e->save();
    }

private:
    mesh_ptrtype M_mesh;
};


FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
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

