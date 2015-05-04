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

#define USE_BOOST_TEST 1
#if defined(USE_BOOST_TEST)
#define BOOST_TEST_MODULE test_worldcomm
#include <testsuite/testsuite.hpp>
#endif

#include <feel/feel.hpp>

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description wcopts( "Test Environment options" );
    wcopts.add_options()
        ( "Npcomm", Feel::po::value<double>()->default_value( 2 ), "Number of communicators (>=2)" )
    ;
    return wcopts.add( Feel::feel_options() );
}

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "test_worldcomm" ,
                           "test_worldcomm" ,
                           "0.1",
                           "MPI communicators test",
                           Feel::AboutData::License_LGPL,
                           "Copyright (c) 2015 Feel++ Consortium" );

    about.addAuthor( "Guillaume Dolle", "developer", "gdolle@unistra.fr", "" );
    return about;
}

#if defined(USE_BOOST_TEST)
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
BOOST_AUTO_TEST_SUITE( worldcomm )

// Create Npcomm communicators and solve in parallel two laplacians with
// two different boundary strong conditions.
BOOST_AUTO_TEST_CASE( test_0 )
{
    using namespace Feel;
    namespace fs = boost::filesystem;

    Environment env(_argc=boost::unit_test::framework::master_test_suite().argc,
                    _argv=boost::unit_test::framework::master_test_suite().argv,
                    _options=makeOptions(),
                    _about=makeAbout() );

    auto world = Environment::worldComm();
    // Number of communicators to create (Np=k*Ncomm).
    const int Npcomm = 2;
    const int Np = world.size();

    BOOST_CHECK( Npcomm >= 2 );

    if( Np>=Npcomm && !(Np%Npcomm) )
    {
        auto color = world.rank() % Npcomm;

        WorldComm w( color );

        auto mesh = loadMesh( _mesh = new Mesh< Simplex<2> >(w),
                              _filename = "test_twodomains.geo",
                              );

        BOOST_WARN( mesh->numberOfPartitions()!= w.localSize() );

        // Create a feel++ worldcomm using local communicators.
        auto bend = backend( _worldcomm=w.localComm() );

        auto Xh = Pch<1>( mesh );
        auto u = Xh->element();
        auto v = Xh->element();

        // Bilinear form2 create the matrix using the backend worldcomm.
        auto a = form2( _test=Xh, _trial=Xh, _backend=bend );

        a+= integrate( _range=markedelements(mesh,"omega"),
                       _expr=gradt(u)*trans(grad(v) ) );

        // Linear form1 create the vector using the backend worldcomm.
        auto l = form1( _test=Xh, _backend=bend );
        l+= integrate( _range=markedelements(mesh,"omega"),
                       _expr=id(v) );

        a += on( _range=markedfaces( mesh, "wall"), _rhs=l, _element=u, _expr=cst(color), _backend=bend );

        // solveb use the backend worldcomm
        a.solveb( _rhs=l, _solution=u, _backend=bend );

        BOOST_MESSAGE( "nMeshParts: " << mesh->numberOfPartitions()
                       << " | color: " << color
                       << " | god rank: " << w.godRank()
                       << " | global rank: " << w.globalRank()
                       << " | local rank: " << w.localRank()
                       //<< " u size: " << u.size()
                       //<< " u[0]: " << u[0] << std::endl;
                       << "\n"
                     );


        // Exporter use the mesh worldcomm
        //w.showMe();

        auto e = exporter( _mesh=mesh, _name=( boost::format("laplacian_%1%") % color ).str() );
        e->add("u",u);
        e->save();

        auto strcase = (boost::format("laplacian_%1%") % color).str()+".case";
        auto caseFileExist = ( (Environment::findFile( strcase )).empty() != 0 );

        // Check if each communicator has created its ensight .case file.
        if( w.localRank() == 0 )
            BOOST_CHECK( caseFileExist );
    }
}
BOOST_AUTO_TEST_SUITE_END()
#endif
