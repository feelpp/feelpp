/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 31 Jan 2016
 
 Copyright (C) 2016 Feel++ Consortium
 
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
#define BOOST_TEST_MODULE test_topetsc
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelalg/topetsc.hpp>
#include <feel/feeldiscr/stencil.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>

/** use Feel namespace */
using namespace Feel;

inline
po::options_description makeOptions()
{
    po::options_description options( "Test toPETSc  Options" );
    return options;
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_topetsc" ,
                     "test_topetsc" ,
                     "0.2",
                     "test topetsc",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2016 Feel++ Consortium" );

    about.addAuthor( "C Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}



FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
BOOST_AUTO_TEST_SUITE( topetsc_suite )

BOOST_AUTO_TEST_CASE( test_backend )
{
    BOOST_MESSAGE( "test_backend");
    auto b = backend();
    auto bp = toPETSc( b );
    BOOST_CHECK( bp );
    BOOST_MESSAGE( "test_backend done");
}

BOOST_AUTO_TEST_CASE( test_backend_not_petsc )
{
    BOOST_MESSAGE( "test_backend_not_petsc");
    auto b = backend(_kind="eigen");
    auto bp = toPETSc( b );
    BOOST_CHECK( bp == 0 );
    BOOST_MESSAGE( "test_backend_not_petsc done");
}

BOOST_AUTO_TEST_CASE( test_vector )
{
    BOOST_MESSAGE( "test_vector");
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Xh = Pch<1>( mesh );
    auto b = backend();
    auto bp = toPETSc( b );
    BOOST_CHECK( bp );
    auto v = b->newVector( Xh );
    auto vp = toPETSc( v );
    BOOST_CHECK( vp );
    Vec vvec = vp->vec();
    BOOST_CHECK( vvec );
    BOOST_MESSAGE( "test_vector done");
}

BOOST_AUTO_TEST_CASE( test_matrix )
{
    BOOST_MESSAGE( "test_matrix");
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Xh = Pch<1>( mesh );
    auto b = backend();
    auto bp = toPETSc( b );
    BOOST_CHECK( bp );
    auto v = b->newMatrix( _test=Xh, _trial=Xh );
    auto vp = toPETSc( v );
    BOOST_CHECK( vp );
    Mat pm = vp->mat();
    BOOST_CHECK( pm );

    BOOST_MESSAGE( "test_matrix done");
}

BOOST_AUTO_TEST_CASE( test_preconditioner )
{
    BOOST_MESSAGE( "test_preconditioner");
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Xh = Pch<1>( mesh );
    auto b = backend();
    auto bp = toPETSc( b );
    BOOST_CHECK( bp );
    auto v = b->newMatrix( _trial=Xh, _test=Xh );
    auto vp = toPETSc( v );
    BOOST_CHECK( vp );
    Mat pm = vp->mat();
    BOOST_CHECK( pm );
    auto p = b->preconditioner();
    BOOST_CHECK( p );
    auto pp = toPETSc( p );
    BOOST_CHECK( pp );
    
    BOOST_MESSAGE( "test_preconditioner done");
}



BOOST_AUTO_TEST_SUITE_END()

