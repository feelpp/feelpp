/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2013-07-05

   Copyright (C) 2013 Université de Strasbourg

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

#define BOOST_TEST_MODULE test_meshmarker
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( mesh_transfinite )

BOOST_AUTO_TEST_CASE( test_load )
{
    using namespace Feel;

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<3>> );
    //auto Xh = Pch<3>( mesh );
    auto area = integrate( _range=boundaryfaces(mesh), _expr=cst(1.) ).evaluate()(0,0);
    auto bottom = integrate( _range=markedfaces(mesh,"Bottom"), _expr=cst(1.) ).evaluate()(0,0);
    auto top = integrate( _range=markedfaces(mesh,"Top"), _expr=cst(1.) ).evaluate()(0,0);
    auto right = integrate( _range=markedfaces(mesh,"Right"), _expr=cst(1.) ).evaluate()(0,0);
    auto left = integrate( _range=markedfaces(mesh,"Left"), _expr=cst(1.) ).evaluate()(0,0);
    auto front = integrate( _range=markedfaces(mesh,"Front"), _expr=cst(1.) ).evaluate()(0,0);
    auto back = integrate( _range=markedfaces(mesh,"Back"), _expr=cst(1.) ).evaluate()(0,0);

    BOOST_CHECK_CLOSE( top, 100, 1e-10 );
    BOOST_CHECK_CLOSE( left, 100, 1e-10 );
    BOOST_CHECK_CLOSE( front, 100, 1e-10 );
    BOOST_CHECK_CLOSE( back, 100, 1e-10 );
    BOOST_CHECK_CLOSE( right, 100, 1e-10 );
    BOOST_CHECK_CLOSE( bottom, 100, 1e-10 );
}
BOOST_AUTO_TEST_SUITE_END()
