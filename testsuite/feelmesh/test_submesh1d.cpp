/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2016-05-17

  Copyright (C) 2016 Feel++ Consortium

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

#define BOOST_TEST_MODULE test_submesh1d

#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/thch.hpp>


FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( submesh1d )

BOOST_AUTO_TEST_CASE( test_submesh1d1 )
{
    using namespace Feel;

    auto mesh3d = loadMesh(_mesh=new Mesh<Simplex<3>>);
    auto Xh3 = THch<1>( mesh3d );
    auto U3 = Xh3->element();
    auto u3 = U3.element<0>();
    auto p3 = U3.element<1>();

    auto mesh1d = createSubmesh(_mesh=mesh3d, _range=markededges(mesh3d,"centerline"));

    size_type nbElements = nelements( elements(mesh1d), true );
    BOOST_CHECK( nbElements > 0 );

}
BOOST_AUTO_TEST_SUITE_END()

