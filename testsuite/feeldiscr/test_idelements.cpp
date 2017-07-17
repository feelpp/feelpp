//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 17 Jul 2017
//! @copyright 2017 Feel++ Consortium
//!
#define BOOST_TEST_MODULE test_idelements
#include <testsuite.hpp>

#include <feel/feelfilters/unitcircle.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelmath/randint.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( testsuite_idelements )

BOOST_AUTO_TEST_CASE( listelements )
{
    auto mesh=unitCircle();
    auto l = randint( 3, mesh->beginElement()->id(), boost::prior(mesh->endElement())->id() );
    LOG(INFO) << "l=" << l;
    auto newmesh = createSubmesh( mesh, idelements( mesh, l ) );
    BOOST_CHECK_EQUAL( newmesh->numElements(), l.size()*Environment::numberOfProcessors() );
}

BOOST_AUTO_TEST_CASE( space_fromlistelements )
{
    auto mesh=unitCircle();
    
    auto l = randint( 3, mesh->beginElement()->id(), boost::prior(mesh->endElement())->id() );
    auto newmesh = createSubmesh( mesh, idelements( mesh, l ) );
    auto Xh=Pdh<0>(newmesh);
    auto v = Xh->element(cst(1.));
    auto s = v.sum();
    BOOST_CHECK_EQUAL( s, l.size()*Environment::numberOfProcessors() );
    
}
    
BOOST_AUTO_TEST_SUITE_END()

