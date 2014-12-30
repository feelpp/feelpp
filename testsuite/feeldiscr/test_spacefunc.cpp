/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-06-20

  Copyright (C) 2014 Feel++ Consortium

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
/**
   \file test_spacefunc.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-06-20
 */
#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/unitcircle.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feelvf/vf.hpp>

int main( int argc, char **argv)
{
    using namespace Feel;
    Environment env(  _argc=argc, _argv=argv );
    auto mesh=unitCircle();
    auto Xh=Pdhv<0>(mesh);
    auto v = Xh->element();
    CHECK( v.nDof() == 2*mesh->numGlobalElements() );
    v.on( _range=elements(mesh), _expr=ones<2,1>() );
    auto s = v.sum();
    CHECK( s == 2*mesh->numGlobalElements() );
    return 0;
}
