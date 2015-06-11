/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-06-20

  Copyright (C) 2014-2015 Feel++ Consortium

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
    CHECK( s == 2*mesh->numGlobalElements() ) << "invalid result " << s << " should be " << 2*mesh->numGlobalElements();

    CHECK( v[Component::X].size() == mesh->numGlobalElements() )
        << "invalid result " << v[Component::Y].size() << " should be " << mesh->numGlobalElements();
    CHECK( v[Component::Y].size() == mesh->numGlobalElements() )
        << "invalid result " << v[Component::Y].size() << " should be " << mesh->numGlobalElements();
    
    v[Component::X].on( _range=elements(mesh), _expr=cst(3.) );
    v[Component::Y].on( _range=elements(mesh), _expr=cst(0) );
    s = v[Component::X].sum();
    CHECK( s == 3*mesh->numGlobalElements() ) << "Invalid sum result = " << s << " should be " << 3*mesh->numGlobalElements();
    v[Component::X].on( _range=elements(mesh), _expr=Px() );
    v[Component::Y].on( _range=elements(mesh), _expr=Py() );
    auto r = normL2( _range=elements(mesh), _expr=idv(v)-vec(Px(),Py()));
    CHECK( r < 1e-12 ) << "invalid L2 error norm " << r << " should be 0";
    v.on(_range=elements(mesh),_expr=vec(3.1*Py(),-2*Px()) );
    r = normL2( _range=elements(mesh), _expr=idv(v[Component::X])-3.1*Py());
    CHECK( r < 1e-12 ) << "invalid L2 error norm " << r << " should be 0";
    r = normL2( _range=elements(mesh), _expr=idv(v[Component::Y])+2*Px());
    CHECK( r < 1e-12 ) << "invalid L2 error norm " << r << " should be 0";
    
    v[Component::X].on( _range=elements(mesh), _expr=cst(0) );
    v[Component::Y].on( _range=elements(mesh), _expr=cst(4.) );
    s = v[Component::Y].sum();
    CHECK( s == 4*mesh->numGlobalElements() ) << "Invalid sum result = " << s << " should be " << 4*mesh->numGlobalElements();
    return 0;
}
