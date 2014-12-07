/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 07 Dec 2014
 
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
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>

using namespace Feel;

struct OpT
{
    typedef Mesh<Simplex<2,1>> mesh_type;
    typedef typename Mesh<Simplex<2,1>>::ptrtype mesh_ptrtype;
    typedef typename Feel::meta::Pch<mesh_type,1>::type space_type;
    typedef typename Feel::meta::Pch<mesh_type,1>::ptrtype space_ptrtype;
    typedef typename Feel::meta::BilinearForm<space_type,space_type>::type form2_type;
    typedef typename Feel::meta::LinearForm<space_type>::type form1_type;

    mesh_ptrtype mesh;
    space_ptrtype Xh;
    form2_type a;
    form1_type l;

    OpT() {}
};


int main(int argc, char** argv)
{

	Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="forms",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    OpT op;
    op.mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    op.Xh = Pch<1>( op.mesh );
    auto u = op.Xh->element("u");
    auto v = op.Xh->element("g");

    op.l = form1( _test=op.Xh );
    op.l = integrate(_range=elements(op.mesh),
                     _expr=id(v));
    
    op.a = form2( _trial=op.Xh, _test=op.Xh);
    op.a = integrate(_range=elements(op.mesh),
                     _expr=gradt(u)*trans(grad(v)) );
    op.a+=on(_range=boundaryfaces(op.mesh), _rhs=op.l, _element=u, _expr=cst(0.) );
    op.a.solve(_rhs=op.l,_solution=u);

    auto e = exporter( _mesh=op.mesh );
    e->addRegions();
    e->add( "u", u );
    e->add( "g", v );
    e->save();
    return 0;
}
