//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
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
//! @author <you>
//! @date 15 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!
//!

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>

using namespace Feel;

int main(int argc, char** argv )
{
    Environment env( _argc=argc, _argv=argv );

    auto mesh = loadMesh( new Mesh<Simplex<2> > );
    auto Xh = Pch<1>(mesh);
    auto u = Xh->element();
    auto v = Xh->element();
    auto kappa = doption("parameters.kappa");
    auto gamma = doption("parameters.gamma");
    auto flux = doption("parameters.f");

    auto a = form2(_test=Xh, _trial=Xh);
    a = integrate( markedelements(mesh, "omega1"),
                   inner(gradt(u),grad(v)) );
    a+= integrate( markedelements(mesh, "omega0"),
                   kappa*inner(gradt(u),grad(v)) );

    auto f = form1(_test=Xh);
    f = integrate( markedfaces(mesh, "base"),
                   flux*id(v) );

    a+= on( _range=markedfaces(mesh, "top"), _element=u, _rhs=f, _expr=cst(0.) );

    a.solve( _rhs=f, _solution=u );

    auto e = exporter(_mesh=mesh);
    e->add("u", u);
    e->save();

    return 0;
}
