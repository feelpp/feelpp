//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-   -*- vim:set fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
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
//! @date 25 Feb 2014
//! @copyright 2014-2017 Feel++ Consortium
//!
//!
// tag::all[]
#include <feel/feel.hpp>

int main(int argc, char**argv )
{
    // tag::env[]
    using namespace Feel;
    using Feel::cout;
	po::options_description laplacianoptions( "Laplacian options" );
	laplacianoptions.add_options()
        ( "mu", po::value<double>()->default_value( 1.0 ), "coeff" )
        ( "no-solve", po::value<bool>()->default_value( false ), "No solve" )
		;

	Environment env( _argc=argc, _argv=argv,
                   _desc=laplacianoptions,
                   _about=about(_name="qs_laplacian",
                                _author="Feel++ Consortium",
                                _email="feelpp-devel@feelpp.org"));
    // end::env[]

    tic();
    // tag::mesh[]
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<FEELPP_DIM,1>>);
    // end::mesh[]
    toc("loadMesh");

    tic();
    // tag::discr[]
    auto Vh = Pch<2>( mesh );
    auto u = Vh->element("u");
    auto mu = doption(_name="mu");
    auto f = expr( soption(_name="functions.f"), "f" );
    auto r_1 = expr( soption(_name="functions.a"), "a" ); // Robin left hand side expression
    auto r_2 = expr( soption(_name="functions.b"), "b" ); // Robin right hand side expression
    auto n = expr( soption(_name="functions.c"), "c" ); // Neumann expression
    auto g = expr( soption(_name="functions.g"), "g" );
    // tag::v[]
    auto v = Vh->element( g, "g" );
    // end::v[]
    // end::discr[]
    toc("Vh");



    tic();
    // tag::vf[]
    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),
                  _expr=f*id(v));
    l+=integrate(_range=markedfaces(mesh,"Robin"), _expr=r_2*id(v));
    l+=integrate(_range=markedfaces(mesh,"Neumann"), _expr=n*id(v));
    // tag::vf[]
    toc("l");

    tic();
    // tag::vf[]
    auto a = form2( _trial=Vh, _test=Vh);
    a = integrate(_range=elements(mesh),
                  _expr=mu*gradt(u)*trans(grad(v)) );
    a+=integrate(_range=markedfaces(mesh,"Robin"), _expr=r_1*idt(u)*id(v));
    a+=on(_range=markedfaces(mesh,"Dirichlet"), _rhs=l, _element=u, _expr=g );
    //! if no markers Robin Neumann or Dirichlet are present in the mesh then
    //! impose Dirichlet boundary conditions over the entire boundary
    if ( !mesh->hasAnyMarker({"Robin", "Neumann","Dirichlet"}) )
        a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=g );
    // end::vf[]
    toc("a");

    tic();
    // tag::solve[]
    //! solve the linear system, find u s.t. a(u,v)=l(v) for all v
    if ( !boption( "no-solve" ) )
        a.solve(_rhs=l,_solution=u);
    // end::solve[]
    toc("a.solve");

    // tag::ug[]
    cout << "||u_h-g||_L2=" << normL2(_range=elements(mesh), _expr=idv(u)-g) << std::endl;
    // end::ug[]

    tic();
    // tag::export[]
    auto e = exporter( _mesh=mesh );
    e->addRegions();
    e->add( "u", u );
    e->add( "g", v );
    e->save();
    // end::export[]

    toc("Exporter");
    return 0;

}
// end::all[]
