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
//! @date 10 Apr 2017
//! @copyright 2017 Feel++ Consortium
//!
// tag::global[]
#include <feel/feel.hpp>

int main(int argc, char**argv )
{
    // tag::env[]
    using namespace Feel;
    using Feel::cout;
	po::options_description laplacianoptions( "Laplacian options" );
	laplacianoptions.add_options()
        ( "no-solve", po::value<bool>()->default_value( false ), "No solve" )
		;

	Environment env( _argc=argc, _argv=argv,
                   _desc=laplacianoptions,
                   _about=about(_name="qs_laplacian",
                                _author="Feel++ Consortium",
                                _email="feelpp-devel@feelpp.org"));
    // end::env[]

    // tag::mesh_space[]
    tic();
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<FEELPP_DIM,1>>);
    toc("loadMesh");

    tic();
    auto Vh = Pch<2>( mesh );
    auto u = Vh->element("u");
    auto mu = expr(soption(_name="functions.mu")); // diffusion term
    auto f = expr( soption(_name="functions.f"), "f" );
    auto r_1 = expr( soption(_name="functions.a"), "a" ); // Robin left hand side expression
    auto r_2 = expr( soption(_name="functions.b"), "b" ); // Robin right hand side expression
    auto n = expr( soption(_name="functions.c"), "c" ); // Neumann expression
    auto solution = expr( checker().solution(), "solution" );
    auto g = checker().check()?solution:expr( soption(_name="functions.g"), "g" );
    // tag::v[]
    auto v = Vh->element( g, "g" );
    // end::v[]
    toc("Vh");
    // end::mesh_space[]

    // tag::forms[]
    tic();
    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),
                  _expr=f*id(v));
    l+=integrate(_range=markedfaces(mesh,"Robin"), _expr=r_2*id(v));
    l+=integrate(_range=markedfaces(mesh,"Neumann"), _expr=n*id(v));
    toc("l");

    tic();
    auto a = form2( _trial=Vh, _test=Vh);
    tic();
    a = integrate(_range=elements(mesh),
                  _expr=mu*inner(gradt(u),grad(v)) );
    toc("a");
    a+=integrate(_range=markedfaces(mesh,"Robin"), _expr=r_1*idt(u)*id(v));
    a+=on(_range=markedfaces(mesh,"Dirichlet"), _rhs=l, _element=u, _expr=g );
    //! if no markers Robin Neumann or Dirichlet are present in the mesh then
    //! impose Dirichlet boundary conditions over the entire boundary
    if ( !mesh->hasAnyMarker({"Robin", "Neumann","Dirichlet"}) )
        a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=g );
    toc("a");
    // end::forms[]

    // tag::solve[]
    tic();
    //! solve the linear system, find u s.t. a(u,v)=l(v) for all v
    if ( !boption( "no-solve" ) )
        a.solve(_rhs=l,_solution=u);
    toc("a.solve");
    // end::solve[]

    // tag::export[]
    tic();
    auto e = exporter( _mesh=mesh );
    e->addRegions();
    e->add( "uh", u );
    if ( checker().check() )
    {
        v.on(_range=elements(mesh), _expr=solution );
        e->add( "solution", v );
    }
    e->save();
    toc("Exporter");
    // end::export[]

    // tag::check[]
    // compute l2 and h1 norm of u-u_h where u=solution
    auto norms = [=]( std::string const& solution ) ->std::map<std::string,double>
        {
            tic();
            double l2 = normL2(_range=elements(mesh), _expr=idv(u)-expr(solution) );
            toc("L2 error norm");
            tic();
            double h1 = normH1(_range=elements(mesh), _expr=idv(u)-expr(solution), _grad_expr=gradv(u)-grad<FEELPP_DIM>(expr(solution)) );
            toc("H1 error norm");
            return { { "L2", l2 }, {  "H1", h1 } };
        };

    int status = checker().runOnce( norms, rate::hp( mesh->hMax(), Vh->fe()->order() ) );
    // end::check[]

    // exit status = 0 means no error
    return !status;

}
// end::global[]
