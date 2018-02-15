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
//! @author Vincent Chabannes <vincent.chabannes@feelpp.org>
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 03 May 2017
//! @copyright 2017 Feel++ Consortium
//!
#include <feel/feel.hpp>

template<typename SpacePtrType>
int
stokes(SpacePtrType Vh)
{
    using namespace Feel;
        
    // tag::mesh_space[]
    tic();
    auto mesh = Vh->mesh();
    auto U = Vh->element("U");
    auto u = U.template element<0>();
    auto p = U.template element<1>();
    auto v = U.template element<0>();
    auto q = U.template element<1>();
    auto mu = doption(_name="mu");
    auto f = expr<FEELPP_DIM,1>( soption(_name="functions.f"), "f" );
    auto solution = expr<FEELPP_DIM,1>( checker().solution(), "solution" );
    auto g = checker().check()?solution:expr<FEELPP_DIM,1>( soption(_name="functions.g"), "g" );
    toc("Vh");
    // end::mesh_space[]
    
    // tag::forms[]
    tic();
    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),
                  _expr=inner(f,id(v)));
    toc("l");

    tic();
    auto a = form2( _trial=Vh, _test=Vh);
    auto Id = eye<FEELPP_DIM,FEELPP_DIM>();
    auto deft = sym(gradt(u));
    auto sigmat = -idt(p)*Id + 2*mu*deft;
    tic();
    a = integrate(_range=elements(mesh),
                  _expr=inner( deft, grad(v) ) );
    a += integrate(_range=elements(mesh),
                   _expr=-idt(p)*div(v) );
    a += integrate(_range=elements(mesh),
                   _expr=id(q)*divt(u) );
    toc("a");

    if ( mesh->hasAnyMarker({"inlet","Dirichlet"}) )
        a+=on(_range=markedfaces(mesh,{"inlet","Dirichlet"}), _rhs=l, _element=u, _expr=g );
    if ( mesh->hasAnyMarker({"wall","letters"}) )
        a+=on(_range=markedfaces(mesh,{"wall","letters"}), _rhs=l, _element=u, _expr=zero<FEELPP_DIM,1>() );
    toc("a");

    tic();
    //! solve the linear system, find u s.t. a(u,v)=l(v) for all v
    if ( !boption( "no-solve" ) )
        a.solve(_rhs=l,_solution=U);
    toc("a.solve");
    // end::forms[]
    
    // tag::export[]
    tic();
    auto e = exporter( _mesh=mesh );
    e->addRegions();
    e->add( "uh", u );
    e->add( "ph", p );
    if ( checker().check() )
    {
        v.on(_range=elements(mesh), _expr=solution );
        e->add( "u", v );
    }
    e->save();
    toc("Exporter");
    // end::export[]
    
    // tag::check[]
    // compute l2 and h1 norm of u-u_h where u=solution
    auto norms = [=]( std::string const& solution ) ->std::map<std::string,double>
        {
            tic();
            auto s = expr<FEELPP_DIM,1>(solution);
            double l2 = normL2(_range=elements(mesh), _expr=idv(u)-s );
            toc("L2 error norm");
            tic();
            double h1 = normH1(_range=elements(mesh), _expr=idv(u)-s, _grad_expr=gradv(u)-grad(s) );
            toc("H1 error norm");
            return { { "L2", l2 }, {  "H1", h1 } };
        };
    int status = checker().runOnce( norms, rate::hp( mesh->hMax(), Vh->template functionSpace<0>()->fe()->order() ) );
    // end::check[]

    return status;
}
int main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description laplacianoptions( "Stokes options" );
	laplacianoptions.add_options()
        ( "mu", po::value<double>()->default_value( 1.0 ), "viscosity" )
        ( "space", po::value<std::string>()->default_value( "P2P1" ), "space type: P2P1(default), P1P1, P1P0" )
        ( "no-solve", po::value<bool>()->default_value( false ), "No solve" )
		;

	Environment env( _argc=argc, _argv=argv,
                   _desc=laplacianoptions,
                   _about=about(_name="qs_stokes",
                                _author="Feel++ Consortium",
                                _email="feelpp-devel@feelpp.org"));

    tic();
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<FEELPP_DIM,1>>);
    toc("loadMesh");

    int status;
    if ( soption("space") == "P1P0" )
        status = stokes( P2ch<Lagrange<1,Vectorial>,Lagrange<0,Scalar,Discontinuous>>( mesh ) );
    else if ( soption("space") == "P1P1" )
        status = stokes( P2ch<Lagrange<1,Vectorial>,Lagrange<1,Scalar>>( mesh ) );
    else
    {
        // default P2P1: good space
        status = stokes( THch<1>( mesh ) );
    }
    return !status;
}
