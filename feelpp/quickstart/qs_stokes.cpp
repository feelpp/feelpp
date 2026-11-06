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
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/checker.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/thch.hpp>
#include <feel/feeldiscr/p2ch.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feeldiscr/check.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelpython/pyexpr.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/vf.hpp>

///[stokes]
template<typename MeshPtrType>
int
stokes(MeshPtrType mesh)
{
    using namespace Feel;

    // tag::mesh_space[]
    tic();
    auto Uh = Pchv<2>(mesh);
    auto Qh = Pch<1>(mesh);
    auto ps = product( Uh, Qh );
    auto u = trial( Uh, "u" );
    auto v = test( Uh, "v" );
    auto p = trial( Qh, "p" );
    auto q = test( Qh, "q" );
    auto W = ps.element();
    auto uh = W( 0_c );
    auto ph = W( 1_c );
    auto mu = doption(_name="mu");
    auto f = expr<FEELPP_DIM,1>( soption(_name="functions.f") );
    auto thechecker = checker(_name="qs_stokes",_solution_key=soption("checker.solution"));
    auto solution = expr<FEELPP_DIM,1>( thechecker.check()? thechecker.solution() : soption(_name="functions.g") );
    auto g = solution;
    toc("Vh");
    // end::mesh_space[]

    // tag::forms[]
    tic();
    auto l = blockform1( ps, solve::strategy::monolithic, backend() );
    l( 0_c ) = integrate(_range=elements(mesh),
                         _expr=inner(f,v));
    toc("l");

    tic();
    auto a = blockform2( ps, solve::strategy::monolithic, backend() );
    a( 0_c, 0_c ) += integrate(_range=elements(mesh),
                               _expr=2*mu*inner( symm_grad(u), symm_grad(v) ) );
    a( 0_c, 1_c ) += integrate(_range=elements(mesh),
                               _expr=-p*div(v) );
    a( 1_c, 0_c ) += integrate(_range=elements(mesh),
                               _expr=q*div(u) );
    toc("a");

    l.close();
    a.close();

    if ( mesh->hasAnyMarker({"inlet","Dirichlet"}) )
        a.row( 0_c ) += on(_range=markedfaces(mesh,{"inlet","Dirichlet"}), _rhs=l( 0_c ), _element=uh, _expr=g, _type="elimination" );
    if ( mesh->hasAnyMarker({"wall","letters"}) )
        a.row( 0_c ) += on(_range=markedfaces(mesh,{"wall","letters"}), _rhs=l( 0_c ), _element=uh, _expr=zero<FEELPP_DIM,1>(), _type="elimination" );
    toc("a");

    tic();
    //! solve the linear system, find u s.t. a(u,v)=l(v) for all v
    if ( !boption( "no-solve" ) )
        a.solve(_rhs=l,_solution=W);
    toc("a.solve");
    // end::forms[]

    // tag::export[]
    tic();
    auto e = exporter( _mesh=mesh );
    e->addRegions();
    e->add( "uh", uh );
    e->add( "ph", ph );
    if ( thechecker.check() )
    {
        e->add( "u", solution );
    }
    e->save();
    toc("Exporter");
    // end::export[]

    return check( thechecker, uh );
}
///[stokes]
int main(int argc, char**argv )
{
    ///[stokes-env]
    using namespace Feel;
    try
    {
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
        ///[stokes-env]

        tic();
        ///[stokes-mesh]
        auto mesh = loadMesh(_mesh=new Mesh<Simplex<FEELPP_DIM,1>>);
        ///[stokes-mesh]
        toc("loadMesh");

        ///[stokes-space]
        return stokes(mesh);
    }
    catch( ... ) 
    {
        handleExceptions();
    }
    return 1;
}
