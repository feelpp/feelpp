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
#include <feel/feelcore/checker.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/check.hpp>
#include <feel/feeldiscr/p2ch.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feeldiscr/thch.hpp>
#include <feel/feeldiscr/traits.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelpython/pyexpr.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/print.hpp>
#include <feel/feelvf/vf.hpp>

///[stokes]
template <typename SpaceT>
int stokes( std::shared_ptr<SpaceT>  Vh )
{
    using namespace Feel;

    // tag::mesh_space[]
    tic();
    auto mesh = Vh->mesh();
    auto U = Vh->element( "U" );
    auto u = U( 0_c );
    auto p = U( 1_c );
    auto v = U( 0_c );
    auto q = U( 1_c );
    auto mu = doption( _name = "mu" );
    auto f = expr<FEELPP_DIM, 1>( soption( _name = "functions.f" ) );
    auto thechecker = checker( _name = "qs_stokes", _solution_key = soption( "checker.solution" ) );
    auto solution = expr<FEELPP_DIM, 1>( thechecker.check() ? thechecker.solution() : soption( _name = "functions.g" ) );
    auto g = solution;
    toc( "Vh" );
    // end::mesh_space[]

    // tag::forms[]
    tic();
    auto l = blockform1( _test=Vh );
    l = integrate( _range = elements( mesh ), _expr = inner( f, id( v ) ) );
    toc( "l" );

    tic();
    auto a = blockform2( _test=Vh );
    auto Id = eye<FEELPP_DIM, FEELPP_DIM>();
    auto deft = sym( gradt( u ) );
    auto def = sym( grad( v ) );
    auto sigmat = -idt( p ) * Id + 2 * mu * deft;
    auto sigma = -id( q ) * Id + 2 * mu * def;
    tic();
    a += integrate( _range = elements( mesh ), _expr = inner( sigmat, grad( v ) ) + id( q ) * divt( u ) );
    toc( "a" );

    tic();
    LOG( INFO ) << fmt::format( "adding Dirichlet boundary conditions" );

    auto addDirichlet = [&]( std::vector<std::string> const& markers, auto const& theexpr )
    {
        if ( mesh->hasAnyMarker( markers ) )
        {
            if ( boption( "weak-bc.use" ) )
            {
                double penalty=doption( "weak-bc.penalty" );
                a += integrate( _range = markedfaces( mesh, markers ), _expr = inner( -sigmat * N(), id( v ) ) - inner( sigma * N(), idt( v ) ) + inner( penalty * idt( u ) / hFace(), id( v ) ) );
                l += integrate( _range = markedfaces( mesh, markers ), _expr = inner( -sigma * N() + penalty * id( v ) / hFace(), theexpr ) );
            }
            else
            {
                a.on( _range = markedfaces( mesh, markers ), _test = 0_c, _rhs = l, _element = u, _expr = theexpr );
            }
        }
    };
    auto z = zero<FEELPP_DIM, 1>();
    addDirichlet( std::vector<std::string>{ "inlet", "Dirichlet" }, g );
    addDirichlet( std::vector<std::string>{ "wall", "letters" }, z );

    LOG( INFO ) << fmt::format( "adding Dirichlet boundary conditions done" );

    toc( "a dirichlet" );
    U.vector()->setOnes();
    tic();
    //! solve the linear system, find u s.t. a(u,v)=l(v) for all v
    if ( !boption( "no-solve" ) )
    {
        LOG( INFO ) << fmt::format( "solving linear system" );
        a.solve( _rhs = l, _solution = U );
        LOG( INFO ) << fmt::format( "solving linear system done" );
    }
    toc( "a.solve" );
    // end::forms[]

    // tag::export[]
    tic();
    auto e = exporter( _mesh = mesh );
    e->addRegions();
    e->add( "uh", u );
    e->add( "Uh0", U( 0_c ) );
    e->add( "ph", p );
    e->add( "Uh1", U( 1_c ) );
    if ( thechecker.check() )
    {
        e->add( "u", solution );
    }
    e->save();
    toc( "Exporter" );
    // end::export[]

    return check( thechecker, U( 0_c ) );
}
///[stokes]
int main( int argc, char** argv )
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
            ( "weak-bc.use", po::value<bool>()->default_value( false ), "Use weak boundary conditions (Nitsche)" )
            ( "weak-bc.penalty", po::value<double>()->default_value( 10 ), "Penalty parameter for weak boundary conditions" )
            ;

        Environment env( _argc = argc, _argv = argv,
                         _desc = laplacianoptions,
                         _about = about( _name = "qs_stokes",
                                         _author = "Feel++ Consortium",
                                         _email = "feelpp-devel@feelpp.org" ) );
        ///[stokes-env]

        tic();
        ///[stokes-mesh]
        auto mesh = loadMesh( _mesh = new Mesh<Simplex<FEELPP_DIM, 1>> );
        ///[stokes-mesh]
        toc( "loadMesh" );

        ///[stokes-space]
        int status;
        if ( soption( "space" ) == "P1P0" )
        {
            auto P1hv = Pchv<1>( mesh );
            auto P0h = Pdh<0>( mesh );

            status = stokes( productPtr( P1hv, P0h ) );
        }
        else if ( soption( "space" ) == "P1P1" )
        {
            auto P1hv = Pchv<1>( mesh );
            auto P1h = Pch<1>( mesh );
            status = stokes( productPtr( P1hv, P1h ) );
        }
        else
        {
            // default P2P1: good space
            auto P2hv = Pchv<2>( mesh );
            auto P1h = Pch<1>( mesh );
            status = stokes( productPtr( P2hv, P1h ) );
        }
        ///[stokes-space]
        return status;
    }
    catch ( ... )
    {
        handleExceptions();
    }
    return 1;
}
