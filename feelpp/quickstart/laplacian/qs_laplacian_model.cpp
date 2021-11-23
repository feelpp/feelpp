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
//! @date June 15 2020
//! @copyright 2020 Feel++ Consortium
//!
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/checker.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/traits.hpp>
#include <feel/feeldiscr/check.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelmodels/modelproperties.hpp>
#include <feel/feelpython/pyexpr.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/print.hpp>
#include <feel/feelpde/cg_laplacian.hpp>


namespace Feel 
{
template<int Dim, int Order>
int cg_laplacian_model()
{
    using Feel::cout;
    // tag::mesh_space[]
    tic();
#if FEELPP_HYPERCUBE == 1
    using mesh_t = Mesh<Hypercube<FEELPP_DIM, 1>>;
    auto mesh = loadMesh( _mesh = new mesh_t );
#else
    using mesh_t = Mesh<Simplex<Dim, 1>>;
    auto mesh = loadMesh( _mesh = new mesh_t );
#endif
    toc( "loadMesh" );

    tic();
    Pch_ptrtype<mesh_t,Order> Vh;
    if ( Environment::vm().count("marker.name") )
        Vh = Pch<Order>( mesh, markedelements(mesh, soption("marker.name") ) );
    else if ( Environment::vm().count("marker.levelset") )
        Vh = Pch<Order>( mesh, elements(mesh, expr(soption("marker.levelset")) ) );
    else
        Vh = Pch<Order>( mesh );

    // cgLaplacian may not solve the problem, hence u is std::optional
    ModelProperties props( Environment::expand(soption("json_filename")) );
    auto to_export = cgLaplacianModel( Vh, props );

    // tag::export[]
    tic();
    auto e = exporter( _mesh = mesh );
    e->addRegions();
    for( auto const& [name,f]: to_export )
        e->add( name, f );
    e->save();
    toc( "Exporter" );
    // end::export[]

    return true; 
}
}
int main( int argc, char** argv )
{
    // tag::env[]
    using namespace Feel;
    
    po::options_description laplacianoptions( "Laplacian options" );
    
    laplacianoptions.add_options()( "no-solve", po::value<bool>()->default_value( false ), "No solve" )
        ( "json_filename", po::value<std::string>()->default_value( "$cfgdir/model.json" ), "json files" )
        ( "pyexpr.filename", po::value<std::string>()->default_value( "$cfgdir/../python/laplacian.py" ), "python filename to execute" );
    laplacianoptions.add( case_options( FEELPP_DIM, "P1" ) );
    laplacianoptions.add_options()( "marker.name", po::value<std::string>(), "marker on which to solve problem" );
    laplacianoptions.add_options()( "marker.levelset", po::value<std::string>(), "marker on which to solve problem" );

    Environment env( _argc = argc, _argv = argv,
                     _desc = laplacianoptions,
                     _about = about( _name = "qs_laplacian_model",
                                     _author = "Feel++ Consortium",
                                     _email = "feelpp-devel@feelpp.org" ) );

    if ( soption( "case.discretization" ) == "P1" )
        return cg_laplacian_model<FEELPP_DIM,1>();
    
}
