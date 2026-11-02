/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2020-2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#include "../discretization.hpp"
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/checker.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/traits.hpp>
#include <feel/feeldiscr/check.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelmodels/modelproperties.hpp>
#include <feel/feelpython/pyexpr.hpp>
#include <feel/feelvf/vf_eval.hpp>
#include <feel/feelvf/print.hpp>
#include <feel/feelpde/cg_laplacian.hpp>


namespace Feel 
{
template<int Dim>
int cg_laplacian_model( RuntimeOrder runtimeOrder )
{
    using Feel::cout;
    // tag::mesh_space[]
    tic();
#if FEELPP_HYPERCUBE == 1
    using mesh_t = Mesh<Hypercube<FEELPP_DIM, Dynamic>>;
    auto mesh = loadMesh( _mesh = new mesh_t( quickstart::geometryOrder() ) );
    quickstart::checkMeshOrder( mesh );
#else
    using mesh_t = Mesh<Simplex<Dim, Dynamic>>;
    auto mesh = loadMesh( _mesh = new mesh_t( quickstart::geometryOrder() ) );
    quickstart::checkMeshOrder( mesh );
#endif
    toc( "loadMesh" );

    tic();
    Pch_ptrtype<mesh_t,Dynamic> Vh;
    if ( Environment::vm().count("marker.name") )
        Vh = Pch<Dynamic>( mesh, markedelements(mesh, soption("marker.name") ), runtimeOrder );
    else if ( Environment::vm().count("marker.levelset") )
        Vh = Pch<Dynamic>( mesh, elements(mesh, expr(soption("marker.levelset")) ), runtimeOrder );
    else
        Vh = Pch<Dynamic>( mesh, runtimeOrder );

    // cgLaplacian may not solve the problem, hence u is std::optional
    ModelProperties props;
    props.setupFromFilenameOption();
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

    return 0;
}
}
int main( int argc, char** argv )
{
    // tag::env[]
    using namespace Feel;
    try
    {
        po::options_description laplacianoptions( "Laplacian options" );
        
        laplacianoptions.add_options()
            ( "order", po::value<int>(), "approximation order (alias for case.discretization=Pk)" )
            ( "geometry.order", po::value<int>()->default_value( 1 ), "mesh geometric transformation order" )
            ( "no-solve", po::value<bool>()->default_value( false ), "No solve" )
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

        return cg_laplacian_model<FEELPP_DIM>( quickstart::lagrangeOrder() );
    }
    catch( ... )
    {
        handleExceptions();
    }
}
