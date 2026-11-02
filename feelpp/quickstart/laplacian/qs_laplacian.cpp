/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2013-2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
// tag::global[]
#include "../discretization.hpp"
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/checker.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/traits.hpp>
#include <feel/feeldiscr/check.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelcore/table.hpp>
#include <feel/feelpython/pyexpr.hpp>
#include <feel/feelvf/vf_dsl.hpp>
#include <feel/feelvf/print.hpp>
#include <feel/feelpde/cg_laplacian.hpp>

using namespace Feel;

template<int Dim>
int cg_laplacian_app( RuntimeOrder runtimeOrder )
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
        Vh = Pch<Dynamic>( mesh, elements(mesh, vf::expr(soption("marker.levelset")) ), runtimeOrder );
    else
        Vh = Pch<Dynamic>( mesh, runtimeOrder );

    Feel::cout << fmt::format( "MESH_STATS elements={} points={} dofs={} order={} geometry_order={} mpi_ranks={}\n",
                              mesh->numGlobalElements(), mesh->numGlobalPoints(), Vh->nDof(),
                              runtimeOrder.value, mesh->order(), Environment::worldComm().globalSize() );

    std::map<std::string,std::string> inputs{{"dim",std::to_string(dimension(mesh))},
                                            {"k",soption("k")},{"r_1",soption("r_1")},{"u",""},{"un",soption("un")},{"f",soption("f")},{"g",soption("g")},{"r_2",soption("r_2")}};
    // if we do not check the results with a manufactured solution,
    // the right hand side is given by functions.f otherwise it is computed by the python script
    auto thechecker = checker( _name= "L1/H1 convergence",
                               _solution_key="p",
                               _gradient_key="grad_p",
                               _inputs=inputs
                               );
    auto check_data = [&thechecker]( std::string data, std::string helper, bool defined = true ) {
                          using namespace std::string_literals;
                          if ( defined )
                          {
                              if ( !thechecker.check() && data.empty() )
                                  throw std::invalid_argument( "Invalid data "s + data + " "s + helper);
                          }
                      };
    check_data( soption( _name = "k" ), " set diffusion coefficient using option --k (eg. --k=1)");
    check_data( soption( _name = "f" ), " set source term using option --f (eg. --f=1)");
    check_data( soption( _name = "g" ), " set Dirichlet condition using option --g (eg. --g=0)",  support( Vh )->hasAnyMarker( {"Dirichlet"} ) );
    check_data( soption( _name = "un" ), " set Neumann condition using option --un (eg. --un=0)", support( Vh )->hasAnyMarker( {"Neumann"} ) );
    check_data( soption( _name = "r_2" ), " set Robin right hand side condition using option --r_2 (eg. --r_2=0)", support( Vh )->hasAnyMarker( {"Robin"} ) );

    auto locals = thechecker.runScript();

    std::string p_exact_str = locals.at("p");
    std::string u_exact_str = locals.at("u");
    auto p_exact = vf::expr( p_exact_str );
    auto u_exact = vf::expr<FEELPP_DIM,1>( u_exact_str );
    auto k = vf::expr( locals.at("k") );
    auto un = vf::expr( locals.at("un") );
    auto f = vf::expr( locals.at("f") );
    auto g = vf::expr( locals.at("g") );
    auto r_1 = vf::expr( locals.at("r_1") );
    auto r_2 = vf::expr( locals.at("r_2") );
    Table summary;
    summary.add_row({"Solving -div(( k grad p ) = f with the following boundary conditions"});
    summary(0,0).format().setFontAlign(Font::Align::center);
    Table data;
    //data.format().hide_border();
    data.add_row({"k",locals.at("k")});
    data.add_row({"f",locals.at("f")});

    if ( support( Vh )->hasAnyMarker( {"Dirichlet"} ) )
    {
        data.add_row({"Dirichlet BC",std::to_string(nelements( markedfaces( support( Vh ), "Dirichlet" ) ))});
        data.add_row({"p", locals.at("g")});
    }
    if ( support( Vh )->hasAnyMarker( {"Neumann"} ) )
    {
        data.add_row({"Neumann BC",std::to_string(nelements( markedfaces( support( Vh ), "Neumann" ) ))});
        data.add_row({"-k*dn(p)", locals.at("un")});
    }
    if ( support( Vh )->hasAnyMarker( {"Robin"} ) )
    {
        data.add_row({"Robin BC",std::to_string(nelements( markedfaces( support( Vh ), "Robin" ) ))});
        data.add_row({"-k*dn(p)+"+locals.at("r_1")+"*u", locals.at("r_2")});
    }
    for ( int i=0;i<data.nRow();++i )
        data(i,0).format().setFontAlign(Font::Align::right);
    summary.add_row({data});
    if ( thechecker.check() )
    {
        Table ex;
        //ex.format().hide_border();
        summary.add_row({"Exact solution to be checked with L2/H1 norms"});
        ex.add_row({"p",locals.at("p")});
        if ( locals.count( "grad_p" ) )
            ex.add_row({"grad(p)",locals.at("grad_p")});
        ex.add_row({"-k*grad(p)",locals.at("u")});
        for ( int i=0;i<ex.nRow();++i )
            ex(i,0).format().setFontAlign(Font::Align::right);
        summary.add_row({ex});
    }
    Feel::cout << summary << std::endl;
    if ( Environment::isMasterRank() )
    {
        std::ofstream ofs("model.txt" );
        ofs << summary << std::endl;
    }

    // cgLaplacian may not solve the problem, hence u is std::optional
    auto opt_u = cgLaplacian( _space=Vh,
                              _data=std::tuple{k,f,g,un,r_1,r_2},
                              _solve=!boption( "no-solve" ),
                              _weak_dirichlet=boption( "weakbc" ) );

    if ( !boption( "no-export" ) )
    {
        // tag::export[]
        tic();
        auto e = exporter( _mesh = mesh );
        e->addRegions();
        if ( opt_u )
        {
            e->add( "p", *opt_u );
            e->add( "u", -k * vf::gradv( *opt_u ), elements( support( Vh ) ), "element" );
        }
        e->add( "k", k );
        e->add( "f", f );
        e->add( "g", g );
        if ( support( Vh )->hasAnyMarker( {"Robin"} ) )
        {
            auto rangeFacesRobin = markedfaces( support( Vh ), "Robin" );
            e->add( "r_1", r_1, rangeFacesRobin );
            e->add( "r_2", r_2, rangeFacesRobin );
        }
        if ( thechecker.check() )
        {
            e->add( "solution", p_exact );
            e->add( "flux", u_exact );
        }
        e->save();
        toc( "Exporter" );
        // end::export[]
    }

    if ( opt_u )
        return check( thechecker, *opt_u );
    return 0;
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
            ( "no-export", po::value<bool>()->default_value( false )->implicit_value( true ), "No export" )
            ( "weakbc", po::value<bool>()->default_value( false ), "Weak Dirichlet conditions" )
            ( "k", po::value<std::string>()->default_value( "1" ), "diffusion coefficient" )
            ( "f", po::value<std::string>()->default_value( "" ), "right hand side" )
            ( "g", po::value<std::string>()->default_value( "" ), "Dirichlet boundary condition" )
            ( "un", po::value<std::string>()->default_value( "" ), "Neumann boundary condition" )
            ( "r_1", po::value<std::string>()->default_value( "1" ), "Robin left hand side coefficient" )
            ( "r_2", po::value<std::string>()->default_value( "" ), "Robin right hand side  coefficient" )
            ( "pyexpr.filename", po::value<std::string>()->default_value( "$cfgdir/../python/laplacian.py" ), "python filename to execute" );
        laplacianoptions.add( case_options( FEELPP_DIM, "P1" ) );
        laplacianoptions.add_options()( "marker.name", po::value<std::string>(), "marker on which to solve problem" );
        laplacianoptions.add_options()( "marker.levelset", po::value<std::string>(), "marker on which to solve problem" );

        Environment env( _argc = argc, _argv = argv,
                         _desc = laplacianoptions,
                         _about = about( _name = "qs_laplacian",
                                         _author = "Feel++ Consortium",
                                         _email = "feelpp-devel@feelpp.org" ) );

        // end::env[]
        return cg_laplacian_app<FEELPP_DIM>( quickstart::lagrangeOrder() );
    }
    catch(...)
    {
        handleExceptions();
    }
    return 1;
}
// end::global[]
