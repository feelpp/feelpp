/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#include "../discretization.hpp"
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/table.hpp>
#include <feel/feeldiscr/brokenh.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelpoly/raviartthomas.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelsystem/system.hpp>
#include <feel/feeltiming/tic.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/options.hpp>

#include <iomanip>
#include <sstream>

namespace Feel
{

inline std::string
formatScientific( double value )
{
    std::ostringstream os;
    os << std::scientific << std::setprecision( 6 ) << value;
    return os.str();
}

inline po::options_description
makeOptions()
{
    po::options_description options( "FeelSystem hybrid Darcy formulation options" );
    options.add( backend_options( "sys-darcy-hybrid" ) );
    options.add( backend_options( prefixvm( "sys-darcy-hybrid", "sc" ) ) );
    options.add_options()
        ( "darcy.order", po::value<int>()->default_value( 0 ), "broken RT flux, pressure, and trace order k >= 0" )
        ( "geometry.order", po::value<int>()->default_value( 1 ), "mesh geometric transformation order" )
        ( "darcy.quad", po::value<int>()->default_value( 0 ), "quadrature order (0: automatic FE/geometry estimate)" )
        ( "no-export", po::value<bool>()->default_value( false )->implicit_value( true ), "disable exporter output" )
        ( "rho", po::value<std::string>()->default_value( "1" ), "electric resistivity in rho*j + grad(V) = f" )
        ( "functions.f", po::value<std::string>()->default_value( "{0,0,0}:x:y:z" ), "right hand side in rho*j + grad(V) = f" )
        ( "functions.g", po::value<std::string>()->default_value( "0" ), "source term in div(j) = g" )
        ( "functions.jn", po::value<std::string>()->default_value( "0" ), "normal current j.n on Neumann faces" )
        ( "solution.V", po::value<std::string>()->default_value( "5*2/pi*atan2(y,x):x:y:z" ), "potential used on Dirichlet faces and for diagnostics" )
        ( "solution.j", po::value<std::string>()->default_value( "{10/pi*y/(x*x+y*y),-10/pi*x/(x*x+y*y),0}:x:y:z" ), "exact current used for diagnostics" )
        ( "checks.tolerance", po::value<double>()->default_value( 1e-8 ), "tolerance for conservation and boundary checks" );
    return options;
}

inline AboutData
makeAbout()
{
    AboutData about( "qs_sys_darcy_hybrid",
                     "qs_sys_darcy_hybrid",
                     "0.1",
                     "FeelSystem quickstart 3D hybrid-mixed Darcy/electrostatic problem",
                     AboutData::License_GPL,
                     "Copyright (C) 2026 University of Strasbourg" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}

int
darcyHybridProblem3D()
{
    using namespace Feel::vf;
    namespace role = Feel::systems::roles;

    auto const fluxOrder = quickstart::checkedOrder( "darcy.order", 0 );
    tic();
    auto mesh = loadMesh( _mesh = new Mesh<Simplex<3, Dynamic>>( quickstart::geometryOrder() ) );
    quickstart::checkMeshOrder( mesh );
    double const meshTime = toc( "mesh generation/loading" );

    tic();
    auto Jh = brokenh<RaviartThomas<Dynamic>>( mesh, fluxOrder );
    auto Vh = Pdh<Dynamic>( mesh, fluxOrder );
    auto faceMesh = createSubmesh( _mesh = mesh, _range = faces( mesh ), _update = 0 );
    quickstart::checkMeshOrder( faceMesh );
    auto Mh = Pdh<Dynamic>( faceMesh, fluxOrder );
    auto H = systems::makeHdgSpaces<role::flux_t, role::potential_t, role::trace_t>(
        mesh, faceMesh, product( Jh, Vh, Mh ) );
    double const spaceTime = toc( "space construction" );

    auto dirichletFaces = markedfaces( mesh, "Dirichlet" );
    auto neumannFaces = markedfaces( mesh, "Neumann" );
    size_type const nBoundaryFaces = nelements( boundaryfaces( mesh ), true );
    CHECK_EQ( nelements( dirichletFaces, true ) + nelements( neumannFaces, true ), nBoundaryFaces )
        << "Dirichlet and Neumann markers must form a disjoint partition of the boundary";

    auto rho = expr( soption( "rho" ) );
    auto f = expr<3, 1>( soption( "functions.f" ) );
    auto g = expr( soption( "functions.g" ) );
    auto jn = expr( soption( "functions.jn" ) );
    auto Vexact = expr( soption( "solution.V" ) );
    auto jexact = expr<3, 1>( soption( "solution.j" ) );

    int const requestedQuadratureOrder = ioption( "darcy.quad" );
    if ( requestedQuadratureOrder < 0 )
        throw std::invalid_argument( "darcy.quad must be non-negative" );
    int const quadratureOrder = requestedQuadratureOrder > 0
        ? requestedQuadratureOrder : 2*fluxOrder.value + 4 + 2*( mesh->order() - 1 );

    auto U = H.element();
    auto j = H.flux().trial();
    auto w = H.flux().test();
    auto V = H.potential().trial();
    auto q = H.potential().test();
    auto lambda = H.trace().trial();
    auto mu = H.trace().test();

    tic();
    auto problem =
        hdgProblem( _spaces = H, _name = "sys-darcy-hybrid" )
            .condensation( hdg::condensation()
                               .enabled( boption( "sc.condense" ) )
                               .condenser( condenser_poisson() ) )
            .block( role::flux, role::flux,
                    integrate( _range = elements( mesh ),
                               _expr = rho * inner( j, w ), _quad = quadratureOrder ) )
            .block( role::flux, role::potential,
                    integrate( _range = elements( mesh ),
                               _expr = -V * div( w ), _quad = quadratureOrder ) )
            .block( role::flux, role::trace,
                    integrate( _range = internalfaces( mesh ),
                               _expr = lambda * ( leftface( normal( w ) ) + rightface( normal( w ) ) ), _quad = quadratureOrder ) )
            .block( role::flux, role::trace,
                    integrate( _range = boundaryfaces( mesh ),
                               _expr = lambda * normal( w ), _quad = quadratureOrder ) )
            .block( role::potential, role::flux,
                    integrate( _range = elements( mesh ),
                               _expr = div( j ) * q, _quad = quadratureOrder ) )
            .block( role::potential, role::potential,
                    integrate( _range = elements( mesh ),
                               _expr = cst( 0. ) * V * q, _quad = quadratureOrder ) )
            .block( role::trace, role::flux,
                    integrate( _range = internalfaces( mesh ),
                               _expr = mu * ( leftfacet( normal( j ) ) + rightfacet( normal( j ) ) ), _quad = quadratureOrder ) )
            .block( role::trace, role::flux,
                    integrate( _range = neumannFaces,
                               _expr = mu * normal( j ), _quad = quadratureOrder ) )
            .block( role::trace, role::trace,
                    integrate( _range = dirichletFaces,
                               _expr = lambda * mu, _quad = quadratureOrder ) )
            .row( role::flux,
                  integrate( _range = elements( mesh ),
                             _expr = inner( f, w ), _quad = quadratureOrder ) )
            .row( role::potential,
                  integrate( _range = elements( mesh ),
                             _expr = g * q, _quad = quadratureOrder ) )
            .row( role::trace,
                  integrate( _range = neumannFaces,
                             _expr = mu * jn, _quad = quadratureOrder ) )
            .row( role::trace,
                  integrate( _range = dirichletFaces,
                             _expr = mu * Vexact, _quad = quadratureOrder ) );

    auto system = problem.assemble();
    double const assemblyTime = toc( "problem assembly" );

    tic();
    auto result = system.solve( _solution = U );
    double const solveTime = toc( "linear solve" );
    CHECK( result.converged )
        << "hybrid Darcy solve did not converge; iterations=" << result.iterations
        << ", residual=" << result.residual;

    auto jh = U( H.index( role::flux ) );
    auto VhSol = U( H.index( role::potential ) );
    double const errorJ = normL2( _range = elements( mesh ), _expr = idv( jh ) - jexact, _quad = quadratureOrder );
    double const errorDivJ = normL2( _range = elements( mesh ), _expr = divv( jh ) - div( jexact ), _quad = quadratureOrder );
    double const errorV = normL2( _range = elements( mesh ), _expr = idv( VhSol ) - Vexact, _quad = quadratureOrder );
    double const darcyError = normL2( _range = elements( mesh ),
                                      _expr = rho * idv( jh ) + trans( grad<3>( Vexact ) ) - f, _quad = quadratureOrder );
    double const sideFluxResidual = normL2( _range = neumannFaces,
                                            _expr = inner( idv( jh ), N() ) - jn, _quad = quadratureOrder );
    // Curved-cell Piola divergence is not generally a plain Pdh function.
    // Check conservation in the pressure test space with assembly quadrature.
    auto massResidualForm = form1( _test = Vh );
    massResidualForm = integrate( _range = elements( mesh ),
                                  _expr = ( divv( jh ) - g ) * q, _quad = quadratureOrder );
    massResidualForm.vectorPtr()->close();
    double const massResidual = massResidualForm.vectorPtr()->l2Norm();

    double const tolerance = doption( "checks.tolerance" );
    CHECK_LT( massResidual, tolerance ) << "the weak mass residual must vanish for every pressure test function";
    CHECK_LT( sideFluxResidual, tolerance ) << "j_h.n must satisfy the prescribed Neumann data";

    Table report;
    report.add_row( { "flux/pressure/trace order", std::to_string( fluxOrder.value ) } );
    report.add_row( { "geometry order", std::to_string( mesh->order() ) } );
    report.add_row( { "mesh generation/loading (s)", formatScientific( meshTime ) } );
    report.add_row( { "space construction (s)", formatScientific( spaceTime ) } );
    report.add_row( { "problem assembly (s)", formatScientific( assemblyTime ) } );
    report.add_row( { "solve (s)", formatScientific( solveTime ) } );
    report.add_row( { "J broken global dof", std::to_string( Jh->nDof() ) } );
    report.add_row( { "V local global dof", std::to_string( Vh->nDof() ) } );
    report.add_row( { "lambda face global dof", std::to_string( Mh->nDof() ) } );
    report.add_row( { "hybrid global dof", std::to_string( H.product().nDof() ) } );
    Feel::cout << report << std::endl;

    Table diagnostics;
    diagnostics.add_row( { "||j-j_h||_L2", formatScientific( errorJ ) } );
    diagnostics.add_row( { "||div(j-j_h)||_L2", formatScientific( errorDivJ ) } );
    diagnostics.add_row( { "||V-V_h||_L2", formatScientific( errorV ) } );
    diagnostics.add_row( { "||rho*j_h+grad(V)-f||_L2", formatScientific( darcyError ) } );
    diagnostics.add_row( { "weak mass residual (Euclidean norm)", formatScientific( massResidual ) } );
    diagnostics.add_row( { "||j_h.n-j_N||_L2(Neumann)", formatScientific( sideFluxResidual ) } );
    Feel::cout << diagnostics << std::endl;

    if ( !boption( "no-export" ) )
    {
        auto e = exporter( _mesh = mesh, _name = "qs_sys_darcy_hybrid" );
        e->addRegions();
        std::set<std::string> const representations{ "element", "nodal" };
        e->add( "j", jh, representations );
        e->add( "V", VhSol, representations );
        e->add( "j_exact", jexact );
        e->add( "V_exact", Vexact );
        e->save();
    }

    return result.converged ? 0 : 1;
}

} // namespace Feel

int
main( int argc, char** argv )
{
    using namespace Feel;

    try
    {
        po::options_description options = makeOptions();
        options.add( case_options( 3, "broken-RT0-P0-P0face" ) );

        Environment env( _argc = argc,
                         _argv = argv,
                         _desc = options,
                         _about = makeAbout() );
        return darcyHybridProblem3D();
    }
    catch ( ... )
    {
        handleExceptions();
    }
    return 1;
}
