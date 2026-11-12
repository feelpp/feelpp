/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/table.hpp>
#include <feel/feeldiscr/brokenh.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelpoly/brezzidouglasmarini.hpp>
#include <feel/feelpoly/raviartthomas.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeltiming/tic.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/vf.hpp>

#include <iomanip>
#include <sstream>

namespace Feel
{
struct DarcySpaceDescriptor
{
    std::string fluxFamily;
    int fluxOrder;
    int pressureOrder;
    int normalTraceOrder;
    int divergenceOrder;
};

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
    po::options_description options( "Hybrid Darcy formulation options" );
    options.add_options()
        ( "darcy.flux-family", po::value<std::string>()->default_value( "rt" ), "H(div) flux family: rt or bdm" )
        ( "darcy.order", po::value<int>()->default_value( 0 ), "mathematical flux order (RT k >= 0, BDM k >= 1)" )
        ( "geometry.order", po::value<int>()->default_value( 1 ), "mesh geometric transformation order" )
        ( "darcy.quad", po::value<int>()->default_value( 0 ), "quadrature order (0: automatic FE/geometry estimate)" )
        ( "darcy.convergence-case", po::value<bool>()->default_value( false )->implicit_value( true ), "use the built-in smooth manufactured convergence case" )
        ( "no-export", po::value<bool>()->default_value( false )->implicit_value( true ), "disable exporter output" )
        ( "rho", po::value<std::string>()->default_value( "1" ), "electric resistivity in rho*j + grad(V) = f" )
        ( "functions.jn", po::value<std::string>()->default_value( "0" ), "normal current j.n on Neumann faces" )
        ( "solution.V", po::value<std::string>()->default_value( "5*2/pi*atan2(y,x):x:y:z" ), "potential used on Dirichlet faces and for diagnostics" )
        ( "solution.j", po::value<std::string>()->default_value( "{10/pi*y/(x*x+y*y),-10/pi*x/(x*x+y*y),0}:x:y:z" ), "exact current used for diagnostics" )
        ( "checks.tolerance", po::value<double>()->default_value( 1e-8 ), "tolerance for conservation and boundary checks" );
    return options;
}

inline AboutData
makeAbout()
{
    AboutData about( "qs_darcy_hybrid",
                     "qs_darcy_hybrid",
                     "0.1",
                     "Quickstart 3D hybrid-mixed Darcy/electrostatic problem",
                     AboutData::License_GPL,
                     "Copyright (C) 2026 University of Strasbourg" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}

template<typename MeshType, typename FluxSpacePtr, typename PressureSpacePtr, typename TraceSpacePtr>
int
darcyHybrid3D( std::shared_ptr<MeshType> const& mesh,
               FluxSpacePtr const& Jh,
               PressureSpacePtr const& Vh,
               TraceSpacePtr const& Mh,
               DarcySpaceDescriptor const& spaceDescriptor,
               double meshTime )
{
    using namespace boost::hana::literals;

    tic();
    auto Xh = product( Jh, Vh, Mh );
    double const spaceTime = toc( "space construction" );

    auto dirichletFaces = markedfaces( mesh, "Dirichlet" );
    auto neumannFaces = markedfaces( mesh, "Neumann" );
    size_type const nBoundaryFaces = nelements( boundaryfaces( mesh ), true );
    CHECK_EQ( nelements( dirichletFaces, true ) + nelements( neumannFaces, true ), nBoundaryFaces )
        << "Dirichlet and Neumann markers must form a disjoint partition of the boundary";

    auto rho = expr( soption( "rho" ) );
    bool const useConvergenceCase = boption( "darcy.convergence-case" );
    // functions.f/g are registered by feel_options(); keep the vector default
    // here without creating ambiguous duplicate command-line options.
    std::string const fExpression = useConvergenceCase || Environment::vm().at( "functions.f" ).defaulted()
                                       ? "{0,0,0}:x:y:z" : soption( "functions.f" );
    std::string const gExpression = useConvergenceCase ? "-30*(x^4+y^4+z^4):x:y:z" : soption( "functions.g" );
    std::string const potentialExpression = useConvergenceCase ? "x^6+y^6+z^6:x:y:z" : soption( "solution.V" );
    std::string const fluxExpression = useConvergenceCase
                                           ? "{-6*x^5,-6*y^5,-6*z^5}:x:y:z"
                                           : soption( "solution.j" );
    auto f = expr<3, 1>( fExpression );
    auto g = expr( gExpression );
    auto jn = expr( soption( "functions.jn" ) );
    auto Vexact = expr( potentialExpression );
    auto jexact = expr<3, 1>( fluxExpression );

    auto U = Xh.element();
    auto T = Xh.element();
    auto j = U( 0_c );
    auto V = U( 1_c );
    auto lambda = U( 2_c );
    auto w = T( 0_c );
    auto q = T( 1_c );
    auto mu = T( 2_c );

    // Curved Piola mass terms are generally nonpolynomial: this estimate
    // increases integration accuracy with geometry order, without claiming
    // exactness. darcy.quad permits quadrature-refinement checks.
    int const automaticQuadratureOrder =
        2*spaceDescriptor.fluxOrder + 4 + 2*( mesh->order() - 1 );
    int const requestedQuadratureOrder = ioption( "darcy.quad" );
    CHECK_GE( requestedQuadratureOrder, 0 ) << "darcy.quad must be non-negative";
    int const quadratureOrder = requestedQuadratureOrder > 0
                                    ? requestedQuadratureOrder
                                    : ( useConvergenceCase
                                            ? std::max( automaticQuadratureOrder, 12 )
                                            : automaticQuadratureOrder );

    tic();
    backend( _rebuild = true );
    solve::strategy const strategy = boption( "sc.condense" )
                                        ? solve::strategy::static_condensation
                                        : solve::strategy::monolithic;
    auto a = blockform2( Xh, strategy, backend(), Pattern::EXTENDED );
    auto l = blockform1( Xh, strategy, backend() );

    // Local compatible mixed equations. The face potential lambda couples cells.
    l( 0_c ) += integrate( _range = elements( mesh ), _expr = inner( f, id( w ) ), _quad = quadratureOrder );
    l( 1_c ) += integrate( _range = elements( mesh ), _expr = g * id( q ), _quad = quadratureOrder );
    l( 2_c ) += integrate( _range = neumannFaces, _expr = id( mu ) * jn, _quad = quadratureOrder );
    l( 2_c ) += integrate( _range = dirichletFaces, _expr = id( mu ) * Vexact, _quad = quadratureOrder );

    a( 0_c, 0_c ) += integrate( _range = elements( mesh ), _expr = rho * inner( idt( j ), id( w ) ), _quad = quadratureOrder );
    a( 0_c, 1_c ) += integrate( _range = elements( mesh ), _expr = -idt( V ) * div( w ), _quad = quadratureOrder );
    a( 0_c, 2_c ) += integrate( _range = internalfaces( mesh ),
                                 _expr = idt( lambda ) * ( leftface( normal( w ) ) + rightface( normal( w ) ) ),
                                 _quad = quadratureOrder );
    a( 0_c, 2_c ) += integrate( _range = boundaryfaces( mesh ), _expr = idt( lambda ) * normal( w ), _quad = quadratureOrder );

    a( 1_c, 0_c ) += integrate( _range = elements( mesh ), _expr = divt( j ) * id( q ), _quad = quadratureOrder );
    // The condenser expects every local block to be registered, including this exact zero.
    a( 1_c, 1_c ) += integrate( _range = elements( mesh ), _expr = cst( 0. ) * idt( V ) * id( q ), _quad = quadratureOrder );

    a( 2_c, 0_c ) += integrate( _range = internalfaces( mesh ),
                                 _expr = id( mu ) * ( leftfacet( normalt( j ) ) + rightfacet( normalt( j ) ) ),
                                 _quad = quadratureOrder );
    a( 2_c, 0_c ) += integrate( _range = neumannFaces, _expr = id( mu ) * normalt( j ), _quad = quadratureOrder );
    a( 2_c, 2_c ) += integrate( _range = dirichletFaces, _expr = idt( lambda ) * id( mu ), _quad = quadratureOrder );
    double const assemblyTime = toc( "space assembly" );

    tic();
    auto solveStatus = a.solve( _rhs = l, _solution = U, _condense = boption( "sc.condense" ) );
    double const solveTime = toc( "linear solve" );
    CHECK( solveStatus.isConverged() )
        << "hybrid Darcy solve did not converge; iterations=" << solveStatus.nIterations()
        << ", residual=" << solveStatus.residual();

    auto jh = U( 0_c );
    auto VhSol = U( 1_c );
    double const errorJ = normL2( _range = elements( mesh ), _expr = idv( jh ) - jexact, _quad = quadratureOrder );
    double const errorDivJ = normL2( _range = elements( mesh ), _expr = divv( jh ) - g, _quad = quadratureOrder );
    double const errorV = normL2( _range = elements( mesh ), _expr = idv( VhSol ) - Vexact, _quad = quadratureOrder );
    double const darcyError = normL2( _range = elements( mesh ),
                                      _expr = rho * idv( jh ) + trans( grad<3>( Vexact ) ) - f,
                                      _quad = quadratureOrder );
    double const sideFluxResidual = normL2( _range = neumannFaces,
                                            _expr = inner( idv( jh ), N() ) - jn,
                                            _quad = quadratureOrder );
    // Test conservation in the pressure space using the assembly quadrature.
    // On curved cells, Piola divergence need not itself lie in plain Pdh.
    auto massResidualForm = form1( _test = Vh );
    massResidualForm = integrate( _range = elements( mesh ),
                                  _expr = ( divv( jh ) - g ) * id( q ),
                                  _quad = quadratureOrder );
    massResidualForm.vectorPtr()->close();
    double const massResidual = massResidualForm.vectorPtr()->l2Norm();

    double const tolerance = doption( "checks.tolerance" );
    CHECK_LT( massResidual, tolerance ) << "the weak mass residual must vanish for every pressure test function";
    CHECK_LT( sideFluxResidual, tolerance ) << "j_h.n must satisfy the prescribed Neumann data";

    Table report;
    report.add_row( { "mesh generation/loading (s)", formatScientific( meshTime ) } );
    report.add_row( { "method", "hybrid" } );
    report.add_row( { "solver strategy", boption( "sc.condense" ) ? "static condensation" : "monolithic" } );
    report.add_row( { "flux family", spaceDescriptor.fluxFamily } );
    report.add_row( { "flux order", std::to_string( spaceDescriptor.fluxOrder ) } );
    report.add_row( { "pressure order", std::to_string( spaceDescriptor.pressureOrder ) } );
    report.add_row( { "normal-trace order", std::to_string( spaceDescriptor.normalTraceOrder ) } );
    report.add_row( { "divergence-image order", std::to_string( spaceDescriptor.divergenceOrder ) } );
    report.add_row( { "space construction (s)", formatScientific( spaceTime ) } );
    report.add_row( { "assembly (s)", formatScientific( assemblyTime ) } );
    report.add_row( { "solve (s)", formatScientific( solveTime ) } );
    report.add_row( { "J broken global dof", std::to_string( Jh->nDof() ) } );
    report.add_row( { "V local global dof", std::to_string( Vh->nDof() ) } );
    report.add_row( { "lambda face global dof", std::to_string( Mh->nDof() ) } );
    report.add_row( { "hybrid global dof", std::to_string( Xh.nDof() ) } );
    Feel::cout << report << std::endl;

    Table diagnostics;
    diagnostics.add_row( { "||j-j_h||_L2", formatScientific( errorJ ) } );
    diagnostics.add_row( { "||div(j_h)-g||_L2", formatScientific( errorDivJ ) } );
    diagnostics.add_row( { "||V-V_h||_L2", formatScientific( errorV ) } );
    diagnostics.add_row( { "||rho*j_h+grad(V)-f||_L2", formatScientific( darcyError ) } );
    diagnostics.add_row( { "weak mass residual (Euclidean norm)", formatScientific( massResidual ) } );
    diagnostics.add_row( { "||j_h.n-j_N||_L2(Neumann)", formatScientific( sideFluxResidual ) } );
    Feel::cout << diagnostics << std::endl;

    if ( !boption( "no-export" ) )
    {
        auto e = exporter( _mesh = mesh, _name = "qs_darcy_hybrid" );
        e->addRegions();
        std::set<std::string> const representations{ "element", "nodal" };
        e->add( "j", jh, representations );
        e->add( "V", VhSol, representations );
        e->add( "j_exact", jexact );
        e->add( "V_exact", Vexact );
        e->save();
    }

    return 0;
}

int
darcyHybrid3D()
{
    int const fluxOrder = ioption( "darcy.order" );
    int const geometryOrder = ioption( "geometry.order" );
    std::string const fluxFamily = soption( "darcy.flux-family" );

    CHECK_GE( fluxOrder, 0 ) << "darcy.order must be non-negative";
    CHECK( fluxFamily == "rt" || fluxFamily == "bdm" )
        << "unsupported darcy.flux-family='" << fluxFamily << "'; expected 'rt' or 'bdm'";
    CHECK( fluxFamily != "bdm" || fluxOrder >= 1 )
        << "BDM uses mathematical orders k >= 1; darcy.order=0 is invalid for darcy.flux-family=bdm";
    CHECK( fluxFamily != "bdm" || fluxOrder <= 2 || boption( "sc.condense" ) )
        << "BDM orders >= 3 currently require --sc.condense=true: the BDM3 monolithic trace system has an unresolved null mode";
    CHECK_GE( geometryOrder, 1 ) << "geometry.order must be positive";
    auto const runtimeFluxOrder = RuntimeOrder::checked( fluxOrder );
    auto const runtimeGeometryOrder = RuntimeOrder::checked( geometryOrder );

    tic();
    auto mesh = loadMesh( _mesh = new Mesh<Simplex<3, Dynamic>>( runtimeGeometryOrder ) );
    CHECK_EQ( mesh->order(), geometryOrder ) << "loaded mesh order must match geometry.order";
    double const meshTime = toc( "mesh generation/loading" );
    auto faceMesh = createSubmesh( _mesh = mesh, _range = faces( mesh ), _update = 0 );
    CHECK_EQ( faceMesh->order(), geometryOrder ) << "face mesh must preserve geometry.order";
    auto Mh = Pdh<Dynamic>( faceMesh, runtimeFluxOrder );

    if ( fluxFamily == "rt" )
    {
        auto Jh = brokenh<RaviartThomas<Dynamic>>( mesh, runtimeFluxOrder );
        auto Vh = Pdh<Dynamic>( mesh, runtimeFluxOrder );
        DarcySpaceDescriptor const descriptor{ "RT", fluxOrder, fluxOrder, fluxOrder, fluxOrder };
        return darcyHybrid3D( mesh, Jh, Vh, Mh, descriptor, meshTime );
    }

    // BDM's implementation order and its divergence image are k-1;
    // its normal trace, constructed above, has mathematical order k.
    auto const pressureOrder = RuntimeOrder::checked( fluxOrder - 1 );
    auto Jh = brokenh<BrezziDouglasMarini<Dynamic>>( mesh, pressureOrder );
    auto Vh = Pdh<Dynamic>( mesh, pressureOrder );
    DarcySpaceDescriptor const descriptor{ "BDM", fluxOrder, fluxOrder - 1,
                                            fluxOrder, fluxOrder - 1 };
    return darcyHybrid3D( mesh, Jh, Vh, Mh, descriptor, meshTime );
}

} // namespace Feel

int
main( int argc, char** argv )
{
    using namespace Feel;

    try
    {
        po::options_description options = makeOptions();
        options.add( case_options( 3, "hybrid-compatible-Hdiv-Pdisc-trace" ) );

        Environment env( _argc = argc,
                         _argv = argv,
                         _desc = options,
                         _about = makeAbout() );
        return darcyHybrid3D();
    }
    catch ( ... )
    {
        handleExceptions();
    }
    return 1;
}
