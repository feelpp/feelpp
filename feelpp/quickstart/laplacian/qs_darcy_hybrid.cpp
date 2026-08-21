/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/table.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelpoly/discontinuous.hpp>
#include <feel/feelpoly/raviartthomas.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeltiming/tic.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/vf.hpp>

#include <iomanip>
#include <sstream>
#include <type_traits>

namespace Feel
{
namespace fem
{

// Static condensation requires element-local flux degrees of freedom.
template<uint16_type N, uint16_type O, typename T = double,
         template<uint16_type, uint16_type, uint16_type> class Convex = Simplex,
         uint16_type TheTAG = 0>
class BrokenRaviartThomas : public RaviartThomas<N, O, T, Convex, TheTAG>
{
public:
    static inline const bool isContinuous = false;
    using continuity_type = Discontinuous;
};

} // namespace fem

template<uint16_type Order, uint16_type TheTAG = 0>
class BrokenRaviartThomas
{
public:
    template<uint16_type N, uint16_type R = N, typename T = double, typename Convex = Simplex<N>>
    struct apply
    {
        using result_type = std::conditional_t<
            Convex::is_simplex,
            fem::BrokenRaviartThomas<N, Order, T, Simplex, TheTAG>,
            fem::BrokenRaviartThomas<N, Order, T, Hypercube, TheTAG>>;
        using type = result_type;
    };

    template<uint16_type NewTAG>
    struct ChangeTag
    {
        using type = BrokenRaviartThomas<Order, NewTAG>;
    };

    using component_basis_type = Lagrange<Order, Scalar>;

    static const uint16_type TAG = TheTAG;
};

template<int Order, typename MeshType, typename T = double>
using broken_dh_type = FunctionSpace<MeshType, bases<BrokenRaviartThomas<Order>>, T,
                                     Periodicity<NoPeriodicity>>;

template<int Order, typename MeshType>
inline std::shared_ptr<broken_dh_type<Order, MeshType>>
BrokenDh( std::shared_ptr<MeshType> const& mesh )
{
    return broken_dh_type<Order, MeshType>::New(
        _mesh = mesh, _worldscomm = makeWorldsComm( 1, mesh->worldComm() ) );
}

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
    AboutData about( "qs_darcy_hybrid",
                     "qs_darcy_hybrid",
                     "0.1",
                     "Quickstart 3D hybrid-mixed Darcy/electrostatic problem",
                     AboutData::License_GPL,
                     "Copyright (C) 2026 University of Strasbourg" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}

int
darcyHybrid3D()
{
    using namespace boost::hana::literals;

    tic();
    auto mesh = loadMesh( _mesh = new Mesh<Simplex<3>> );
    double const meshTime = toc( "mesh generation/loading" );

    tic();
    auto Jh = BrokenDh<0>( mesh );
    auto Vh = Pdh<0>( mesh );
    auto faceMesh = createSubmesh( _mesh = mesh, _range = faces( mesh ), _update = 0 );
    auto Mh = Pdh<0>( faceMesh );
    auto Xh = product( Jh, Vh, Mh );
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

    auto U = Xh.element();
    auto T = Xh.element();
    auto j = U( 0_c );
    auto V = U( 1_c );
    auto lambda = U( 2_c );
    auto w = T( 0_c );
    auto q = T( 1_c );
    auto mu = T( 2_c );

    tic();
    backend( _rebuild = true );
    solve::strategy const strategy = boption( "sc.condense" )
                                        ? solve::strategy::static_condensation
                                        : solve::strategy::monolithic;
    auto a = blockform2( Xh, strategy, backend() );
    auto l = blockform1( Xh, strategy, backend() );

    // Local RT0-P0 mixed equations. The face potential lambda couples cells.
    l( 0_c ) += integrate( _range = elements( mesh ), _expr = inner( f, id( w ) ) );
    l( 1_c ) += integrate( _range = elements( mesh ), _expr = g * id( q ) );
    l( 2_c ) += integrate( _range = neumannFaces, _expr = id( mu ) * jn );
    l( 2_c ) += integrate( _range = dirichletFaces, _expr = id( mu ) * Vexact );

    a( 0_c, 0_c ) += integrate( _range = elements( mesh ), _expr = rho * inner( idt( j ), id( w ) ) );
    a( 0_c, 1_c ) += integrate( _range = elements( mesh ), _expr = -idt( V ) * div( w ) );
    a( 0_c, 2_c ) += integrate( _range = internalfaces( mesh ),
                                 _expr = idt( lambda ) * ( leftface( normal( w ) ) + rightface( normal( w ) ) ) );
    a( 0_c, 2_c ) += integrate( _range = boundaryfaces( mesh ), _expr = idt( lambda ) * normal( w ) );

    a( 1_c, 0_c ) += integrate( _range = elements( mesh ), _expr = divt( j ) * id( q ) );
    // The condenser expects every local block to be registered, including this exact zero.
    a( 1_c, 1_c ) += integrate( _range = elements( mesh ), _expr = cst( 0. ) * idt( V ) * id( q ) );

    a( 2_c, 0_c ) += integrate( _range = internalfaces( mesh ),
                                 _expr = id( mu ) * ( leftfacet( normalt( j ) ) + rightfacet( normalt( j ) ) ) );
    a( 2_c, 0_c ) += integrate( _range = neumannFaces, _expr = id( mu ) * normalt( j ) );
    a( 2_c, 2_c ) += integrate( _range = dirichletFaces, _expr = idt( lambda ) * id( mu ) );
    double const assemblyTime = toc( "space assembly" );

    tic();
    auto solveStatus = a.solve( _rhs = l, _solution = U, _condense = boption( "sc.condense" ) );
    double const solveTime = toc( "linear solve" );
    CHECK( solveStatus.isConverged() )
        << "hybrid Darcy solve did not converge; iterations=" << solveStatus.nIterations()
        << ", residual=" << solveStatus.residual();

    auto jh = U( 0_c );
    auto VhSol = U( 1_c );
    double const errorJ = normL2( _range = elements( mesh ), _expr = idv( jh ) - jexact );
    double const errorDivJ = normL2( _range = elements( mesh ), _expr = divv( jh ) - div( jexact ) );
    double const errorV = normL2( _range = elements( mesh ), _expr = idv( VhSol ) - Vexact );
    double const darcyError = normL2( _range = elements( mesh ),
                                      _expr = rho * idv( jh ) + trans( grad<3>( Vexact ) ) - f );
    double const sideFluxResidual = normL2( _range = neumannFaces,
                                            _expr = inner( idv( jh ), N() ) - jn );
    auto massResidualField = project( _space = Vh, _range = elements( mesh ),
                                      _expr = divv( jh ) - g );
    double const massResidual = normL2( _range = elements( mesh ),
                                        _expr = idv( massResidualField ) );

    double const tolerance = doption( "checks.tolerance" );
    CHECK_LT( massResidual, tolerance ) << "div(j_h) must equal the P0 projection of g cellwise";
    CHECK_LT( sideFluxResidual, tolerance ) << "j_h.n must satisfy the prescribed Neumann data";

    Table report;
    report.add_row( { "mesh generation/loading (s)", formatScientific( meshTime ) } );
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
    diagnostics.add_row( { "||div(j-j_h)||_L2", formatScientific( errorDivJ ) } );
    diagnostics.add_row( { "||V-V_h||_L2", formatScientific( errorV ) } );
    diagnostics.add_row( { "||rho*j_h+grad(V)-f||_L2", formatScientific( darcyError ) } );
    diagnostics.add_row( { "||P0(div(j_h)-g)||_L2", formatScientific( massResidual ) } );
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
        return darcyHybrid3D();
    }
    catch ( ... )
    {
        handleExceptions();
    }
    return 1;
}
