// -*- coding: utf-8; mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
#include <feel/feel.hpp>

#include <cmath>
#include <functional>
#include <iomanip>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

struct RunConfig
{
    int order = 1;
    std::string mode = "runtime";
    bool convergence = false;
    int nRefine = 1;
    double refineRatio = 2.0;
};

struct RunResult
{
    double h = 0;
    std::size_t nDof = 0;
    int runtimeOrder = 0;
    double l2 = 0;
    double h1 = 0;
};

template<int Dim>
Feel::Symbols makeSymbols()
{
    if constexpr ( Dim == 2 )
        return Feel::Symbols{ "x", "y" };
    else
        return Feel::Symbols{ "x", "y", "z" };
}

template<int Dim>
std::string defaultExactFunction()
{
    if constexpr ( Dim == 2 )
        return "x+y";
    else
        return "x+y+z";
}

template<int Dim>
std::string exactFunctionExpression()
{
    auto g = Feel::soption( Feel::_name = "functions.g" );
    if ( g.empty() )
        g = defaultExactFunction<Dim>();
    return g;
}

double roundSmallToZero( double value )
{
    return ( std::abs( value ) < 1e-12 ) ? 0.0 : value;
}

std::optional<double> computeConvergenceRate( double errCur, double errPrev, double hCur, double hPrev )
{
    if ( errCur <= 0.0 || errPrev <= 0.0 || hCur <= 0.0 || hPrev <= 0.0 )
        return std::nullopt;

    auto denom = std::log( hCur / hPrev );
    if ( std::abs( denom ) < 1e-14 )
        return std::nullopt;

    return std::log( errCur / errPrev ) / denom;
}

std::string scientific( double value, int precision = 6 )
{
    std::ostringstream os;
    os << std::scientific << std::setprecision( precision ) << value;
    return os.str();
}

void printConvergenceTable( std::vector<RunResult> const& results, std::string const& g )
{
    if ( !Feel::Environment::isMasterRank() || results.empty() )
        return;

    std::cout << "Convergence table (g = " << g << ")\n";
    std::cout << std::left
              << std::setw( 7 ) << "level"
              << std::setw( 14 ) << "h"
              << std::setw( 10 ) << "runtimeP"
              << std::setw( 12 ) << "nDof"
              << std::setw( 16 ) << "L2"
              << std::setw( 14 ) << "rate(L2)"
              << std::setw( 16 ) << "H1"
              << std::setw( 14 ) << "rate(H1)" << "\n";

    for ( std::size_t i = 0; i < results.size(); ++i )
    {
        auto const& r = results[i];
        auto rateL2 = ( i == 0 ) ? std::optional<double>{}
                                 : computeConvergenceRate( r.l2, results[i - 1].l2, r.h, results[i - 1].h );
        auto rateH1 = ( i == 0 ) ? std::optional<double>{}
                                 : computeConvergenceRate( r.h1, results[i - 1].h1, r.h, results[i - 1].h );

        auto rateL2Str = rateL2 ? scientific( *rateL2, 4 ) : std::string( "-" );
        auto rateH1Str = rateH1 ? scientific( *rateH1, 4 ) : std::string( "-" );
        auto runtimeOrder = "P" + std::to_string( r.runtimeOrder );

        std::cout << std::left
                  << std::setw( 7 ) << i
                  << std::setw( 14 ) << scientific( r.h, 4 )
                  << std::setw( 10 ) << runtimeOrder
                  << std::setw( 12 ) << r.nDof
                  << std::setw( 16 ) << scientific( r.l2, 6 )
                  << std::setw( 14 ) << rateL2Str
                  << std::setw( 16 ) << scientific( r.h1, 6 )
                  << std::setw( 14 ) << rateH1Str << "\n";
    }
}

template<int Dim, typename SpacePtrType>
RunResult solveSimplexLaplacian( int order,
                                 std::string const& modeTag,
                                 SpacePtrType const& Vh,
                                 bool printDetails,
                                 bool exportResults )
{
    using namespace Feel;

    auto mesh = Vh->mesh();
    auto u = Vh->element( "u" );
    auto v = Vh->element( "v" );
    if ( printDetails && Environment::isMasterRank() )
    {
        auto const& fe = Vh->dof()->fe();
        std::cout << "dof table local dof/component : " << Vh->dof()->nLocalDof( true ) << "\n";
        std::cout << "dof table local dof/element   : " << Vh->dof()->nLocalDof() << "\n";
        if constexpr ( requires { fe.runtimeLocalDof(); } )
            std::cout << "fe runtimeLocalDof            : " << fe.runtimeLocalDof() << "\n";
        if constexpr ( requires { fe.nbDof(); } )
            std::cout << "fe nbDof                      : " << fe.nbDof() << "\n";
    }

    auto g = exactFunctionExpression<Dim>();
    auto vars = makeSymbols<Dim>();
    auto lap = GiNaC::laplacian( g, vars );
    auto gradg = GiNaC::grad( g, vars );

    auto cacheSeed = g + "_d" + std::to_string( Dim ) + "_p" + std::to_string( order );
    auto cacheId = std::to_string( std::hash<std::string>{}( cacheSeed ) );
    auto gexpr = expr<8>( g, vars, "gexpr_" + cacheId );
    auto lapexpr = expr<1, 1, 8>( lap, vars, "lapg", "lapexpr_" + cacheId );
    auto gradexpr = expr<1, Dim, 8>( gradg, vars, "gradg", "gradexpr_" + cacheId );

    auto elementRange = elements( mesh );
    auto rhsExpr = -lapexpr * id( v );
    auto lhsExpr = gradt( u ) * trans( grad( v ) );
    auto l2ErrExpr = idv( u ) - gexpr;
    auto h1ErrExpr = gradv( u ) - gradexpr;

    auto autoQuadOrder = [&]<typename ExprType>( ExprType const& ex ) {
        using range_type = std::decay_t<decltype( elementRange )>;
        using expr_type = std::decay_t<ExprType>;
        return static_cast<int>( vf::ExpressionOrder<range_type, expr_type>::value( ex ) );
    };

    int const quadOverride = ioption( "case.quad" );
    int const quadRhs = ( quadOverride > 0 ) ? quadOverride : autoQuadOrder( rhsExpr );
    int const quadLhs = ( quadOverride > 0 ) ? quadOverride : autoQuadOrder( lhsExpr );
    int const quadL2 = ( quadOverride > 0 ) ? quadOverride : autoQuadOrder( inner( l2ErrExpr ) );
    int const quadH1 = ( quadOverride > 0 ) ? quadOverride : autoQuadOrder( inner( h1ErrExpr ) );

    auto l = form1( _test = Vh );
    l = integrate( _range = elementRange, _expr = rhsExpr, _quad = quadRhs );

    auto a = form2( _trial = Vh, _test = Vh );
    a = integrate( _range = elementRange, _expr = lhsExpr, _quad = quadLhs );

    // Export the unconstrained operator before Dirichlet elimination when requested.
    if ( !soption( _name = "export-matlab" ).empty() )
    {
        a.matrixPtr()->printMatlab( "A_before_dirichlet.m" );
        l.vectorPtr()->printMatlab( "b_before_dirichlet.m" );
    }

    a += on( _range = boundaryfaces( mesh ), _rhs = l, _element = u, _expr = gexpr );
    a.solve( _rhs = l, _solution = u, _rebuild = true );

    auto l2 = normL2( _range = elementRange, _expr = l2ErrExpr, _quad = quadL2 );
    auto semih1 = normL2( _range = elementRange, _expr = h1ErrExpr, _quad = quadH1 );
    auto h1 = math::sqrt( l2 * l2 + semih1 * semih1 );
    auto l2Display = roundSmallToZero( l2 );
    auto h1Display = roundSmallToZero( h1 );
    if ( printDetails && Environment::isMasterRank() )
    {
        std::cout << "simplex dim/order      : " << Dim << "D/P" << order << " (" << modeTag << ")\n";
        if ( quadOverride > 0 )
            std::cout << "quadrature order       : " << quadOverride << " (override)\n";
        else
            std::cout << "quadrature order       : rhs=" << quadRhs
                      << " lhs=" << quadLhs
                      << " l2=" << quadL2
                      << " h1=" << quadH1 << " (auto)\n";
        std::cout << "Exact solution g       : " << g << "\n";
        std::cout << "lap(g)                 : " << lap << "\n";
        std::cout << "grad(g)                : " << gradg << "\n";
        std::cout << "L2 error  ||u-g||_L2   = " << l2Display << "\n";
        std::cout << "H1 error  ||u-g||_H1   = " << h1Display << "\n";
    }

    if ( exportResults )
    {
        auto e = exporter( _mesh = mesh );
        e->addRegions();
        auto lapg = Vh->element( "lapg" );
        lapg.on( _range = elements( mesh ), _expr = lapexpr );
        e->add( "u", u );
        v = vf::project( _space = Vh, _range = elements( mesh ), _expr = gexpr );
        e->add( "uexact", v );
        e->add( "lapg", lapg );
        e->add( "gradg", gradexpr, "element" );
        e->save();
    }

    int runtimeOrder = order;
    if constexpr ( requires { Vh->runtimeOrder(); } )
        runtimeOrder = static_cast<int>( Vh->runtimeOrder() );

    RunResult result;
    result.h = mesh->hAverage();
    result.nDof = Vh->nDof();
    result.runtimeOrder = runtimeOrder;
    result.l2 = l2Display;
    result.h1 = h1Display;
    return result;
}

template<int Dim>
    requires ( Dim >= 1 && Dim <= 3 )
RunResult runSimplexLaplacianOnMesh( std::shared_ptr<Feel::Mesh<Feel::Simplex<Dim>>> const& mesh,
                                     RunConfig const& config,
                                     bool printDetails,
                                     bool exportResults )
{
    using namespace Feel;

    auto runRuntime = [&]() {
        auto Vh = Pch<Dynamic>( mesh, RuntimeOrder( config.order ) );
        return solveSimplexLaplacian<Dim>( config.order, "runtime", Vh, printDetails, exportResults );
    };
    auto runCompile = [&]<int P>() {
        auto Vh = Pch<P>( mesh );
        return solveSimplexLaplacian<Dim>( config.order, "compile", Vh, printDetails, exportResults );
    };

    //if ( config.mode == "runtime" )
    return runRuntime();
    //if ( config.mode == "compile" )
    //{
    //    switch ( config.order )
    //    {
    //    case 1: return runCompile.template operator()<1>();
    //    case 2: return runCompile.template operator()<2>();
    //    case 3: return runCompile.template operator()<3>();
    //    //case 4: return runCompile.template operator()<4>();
    //    default:
    //        CHECK( false ) << "compile mode supports case.order in {1,2,3}, got " << config.order;
    //    }
    //}
    CHECK( false ) << "invalid case.mode=" << config.mode << ", expected runtime or compile";

    return {};
}

template<int Dim>
    requires ( Dim >= 1 && Dim <= 3 )
int runSimplexLaplacian( RunConfig const& config )
{
    using namespace Feel;

    if ( !config.convergence )
    {
        auto mesh = loadMesh( _mesh = new Mesh<Simplex<Dim>> );
        runSimplexLaplacianOnMesh<Dim>( mesh, config, true, true );
        return 0;
    }

    CHECK( config.nRefine >= 1 ) << "invalid case.nrefine=" << config.nRefine << ", expected >= 1";
    CHECK( config.refineRatio > 1.0 ) << "invalid case.refine-ratio=" << config.refineRatio << ", expected > 1";

    std::vector<RunResult> results;
    results.reserve( config.nRefine );

    auto h0 = doption( "gmsh.hsize" );
    for ( int level = 0; level < config.nRefine; ++level )
    {
        auto hLevel = h0 / std::pow( config.refineRatio, static_cast<double>( level ) );
        auto mesh = loadMesh( _mesh = new Mesh<Simplex<Dim>>, _h = hLevel );
        auto result = runSimplexLaplacianOnMesh<Dim>( mesh, config, false, false );
        results.push_back( result );
    }

    printConvergenceTable( results, exactFunctionExpression<Dim>() );
    return 0;
}

int main( int argc, char** argv )
{
    using namespace Feel;

    try
    {
        po::options_description desc( "simple_lap (simplex, dynamic order)" );
        desc.add( case_options( FEELPP_DIM, "P1" ) ).add_options()
            ( "case.order", po::value<int>()->default_value( 1 ), "polynomial order (P_k), k>=1" )
            ( "case.quad", po::value<int>()->default_value( -1 ), "quadrature order (>0: override, <=0: auto from expression order)" )
            ( "case.mode", po::value<std::string>()->default_value( "runtime" ), "space mode: runtime or compile" )
            ( "case.convergence", po::value<bool>()->default_value( false ), "run convergence mode over mesh refinements" )
            ( "case.nrefine", po::value<int>()->default_value( 1 ), "number of refinement levels in convergence mode" )
            ( "case.refine-ratio", po::value<double>()->default_value( 2.0 ), "mesh refinement ratio (>1) between levels" );

        Environment env( _argc = argc, _argv = argv,
                         _desc = desc,
                         _about = about( _name = "simple_lap",
                                         _author = "Feel++ Consortium",
                                         _email = "feelpp-devel@feelpp.org" ) );

        RunConfig config;
        config.order = ioption( "case.order" );
        config.mode = soption( "case.mode" );
        config.convergence = boption( "case.convergence" );
        config.nRefine = ioption( "case.nrefine" );
        config.refineRatio = doption( "case.refine-ratio" );

        CHECK( config.order >= 1 ) << "invalid case.order=" << config.order << ", expected >= 1";
        return runSimplexLaplacian<FEELPP_DIM>( config );
    }
    catch ( ... )
    {
        handleExceptions();
    }
    return 1;
}
