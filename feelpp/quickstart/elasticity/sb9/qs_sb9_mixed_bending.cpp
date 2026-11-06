/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/tensorformat.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhm.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeltiming/tic.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/sb9_bending.hpp>
#include <feel/feelvf/sb9_pinching.hpp>
#include <feel/feelvf/sb9_quadrature.hpp>
#include <feel/feelvf/sb9_shear.hpp>
#include <feel/feelvf/sb9_stabilization.hpp>
#include <feel/feelvf/sb9_strain.hpp>
#include <feel/feelvf/vf.hpp>

#include "qs_shell_benchmark_framework.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <numbers>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <set>
#include <vector>

using namespace Feel;

namespace
{
namespace qsb = Feel::Quickstart::ShellBenchmark;
using BenchmarkConfig = qsb::BenchmarkConfig;
using json = qsb::json;
using mesh_type = qsb::mesh_type;
using qsb::benchmarkPreset;
using qsb::buildMesh;
using qsb::loadBenchmarkSpecs;
using qsb::shellVectorExpression;

SymmetricTensorFormat
matlabEpsilonFormat()
{
    // validation.pdf / JSON epsilon vectors are xx,yy,zz,xy,xz,yz with
    // engineering shear entries.
    return { SymmetricTensorOrder::DiagonalFirst,
             SymmetricTensorScaling::EngineeringShear };
}

SymmetricTensorFormat
matlabSigmaFormat()
{
    // validation.pdf / JSON sigma vectors use the same order but keep tensor
    // shear entries.
    return { SymmetricTensorOrder::DiagonalFirst,
             SymmetricTensorScaling::Tensor };
}

node_type
toNode( qsb::shell_vec const& p )
{
    node_type n( 3 );
    n( 0 ) = p( 0 );
    n( 1 ) = p( 1 );
    n( 2 ) = p( 2 );
    return n;
}

template<typename PointRangeType>
double
pointLoadScale( PointRangeType const& pointRange, std::string const& quantity, std::string const& label )
{
    if ( quantity == "per-point" || quantity == "per_point" || quantity == "point" )
        return 1.0;

    CHECK( quantity == "total" ) << label << " quantity must be 'per-point' or 'total'";
    size_type const nMarkedPoints = nelements( pointRange, true );
    CHECK( nMarkedPoints > 0 ) << label << " uses quantity=total but selects no marked point";
    return 1.0 / nMarkedPoints;
}

template <auto Kind, typename CacheType>
Eigen::Matrix<double, 6, 24>
sb9CoefficientMatrix( CacheType const& cache )
{
    Eigen::Matrix<double, 6, 24> matrix;
    matrix.setZero();

    for ( uint16_type component = 0; component < 3; ++component )
    {
        for ( uint16_type node = 0; node < 8; ++node )
        {
            Eigen::Matrix<double, 6, 1> coeff;
            coeff.setZero();
            cache.template fillVectorCoefficients<Kind>( coeff, node, component );
            matrix.col( component*8 + node ) = coeff;
        }
    }

    return matrix;
}

template <typename SpaceType, typename FieldElementType>
qsb::shell_voigt
evaluateSymmetricFieldReference( SpaceType const& Sh,
                                 FieldElementType const& field,
                                 qsb::shell_vec const& point,
                                 SymmetricTensorFormat const& format )
{
    auto ctx = Sh->context();
    ctx.add( toNode( point ) );
    auto values = field.evaluateSymmetric( ctx, format, true );
    CHECK( values.rows() >= 1 ) << "expected at least one tensor-field evaluation point";
    CHECK( values.cols() == 6 ) << "SB9 validation expects six tensor components";

    qsb::shell_voigt result;
    for ( int c = 0; c < 6; ++c )
        result( c ) = values( 0, c );
    return result;
}

struct SymmetricTensorComponents
{
    double xx = 0.0;
    double yy = 0.0;
    double zz = 0.0;
    double xy = 0.0;
    double xz = 0.0;
    double yz = 0.0;
};

SymmetricTensorComponents
toSymmetricTensorComponents( qsb::shell_voigt const& values, bool engineeringShear )
{
    double const shearScale = engineeringShear ? 0.5 : 1.0;
    return { values( 0 ),
             values( 1 ),
             values( 2 ),
             shearScale * values( 3 ),
             shearScale * values( 4 ),
             shearScale * values( 5 ) };
}

template <typename FieldElementType, typename MeshType>
void
assignSymmetricTensorOnElement( FieldElementType& field,
                                std::shared_ptr<MeshType> const& mesh,
                                index_type elementId,
                                qsb::shell_voigt const& values,
                                bool engineeringShear )
{
    auto const tensor = toSymmetricTensorComponents( values, engineeringShear );
    field.tensorComponent( Component::X, Component::X ).on( _range=idedelements( mesh, elementId ), _expr=cst( tensor.xx ) );
    field.tensorComponent( Component::Y, Component::Y ).on( _range=idedelements( mesh, elementId ), _expr=cst( tensor.yy ) );
    field.tensorComponent( Component::Z, Component::Z ).on( _range=idedelements( mesh, elementId ), _expr=cst( tensor.zz ) );
    field.tensorComponent( Component::X, Component::Y ).on( _range=idedelements( mesh, elementId ), _expr=cst( tensor.xy ) );
    field.tensorComponent( Component::X, Component::Z ).on( _range=idedelements( mesh, elementId ), _expr=cst( tensor.xz ) );
    field.tensorComponent( Component::Y, Component::Z ).on( _range=idedelements( mesh, elementId ), _expr=cst( tensor.yz ) );
}

template <typename GmcPtrType, typename DisplacementElementType, typename AlphaElementType>
std::pair<qsb::shell_voigt, qsb::shell_voigt>
evaluateSb9AtGmc( GmcPtrType const& gmc,
                  DisplacementElementType const& uh,
                  AlphaElementType const& alphah,
                  double lambda,
                  double mu,
                  double alphaScale,
                  double pinchingBpzScale,
                  double shellShearFactor )
{
    auto const data = Feel::vf::detail::computeShellCellGeometry( gmc.get() );
    using geometry_data_type = std::decay_t<decltype( data )>;
    Feel::vf::detail::SB9BendingKernelCache<geometry_data_type> const bendingCache( data );
    Feel::vf::detail::SB9PinchingKernelCache<geometry_data_type> const pinchingCache( data );
    Feel::vf::detail::SB9ShearKernelCache<geometry_data_type> const shearCache( data );

    Eigen::Matrix<double, 24, 1> uLocal;
    uLocal.setZero();
    index_type const elementId = gmc->element().id();
    for ( uint16_type component = 0; component < 3; ++component )
        for ( uint16_type node = 0; node < 8; ++node )
            uLocal( component*8 + node ) = uh.localToGlobal( elementId, node, component );

    auto const bm0 = sb9CoefficientMatrix<Feel::vf::detail::SB9BendingKind::Bm0>( bendingCache );
    auto const bb0 = sb9CoefficientMatrix<Feel::vf::detail::SB9BendingKind::Bb0>( bendingCache );
    auto const bpc = sb9CoefficientMatrix<Feel::vf::detail::SB9PinchingKind::Bpc>( pinchingCache );
    auto const bpz = sb9CoefficientMatrix<Feel::vf::detail::SB9PinchingKind::Bpz>( pinchingCache );
    auto const bc0 = sb9CoefficientMatrix<Feel::vf::detail::SB9ShearKind::Bc0>( shearCache );

    double const zetaValue = gmc->xRefs()( 2, 0 );
    double const shearWeight = shellShearFactor * ( 1.0 - zetaValue*zetaValue );
    auto const membrane = ( bm0 + zetaValue*bb0 ) * uLocal;
    auto const pinching = ( bpc + pinchingBpzScale*zetaValue*bpz ) * uLocal;
    auto const shear = bc0 * uLocal;

    Eigen::Matrix<double, 6, 1> epsilonMandel;
    epsilonMandel.setZero();
    epsilonMandel( 0 ) = membrane( 0 );
    epsilonMandel( 1 ) = std::numbers::sqrt2_v<double> * membrane( 3 );
    epsilonMandel( 2 ) = std::numbers::sqrt2_v<double> * shearWeight * shear( 4 );
    epsilonMandel( 3 ) = membrane( 1 );
    epsilonMandel( 4 ) = std::numbers::sqrt2_v<double> * shearWeight * shear( 5 );
    epsilonMandel( 5 ) = pinching( 2 ) +
                         alphaScale * ( -4.0*zetaValue/data.thickness ) * alphah.localToGlobal( elementId, 0, 0 );

    qsb::shell_voigt epsilon;
    epsilon << epsilonMandel( 0 ),
               epsilonMandel( 3 ),
               epsilonMandel( 5 ),
               std::numbers::sqrt2_v<double> * epsilonMandel( 1 ),
               std::numbers::sqrt2_v<double> * epsilonMandel( 2 ),
               std::numbers::sqrt2_v<double> * epsilonMandel( 4 );

    qsb::shell_voigt sigma;
    double const traceEpsilon = epsilon( 0 ) + epsilon( 1 ) + epsilon( 2 );
    sigma << lambda*traceEpsilon + 2.0*mu*epsilon( 0 ),
             lambda*traceEpsilon + 2.0*mu*epsilon( 1 ),
             lambda*traceEpsilon + 2.0*mu*epsilon( 2 ),
             mu*epsilon( 3 ),
             mu*epsilon( 4 ),
             mu*epsilon( 5 );

    return { epsilon, sigma };
}

template <typename SpaceType, typename DisplacementElementType, typename AlphaElementType>
std::pair<qsb::shell_voigt, qsb::shell_voigt>
evaluateSb9FieldReference( SpaceType const& Uh,
                           DisplacementElementType const& uh,
                           AlphaElementType const& alphah,
                           qsb::shell_vec const& point,
                           double lambda,
                           double mu,
                           double alphaScale,
                           double pinchingBpzScale,
                           double shellShearFactor )
{
    Eigen::Matrix<double, 6, 1> localEpsilon;
    Eigen::Matrix<double, 6, 1> localSigma;
    localEpsilon.setZero();
    localSigma.setZero();

    auto ctx = Uh->context();
    ctx.add( toNode( point ) );

    for ( auto const& [ctxId, ctxEntry] : ctx )
    {
        (void)ctxId;
        auto const& basisContext = std::get<0>( ctxEntry );
        auto const& gmc = basisContext->gmContext();
        auto const [epsilonValue, sigmaValue] =
            evaluateSb9AtGmc( gmc, uh, alphah, lambda, mu, alphaScale, pinchingBpzScale, shellShearFactor );
        localEpsilon += epsilonValue;
        localSigma += sigmaValue;
    }

    qsb::shell_voigt epsilon;
    qsb::shell_voigt sigma;
    if ( Uh->worldComm().globalSize() > 1 )
    {
        mpi::all_reduce( Uh->worldComm(), localEpsilon, epsilon,
                         []( qsb::shell_voigt const& x, qsb::shell_voigt const& y )
                         {
                             return x + y;
                         } );
        mpi::all_reduce( Uh->worldComm(), localSigma, sigma,
                         []( qsb::shell_voigt const& x, qsb::shell_voigt const& y )
                         {
                             return x + y;
                         } );
    }
    else
    {
        epsilon = localEpsilon;
        sigma = localSigma;
    }

    return { epsilon, sigma };
}

template <typename MeshPtrType, typename DisplacementElementType, typename AlphaElementType,
          typename EpsilonElementType, typename SigmaElementType>
void
fillSb9SymmetricFields( MeshPtrType const& mesh,
                        DisplacementElementType const& uh,
                        AlphaElementType const& alphah,
                        EpsilonElementType& epsilonh,
                        SigmaElementType& sigmah,
                        double lambda,
                        double mu,
                        double alphaScale,
                        double pinchingBpzScale,
                        double shellShearFactor )
{
    typename mesh_type::gm_type::matrix_node_t_type referenceCenter( mesh_type::nDim, 1 );
    referenceCenter.clear();
    auto gmpc = mesh->gm()->preCompute( mesh->gm(), referenceCenter );

    for ( auto const& elementRef : elements( mesh ) )
    {
        auto const& element = boost::unwrap_ref( elementRef );
        auto gmc = mesh->gm()->template context<vm::POINT|vm::JACOBIAN|vm::KB>( element, gmpc );
        auto const [epsilonValue, sigmaValue] =
            evaluateSb9AtGmc( gmc, uh, alphah, lambda, mu, alphaScale, pinchingBpzScale, shellShearFactor );

        assignSymmetricTensorOnElement( epsilonh, mesh, element.id(), epsilonValue, true );
        assignSymmetricTensorOnElement( sigmah, mesh, element.id(), sigmaValue, false );
    }
}

template <typename SpaceType, typename EpsilonElementType, typename SigmaElementType>
bool
checkFieldReferences( qsb::BenchmarkConfig const& cfg,
                      SpaceType const& Sh,
                      EpsilonElementType const& epsilonh,
                      SigmaElementType const& sigmah )
{
    bool ok = true;

    for ( auto const& fieldReference : cfg.fieldReferences )
    {
        std::string target;
        if ( fieldReference.targets.find( "sb9g25" ) != fieldReference.targets.end() )
            target = "sb9g25";
        else if ( fieldReference.targets.find( "sb9" ) != fieldReference.targets.end() )
            target = "sb9";
        else
            continue;

        auto const [epsilonValue, sigmaValue] =
            std::make_pair( evaluateSymmetricFieldReference( Sh, epsilonh, fieldReference.point, matlabEpsilonFormat() ),
                            evaluateSymmetricFieldReference( Sh, sigmah, fieldReference.point, matlabSigmaFormat() ) );
        ok = qsb::reportFieldReferenceComparison( fieldReference, target, epsilonValue, sigmaValue,
                                                  boption( "check-reference" ), std::cout ) && ok;
    }

    return ok;
}

template <typename SpaceType>
auto
matlabFieldContext( SpaceType const& Sh,
                    qsb::BenchmarkConfig const& cfg,
                    qsb::shell_vec const& fallbackPoint )
{
    auto ctx = Sh->context();
    if ( cfg.fieldReferences.empty() )
    {
        ctx.add( toNode( fallbackPoint ) );
        return ctx;
    }

    for ( auto const& fieldReference : cfg.fieldReferences )
        ctx.add( toNode( fieldReference.point ) );
    return ctx;
}

po::options_description
makeOptions()
{
    po::options_description options( "qs_sb9_mixed_bending options" );
    options.add_options()
        ( "specs", po::value<std::string>(),
          "json benchmark specification file (typically provided by qs_sb9_mixed_bending.cfg)" )
        ( "benchmark", po::value<std::string>()->default_value( "square-plate" ),
          "benchmark name in the JSON specification file" )
        ( "E", po::value<double>()->default_value( -1.0 ), "override Young modulus from the benchmark preset" )
        ( "nu", po::value<double>()->default_value( -1.0 ), "override Poisson ratio from the benchmark preset" )
        ( "pressure", po::value<double>()->default_value( -1.0 ),
          "override every pressure load magnitude; positive pressure acts inward as -p*N()" )
        ( "alpha-scale", po::value<double>()->default_value( 1.0 ),
          "scale applied to the shell-normal SB9 scalar correction" )
        ( "pinching-bpz-scale", po::value<double>()->default_value( 1.0 ),
          "scale applied to the zeta*Bpz pinching term; 1 enables the full SB9 pinching formulation" )
        ( "shell-shear-factor", po::value<double>()->default_value( 1.25 ),
          "Hallquist-like transverse shear shape factor coefficient" )
        ( "shell-shear-stab", po::value<double>()->default_value( 0.25 ),
          "scale applied to the Hallquist Bc1/Bc2 stabilization block" )
        ( "shell-bs-stab", po::value<double>()->default_value( 1.0 ),
          "scale applied to the additional SB9 Bs1..Bs4 stabilization block; 1 matches the MATLAB Ds stabilization coefficients" )
        ( "functions.f", po::value<std::string>()->default_value( "" ), "override body force expression" )
        ( "no-solve", po::value<bool>()->default_value( false ), "assemble only" )
        ( "show-timings", po::value<bool>()->default_value( true ),
          "print per-term assembly and solve timings" )
        ( "check-reference", po::value<bool>()->default_value( true ),
          "fail when a benchmark reference value exceeds its configured tolerance" )
        ( "print-matlab-fields", po::value<bool>()->default_value( false ),
          "write u, alpha, epsilon, and sigma fields; tensor fields use the MATLAB validation order" )
        ( "export-thickness", po::value<bool>()->default_value( true ), "export shell thickness projected to P0" )
        ( "export-normal-displacement", po::value<bool>()->default_value( true ),
          "export inner(u,shellNormal()) projected to a scalar field" );
    return options.add( feel_options() );
}

AboutData
makeAbout()
{
    AboutData about( "qs_sb9_mixed_bending",
                     "qs_sb9_mixed_bending",
                     "0.1",
                     "SB9-inspired mixed Q1 x P0 shell bending prototype",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) Feel++ Consortium" );
    about.addAuthor( "Feel++ Consortium", "developer", "feelpp-devel@feelpp.org", "" );
    return about;
}
} // namespace

int
main( int argc, char** argv )
{
    using namespace Feel;
    using namespace vf;

    try
    {
        Environment env( _argc=argc, _argv=argv, _desc=makeOptions(), _about=makeAbout() );

        if ( !Environment::vm().count( "specs" ) || soption( "specs" ).empty() )
            throw std::invalid_argument( "missing required option 'specs'; set it in qs_sb9_mixed_bending.cfg or pass --specs" );

        std::string const specsPath = Environment::expand( soption( "specs" ) );
        auto const specs = loadBenchmarkSpecs( specsPath );
        auto cfg = benchmarkPreset( specs, soption( "benchmark" ) );
        auto meshData = buildMesh( cfg );
        auto mesh = meshData.mesh;

        double const E = ( doption( "E" ) > 0.0 ) ? doption( "E" ) : cfg.E;
        double const nu = ( doption( "nu" ) >= 0.0 ) ? doption( "nu" ) : cfg.nu;
        double const pressureOverride = doption( "pressure" );
        double const alphaScale = doption( "alpha-scale" );
        double const pinchingBpzScale = doption( "pinching-bpz-scale" );
        double const shellShearFactor = doption( "shell-shear-factor" );
        double const shellShearStab = doption( "shell-shear-stab" );
        double const shellBsStab = doption( "shell-bs-stab" );
        double const lambda = E * nu / ( ( 1.0 + nu ) * ( 1.0 - 2.0 * nu ) );
        double const mu = E / ( 2.0 * ( 1.0 + nu ) );
        std::string const bodyForceExpression =
            soption( "functions.f" ).empty() ? cfg.bodyForceExpression : soption( "functions.f" );

        std::cout << "benchmark: " << cfg.name << "\n"
                  << "description: " << cfg.description << "\n"
                  << "specs: " << specsPath << "\n"
                  << "E=" << E << ", nu=" << nu << "\n"
                  << "pinching-bpz-scale=" << pinchingBpzScale << "\n"
                  << "solver strategy: " << ( boption( "sc.condense" ) ? "static_condensation" : "monolithic" ) << "\n"
                  << "loads: pressures=" << cfg.pressureLoads.size()
                  << ", tractions=" << cfg.tractionLoads.size()
                  << ", total-forces=" << cfg.totalForceLoads.size()
                  << ", point-loads=" << cfg.pointLoads.size()
                  << ", point-moments=" << cfg.pointMoments.size() << "\n"
                  << "constraints: face-clamps=" << cfg.clampMarkers.size()
                  << ", point-constraints=" << cfg.pointConstraints.size() << "\n";

        bool const showTimings = boption( "show-timings" );
        auto ticIf = [showTimings]()
        {
            if ( showTimings )
                tic();
        };
        auto tocIf = [showTimings]( char const* label )
        {
            if ( showTimings )
                toc( label, true );
        };

        ticIf();
        auto Uh = Pchv<1>( mesh );
        auto Ah = Pdh<0>( mesh );
        auto Xh = productPtr( Uh, Ah );
        auto algebraBackend = backend( _rebuild=true, _worldcomm=Uh->worldCommPtr() );
        tocIf( "sb9.spaces" );
        solve::strategy const strategy =
            boption( "sc.condense" ) ? solve::strategy::static_condensation : solve::strategy::monolithic;

        auto u = trial( Uh, "u" );
        auto v = test( Uh, "v" );
        auto alpha = trial( Ah, "alpha" );
        auto beta = test( Ah, "beta" );
        auto X = Xh->element();
        auto uBc = X( 0_c );

        auto a = blockform2( *Xh, strategy, algebraBackend );
        auto zt = zeta();
        auto C = isotropic_stiffness<3>( lambda, mu );
        auto invJ0 = inv( shellJacobian0() );
        auto shellShearWeight = cst( shellShearFactor ) * ( cst( 1.0 ) - zt * zt );
        auto membraneTrial = sb9MembraneBending( u, zt );
        auto membraneTest = sb9MembraneBending( v, zt );
        auto pinchingTrial = sb9Pinching( u, zt, cst( pinchingBpzScale ) );
        auto pinchingTest = sb9Pinching( v, zt, cst( pinchingBpzScale ) );
        auto shearTrial = sb9Shear( u, shellShearWeight );
        auto shearTest = sb9Shear( v, shellShearWeight );
        auto epsShellTrialMandel = vec( component<0, 0>( membraneTrial ),
                                        component<1, 0>( membraneTrial ),
                                        component<2, 0>( shearTrial ),
                                        component<3, 0>( membraneTrial ),
                                        component<4, 0>( shearTrial ),
                                        component<5, 0>( pinchingTrial ) );
        auto epsShellTestMandel = vec( component<0, 0>( membraneTest ),
                                       component<1, 0>( membraneTest ),
                                       component<2, 0>( shearTest ),
                                       component<3, 0>( membraneTest ),
                                       component<4, 0>( shearTest ),
                                       component<5, 0>( pinchingTest ) );
        auto epsAlphaTrialMandel = sb9W9( alpha, cst( alphaScale ) );
        auto epsAlphaTestMandel = sb9W9( beta, cst( alphaScale ) );
        auto bs1Trial = sb9Bs1( u );
        auto bs1Test = sb9Bs1( v );
        auto bs2Trial = sb9Bs2( u );
        auto bs2Test = sb9Bs2( v );
        auto bs3Trial = sb9Bs3( u );
        auto bs3Test = sb9Bs3( v );
        auto bs4Trial = sb9Bs4( u );
        auto bs4Test = sb9Bs4( v );

        // Hallquist transverse-shear block through Bc1/Bc2, plus the Bs1..Bs4
        // terms from the MATLAB element. The default Bs scale is one: the
        // MATLAB Ds coefficients already include the small 1e-4 factor only on
        // the bending-like Ds4(1:2) terms.
        auto shearStabWeight = cst( shellShearStab * mu * 5.0 / 18.0 );
        double const dNormal = lambda + 2.0 * mu;
        auto invJ011 = component<0, 0>( invJ0 );
        auto invJ021 = component<1, 0>( invJ0 );
        auto invJ022 = component<1, 1>( invJ0 );
        auto invJ033 = component<2, 2>( invJ0 );
        auto ds1 = cst( shellBsStab * dNormal / 3.0 ) * invJ033 * invJ033;
        auto ds2 = ds1;
        auto ds3x = cst( shellBsStab * dNormal / 3.0 ) * invJ011 * invJ011;
        auto ds3y = cst( shellBsStab * dNormal / 3.0 ) * ( invJ021 * invJ021 + invJ022 * invJ022 );
        auto ds4x = cst( shellBsStab * dNormal * 1.0e-4 / 9.0 ) * invJ011 * invJ011;
        auto ds4y = cst( shellBsStab * dNormal * 1.0e-4 / 9.0 ) * ( invJ021 * invJ021 + invJ022 * invJ022 );
        auto ds4z = cst( shellBsStab * dNormal / 9.0 ) * invJ033 * invJ033;

        auto sb9Quad = sb9ThroughThicknessLobatto5();

        ticIf();
        a( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                    _quad=sb9Quad,
                                    _quad1=sb9Quad,
                                    _expr=ddot( C, epsShellTrialMandel, epsShellTestMandel ) );
        tocIf( "sb9.a(0,0).elastic" );

        ticIf();
        a( 0_c, 1_c ) += integrate( _range=elements( mesh ),
                                    _quad=sb9Quad,
                                    _quad1=sb9Quad,
                                    _expr=ddot( C, epsAlphaTrialMandel, epsShellTestMandel ) );
        tocIf( "sb9.a(0,1)" );

        ticIf();
        a( 1_c, 0_c ) += integrate( _range=elements( mesh ),
                                    _quad=sb9Quad,
                                    _quad1=sb9Quad,
                                    _expr=ddot( C, epsShellTrialMandel, epsAlphaTestMandel ) );
        tocIf( "sb9.a(1,0)" );

        ticIf();
        a( 1_c, 1_c ) += integrate( _range=elements( mesh ),
                                    _quad=sb9Quad,
                                    _quad1=sb9Quad,
                                    _expr=ddot( C, epsAlphaTrialMandel, epsAlphaTestMandel ) );
        tocIf( "sb9.a(1,1)" );

        ticIf();
        a( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                    _quad=sb9Quad,
                                    _quad1=sb9Quad,
                                    _expr=shearStabWeight *
                                          inner( sb9Bc1( u ), sb9Bc1( v ) ) );
        tocIf( "sb9.a(0,0).shear-stab13" );

        ticIf();
        a( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                    _quad=sb9Quad,
                                    _quad1=sb9Quad,
                                    _expr=shearStabWeight *
                                          inner( sb9Bc2( u ), sb9Bc2( v ) ) );
        tocIf( "sb9.a(0,0).shear-stab23" );

        ticIf();
        a( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                    _quad=sb9Quad,
                                    _quad1=sb9Quad,
                                    _expr=ds1 *
                                          component<0, 0>( bs1Trial ) *
                                          component<0, 0>( bs1Test ) );
        a( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                    _quad=sb9Quad,
                                    _quad1=sb9Quad,
                                    _expr=ds2 *
                                          component<0, 0>( bs2Trial ) *
                                          component<0, 0>( bs2Test ) );
        tocIf( "sb9.a(0,0).bs12-stab" );

        ticIf();
        a( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                    _quad=sb9Quad,
                                    _quad1=sb9Quad,
                                    _expr=ds3x *
                                          component<0, 0>( bs3Trial ) *
                                          component<0, 0>( bs3Test ) +
                                          ds3y *
                                          component<1, 0>( bs3Trial ) *
                                          component<1, 0>( bs3Test ) );
        tocIf( "sb9.a(0,0).bs3-stab" );

        ticIf();
        a( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                    _quad=sb9Quad,
                                    _quad1=sb9Quad,
                                    _expr=ds4x *
                                          component<0, 0>( bs4Trial ) *
                                          component<0, 0>( bs4Test ) +
                                          ds4y *
                                          component<1, 0>( bs4Trial ) *
                                          component<1, 0>( bs4Test ) +
                                          ds4z *
                                          component<2, 0>( bs4Trial ) *
                                          component<2, 0>( bs4Test ) );
        tocIf( "sb9.a(0,0).bs4-stab" );

        auto l = blockform1( *Xh, strategy, algebraBackend );
        auto f = expr<3, 1>( bodyForceExpression, "f" );
        auto lDisplacement = l( 0_c );

        // In monolithic blockform1, each block assembles into the same global vector.
        // Use += to accumulate block contributions without reinitializing the other blocks.
        // The same pattern is safe in static condensation because the local SB9
        // condenser reads the assembled block rhs stored in the condensed vector.
        ticIf();
        lDisplacement += integrate( _range=elements( mesh ),
                                    _expr=inner( f, v ) );
        tocIf( "sb9.l(0).body-force" );

        for ( auto const& pressureLoad : cfg.pressureLoads )
        {
            double const pressureValue = ( pressureOverride >= 0.0 ) ? pressureOverride : pressureLoad.value;
            ticIf();
            lDisplacement += integrate( _range=markedfaces( mesh, pressureLoad.marker ),
                                        _expr=inner( -cst( pressureValue ) * N(), v ) );
            tocIf( "sb9.l(0).pressure" );
        }

        for ( auto const& tractionLoad : cfg.tractionLoads )
        {
            auto traction = expr<3, 1>( tractionLoad.expression, "traction" );
            ticIf();
            lDisplacement += integrate( _range=markedfaces( mesh, tractionLoad.marker ),
                                        _expr=inner( traction, v ) );
            tocIf( "sb9.l(0).traction" );
        }

        for ( auto const& totalForceLoad : cfg.totalForceLoads )
        {
            double const loadedArea = measure( _range=markedfaces( mesh, totalForceLoad.marker ) );
            CHECK( loadedArea > 0.0 ) << "marked face '" << totalForceLoad.marker
                                      << "' has zero measure, cannot distribute a total force";
            auto traction = expr<3, 1>( shellVectorExpression( totalForceLoad.value / loadedArea ),
                                        "traction" );
            ticIf();
            lDisplacement += integrate( _range=markedfaces( mesh, totalForceLoad.marker ),
                                        _expr=inner( traction, v ) );
            tocIf( "sb9.l(0).total-force" );
        }

        for ( auto const& pointLoad : cfg.pointLoads )
        {
            auto pointRange = markedpoints( mesh, pointLoad.marker );
            double const scale = pointLoadScale( pointRange, pointLoad.quantity, "point load '" + pointLoad.marker + "'" );
            auto pointLoadExpr = expr<3, 1>( shellVectorExpression( scale * pointLoad.value ), "point_load" );
            ticIf();
            lDisplacement += integrate( _range=pointRange,
                                        _expr=inner( pointLoadExpr, id( v ) ) );
            tocIf( "sb9.l(0).point-load" );
        }

        for ( auto const& pointMoment : cfg.pointMoments )
        {
            auto pointRange = markedpoints( mesh, pointMoment.marker );
            double const scale = pointLoadScale( pointRange, pointMoment.quantity, "point moment '" + pointMoment.marker + "'" );
            auto pointMomentExpr = expr<3, 1>( shellVectorExpression( scale * pointMoment.value ), "point_moment" );
            ticIf();
            lDisplacement += integrate( _range=pointRange,
                                        _expr=inner( pointMomentExpr, omega( v ) ) );
            tocIf( "sb9.l(0).point-moment" );
        }

        ticIf();
        l.close();
        tocIf( "sb9.l.close" );

        // Strong displacement Dirichlet conditions act on the whole mixed row
        // block associated with u, not only on the (0,0) matrix sub-block.
        // In monolithic mode the elimination is applied immediately to the
        // assembled matrix. In static-condensation mode the same row operation
        // is deferred and later applied to the condensed displacement system.
        ticIf();
        a.close();
        tocIf( "sb9.a.close" );
        // The RHS can be passed as l(0_c): this block view aliases the
        // displacement part of the global rhs in monolithic mode and carries
        // the same prescribed values for the condensed displacement solve.
        for ( auto const& clampMarker : cfg.clampMarkers )
        {
            ticIf();
            a.row( 0_c ) += on( _range=markedfaces( mesh, clampMarker ),
                                _rhs=l( 0_c ),
                                _element=uBc,
                                _expr=zero<3, 1>(),
                                _type="elimination" );
            tocIf( "sb9.on(0)" );
        }

        auto uBcX = uBc[ComponentType::X];
        auto uBcY = uBc[ComponentType::Y];
        auto uBcZ = uBc[ComponentType::Z];
        auto pointConstraintComponent = [&]( BenchmarkConfig::PointConstraint const& constraint, int c )
        {
            if ( !constraint.components[c] )
            {
                return;
            }

            ticIf();
            switch ( c )
            {
            case 0:
                a.row( 0_c ) += on( _range=markedpoints( mesh, constraint.marker ),
                                    _rhs=l( 0_c ),
                                    _element=uBcX,
                                    _expr=cst( constraint.value( 0 ) ),
                                    _type="elimination" );
                break;
            case 1:
                a.row( 0_c ) += on( _range=markedpoints( mesh, constraint.marker ),
                                    _rhs=l( 0_c ),
                                    _element=uBcY,
                                    _expr=cst( constraint.value( 1 ) ),
                                    _type="elimination" );
                break;
            case 2:
                a.row( 0_c ) += on( _range=markedpoints( mesh, constraint.marker ),
                                    _rhs=l( 0_c ),
                                    _element=uBcZ,
                                    _expr=cst( constraint.value( 2 ) ),
                                    _type="elimination" );
                break;
            default:
                CHECK( false ) << "invalid point constraint component " << c;
            }
            tocIf( "sb9.on(point)" );
        };
        for ( auto const& pointConstraint : cfg.pointConstraints )
        {
            for ( int c = 0; c < 3; ++c )
            {
                pointConstraintComponent( pointConstraint, c );
            }
        }

        if ( !boption( "no-solve" ) )
        {
            ticIf();
            a.solve( _rhs=l, _solution=X,
                     _condense=boption( "sc.condense" ),
                     _condenser=condenser_sb9() );
            tocIf( "sb9.solve" );
        }

        // Re-extract the product-space components after the solve. The
        // sub-elements returned by X(i_c) are value objects, so post-processing
        // must read fresh copies from the solved product element.
        auto uh = X( 0_c );
        auto alphah = X( 1_c );

        auto scalarSpace = Pch<1>( mesh );
        auto thicknessSpace = Pdh<0>( mesh );
        auto symmetricSpace = Pdhms<1>( mesh );
        auto shellThicknessField = vf::project( _space=thicknessSpace, _range=elements( mesh ), _expr=shellThickness() );
        auto shellAreaField = vf::project( _space=thicknessSpace, _range=elements( mesh ), _expr=shellArea0() );
        auto normalDispField = vf::project( _space=scalarSpace,
                                            _range=elements( mesh ),
                                            _expr=inner( idv( uh ), shellNormal() ) );
        auto epsilonh = symmetricSpace->element( "epsilon" );
        auto sigmah = symmetricSpace->element( "sigma" );
        fillSb9SymmetricFields( mesh, uh, alphah, epsilonh, sigmah,
                                lambda, mu, alphaScale, pinchingBpzScale, shellShearFactor );
        if ( boption( "print-matlab-fields" ) )
        {
            uh.printMatlab( "u.m" );
            alphah.printMatlab( "alpha.m" );

            auto fieldReferenceContext = matlabFieldContext( symmetricSpace, cfg, meshData.probe );
            epsilonh.printMatlab( "epsilon_matlab",
                                  fieldReferenceContext,
                                  matlabEpsilonFormat(),
                                  true,
                                  "epsilon" );
            sigmah.printMatlab( "sigma_matlab",
                                fieldReferenceContext,
                                matlabSigmaFormat(),
                                true,
                                "sigma" );
        }

        auto dispNormMax = normLinf( _range=elements( mesh ), _pset=_Q<2>(), _expr=norm2( idv( uh ) ) );
        auto probeCtx = Uh->context();
        probeCtx.add( toNode( meshData.probe ) );
        auto const probeValue = uh.evaluate( probeCtx, false );
        CHECK( probeValue.size() >= 3 )
            << "expected a 3-component displacement probe value";
        double const probeDisplacement =
            cfg.probeDirection( 0 ) * probeValue( 0 ) +
            cfg.probeDirection( 1 ) * probeValue( 1 ) +
            cfg.probeDirection( 2 ) * probeValue( 2 );

        std::cout << "max displacement norm = " << dispNormMax.value()
                  << " at " << dispNormMax.arg().transpose() << "\n"
                  << cfg.probeLabel << " = " << probeDisplacement << "\n";
        if ( cfg.reference.hasProbeValue )
        {
            double const absoluteError = probeDisplacement - cfg.reference.probeValue;
            double const relativeError = ( std::abs( cfg.reference.probeValue ) > 0.0 ) ?
                                             std::abs( absoluteError / cfg.reference.probeValue ) :
                                             std::abs( absoluteError );
            std::cout << "reference " << cfg.probeLabel << " = " << cfg.reference.probeValue;
            if ( !cfg.reference.description.empty() )
            {
                std::cout << " (" << cfg.reference.description << ")";
            }
            std::cout << "\n"
                      << "absolute error = " << absoluteError << "\n"
                      << "relative error = " << relativeError << "\n";

            if ( boption( "check-reference" ) &&
                 ( cfg.reference.hasAbsoluteTolerance || cfg.reference.hasRelativeTolerance ) )
            {
                bool referenceOk = true;
                if ( cfg.reference.hasAbsoluteTolerance )
                {
                    std::cout << "absolute tolerance = " << cfg.reference.absoluteTolerance << "\n";
                    referenceOk = referenceOk && ( std::abs( absoluteError ) <= cfg.reference.absoluteTolerance );
                }
                if ( cfg.reference.hasRelativeTolerance )
                {
                    std::cout << "relative tolerance = " << cfg.reference.relativeTolerance << "\n";
                    referenceOk = referenceOk && ( relativeError <= cfg.reference.relativeTolerance );
                }
                std::cout << "reference check = " << ( referenceOk ? "PASS" : "FAIL" ) << "\n";
                if ( !referenceOk )
                {
                    throw std::runtime_error( "benchmark reference tolerance check failed" );
                }
            }
        }

        bool const fieldReferenceOk = checkFieldReferences( cfg, symmetricSpace, epsilonh, sigmah );
        if ( !fieldReferenceOk )
            throw std::runtime_error( "benchmark field reference tolerance check failed" );

        auto e = exporter( _mesh=mesh );
        e->addRegions();
        e->add( "u", uh );
        e->add( "alpha", alphah );
        e->add( "epsilon", epsilonh );
        e->add( "sigma", sigmah );
        if ( boption( "export-thickness" ) )
        {
            e->add( "shell_thickness", shellThicknessField );
            e->add( "shell_area0", shellAreaField );
        }
        if ( boption( "export-normal-displacement" ) )
        {
            e->add( "u_normal", normalDispField );
        }
        e->save();
    }
    catch ( ... )
    {
        handleExceptions();
        return 1;
    }

    return 0;
}
