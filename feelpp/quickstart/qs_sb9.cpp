/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhm.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/geo.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelmesh/hypercube.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/sb9_bending.hpp>
#include <feel/feelvf/sb9_pinching.hpp>
#include <feel/feelvf/sb9_quadrature.hpp>
#include <feel/feelvf/sb9_shear.hpp>
#include <feel/feelvf/sb9_strain.hpp>
#include <feel/feelvf/sb9_stabilization.hpp>
#include <feel/feelvf/shellgeometric.hpp>
#include <feel/feelvf/vf.hpp>

#include "qs_elasticity_case.hpp"
#include "qs_elasticity_checks.hpp"
#include "qs_sb9_postprocess.hpp"

#include <boost/format.hpp>

#include <iostream>
#include <set>
#include <stdexcept>
#include <string>

using namespace Feel;

namespace
{
namespace qsec = Feel::Quickstart::ElasticityChecks;
namespace qsecase = Feel::Quickstart::ElasticityCase;
namespace qssb9 = Feel::Quickstart::SB9;

using mesh_type = Mesh<Hypercube<3>>;
using PointConstraintConfig = qsecase::PointConstraintConfig;

po::options_description
makeOptions()
{
    po::options_description options( "qs_sb9 options" );
    qsecase::addOptions( options, "{0,0,0}", "{0,0,0}" );
    options.add_options()
        ( "load", po::value<double>()->default_value( -1.0 ),
          "fallback constant z traction on XPlus when no elasticity.case is provided" )
        ( "alpha-scale", po::value<double>()->default_value( 1.0 ),
          "scale applied to the shell-normal SB9 scalar correction" )
        ( "pinching-bpz-scale", po::value<double>()->default_value( 1.0 ),
          "scale applied to the zeta*Bpz pinching term; 1 enables the full SB9 pinching formulation" )
        ( "shell-shear-factor", po::value<double>()->default_value( 1.25 ),
          "Hallquist-like transverse shear shape factor coefficient" )
        ( "shell-shear-stab", po::value<double>()->default_value( 2.0 ),
          "scale applied to the Hallquist Bc1/Bc2 stabilization block" )
        ( "shell-bs-stab", po::value<double>()->default_value( 1.0 ),
          "scale applied to the additional SB9 Bs1..Bs4 stabilization block" )
        ( "shell-membrane-stab", po::value<double>()->default_value( 1.0 ),
          "membrane stabilization coefficient applied to the Bs3 block" )
        ( "shell-bending-stab", po::value<double>()->default_value( 1.0e-4 ),
          "bending stabilization coefficient applied to the in-plane Bs4 block" )
        ( "shell-pinching-stab", po::value<double>()->default_value( 1.0 ),
          "pinching stabilization coefficient applied to the Bs1/Bs2 and normal Bs4 blocks" )
        ( "monolithic", po::value<bool>()->default_value( false ),
          "use the monolithic mixed solve instead of the default SB9 static condensation" )
        ( "check-reference", po::value<bool>()->default_value( true ),
          "fail when a JSON reference value exceeds its configured tolerance" )
        ( "print-matlab-fields", po::value<bool>()->default_value( false ),
          "write u, alpha, epsilon, and sigma fields; tensor fields use the MATLAB validation order" )
        ( "export-thickness", po::value<bool>()->default_value( true ),
          "export shell thickness and mid-surface area diagnostics" )
        ( "export-normal-displacement", po::value<bool>()->default_value( true ),
          "export inner(u,shellNormal()) projected to a scalar field" )
        ( "show-timings", po::value<bool>()->default_value( false ),
          "print coarse assembly, solve, and postprocessing timings" );
    return options.add( feel_options() );
}

AboutData
makeAbout()
{
    AboutData about( "qs_sb9",
                     "qs_sb9",
                     "0.1",
                     "SB9 mixed shell formulation",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) Feel++ Consortium" );
    about.addAuthor( "Feel++ Consortium", "developer", "feelpp-devel@feelpp.org", "" );
    return about;
}

std::string
unitShellGeo()
{
    return R"(
Mesh.RecombineAll = 1;

Point(1) = {0, 0, -0.05, 1};
Point(2) = {1, 0, -0.05, 1};
Point(3) = {0, 1, -0.05, 1};
Point(4) = {1, 1, -0.05, 1};
Point(5) = {0, 0,  0.05, 1};
Point(6) = {1, 0,  0.05, 1};
Point(7) = {0, 1,  0.05, 1};
Point(8) = {1, 1,  0.05, 1};

Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};
Line(5) = {5, 6};
Line(6) = {6, 8};
Line(7) = {8, 7};
Line(8) = {7, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {4, 8};
Line(12) = {3, 7};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Line Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};
Line Loop(3) = {1, 10, -5, -9};
Plane Surface(3) = {3};
Line Loop(4) = {2, 11, -6, -10};
Plane Surface(4) = {4};
Line Loop(5) = {-3, 11, 7, -12};
Plane Surface(5) = {5};
Line Loop(6) = {-4, 12, 8, -9};
Plane Surface(6) = {6};

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

Transfinite Line {1, 3, 5, 7} = 2;
Transfinite Line {2, 4, 6, 8} = 2;
Transfinite Line {9, 10, 11, 12} = 2;
Transfinite Surface {1} = {1, 2, 4, 3};
Transfinite Surface {2} = {5, 6, 8, 7};
Transfinite Surface {3} = {1, 2, 6, 5};
Transfinite Surface {4} = {2, 4, 8, 6};
Transfinite Surface {5} = {3, 4, 8, 7};
Transfinite Surface {6} = {1, 3, 7, 5};
Transfinite Volume {1} = {1, 2, 4, 3, 5, 6, 8, 7};
Recombine Surface {1, 2, 3, 4, 5, 6};
Recombine Volume {1};

Physical Surface("ZMoins") = {1};
Physical Surface("ZPlus") = {2};
Physical Surface("YMoins") = {3};
Physical Surface("XPlus") = {4};
Physical Surface("YPlus") = {5};
Physical Surface("XMoins") = {6};
Physical Volume("Shell") = {1};
)";
}

std::shared_ptr<mesh_type>
createUnitShellMesh()
{
    Environment::changeRepository( _directory=boost::format( "quickstart/%1%/" ) % Environment::about().appName() );
    return createGMSHMesh( _mesh=new mesh_type,
                           _desc=geo( _filename="qs_sb9_unit_shell.geo",
                                      _desc=unitShellGeo(),
                                      _dim=3,
                                      _order=1,
                                      _h=1.0 ),
                           _force_rebuild=boption( "gmsh.rebuild" ) );
}

std::shared_ptr<mesh_type>
createMesh( qsecase::Config const& cfg )
{
    if ( cfg.meshFilename.empty() )
        return createUnitShellMesh();
    return loadMesh( _mesh=new mesh_type, _filename=cfg.meshFilename );
}

bool
hasElasticityCase()
{
    return qsecase::optionExplicitlySet( "elasticity.case" ) && !soption( "elasticity.case" ).empty();
}

qsecase::Config
caseConfigFromInput()
{
    auto cfg = qsecase::fromEnvironment( "{0,0,0}", "{0,0,0}" );
    if ( !hasElasticityCase() )
    {
        cfg.dirichletMarkers = { "XMoins" };
        cfg.faceTractions.push_back(
            { { "XPlus" },
              "{0,0," + qsecase::formatDouble( doption( "load" ) ) + "}" } );
    }
    cfg.referenceChecks.target =
        qssb9::preferredTarget( cfg.referenceChecks,
                                qsecase::optionExplicitlySet( "checks.target" ) &&
                                !soption( "checks.target" ).empty() );
    return cfg;
}

template<typename PointRangeType>
double
pointLoadScale( PointRangeType const& pointRange, std::string const& quantity, std::string const& label )
{
    if ( quantity == "per-point" || quantity == "per_point" || quantity == "point" )
        return 1.0;
    if ( quantity != "total" )
        throw std::invalid_argument( label + " quantity must be 'per-point' or 'total'" );

    size_type const nMarkedPoints = nelements( pointRange, true );
    if ( nMarkedPoints == 0 )
        throw std::invalid_argument( label + " uses quantity=total but selects no marked point" );
    return 1.0 / nMarkedPoints;
}

qsec::point_type<3>
fallbackPostprocessPoint( qsec::Config<3> const& checks )
{
    if ( !checks.fieldReferences.empty() )
        return checks.fieldReferences.front().point;
    for ( auto const& probe : checks.probes )
        if ( probe.hasProbe )
            return probe.point;
    qsec::point_type<3> point;
    point << 0.5, 0.5, 0.0;
    return point;
}

void
printMeshInfo( std::shared_ptr<mesh_type> const& mesh )
{
    if ( !Environment::isMasterRank() )
        return;

    std::cout << "mesh: elements=" << mesh->numGlobalElements()
              << ", faces=" << mesh->numGlobalFaces()
              << ", nodes=" << mesh->numGlobalPoints() << "\n";
}

template<typename UhType, typename AhType, typename XhType>
void
printUnknownSpaceInfo( UhType const& Uh, AhType const& Ah, XhType const& Xh )
{
    if ( !Environment::isMasterRank() )
        return;

    std::cout << "spaces: Uh=Pchv<1> ndof=" << Uh->nDof()
              << ", Ah=Pdh<0> ndof=" << Ah->nDof()
              << ", Xh=Uh x Ah ndof=" << Xh->nDof() << "\n";
}

template<typename ScalarSpaceType, typename CellSpaceType, typename TensorSpaceType>
void
printPostprocessSpaceInfo( ScalarSpaceType const& scalarSpace,
                           CellSpaceType const& cellSpace,
                           TensorSpaceType const& tensorSpace )
{
    if ( !Environment::isMasterRank() )
        return;

    std::cout << "postprocess spaces: Pch<1> ndof=" << scalarSpace->nDof()
              << ", Pdh<0> ndof=" << cellSpace->nDof()
              << ", Pdhms<1> ndof=" << tensorSpace->nDof() << "\n";
}
} // namespace

int
main( int argc, char** argv )
{
    using namespace vf;

    try
    {
        Environment env( _argc=argc, _argv=argv, _desc=makeOptions(), _about=makeAbout() );

        auto cfg = caseConfigFromInput();
        auto mesh = createMesh( cfg );
        printMeshInfo( mesh );
        if ( cfg.expectedElements )
        {
            auto const nElements = nelements( elements( mesh ), true );
            if ( nElements != *cfg.expectedElements )
                throw std::runtime_error( "mesh element count check failed: expected " +
                                          std::to_string( *cfg.expectedElements ) +
                                          ", got " + std::to_string( nElements ) );
        }

        double const E = cfg.E;
        double const nu = cfg.nu;
        double const lambda = E * nu / ( ( 1.0 + nu ) * ( 1.0 - 2.0 * nu ) );
        double const mu = E / ( 2.0 * ( 1.0 + nu ) );
        double const alphaScale = doption( "alpha-scale" );
        double const pinchingBpzScale = doption( "pinching-bpz-scale" );
        double const shellShearFactor = doption( "shell-shear-factor" );
        double const shellShearStab = doption( "shell-shear-stab" );
        double const shellBsStab = doption( "shell-bs-stab" );
        double const shellMembraneStab = doption( "shell-membrane-stab" );
        double const shellBendingStab = doption( "shell-bending-stab" );
        double const shellPinchingStab = doption( "shell-pinching-stab" );
        bool const useStaticCondensation = !boption( "monolithic" );
        bool const solved = !boption( "no-solve" );

        std::cout << "qs_sb9: E=" << E
                  << ", nu=" << nu
                  << ", pinching-bpz-scale=" << pinchingBpzScale
                  << ", strategy=" << ( useStaticCondensation ? "static_condensation" : "monolithic" )
                  << ", check-target=" << cfg.referenceChecks.target << "\n";

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
        printUnknownSpaceInfo( Uh, Ah, Xh );
        tocIf( "sb9.spaces" );

        solve::strategy const strategy = useStaticCondensation
                                             ? solve::strategy::static_condensation
                                             : solve::strategy::monolithic;
        auto a = blockform2( *Xh, strategy, algebraBackend );
        auto l = blockform1( *Xh, strategy, algebraBackend );

        auto u = trial( Uh, "u" );
        auto v = test( Uh, "v" );
        auto alpha = trial( Ah, "alpha" );
        auto beta = test( Ah, "beta" );
        auto X = Xh->element();
        auto uBc = X( 0_c );

        auto zt = zeta();
        auto C = isotropic_stiffness<3>( lambda, mu );
        auto shellShearWeight = cst( shellShearFactor ) * ( cst( 1.0 ) - zt * zt );

        auto pinchingScale = cst( pinchingBpzScale );
        auto epsShellTrialMandel = sb9ShellStrain( u, zt, shellShearWeight, pinchingScale );
        auto epsShellTestMandel = sb9ShellStrain( v, zt, shellShearWeight, pinchingScale );
        auto alphaTrial = sb9W9( alpha, cst( alphaScale ) );
        auto alphaTest = sb9W9( beta, cst( alphaScale ) );

        auto sb9Quad = sb9ThroughThicknessLobatto5();

        ticIf();
        a( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                    _quad=sb9Quad,
                                    _expr=ddot( C, epsShellTrialMandel, epsShellTestMandel ) );
        a( 0_c, 1_c ) += integrate( _range=elements( mesh ),
                                    _quad=sb9Quad,
                                    _expr=ddot( C, alphaTrial, epsShellTestMandel ) );
        a( 1_c, 0_c ) += integrate( _range=elements( mesh ),
                                    _quad=sb9Quad,
                                    _expr=ddot( C, epsShellTrialMandel, alphaTest ) );
        a( 1_c, 1_c ) += integrate( _range=elements( mesh ),
                                    _quad=sb9Quad,
                                    _expr=ddot( C, alphaTrial, alphaTest ) );
        tocIf( "sb9.elastic" );

        auto sb9Stabilization = [&]( auto const& trialField, auto const& testField )
        {
            auto invJ0 = inv( shellJacobian0() );
            auto j11 = component<0, 0>( invJ0 );
            auto j21 = component<1, 0>( invJ0 );
            auto j22 = component<1, 1>( invJ0 );
            auto j33 = component<2, 2>( invJ0 );

            auto normalStiffness = cst( shellBsStab * ( lambda + 2.0 * mu ) );
            auto dz = normalStiffness * j33 * j33;
            auto dx = normalStiffness * j11 * j11;
            auto dy = normalStiffness * ( j21 * j21 + j22 * j22 );

            auto c0 = []( auto const& Bu, auto const& Bv )
            {
                return component<0, 0>( Bu ) * component<0, 0>( Bv );
            };
            auto c1 = []( auto const& Bu, auto const& Bv )
            {
                return component<1, 0>( Bu ) * component<1, 0>( Bv );
            };
            auto c2 = []( auto const& Bu, auto const& Bv )
            {
                return component<2, 0>( Bu ) * component<2, 0>( Bv );
            };

            auto bs1u = sb9Bs1( trialField );
            auto bs1v = sb9Bs1( testField );
            auto bs2u = sb9Bs2( trialField );
            auto bs2v = sb9Bs2( testField );
            auto bs3u = sb9Bs3( trialField );
            auto bs3v = sb9Bs3( testField );
            auto bs4u = sb9Bs4( trialField );
            auto bs4v = sb9Bs4( testField );

            auto shearPart = cst( shellShearStab * mu * 5.0 / 18.0 ) *
                             ( inner( sb9Bc1( trialField ), sb9Bc1( testField ) ) +
                               inner( sb9Bc2( trialField ), sb9Bc2( testField ) ) );

            auto bsPart = cst( shellPinchingStab / 3.0 ) * dz *
                              ( c0( bs1u, bs1v ) + c0( bs2u, bs2v ) ) +
                          cst( shellMembraneStab / 3.0 ) *
                              ( dx * c0( bs3u, bs3v ) + dy * c1( bs3u, bs3v ) ) +
                          cst( 1.0 / 9.0 ) *
                              ( cst( shellBendingStab ) *
                                    ( dx * c0( bs4u, bs4v ) + dy * c1( bs4u, bs4v ) ) +
                                cst( shellPinchingStab ) * dz * c2( bs4u, bs4v ) );

            return shearPart + bsPart;
        };

        ticIf();
        a( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                    _quad=sb9Quad,
                                    _expr=sb9Stabilization( u, v ) );
        tocIf( "sb9.stabilization" );

        auto lDisplacement = l( 0_c );
        auto f = expr<3, 1>( cfg.bodyForceExpression, "f" );
        lDisplacement += integrate( _range=elements( mesh ),
                                    _expr=inner( f, v ) );

        for ( auto const& load : cfg.facePressures )
        {
            auto pressure = expr( load.expression, "pressure" );
            lDisplacement += integrate( _range=markedfaces( mesh, load.marker ),
                                        _expr=inner( -pressure * N(), v ) );
        }

        for ( auto const& load : cfg.faceTractions )
        {
            auto traction = expr<3, 1>( load.expression, "traction" );
            lDisplacement += integrate( _range=markedfaces( mesh, load.markers ),
                                        _expr=inner( traction, v ) );
        }

        for ( auto const& load : cfg.faceTotalForces )
        {
            auto forceRange = markedfaces( mesh, load.markers );
            double const area = integrate( _range=forceRange, _expr=cst( 1.0 ) ).evaluate()( 0, 0 );
            if ( area <= 0.0 )
                throw std::invalid_argument( "face total force selects no marked face" );

            auto totalForce = expr<3, 1>( load.expression, "face_total_force" );
            lDisplacement += integrate( _range=forceRange,
                                        _expr=inner( cst( 1.0/area ) * totalForce, v ) );
        }

        for ( auto const& load : cfg.pointForces )
        {
            auto pointRange = markedpoints( mesh, load.markers );
            double const scale = pointLoadScale( pointRange, load.quantity, "point force" );
            auto pointForce = expr<3, 1>( load.expression, "point_force" );
            lDisplacement += integrate( _range=pointRange,
                                        _expr=cst( scale ) * inner( pointForce, id( v ) ) );
        }

        for ( auto const& load : cfg.pointMoments )
        {
            auto pointRange = markedpoints( mesh, load.markers );
            double const scale = pointLoadScale( pointRange, load.quantity, "point moment" );
            auto pointMoment = expr<3, 1>( load.expression, "point_moment" );
            lDisplacement += integrate( _range=pointRange,
                                        _expr=cst( scale ) * inner( pointMoment, omega( v ) ) );
        }

        l.close();
        a.close();

        auto g = expr<3, 1>( cfg.dirichletExpression, "g" );
        for ( auto const& clampMarker : cfg.dirichletMarkers )
        {
            a.row( 0_c ) += on( _range=markedfaces( mesh, clampMarker ),
                                _rhs=l( 0_c ),
                                _element=uBc,
                                _expr=g,
                                _type="elimination" );
        }

        auto uBcX = uBc[ComponentType::X];
        auto uBcY = uBc[ComponentType::Y];
        auto uBcZ = uBc[ComponentType::Z];
        auto applyPointConstraintComponent = [&]( PointConstraintConfig const& constraint, int component )
        {
            if ( !constraint.components[component] )
                return;

            auto value = expr( constraint.values[component], "point_constraint" );
            switch ( component )
            {
            case 0:
                a.row( 0_c ) += on( _range=markedpoints( mesh, constraint.marker ),
                                    _rhs=l( 0_c ),
                                    _element=uBcX,
                                    _expr=value,
                                    _type="elimination" );
                break;
            case 1:
                a.row( 0_c ) += on( _range=markedpoints( mesh, constraint.marker ),
                                    _rhs=l( 0_c ),
                                    _element=uBcY,
                                    _expr=value,
                                    _type="elimination" );
                break;
            case 2:
                a.row( 0_c ) += on( _range=markedpoints( mesh, constraint.marker ),
                                    _rhs=l( 0_c ),
                                    _element=uBcZ,
                                    _expr=value,
                                    _type="elimination" );
                break;
            default:
                throw std::invalid_argument( "invalid point constraint component" );
            }
        };
        for ( auto const& constraint : cfg.pointConstraints )
            for ( int component = 0; component < 3; ++component )
                applyPointConstraintComponent( constraint, component );

        auto applyFaceConstraintComponent = [&]( PointConstraintConfig const& constraint, int component )
        {
            if ( !constraint.components[component] )
                return;

            auto value = expr( constraint.values[component], "face_constraint" );
            switch ( component )
            {
            case 0:
                a.row( 0_c ) += on( _range=markedfaces( mesh, constraint.marker ),
                                    _rhs=l( 0_c ),
                                    _element=uBcX,
                                    _expr=value,
                                    _type="elimination" );
                break;
            case 1:
                a.row( 0_c ) += on( _range=markedfaces( mesh, constraint.marker ),
                                    _rhs=l( 0_c ),
                                    _element=uBcY,
                                    _expr=value,
                                    _type="elimination" );
                break;
            case 2:
                a.row( 0_c ) += on( _range=markedfaces( mesh, constraint.marker ),
                                    _rhs=l( 0_c ),
                                    _element=uBcZ,
                                    _expr=value,
                                    _type="elimination" );
                break;
            default:
                throw std::invalid_argument( "invalid face constraint component" );
            }
        };
        for ( auto const& constraint : cfg.faceConstraints )
            for ( int component = 0; component < 3; ++component )
                applyFaceConstraintComponent( constraint, component );

        if ( solved )
        {
            ticIf();
            a.solve( _rhs=l, _solution=X,
                     _condense=useStaticCondensation,
                     _condenser=condenser_sb9() );
            tocIf( "sb9.solve" );
        }

        auto uh = X( 0_c );
        auto alphah = X( 1_c );

        ticIf();
        auto scalarSpace = Pch<1>( mesh );
        auto cellSpace = Pdh<0>( mesh );
        auto tensorSpace = Pdhms<1>( mesh );
        printPostprocessSpaceInfo( scalarSpace, cellSpace, tensorSpace );
        auto thickness = vf::project( _space=cellSpace, _range=elements( mesh ), _expr=shellThickness() );
        auto area0 = vf::project( _space=cellSpace, _range=elements( mesh ), _expr=shellArea0() );
        auto normalDisplacement = vf::project( _space=scalarSpace,
                                               _range=elements( mesh ),
                                               _expr=inner( idv( uh ), shellNormal() ) );
        auto epsilon = tensorSpace->element( "epsilon" );
        auto sigma = tensorSpace->element( "sigma" );
        qssb9::fillSymmetricFields( mesh, uh, alphah, epsilon, sigma,
                                    lambda, mu, alphaScale, pinchingBpzScale, shellShearFactor );
        tocIf( "sb9.postprocess" );

        if ( boption( "print-matlab-fields" ) )
        {
            uh.printMatlab( "u.m" );
            alphah.printMatlab( "alpha.m" );

            auto fieldContext = qssb9::matlabFieldContext( tensorSpace,
                                                           cfg.referenceChecks,
                                                           fallbackPostprocessPoint( cfg.referenceChecks ) );
            epsilon.printMatlab( "epsilon_matlab",
                                 fieldContext,
                                 qssb9::matlabEpsilonFormat(),
                                 true,
                                 "epsilon" );
            sigma.printMatlab( "sigma_matlab",
                               fieldContext,
                               qssb9::matlabSigmaFormat(),
                               true,
                               "sigma" );
        }

        auto dispNormMax = normLinf( _range=elements( mesh ), _pset=_Q<2>(), _expr=norm2( idv( uh ) ) );
        std::cout << "max displacement norm = " << dispNormMax.value()
                  << " at " << dispNormMax.arg().transpose() << "\n";

        bool const runReferenceChecks = solved && boption( "check-reference" );
        qsec::ElasticityReferenceChecker<3> checks( "qs_sb9" );
        checks.setConfig( cfg.referenceChecks )
              .setTarget( cfg.referenceChecks.target )
              .addDisplacementProbe( Uh, uh, runReferenceChecks )
              .addCantileverChecks( mesh, Uh, uh, E, runReferenceChecks );
        checks.add( "SB9 field reference values",
                    runReferenceChecks && cfg.referenceChecks.hasFieldReferences(),
                    [&cfg,&tensorSpace,&epsilon,&sigma]() {
                        return qssb9::checkFieldReferences( cfg.referenceChecks,
                                                            tensorSpace,
                                                            epsilon,
                                                            sigma,
                                                            cfg.referenceChecks.target,
                                                            true );
                    } );

        auto e = exporter( _mesh=mesh );
        e->addRegions();
        e->add( "u", uh );
        e->add( "alpha", alphah );
        e->add( "epsilon", epsilon );
        e->add( "sigma", sigma );
        if ( boption( "export-thickness" ) )
        {
            e->add( "shell_thickness", thickness );
            e->add( "shell_area0", area0 );
        }
        if ( boption( "export-normal-displacement" ) )
            e->add( "u_normal", normalDisplacement );
        e->save();

        return checks.run();
    }
    catch ( ... )
    {
        handleExceptions();
    }

    return 1;
}
