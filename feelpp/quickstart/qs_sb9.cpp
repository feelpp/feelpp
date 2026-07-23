/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/tensorformat.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhm.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/geo.hpp>
#include <feel/feelmesh/hypercube.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/sb9_bending.hpp>
#include <feel/feelvf/sb9_pinching.hpp>
#include <feel/feelvf/sb9_quadrature.hpp>
#include <feel/feelvf/sb9_shear.hpp>
#include <feel/feelvf/sb9_strain.hpp>
#include <feel/feelvf/sb9_stabilization.hpp>
#include <feel/feelvf/vf.hpp>

#include <boost/format.hpp>

using namespace Feel;

namespace
{
using mesh_type = Mesh<Hypercube<3>>;

po::options_description
makeOptions()
{
    po::options_description options( "qs_sb9 options" );
    options.add_options()
        ( "E", po::value<double>()->default_value( 10.0 ), "Young modulus" )
        ( "nu", po::value<double>()->default_value( 0.3 ), "Poisson ratio" )
        ( "load", po::value<double>()->default_value( -1.0 ), "constant z traction on XPlus" )
        ( "alpha-scale", po::value<double>()->default_value( 1.0 ), "scale applied to the internal SB9 scalar mode" )
        ( "monolithic", po::value<bool>()->default_value( false ), "use the monolithic mixed solve instead of SB9 static condensation" )
        ( "print-matlab-fields", po::value<bool>()->default_value( false ),
          "print epsilon/sigma tensor fields in MATLAB validation order" )
        ( "no-solve", po::value<bool>()->default_value( false ), "assemble only" );
    return options.add( feel_options() );
}

AboutData
makeAbout()
{
    AboutData about( "qs_sb9",
                     "qs_sb9",
                     "0.1",
                     "Minimal SB9 mixed shell formulation",
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
} // namespace

int
main( int argc, char** argv )
{
    using namespace vf;

    Environment env( _argc=argc, _argv=argv, _desc=makeOptions(), _about=makeAbout() );

    auto mesh = createUnitShellMesh();

    double const E = doption( "E" );
    double const nu = doption( "nu" );
    double const lambda = E * nu / ( ( 1.0 + nu ) * ( 1.0 - 2.0 * nu ) );
    double const mu = E / ( 2.0 * ( 1.0 + nu ) );
    double const alphaScale = doption( "alpha-scale" );

    auto Uh = Pchv<1>( mesh );
    auto Ah = Pdh<0>( mesh );
    auto Xh = productPtr( Uh, Ah );
    auto Th = Pdhms<1>( mesh );

    auto epsilon = Th->element( "epsilon" );
    auto sigma = Th->element( "sigma" );
    [[maybe_unused]] auto epsilonXY = epsilon.tensorComponent( Component::X, Component::Y );
    [[maybe_unused]] auto sigmaZZ = sigma.tensorComponent( Component::Z, Component::Z );

    auto u = trial( Uh, "u" );
    auto v = test( Uh, "v" );
    auto alpha = trial( Ah, "alpha" );
    auto beta = test( Ah, "beta" );

    auto X = Xh->element();
    auto uBc = X( 0_c );

    auto backend = Feel::backend( _rebuild=true, _worldcomm=Uh->worldCommPtr() );
    bool const useStaticCondensation = !boption( "monolithic" );
    solve::strategy const strategy = useStaticCondensation
                                         ? solve::strategy::static_condensation
                                         : solve::strategy::monolithic;
    auto a = blockform2( *Xh, strategy, backend );
    auto l = blockform1( *Xh, strategy, backend );

    auto zt = zeta();
    auto C = isotropic_stiffness<3>( lambda, mu );
    auto shearWeight = cst( 5.0 / 4.0 ) * ( cst( 1.0 ) - zt * zt );

    auto epsU = sb9MembraneBending( u, zt );
    auto epsV = sb9MembraneBending( v, zt );
    auto epsP = sb9Pinching( u, zt );
    auto epsQ = sb9Pinching( v, zt );
    auto epsS = sb9Shear( u, shearWeight );
    auto epsT = sb9Shear( v, shearWeight );
    auto epsW = sb9W9( alpha, cst( alphaScale ) );
    auto epsZ = sb9W9( beta, cst( alphaScale ) );

    auto epsShellTrial = vec( component<0, 0>( epsU ),
                              component<1, 0>( epsU ),
                              component<2, 0>( epsS ),
                              component<3, 0>( epsU ),
                              component<4, 0>( epsS ),
                              component<5, 0>( epsP ) );
    auto epsShellTest = vec( component<0, 0>( epsV ),
                             component<1, 0>( epsV ),
                             component<2, 0>( epsT ),
                             component<3, 0>( epsV ),
                             component<4, 0>( epsT ),
                             component<5, 0>( epsQ ) );

    // SB9 assembly uses Feel++ compact symmetric storage in Mandel scaling:
    //   Storage + Mandel = xx,xy,xz,yy,yz,zz with sqrt(2)-scaled shear slots.
    // MATLAB validation tables use:
    //   epsilon: xx,yy,zz,xy,xz,yz with engineering shear components,
    //   sigma:   xx,yy,zz,xy,xz,yz with tensor shear components.
    [[maybe_unused]] SymmetricTensorFormat const sb9AssemblyFormat{
        SymmetricTensorOrder::Storage, SymmetricTensorScaling::Mandel };
    [[maybe_unused]] SymmetricTensorFormat const matlabEpsilonFormat{
        SymmetricTensorOrder::DiagonalFirst, SymmetricTensorScaling::EngineeringShear };
    [[maybe_unused]] SymmetricTensorFormat const matlabSigmaFormat{
        SymmetricTensorOrder::DiagonalFirst, SymmetricTensorScaling::Tensor };

    // The tensor fields above are the semantic storage target for postprocessed
    // SB9 strain/stress values. Use tensorComponent(i,j) to assign physical
    // tensor entries, then evaluate/print with matlabEpsilonFormat or
    // matlabSigmaFormat when comparing against the MATLAB validation order.

    auto sb9Quad = sb9ThroughThicknessLobatto5();

    // Minimal SB9 mixed elastic formulation:
    //
    // [ u     ]  displacement Q1 vector field, 24 element dofs
    // [ alpha ]  internal P0 scalar mode, one element dof
    //
    a( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                _quad=sb9Quad,
                                _expr=ddot( C, epsShellTrial, epsShellTest ) );
    a( 0_c, 1_c ) += integrate( _range=elements( mesh ),
                                _quad=sb9Quad,
                                _expr=ddot( C, epsW, epsShellTest ) );
    a( 1_c, 0_c ) += integrate( _range=elements( mesh ),
                                _quad=sb9Quad,
                                _expr=ddot( C, epsShellTrial, epsZ ) );
    a( 1_c, 1_c ) += integrate( _range=elements( mesh ),
                                _quad=sb9Quad,
                                _expr=ddot( C, epsW, epsZ ) );

#if 0
    // Transient extension: mass belongs only to the physical Q1 displacement
    // field. The internal SB9 scalar alpha is an assumed-strain/static
    // condensation variable and is intentionally not included in inertia.
    double constexpr rho = 1.0;
    auto massQuad = sb9LumpedMassLobatto();
    auto m = form2( _trial=Uh, _test=Uh );
    m = integrate( _range=elements( mesh ),
                   _quad=massQuad,
                   _expr=cst( rho ) * inner( u, v ) );
#endif

    // Compact equivalent of the stabilization block used in
    // qs_sb9_mixed_bending.cpp.
    auto sb9Stabilization = [&]( auto const& trialField, auto const& testField )
    {
        double constexpr shearStabilizationScale = 1.0;
        double constexpr bsStabilizationScale = 1.0;

        auto invJ0 = inv( shellJacobian0() );
        auto j11 = component<0, 0>( invJ0 );
        auto j21 = component<1, 0>( invJ0 );
        auto j22 = component<1, 1>( invJ0 );
        auto j33 = component<2, 2>( invJ0 );

        auto normalStiffness = cst( bsStabilizationScale * ( lambda + 2.0 * mu ) );
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

        auto shearPart = cst( shearStabilizationScale * mu * 5.0 / 18.0 ) *
                         ( inner( sb9Bc1( trialField ), sb9Bc1( testField ) ) +
                           inner( sb9Bc2( trialField ), sb9Bc2( testField ) ) );

        auto bsPart = cst( 1.0 / 3.0 ) * dz * ( c0( bs1u, bs1v ) + c0( bs2u, bs2v ) ) +
                      cst( 1.0 / 3.0 ) * ( dx * c0( bs3u, bs3v ) + dy * c1( bs3u, bs3v ) ) +
                      cst( 1.0 / 9.0 ) * ( cst( 1.0e-4 ) *
                                           ( dx * c0( bs4u, bs4v ) + dy * c1( bs4u, bs4v ) ) +
                                           dz * c2( bs4u, bs4v ) );

        return shearPart + bsPart;
    };

    a( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                _quad=sb9Quad,
                                _expr=sb9Stabilization( u, v ) );

    l( 0_c ) += integrate( _range=markedfaces( mesh, "XPlus" ),
                           _expr=inner( vec( cst( 0.0 ), cst( 0.0 ), cst( doption( "load" ) ) ), v ) );

    l.close();
    a.close();

    a.row( 0_c ) += on( _range=markedfaces( mesh, "XMoins" ),
                        _rhs=l( 0_c ),
                        _element=uBc,
                        _expr=zero<3, 1>(),
                        _type="elimination" );

    if ( !boption( "no-solve" ) )
        a.solve( _rhs=l, _solution=X,
                 _condense=useStaticCondensation,
                 _condenser=condenser_sb9() );

    auto uh = X( 0_c );
    auto alphah = X( 1_c );

    if ( boption( "print-matlab-fields" ) )
    {
        auto matlabContext = Th->context();
        node_type center( 3 );
        center( 0 ) = 0.5;
        center( 1 ) = 0.5;
        center( 2 ) = 0.0;
        matlabContext.add( center );

        epsilon.printMatlab( "epsilon_matlab",
                             matlabContext,
                             matlabEpsilonFormat,
                             true,
                             "epsilon" );
        sigma.printMatlab( "sigma_matlab",
                           matlabContext,
                           matlabSigmaFormat,
                           true,
                           "sigma" );
    }

    std::cout << "qs_sb9 minimal formulation assembled"
              << " with E=" << E
              << ", nu=" << nu
              << ", load=" << doption( "load" )
              << ", strategy=" << ( useStaticCondensation ? "static_condensation" : "monolithic" ) << "\n";

    auto e = exporter( _mesh=mesh );
    e->addRegions();
    e->add( "u", uh );
    e->add( "alpha", alphah );
    e->save();

    return 0;
}
