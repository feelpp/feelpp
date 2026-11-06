/*
 * SPDX-FileContributor: Hanna Chetouane <h.chetouane@unistra.fr>
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
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
#include <feel/feelvf/sb9.hpp>
#include <feel/feelvf/sb9_quadrature.hpp>
#include <feel/feelvf/sb9_shear.hpp>
#include <feel/feelvf/vf.hpp>

#include <boost/format.hpp>

using namespace Feel;

namespace
{
using mesh_type = Mesh<Hypercube<3,Dynamic>>;

po::options_description
makeOptions()
{
    po::options_description options( "qs_sb9_minimal options" );
    options.add_options()
        ( "E", po::value<double>()->default_value( 1.0 ), "Young modulus" )
        ( "nu", po::value<double>()->default_value( 0.3 ), "Poisson ratio" )
        ( "load", po::value<double>()->default_value( 1.0 ), "constant x traction on XPlus" )
        ( "monolithic", po::value<bool>()->default_value( true ), "use the monolithic mixed solve instead of SB9 static condensation" )
        ( "print-matlab-fields", po::value<bool>()->default_value( false ),
          "print epsilon/sigma tensor fields in MATLAB validation order" )
        ( "no-solve", po::value<bool>()->default_value( false ), "assemble only" );
    return options.add( feel_options() );
}

AboutData
makeAbout()
{
    AboutData about( "qs_sb9_minimal",
                     "qs_sb9_minimal",
                     "0.1",
                     "Minimal SB9 mixed shell formulation with stabilization terms",
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

Point(1) = {0, 0, 0, 1};
Point(2) = {10, 0, 0, 1};
Point(3) = {0, 10, 0, 1};
Point(4) = {10, 10, 0, 1};
Point(5) = {0, 0, 1, 1};
Point(6) = {10, 0, 1, 1};
Point(7) = {0, 10, 1, 1};
Point(8) = {10, 10, 1, 1};

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
Physical Point("Point1") = {1};
Physical Point("Point2") = {2};
Physical Point("Point3") = {3};
Physical Point("Point4") = {4};
Physical Point("Point5") = {5};
Physical Point("Point6") = {6};
Physical Point("Point7") = {7};
Physical Point("Point8") = {8};
)";
}

std::shared_ptr<mesh_type>
createUnitShellMesh()
{
    Environment::changeRepository( _directory=boost::format( "quickstart/%1%/" ) % Environment::about().appName() );
    return createGMSHMesh( _mesh=new mesh_type( RuntimeOrder{1} ),
                           _desc=geo( _filename="qs_sb9_minimal_unit_shell.geo",
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

    auto Uh = Pchv<Dynamic>( mesh, RuntimeOrder{1} );
    auto Ah = Pdh<Dynamic>( mesh, RuntimeOrder{0} );
    auto Xh = productPtr( Uh, Ah );
    auto Th = Pdhms<Dynamic>( mesh, RuntimeOrder{1} );

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

    bool const useStaticCondensation = !boption( "monolithic" );
    solve::strategy const strategy = useStaticCondensation
                                         ? solve::strategy::static_condensation
                                         : solve::strategy::monolithic;
    auto a = blockform2( *Xh, strategy, backend() );
    auto l = blockform1( *Xh, strategy, backend() );

    auto zt = zeta();
    auto C = isotropic_stiffness<3>( lambda, mu );

    auto epsU = sb9MembraneBending( u, zt );
    auto epsV = sb9MembraneBending( v, zt );
    auto epsP = sb9Pinching( u, zt );
    auto epsQ = sb9Pinching( v, zt );
    auto epsS = sb9Shearing( u, zt );
    auto epsT = sb9Shearing( v, zt );
    auto epsW = sb9PinchingW9( alpha );
    auto epsZ = sb9PinchingW9( beta );

    double coef = 1/std::sqrt(2);
    auto epsShellTrial = vec( component<0, 0>( epsU ),           // E11   
                              coef * component<1, 0>( epsU ),    // E12
                              coef * component<2, 0>( epsS ),    // E13
                              component<3, 0>( epsU ),           // E22
                              coef * component<4, 0>( epsS ),    // E23
                              component<5, 0>( epsP ) );         // E33
    auto epsShellTest = vec( component<0, 0>( epsV ),            // E11
                             coef * component<1, 0>( epsV ),     // E12
                             coef * component<2, 0>( epsT ),     // E13
                             component<3, 0>( epsV ),            // E22
                             coef * component<4, 0>( epsT ),     // E23
                             component<5, 0>( epsQ ) );          // E33


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


    a( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                _quad=sb9Quad,
                                _expr=ddot( C, epsShellTrial, epsShellTest ) );

    // Compose the former sb9Stabilization helper with the current compact
    // Bs1/Bs2 (one row), Bs3 (two rows), and Bs4 (three rows) operators.
    // Retain the original minimal example's element-volume weighting.
    auto volume = shellThickness() * shellArea0();
    auto normalStiffness = cst( lambda + 2.0 * mu ) * volume;
    auto ds1 = normalStiffness * shellInvJ0_00() * shellInvJ0_00() / cst( 3.0 );
    auto ds2 = normalStiffness * ( shellInvJ0_10() * shellInvJ0_10() +
                                  shellInvJ0_11() * shellInvJ0_11() ) / cst( 3.0 );
    auto ds3 = normalStiffness * shellInvJ0_22() * shellInvJ0_22() / cst( 3.0 );
    auto bs1u = sb9Bs1( u );
    auto bs1v = sb9Bs1( v );
    auto bs2u = sb9Bs2( u );
    auto bs2v = sb9Bs2( v );
    auto bs3u = sb9Bs3( u );
    auto bs3v = sb9Bs3( v );
    auto bs4u = sb9Bs4( u );
    auto bs4v = sb9Bs4( v );
    auto modeStabilization =
        ds1 * ( component<0, 0>( bs3u ) * component<0, 0>( bs3v ) +
                component<0, 0>( bs4u ) * component<0, 0>( bs4v ) / cst( 3.0 ) ) +
        ds2 * ( component<1, 0>( bs3u ) * component<1, 0>( bs3v ) +
                component<1, 0>( bs4u ) * component<1, 0>( bs4v ) / cst( 3.0 ) ) +
        ds3 * ( component<0, 0>( bs1u ) * component<0, 0>( bs1v ) +
                component<0, 0>( bs2u ) * component<0, 0>( bs2v ) +
                component<2, 0>( bs4u ) * component<2, 0>( bs4v ) / cst( 3.0 ) );
    auto shearStabilization = cst( 5.0 * mu / 18.0 ) * volume *
                             ( inner( sb9Bc1( u ), sb9Bc1( v ) ) +
                               inner( sb9Bc2( u ), sb9Bc2( v ) ) );
    a( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                _quad=sb9Quad,
                                _expr=modeStabilization + shearStabilization );

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

    l( 0_c ) += integrate( _range=markedfaces( mesh, "XPlus" ),
                           _expr=inner( vec( cst( 1.0 ), cst( 0.0 ), cst( 0.0 ) ), v ) );
    
    a.matrix().printMatlab( "KSB9.m" );
    // l.vector().printMatlab( "FSB9.m" );

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
        center( 2 ) = 0.5;
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

    std::cout << "qs_sb9_minimal formulation assembled"
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
