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
#include <feel/feelts/newmark.hpp>
#include <feel/feelfit/fit.hpp>

#include <boost/format.hpp>

#include <iostream>
#include <set>
#include <stdexcept>
#include <string>

using namespace Feel;

namespace
{

using sb9_mesh_type = Mesh<Hypercube<3>>;

po::options_description
makeOptions()
{
    po::options_description options( "sb9 options" );
    options.add_options()
        ( "shell-membrane-stab", po::value<double>()->default_value( 1.0 ),
          "membrane stablization coefficient applied to the additional SB9 Bs1..Bs4 stabilization block" )
        ( "shell-bending-stab", po::value<double>()->default_value( 1.0 ),
          "bending stablization coefficient applied to the additional SB9 Bs1..Bs4 stabilization block" )
        ( "shell-pinching-stab", po::value<double>()->default_value( 1.0 ),
          "pinching applied to the additional SB9 Bs1..Bs4 stabilization block" )
        ( "monolithic", po::value<bool>()->default_value( false ),
          "use the monolithic mixed solve instead of the default SB9 static condensation" )
        ( "E", po::value<double>()->default_value( 1.0 ),
          "Young coefficient" )
        ( "nu", po::value<double>()->default_value( 0.0 ),
          "Poisson ratio" )
        ( "sb9-ts.dynamic", po::value<bool>()->default_value( false ),
          "dynamic mode" )
        ( "sb9-ts.rho", po::value<double>()->default_value( 0.0 ),
          "density of the solid" )
        ( "sb9-ts.final-time", po::value<double>()->default_value( 1.0 ),
          "final time of the simulation" )
        ( "sb9-ts.init-time", po::value<double>()->default_value( 0.0 ),
          "init time of the simulation" )
        ( "sb9-ts.time-step", po::value<double>()->default_value( 0.1 ),
          "time step of the simulation" )
        ( "csv.filename", po::value<std::string>()->default_value( "" ),
          "csv file" )
        ;
    return options.add( feel_options() );
}

AboutData
makeAbout()
{
    AboutData about( "sb9",
                     "sb9",
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

std::shared_ptr<sb9_mesh_type>
createUnitShellMesh()
{
  Environment::changeRepository( _directory=boost::format( "quickstart/%1%/" ) % Environment::about().appName() );
  return createGMSHMesh( _mesh=new sb9_mesh_type,
                           _desc=geo( _filename="sb9_unit_shell.geo",
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

    try
    {
        Environment env( _argc=argc, _argv=argv, _desc=makeOptions(), _about=makeAbout() );

        std::shared_ptr<sb9_mesh_type> mesh;
        if ( Environment::vm().count( "gmsh.filename" ) )
        {
          auto const mesh_file_opt = soption( _name = "gmsh.filename" );
          if ( !mesh_file_opt.empty() )
          {
            auto mesh_file = Environment::expand( mesh_file_opt );
            mesh = loadMesh( _mesh = new sb9_mesh_type(), _filename = mesh_file );
          }
        }
        if ( !mesh )
        {
          std::cout << "No gmsh filename provided, use createUnitShellMesh()" << std::endl;
          mesh = createUnitShellMesh();
        }

        auto const csv_file = Environment::expand( soption( _name = "csv.filename" ) );

        double const E = doption( "E" );
        double const nu = doption( "nu" );
        std::cout << "E = " << E << "; nu = " << nu << std::endl;

        double const lambda = E * nu / ( ( 1.0 + nu ) * ( 1.0 - 2.0 * nu ) );
        double const mu = E / ( 2.0 * ( 1.0 + nu ) );

        double const shellMembraneStab = doption( "shell-membrane-stab" );
        double const shellBendingStab = doption( "shell-bending-stab" );
        double const shellPinchingStab = doption( "shell-pinching-stab" );
        bool const useStaticCondensation = !boption( "monolithic" );

        std::cout << "sb9: E=" << E
                  << ", nu=" << nu
                  << ", strategy=" << ( useStaticCondensation ? "static_condensation" : "monolithic" ) << "\n";


        bool const dynamic = boption( "sb9-ts.dynamic" );
        double const rho = doption( "sb9-ts.rho" );    
        double const t0 = doption( "sb9-ts.init-time" ); 
        double const tf = doption( "sb9-ts.final-time" );
        double const dt = doption( "sb9-ts.time-step" );
        
        if( dynamic ) {
            std::cout << "dynamic : " << dynamic
                      << ", rho = " << rho
                      << ", init time = " << t0
                      << ", final time = " << tf
                      << ", time step = " << dt << "\n";
        }

        auto Uh = Pchv<1>( mesh );
        auto Ah = Pdh<0>( mesh );
        auto Xh = productPtr( Uh, Ah );
        auto algebraBackend = backend( _rebuild=true, _worldcomm=Uh->worldCommPtr() );

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
        auto shellShearWeight = cst( 1.25 ) * ( cst( 1.0 ) - zt * zt );  

        auto membraneTrial = sb9MembraneBending( u, zt );
        auto membraneTest = sb9MembraneBending( v, zt );
        auto pinchingTrial = sb9Pinching( u, zt, cst( 1.0 ) );
        auto pinchingTest = sb9Pinching( v, zt, cst( 1.0 ) );
        auto shearTrial = sb9Shear( u, shellShearWeight );
        auto shearTest = sb9Shear( v, shellShearWeight );
        auto alphaTrial = sb9W9( alpha, cst( 1.0 ) );
        auto alphaTest = sb9W9( beta, cst( 1.0 ) );

        double coef = 1/std::sqrt(2);
        auto epsShellTrialMandel = vec( component<0, 0>( membraneTrial ),
                                        cst( coef ) * component<1, 0>( membraneTrial ),
                                        cst( coef ) * component<2, 0>( shearTrial ),
                                        component<3, 0>( membraneTrial ),
                                        cst( coef ) * component<4, 0>( shearTrial ),
                                        component<5, 0>( pinchingTrial ) );
        auto epsShellTestMandel = vec( component<0, 0>( membraneTest ),
                                       cst( coef ) * component<1, 0>( membraneTest ),
                                       cst( coef ) * component<2, 0>( shearTest ),
                                       component<3, 0>( membraneTest ),
                                       cst( coef ) * component<4, 0>( shearTest ),
                                       component<5, 0>( pinchingTest ) );


        auto sb9Quad = sb9ThroughThicknessLobatto5();

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

        auto sb9Stabilization = [&]( auto const& trialField, auto const& testField )
        {
            auto invJ0 = inv( shellJacobian0() );
            auto j11 = component<0, 0>( invJ0 );
            auto j21 = component<1, 0>( invJ0 );
            auto j22 = component<1, 1>( invJ0 );
            auto j33 = component<2, 2>( invJ0 );

            auto normalStiffness = cst( ( lambda + 2.0 * mu ) );
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

            auto shearPart = cst( 2.0 * 5.0/18.0 * mu ) * 
                             ( inner( sb9Bc1( trialField ), sb9Bc1( testField ) ) +
                               inner( sb9Bc2( trialField ), sb9Bc2( testField ) ) );

            auto bsPart = cst( 1.0 / 3.0 ) * shellPinchingStab * dz * ( c0( bs1u, bs1v ) + c0( bs2u, bs2v ) ) +
                          cst( 1.0 / 3.0 ) * shellMembraneStab * ( dx * c0( bs3u, bs3v ) + dy * c1( bs3u, bs3v ) ) +
                          cst( 1.0 / 9.0 ) * (  
                                               ( shellBendingStab * dx * c0( bs4u, bs4v ) + shellBendingStab * dy * c1( bs4u, bs4v ) ) +
                                               shellPinchingStab * dz * c2( bs4u, bs4v ) );

            return shearPart + bsPart;
        };

        a( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                    _quad=sb9Quad,
                                    _expr=sb9Stabilization( u, v ) );


        auto e = exporter( _mesh=mesh );
        e->addRegions();

        if ( !dynamic )
        {
          auto lDisplacement = l( 0_c );
          lDisplacement += integrate( _range=markedfaces( mesh, "XPlus" ),
                                          _expr=inner( vec( cst(1.0), cst(0.0), cst(0.0) ), v ) );
          l.close();
          a.close();

          a.row( 0_c ) += on( _range=markedfaces( mesh, "XMoins" ),
                                  _rhs=l( 0_c ),
                                  _element=uBc,
                                  _expr=vec( cst(0.0), cst(0.0), cst(0.0) ),
                                  _type="elimination" );

          a.solve( _rhs=l, _solution=X,
                        _condense=useStaticCondensation,
                        _condenser=condenser_sb9() );

          auto uh = X(0_c);
          e->add( "u", uh );
        }

        else 
        {
          auto uh = X(0_c);
          uh.zero();

          auto newmark_scheme = newmark(_space = Uh, _name = "sb9", _initial_time = t0, _final_time = tf,
                          _time_step = dt, _rank_proc_in_files_name = true, _restart_at_last_save = true);
          newmark_scheme->start();
          newmark_scheme->initialize( uh ); 

          a( 0_c, 0_c ) += integrate( _range=elements( mesh ), _expr= rho * newmark_scheme->polyDerivCoefficient() * inner( id(u), id(v) ) );

          auto time = newmark_scheme->time();
          while ( !newmark_scheme->isFinished() )
          {
            time = newmark_scheme->time();
            l( 0_c ).zero();
            l( 0_c ) += integrate( _range = elements( mesh ), _expr = rho * inner( idv( newmark_scheme->polyDeriv() ), id(v) ) );
            
            // plate (plaque_percee or plaque_pas_percee)
            l( 0_c ) += integrate( _range = markedpoints( mesh, "force_apply" ),   
                                   _expr = inner( vec( cst(0.0), cst(0.0), - fit( expr( time ), csv_file, "Temps", "Marteau", "P1") ), id(v) ) );
            // tank (cuve)
            // l( 0_c ) += integrate( _range = markedpoints( mesh, "pt_P1" ),   
            //                        _expr = inner( vec( cst(0.0), cst(0.0), - fit( evaluation_time, csv_file, "Temps", "Marteau", "P1") ), id(v) ) );

            l.close();
            a.close();

            a.solve( _rhs=l, _solution=X,
                          _condense=useStaticCondensation,
                          _condenser=condenser_sb9() );

            uh = X(0_c);
            newmark_scheme->next( uh );

            e->step(time)->add( "u", uh );
            e->step(time)->add("acceleration", newmark_scheme->currentAcceleration() );
          }
        }

        e->save();
       
        auto uh = X( 0_c );
        auto alphah = X( 1_c );
        auto dispNormMax = normLinf( _range=elements( mesh ), _pset=_Q<2>(), _expr=norm2( idv( uh ) ) );
        std::cout << "max displacement norm = " << dispNormMax.value()
                  << " at " << dispNormMax.arg().transpose() << "\n";
    }
    catch ( ... )
    {
        handleExceptions();
    }

    return 0;
}
