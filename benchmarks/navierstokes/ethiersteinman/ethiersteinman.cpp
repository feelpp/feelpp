/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2007-04-03

  Copyright (C) 2007 EPFL

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file ethiersteinman.cpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2007-04-03
 */

#include <sstream>

#include "ethiersteinman.hpp"

#include <feel/feelcore/pslogger.hpp>

#include <feel/feelfilters/importergmsh.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>

namespace Feel
{

const uint16_type EthierSteinman::uOrder;
const uint16_type EthierSteinman::pOrder;
const int EthierSteinman::Dim;
const uint16_type EthierSteinman::imOrder;

EthierSteinman::EthierSteinman( int argc, char** argv, AboutData const& ad )
    :
    super( argc, argv, ad ),
    M_meshSize( this->vm()["hsize"].as<double>() ),
    M_bcCoeffDiff( this->vm()["bccoeffdiff"].as<double>() ),
    M_bcCoeffConv( this->vm()["bccoeffconv"].as<double>() ),
    M_mu ( this->vm()["mu"].as<double>() ),
    M_exporter( new ExporterEnsight<mesh_type>( "ethiersteinman" ) ),
    M_timeSet( new timeset_type( "ethiersteinman" ) ),
    M_timers(),
    M_im(),
    M_uErrorL2( -1.0 ),
    M_uErrorH1( -1.0 ),
    M_pErrorL2( -1.0 ),
    M_divError( -1.0 )
{
    std::stringstream cmdline;

    for ( int i=0; i<argc; ++i )
    {
        cmdline << argv[i] << " ";
    }

    VLOG(1) << cmdline.str() << "\n";

    VLOG(1) << "[EthierSteinman] hsize   = " << M_meshSize << "\n";
    VLOG(1) << "[EthierSteinman] bcCdiff = " << M_bcCoeffDiff << "\n";
    VLOG(1) << "[EthierSteinman] bcCconv = " << M_bcCoeffConv << "\n";
    VLOG(1) << "[EthierSteinman] mu      = " << M_mu << "\n";
    VLOG(1) << "[EthierSteinman] export  = "
            << this->vm()["export"].as<int>() << "\n";

    M_timeSet->setTimeIncrement( this->vm()["dt"].as<double>() );
    M_exporter->addTimeSet( M_timeSet );
}

EthierSteinman::EthierSteinman( int argc,
                                char** argv,
                                AboutData const& ad,
                                po::options_description const& od )
    :
    super( argc, argv, ad, od ),
    M_meshSize( this->vm()["hsize"].as<double>() ),
    M_bcCoeffDiff( this->vm()["bccoeffdiff"].as<double>() ),
    M_bcCoeffConv( this->vm()["bccoeffconv"].as<double>() ),
    M_mu ( this->vm()["mu"].as<double>() ),
    M_exporter( new ExporterEnsight<mesh_type>( "ethiersteinman" ) ),
    M_timeSet( new timeset_type( "ethiersteinman" ) ),
    M_timers(),
    M_im(),
    M_uErrorL2( -1.0 ),
    M_uErrorH1( -1.0 ),
    M_pErrorL2( -1.0 ),
    M_divError( -1.0 )
{
    std::stringstream cmdline;

    for ( int i=0; i<argc; ++i )
    {
        cmdline << argv[i] << " ";
    }

    VLOG(1) << cmdline.str() << "\n";

    VLOG(1) << "[EthierSteinman] hsize   = " << M_meshSize << "\n";
    VLOG(1) << "[EthierSteinman] bcCdiff = " << M_bcCoeffDiff << "\n";
    VLOG(1) << "[EthierSteinman] bcCconv = " << M_bcCoeffConv << "\n";
    VLOG(1) << "[EthierSteinman] mu      = " << M_mu << "\n";
    VLOG(1) << "[EthierSteinman] export  = "
            << this->vm()["export"].as<int>() << "\n";

    M_timeSet->setTimeIncrement( this->vm()["dt"].as<double>() );
    M_exporter->addTimeSet( M_timeSet );
}

EthierSteinman::EthierSteinman( EthierSteinman const& tc )
    :
    super( tc ),
    M_meshSize( tc.M_meshSize ),
    M_bcCoeffDiff( tc.M_bcCoeffDiff ),
    M_bcCoeffConv( tc.M_bcCoeffConv ),
    M_mu( tc.M_mu ),
    M_exporter( new ExporterEnsight<mesh_type>( "ethiersteinman" ) ),
    M_timeSet( new timeset_type( "ethiersteinman" ) ),
    M_timers( tc.M_timers ),
    M_im(),
    M_uErrorL2( tc.M_uErrorL2 ),
    M_uErrorH1( tc.M_uErrorH1 ),
    M_pErrorL2( tc.M_pErrorL2 ),
    M_divError( tc.M_divError )
{
    VLOG(1) << "[EthierSteinman] hsize   = " << M_meshSize << "\n";
    VLOG(1) << "[EthierSteinman] bcCdiff = " << M_bcCoeffDiff << "\n";
    VLOG(1) << "[EthierSteinman] bcCconv = " << M_bcCoeffConv << "\n";
    VLOG(1) << "[EthierSteinman] mu      = " << M_mu << "\n";
    VLOG(1) << "[EthierSteinman] export  = "
            << this->vm()["export"].as<int>() << "\n";

    M_timeSet->setTimeIncrement( this->vm()["dt"].as<double>() );
    M_exporter->addTimeSet( M_timeSet );
}

EthierSteinman::mesh_ptr_type
EthierSteinman::createMesh( double meshSize, double R )
{
    M_timers["mesh"].first.restart();
    mesh_ptr_type mesh( new mesh_type );

    if ( R == 0 )
    {
        GmshHypercubeDomain<Dim,1,ENTITY> td;
        td.setCharacteristicLength( meshSize );
        td.setX( std::make_pair( -1.0, 1.0 ) );
        td.setY( std::make_pair( -1.0, 1.0 ) );
        //         td.setZ( std::make_pair( -1.0, 1.0 ) );
        ImporterGmsh<mesh_type>
        import( td.generate( ENTITY<Dim,1,Dim>::name().c_str() ) );
        mesh->accept( import );
    }

    else
    {
        Gmsh gmsh;
        gmsh.setOrder( GMSH_ORDER_ONE );

        std::ostringstream mesh_desc;
        std::ostringstream mesh_name;
        mesh_desc << "h=" << meshSize << ";\n"
                  << "R=" << R << ";\n"
                  << "Point(1) = { R, 0,0,h};\n"
                  << "Point(2) = { 0, R,0,h};\n"
                  << "Point(3) = {-R, 0,0,h};\n"
                  << "Point(4) = { 0,-R,0,h};\n"
                  << "Point(5) = { 0, 0,0,h};\n"
                  << "Circle(1) = {1,5,2};\n"
                  << "Circle(2) = {2,5,3};\n"
                  << "Circle(3) = {3,5,4};\n"
                  << "Circle(4) = {4,5,1};\n"
                  << "Line Loop(5) = {1,2,3,4};\n"
                  << "Plane Surface(6) = {5};\n"
                  << "Extrude Surface {6, {0,0,1}};\n"
                  << "Physical Surface(29) = {15,19,23,27};\n"
                  << "Physical Surface(30) = {28};\n"
                  << "Physical Surface(31) = {6};\n"
                  << "Surface Loop(32) = {6,-15,-19,-23,-27,-28};\n"
                  << "Volume(1) = {32};\n"
                  << "Physical Volume(2) = {1};\n";
        mesh_name << "poiseuille." << meshSize;

        std::string fname = gmsh.generate( mesh_name.str(), mesh_desc.str() );
        ImporterGmsh<mesh_type> import( fname );
        mesh->accept( import );
    }

    M_timers["mesh"].second = M_timers["mesh"].first.elapsed();
    VLOG(1) << "[timer] createMesh(): " << M_timers["mesh"].second << "\n";
    return mesh;
} // EthierSteinman::createMesh

void
EthierSteinman::run()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    double dt = this->vm()["dt"].as<double>();

    this->changeRepository( boost::format( "%1%/Re_%2%/P%3%P%4%/h_%5%/dt_%6%" )
                            % this->about().appName()
                            % ( 1.0/M_mu )
                            % uOrder
                            % pOrder
                            % M_meshSize
                            % dt );
    this->setLogs();

    PsLogger psLogger( "ps.log" );
    psLogger.log( "t=0, start" );

    using namespace Feel::vf;

    /*
     * First we create the mesh
     */
    std::cout << "[EthierSteinman] mesh\n" << std::flush;
    //     value_type R = 0.5;
    mesh_ptr_type mesh = createMesh( M_meshSize );
    psLogger.log( "t=0, meshed" );

    /*
     * The function spaces and some associate elements are then defined
     */
    M_timers["init"].first.restart();
    std::cout << "[EthierSteinman] spaces\n" << std::flush;
    space_U_ptrtype space_U = space_U_type::New( mesh );
    space_u_ptrtype space_u = space_u_type::New( mesh );
    space_p_ptrtype space_p = space_p_type::New( mesh );
    space_i_ptrtype space_i = space_i_type::New( mesh );
    //space_u->dof()->showMe();
    psLogger.log( "t=0, spaces" );

    VLOG(1) << "[EthierSteinman] velocity dofs per component "
            << space_u->nbDof() << "\n";
    VLOG(1) << "[EthierSteinman] velocity dofs total         "
            << space_U->nbDof() << "\n";
    VLOG(1) << "[EthierSteinman] pressure dofs               "
            << space_p->nbDof() << "\n";
    VLOG(1) << "[EthierSteinman] total    dofs               "
            << space_U->nbDof() + space_p->nbDof() << "\n";

    std::cout << "[EthierSteinman] elements\n" << std::flush;
    element_U_type U    ( space_U, "U" );
    element_u_type ux   ( space_u, "ux" );
    element_u_type uy   ( space_u, "uy" );
    element_u_type uz   ( space_u, "uz" );
    element_u_type uxn  ( space_u, "uxn" );
    element_u_type uyn  ( space_u, "uyn" );
    element_u_type uzn  ( space_u, "uzn" );
    element_u_type uxo  ( space_u, "uxo" );
    element_u_type uyo  ( space_u, "uyo" );
    element_u_type uzo  ( space_u, "uzo" );
    element_u_type ux0  ( space_u, "ux0" );
    element_u_type uy0  ( space_u, "uy0" );
    element_u_type uz0  ( space_u, "uz0" );
    element_p_type p    ( space_p, "p" );
    element_p_type pn   ( space_p, "pn" );
    element_p_type po   ( space_p, "po" );
    element_p_type p0   ( space_p, "p0" );
    element_i_type indic( space_i, "indic" );
    psLogger.log( "t=0, elements" );

    // -- initial condition
    std::cout << "[EthierSteinman] initial conditions\n" << std::flush;

    // ethier steinman solution
    //     value_type pi = 4.0 * math::atan( value_type( 1.0 ) );
    //     value_type a = pi/2;
    //     value_type d = pi/4;
    //     value_type timeFactor = std::exp(-d*d*M_mu*dt);
    //     ux0 = vf::project( space_u, elements(*mesh),
    //                        -a *
    //                        ( exp(a*Px()) * sin(a*Py()+d*Pz()) +
    //                          exp(a*Pz()) * cos(a*Px()+d*Py()) ) );
    //     uy0 = vf::project( space_u, elements(*mesh),
    //                        -a *
    //                        ( exp(a*Py()) * sin(a*Pz()+d*Px()) +
    //                          exp(a*Px()) * cos(a*Py()+d*Pz()) ) );
    //     uz0 = vf::project( space_u, elements(*mesh),
    //                        -a *
    //                        ( exp(a*Pz()) * sin(a*Px()+d*Py()) +
    //                          exp(a*Py()) * cos(a*Pz()+d*Px()) ) );
    //     p0 = vf::project( space_p, elements(*mesh),
    //                       -a*a / 2 *
    //                       ( exp(2*a*Px()) + exp(2*a*Py()) + exp(2*a*Pz()) +
    //                         2 * sin(a*Px()+d*Py()) * cos(a*Pz()+d*Px()) * exp(a*(Py()+Pz())) +
    //                         2 * sin(a*Py()+d*Pz()) * cos(a*Px()+d*Py()) * exp(a*(Pz()+Px())) +
    //                         2 * sin(a*Pz()+d*Px()) * cos(a*Py()+d*Pz()) * exp(a*(Px()+Py())) ) );

    // zero solution
    //     value_type timeFactor = 1.0;
    //     ux0 = vf::project( space_u, elements(*mesh), constant(0.0) );
    //     uy0 = vf::project( space_u, elements(*mesh), constant(0.0) );
    //     uz0 = vf::project( space_u, elements(*mesh), constant(0.0) );
    //     p0  = vf::project( space_p, elements(*mesh), constant(0.0) );

    // constant 1-2-3 solution
    //     value_type timeFactor = 1.0;
    //     ux0 = vf::project( space_u, elements(*mesh), constant(1.0) );
    //     uy0 = vf::project( space_u, elements(*mesh), constant(2.0) );
    //     uz0 = vf::project( space_u, elements(*mesh), constant(3.0) );
    //     p0  = vf::project( space_p, elements(*mesh), constant(0.0) );

    // Poiseuille flow
    //     value_type timeFactor = 1.0;
    //     ux0 = vf::project( space_u, elements(*mesh), constant(0.0) );
    //     uy0 = vf::project( space_u, elements(*mesh), constant(0.0) );
    //     uz0 = vf::project( space_u, elements(*mesh),
    //                        1.0 - ( Px()*Px() + Py()*Py() ) / (R*R) );
    //     p0  = vf::project( space_p, elements(*mesh), -4.0/(R*R)*Pz() );

    // Kim-Moin solution
    value_type pi = 4.0 * math::atan( value_type( 1.0 ) );
    value_type timeFactor = std::exp( -2. * pi * pi * M_mu * dt );
    ux0 = vf::project( space_u, elements( *mesh ), -cos( pi*Px() )*sin( pi*Py() ) );
    uy0 = vf::project( space_u, elements( *mesh ),  sin( pi*Px() )*cos( pi*Py() ) );
    p0  = vf::project( space_p, elements( *mesh ),
                       -0.25*( cos( 2*pi*Px() )+cos( 2*pi*Py() ) ) );

    ux = ux0;
    uy = uy0;
    uz = uz0;
    p  = p0;
    uxn = ux;
    uyn = uy;
    uzn = uz;
    pn  = p;
    U.comp( X ) = ux;
    U.comp( Y ) = uy;
    //     U.comp(Z) = uz;

    indic = vf::project( space_i, elements( *mesh ), constant( 0.0 ) );
    mesh->updateMarker2( indic );

    std::cout << "[EthierSteinman] backends\n" << std::flush;
    backend_ptrtype backendNS( new backend_type );
    backendNS->set_noisy( this->vm()["verbose"].as<int>() );
    backendNS->set_maxiter( this->vm()["maxiter"].as<int>() );
    backendNS->set_fillin( this->vm()["fillin"].as<int>() );
    backendNS->set_threshold( this->vm()["threshold"].as<double>() );
    backendNS->set_tol( this->vm()["tolerance"].as<double>() );
    backendNS->set_solver_type( this->vm()["solver"].as<std::string>() );
    backendNS->set_preconditioner_type
    ( this->vm()["precond"].as<std::string>() );
    backendNS->set_restart( this->vm()["restart"].as<int>() );
    backendNS->set_direct( this->vm().count( "direct" ) );

    backend_ptrtype backendSymm( new backend_type );
    backendSymm->set_noisy( this->vm()["verbose"].as<int>() );
    backendSymm->set_maxiter( this->vm()["maxiter"].as<int>() );
    backendSymm->set_fillin( this->vm()["fillin"].as<int>() );
    backendSymm->set_threshold( this->vm()["threshold"].as<double>() );
    backendSymm->set_tol( this->vm()["tolerance"].as<double>() );
    backendSymm->set_preconditioner_type( "id" );
    backendSymm->set_symmetric( true );
    backendSymm->set_direct( this->vm().count( "direct" ) );
    psLogger.log( "t=0, backends" );

    std::set<flag_type> dirichletFlags;
    std::set<flag_type> neumannFlags;

    //     dirichletFlags.insert(6);
    //     dirichletFlags.insert(15);
    //     dirichletFlags.insert(19);
    //     dirichletFlags.insert(23);
    //     dirichletFlags.insert(27);
    //     dirichletFlags.insert(28);
    //     neumannFlags.insert(28);

    //     dirichletFlags.insert(29);
    //     dirichletFlags.insert(31);
    //     dirichletFlags.insert(30);

    dirichletFlags.insert( 1 );
    dirichletFlags.insert( 2 );
    dirichletFlags.insert( 3 );
    dirichletFlags.insert( 4 );

    std::cout << "[EthierSteinman] oseen creation\n" << std::flush;
    Oseen<space_u_type, space_p_type, imOrder, ENTITY>
    oseen( space_u, space_p, backendNS, dirichletFlags, neumannFlags );
    oseen.set_bccoeffdiff( M_bcCoeffDiff );
    oseen.set_bccoeffconv( M_bcCoeffConv );
    oseen.set_stabcoeffdiv( this->vm()["stabcoeff-div"].as<double>() );
    oseen.set_stabcoeffp( this->vm()["stabcoeff-p"].as<double>() );
    double stabTheta = this->vm()["stabtheta"].as<double>();

    if ( stabTheta > 0.0 )
    {
        oseen.decouplePstab( pn, stabTheta );
    }

    oseen.set_epscompress( this->vm()["epscompress"].as<double>() );
    oseen.set_divdivcoeff( this->vm()["divdivcoeff"].as<double>() );
    psLogger.log( "t=0, oseen" );

    // mass matrix on velocity space, for L2 norms
    std::cout << "[EthierSteinman] mass matrices\n" << std::flush;
    OperatorLinear<space_u_type, space_u_type, backend_type>
    massU( space_u, space_u, backendSymm );
    massU = integrate( elements( *mesh ), M_im, id( ux )*idt( ux ) );

    // mass matrix on pressure space, for L2 norms
    OperatorLinear<space_p_type, space_p_type, backend_type>
    massP( space_p, space_p, backendSymm );
    massP = integrate( elements( *mesh ), M_im, id( p )*idt( p ) );

    // laplace matrix on velocity space, for H1 norms
    OperatorLinear<space_u_type, space_u_type, backend_type>
    laplaceU( space_u, space_u, backendSymm );
    laplaceU = integrate( elements( *mesh ), M_im,
                          dx( ux )*dxt( ux ) + dy( ux )*dyt( ux ) );

    psLogger.log( "t=0, mass matrices" );

    double time = 0;
    std::cout << "[EthierSteinman] export\n" << std::flush;
    this->exportResults( 0, time, U, ux, uy, uz, p,
                         massU, laplaceU, massP,
                         ux0, uy0, uz0, p0 );

    time               = dt;
    double fixpointTol = this->vm()["fixpointtol"].as<double>();
    int    maxSubIter  = this->vm()["maxsubiter" ].as<int>();

    value_type eps = type_traits<value_type>::epsilon();
    value_type epsCompress = this->vm()["epscompress"].as<double>();

    M_timers["init"].second = M_timers["init"].first.elapsed();

    // basic initialization of oseen operator
    std::cout << "[EthierSteinman] update oseen\n" << std::flush;
    //     oseen.update( /* itRan = */ elements(*mesh),
    //                   /* sigma = */ 1.0/dt,
    //                   /* nuInc = */ M_mu,
    //                   /* nuAbs = */ 0.0,
    //                   /* beta  = */ 0.0*oneX(),
    //                   /* f     = */ 0.0*oneX(),
    //                   /* c     = */ 0.0,
    //                   /* g     = */ 0.0*oneX(),
    //                   /* noSlip= */ 0.0,
    //                   /* updtJ = */ false );
    oseenUpdateInit( oseen, mesh, dt );
    psLogger.log( "t=0, oseen update" );
    std::stringstream msg;

    // --- Time loop
    for ( int iter = 0;
            time-dt/2 < this->vm()["ft"].as<double>();
            ++iter, time += dt )
    {
        VLOG(1) << "[EthierSteinman] t = " << time << "\n";
        std::cout << "[EthierSteinman] t = " << time << "\n";

        ux0 *= timeFactor;
        uy0 *= timeFactor;
        uz0 *= timeFactor;
        p0 *= timeFactor*timeFactor;

        double fixpointErr = 2*fixpointTol+1.0;
        uint32_type subiter;

        for ( subiter = 0;
                ( fixpointErr>fixpointTol ) && ( subiter < maxSubIter );
                ++subiter )
        {
            std::cout << "[EthierSteinman] update oseen"
                      << ( backendNS->reusePC() ? "" : " (rebuild ip)" )
                      << "\n" << std::flush;

            M_timers["updateNS"].first.restart();

            if ( this->vm().count( "bdf1" ) || ( iter == 0 ) )
            {
                // BDF1
                oseenUpdateBdf1( oseen, dt, uxn, uyn, uzn, ux, uy, uz,
                                 ux0, uy0, uz0,
                                 !backendNS->reusePC(),
                                 pn );
            }

            else if ( ( iter == 1 ) && ( subiter == 0 ) )
            {
                // transition from BDF1 to BDF2
                oseenUpdateBdf2Trans( oseen, dt, uxn, uyn, uzn, ux, uy, uz,
                                      uxo, uyo, uzo, ux0, uy0, uz0,
                                      !backendNS->reusePC(),
                                      pn );
            }

            else
            {
                // BDF2
                oseenUpdateBdf2( oseen, dt, uxn, uyn, uzn, ux, uy, uz,
                                 uxo, uyo, uzo, ux0, uy0, uz0,
                                 !backendNS->reusePC(),
                                 pn );
            }

            M_timers["updateNS"].second += M_timers["updateNS"].first.elapsed();
            VLOG(1) << "[EthierSteinman] NS assembly time: "
                    << M_timers["updateNS"].first.elapsed() << "\n";
            msg.str( "" );
            msg << "t=" << time << ", subiter " << subiter << " oseen update";
            psLogger.log( msg.str() );

            // --- solve oseen
            std::cout << "[EthierSteinman] solve  oseen"
                      << ( backendNS->reusePC() ? "" : " (rebuild pc)" )
                      << "\n" << std::flush;

            if ( !backendNS->reusePC() )
                VLOG(1) << "[EthierSteinman] NS solving: rebuild pc\n";

            M_timers["solverNS"].first.restart();
            oseen.solve();
            M_timers["solverNS"].second += M_timers["solverNS"].first.elapsed();
            VLOG(1) << "[EthierSteinman] NS solving  time: "
                    << M_timers["solverNS"].first.elapsed() << "\n";
            VLOG(1) << "[EthierSteinman] NS solving  iterations: "
                    << backendNS->get_iteration() <<"\n";

            if ( backendNS->reuseFailed() )
            {
                VLOG(1) << "[EthierSteinman] NS solving: pc reuse failed\n";
                std::cout << "[EthierSteinman]    pc reuse failed\n";
            }

            if ( !backendNS->converged() )
                VLOG(1) << "[EthierSteinman] NS solving: didn't converge\n";

            msg.str( "" );
            msg << "t=" << time << ", subiter " << subiter << " oseen solve";
            psLogger.log( msg.str() );

            const element_u_type& uxnn = oseen.velocityX();
            const element_u_type& uynn = oseen.velocityY();
            //             const element_u_type& uznn = oseen.velocityZ();
            const element_p_type& pnn  = oseen.pressure();

            // --- fixpoint error calculation
            std::cout << "[EthierSteinman] fixpoint computations\n" << std::flush;
            element_u_type uxIncr( uxn );
            uxIncr -= uxnn;
            element_u_type uyIncr( uyn );
            uyIncr -= uynn;
            //             element_u_type uzIncr( uzn );
            //             uzIncr -= uznn;
            double fixpointIncU = std::sqrt( ( massU( uxIncr )( uxIncr )
                                               + massU( uyIncr )( uyIncr )
                                               //                                                + massU(uzIncr)(uzIncr)
                                             ) /
                                             ( massU( uxnn )  ( uxnn )
                                               + massU( uynn )  ( uynn )
                                               //                                                + massU(uznn)  (uznn)
                                             ) );
            element_p_type pIncr( pn );
            pIncr -= pnn;
            double fixpointIncP = std::sqrt( massP( pIncr )( pIncr ) /
                                             massP( pnn )( pnn ) );
            fixpointErr = std::sqrt( fixpointIncU*fixpointIncU /
                                     std::pow( M_meshSize, 2*( 1+uOrder ) ) +
                                     fixpointIncP*fixpointIncP /
                                     std::pow( M_meshSize, 2*( 1+pOrder ) ) );

            VLOG(1) << "[EthierSteinman] fixpoint iteration   "
                    << subiter  << "\n";
            VLOG(1) << "[EthierSteinman] fixpoint increm. u = "
                    << fixpointIncU << "\n";
            VLOG(1) << "[EthierSteinman] fixpoint increm. p = "
                    << fixpointIncP << "\n";
            VLOG(1) << "[EthierSteinman] fixpoint error     = "
                    << fixpointErr  << "\n";

            uxn = uxnn;
            uyn = uynn;
            //             uzn = uznn;
            pn  = pnn;

        } // nonlinear/subiteration loop

        // --- post processing
        std::cout << "[EthierSteinman] post processing\n";
        VLOG(1) << "[EthierSteinman] #subiter = " << subiter << "\n";

        double divError =
            std::sqrt( integrate( elements( *mesh ), M_im,
                                  vf::pow( dxv( uxn )
                                           + dyv( uyn )
                                           //                                           + dzv(uzn)
                                           ,
                                           2.0 )
                                ).evaluate()( 0,0 ) );
        VLOG(1) << "[EthierSteinman] ||div u||_2 = " << divError << "\n";

        VLOG(1) << "[EthierSteinman] p stabil. energy  = "
                << oseen.stabilizationEnergyP() << "\n";
        VLOG(1) << "[EthierSteinman] u stabil. energy  = "
                << oseen.stabilizationEnergyU() << "\n";

        // --- time shift
        uxo = ux;
        uyo = uy;
        uzo = uz;
        po  = p;

        ux = uxn;
        uy = uyn;
        uz = uzn;
        p  = pn;

        uxn = 2*ux - uxo; // extrapolation
        uyn = 2*uy - uyo;
        uzn = 2*uz - uzo;
        pn  = 2*p  - po;

        U.comp( X ) = ux;
        U.comp( Y ) = uy;
        //         U.comp(Z) = uz;

        std::cout << "[EthierSteinman] export\n" << std::flush;
        this->exportResults( iter+1, time, U, ux, uy, uz, p,
                             massU, laplaceU, massP,
                             ux0, uy0, uz0, p0 );
        msg.str( "" );
        msg << "t=" << time << " export";
        psLogger.log( msg.str() );

    } // time loop

    VLOG(1) << "[EthierSteinman] total timings:\n";
    std::map<std::string,std::pair<boost::timer,double> >::iterator it;

    for ( it=M_timers.begin(); it!=M_timers.end(); ++it )
    {
        VLOG(1) << "[EthierSteinman]   " << it->first << ": "
                << it->second.second << "\n";
    }

} // EthierSteinman::run

void
EthierSteinman::exportResults( int iter,
                               double time,
                               element_U_type& U,
                               element_u_type& ux,
                               element_u_type& uy,
                               element_u_type& uz,
                               element_p_type& p,
                               OperatorLinear<space_u_type, space_u_type, backend_type>& massU,
                               OperatorLinear<space_u_type, space_u_type, backend_type>& laplaceU,
                               OperatorLinear<space_p_type, space_p_type, backend_type>& massP,
                               element_u_type& ux0,
                               element_u_type& uy0,
                               element_u_type& uz0,
                               element_p_type& p0 )
{
    using namespace Feel::vf;

    M_timers["errors"].first.restart();
    mesh_ptr_type mesh = ux.functionSpace()->mesh();

    element_u_type psi( ux );
    FsFunctionalLinear<space_u_type, backend_type> rhsPsi( ux.functionSpace() );
    rhsPsi =
        integrate( elements( *mesh ), M_im,
                   ( dxv( uy ) - dyv( ux ) ) * id( ux )
                 ) +
        integrate( boundaryfaces( *mesh ), M_im,
                   ( -idv( uy )*Nx()+idv( ux )*Ny() ) * id( ux )
                 );
    laplaceU.applyInverse( psi, rhsPsi );

    element_u_type dux( ux );
    dux -= ux0;
    value_type uErrorL2 = massU( dux )( dux );
    value_type uErrorH1 = laplaceU( dux )( dux );

    element_u_type duy( uy );
    duy -= uy0;
    uErrorL2 += massU( duy )( duy );
    uErrorH1 += laplaceU( duy )( duy );

    element_u_type duz( uz );
    duz -= uz0;
    //     uErrorL2 += massU(duz)(duz);
    //     uErrorH1 += laplaceU(duz)(duz);

    uErrorH1 += uErrorL2;
    uErrorL2 = std::sqrt( uErrorL2 );
    uErrorH1 = std::sqrt( uErrorH1 );

    if ( M_uErrorL2 < 0.0 )
        M_uErrorL2 = uErrorL2;

    if ( M_uErrorH1 < 0.0 )
        M_uErrorH1 = uErrorH1;

    element_p_type pOne( p.functionSpace(), "pOne" );
    pOne = vf::project( p.functionSpace(), elements( *mesh ), constant( 1.0 ) );
    element_p_type dp( p );
    dp -= p0;
    double int0 = massP( pOne )( pOne );
    double int1 = massP( dp )( pOne );
    pOne *= int1/int0;
    dp -= pOne;
    double int2 = massP( dp )( dp );
    double pErrorL2 = std::sqrt( int2 ); //- int1*int1/int0 );

    if ( M_pErrorL2 < 0.0 )
        M_pErrorL2 = pErrorL2;

    double divError =
        std::sqrt( integrate( elements( *mesh ), M_im,
                              vf::pow( dxv( ux )
                                       + dyv( uy )
                                       //                                       + dzv(uz)
                                       , 2.0 )
                            ).evaluate()( 0,0 ) );

    if ( M_divError < 0.0 )
        M_divError = divError;

    M_timers["errors"].second += M_timers["errors"].first.elapsed();
    VLOG(1) << "[EthierSteinman] error comp. t.: "
            << M_timers["export"].first.elapsed()
            << "\n";

    VLOG(1) << "[EthierSteinman] ||u-u_ex||_L2 = " << uErrorL2 << "\n";
    VLOG(1) << "[EthierSteinman] ||u-u_ex||_H1 = " << uErrorH1 << "\n";
    VLOG(1) << "[EthierSteinman] ||p-p_ex||_L2 = " << pErrorL2 << "\n";
    VLOG(1) << "[EthierSteinman] ||div u ||_L2 = " << divError << "\n";

    VLOG(1) << "[EthierSteinman] ||u-u_ex||_L2 / ||u_0-u_ex||_L2 = "
            << uErrorL2/M_uErrorL2 << "\n";
    VLOG(1) << "[EthierSteinman] ||u-u_ex||_H1 / ||u_0-u_ex||_H1 = "
            << uErrorH1/M_uErrorH1 << "\n";
    VLOG(1) << "[EthierSteinman] ||p-p_ex||_L2 / ||p_0-p_ex||_L2 = "
            << pErrorL2/M_pErrorL2 << "\n";
    VLOG(1) << "[EthierSteinman] ||div u ||_L2 / ||div u_0 ||_L2 = "
            << divError/M_divError << "\n";

    M_timers["export"].first.restart();

    // -- EXPORT --
    if ( this->vm()["export"].as<int>() > 0 &&
            iter % this->vm()["export"].as<int>() == 0 )
    {
        timeset_type::step_ptrtype timeStep = M_timeSet->step( time );
        timeStep->setMesh( ux.functionSpace()->mesh() );
        timeStep->add( "U", U );
        timeStep->add( "ux", ux );
        timeStep->add( "uy", uy );
        timeStep->add( "uz", uz );
        timeStep->add( "p", p );
        timeStep->add( "uxError", dux );
        timeStep->add( "uyError", duy );
        timeStep->add( "uzError", duz );
        timeStep->add( "pError", dp );
        timeStep->add( "psi", psi );
        M_exporter->save();
    } // export

    M_timers["export"].second += M_timers["export"].first.elapsed();
    VLOG(1) << "[EthierSteinman] exporting time: "
            << M_timers["export"].first.elapsed()
            << "\n";
} // EthierSteinman::exportResults

} // Feel
