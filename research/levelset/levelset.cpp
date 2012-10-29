/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2007-06-19

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
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file levelset.cpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2007-06-19
 */

#include "levelset.hpp"

#include <feel/feelcore/pslogger.hpp>

namespace Feel
{

LevelSet::LevelSet( int argc, char** argv, AboutData const& ad )
    :
    super( argc, argv, ad ),
    M_meshSize( this->vm()["hsize"].as<double>() ),
    M_domainSize( 0.0 ),
    M_exporter( new ExporterEnsight<mesh_type>( "levelset" ) ),
    M_timeSet( new timeset_type( "levelset" ) ),
    M_timers(),
    M_im()
{
    VLOG(1) << "[LevelSet] hsize = " << M_meshSize << "\n";
    VLOG(1) << "[LevelSet] export = " << this->vm()["export"].as<int>()
            << "\n";

    M_timeSet->setTimeIncrement( this->vm()["dt"].as<double>() );
    M_exporter->addTimeSet( M_timeSet );
}

LevelSet::LevelSet( int argc,
                    char** argv,
                    AboutData const& ad,
                    po::options_description const& od )
    :
    super( argc, argv, ad, od ),
    M_meshSize( this->vm()["hsize"].as<double>() ),
    M_domainSize( 0.0 ),
    M_exporter( new ExporterEnsight<mesh_type>( "levelset" ) ),
    M_timeSet( new timeset_type( "levelset" ) ),
    M_timers(),
    M_im()
{
    VLOG(1) << "[LevelSet] hsize = " << M_meshSize << "\n";
    VLOG(1) << "[LevelSet] export = " << this->vm()["export"].as<int>()
            << "\n";

    M_timeSet->setTimeIncrement( this->vm()["dt"].as<double>() );
    M_exporter->addTimeSet( M_timeSet );
}

LevelSet::LevelSet( LevelSet const& tc )
    :
    super( tc ),
    M_meshSize( tc.M_meshSize ),
    M_domainSize( 0.0 ),
    M_exporter( new ExporterEnsight<mesh_type>( "levelset" ) ),
    M_timeSet( new timeset_type( "levelset" ) ),
    M_timers( tc.M_timers ),
    M_im()
{
    VLOG(1) << "[LevelSet] hsize = " << M_meshSize << "\n";
    VLOG(1) << "[LevelSet] export = " << this->vm()["export"].as<int>()
            << "\n";

    M_timeSet->setTimeIncrement( this->vm()["dt"].as<double>() );
    M_exporter->addTimeSet( M_timeSet );
}

void
LevelSet::operator()()
{
    run();
}

LevelSet::mesh_ptr_type
LevelSet::createMesh( double meshSize )
{
    M_timers["mesh"].first.restart();
    mesh_ptr_type mesh( new mesh_type );

    GmshHypercubeDomain<Dim,1,Dim,ENTITY> td;
    td.setCharacteristicLength( meshSize );
    td.setX( std::make_pair( -1.0, 1.0 ) );
    td.setY( std::make_pair( -1.0, 1.0 ) );

    if ( Dim > 2 )
    {
        td.setZ( std::make_pair( -1.0, 1.0 ) );
    }

    std::string fname = td.generate( ENTITY<Dim,1,Dim>::name() );

    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );

    M_timers["mesh"].second = M_timers["mesh"].first.elapsed();
    VLOG(1) << "[LevelSet] createMesh(): "
            << M_timers["mesh"].second << "\n";

    return mesh;

} // LevelSet::createMesh


void
LevelSet::run()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    this->changeRepository( boost::format( "%1%/h_%2%/dt_%3%" )
                            % this->about().appName()
                            % M_meshSize
                            % this->vm()["dt"].as<double>() );
    this->setLogs();

    PsLogger psLogger( "ps.log" );
    psLogger.log( "t=0, start" );

    using namespace Feel::vf;

    /*
     * First we create the mesh
     */
    mesh_ptr_type mesh = createMesh( M_meshSize );
    M_domainSize = integrate( elements( mesh ), M_im,
                              constant( 1.0 ) ).evaluate()( 0,0 );
    psLogger.log( "t=0, meshed" );

    /*
     * The function spaces and some associate elements are then defined
     */
    M_timers["init"].first.restart();
    space_p_ptrtype space_p = space_p_type::New( mesh );
    space_i_ptrtype space_i = space_i_type::New( mesh );
    //space_i->dof()->showMe();
    psLogger.log( "t=0, spaces" );

    element_p_type phi  ( space_p, "phi" );
    element_p_type phio ( space_p, "phio" );
    element_p_type psi  ( space_p, "psi" );
    element_i_type kappa( space_i, "kappa" );
    psLogger.log( "t=0, elements" );

    VLOG(1) << "[LevelSet] h = " << M_meshSize << "\n";
    VLOG(1) << "[LevelSet] N = " << space_p->nDof() << "\n";

    backend_ptrtype backend( backend_type::build( this->vm() ) );
#warning TODO
#if 0
    backend->set_solver_type( this->vm()["solver"].as<std::string>() );
    backend->set_noisy( this->vm()["verbose"].as<int>() );
    backend->set_maxiter( this->vm()["maxiter"].as<int>() );
    backend->set_tol( this->vm()["tolerance"].as<double>() );
    backend->set_fillin( this->vm()["fillin"].as<int>() );
    backend->set_threshold( this->vm()["threshold"].as<double>() );
#endif
    backendS_ptrtype backendSymmP1( backend_type::build( this->vm() ) );
#warning TODO
#if 0
    backendSymmP1->set_solver_type( "cg" );
    backendSymmP1->set_noisy( this->vm()["verbose"].as<int>() );
    backendSymmP1->set_maxiter( this->vm()["maxiter"].as<int>() );
    backendSymmP1->set_tol( this->vm()["tolerance"].as<double>() );
    backendSymmP1->set_fillin( this->vm()["fillin"].as<int>() );
    backendSymmP1->set_threshold( this->vm()["threshold"].as<double>() );
    backendSymmP1->set_symmetric( true );
#endif
    backendS_ptrtype backendSymmP0( backend_type::build( this->vm() ) );
#warning TODO
#if 0
    backendSymmP0->set_solver_type( "cg" );
    backendSymmP0->set_noisy( this->vm()["verbose"].as<int>() );
    backendSymmP0->set_maxiter( this->vm()["maxiter"].as<int>() );
    backendSymmP0->set_tol( this->vm()["tolerance"].as<double>() );
    backendSymmP0->set_fillin( this->vm()["fillin"].as<int>() );
    backendSymmP0->set_threshold( this->vm()["threshold"].as<double>() );
    backendSymmP0->set_symmetric( true );
#endif
    psLogger.log( "t=0, backends" );

    typedef __typeof__( elements( mesh ) ) IteratorRange;
    ReinitializerFMS<space_p_type, IteratorRange> reinitializerFMS( space_p, elements( mesh ) );
    ReinitializerILP<space_p_type, ENTITY> reinitializerILP( space_p, backendSymmP1 );
    Indicator<space_p_type, ENTITY> indicator( space_p, space_i, backendSymmP0 );
    psLogger.log( "t=0, reinitializers" );

    // -- initial condition
    value_type x0 = 0.5;
    value_type y0 = 0.5; // 0.2
    value_type radius = 0.25; // 0.1
    value_type slotWidth = 0.1;
    value_type slotDepth = 0.0;
    value_type pi = 4.0 * math::atan( value_type( 1.0 ) );

    // bubble
    //     AUTO( phi0, pow(radius,2.0) - pow(Px()-x0,2.0) - pow(Py()-y0,2.0) );
    //     AUTO( phi1, radius - sqrt( pow( Px()-x0, 2.0 ) + pow( Py()-y0, 2.0 ) ) );
    //     value_type mass0 = pi*std::pow(radius,2.0);

    // Zalesak slotted disk
    //     AUTO( phi0, max( ( min( slotWidth+(Px()-x0),
    //                             (min( slotWidth-(Px()-x0),
    //                                   slotDepth-(Py()-y0) ) ) ) ),
    //                      sqrt( pow(Px()-x0,2.0) + pow(Py()-y0,2.0) ) - radius) );
    //     value_type mass0 = (pi-std::asin(slotWidth/radius))*std::pow(radius,2.0)
    //         - slotWidth *
    //         ( std::sqrt(std::pow(radius,2.0)-std::pow(slotWidth,2.0)) +
    //           2.0*slotDepth );

    // circle
    //     AUTO( phi0, ( pow( pow(Px()-x0,2.0) + pow(Py()-y0,2.0), 2.0*Px()+0.5 ) -
    //                   pow( pow(radius ,2.0)                   , 2.0*Px()+0.5 ) ) );
    //     AUTO( phi1, 2*(sqrt(pow(Px()-x0,2.0)+pow(Py()-y0,2.0))-radius) );
    //     value_type mass0 = pi*std::pow(radius,2.0);

    // circle Fujima-Ohmori on (-1,1)^2
    AUTO( phi0,
          ( exp( -1.0 ) -
            vf::exp( -9.0*trans( P()-1./3.*oneY() )*( P()-1./3.*oneY() ) ) )
          / exp( -1.0 ) );
    AUTO( phi1, vf::sqrt( trans( P()-1./3.*oneY() )*( P()-1./3.*oneY() ) ) - 1./3. );
    value_type mass0 = pi*std::pow( 1./3.,Dim )*( 1.+Dim )/3.;

    // straight Fujima-Ohmori line y=0 on (-1,1)^2
    //     AUTO( phi0, Py() );
    //     AUTO( phi1, Py() );
    //     value_type mass0 = 2.0;

    // straight line
    //     AUTO( phi0, (Px()-x0)*exp(10*Py()) );
    //     AUTO( phi0, 2*(Px()-x0) );
    //     AUTO( phi1, Px()-x0 );
    //     value_type mass0 = x0;

    // square
    //     AUTO( phi0, 2*max( (max(-radius-(Px()-x0),-radius+(Px()-x0))),
    //                        (max(-radius-(Py()-y0),-radius+(Py()-y0))) ) );
    //     AUTO( phi1, max(max(-radius-(Px()-x0),-radius+(Px()-x0)),
    //                     max(-radius-(Py()-y0),-radius+(Py()-y0)));
    //     value_type mass0 = 4*pow(radius,2.0);

    // two circles intersecting
    //     AUTO( phi0, min(sqrt( pow(Px()-x0+radius/2,2.0) + pow(Py()-y0,2.0) )
    //                     -radius,
    //                     sqrt( pow(Px()-x0-radius/2,2.0) + pow(Py()-y0,2.0) )
    //                     -radius) );
    //     AUTO( phi1, min(sqrt( pow(Px()-x0+radius/2,2.0) + pow(Py()-y0,2.0) )
    //                     -radius,
    //                     sqrt( pow(Px()-x0-radius/2,2.0) + pow(Py()-y0.2.0) )
    //                     -radius);
    //     value_type mass0 = pow(radius,2.0)*(4*pi/3+sqrt(3.)/2);

    psi = vf::project( space_p, elements( mesh ), phi0 );

    double massBefore = mass( psi );
    VLOG(1) << "[LevelSet] mass before reinit = " << massBefore << "\n";
    VLOG(1) << "[LevelSet]   rel. mass error  = " << massBefore/mass0-1.0
            << "\n";
    M_timers["init"].second = M_timers["init"].first.elapsed();

    M_timers["reinit"].first.restart();
    indicator.update( psi );
    VLOG(1) << "[LevelSet] indicator update : "
            << M_timers["reinit"].first.elapsed() << "\n";
    kappa = indicator.indicatorGamma();
    phi = reinitializerILP( psi, kappa );
    VLOG(1) << "[LevelSet] + ILP            : "
            << M_timers["reinit"].first.elapsed() << "\n";
    phi = reinitializerFMS( phi );
    VLOG(1) << "[LevelSet] + FMS            : "
            << M_timers["reinit"].first.elapsed() << "\n";
    M_timers["reinit"].second = M_timers["reinit"].first.elapsed();

    M_timers["init"].first.restart();

    statsAfterReinit( psi, phi, massBefore, mass0 );

    AdvReact<space_p_type, imOrder, ENTITY> advreact( space_p, backend );
    advreact.set_stabcoeff( this->vm()["stabcoeff"].as<double>() );
    psLogger.log( "t=0, advreact" );

    //     AUTO( betax, 2*pi*(0.5-Py()) );
    //     AUTO( betay, 2*pi*(Px()-0.5) );
    //     AUTO( betax,  sin(2.0*pi*Py())*sin(2.0*pi*Px()) );
    //     AUTO( betay, -sin(3.0*pi*Px())*sin(1.0*pi*Py()) );
    AUTO( betax, -2.0*pi*Py()*( 1.0-vf::pow( Px(),2.0 ) ) );
    AUTO( betay,  2.0*pi*Px()*( 1.0-vf::pow( Py(),2.0 ) ) );

    //     AUTO( beta, betax*oneX()+betay*oneY() );

    element_p_type vx = project( space_p, elements( mesh ), betax );
    element_p_type vy = project( space_p, elements( mesh ), betay );
    //     AUTO( beta, idv(vx)*oneX()+idv(vy)*oneY() );

    double time        = 0;
    this->exportResults( 0, time, phi, psi, kappa, vx, vy );
    double dt          = this->vm()["dt"].as<double>();
    time               = dt;
    phio               = phi;

    const value_type theta = this->vm()["theta"].as<double>();
    M_timers["init"].second += M_timers["init"].first.elapsed();
    const double finalTime = this->vm()["ft"].as<double>();

    // square of minimal gradient criterion for reinitialization
    double mingrad2 = this->vm()["mingrad"].as<double>();
    mingrad2 *= mingrad2;

    int reinitEvery = this->vm()["reinitevery"].as<int>();
    bool doReinit = true;
    std::stringstream msg;

    // --- Time loop
    for ( int iter = 0; time-dt/2 < finalTime; ++iter, time += dt )
    {
        // --- update advreact
        std::cout << "[LevelSet] update advreact\n" << std::flush;
        M_timers["update"].first.restart();
        double sign = time-dt/2 > finalTime/2 ? -1 : 1;

        if ( this->vm().count( "bdf2" ) && ( doReinit ) )
        {
            advReactUpdateBdf2( advreact, dt, sign, vx, vy, phi, phio, false );//!backend->reusePC() );
        }

        else
        {
            advReactUpdateCN( advreact, dt, theta, sign, vx, vy, phi, false );// !backend->reusePC() );
        }

        M_timers["update"].second += M_timers["update"].first.elapsed();
        VLOG(1) << "[LevelSet] assembly time: "
                << M_timers["update"].first.elapsed() << "\n";
        msg.str( "" );
        msg << "t=" << time << " advreact update";
        psLogger.log( msg.str() );

        // --- solve advreact
        std::cout << "[LevelSet] solve  advreact\n" << std::flush;
        M_timers["solve"].first.restart();
        advreact.solve();
        M_timers["solve"].second += M_timers["solve"].first.elapsed();
        VLOG(1) << "[LevelSet] solving  time: "
                << M_timers["solve"].first.elapsed() << "\n";
        msg.str( "" );
        msg << "t=" << time << " advreact solve";
        psLogger.log( msg.str() );

        psi = advreact.phi();

        phio = phi;

        VLOG(1) << "[LevelSet] t = " << time << "\n";
        std::cout << "[LevelSet] t = " << time << "\n";

        massBefore = mass( psi );
        VLOG(1) << "[LevelSet] mass before reinit = " << massBefore << "\n";
        VLOG(1) << "[LevelSet]   rel. mass error  = "
                << massBefore/mass0-1.0 << "\n";

        // --- reinitialize
        if ( reinitEvery == 0 )
        {
            double reinitIndicator =
                integrate( elements( mesh ), M_im,
                           chi( ( gradv( psi )*trans( gradv( psi ) ) <
                                  mingrad2 )
                                && ( pow( idv( psi ),2.0 ) < pow( h(),2.0 )*
                                     ( gradv( psi )*trans( gradv( psi ) ) ) )
                              ) ).evaluate()( 0,0 );
            doReinit = ( reinitIndicator > 0.0 );
        }

        else
        {
            doReinit = ( ( iter+1 ) % reinitEvery == 0 );
        }

        if ( doReinit )
        {
            std::cout << "[LevelSet] reinitialize\n" << std::flush;
            M_timers["reinit"].first.restart();
            indicator.update( psi );
            kappa = indicator.indicatorGamma();
            phi = reinitializerILP( psi, kappa );
            phi = reinitializerFMS( phi );
            M_timers["reinit"].second +=
                M_timers["reinit"].first.elapsed();

            statsAfterReinit( psi, phi, massBefore, mass0 );
            msg.str( "" );
            msg << "t=" << time << " reinit";
            psLogger.log( msg.str() );
        }

        else
        {
            phi = psi;
        }

        this->exportResults( iter+1, time, phi, psi, kappa, vx, vy );
        msg.str( "" );
        msg << "t=" << time << " export";
        psLogger.log( msg.str() );

    } // time loop

    double forthBackL2 =
        std::sqrt( integrate( elements( mesh ), M_im,
                              vf::pow( idv( phi ) - phi1, 2.0 ) ).evaluate()( 0,0 ) );
    VLOG(1) << "[LevelSet] forth-back L2 error = " << forthBackL2 << "\n";

    double signChgErr = integrate( elements( mesh ), M_im,
                                   chi( idv( phi )*phi1<0 )
                                 ).evaluate()( 0,0 );
    VLOG(1) << "[LevelSet] global sgn chg err = " << signChgErr << "\n";
    VLOG(1) << "[LevelSet] gl.  sgn chg err/h = "
            << signChgErr/M_meshSize << "\n";

    VLOG(1) << "[LevelSet] total timings:\n";
    VLOG(1) << "[LevelSet]   init:     " << M_timers["init"].second << "\n";
    VLOG(1) << "[LevelSet]   mesh:     " << M_timers["mesh"].second << "\n";
    VLOG(1) << "[LevelSet]   assembly: " << M_timers["update"].second << "\n";
    VLOG(1) << "[LevelSet]   solving:  " << M_timers["solve"].second << "\n";
    VLOG(1) << "[LevelSet]   reinit.:  " << M_timers["reinit"].second << "\n";

} // LevelSet::run

double
LevelSet::mass( const element_p_type& psi )
{
    using namespace Feel::vf;
    return integrate( elements( *( psi.functionSpace()->mesh() ) ), M_im,
                      chi( idv( psi )<0 )
                    ).evaluate()( 0,0 );
}

void
LevelSet::exportResults( int iter,
                         double time,
                         element_p_type& phi,
                         element_p_type& psi,
                         element_i_type& kappa,
                         element_p_type& vx,
                         element_p_type& vy )
{
    M_timers["export"].first.restart();

    // -- EXPORT --
    if ( this->vm()["export"].as<int>() > 0 &&
            iter % this->vm()["export"].as<int>() == 0 )
    {
        timeset_type::step_ptrtype timeStep = M_timeSet->step( time );
        timeStep->setMesh( phi.functionSpace()->mesh() );
        timeStep->add( "phi", phi );
        timeStep->add( "psi", psi );
        timeStep->add( "kappa", kappa );
        timeStep->add( "vx", vx );
        timeStep->add( "vy", vy );
        M_exporter->save();
    } // export

    M_timers["export"].second += M_timers["export"].first.elapsed();
    VLOG(1) << "[LevelSet] exporting time: "
            << M_timers["export"].second << "\n";
} // LevelSet::export

void
LevelSet::statsAfterReinit( const element_p_type& psi,
                            const element_p_type& phi,
                            double massBefore,
                            double mass0 )
{
    using namespace Feel::vf;

    double massAfter = mass( phi );
    VLOG(1) << "[LevelSet] mass after  reinit = " << massAfter << "\n";
    VLOG(1) << "[LevelSet]   rel. mass error  = " << massAfter/mass0-1.0
            << "\n";
    VLOG(1) << "[LevelSet]   rel. reinit err. = " << massAfter/massBefore-1.0
            << "\n";

    double signChgErr =
        integrate( elements( *( psi.functionSpace()->mesh() ) ), M_im,
                   chi( idv( phi )*idv( psi )<0 )
                 ).evaluate()( 0,0 );
    VLOG(1) << "[LevelSet] sign change error  = " << signChgErr << "\n";
    VLOG(1) << "[LevelSet] sign chg. error/h  = " << signChgErr/M_meshSize
            << "\n";

    double gradMag =
        integrate( elements( *( psi.functionSpace()->mesh() ) ), M_im,
                   sqrt( gradv( phi )*trans( gradv( phi ) ) )
                 ).evaluate()( 0,0 ) / M_domainSize;
    VLOG(1) << "[LevelSet] mean gradient magnitude = " << gradMag
            << "\n";

    double distDist =
        std::sqrt( integrate( elements( *( psi.functionSpace()->mesh() ) ), M_im,
                              pow( sqrt( gradv( phi )*trans( gradv( phi ) ) ) - 1.0,
                                   2.0 )
                            ).evaluate()( 0,0 ) );
    VLOG(1) << "[LevelSet] distance from dist = " << distDist << "\n";

} // LevelSet::statsAfterReinit

} // Feel
