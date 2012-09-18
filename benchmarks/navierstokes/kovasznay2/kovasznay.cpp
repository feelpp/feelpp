/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2007-01-25

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
   \file kovasznay.cpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2007-01-25
 */

#include <sstream>

#include "kovasznay.hpp"

#include <feel/feelfilters/importergmsh.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>

#include <feel/feelcore/pslogger.hpp>
#include <feel/feeldiscr/operatorlinear.hpp>

namespace Feel
{

Kovasznay::Kovasznay( int argc, char** argv, AboutData const& ad )
    :
    super( argc, argv, ad ),
    M_meshSize( this->vm()["hsize"].as<double>() ),
    M_nu( this->vm()["nu"].as<double>() ),
    M_exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) ),
    M_timers(),
    M_im(),
    M_uErrorL2( -1.0 ),
    M_uErrorH1( -1.0 ),
    M_pErrorL2( -1.0 ),
    M_divError( -1.0 )
{
    VLOG(1) << "[Kovasznay] hsize       = " << M_meshSize << "\n";
    VLOG(1) << "[Kovasznay] nu          = " << M_nu << "\n";
    VLOG(1) << "[Kovasznay] export      = "
            << this->vm()["doexport"].as<int>() << "\n";

}

Kovasznay::Kovasznay( int argc,
                      char** argv,
                      AboutData const& ad,
                      po::options_description const& od )
    :
    super( argc, argv, ad, od ),
    M_meshSize( this->vm()["hsize"].as<double>() ),
    M_nu( this->vm()["nu"].as<double>() ),
    M_exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) ),
    M_timers(),
    M_im(),
    M_uErrorL2( -1.0 ),
    M_uErrorH1( -1.0 ),
    M_pErrorL2( -1.0 ),
    M_divError( -1.0 )
{
    VLOG(1) << "[Kovasznay] hsize   = " << M_meshSize << "\n";
    VLOG(1) << "[Kovasznay] nu      = " << M_nu << "\n";
    VLOG(1) << "[Kovasznay] export  = "
            << this->vm()["doexport"].as<int>() << "\n";

}

Kovasznay::Kovasznay( Kovasznay const& tc )
    :
    super( tc ),
    M_meshSize( tc.M_meshSize ),
    M_nu( tc.M_nu ),
    M_exporter( Exporter<mesh_type>::New( "kovasznay" ) ),
    M_timers( tc.M_timers ),
    M_im(),
    M_uErrorL2( -1.0 ),
    M_uErrorH1( -1.0 ),
    M_pErrorL2( -1.0 ),
    M_divError( -1.0 )
{
    VLOG(1) << "[Kovasznay] hsize   = " << M_meshSize << "\n";
    VLOG(1) << "[Kovasznay] nu      = " << M_nu << "\n";
    VLOG(1) << "[Kovasznay] export  = "
            << this->vm()["doexport"].as<int>() << "\n";

}

Kovasznay::mesh_ptr_type
Kovasznay::createMesh( double meshSize )
{
    M_timers["mesh"].first.restart();
#if 0
    mesh_ptr_type mesh( new mesh_type );

    GmshHypercubeDomain<Dim,1,ENTITY> td;
    td.setCharacteristicLength( meshSize );
    td.setX( std::make_pair( -0.5, 1. ) );
    td.setY( std::make_pair( -0.5, 1.5 ) );
    ImporterGmsh<mesh_type>
    import( td.generate( ENTITY<Dim,1,Dim>::name().c_str() ) );
    mesh->accept( import );
#else
    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _desc=domain( _name="square",
                                        _shape="hypercube",
                                        _usenames=false,
                                        _dim=2,
                                        _h=meshSize,
                                        _xmin=-0.5,_xmax=1.,
                                        _ymin=-0.5,_ymax=1.5 ) );

#endif
    M_timers["mesh"].second = M_timers["mesh"].first.elapsed();
    VLOG(1) << "[timer] createMesh(): " << M_timers["mesh"].second << "\n";
    return mesh;
} // Kovasznay::createMesh


void
Kovasznay::run()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    this->changeRepository( boost::format( "%1%/h_%2%/Re_%3%" )
                            % this->about().appName()
                            % M_meshSize
                            % ( 1./M_nu ) );
    this->setLogs();

    PsLogger psLogger( "ps.log" );
    psLogger.log( "t=0, start" );

    using namespace Feel::vf;

    /*
     * First we create the mesh
     */
    mesh_ptr_type mesh = createMesh( M_meshSize );
    psLogger.log( "t=0, meshed" );

    /*
     * The function spaces and some associate elements are then defined
     */
    M_timers["init"].first.restart();
    space_ptrtype space = space_type::New( mesh );
    space_U_ptrtype space_U = space->functionSpace<0>();
    space_p_ptrtype space_p = space->functionSpace<1>();
    space_i_ptrtype space_i = space_i_type::New( mesh );
    //space_u->dof()->showMe();
    psLogger.log( "t=0, spaces" );

    VLOG(1) << "[Kovasznay] velocity dofs total         "
            << space_U->nDof() << "\n";
    VLOG(1) << "[Kovasznay] pressure dofs               "
            << space_p->nDof() << "\n";
    VLOG(1) << "[Kovasznay] total    dofs               "
            << space_U->nDof() + space_p->nDof() << "\n";
#if 0
    element_u_type ux   ( space_u, "ux" );
    element_u_type uy   ( space_u, "uy" );
    element_u_type rx   ( space_u, "rx" );
    element_u_type ry   ( space_u, "ry" );
    element_u_type uxo  ( space_u, "uxo" );
    element_u_type uyo  ( space_u, "uyo" );
    element_u_type rxo  ( space_u, "rxo" );
    element_u_type ryo  ( space_u, "ryo" );
    element_u_type uxn  ( space_u, "uxn" );
    element_u_type uyn  ( space_u, "uyn" );
    element_u_type dux  ( space_u, "dux" );
    element_u_type duy  ( space_u, "duy" );
    element_U_type U    ( space_U, "U" );
    element_U_type Un   ( space_U, "Un" );
    element_p_type p    ( space_p, "p" );
    element_p_type rp   ( space_p, "rp" );
    element_p_type po   ( space_p, "po" );
    element_p_type rpo  ( space_p, "rpo" );
    element_p_type pn   ( space_p, "pn" );
    element_p_type dp   ( space_p, "dp" );
    //     element_p_type pl   ( space_p, "pl" );
    element_i_type phi  ( space_i, "phi" );
    psLogger.log( "t=0, elements" );

    // -- exact solution
    value_type pi = 4.0 * math::atan( value_type( 1.0 ) );
    value_type lambda = 1./( 2.*M_nu ) - std::sqrt( 1./( 4.*M_nu*M_nu ) + 4.*pi*pi );
    AUTO( uxe, 1. - exp( lambda * Px() ) * cos( 2.*pi*Py() ) );
    AUTO( uye, lambda/( 2.*pi ) * exp( lambda * Px() ) * sin( 2.*pi*Py() ) );
    AUTO( pe, 0.5*( 1.-exp( 2.*lambda*Px() ) ) );
    //     value_type eps = type_traits<value_type>::epsilon();
    //     AUTO( uxe, chi(Py() > 1.-eps ) );
    //     AUTO( uye, constant(0.) );
    //     AUTO( pe, constant(0.) );

    // -- use projection of exact solution as initial guess!
    U = vf::project( space_U, elements( *mesh ), uxe*oneX()+uye*oneY() );
    p  = vf::project( space_p, elements( *mesh ), pe );

    // -- dummy function to define empty iterator range
    phi = project( space_i, elements( *mesh ), constant( 0. ) );
    mesh->updateMarker2( phi );

    this->exportResults( -1, U, p );

    backend_ptrtype backendNS  ( new backend_type( this->vm(), "oseen" ) );
    backend_ptrtype backendSymm( new backend_type( this->vm(), "symm"  ) );
    psLogger.log( "t=0, backends" );

    std::set<flag_type> dirichletFlags;
    std::set<flag_type> neumannFlags;

    for ( flag_type flag=1; flag<=4; ++flag )
    {
        dirichletFlags.insert( flag );
    }

    Oseen<space_type, imOrder, ENTITY>
    oseen( space, backendNS, dirichletFlags, neumannFlags, this->vm() );
    //     oseen.decouplePstab( p, this->vm()["stabtheta"].as<double>() );
    psLogger.log( "t=0, oseen" );

    double fixpointTol = this->vm()["fixpointtol"].as<double>();
    int    maxSubIter  = this->vm()["maxsubiter" ].as<int>();

    // mass matrix on velocity space, for L2 norms
    OperatorLinear<space_u_type, space_u_type, backend_type>
    massU( space_u, space_u, backendSymm );
    massU = integrate( elements( *mesh ), M_im, id( ux )*idt( ux ) );

    // mass matrix on pressure space, for L2 norms
    OperatorLinear<space_p_type, space_p_type, backend_type>
    massP( space_p, space_p, backendSymm );
    massP = integrate( elements( *mesh ), M_im, id( p )*idt( p ) );
    psLogger.log( "t=0, mass matrices" );

    M_timers["init"].second = M_timers["init"].first.elapsed();

    // -- first nonlinear iteration with fixed omega and initialization of
    //    oseen operator

    std::stringstream msg;
    uint32_type subiter = 0;

    std::cout << "[Kovasznay] update oseen (rebuild ip)\n";
    M_timers["updateNS"].first.restart();
    std::cout << "[Kovasznay] update oseen (rebuild ip)\n" << std::flush;

    oseenUpdateInit( oseen, U );

    M_timers["updateNS"].second += M_timers["updateNS"].first.elapsed();
    VLOG(1) << "[Kovasznay] NS assembly time: "
            << M_timers["updateNS"].first.elapsed() << "\n";
    msg.str( "" );
    msg << "subiter " << subiter << " oseen update";
    psLogger.log( msg.str() );

    // --- solve oseen
    std::cout << "[Kovasznay] solve  oseen"
              << ( backendNS->reusePC() ? "" : " (rebuild pc)" )
              << "\n" << std::flush;

    if ( !backendNS->reusePC() )
        VLOG(1) << "[Kovasznay] NS solving: rebuild pc\n";

    M_timers["solverNS"].first.restart();
    oseen.solve();
    M_timers["solverNS"].second += M_timers["solverNS"].first.elapsed();
    VLOG(1) << "[Kovasznay] NS solving  time: "
            << M_timers["solverNS"].first.elapsed() << "\n";
    VLOG(1) << "[Kovasznay] NS solving  iterations: "
            << backendNS->get_iteration() <<"\n";

    if ( backendNS->reuseFailed() )
    {
        VLOG(1) << "[Kovasznay] NS solving: pc reuse failed\n";
        std::cout << "[Kovasznay]    pc reuse failed\n";
    }

    if ( !backendNS->converged() )
        VLOG(1) << "[Kovasznay] NS solving: didn't converge\n";

    msg.str( "" );
    msg << "subiter " << subiter << " oseen solve";
    psLogger.log( msg.str() );

    Un  = oseen.velocity();
    uxn = Un.comp( X );
    uyn = Un.comp( Y );
    pn  = oseen.pressure();
    //     const element_U_type& Un  = oseen.velocity();
    //     const element_u_type& uxn = Un.comp(X);
    //     const element_u_type& uyn = Un.comp(Y);
    //     const element_p_type& pn  = oseen.pressure();

    //         dux  = ux;          duy  = uy;          dp  = p;
    //         ux   = uxn;         uy   = uyn;         p   = pn;
    //         dux -= ux;          duy -= uy;          dp -= p;

    value_type omegaMin = 1.e-3;
    value_type omegaLow = 1.e-2;
    value_type omegaHigh= 1.e+0;
    value_type omegaMax = 1.e+1;
    value_type omega = omegaLow;

    rx   = ux;
    ry   = uy;
    rp  = p;
    rx  -= uxn;
    ry  -= uyn;
    rp -= pn;
    uxo  = ux;
    uyo  = uy;
    po  = p;
    rxo  = rx;
    ryo  = ry;
    rpo = rp;
    dux  = rx;
    duy  = ry;
    dp  = rp;
    dux *= -omega;
    duy *= -omega;
    dp *= -omega;
    ux  += dux;
    uy  += duy;
    p  += dp;

    double fixpointIncU = std::sqrt( ( massU( dux )( dux ) + massU( duy )( duy ) ) /
                                     ( massU( ux )( ux ) + massU( uy )( uy ) ) );
    double fixpointIncP = std::sqrt( massP( dp )( dp ) / massP( p )( p ) );
    double fixpointErr  = std::sqrt( fixpointIncU*fixpointIncU /
                                     std::pow( M_meshSize, 2*( 1+uOrder ) ) +
                                     fixpointIncP*fixpointIncP /
                                     std::pow( M_meshSize, 2*( 1+pOrder ) ) );

    VLOG(1) << "[Kovasznay] fixpoint iteration   " << subiter << "\n";
    VLOG(1) << "[Kovasznay] fixpoint increm. u = " << fixpointIncU << "\n";
    VLOG(1) << "[Kovasznay] fixpoint increm. p = " << fixpointIncP << "\n";
    VLOG(1) << "[Kovasznay] fixpoint error     = " << fixpointErr  << "\n";

    this->exportResults( subiter, U, p );
    psLogger.log( "doexport" );

    //     pl = p;
    //     oseen.decouplePstab( pl, this->vm()["stabtheta"].as<double>() );

    for ( subiter = 1;
            ( fixpointErr>fixpointTol ) && ( subiter < maxSubIter );
            ++subiter )
    {
        // --- update oseen
        std::cout << "[Kovasznay] update oseen"
                  << ( backendNS->reusePC() ? "" : " (rebuild ip)" )
                  << "\n" << std::flush;

        M_timers["updateNS"].first.restart();
        oseenUpdateIter( oseen, U, !backendNS->reusePC() );

        M_timers["updateNS"].second += M_timers["updateNS"].first.elapsed();
        VLOG(1) << "[Kovasznay] NS assembly time: "
                << M_timers["updateNS"].first.elapsed() << "\n";
        msg.str( "" );
        msg << "subiter " << subiter << " oseen update";
        psLogger.log( msg.str() );

        // --- solve oseen
        std::cout << "[Kovasznay] solve  oseen"
                  << ( backendNS->reusePC() ? "" : " (rebuild pc)" )
                  << "\n" << std::flush;

        if ( !backendNS->reusePC() )
            VLOG(1) << "[Kovasznay] NS solving: rebuild pc\n";

        M_timers["solverNS"].first.restart();
        oseen.solve();
        M_timers["solverNS"].second += M_timers["solverNS"].first.elapsed();
        VLOG(1) << "[Kovasznay] NS solving  time: "
                << M_timers["solverNS"].first.elapsed() << "\n";
        VLOG(1) << "[Kovasznay] NS solving  iterations: "
                << backendNS->get_iteration() <<"\n";

        if ( backendNS->reuseFailed() )
        {
            VLOG(1) << "[Kovasznay] NS solving: pc reuse failed\n";
            std::cout << "[Kovasznay]    pc reuse failed\n";
        }

        if ( !backendNS->converged() )
            VLOG(1) << "[Kovasznay] NS solving: didn't converge\n";

        msg.str( "" );
        msg << "subiter " << subiter << " oseen solve";
        psLogger.log( msg.str() );

        Un  = oseen.velocity();
        pn  = oseen.pressure();

        //         dux  = ux;          duy  = uy;          dp  = p;
        //         ux   = uxn;         uy   = uyn;         p   = pn;
        //         dux -= ux;          duy -= uy;          dp -= p;

        rx   = ux;
        ry   = uy;
        rp   = p;
        rx  -= uxn;
        ry  -= uyn;
        rp  -= pn;

        element_u_type drx( rx );
        element_u_type dry( ry );
        element_p_type drp( rp );
        drx -= rxo;
        dry -= ryo;
        drp -= rpo;

        double denom = massU( drx )( drx ) + massU( dry )( dry ) + massP( drp )( drp );

        if ( denom == 0.0 )
        {
            omega = omegaLow;
            VLOG(1) << "[Kovasznay] omega = Inf -> " << omegaLow << "\n";
        }

        else
        {
            double num = massU( drx )( dux ) + massU( dry )( duy ) + massP( drp )( dp );
            omega = num / denom;

            if ( omega < omegaMin )
            {
                VLOG(1) << "[Kovasznay] omega = " << omega << " -> "
                        << omegaLow << "\n";
                omega = omegaLow;
            }

            else if ( omega > omegaMax )
            {
                VLOG(1) << "[Kovasznay] omega = " << omega << " -> "
                        << omegaHigh << "\n";
                omega = omegaHigh;
            }

            else
            {
                VLOG(1) << "[Kovasznay] omega = " << omega << "\n";
            }
        }

        uxo  = ux;
        uyo  = uy;
        po   = p;
        rxo  = rx;
        ryo  = ry;
        rpo  = rp;
        dux  = rx;
        duy  = ry;
        dp   = rp;
        dux *= -omega;
        duy *= -omega;
        dp  *= -omega;
        ux  += dux;
        uy  += duy;
        p   += dp;

        fixpointIncU = std::sqrt( ( massU( drx )( drx ) + massU( dry )( dry ) ) /
                                  ( massU( ux )( ux ) + massU( uy )( uy ) ) );
        fixpointIncP = std::sqrt( massP( rp )( rp ) / massP( p )( p ) );
        fixpointErr  = std::sqrt( fixpointIncU*fixpointIncU /
                                  std::pow( M_meshSize, 2*( 1+uOrder ) ) +
                                  fixpointIncP*fixpointIncP /
                                  std::pow( M_meshSize, 2*( 1+pOrder ) ) );

        VLOG(1) << "[Kovasznay] fixpoint iteration   " << subiter << "\n";
        VLOG(1) << "[Kovasznay] fixpoint increm. u = " << fixpointIncU << "\n";
        VLOG(1) << "[Kovasznay] fixpoint increm. p = " << fixpointIncP << "\n";
        VLOG(1) << "[Kovasznay] fixpoint error     = " << fixpointErr  << "\n";

        p = oseen.pressure();
        this->exportResults( subiter, U, p );
        psLogger.log( "doexport" );

    } // nonlinear/subiteration loop

    // --- post processing
    VLOG(1) << "[Kovasznay] #subiter = " << subiter << "\n";
#endif
    VLOG(1) << "[Kovasznay] total timings:\n";
    std::map<std::string,std::pair<boost::timer,double> >::iterator it;

    for ( it=M_timers.begin(); it!=M_timers.end(); ++it )
    {
        VLOG(1) << "[Kovasznay]   " << it->first << ": "
                << it->second.second << "\n";
    }

} // Kovasznay::run

void
Kovasznay::exportResults( int iter,
                          element_U_type& U,
                          element_p_type& p )
{
    // --- error calculations
    using namespace Feel::vf;
    value_type pi = 4.0 * math::atan( value_type( 1.0 ) );
    value_type lambda = 1./( 2.*M_nu ) - std::sqrt( 1./( 4.*M_nu*M_nu ) + 4.*pi*pi );
    AUTO( uxe, 1. - exp( lambda * Px() ) * cos( 2.*pi*Py() ) );
    AUTO( uye, lambda/( 2.*pi ) * exp( lambda * Px() ) * sin( 2.*pi*Py() ) );
    AUTO( pe, 0.5*( 1.-exp( 2.*lambda*Px() ) ) );
    AUTO( uxedx, -lambda * exp( lambda * Px() ) * cos( 2.*pi*Py() ) );
    AUTO( uyedx, lambda*lambda/( 2.*pi ) * exp( lambda * Px() ) * sin( 2.*pi*Py() ) );
    AUTO( uxedy, 2.*pi* exp( lambda * Px() ) * sin( 2.*pi*Py() ) );
    AUTO( uyedy, lambda * exp( lambda * Px() ) * cos( 2.*pi*Py() ) );
    mesh_ptr_type mesh = U.functionSpace()->mesh();

    value_type uErrorL2 =
        std::sqrt( integrate( elements( *mesh ), M_im,
                              vf::pow( uxe*oneX()+uye*oneY()-idv( U ) , 2.0 )
                            ).evaluate()( 0,0 ) );

    if ( M_uErrorL2 < 0.0 )
        M_uErrorL2 = uErrorL2;

    value_type uErrorH1 =
        std::sqrt( uErrorL2 * uErrorL2 +
                   integrate( elements( *mesh ), M_im,
                              vf::pow( ( uxedx*oneX()+uxedy*oneY() )-trans( gradv( U.comp( X ) ) ), 2.0 ) +
                              vf::pow( ( uyedx*oneX()+uyedy*oneY() )-trans( gradv( U.comp( Y ) ) ), 2.0 )
                            ).evaluate()( 0,0 ) );

    if ( M_uErrorH1 < 0.0 )
        M_uErrorH1 = uErrorH1;

    double int2 = integrate( elements( *mesh ), M_im, vf::pow( pe - idv( p ), 2.0 ) ).evaluate()( 0,0 );
    double int1 = integrate( elements( *mesh ), M_im, pe - idv( p ) ).evaluate()( 0,0 );
    double int0 = integrate( elements( *mesh ), M_im, constant( 1.0 ) ).evaluate()( 0,0 );
    double pErrorL2 = std::sqrt( int2 - int1*int1/int0 );

    if ( M_pErrorL2 < 0.0 )
        M_pErrorL2 = pErrorL2;

    double divError = std::sqrt( integrate( elements( *mesh ), M_im,
                                            vf::pow( divv( U ), 2.0 )
                                          ).evaluate()( 0,0 ) );

    if ( M_divError < 0.0 )
        M_divError = divError;

    VLOG(1) << "[Kovasznay] ||u-u_ex||_L2 = " << uErrorL2 << "\n";
    VLOG(1) << "[Kovasznay] ||u-u_ex||_H1 = " << uErrorH1 << "\n";
    VLOG(1) << "[Kovasznay] ||p-p_ex||_L2 = " << pErrorL2 << "\n";
    VLOG(1) << "[Kovasznay] ||div u ||_L2 = " << divError << "\n";

    VLOG(1) << "[Kovasznay] ||u-u_ex||_L2 / ||u_0-u_ex||_L2 = "
            << uErrorL2/M_uErrorL2 << "\n";
    VLOG(1) << "[Kovasznay] ||u-u_ex||_H1 / ||u_0-u_ex||_H1 = "
            << uErrorH1/M_uErrorH1 << "\n";
    VLOG(1) << "[Kovasznay] ||p-p_ex||_L2 / ||p_0-p_ex||_L2 = "
            << pErrorL2/M_pErrorL2 << "\n";
    VLOG(1) << "[Kovasznay] ||div u ||_L2 / ||div u_0 ||_L2 = "
            << divError/M_divError << "\n";

    M_timers["doexport"].first.restart();

    // -- EXPORT --
    if ( this->vm()["doexport"].as<int>() > 0 &&
            iter % this->vm()["doexport"].as<int>() == 0 )
    {
        element_U_type UErr( U );
        element_p_type pErr ( p  );
        double dpm = -int1/int0;
        UErr = vf::project( U.functionSpace(), elements( *mesh ), idv( U ) - uxe*oneX()-uye*oneY() );
        pErr  = vf::project( p.functionSpace(),  elements( *mesh ), idv( p )  - pe  - dpm );

        M_exporter->step( iter )->setMesh( mesh );
        M_exporter->step( iter )->add( "U", U );
        M_exporter->step( iter )->add( "p", p );
        M_exporter->step( iter )->add( "UErr", UErr );
        M_exporter->step( iter )->add( "pErr", pErr );
        M_exporter->save();
    } // export

    M_timers["doexport"].second += M_timers["doexport"].first.elapsed();
    VLOG(1) << "[Kovasznay] exporting time: "
            << M_timers["doexport"].first.elapsed()
            << "\n";
} // Kovasznay::exportResults

} // Feel
