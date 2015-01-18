/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-03-04

  Copyright (C) 2009 Universite Joseph Fourier (Grenoble I)

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
   \file convection_residual.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-03-04
 */
#include "convection.hpp"
#include"functionSup.cpp"
// variational formulation language
#include <feel/feelvf/vf.hpp>

//<int Order_s, int Order_p, int Order_t>
void
Convection::updateResidual( const vector_ptrtype& X, vector_ptrtype& R )
{

    boost::timer ti;
    LOG(INFO) << "[updateResidual] start\n";

    mesh_ptrtype mesh = Xh->mesh();

    element_type U( Xh, "u" );
    U = *X;
    element_type Un( Xh, "un" );
    element_type V( Xh, "v" );
    element_0_type u = U. element<0>(); // fonction vitesse
    element_0_type un = Un. element<0>(); // fonction vitesse
    element_0_type v = V. element<0>(); // fonction test vitesse
    element_1_type p = U. element<1>(); // fonction pression
    element_1_type pn = Un. element<1>(); // fonction pression
    element_1_type q = V. element<1>(); // fonction test pression
    element_2_type t = U. element<2>(); // fonction temperature
    element_2_type tn = Un. element<2>(); // fonction temperature
    element_2_type s = V. element<2>(); // fonction test temperature
#if defined( FEELPP_USE_LM )
    element_3_type xi = U. element<3>(); // fonction multipliers
    element_3_type eta = V. element<3>(); // fonction test multipliers
#endif

    double gr( M_current_Grashofs );
    double sqgr( 1/math::sqrt( gr ) );
    double pr = M_current_Prandtl;
    double sqgrpr( 1/( pr*math::sqrt( gr ) ) );
    double gamma( this->vm()["penalbc"]. as<double>() );
    double k=this->vm()["k"]. as<double>();
    double nu=this->vm()["nu"]. as<double>();
    double rho=this->vm()["rho"]. as<double>();
    double dt=this->vm()["dt"]. as<double>();
    int adim=this->vm()["adim"]. as<int>();

    double T0=this->vm()["T0"]. as<double>();
    double neum=this->vm()["neum"]. as<double>();
    double pC=1;

    if ( adim == 0 ) pC = this->vm()["pC"]. as<double>();

    LOG(INFO) << "residual: residual: gr = " << gr << "\n";
    LOG(INFO) << "residual: pr = " << pr << "\n";
    LOG(INFO) << "residual: sqgr = " << sqgr << "\n";
    LOG(INFO) << "residual: sqgrpr = " << sqgrpr << "\n";
    LOG(INFO) << "residual: gamma = " << gamma << "\n";
    int weakdir=this->vm()["weakdir"]. as<int>();



    //choix de la valeur des paramètres dimensionnés ou adimensionnés
    double a=0.0,b=0.0,c=0.0;

    if ( adim==1 )
    {
        a=1;
        b=sqgr;
        c=sqgrpr;
    }

    else
    {
        a=rho;
        b=nu;
        c=k;
    }

    double expansion = 1;

    if ( adim == 0 ) expansion=3.7e-3;


    // -- Partie Navier-Stokes --

    // gravity

    form1( Xh, _vector=R ) =
        integrate ( elements( mesh ),
                    // convection
                    cst( a )*trans( gradv( u )*idv( u ) )*id( v ) );


    form1( Xh, _vector=R ) +=
        integrate ( elements( mesh ),
                    // heat diffusion
                    cst( b ) * trace( gradv( u ) * trans( grad( v ) ) ) );
    form1( Xh, _vector=R ) +=
        integrate ( elements( mesh ),
                    // pressure-velocity terms
                    +divv( u ) * id( q ) - idv( p ) * div( v )  );

#if defined( FEELPP_USE_LM )
    // only when using pressure is up to a constant
    form1( Xh, _vector=R ) +=
        integrate ( elements( mesh ),
                    // multipliers for zero-mean pressure
                    +id( q )*idv( xi )+idv( p )*id( eta ) );
#endif

    auto SigmaNv = ( -idv( p )*N()+cst( b )*gradv( u )*N() );
    auto SigmaN = ( -id( q )*N()+cst( b )*grad( v )*N() );

#if CONVECTION_DIM == 2
    // right hand side
    form1( Xh, _vector=R ) +=
        integrate ( elements( mesh ),
                    // buyoancy force
                    -expansion*trans( vec( constant( 0. ),idv( t ) ) )*id( v ) );
#else
    // right hand side
    form1( Xh, _vector=R ) +=
        integrate ( elements( mesh ),
                    // buyoancy force
                    -expansion*trans( vec( cst(0.), constant( 0. ),idv( t ) ) )*id( v ) );
#endif

    // -- Partie Chaleur --
    form1( Xh, _vector=R ) +=
        integrate ( elements( mesh ),
                    // heat convection by the fluid
                    pC*grad( s )*( idv( t ) * idv( u ) ) );
    form1( Xh, _vector=R ) +=
        integrate ( boundaryfaces( mesh ),
                    pC*( trans( idv( u ) )*N() )*id( s )*idv( t ) );



    form1( Xh, _vector=R ) +=
        integrate ( elements( mesh ),
                    // heat diffusion
                    cst( c ) * gradv( t ) * trans( grad( s ) ) );

    if ( adim==1 )
        neum=1;
    form1( Xh, _vector=R ) +=
        integrate ( markedfaces( mesh, "Tflux"),
                    // heat flux on the right side
                    - id( s )*cst( c )*cst( neum )  );

    if ( weakdir == 1 )
    {
        auto SigmaNv = ( -idv( p )*N()+cst_ref( sqgr )*gradv( u )*N() );
        auto SigmaN = ( -id( q )*N()+cst_ref( sqgr )*grad( v )*N() );
        form1( Xh, _vector=R ) +=
            integrate ( marked2faces( mesh, "F.wall" ),
                        // weak Dirichlet condition at the walls (u=0)
                        -trans( SigmaNv )*id( v )
                        -trans( SigmaN )*idv( u )
                        +gamma*trans( idv( u ) )*id( v )/hFace() );

        form1( Xh, _vector=R ) +=
            integrate ( markedfaces( mesh, "Tfixed" ),
                        // weak dirichlet condition T=T_0 | left side
                        -gradv( t )*N()*id( s )*cst_ref( sqgrpr )
                        -grad( s )*N()*idv( t )*cst_ref( sqgrpr ) );
        form1( Xh, _vector=R ) 	+=
            integrate ( markedfaces( mesh, "Tfixed"  ),
                        cst_ref( gamma )*idv( t )*id( s )/hFace() );

    }
    R->close();

    if ( weakdir == 0 )
    {
        if ( adim==1 )
        {
            LOG(INFO) << "dirichlet marker Tfixed : " <<  t.start() << "  "  << mesh->markerName("Tfixed") << "\n";
            modifVec( markedfaces( mesh,"Tfixed" ),t,R,cst( -1.0 ) );
        }
        else
        {
            modifVec( markedfaces( mesh,"Tfixed" ),t,R,cst( T0 ) );
        }
    }
    LOG(INFO) << "dirichlet boundary faces(velocity) : " <<  u.start()  << "\n";
    modifVec( boundaryfaces( mesh ),u,R,one()*0 );


    LOG(INFO) << "[updateResidual] done in " << ti.elapsed() << "s\n";

}


// instantiation
// class Convection<2,1,2>;
