/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2012-03-22

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2012-03-22
 */
#include "convection.hpp"


template<int Order_s, int Order_p, int Order_t>
void
Convection<Order_s,Order_p,Order_t>::updateResidual( const vector_ptrtype& X, vector_ptrtype& R )
{
    boost::timer ti;
    Log() << "[updateResidual] start\n";

    mesh_ptrtype mesh = Xh->mesh();

    element_type U( Xh, "u" );
    U = *X;
    element_type V( Xh, "v" );
    element_type W( Xh, "v" );
    element_0_type u = U.template element<0>(); // fonction vitesse
    element_0_type v = V.template element<0>(); // fonction test vitesse
    element_1_type p = U.template element<1>(); // fonction pression
    element_1_type q = V.template element<1>(); // fonction test pression
    element_2_type t = U.template element<2>(); // fonction temperature
    element_2_type s = V.template element<2>(); // fonction test temperature
    element_3_type xi = U.template element<3>(); // fonction multipliers
    element_3_type eta = V.template element<3>(); // fonction test multipliers

    double gr( M_current_Grashofs );
    double sqgr( 1/math::sqrt( gr ) );
    double pr = M_current_Prandtl;
    double sqgrpr( 1/( pr*math::sqrt( gr ) ) );
    double gamma( this->vm()["penalbc"].template as<double>() );

    Log() << "gr = " << gr << "\n";
    Log() << "pr = " << pr << "\n";
    Log() << "sqgr = " << sqgr << "\n";
    Log() << "sqgrpr = " << sqgrpr << "\n";
    Log() << "gamma = " << gamma << "\n";

    form1( Xh, R, _init=true ) =
        integrate ( elements( mesh ),
                    // convection
                    1*trans( gradv( u )*idv( u ) )*id( v ) );


    form1( Xh, R ) +=
        integrate ( elements( mesh ),
                    // heat diffusion
                    cst_ref( sqgr ) * trace( gradv( u ) * trans( grad( v ) ) ) );
    form1( Xh, R ) +=
        integrate ( elements( mesh ),
                    // pressure-velocity terms
                    +divv( u ) * id( q ) - idv( p ) * div( v ) );
    form1( Xh, R ) +=
        integrate ( elements( mesh ),
                    // multipliers for zero-mean pressure
                    +id( q )*idv( xi )+idv( p )*id( eta )  );

    auto SigmaNv = ( -idv( p )*N()+cst_ref( sqgr )*gradv( u )*N() );
    auto SigmaN = ( -id( q )*N()+cst_ref( sqgr )*grad( v )*N() );
    form1( Xh, R ) +=
        integrate ( boundaryfaces( mesh ),
                    // weak Dirichlet condition at the walls (u=0)
                    -trans( SigmaNv )*id( v )
                    -trans( SigmaN )*idv( u )
                    +gamma*trans( idv( u ) )*id( v )/hFace() );

    // right hand side
    form1( Xh, R ) +=
        integrate ( elements( mesh ),
                    // buyoancy force
                    -trans( vec( constant( 0. ),idv( t ) ) )*id( v ) );
    //-trans(idv(t)*oneY())*id(v) );

    // -- Partie Chaleur --
    form1( Xh, R ) +=
        integrate ( elements( mesh ),
                    // heat convection by the fluid
                    grad( s )*( idv( t ) * idv( u ) ) );
    form1( Xh, R ) +=
        integrate ( boundaryfaces( mesh ),
                    ( trans( idv( u ) )*N() )*id( s )*idv( t ) );

    form1( Xh, R ) +=
        integrate ( boundaryfaces( mesh ),
                    ( trans( idv( u ) )*N() )*id( s )*idv( t ) );

    form1( Xh, R ) +=
        integrate ( elements( mesh ),
                    // heat diffusion
                    cst_ref( sqgrpr ) * gradv( t ) * trans( grad( s ) ) );

    form1( Xh, R ) +=
        integrate ( markedfaces( mesh,mesh->markerName( "Tflux" ) ),
                    // heat flux on the right side
                    - id( s )*cst_ref( sqgrpr ) );

    form1( Xh, R ) +=
        integrate ( markedfaces( mesh,mesh->markerName( "Tfixed" ) ),
                    // weak dirichlet condition T=T_0 | left side
                    -gradv( t )*N()*id( s )*cst_ref( sqgrpr )
                    -grad( s )*N()*idv( t )*cst_ref( sqgrpr ) );
    form1( Xh, R ) +=
        integrate ( boundaryfaces( mesh ),
                    ( trans( idv( u ) )*N() )*idv( t )*id( s ) );

    // -- weak Dirichlet conditions : 0 for temperature and velocity
    form1( Xh,R ) 	+=
        integrate ( markedfaces( mesh,mesh->markerName( "Tfixed" ) ),
                    cst_ref( gamma )*idv( t )*id( s )/hFace() );

    R->close();
    Log() << "[updateResidual] done in " << ti.elapsed() << "s\n";

}

template class Convection<2,1,2>;
