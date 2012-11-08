/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-03-05

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
   \file convection_jacobian1.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-03-05
 */
#include "convection.hpp"

// variational formulation language
#include <feel/feelvf/vf.hpp>

// <int Order_s, int Order_p, int Order_t>
void Convection ::updateJacobian1( const vector_ptrtype& X,
        sparse_matrix_ptrtype& D )
{
    auto mesh = Xh->mesh();

    auto U = Xh->element( "u" );
    U = *X;
    auto V = Xh->element( "v" );
    auto W = Xh->element( "v" );
    auto u = U. element<0>(); // fonction vitesse
    auto v = V. element<0>(); // fonction test vitesse
    auto p = U. element<1>(); // fonction pression
    auto q = V. element<1>(); // fonction test pression
    auto t = U. element<2>(); // fonction temperature
    auto s = V. element<2>(); // fonction test temperature
#if defined( FEELPP_USE_LM )
    auto xi = U. element<3>(); // fonction multipliers
    auto eta = V. element<3>(); // fonction test multipliers
#endif

    LOG(INFO) << "[updateJacobian1] ||U|| = " << U.l2Norm() << "\n";
    LOG(INFO) << "[updateJacobian1] ||u|| = " << u.l2Norm() << "\n";
    LOG(INFO) << "[updateJacobian1] ||p|| = " << p.l2Norm() << "\n";
    LOG(INFO) << "[updateJacobian1] ||t|| = " << t.l2Norm() << "\n";
#if defined( FEELPP_USE_LM )
    LOG(INFO) << "[updateJacobian1] ||xi|| = " << xi.l2Norm() << "\n";
#endif

    double gr= M_current_Grashofs;
    double sqgr( 1/math::sqrt( gr ) );
    double pr=M_current_Prandtl;
    double sqgrpr( 1/( pr*math::sqrt( gr ) ) );
    double gamma( this->vm()["penalbc"]. as<double>() );
    int steady( this->vm()["steady"]. as<int>() );
    double dt=this->vm()["dt"]. as<double>();

    int adim=this->vm()["adim"]. as<int>();

    double k=this->vm()["k"]. as<double>();
    double nu=this->vm()["nu"]. as<double>();
    double rho=this->vm()["rho"]. as<double>();

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


    // Fluid-NS
    // fluid convection derivatives: attention 2 terms

    form2( _test=Xh,_trial=Xh, _matrix=D )  += integrate ( _range=elements( mesh ),_expr=cst( a )*trans( id( v ) )*( gradv( u ) )*idt( u ) );
    form2( _test=Xh,_trial=Xh, _matrix=D )  +=integrate ( _range=elements( mesh ), _expr=cst( a )*trans( id( v ) )*( gradt( u )*idv( u ) ) );

    LOG(INFO) << "[updateJacobian1] Convection terms done\n";


}


// instantiation
// class Convection<2,1,2>;
