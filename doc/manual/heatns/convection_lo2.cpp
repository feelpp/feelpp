/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-03-12

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
   \file convection_lo.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-03-12
 */
#include "convection.hpp"

// variational formulation language
#include <feel/feelvf/vf.hpp>

// <int Order_s, int Order_p, int Order_t>
void Convection ::initLinearOperator2( sparse_matrix_ptrtype& L )
{
    boost::timer ti;
    LOG(INFO) << "[initLinearOperator2] start\n";

    mesh_ptrtype mesh = Xh->mesh();
    element_type U( Xh, "u" );
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

    double gr= M_current_Grashofs;
    double sqgr( 1/math::sqrt( gr ) );
    double pr = M_current_Prandtl;
    double sqgrpr( 1/( pr*math::sqrt( gr ) ) );
    double gamma( this->vm()["penalbc"]. as<double>() );

    double k=this->vm()["k"]. as<double>();
    double nu=this->vm()["nu"]. as<double>();
    double rho=this->vm()["rho"]. as<double>();
    //double dt=this->vm()["dt"]. as<double>();
    int adim=this->vm()["adim"]. as<int>();
    int weakdir=this->vm()["weakdir"]. as<int>();
    //choix de la valeur des paramètres dimensionnés ou adimensionnés
    double a=0.0,b=0.0,c=0.0;
    double pC=1;

    if ( adim == 0 ) pC = this->vm()["pC"]. as<double>();

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

    auto bf = form2( _test=Xh, _trial=Xh, _matrix=L );
    // Temperature
#if CONVECTION_DIM==2
    // buyoancy forces c(theta,v)
    bf +=integrate( _range=elements( mesh ),
                    _expr=-expansion*idt( t )*( trans( vec( constant( 0. ),constant( 1.0 ) ) )*id( v ) ) );
#else
    bf +=integrate( _range=elements( mesh ),
                    _expr=-expansion*idt( t )*( trans( vec( cst(0.), constant( 0. ),constant( 1.0 ) ) )*id( v ) ) );
#endif

    LOG(INFO) << "[initLinearOperator] temperature Force terms done\n";
    // heat conduction/diffusion: e(beta1,theta,chi)+f(theta,chi)
    bf  += integrate( _range=elements( mesh ),
                      _expr=cst( c )*gradt( t )*trans( grad( s ) ) );
    LOG(INFO) << "[initLinearOperator] Temperature Diffusion terms done\n";

    if ( weakdir == 1 )
    {
        // weak Dirichlet on temperature (T=0|left wall)
        bf  += integrate ( markedfaces( mesh, "Tfixed"  ),
                           - gradt( t )*N()*id( s )*cst_ref( sqgrpr ) );
        bf  += integrate ( markedfaces( mesh, "Tfixed"  ),
                           - grad( s )*N()*idt( t )*cst_ref( sqgrpr ) );
        bf  += integrate ( markedfaces( mesh, "Tfixed"  ),
                           gamma*idt( t )*id( s )/hFace() );
    }

    LOG(INFO) << "[initLinearOperator2] done in " << ti.elapsed() << "s\n";
}

// instantiation
// class Convection<2,1,2>;
