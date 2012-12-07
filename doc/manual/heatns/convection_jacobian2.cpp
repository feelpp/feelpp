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
   \file convection_jacobian2.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-03-05
 */
#include "convection.hpp"

// variational formulation language
#include <feel/feelvf/vf.hpp>

// <int Order_s, int Order_p, int Order_t>
void Convection::updateJacobian2( const vector_ptrtype& X,
        sparse_matrix_ptrtype& D )
{
    mesh_ptrtype mesh = Xh->mesh();
    element_type U( Xh, "u" );
    U = *X;
    LOG(INFO) << "[updateJacobian] ||X|| = " << U.l2Norm() << "\n";
    element_type V( Xh, "v" );
    element_type W( Xh, "v" );
    element_0_type u = U. element<0>(); // fonction vitesse
    element_0_type v = V. element<0>(); // fonction test vitesse
    element_1_type p = U. element<1>(); // fonction pression
    element_1_type q = V. element<1>(); // fonction test pression
    element_2_type t = U. element<2>(); // fonction temperature
    element_2_type s = V. element<2>(); // fonction test temperature
#if defined( FEELPP_USE_LM )
    element_3_type xi = U. element<3>(); // fonction multipliers
    element_3_type eta = V. element<3>(); // fonction test multipliers
#endif

    double gr= M_current_Grashofs;
    double sqgr( 1/math::sqrt( gr ) );
    double pr = M_current_Prandtl;
    double sqgrpr=( 1/( pr*math::sqrt( gr ) ) );
    double gamma( this->vm()["penalbc"]. as<double>() );
    double dt=this->vm()["dt"]. as<double>();
    int adim=this->vm()["adim"]. as<int>();
    double pC=1;

    if ( adim == 0 ) pC = this->vm()["pC"]. as<double>();

    //
    // temperature derivatives
    //
    // heat convection by the fluid: attention 2 terms
    form2( _test=Xh, _trial=Xh, _matrix=D ) +=
        integrate ( elements( mesh ),
                    pC*grad( s )*( idv( t )*idt( u ) ) );

    form2( _test=Xh, _trial=Xh, _matrix=D ) +=
        integrate ( elements( mesh ),
                    pC*grad( s )*( idt( t )*idv( u ) ) );

    form2( _test=Xh, _trial=Xh, _matrix=D ) +=
        integrate ( boundaryfaces( mesh ),
                    pC*( trans( idv( u ) )*N() )*id( s )*idt( t ) );

    form2( _test=Xh, _trial=Xh, _matrix=D ) +=
        integrate ( boundaryfaces( mesh ),
                    pC*( trans( idt( u ) )*N() )*id( s )*idv( t ) );

    LOG(INFO) << "[updateJacobian2] Temperature convection terms done\n";


    LOG(INFO) << "[updateJacobian2] Temperature weak Dirichlet BC terms done\n";


}

// instantiation
// class Convection<2,1,2>;
