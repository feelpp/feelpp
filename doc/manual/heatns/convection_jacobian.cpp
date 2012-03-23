/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

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
   \file convection_jacobian.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2012-03-22
 */
#include "convection.hpp"

template <int Order_s, int Order_p, int Order_t>
void Convection<Order_s,Order_p,Order_t> ::updateJacobian( const vector_ptrtype& X,
                                                           sparse_matrix_ptrtype& J)
{
    boost::timer ti;
    Log() << "[updateJacobian] start\n";

    //D->zero();
    //L->zero();
    this->initLinearOperator(J);

    this->updateJacobian1( X, J );

    //L->close();

    //J->mat() = MatDuplicate( L->mat() );
    //*J = L;
    //J = L;
    J->close();
    //L->printMatlab( "L.m" );
    //J->printMatlab( "J1.m" );
    //D->printMatlab( "D1.m" );
    //J->addMatrix( 1.0, D );
    //J->printMatlab( "J.m" );
    Log() << "[updateJacobian] done in " << ti.elapsed() << "s\n";

}

template <int Order_s, int Order_p, int Order_t>
void Convection<Order_s,Order_p,Order_t> ::updateJacobian1( const vector_ptrtype& X,
                                                            sparse_matrix_ptrtype& D)
{
#if 1
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

    Log() << "[updateJacobian1] ||U|| = " << U.l2Norm() << "\n";
    Log() << "[updateJacobian1] ||u|| = " << u.l2Norm() << "\n";
    Log() << "[updateJacobian1] ||p|| = " << p.l2Norm() << "\n";
    Log() << "[updateJacobian1] ||t|| = " << t.l2Norm() << "\n";
    Log() << "[updateJacobian1] ||xi|| = " << xi.l2Norm() << "\n";

    double gr= M_current_Grashofs;
    double sqgr(1/math::sqrt(gr));
    double pr=M_current_Prandtl;
    double sqgrpr(1/(pr*math::sqrt(gr)));
    double gamma(this->vm()["penalbc"].template as<double>());

    // Fluid-NS
    // fluid convection derivatives: attention 2 terms
    form2( Xh,Xh, D )  +=
        integrate (elements(mesh),
                   trans(id(v))*(gradv(u))*idt(u));
    form2( Xh,Xh, D )  +=
        integrate (elements(mesh),
                   trans(id(v))*(gradt(u)*idv(u)) );
    Log() << "[updateJacobian1] Convection terms done\n";


    //
    // temperature derivatives
    //
    // heat convection by the fluid: attention 2 terms
    form2( Xh,Xh, D ) +=
        integrate ( elements(mesh),
                    grad(s)*(idv(t)*idt(u)) );
    form2( Xh,Xh, D ) +=
        integrate ( elements(mesh),
                    grad(s)*(idt(t)*idv(u)) );
    form2( Xh,Xh, D ) +=
        integrate ( boundaryfaces(mesh),
                    (trans(idv(u))*N())*id(s)*idt(t) );
    form2( Xh,Xh, D ) +=
        integrate ( boundaryfaces(mesh),
                    (trans(idt(u))*N())*id(s)*idv(t) );
    Log() << "[updateJacobian] Temperature convection terms done\n";
#endif
}


template class Convection<2,1,2>;
