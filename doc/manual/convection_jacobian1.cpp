/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-03-05

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-03-05
 */
#include "convection.hpp"

// variational formulation language
#include <feel/feelvf/vf.hpp>

template <int Order_s, int Order_p, int Order_t>
void Convection<Order_s,Order_p,Order_t> ::updateJacobian1( const vector_ptrtype& X,
                                                            sparse_matrix_ptrtype& D)
{
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
                   _Q<3*Order_s-1>(),
                   trans(id(v))*(gradv(u))*idt(u));
    form2( Xh,Xh, D )  +=
        integrate (elements(mesh),
                   _Q<3*Order_s-1>(),
                   trans(id(v))*(gradt(u)*idv(u)) );
    Log() << "[updateJacobian1] Convection terms done\n";


}


// instantiation
template class Convection<2,1,2>;
