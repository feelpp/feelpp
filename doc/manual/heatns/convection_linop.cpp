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
   \file convection_linop.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2012-03-22
 */
#include "convection.hpp"




template <int Order_s, int Order_p, int Order_t>
void Convection<Order_s,Order_p,Order_t> ::initLinearOperator( sparse_matrix_ptrtype& L )
{
#if 1
    boost::timer ti;
    Log() << "[initLinearOperator] start\n";

    mesh_ptrtype mesh = Xh->mesh();
    element_type U( Xh, "u" );
    element_type V( Xh, "v" );
	element_0_type u = U.template element<0>(); // fonction vitesse
    element_0_type v = V.template element<0>(); // fonction test vitesse
	element_1_type p = U.template element<1>(); // fonction pression
	element_1_type q = V.template element<1>(); // fonction test pression
	element_2_type t = U.template element<2>(); // fonction temperature
	element_2_type s = V.template element<2>(); // fonction test temperature
    element_3_type xi = U.template element<3>(); // fonction multipliers
	element_3_type eta = V.template element<3>(); // fonction test multipliers

    //M_oplin = oplin_ptrtype( new oplin_type( Xh, Xh, M_backend ) );
    //M_oplin->mat().zero();

    double gr= M_current_Grashofs;
    double sqgr(1/math::sqrt(gr));
    double pr = M_current_Prandtl;
    double sqgrpr(1/(pr*math::sqrt(gr)));
    double gamma(this->vm()["penalbc"].template as<double>());

    auto bf = form2( _test=Xh, _trial=Xh, _matrix=L );
    // Fluid
    // diffusion
    bf =integrate(elements(mesh), cst_ref(sqgr)*trace(gradt(u)*trans(grad(v))));
    Log() << "[initLinearOperator] Fluid Diffusion terms done\n";
    // pressure-velocity terms
    bf  += integrate ( elements(mesh),  - idt(p) * div(v) );
    bf  += integrate ( elements(mesh),    divt(u) * id(q) );

    Log() << "[initLinearOperator] Fluid Pressure-Velocity terms done\n";
    // multipliers for zero-mean pressure
    bf  += integrate ( elements(mesh),  id(q)*idt(xi) );
    bf  += integrate ( elements(mesh),  idt(p)*id(eta) );
    Log() << "[initLinearOperator] Fluid Pressure-Multipliers terms done\n";

    // weak Dirichlet condition at the walls (u=0)
    auto SigmaNt = (-idt(p)*N()+cst_ref(sqgr)*gradt(u)*N());
    auto SigmaN = (-id(q)*N()+cst_ref(sqgr)*grad(v)*N());
    bf  += integrate ( boundaryfaces(mesh),-trans(SigmaNt)*id(v) );
    bf  += integrate ( boundaryfaces(mesh),-trans(SigmaN)*idt(u) );
    bf  += integrate ( boundaryfaces(mesh),+gamma*trans(idt(u))*id(v)/hFace());
    Log() << "[initLinearOperator] Fluid Dirichlet weak BC terms done\n";

    Log() << "[initLinearOperator] done in " << ti.elapsed() << "s\n";


    // Temperature
    // buyoancy forces c(theta,v)
    bf +=integrate(elements(mesh),-idt(t)*(trans(vec(constant(0.),constant(1.0)))*id(v)));

    Log() << "[initLinearOperator] temperature Force terms done\n";
    // heat conduction/diffusion: e(beta1,theta,chi)+f(theta,chi)
    bf  += integrate(elements(mesh), cst_ref(sqgrpr)*gradt(t)*trans(grad(s)));
    Log() << "[initLinearOperator] Temperature Diffusion terms done\n";

    // weak Dirichlet on temperature (T=0|left wall)
    bf  += integrate ( markedfaces(mesh,mesh->markerName( "Tfixed" )),
                            - gradt(t)*N()*id(s)*cst_ref(sqgrpr) );
    bf  += integrate ( markedfaces(mesh,mesh->markerName( "Tfixed" )),
                            - grad(s)*N()*idt(t)*cst_ref(sqgrpr) );
    bf  += integrate ( markedfaces(mesh,mesh->markerName( "Tfixed" )),
                            gamma*idt(t)*id(s)/hFace());
    Log() << "[initLinearOperator2] done in " << ti.elapsed() << "s\n";
#endif
}

template class Convection<2,1,2>;
