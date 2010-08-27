/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-03-12

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
   \file convection_lo.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-03-12
 */
#include "convection.hpp"

// variational formulation language
#include <feel/feelvf/vf.hpp>

template <int Order_s, int Order_p, int Order_t>
void Convection<Order_s,Order_p,Order_t> ::initLinearOperator( sparse_matrix_ptrtype& L )
{
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

    // Fluid
    // diffusion
    form2( Xh, Xh, L, _init=true ) =integrate(elements(mesh),_Q<2*Order_s-2>(), cst_ref(sqgr)*trace(gradt(u)*trans(grad(v))));
    Log() << "[initLinearOperator] Fluid Diffusion terms done\n";
    // pressure-velocity terms
    form2( Xh, Xh, L )  += integrate ( elements(mesh), _Q<Order_s+Order_p>(), - idt(p) * div(v) );
    form2( Xh, Xh, L )  += integrate ( elements(mesh), _Q<Order_s+Order_p>(),   divt(u) * id(q) );

    Log() << "[initLinearOperator] Fluid Pressure-Velocity terms done\n";
    // multipliers for zero-mean pressure
    form2( Xh, Xh, L )  += integrate ( elements(mesh), _Q<Order_p>(), id(q)*idt(xi) );
    form2( Xh, Xh, L )  += integrate ( elements(mesh), _Q<Order_p>(), idt(p)*id(eta) );
    form2( Xh, Xh, L )  += integrate ( elements(mesh), _Q<0>(), 0*idt(xi)*id(eta) );
    Log() << "[initLinearOperator] Fluid Pressure-Multipliers terms done\n";

    // weak Dirichlet condition at the walls (u=0)
    AUTO( SigmaNt, (-idt(p)*N()+cst_ref(sqgr)*gradt(u)*N()) );
    AUTO( SigmaN, (-id(q)*N()+cst_ref(sqgr)*grad(v)*N()) );
    form2( Xh, Xh, L )  += integrate ( boundaryfaces(mesh),_Q<2*Order_s-1>(), -trans(SigmaNt)*id(v) );
    form2( Xh, Xh, L )  += integrate ( boundaryfaces(mesh),_Q<2*Order_s-1>(), -trans(SigmaN)*idt(u) );
    form2( Xh, Xh, L )  += integrate ( boundaryfaces(mesh), _Q<2*Order_s>(), +gamma*trans(idt(u))*id(v)/hFace());
    Log() << "[initLinearOperator] Fluid Dirichlet weak BC terms done\n";

    Log() << "[initLinearOperator] done in " << ti.elapsed() << "s\n";

}

// instantiation
template class Convection<2,1,2>;

