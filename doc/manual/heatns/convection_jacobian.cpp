/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-03-04

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
   \file convection_jacobian.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-03-04
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
    this->initLinearOperator2(J);


    this->updateJacobian1( X, J );
    this->updateJacobian2( X, J );

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

// instantiation
template class Convection<2,1,2>;
