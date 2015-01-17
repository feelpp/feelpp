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
   \file convection_jacobian.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-03-04
 */
#include <convection.hpp>

// <int Order_s, int Order_p, int Order_t>
void Convection ::updateJacobian( const vector_ptrtype& X,
        sparse_matrix_ptrtype& J )
{
    boost::timer ti;
    LOG(INFO) << "[updateJacobian] start\n";

    if ( !J )
    {
        J =  M_backend->newMatrix( _test=Xh, _trial=Xh );
        M_D =  M_backend->newMatrix( _test=Xh, _trial=Xh );
        M_L =  M_backend->newMatrix( _test=Xh, _trial=Xh );

        this->initLinearOperator( M_L );
        this->initLinearOperator2( M_L );

    }
    else
    {
        M_D->zero();
        J->zero();
    }

    this->updateJacobian1( X, M_D );
    this->updateJacobian2( X, M_D );

    M_D->close();
    M_L->close();
    J->close();
    
    //conditions fortes de dir
    auto mesh = Xh->mesh();
    auto U =  Xh->element( "u" );
    auto u = U. element<0>(); // fonction vitesse
    auto t= U. element<2>(); // fonction temperature
    auto Rtemp =  M_backend->newVector( Xh );
    int weakdir( this->vm()["weakdir"]. as<int>() );
    int adim = this->vm()["adim"]. as<int>();
    double T0 = this->vm()["T0"]. as<double>();

    if ( weakdir == 0 )
    {
        //vitesse
        form2( Xh, Xh, _matrix=M_D )  += on( boundaryfaces( mesh ),u, Rtemp,one()*0. );

        if ( adim==1 )
            //temperature
            form2( Xh, Xh, _matrix=M_D )  += on ( markedfaces( mesh, "Tfixed" ),t,Rtemp,cst( 0.0 ) );

        else
            form2( Xh, Xh, _matrix=M_D )  += on ( markedfaces( mesh, "Tfixed" ),t,Rtemp,cst( T0 ) );
    }
    LOG(INFO) << "[updateJacobian] done in " << ti.elapsed() << "s\n";

    J->addMatrix( 1.0, M_L);
    J->addMatrix( 1.0, M_D);
    
}

// instantiation
//emplate class Convection<2,1,2>;
