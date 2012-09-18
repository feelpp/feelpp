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
   \file convection_other.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Elisa Schenone
   \date 2012-08-13
 */
#include <boost/lexical_cast.hpp>

#include "convection_crb.hpp"

typedef Eigen::VectorXd theta_vector_type;

void Convection_crb::init()
{
    
    mesh_ptrtype mesh;
    mesh = createGMSHMesh( _mesh=new mesh_type,
                          _desc=createMesh(),
                          _update=MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
    Xh = space_type::New( mesh );

    pT = element_ptrtype( new element_type( Xh ) );

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
    
    using namespace Feel::vf;
    Feel::ParameterSpace<2>::Element mu_min( M_Dmu );
    double Grmin( this->vm()["Grmin"]. as<double>() );
    double Prmin( this->vm()["Prmin"]. as<double>() );
    mu_min << Grmin, Prmin;
    M_Dmu->setMin( mu_min );
    Feel::ParameterSpace<2>::Element mu_max( M_Dmu );
    double Grmax( this->vm()["Grmax"]. as<double>() );
    double Prmax( this->vm()["Prmax"]. as<double>() );
    mu_max << Grmax, Prmax;
    M_Dmu->setMax( mu_max );
    
    
    //  initialisation de A0, A1, A2
    M_Aq.resize( Qa() );
    for (int ii=0; ii<Qa(); ii++)
    {
        M_Aq[ii] = M_backend->newMatrix( Xh, Xh );
    }
        
    M_Fq.resize( Nl() );
    for (int ii=0; ii<Nl(); ii++)
    {
        M_Fq[ii].resize( Ql(ii) );
        for (int jj=0; jj<Ql(ii); jj++)
        {
            M_Fq[ii][jj] = M_backend->newVector( Xh );
        }        
    }
        
    D = M_backend->newMatrix( Xh, Xh );
    F = M_backend->newVector( Xh );
    
    double gamma( this->vm()["penalbc"]. as<double>() );
    double k=this->vm()["k"]. as<double>();
    double nu=this->vm()["nu"]. as<double>();
    double rho=this->vm()["rho"]. as<double>();
    int weakdir=this->vm()["weakdir"]. as<int>();
    
    int adim=this->vm()["adim"]. as<int>();    
    
    // Fluid
    // diffusion
    
    // M_Aq[0] = C = grad(u)*grav(v)
    form2( _test=Xh, _trial=Xh, _matrix=M_Aq[0], _init=true ) = integrate( _range=elements( mesh ),
                                                                          _expr=trace( gradt( u )*trans( grad( v ) ) )  );
    
    // heat diffusion: M_Aq[1] = G = grad(t)*grad(s)
    form2( _test=Xh, _trial=Xh, _matrix=M_Aq[1], _init=true ) = integrate( _range=elements( mesh ),
                                                                          _expr=gradt( t )*trans( grad( s ) ) );

    // pressure-velocity terms
    // M_Aq[2] = -B = -p*div(v)
    form2( _test=Xh, _trial=Xh, _matrix=M_Aq[2], _init=true ) = integrate ( _range=elements( mesh ), _expr=- idt( p ) * div( v ) );
    
    // M_Aq[2] += B^t = q*div(u)
    form2( _test=Xh, _trial=Xh, _matrix=M_Aq[2] ) += integrate ( _range=elements( mesh ), _expr=divt( u ) * id( q ) );
    
    if ( weakdir == 1 )
    {
        // weak Dirichlet condition at the walls (u=0)
        form2( _test=Xh, _trial=Xh, _matrix=M_Aq[0] ) += integrate ( marked2faces( mesh, "F.wall" ), -trans( gradt( u )*N() )*id( v ) );
        form2( _test=Xh, _trial=Xh, _matrix=M_Aq[2] ) += integrate ( marked2faces( mesh, "F.wall" ), -trans( -idt( p )*N() )*id( v ) );
        form2( _test=Xh, _trial=Xh, _matrix=M_Aq[0] )  += integrate ( marked2faces( mesh, "F.wall" ), -trans( grad( v )*N() )*idt( u ) );
        form2( _test=Xh, _trial=Xh, _matrix=M_Aq[2] )  += integrate ( marked2faces( mesh, "F.wall" ), -trans( -id( q )*N() )*idt( u ) );
        form2( _test=Xh, _trial=Xh, _matrix=M_Aq[2] )  += integrate ( marked2faces( mesh, "F.wall" ), +gamma*trans( idt( u ) )*id( v )/hFace() );
    }
    
#if defined( FEELPP_USE_LM )
    // multipliers for zero-mean pressure
    form2( _test=Xh, _trial=Xh, _matrix=M_Aq[2] ) += integrate ( _range=elements( mesh ), _expr=id( q )*idt( xi ) );
    form2( _test=Xh, _trial=Xh, _matrix=M_Aq[2] ) += integrate ( _range=elements( mesh ), _expr=idt( p )*id( eta ) );
#endif
    
    
    double expansion = 1;
    
    if ( adim == 0 ) expansion=3.7e-3;
    
    // Temperature
    
#if CONVECTION_DIM==2
    // buyoancy forces c(theta,v)
    form2( _test=Xh, _trial=Xh, _matrix=M_Aq[2] ) +=integrate( _range=elements( mesh ),
                                                              _expr=-expansion*idt( t )*( trans( vec( constant( 0. ),constant( 1.0 ) ) )*id( v ) ) );
#else
    form2( _test=Xh, _trial=Xh, _matrix=M_Aq[2] ) +=integrate( _range=elements( mesh ),
                                                              _expr=-expansion*idt( t )*( trans( vec( cst(0.), constant( 0. ),constant( 1.0 ) ) )*id( v ) ) );
#endif
    
    
    if ( weakdir == 1 )
    {
        // weak Dirichlet on temperature (T=0|left wall)
        form2( _test=Xh, _trial=Xh, _matrix=M_Aq[1] ) += integrate ( markedfaces( mesh,mesh->markerName( "Tfixed" ) ),
                                                                    - gradt( t )*N()*id( s ) );
        form2( _test=Xh, _trial=Xh, _matrix=M_Aq[1] ) += integrate ( markedfaces( mesh,mesh->markerName( "Tfixed" ) ),
                                                                    - grad( s )*N()*idt( t ) );
        form2( _test=Xh, _trial=Xh, _matrix=M_Aq[2] ) += integrate ( markedfaces( mesh,mesh->markerName( "Tfixed" ) ),
                                                                    gamma*idt( t )*id( s )/hFace() );
    }
        
    M = M_backend->newMatrix( Xh, Xh );
    
    form2( Xh, Xh, M, _init=true ) = integrate( elements( mesh ), trans( id( v ) )*idt( u ) + trace( grad( v )*trans( gradt( u ) ) ) );
    
    form2( Xh, Xh, M ) += integrate( _range=elements( mesh ), _expr=id( q )*idt( p ) + grad( q )*trans( gradt( p ) ) );
    
    form2( Xh, Xh, M ) += integrate( _range=elements( mesh ), _expr=id( s )*idt( t ) + grad( s )*trans( gradt( t ) ) );
    
    M->close();
        
}

// \return the number of terms in affine decomposition of left hand
// side bilinear form
int Convection_crb::Qa() const
{
    return 4;
}

/**
 * there is at least one output which is the right hand side of the
 * primal problem
 *
 * \return number of outputs associated to the model
 */
int Convection_crb::Nl() const
{
    return 2;
}

/**
 * \param l the index of output
 * \return number of terms  in affine decomposition of the \p q th output term
 */
int Convection_crb::Ql( int l ) const
{
//    if ( l == 0 ) return 2;
    
    return 1;
}

/**
 * \brief compute the theta coefficient for both bilinear and linear form
 * \param mu parameter to evaluate the coefficients
 */
boost::tuple<theta_vector_type, std::vector<theta_vector_type> >
Convection_crb::computeThetaq( parameter_type const& mu, double )
{    
    M_thetaAq.resize( Qa() );
    M_thetaAq( 0 ) = - 1/math::sqrt( mu( 0 ) ); // k_1
    M_thetaAq( 1 ) = 1/( mu( 1 )*math::sqrt( mu ( 0 ) ) ); // k_2
    M_thetaAq( 2 ) = 1;
    M_thetaAq( 3 ) = 1;
    
    M_thetaFq.resize( Nl() );        
    M_thetaFq[0].resize( Ql( 0 ) );
    M_thetaFq[0]( 0 ) = 1. ; //mu( 2 ); // delta
    //    M_thetaFq[0]( 1 ) = mu( 3 ); // phi
    
    M_thetaFq[1].resize( Ql( 1 ) );
    M_thetaFq[1]( 0 ) = 1;

    return boost::make_tuple( M_thetaAq, M_thetaFq );
}

void
Convection_crb::update( parameter_type const& mu )
{    
    *D = *M_Aq[0];
    
    for ( size_type q = 1; q < M_Aq.size(); ++q )
    {
        //std::cout << "[affine decomp] scale q=" << q << " with " << M_thetaAq[q] << "\n";
        D->addMatrix( M_thetaAq[q], M_Aq[q] );
    }
    
    D->close();
    F->zero();
    
    for ( size_type q = 0; q < M_Fq[0].size(); ++q )
    {
        //std::cout << "[affine decomp] scale q=" << q << " with " << M_thetaFq[0][q] << "\n";
        F->add( M_thetaFq[0][q], M_Fq[0][q] );
    }
    F->close();
    
}

void Convection_crb ::updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J)
{
    Log() << "[updateJacobian] start\n";
    
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

    
    double gamma( this->vm()["penalbc"]. as<double>() );
    
    int adim=this->vm()["adim"]. as<int>();
    
    //conditions fortes de dir
    auto Rtemp =  M_backend->newVector( Xh );
    int weakdir( this->vm()["weakdir"]. as<int>() );
    double T0 = this->vm()["T0"]. as<double>();
        
    double pC = 1.0;
    M_Aq[3]->zero();

    // Fluid-NS
    // fluid convection derivatives: attention 2 terms
    
    form2( _test=Xh,_trial=Xh, _matrix=M_Aq[3] )  += integrate ( _range=elements( mesh ), _expr=cst( 1. )*trans( id( v ) )*( gradv( u ) )*idt( u ) );
    form2( _test=Xh,_trial=Xh, _matrix=M_Aq[3] )  += integrate ( _range=elements( mesh ), _expr=cst( 1. )*trans( id( v ) )*( gradt( u ) )*idv( u ) );

        
    //    // temperature derivatives
    //
    // heat convection by the fluid: attention 2 terms
    form2( _test=Xh, _trial=Xh, _matrix=M_Aq[3] ) +=
    integrate ( elements( mesh ),
               pC*grad( s )*( idv( t )*idt( u ) ) );
    
    form2( _test=Xh, _trial=Xh, _matrix=M_Aq[3] ) +=
    integrate ( elements( mesh ),
               pC*grad( s )*( idt( t )*idv( u ) ) );
    
    form2( _test=Xh, _trial=Xh, _matrix=M_Aq[3] ) +=
    integrate ( boundaryfaces( mesh ),
               pC*( trans( idv( u ) )*N() )*id( s )*idt( t ) );
    
    form2( _test=Xh, _trial=Xh, _matrix=M_Aq[3] ) +=
    integrate ( boundaryfaces( mesh ),
               pC*( trans( idt( u ) )*N() )*id( s )*idv( t ) );
    
    if ( weakdir == 0 )
    {
        //vitesse
        form2( Xh, Xh, M_Aq[3] )  += on( boundaryfaces( mesh ),u, Rtemp,one()*0. );
        
        if ( adim==1 )
            //temperature
            form2( Xh, Xh, M_Aq[3] )  += on ( markedfaces( mesh, "Tfixed" ),t,Rtemp,cst( 0.0 ) );
        
        else
            form2( Xh, Xh, M_Aq[3] )  += on ( markedfaces( mesh, "Tfixed" ),t,Rtemp,cst( T0 ) );
    }

    M_Aq[3]->close();

    J->zero();
    J->addMatrix(1.,D);
    J->addMatrix(1.,M_Aq[3]);
    
}


// instantiation
// class Convection_crb<2,1,2>;
