/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-11-22

  Copyright (C) 2013 Université de Strasbourg

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
   \file drivencavity_impl.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-11-22
 */
#ifndef FEELPP_DRIVENCAVITY_IMPL_HPP_H
#define FEELPP_DRIVENCAVITY_IMPL_HPP_H 1

#include "drivencavity.hpp"

namespace Feel
{
template<int Dim>
DrivenCavity<Dim>::DrivenCavity( )
    :
    super( ),
    Re(      doption("Re") ),
    penalbc( doption("bccoeff") ),
    exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
{

}

template<int Dim>
void DrivenCavity<Dim>::init()
{
    if ( this->vm().count( "help" ) )
    {
        if ( Environment::worldComm().isMasterRank() )
            std::cout << this->optionsDescription() << "\n";
        return;
    }
    mesh = loadMesh(_mesh=new Mesh<Simplex<Dim>>);
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "number of elements of " << Dim << "D: " << mesh->numElements() << "\n";

    Vh = space_type::New( mesh );
}

template<int Dim>
void
DrivenCavity<Dim>::Jacobian(const vector_ptrtype& X, sparse_matrix_ptrtype& J)
{
    auto U = Vh->element( "(u,p)" );
    auto V = Vh->element( "(v,q)" );
    auto u = U.template element<0>( "u" );
    auto v = V.template element<0>( "u" );
    auto p = U.template element<1>( "p" );
    auto q = V.template element<1>( "p" );
    //#if defined( FEELPP_USE_LM )
    auto lambda = U.template element<2>();
    auto nu = V.template element<2>();

    //#endif

    if (!J) J = backend(_name="newtonns")->newMatrix( _test=Vh, _trial=Vh );
    auto a = form2( _test=Vh, _trial=Vh, _matrix=J );
    a = integrate( _range=elements( mesh ), _expr=inner(gradt( u ),grad( v ))/Re );
    a += integrate( _range=elements( mesh ), _expr=id(q)*divt(u) -idt(p)*div(v) );
    // Convective terms
    a += integrate( _range=elements( mesh ), _expr=trans(id(v))*gradv(u)*idt(u));
    a += integrate( _range=elements( mesh ), _expr=trans(id(v))*gradt(u)*idv(u));

    //#if defined( FEELPP_USE_LM )
    a += integrate(_range=elements(mesh), _expr=id(q)*idt(lambda)+idt(p)*id(nu));
    //#elif
    //a += integrate(_range=elements(mesh), _expr=idt(p)*id(nu));

    //Weak Dirichlet conditions
    a += integrate( _range=boundaryfaces( mesh ),_expr= -trans( -idt(p)*N()+gradt(u)*N()/Re )*id( v ));//
    a += integrate( _range=boundaryfaces( mesh ),_expr= -trans( -id(p)*N()+grad(u)*N()/Re )*idt( u ));//
    a += integrate( _range=boundaryfaces( mesh ),_expr= +penalbc*inner( idt( u ),id( v ) )/hFace() );


}

template<int Dim>
void
DrivenCavity<Dim>::Residual(const vector_ptrtype& X, vector_ptrtype& R)
{
    auto U = Vh->element( "(u,p)" );
    auto V = Vh->element( "(v,q)" );
    auto u = U.template element<0>( "u" );
    //auto u_exact = U.template element<0>( "u_exact" );
    auto v = V.template element<0>( "u" );
    auto p = U.template element<1>( "p" );
    auto q = V.template element<1>( "p" );
    //#if defined( FEELPP_USE_LM )
    auto lambda = U.template element<2>();
    auto nu = V.template element<2>();
    //#endif

    auto uex=unitX();
    auto u_exact=vf::project(_space=Vh->template functionSpace<0>(), _range=markedfaces(mesh, "wall2"), _expr=uex );


    U=*X;
    auto r = form1( _test=Vh, _vector=R );
    //r += integrate( elements( mesh ),-inner( f,id( v ) ) );
    r = integrate( _range=elements( mesh ), _expr= trans(gradv( u )*idv(u))*id(v));//convective term
    r += integrate( _range=elements( mesh ), _expr= inner(gradv( u ),grad( v ))/Re );
    r +=  integrate( _range=elements( mesh ),_expr= -idv(p)*div(v) + id(q)*divv(u));
    //#if defined( FEELPP_USE_LM )
    r += integrate ( _range=elements( mesh ), _expr= +id( q )*idv( lambda )+idv( p )*id( nu ) );
    //#endif


    //Weak Dirichlet
    auto SigmaNv = ( -idv( p )*N() + gradv( u )*N()/Re );
    auto SigmaN = ( -id( q )*N() + grad( v )*N()/Re );
    r +=integrate ( _range=boundaryfaces(mesh), _expr= - trans( SigmaNv )*id( v ) - trans( SigmaN )*( idv( u ) -idv(u_exact) ) + penalbc*trans( idv( u ) - idv(u_exact) )*id( v )/hFace() );


}

template<int Dim>
void DrivenCavity<Dim>::exportResults( element_type const& U )
{
    auto uex=unitX();
    auto u_exact=vf::project(_space=Vh->template functionSpace<0>(), _range=markedfaces(mesh, "wall2"), _expr=uex );

    if ( exporter->doExport() )
    {
        exporter->step( 0 )->setMesh( U.functionSpace()->mesh() );
        exporter->step( 0 )->addRegions();
        exporter->step( 0 )->add( "u", U.template element<0>() );
        exporter->step( 0 )->add( "p", U.template element<1>() );
        exporter->step( 0 )->add( "uex", u_exact);
        exporter->save();
    }

}
template<int Dim>
void DrivenCavity<Dim>::run()
{
    this->init();

    auto U = Vh->element( "(u,p)" );
    auto V = Vh->element( "(v,q)" );
    auto u = U.template element<0>( "u" );
    auto v = V.template element<0>( "u" );
    auto p = U.template element<1>( "p" );
    auto q = V.template element<1>( "p" );
    //#if defined( FEELPP_USE_LM )
    auto lambda = U.template element<2>();
    auto nu = V.template element<2>();
    //#endif



    u.on( _range=elements(mesh), _expr=zero<Dim,1>());
    p.on( _range=elements(mesh), _expr=constant(0.0));

    if ( Environment::worldComm().isMasterRank() )
        std::cout << "Initializing residual " << "\n";
    backend(_name="newtonns")->nlSolver()->residual = boost::bind( &DrivenCavity::Residual,
                                                                   boost::ref( *this ), _1, _2 );
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "Initializing the Jacobian matrix" << "\n";
    backend(_name="newtonns")->nlSolver()->jacobian = boost::bind( &DrivenCavity::Jacobian,
                                                                   boost::ref( *this ), _1, _2 );

    if ( boption(_name="continuation" ) )
    {
        double ReTarget = Re;
        int N = std::ceil( std::log( Re ) );
        for( int i  = 0; i <= N; ++i )
        {
            Re = std::exp( std::log(1)+i*(std::log(ReTarget)-std::log(1))/double(N));
            if ( Environment::worldComm().isMasterRank() )
                std::cout << "Start solving for Reynolds = " << Re << "\n";
            backend(_name="newtonns")->nlSolve( _solution=U );
            if ( Environment::worldComm().isMasterRank() )
                std::cout << "Done!" << "\n";
        }
    }
    else
    {
        if ( Environment::worldComm().isMasterRank() )
            std::cout << "Start solving for Reynolds = " << Re << "\n";
        backend(_name="newtonns")->nlSolve( _solution=U );
        if ( Environment::worldComm().isMasterRank() )
            std::cout << "Done!" << "\n";
    }

    this->exportResults( U );
}
} // namespace Feel
#endif /* FEELPP_DRIVENCAVITY_IMPL_HPP_H */
