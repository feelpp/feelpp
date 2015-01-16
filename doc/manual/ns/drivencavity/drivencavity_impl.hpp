/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
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
   \author Christophe Prud'homme <prudhomme@unistra.fr>
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

    if (!J) J = backend(_name="newtonns")->newMatrix( Vh, Vh );
    auto a = form2( _test=Vh, _trial=Vh, _matrix=J );
    a = integrate( elements( mesh ), inner(gradt( u ),grad( v ))/Re );
    a += integrate( elements( mesh ), id(q)*divt(u) -idt(p)*div(v) );
    // Convective terms
    a += integrate( elements( mesh ), trans(id(v))*gradv(u)*idt(u));
    a += integrate( elements( mesh ), trans(id(v))*gradt(u)*idv(u));

    //#if defined( FEELPP_USE_LM )
    a += integrate(elements(mesh), id(q)*idt(lambda)+idt(p)*id(nu));
    //#elif
    //a += integrate(elements(mesh), idt(p)*id(nu));

    //Weak Dirichlet conditions
    a += integrate( boundaryfaces( mesh ),-trans( -idt(p)*N()+gradt(u)*N()/Re )*id( v ));//
    a += integrate( boundaryfaces( mesh ),-trans( -id(p)*N()+grad(u)*N()/Re )*idt( u ));//
    a += integrate( boundaryfaces( mesh ), +penalbc*inner( idt( u ),id( v ) )/hFace() );


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
    auto u_exact=vf::project(Vh->template functionSpace<0>(), markedfaces(mesh, "wall2"), uex );


    U=*X;
    auto r = form1( _test=Vh, _vector=R );
    //r += integrate( elements( mesh ),-inner( f,id( v ) ) );
    r = integrate( elements( mesh ), trans(gradv( u )*idv(u))*id(v));//convective term
    r += integrate( elements( mesh ), inner(gradv( u ),grad( v ))/Re );
    r +=  integrate( elements( mesh ),-idv(p)*div(v) + id(q)*divv(u));
    //#if defined( FEELPP_USE_LM )
    r += integrate ( elements( mesh ), +id( q )*idv( lambda )+idv( p )*id( nu ) );
    //#endif


    //Weak Dirichlet
    auto SigmaNv = ( -idv( p )*N() + gradv( u )*N()/Re );
    auto SigmaN = ( -id( q )*N() + grad( v )*N()/Re );
    r +=integrate ( boundaryfaces(mesh), - trans( SigmaNv )*id( v ) - trans( SigmaN )*( idv( u ) -idv(u_exact) ) + penalbc*trans( idv( u ) - idv(u_exact) )*id( v )/hFace() );


}

template<int Dim>
void DrivenCavity<Dim>::exportResults( element_type const& U )
{
    auto uex=unitX();
    auto u_exact=vf::project(Vh->template functionSpace<0>(), markedfaces(mesh, "wall2"), uex );

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



    u=vf::project(Vh->template functionSpace<0>(), elements(mesh), zero<Dim,1>());
    p=vf::project(Vh->template functionSpace<1>(), elements(mesh), constant(0.0));

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
