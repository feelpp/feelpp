/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-05-27

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file stvenant_kirchhoff_impl.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-05-27
 */
#ifndef __STVENANT_KICHHOFF_IMPL_H
#define __STVENANT_KICHHOFF_IMPL_H 1


#include <stvenant_kirchhoff.hpp>

namespace Feel
{
template<int Dim, int Order>
StVenantKirchhoff<Dim,Order>::StVenantKirchhoff( po::variables_map const& vm )
    :
    super( Dim, vm ),
    M_backend( backend_type::build( this->vm() ) ),
    M_Xh(),
    exporter( Exporter<mesh_type>::New( soption("exporter") )->setOptions( this->vm() ) ),
    timeSet( new timeset_type( "stvenant_kirchhoff" ) )
{
    //timeSet->setTimeIncrement( this->dt()/this->nSubSteps() );
    timeSet->setTimeIncrement( this->dt() );
    exporter->addTimeSet( timeSet );
    exporter->setPrefix( "stvenant_kirchhoff" );

    /**
     * Physical data
     */
    E = 21*1e5;
    sigma = 0.28;
    mu = E/( 2*( 1+sigma ) );
    lambda = E*sigma/( ( 1+sigma )*( 1-2*sigma ) );
    density = 1;
    gravity = -density*0.05;


    mesh_ptrtype mesh = loadMesh();

    M_Xh = functionspace_ptrtype( functionspace_type::New( mesh ) );
    un2 = element_ptrtype( new element_type( M_Xh, "un2" ) );
    un1 = element_ptrtype( new element_type( M_Xh, "un1" ) );
    un = element_ptrtype( new element_type( M_Xh, "un" ) );

    M_oplin = oplin_ptrtype( new oplin_type( M_Xh, M_Xh, M_backend ) );
    M_jac = oplin_ptrtype( new oplin_type( M_Xh, M_Xh, M_backend ) );
    M_residual = funlin_ptrtype( new funlin_type( M_Xh, M_backend ) );

    M_Dh = displacement_functionspace_ptrtype( displacement_functionspace_type::New( mesh ) );
    M_opelas = opdisplacement_ptrtype( new opdisplacement_type( M_Dh, M_Dh, M_backend ) );
    M_lfelas = fundisplacement_ptrtype( new fundisplacement_type( M_Dh, M_backend ) );

    // Operator Lag P1
    M_displacement_oplagp1 = displacement_oplagp1_ptrtype( new displacement_oplagp1_type( M_Dh, M_backend ) );

    LOG(INFO) << "[Turek::initLinearOperators] done\n";

}

template<int Dim, int Order>
typename StVenantKirchhoff<Dim,Order>::mesh_ptrtype
StVenantKirchhoff<Dim,Order>::loadMesh()
{
    mesh_ptrtype mesh( new mesh_type );

    ImporterGmsh<mesh_type> import( this->createMesh() );
    import.setVersion( "2.0" );
    mesh->accept( import );

    return mesh;
}


template<int Dim, int Order>
void
StVenantKirchhoff<Dim, Order>::updateResidual( const vector_ptrtype& X, vector_ptrtype& R )
{
    boost::timer ti;
    LOG(INFO) << "[updateResidual] start\n";

    mesh_ptrtype mesh = M_Xh->mesh();

    element_type U( M_Xh, "U" );
    element_type V( M_Xh, "V" );
    element_u_type u = U.template element<0>();
    element_u_type v = V.template element<0>();
    element_v_type uu = U.template element<1>();
    element_v_type vv = V.template element<1>();
    U = *X;

    im_type im;

    AUTO( g, constant( 0.0 ) );
    AUTO( defv, 0.5*( gradv( u )+trans( gradv( u ) ) ) );
    AUTO( def, 0.5*( grad( v )+trans( grad( v ) ) ) );
    AUTO( Id, ( mat<Dim,Dim>( cst( 1 ), cst( 0 ), cst( 0 ), cst( 1. ) ) ) );
    //std::cout << "u = " << u << "\n";


    *M_residual =
        integrate( elements( mesh ), im,
                   .5*mu*( trace( ( gradv( u )*trans( gradv( u ) ) )*grad( v ) ) )+
                   .25*lambda*trace( gradv( u )*trans( gradv( u ) ) )*div( v ) -
                   trans( gravity*oneY() )*id( v ) );

#if 0
    //AUTO( eta, 0.1*Px()*( Px() -5 )*(Px()-2.5)*sin( omega*M_PI*cst_ref(time)  ) );

    // force applied at the bottom
    *M_residual +=
        integrate( markedfaces( mesh, 2 ), im,
                   -trans( eta*oneY() )*id( v ) );
#endif

    *M_residual +=
        integrate( elements( mesh ), im,
                   - density*trans( idv( M_bdf->derivate( this->timeOrder(), this->dt() ).template element<1>() ) ) *id( v )
                   //-density*trans(2*idv(un->template element<0>())-idv(un1->template element<0>())) *id(v) /(this->dt()*this->dt())
                 );
    *M_residual +=
        integrate( elements( mesh ), im,
                   //+ trans(idv( u ))*id(vv)*M_bdf->derivateCoefficient( this->timeOrder(), 0, this->this->dt()() )
                   - trans( idv( M_bdf->derivate( this->timeOrder(), this->dt() ).template element<0>() ) )*id( vv )
                 );
    FsFunctionalLinear<functionspace_type> flin( M_Xh, M_backend );
    M_oplin->apply( U, flin );





    M_residual->add( flin );
    M_residual->close();
    *R = M_residual->container();
    LOG(INFO) << "[updateResidual] done in " << ti.elapsed() << "s\n";
}
template<int Dim, int Order>
void
StVenantKirchhoff<Dim, Order>::updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J )
{
    boost::timer ti;
    LOG(INFO) << "[updateJacobian] start\n";
    static bool is_init = false;
    mesh_ptrtype mesh = M_Xh->mesh();

    element_type U( M_Xh, "U" );
    element_type V( M_Xh, "V" );
    element_u_type u = U.template element<0>();
    element_u_type v = V.template element<0>();
    element_v_type uu = U.template element<1>();
    element_v_type vv = V.template element<1>();

    U = *X;

    im_type im;

    if ( is_init == false )
    {
        *M_jac = integrate( elements( mesh ), im,
                            .5*mu*( trace( ( gradv( u )*trans( gradt( u ) ) )*grad( v ) ) )+
                            .25*lambda*trace( gradv( u )*trans( gradt( u ) ) )*div( v )

                            + 0*trans( idt( u ) )*id( vv )
                            + 0*trans( idt( uu ) )*id( v )
                            + 0*trans( idt( uu ) )*id( vv )

                          );

        is_init = true;
    }

    else
    {
        M_jac->matPtr()->zero();
        *M_jac += integrate( elements( mesh ), im,
                             .5*mu*( trace( ( gradv( u )*trans( gradt( u ) ) )*grad( v ) ) )+
                             .25*lambda*trace( gradv( u )*trans( gradt( u ) ) )*div( v ) );
    }

    M_jac->close();
    M_jac->matPtr()->addMatrix( 1.0, M_oplin->mat() );
    J = M_jac->matPtr();
    LOG(INFO) << "[updateJacobian] done in " << ti.elapsed() << "s\n";
}
template<int Dim, int Order>
void
StVenantKirchhoff<Dim, Order>::updateResidualJacobian( const vector_ptrtype& X, vector_ptrtype& R, sparse_matrix_ptrtype& J )
{
}

template<int Dim, int Order>
void
StVenantKirchhoff<Dim, Order>::initElastoStaticProblem()
{
    displacement_element_type u( M_Dh, "u" );
    displacement_element_type v( M_Dh, "v" );

    mesh_ptrtype mesh = M_Dh->mesh();
    im_type im;
    *M_lfelas = integrate( elements( mesh ), im, trans( gravity*oneY() )*id( v ) );
    M_lfelas->close();


    AUTO( deft, 0.5*( gradt( u )+trans( gradt( u ) ) ) );
    AUTO( def, 0.5*( grad( v )+trans( grad( v ) ) ) );
    AUTO( Id, ( mat<Dim,Dim>( cst( 1 ), cst( 0 ), cst( 0 ), cst( 1. ) ) ) );

    *M_opelas =
        integrate( elements( mesh ), im,
                   lambda*divt( u )*div( v )  +
                   2*mu*trace( trans( deft )*def ) );

    BOOST_FOREACH( std::string marker, this->dirichletMarkers() )
    {
        LOG(INFO) << "[LinearElasticty::Dirichlet] weakbc boundary " << marker << " id : " << mesh->markerName( marker ) << "\n";
        *M_opelas +=
            integrate( markedfaces( mesh,mesh->markerName( marker ) ), im,
                       - trans( ( 2*mu*deft+lambda*trace( deft )*Id )*N() )*id( v )
                       - trans( ( 2*mu*def+lambda*trace( def )*Id )*N() )*idt( u )
                       + this->gammaBc()*trans( idt( u ) )*id( v )/hFace() );
    }

    M_opelas->close();

}
template<int Dim, int Order>
void
StVenantKirchhoff<Dim, Order>::initLinearPart()
{
    element_type U( M_Xh, "U" );
    element_type V( M_Xh, "V" );
    element_u_type u = U.template element<0>();
    element_u_type v = V.template element<0>();
    element_v_type uu = U.template element<1>();
    element_v_type vv = V.template element<1>();


    AUTO( deft, 0.5*( gradt( u )+trans( gradt( u ) ) ) );
    AUTO( def, 0.5*( grad( v )+trans( grad( v ) ) ) );
    AUTO( Id, ( mat<Dim,Dim>( cst( 1 ), cst( 0 ), cst( 0 ), cst( 1. ) ) ) );

    mesh_ptrtype mesh = M_Xh->mesh();
    im_type im;

    *M_oplin =
        integrate( elements( mesh ), im,
                   density*trans( idt( uu ) )*id( v )*M_bdf->derivateCoefficient( this->timeOrder(), 0, this->dt() ) +

                   lambda*divt( u )*div( v )  +
                   2*mu*trace( trans( deft )*def )



                   + trans( idt( u ) )*id( vv )*M_bdf->derivateCoefficient( this->timeOrder(), 0, this->dt() )

                   - trans( idt( uu ) )*id( vv )
                 );

    BOOST_FOREACH( std::string marker, this->dirichletMarkers() )
    {
        LOG(INFO) << "[LinearPart::Dirichlet] weakbc boundary " << marker << " id : " << mesh->markerName( marker ) << "\n";
        *M_oplin +=
            integrate( markedfaces( mesh,mesh->markerName( marker ) ), im,
                       - trans( ( 2*mu*deft+lambda*trace( deft )*Id )*N() )*id( v )
                       - trans( ( 2*mu*def+lambda*trace( def )*Id )*N() )*idt( u )
                       + this->gammaBc()*trans( idt( u ) )*id( v )/hFace() );
    }

    M_oplin->close();
    //M_oplin->mat().printMatlab( "oplin.m" );

}

template<int Dim, int Order>
void
StVenantKirchhoff<Dim, Order>::run()
{
    using namespace Feel::vf;
    mesh_ptrtype mesh = M_Xh->mesh();

    element_type U( M_Xh, "U" );
    element_type V( M_Xh, "V" );
    element_u_type u = U.template element<0>();
    element_u_type v = V.template element<0>();
    element_v_type uu = U.template element<1>();
    element_v_type vv = V.template element<1>();

    M_bdf = bdf_ptrtype( new bdf_type( M_Xh ) );

    // init linear elastostatic problem
    initElastoStaticProblem();

    // init linear part of the StVenantKirchhoff model
    initLinearPart();



    M_backend->nlSolver()->residual = boost::bind( &self_type::updateResidual, boost::ref( *this ), _1, _2 );
    M_backend->nlSolver()->jacobian = boost::bind( &self_type::updateJacobian, boost::ref( *this ), _1, _2 );

    //u = project( M_Xh, elements(mesh), constant(0.)*one() );
    U.zero();

    un->zero();
    un1->zero();


    M_bdf->initialize( U );
    im_type im;
    // solve the elastostatic problem
    displacement_element_type d( M_Dh, "d" );
    this->linearSolve( M_opelas->matPtr(), d, M_lfelas->containerPtr() );

    double d_norm2 = math::sqrt( integrate( elements( mesh ), im, trans( idv( d ) )*( idv( d ) ) ).evaluate()( 0, 0 ) );

    time = 0;
    exportResults( time, U, d );


    vector_ptrtype Un( M_backend->newVector( U.functionSpace() ) );
    vector_ptrtype R( M_backend->newVector( U.functionSpace() ) );
    sparse_matrix_ptrtype J;


    boost::timer ttotal;
    int iterations = 0;
    double error1 = 1;
    double error2 = 1;
    double error3 = 1;

    //for( time = this->dt(), iterations = 0; error1 >= 1e-6 || iterations <= 2; time +=this->dt(), ++iterations )
    for ( time += this->dt(), iterations = 1; time <= this->T(); time +=this->dt(), ++iterations )
        //while ( error1 >= 1e-6 && iterations <= 2 )
    {
        boost::timer ti;
        LOG(INFO) << "============================================================\n";
        LOG(INFO) << "time: " << time << "s, iteration: " << iterations << "\n";

        *Un = U;
        this->updateResidual( Un, R );
        this->updateJacobian( Un, J );


        solve( J, U, R );

        *un1 = *un;
        *un = U;

        M_bdf->shiftRight( U );

        error1 = math::sqrt( integrate( elements( mesh ), im, trans( idv( un->template element<0>() )-idv( un1->template element<0>() ) )*( idv( un->template element<0>() )-idv( un1->template element<0>() ) ) ).evaluate()( 0, 0 ) )/d_norm2;
        error2 = math::sqrt( integrate( elements( mesh ), im, trans( idv( U.template element<0>() )-idv( d ) )*( idv( U.template element<0>() )-idv( d ) ) ).evaluate()( 0, 0 ) )/d_norm2;


        error3 = math::sqrt( integrate( elements( mesh ), im,
                                        trans( idv( M_bdf->unknown( 0 ).template element<0>() )-idv( M_bdf->unknown( 1 ).template element<0>() ) )*
                                        ( idv( M_bdf->unknown( 0 ).template element<0>() )-idv( M_bdf->unknown( 1 ).template element<0>() ) ) ).evaluate()( 0, 0 ) )/d_norm2;

        LOG(INFO) << "                      ||d||_0 = " << d_norm2 << "\n";
        LOG(INFO) << "                   ||d(t)||_0 = " << math::sqrt( integrate( elements( mesh ), im, trans( idv( un->template element<0>() ) )*idv( un->template element<0>() ) ).evaluate()( 0, 0 ) ) << "\n";
        LOG(INFO) << "||dh(t+1) - dh(t)||_0/||d||_0 = " << error1 << "\n";
        LOG(INFO) << "     ||dh(t) - dh||_0/||d||_0 = " << error2 << "\n";
        LOG(INFO) << "||dh(t+1) - dh(t)||_0/||d||_0 = (bdf) " << error3 << "\n";

        exportResults( time, U, d );

        LOG(INFO) << "time spent in iteration :  " << ti.elapsed() << "s\n";

        if ( error1  <= 1e-8 && iterations > 2 )
            break;
    }

    LOG(INFO) << "total time spent :  " << ttotal.elapsed() << "s\n";
    LOG(INFO) << "total number of iterations :  " << iterations << "\n";


} // StVenantKirchhoff::run

template<int Dim, int Order>
void
StVenantKirchhoff<Dim, Order>::linearSolve( sparse_matrix_ptrtype& D,
        displacement_element_type& u,
        vector_ptrtype& F )
{
    vector_ptrtype U( M_backend->newVector( u.functionSpace() ) );
    *U = u;
    M_backend->solve( D, D, U, F );
    u = *U;


} // StVenantKirchhoff::solve

template<int Dim, int Order>
void
StVenantKirchhoff<Dim, Order>::solve( sparse_matrix_ptrtype& D,
                                      element_type& u,
                                      vector_ptrtype& F )
{
    vector_ptrtype U( M_backend->newVector( u.functionSpace() ) );
    *U = u;
    M_backend->nlSolve( D, U, F, 1e-10, 10 );
    u = *U;


} // StVenantKirchhoff::solve


template<int Dim, int Order>
void
StVenantKirchhoff<Dim, Order>::exportResults( double time, element_type& U, displacement_element_type& d )
{

    LOG(INFO) << "exportResults starts\n";

    element_type S( M_Xh, "S" );
    S = U;
    S.add( -1, M_bdf->unknown( 0 ) );

    double subdt = this->dt()/this->nSubSteps();;
    //for( uint16_type i = 1; i <= this->nSubSteps(); ++i )
    {
        //typename timeset_type::step_ptrtype timeStep = timeSet->step( time-this->dt()+i*subdt );
        typename timeset_type::step_ptrtype timeStep = timeSet->step( time );
        timeStep->setMesh( M_displacement_oplagp1->dualImageSpace()->mesh() );
        timeStep->add( "pid",
                       regionProcess( boost::shared_ptr<p0_space_type>( new p0_space_type( U.functionSpace()->mesh() ) ) ) );


        element_type V( M_Xh, "V" );
        V = U;
        //V.add( double(i)/subdt, S );
        timeStep->add( "displ", U.template element<0>() );
        timeStep->add( "veloc", U.template element<1>() );
        timeStep->add( "displ_static", d );

        mesh_ptrtype mesh = M_Xh->mesh();
        im_type im;
        double length = integrate( markedfaces( mesh,mesh->markerName( "right" ) ), im, constant( 1.0 ) ).evaluate()( 0, 0 );
        double avg_displ = integrate( markedfaces( mesh,mesh->markerName( "right" ) ), im, sqrt( trans( idv( U.template element<0>() ) )*idv( U.template element<0>() ) ) ).evaluate()( 0, 0 )/length;
        double avg_displ_static = integrate( markedfaces( mesh,mesh->markerName( "right" ) ), im, sqrt( trans( idv( d ) )*idv( d ) ) ).evaluate()( 0, 0 )/length;

        LOG(INFO) << "[postProcess]                  length = " << length << "\n";
        LOG(INFO) << "[postProcess]        avg_displacement = " << avg_displ << "\n";
        LOG(INFO) << "[postProcess] avg_displacement_static = " << avg_displ_static << "\n";
        timeStep->addScalar( "avg_displ", avg_displ );
        timeStep->addScalar( "avg_displ_static", avg_displ_static );

        exporter->save();
        timeStep->setState( STEP_ON_DISK );
    }
} // StVenantKirchhoff::export

}
#endif /* __STVENANT_KICHHOFF_IMPL_H */
