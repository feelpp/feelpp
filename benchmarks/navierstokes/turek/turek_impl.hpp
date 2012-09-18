/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-05-29

  Copyright (C) 2008,2011 Université Joseph Fourier (Grenoble I)

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
   \file turek_impl.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-05-29
 */
#ifndef __TurekImpl_H
#define __TurekImpl_H 1

#include <turek.hpp>

#include <feel/feelfilters/importergmsh.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>

#include <feel/feelvf/vf.hpp>

namespace Feel
{
using namespace vf;


template<int Dim, int Order, int GeoOrder>
Turek<Dim, Order, GeoOrder>::Turek( po::variables_map const& vm )
    :
    super( Dim, vm ),
    M_backend( backend_type::build( this->vm() ) ),
    M_backend_symm_v( backend_type::build( this->vm() ) ),
    M_backend_symm_s( backend_type::build( this->vm() ) ),
    M_Xh(),
    exporter( Exporter<mesh_type>::New( this->vm(), Data::makeAbout().appName() ) ),
    M_data( "data.txt" )
{

    // print data
    this->print();

    mesh_ptrtype mesh = createMesh();

    M_Xh = fluid_functionspace_ptrtype( fluid_functionspace_type::New( mesh ) );
    Un1 = fluid_element_ptrtype( new fluid_element_type( M_Xh, "un1" ) );
    Un = fluid_element_ptrtype( new fluid_element_type( M_Xh, "un" ) );
    LOG(INFO) << "[Turek] Constructor done\n";

    M_backend->nlSolver()->residual = boost::bind( &self_type::updateResidual, boost::ref( *this ), _1, _2 );
    M_backend->nlSolver()->jacobian = boost::bind( &self_type::updateJacobian, boost::ref( *this ), _1, _2 );

    // jacobian and residual
    M_jac = oplin_ptrtype( new oplin_type( M_Xh, M_Xh, M_backend ) );
    M_residual = funlin_ptrtype( new funlin_type( M_Xh, M_backend ) );
    M_stokes_rhs = funlin_ptrtype( new funlin_type( M_Xh, M_backend ) );


    // Operator Lag P1
    M_velocity_oplagp1 = velocity_oplagp1_ptrtype( new velocity_oplagp1_type( M_Xh->template functionSpace<0>(), M_backend ) );
    M_pressure_oplagp1 = pressure_oplagp1_ptrtype( new pressure_oplagp1_type( M_Xh->template functionSpace<1>(), M_backend ) );

    this->initLinearOperators();

    M_bdf = bdf_ptrtype( new bdf_type( this->vm(),
                                       M_Xh,
                                       "fluid" ) );
    M_bdf->print();

    this->createMatlabScript();
}

template<int Dim, int Order, int GeoOrder>
typename Turek<Dim, Order, GeoOrder>::mesh_ptrtype
Turek<Dim, Order, GeoOrder>::createMesh()
{
    LOG(INFO) << "[Turek] CreateMesh starts\n";
    mesh_ptrtype mesh( new mesh_type );
    //mesh->setRenumber( false );

    Gmsh gmsh;
    gmsh.setOrder( GeoOrder );
    std::string mesh_name, mesh_desc;
    boost::tie( mesh_name, mesh_desc ) = this->createCylinder();
    std::string fname = gmsh.generate( mesh_name, mesh_desc );

    ImporterGmsh<mesh_type> import( fname );
    import.setVersion( "2.0" );
    mesh->accept( import );

    LOG(INFO) << "[Turek] CreateMesh done\n";
    return mesh;
} // Turek::createMesh

template<int Dim, int Order, int GeoOrder>
void
Turek<Dim, Order, GeoOrder>::nlsolve( sparse_matrix_ptrtype& D,
                                      fluid_element_type& u,
                                      vector_ptrtype& F )
{
    boost::timer ti;
    vector_ptrtype U( M_backend->newVector( u.functionSpace() ) );
    *U = u;

    if ( this->useSamePreconditioner() )
        M_backend->setPrecMatrixStructure( SAME_PRECONDITIONER );

    else
        M_backend->setPrecMatrixStructure( SAME_NONZERO_PATTERN );

    M_backend->nlSolve( D, U, F, 1e-10, 10 );
    //M_backend->solve( D, D, U, F );
    u = *U;
    LOG(INFO) << "time spent in nonlinear solve :  " << ti.elapsed() << "s\n";

} // Turek::solve

template<int Dim, int Order, int GeoOrder>
void
Turek<Dim, Order, GeoOrder>::solve( sparse_matrix_ptrtype& D,
                                    fluid_element_type& u,
                                    vector_ptrtype& F )
{
    boost::timer ti;
    vector_ptrtype U( M_backend->newVector( u.functionSpace() ) );
    *U = u;
    //M_backend->nlSolve( D, U, F, 1e-10, 10 );
    bool conv;
    int its;
    double res;
    boost::tie( conv, its, res ) = M_backend->solve( D, D, U, F );
    u = *U;
    LOG(INFO) << "[linear solve]           converged : " <<  conv << "\n";
    LOG(INFO) << "[linear solve] number of iterations: " <<  its << "\n";
    LOG(INFO) << "[linear solve]             residual: " <<  res << "\n";
    LOG(INFO) << "time spent in linear solve :  " << ti.elapsed() << "s\n";

} // Turek::solve

template<int Dim, int Order, int GeoOrder>
void
Turek<Dim, Order, GeoOrder>::solve( sparse_matrix_ptrtype& D,
                                    velocity_element_type& u,
                                    vector_ptrtype& F )
{
    boost::timer ti;
    vector_ptrtype U( M_backend_symm_v->newVector( u.functionSpace() ) );
    *U = u;
    //M_backend->nlSolve( D, U, F, 1e-10, 10 );
    bool conv;
    int its;
    double res;
    boost::tie( conv, its, res ) = M_backend_symm_v->solve( D, D, U, F );
    u = *U;
    LOG(INFO) << "[linear solve]           converged : " <<  conv << "\n";
    LOG(INFO) << "[linear solve] number of iterations: " <<  its << "\n";
    LOG(INFO) << "[linear solve]             residual: " <<  res << "\n";
    LOG(INFO) << "time spent in linear solve :  " << ti.elapsed() << "s\n";

} // Turek::solve

template<int Dim, int Order, int GeoOrder>
void
Turek<Dim, Order, GeoOrder>::solve( sparse_matrix_ptrtype& D,
                                    pressure_element_type& u,
                                    vector_ptrtype& F )
{
    boost::timer ti;
    vector_ptrtype U( M_backend_symm_s->newVector( u.functionSpace() ) );
    *U = u;
    //M_backend->nlSolve( D, U, F, 1e-10, 10 );
    bool conv;
    int its;
    double res;
    boost::tie( conv, its, res ) = M_backend_symm_s->solve( D, D, U, F );
    u = *U;
    LOG(INFO) << "[linear solve]           converged : " <<  conv << "\n";
    LOG(INFO) << "[linear solve] number of iterations: " <<  its << "\n";
    LOG(INFO) << "[linear solve]             residual: " <<  res << "\n";
    LOG(INFO) << "time spent in linear solve :  " << ti.elapsed() << "s\n";

} // Turek::solve





template<int Dim, int Order, int GeoOrder>
void
Turek<Dim, Order, GeoOrder>::initLinearOperators()
{
    boost::timer ti, total_time;
    LOG(INFO) << "[Turek::initLinearOperators] start\n";
    using namespace Feel::vf;
    mesh_ptrtype mesh = M_Xh->mesh();

    fluid_element_type U( M_Xh, "U" );
    fluid_element_type V( M_Xh, "V" );

    fluid_element_0_type u = U.template element<0>();
    fluid_element_1_type p = U.template element<1>();
    //fluid_element_2_type lambda = U.template element<2>();

    fluid_element_0_type v = V.template element<0>();
    fluid_element_1_type q = V.template element<1>();

    LOG(INFO) << "[Turek::initLinearOperators] space+elements init done in " << ti.elapsed() << "s\n";
    ti.restart();

    M_mass_v = op_vector_ptrtype( new op_vector_type( M_Xh->template functionSpace<0>(), M_Xh->template functionSpace<0>(), M_backend ) );
    M_mass_s = op_scalar_ptrtype( new op_scalar_type( M_Xh->template functionSpace<1>(), M_Xh->template functionSpace<1>(), M_backend ) );
    M_oplin = oplin_ptrtype( new oplin_type( M_Xh, M_Xh, M_backend ) );
    AUTO( deft, 0.5*( gradt( u )+trans( gradt( u ) ) ) );
    AUTO( def, 0.5*( grad( v )+trans( grad( v ) ) ) );
    AUTO( defv, 0.5*( gradv( u )+trans( gradv( u ) ) ) );
    AUTO( Id, ( mat<Dim,Dim>( cst( 1 ), cst( 0 ), cst( 0 ), cst( 1. ) ) ) );
    AUTO( SigmaNt, ( -idt( p )*N()+2*this->nu()*deft*N() ) );
    AUTO( SigmaN, ( -id( p )*N()+2*this->nu()*def*N() ) );
    AUTO( SigmaNv, ( -idv( p )*N()+2*this->nu()*defv*N() ) );

    *M_mass_v = integrate( elements( mesh ), _Q<2*uOrder+2*( GeoOrder-1 )>(), trans( idt( u ) )*id( v ) );
    M_mass_v->close();
    *M_mass_s = integrate( elements( mesh ), _Q<2*uOrder+2*( GeoOrder-1 )>(), trans( idt( p ) )*id( q ) );
    M_mass_s->close();

    // oplin
    LOG(INFO) << "[add element stokes terms] (nu * (nabla u+ nabla^T u)  : nabla v)\n";
    *M_oplin =
        integrate( elements( mesh ), _Q<2*( uOrder-1 )+4*( GeoOrder-1 )>(),
                   2*this->nu()*trace( deft*def ) );
    LOG(INFO) << "[add element stokes terms] (nu * (nabla u+ nabla^T u)  : nabla v) done in " << ti.elapsed() << "s\n";
    ti.restart();

    LOG(INFO) << "[add element stokes terms] ( p, div(v) ) + ( div(u), q )\n";
    *M_oplin +=
        integrate( elements( mesh ), _Q<2*( uOrder-1 )+3*( GeoOrder-1 )>(),
                   - div( v )*idt( p ) + divt( u )*id( q )
                 );
    LOG(INFO) << "[add element stokes terms] ( p, div(v) ) + ( div(u), q ) done in " << ti.elapsed()<< "s\n";
    ti.restart();


#if 1
    *M_oplin +=
        integrate( elements( mesh ), _Q<2*( uOrder-1 )+2*( GeoOrder-1 )>(),
                   this->epsPseudoCompressibility()*idt( p )*id( q )
                 );
    LOG(INFO) << "[add element stokes terms] ( epsilon p, q) ) done in " << ti.elapsed() << "s\n";
    ti.restart();

#endif
#if 0
    *M_oplin +=
        integrate( elements( mesh ), _Q<2*( uOrder-1 )+2*( GeoOrder-1 )>(),
                   idt( p )*id( lambda ) + id( p )*idt( lambda )
                 );
#endif
#if 0

    if ( this->deltaDivDiv() != 0.0 )
    {
        LOG(INFO) << "[add element divdiv term] (h delta div u,div v)\n";
        *M_oplin +=
            integrate( elements( mesh ), _Q<2*uOrder-2+3*( GeoOrder-1 )>(),
                       h()*this->deltaDivDiv()*div( v )*divt( u )
                     );
    }

#endif
    BOOST_FOREACH( std::string marker, this->dirichletVelocityMarkers() )
    {
        LOG(INFO) << "[add weakbc boundary terms velocity] boundary " << marker << " id : " << mesh->markerName( marker ) << "\n";
        *M_oplin +=
            integrate( markedfaces( mesh,mesh->markerName( marker ) ), _Q<( 2*uOrder-1 )+3*( GeoOrder-1 )>(),
                       -trans( SigmaNt )*id( v )
                       -trans( SigmaN )*idt( u )
                     );
        *M_oplin +=
            integrate( markedfaces( mesh,mesh->markerName( marker ) ), _Q<2*uOrder+2*( GeoOrder-1 )>(),
                       +this->gammaBc()*trans( idt( u ) )*id( v )/hFace()


                     );
        LOG(INFO) << "[Turek::initLinearOperators] oplin marked faces with marker " << marker << " integration done in " << ti.elapsed() << "s\n";
        ti.restart();
    }

#if 0
    BOOST_FOREACH( std::string marker, this->dirichletPressureMarkers() )
    {
        LOG(INFO) << "[add weakbc boundary terms pressure] boundary " << marker << " id : " << mesh->markerName( marker ) << "\n";
        *M_oplin +=
            integrate( markedfaces( mesh,mesh->markerName( marker ) ), _Q<( 2*uOrder-1 )+3*( GeoOrder-1 )>(),
                       -trans( -idt( p )*N() )*id( v )
                       -trans( -id( q )*N() )*idt( u )
                     );
        *M_oplin +=
            integrate( markedfaces( mesh,mesh->markerName( marker ) ), _Q<2*( uOrder-1 )+2*( GeoOrder-1 )>(),
                       +this->gammaBc()*idt( p )*id( q )/hFace()


                     );
        LOG(INFO) << "[Turek::initLinearOperators] oplin marked faces with marker " << marker << " integration done in " << ti.elapsed() << "s\n";
        ti.restart();
    }
#endif

    M_oplin->close();
    LOG(INFO) << "[Turek::initLinearOperators] oplin close in " << ti.elapsed() << "s\n";
    ti.restart();

    // stokes linear form
    *M_stokes_rhs =
        integrate( markedfaces( mesh,mesh->markerName( "inflow" ) ), _Q<3*uOrder+2*Dim+2*( GeoOrder-1 )>(),
                   trans( idf( this->inflow( time ) ) )*( -SigmaN+this->gammaBc()*id( v )/hFace() ) );
    LOG(INFO) << "[Turek::initLinearOperators] stokes  in " << ti.elapsed() << "s\n";
    ti.restart();
    M_stokes_rhs->close();
    LOG(INFO) << "[Turek::initLinearOperators] stokes close  in " << ti.elapsed() << "s\n";
    ti.restart();

    if ( this->vm().count( "export-matlab" ) )
        M_stokes_rhs->containerPtr()->printMatlab( "stokes_rhs.m" );

    LOG(INFO) << "[Turek::initLinearOperators] done in " << total_time.elapsed() << "\n";


}

template<int Dim, int Order, int GeoOrder>
void
Turek<Dim, Order, GeoOrder>::updateResidual( const vector_ptrtype& X, vector_ptrtype& R )
{

    boost::timer ti;
    LOG(INFO) << "[updateResidual] start\n";

    mesh_ptrtype mesh = M_Xh->mesh();

    fluid_element_type U( M_Xh, "U" );
    fluid_element_type V( M_Xh, "V" );
    fluid_element_0_type u = U.template element<0>();
    fluid_element_0_type v = V.template element<0>();
    fluid_element_1_type p = U.template element<1>();
    fluid_element_1_type q = V.template element<1>();

    fluid_element_0_type un = Un->template element<0>();
    fluid_element_0_type un1 = Un1->template element<0>();

    U = *X;

    AUTO( deft, 0.5*( gradt( u )+trans( gradt( u ) ) ) );
    AUTO( def, 0.5*( grad( v )+trans( grad( v ) ) ) );
    AUTO( Id, ( mat<Dim,Dim>( cst( 1 ), cst( 0 ), cst( 0 ), cst( 1. ) ) ) );
    AUTO( SigmaNt, ( -idt( p )*N()+2*this->nu()*deft*N() ) );
    AUTO( SigmaN, ( -id( p )*N()+2*this->nu()*def*N() ) );
    AUTO( beta, idv( u ) );
    // add the right hand side contribution from the non-homogeneous
    // Dirichlet contribution
    *M_residual =
        integrate( elements( mesh ), _Q<3*uOrder-1+3*( GeoOrder-1 )>(),
                   this->rho()* ( trans( idv( u ) )*id( v )*M_bdf->polyDerivCoefficient( 0 ) +
                                  trans( gradv( u )*idv( u ) )*id( v )
                                  - trans( idv( M_bdf->polyDeriv().template element<0>() ) ) *id( v ) ) ) +

        integrate( markedfaces( mesh,mesh->markerName( "inflow" ) ), _Q<3*uOrder+2*Dim+2*( GeoOrder-1 )>(),
                   //- this->gammaBc()*max(sqrt(trans(beta)*beta),this->nu()/hFace())*(trans(idf(this->inflow(time)))*N())*(trans(id(v))*N())
                   //-trans(vec( 4*this->Um()*Py()*(this->H()-Py())/math::pow(this->H(),2), constant(0.)))*( -SigmaN+this->gammaBc()*id(v)/hFace() ) );
                   - trans( idf( this->inflow( time ) ) )*( -SigmaN+this->gammaBc()*id( v )/hFace() ) );

    FsFunctionalLinear<fluid_functionspace_type> flin( M_Xh, M_backend );
    M_oplin->apply( U, flin );

    M_residual->add( flin );
    M_residual->close();
    *R = M_residual->container();

    LOG(INFO) << "[updateResidual] done in " << ti.elapsed() << "s\n";
}
template<int Dim, int Order, int GeoOrder>
void
Turek<Dim, Order, GeoOrder>::updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J )
{
    boost::timer ti;
    LOG(INFO) << "[updateJacobian] start\n";

    static bool is_init = false;

    mesh_ptrtype mesh = M_Xh->mesh();

    fluid_element_type U( M_Xh, "U" );
    fluid_element_type V( M_Xh, "V" );
    fluid_element_0_type u = U.template element<0>();
    fluid_element_0_type v = V.template element<0>();
    fluid_element_1_type p = U.template element<1>();
    fluid_element_1_type q = V.template element<1>();
    //fluid_element_2_type lambda = U.template element<2>();
    U = *X;

    if ( is_init == false )
    {
        *M_jac = integrate( elements( mesh ), _Q<3*uOrder-1+3*( GeoOrder-1 )>(),
                            this->rho()*trans( idt( u ) )*id( v )*M_bdf->polyDerivCoefficient( 0 ) +
                            this->rho()*trans( gradt( u )*idv( u ) )*id( v ) );
        *M_jac += integrate( elements( mesh ), _Q<3*uOrder-1+3*( GeoOrder-1 )>(),
                             +0*div( v )*idt( p )+0*divt( u )*id( q ) +0*idt( p )*id( q ) );
#if 0
        *M_jac += integrate( elements( mesh ), _Q<3*uOrder-1+3*( GeoOrder-1 )>(),
                             +0*idt( lambda )*id( p )+0*id( lambda )*idt( p ) );
#endif
        //this->rho()*trans(gradt(u)*idv(u))*id(v) );
        is_init = true;
    }

    else
    {
        M_jac->matPtr()->zero();
        *M_jac += integrate( elements( mesh ), _Q<3*uOrder-1+3*( GeoOrder-1 )>(),
                             this->rho()*trans( idt( u ) )*id( v )*M_bdf->polyDerivCoefficient( 0 ) +
                             this->rho()*trans( gradt( u )*idv( u ) )*id( v ) );
    }

#if 0
    AUTO( beta, idv( u ) );
    *M_jac += integrate( boundaryfaces( mesh ), _Q<4*uOrder+2*( GeoOrder-1 )>(),
                         this->gammaBc()*max( sqrt( trans( beta )*beta ),this->nu()/hFace() )*( trans( idt( u ) )*N() )*( trans( id( v ) )*N() ) );
#endif
    M_jac->close();

    if ( this->vm().count( "export-matlab" ) )
    {
        M_jac->matPtr()->printMatlab( "jac1.m" );
        M_oplin->matPtr()->printMatlab( "stokes1.m" );
    }

    M_jac->matPtr()->addMatrix( 1.0, M_oplin->mat() );
    J = M_jac->matPtr();

    if ( this->vm().count( "export-matlab" ) )
        J->printMatlab( "jac.m" );

    LOG(INFO) << "[updateJacobian] done in " << ti.elapsed() << "s\n";
}
template<int Dim, int Order, int GeoOrder>
void
Turek<Dim, Order, GeoOrder>::updateResidualJacobian( const vector_ptrtype& X, vector_ptrtype& R, sparse_matrix_ptrtype& J )
{
}

template<int Dim, int Order, int GeoOrder>
void
Turek<Dim, Order, GeoOrder>::run()
{
    time = 0.0;

    boost::timer ti;

    M_Xh->printInfo();


    LOG(INFO) << "[Turek] run() starts\n";
    using namespace Feel::vf;
    mesh_ptrtype mesh = M_Xh->mesh();

    fluid_element_type U( M_Xh, "U" );

    fluid_element_0_type u = U.template element<0>();
    fluid_element_1_type p = U.template element<1>();

    LOG(INFO) << "space+elements init done in " << ti.elapsed() << "s\n";
    ti.restart();

    //u = project( M_Xh, elements(mesh), constant(0.)*one() );
    u.zero();
    p.zero();

    vector_ptrtype S( M_backend->newVector( U.functionSpace() ) );
    vector_ptrtype R( M_backend->newVector( U.functionSpace() ) );
    sparse_matrix_ptrtype J;

    if ( ( this->init() == INIT_WITH_STOKES ) && ( M_bdf->timeInitial() == 0.0 ) )
    {
        this->solve( M_oplin->matPtr(), U, M_stokes_rhs->containerPtr() );
        this->normL2Div( U );
        exportResults( M_bdf->time(), U );
    }

    if ( M_bdf->timeInitial() > 0.0 )
    {
        U = M_bdf->unknown( 0 );
    }

    else
    {
        M_bdf->initialize( U );
    }


    boost::timer ttotal;
    LOG(INFO) << "[run] start iterating in time\n";


    for ( M_bdf->start(); M_bdf->isFinished() == false; M_bdf->next() )
    {
        LOG(INFO) << "============================================================\n";
        LOG(INFO) << "time: " << M_bdf->time() << "s, iteration: " << M_bdf->iteration() << "\n";

        *S = U;

        this->updateResidual( S, R );
        this->updateJacobian( S, J );

        //R->scale( -1 );
        nlsolve( J, U, R );
        exportResults( M_bdf->time(), U );

        M_bdf->shiftRight( U );

        this->normL2Div( U );
        LOG(INFO) << "time spent in iteration = " << M_bdf->realTimePerIteration() << "s\n";
    }

    LOG(INFO) << "total time spent :  " << ttotal.elapsed() << "s\n";
} // Turek::run



template<int Dim, int Order, int GeoOrder>
double
Turek<Dim, Order, GeoOrder>::normL2Div( fluid_element_type& U ) const
{
    fluid_element_0_type u = U.template element<0>();
    fluid_element_1_type p = U.template element<1>();
    mesh_ptrtype mesh = M_Xh->mesh();

    double int_one_times_N = integrate( markedfaces( u.functionSpace()->mesh(),mesh->markerName( "cylinder" ) ), _Q<2*( GeoOrder-1 )>(), trans( one() )*N() ).evaluate()( 0, 0 );
    LOG(INFO) << "             int 1.N = " << int_one_times_N << "\n";
    //
    // average pressure
    //
    double measure = integrate( elements( u.functionSpace()->mesh() ), _Q<2*( GeoOrder-1 )>(), constant( 1.0 ) ).evaluate()( 0, 0 );
    LOG(INFO) << "                    area = " << measure << "\n";
    double meanp = integrate( elements( u.functionSpace()->mesh() ), _Q<( uOrder-1 )+2*( GeoOrder-1 )>(), idv( p ) ).evaluate()( 0, 0 )/measure;
    LOG(INFO) << "                mean( p )= " << meanp << "\n";

    //
    // average pressure at outflow
    //
    double measure_outflow = integrate( markedfaces( u.functionSpace()->mesh(),mesh->markerName( "outflow" ) ), _Q<2*( GeoOrder-1 )>(),
                                        constant( 1.0 ) ).evaluate()( 0, 0 );
    LOG(INFO) << "         measure outflow = " << measure_outflow << "\n";
    double meanp_outflow = integrate( markedfaces( u.functionSpace()->mesh(),mesh->markerName( "outflow" ) ), _Q<( uOrder-1 )+2*( GeoOrder-1 )>(),
                                      idv( p ) ).evaluate()( 0, 0 )/measure_outflow;
    LOG(INFO) << "        mean_outflow( p )= " << meanp_outflow << "\n";


    //
    // Divergence
    //
    double intdiv2 = math::sqrt( integrate( elements( u.functionSpace()->mesh() ), _Q<2*( uOrder-1 )+4*( GeoOrder-1 )>(), divv( u )*divv( u ) ).evaluate()( 0, 0 ) );
    LOG(INFO) << "             |div(un)|_2 = " << intdiv2 << "\n";
    double intdiv = integrate( elements( u.functionSpace()->mesh() ), _Q<( uOrder-1 )+3*( GeoOrder-1 )>(), divv( u ) ).evaluate()( 0, 0 );
    LOG(INFO) << "            int( div(un))= " << intdiv << "\n";
    double intun_wall = integrate( markedfaces( u.functionSpace()->mesh(),mesh->markerName( "wall" ) ), _Q<uOrder+2*( GeoOrder-1 )>(), idv( u )*N() ).evaluate()( 0, 0 );
    LOG(INFO) << "          int(u.n)|_wall = " << intun_wall << "\n";
    double intun_inflow = integrate( markedfaces( u.functionSpace()->mesh(),mesh->markerName( "inflow" ) ), _Q<uOrder+2*( GeoOrder-1 )>(), idv( u )*N() ).evaluate()( 0, 0 );
    LOG(INFO) << "        int(u.n)|_inflow = " << intun_inflow << "\n";
    double intun_outflow = integrate( markedfaces( u.functionSpace()->mesh(),mesh->markerName( "outflow" ) ), _Q<uOrder+2*( GeoOrder-1 )>(), idv( u )*N() ).evaluate()( 0, 0 );
    LOG(INFO) << "       int(u.n)|_outflow = " << intun_outflow << "\n";
    double intun_cylinder = integrate( markedfaces( u.functionSpace()->mesh(),mesh->markerName( "cylinder" ) ), _Q<uOrder+2*( GeoOrder-1 )>(), idv( u )*N() ).evaluate()( 0, 0 );
    LOG(INFO) << "      int(u.n)|_cylinder = " << intun_cylinder << "\n";

    return intdiv;
}
template<int Dim, int Order, int GeoOrder>
void
Turek<Dim, Order, GeoOrder>::exportResults( double time, fluid_element_type& U )
{
    boost::timer total_ti;
    boost::timer ti;
    LOG(INFO) << "exportResults starts\n";

    auto mesh = M_Xh->mesh();
    auto u = U.template element<0>();
    auto p = U.template element<1>();

    double DeltaP = p( this->xa() )( 0, 0, 0 ) - p( this->xe() )( 0, 0, 0 );
    LOG(INFO) << "DeltaP=" << DeltaP << "\n";

    LOG(INFO) << "[exportResults] Dp : " << ti.elapsed() << "\n";
    ti.restart();

    auto defv = 0.5*( gradv( u )+trans( gradv( u ) ) );
    auto SigmaNv = ( -idv( p )*N()+2*this->nu()*defv*N() );

    auto Force = integrate( markedfaces( mesh,mesh->markerName( "cylinder" ) ), _Q<uOrder-1+3*( GeoOrder-1 )>(), SigmaNv ).evaluate();
    LOG(INFO) << "Force=" << Force << "\n";
    LOG(INFO) << "Scaling=" << this->scalingForce() << "\n";
    Force *= this->scalingForce();
    LOG(INFO) << "Force after scaling=" << Force << "\n";

    LOG(INFO) << "[exportResults] CD, CL : " << ti.elapsed() << "\n";

    M_data.precision( 8 );
    M_data.setf( std::ios::scientific );
    M_data << std::setw( 16 ) << time
           << std::setw( 16 ) << Force( 0, 0 )
           << std::setw( 16 ) << Force( 1, 0 )
           << std::setw( 16 ) << DeltaP
           << std::setw( 16 ) << u.linftyNorm() << std::endl;


    ti.restart();

    if ( this->doExport() != NO_EXPORT )
    {
        if ( this->doExport() == EXPORT_SAME_MESH )
            exporter->step( time )->setMesh( M_pressure_oplagp1->dualImageSpace()->mesh() );

        else if ( this->doExport() == EXPORT_LAGP1_MESH )
            exporter->step( time )->setMesh( M_velocity_oplagp1->dualImageSpace()->mesh() );

        else
        {
            LOG(INFO) << "invalid export strategy: using EXPORT_SAME_MESH\n";
            exporter->step( time )->setMesh( M_pressure_oplagp1->dualImageSpace()->mesh() );
        }

        LOG(INFO) << "[exportResults] setMesh : " << ti.elapsed() << "\n";
        ti.restart();

        exporter->step( time )->addScalar( "CD", Force( 0,0 ) );
        exporter->step( time )->addScalar( "CL", Force( 1,0 ) );
        exporter->step( time )->addScalar( "DP", DeltaP );
        exporter->step( time )->addScalar( "Uinf", u.linftyNorm() );


        velocity_element_type force( M_Xh->template functionSpace<0>(), "force" );
        force = vf::project( M_Xh->template functionSpace<0>(), boundaryfaces( mesh ), SigmaNv );
        exporter->step( time )->add( "total_force_per_area", force );
        force = vf::project( M_Xh->template functionSpace<0>(), boundaryfaces( mesh ), 2*this->nu()*defv*N() );
        exporter->step( time )->add( "viscous_force_per_area", force );


        //else
        //timeStep->setMesh( M_Xh->mesh() );
        exporter->step( time )->add( "pid",
                                     regionProcess( boost::shared_ptr<p0_space_type>( new p0_space_type( M_velocity_oplagp1->dualImageSpace()->mesh() ) ) ) );


        exporter->step( time )->add( "velocity", u );
        exporter->step( time )->add( "u", u.comp( X ) );
        exporter->step( time )->add( "v", u.comp( Y ) );

        if ( Dim == 3 )
            exporter->step( time )->add( "w", u.comp( Z ) );

        exporter->step( time )->add( "pressure", p );

        LOG(INFO) << "[exportResults] pid, U, u, v, p : " << ti.elapsed() << "\n";
        ti.restart();
        //
        pressure_element_type aux( M_Xh->template functionSpace<1>(), "mag" );
        aux = vf::project( M_Xh->template functionSpace<1>(), elements( mesh ), sqrt( trans( idv( u ) )*idv( u ) ) );

        exporter->step( time )->add( "velocity_mag", aux );

        LOG(INFO) << "[exportResults] ||U|| : " << ti.elapsed() << "\n";
        ti.restart();

        // cell Reynolds number
        aux = vf::project( M_Xh->template functionSpace<1>(), elements( mesh ), this->rho()*sqrt( trans( idv( u ) )*idv( u ) )*h()/this->nu() );

        exporter->step( time )->add( "cellRe", aux );

        LOG(INFO) << "[exportResults] cellRe : " << ti.elapsed() << "\n";
        ti.restart();

        if ( Dim == 3 )
        {
            // vorticity
            velocity_element_type vort( M_Xh->template functionSpace<0>(), "vort" );

            vector_ptrtype F_vort( M_backend->newVector( M_Xh->template functionSpace<0>() ) );
            form1( M_Xh->template functionSpace<0>(), F_vort ) = integrate( elements( mesh ),
                    _Q<uOrder+3*( GeoOrder-1 )>(),
                    trans( curlv( u ) )*id( vort ) );
            F_vort->close();

            this->solve( M_mass_v->matPtr(), vort, F_vort );

            exporter->step( time )->add( "vorticity", vort );
            exporter->step( time )->add( "vorticity_x", vort.comp( X ) );
            exporter->step( time )->add( "vorticity_y", vort.comp( Y ) );
            exporter->step( time )->add( "vorticity_z", vort.comp( Z ) );
        }

        else
        {
            vector_ptrtype F_vort( M_backend->newVector( M_Xh->template functionSpace<1>() ) );
            form1( M_Xh->template functionSpace<1>(), F_vort ) = integrate( elements( mesh ),
                    _Q<uOrder+3*( GeoOrder-1 )>(),
                    ( dxv( u.comp( Y ) )-dyv( u.comp( X ) ) )*id( aux ) );
            F_vort->close();

            this->solve( M_mass_s->matPtr(), aux, F_vort );

            exporter->step( time )->add( "vorticity", aux );
        }

        LOG(INFO) << "[exportResults] vorticity : " << ti.elapsed() << "\n";
        ti.restart();
        exporter->save();
    }


    LOG(INFO) << "[exportResults] save : " << ti.elapsed() << "\n";

    LOG(INFO) << "time spent in exportResults :  " << total_ti.elapsed() << "s\n";
} // Turek::export





}

#endif /* __TurekImpl_H */
