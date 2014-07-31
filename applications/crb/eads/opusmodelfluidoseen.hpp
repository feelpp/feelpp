/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-11-15

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file opusmodelfluid.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-11-15
 */
#ifndef _OPUSMODELFLUIDOSEEN_HPP_
#define _OPUSMODELFLUIDOSEEN_HPP_ 1

#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/foreach.hpp>

#include <feel/options.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feelts/bdf.hpp>

namespace Feel
{
/**
 * \addtogroup Models
 * \\@{
 */
Feel::po::options_description opusModelFluidOseenOptions();

/**
 * \class OpusModelFluidOseen
 * \brief class for the fluid Navier-Stokes model
 *
 * This class implements the fully nonlinear incompressible
 * Navier-Stokes equations.
 *
 * @author Christophe Prud'homme
 */
template<typename SpaceType>
class OpusModelFluidOseen : public OpusModelBase
{
    typedef OpusModelBase super;
public:
#define Entity Simplex
    /**
     * Typedefs  and Constants
     */
    static const uint16_type Dim = SpaceType::nDim;
    static const uint16_type Order = SpaceType::template Basis<0>::type::nOrder;
    static const uint16_type uOrder = SpaceType::template Basis<0>::type::nOrder;
    static const uint16_type pOrder = SpaceType::template Basis<1>::type::nOrder;

    static const uint16_type imOrder = 2*Order;
    static const uint16_type GeoOrder = 1;
    typedef OpusModelFluidOseen<SpaceType> self_type;

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    typedef typename SpaceType::mesh_type mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef SpaceType fluid_functionspace_type;
    typedef boost::shared_ptr<SpaceType> fluid_functionspace_ptrtype;

    typedef typename fluid_functionspace_type::element_type fluid_element_type;
    typedef typename fluid_element_type::template sub_element<0>::type fluid_element_0_type;
    typedef typename fluid_element_type::template sub_element<1>::type fluid_element_1_type;
    typedef typename fluid_element_type::template sub_element<0>::type velocity_element_type;
    typedef typename fluid_element_type::template sub_element<1>::type pressure_element_type;

    typedef typename fluid_functionspace_type::template sub_functionspace<0>::type velocity_functionspace_type;
    typedef typename fluid_functionspace_type::template sub_functionspace<1>::type pressure_functionspace_type;


    /* Operators */
    typedef OperatorLinear<fluid_functionspace_type, fluid_functionspace_type> oplin_type;
    typedef boost::shared_ptr<oplin_type> oplin_ptrtype;
    typedef OperatorLinear<velocity_functionspace_type, velocity_functionspace_type> op_vector_type;
    typedef boost::shared_ptr<op_vector_type> op_vector_ptrtype;
    typedef OperatorLinear<pressure_functionspace_type, pressure_functionspace_type> op_scalar_type;
    typedef boost::shared_ptr<op_scalar_type> op_scalar_ptrtype;
    typedef FsFunctionalLinear<fluid_functionspace_type> funlin_type;
    typedef boost::shared_ptr<funlin_type> funlin_ptrtype;


    /**
     * constructor: Xh space and some space functions are initialized
     */
    OpusModelFluidOseen( po::variables_map const& vm, fluid_functionspace_ptrtype const& Xh );

    ~OpusModelFluidOseen() {}

    // return air density (kg/m^3)
    double rho() const
    {
        return 1.204;
    }
    //! \return air viscosity (kg/(m·s))
    double nu() const
    {
        return 1.78* 1e-5;
    }
    void setFluidFlowRate( double r )
    {
        M_flow_rate = r;
    }
#if 0
    template< typename MassExpr,
              typename DiffExpr,
              typename ConvExpr,
              typename RhsExpr >
    void update( MassExpr const& mass, DiffExpr const& diff, ConvExpr const& conv, RhsExpr const& rhsExpr );
#else
    void update( double time );
#endif

    void solve( fluid_element_type& U );

private:

    void initLinearOperators();
    void updateResidual( const vector_ptrtype& X, vector_ptrtype& R );
    void updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J );
    void updateResidualJacobian( const vector_ptrtype& X, vector_ptrtype& R, sparse_matrix_ptrtype& J );
private:

    backend_ptrtype M_backend;

    double M_time;
    double M_flow_rate, M_current_flow_rate;

    sparse_matrix_ptrtype M_D;
    vector_ptrtype M_F;

    fluid_functionspace_ptrtype M_Xh;

    oplin_ptrtype M_oplin;
    oplin_ptrtype M_jac;

    funlin_ptrtype M_residual;

    op_vector_ptrtype M_mass_v;
    op_scalar_ptrtype M_mass_s;
}; // OpusModelFluidOseen

template<typename SpaceType>
OpusModelFluidOseen<SpaceType>::OpusModelFluidOseen( po::variables_map const& vm, fluid_functionspace_ptrtype const& Xh )
    :
    super( vm ),
    M_time( 100 ),
    M_flow_rate( vm["fluid.flow-rate"].as<double>()  ),
    M_current_flow_rate( vm["fluid.flow-rate"].as<double>()  ),
    M_backend( backend_type::build( vm, "fluid" ) ),
    M_Xh( Xh ),
    M_D(),
    M_F()
{
    LOG(INFO) << "[OpusModelFluidOseen] flow rate = " << M_flow_rate << "\n";
    FEELPP_ASSERT( M_backend != 0 ).error( "[OpusModelFluidOseen] invalid backend" );
    FEELPP_ASSERT( M_Xh != 0 ).error( "[OpusModelFluidOseen] invalid functionspace_ptrtype" );

    M_D = M_backend->newMatrix( M_Xh, M_Xh );
    FEELPP_ASSERT( M_D != 0 ).error( "invalid matrix" );

    M_F = M_backend->newVector( M_Xh );
    FEELPP_ASSERT( M_F != 0 ).error( "invalid vector" );

    M_backend->nlSolver()->residual = boost::bind( &self_type::updateResidual, boost::ref( *this ), _1, _2 );
    M_backend->nlSolver()->jacobian = boost::bind( &self_type::updateJacobian, boost::ref( *this ), _1, _2 );

    // jacobian and residual
    M_jac = oplin_ptrtype( new oplin_type( M_Xh, M_Xh, M_backend ) );
    M_residual = funlin_ptrtype( new funlin_type( M_Xh, M_backend ) );

    this->initLinearOperators();
}

template<typename SpaceType>
void
OpusModelFluidOseen<SpaceType>::initLinearOperators()
{
    using namespace vf;
    boost::timer ti, total_time;
    LOG(INFO) << "[OpusModelFluidOseen::initLinearOperators] start\n";
    using namespace Feel::vf;
    mesh_ptrtype mesh = M_Xh->mesh();

    fluid_element_type U( M_Xh, "U" );
    fluid_element_type V( M_Xh, "V" );

    fluid_element_0_type u = U.template element<0>();
    fluid_element_1_type p = U.template element<1>();
    //fluid_element_2_type lambda = U.template element<2>();

    fluid_element_0_type v = V.template element<0>();
    fluid_element_1_type q = V.template element<1>();

    LOG(INFO) << "[OpusModelFluidOseen::initLinearOperators] space+elements init done in " << ti.elapsed() << "s\n";
    ti.restart();
    std::cout << "perimeter:"  << integrate( _range=boundaryfaces( mesh ), _expr=cst( 1.0 ) ).evaluate()( 0, 0 ) << std::endl;
    //M_mass_v = op_vector_ptrtype( new op_vector_type( M_Xh->template functionSpace<0>(), M_Xh->template functionSpace<0>(), M_backend ) );
    //M_mass_s = op_scalar_ptrtype( new op_scalar_type( M_Xh->template functionSpace<1>(), M_Xh->template functionSpace<1>(), M_backend ) );
    M_oplin = oplin_ptrtype( new oplin_type( M_Xh, M_Xh, M_backend ) );
    auto deft=  sym( gradt( u ) );
    auto def = sym( grad( u ) );
    auto defv = sym( gradv( u ) );
    auto Id = ( mat<Dim,Dim>( cst( 1 ), cst( 0 ), cst( 0 ), cst( 1. ) ) );
    auto Sigmat = -idt( p )*Id + 2*this->nu()*deft;
    auto SigmaNt = ( -idt( p )*N()+2*this->nu()*deft*N() );
    auto SigmaN = ( -id( p )*N()+2*this->nu()*def*N() );
    auto SigmaNv = ( -idv( p )*N()+2*this->nu()*defv*N() );

    //*M_mass_v = integrate( elements(mesh), _Q<2*uOrder+2*(GeoOrder-1)>(), trans(idt(u))*id(v) );
    //M_mass_v->close();
    //*M_mass_s = integrate( elements(mesh), _Q<2*uOrder+2*(GeoOrder-1)>(), trans(idt(p))*id(q) );
    //M_mass_s->close();

    // oplin
    LOG(INFO) << "[OpusModelFluidOseen::add element stokes terms] (nu * (nabla u+ nabla^T u)  : nabla v)\n";
    *M_oplin =
        integrate( _range=elements( mesh ), _expr=trace( Sigmat*trans( grad( v ) ) ) );
    LOG(INFO) << "[OpusModelFluidOseen::add element stokes terms] (nu * (nabla u+ nabla^T u)  : nabla v) done in " << ti.elapsed() << "s\n";
    ti.restart();

    LOG(INFO) << "[OpusModelFluidOseen::add element stokes terms] ( p, div(v) ) + ( div(u), q )\n";
    *M_oplin += integrate( _range=elements( mesh ), _expr=divt( u )*id( q ) );
    LOG(INFO) << "[OpusModelFluidOseen::add element stokes terms] ( p, div(v) ) + ( div(u), q ) done in " << ti.elapsed()<< "s\n";
    ti.restart();

#if 1
    *M_oplin +=
        integrate( _range=elements( mesh ), _expr=this->data()->epsPseudoCompressibility()*idt( p )*id( q ) );
    LOG(INFO) << "[OpusModelFluidOseen::add element stokes terms] ( epsilon p, q) ) done in " << ti.elapsed() << "s\n";
#endif
    ti.restart();

    BOOST_FOREACH( std::string marker, this->data()->dirichletVelocityMarkers() )
    {
        std::cout << "  -- dirichlet marker: "  << marker << "\n";
        std::cout << "  -- dirichlet perimeter:"  << integrate( _range=markedfaces( mesh,marker ), _expr=cst( 1.0 ) ).evaluate()( 0, 0 ) << std::endl;
    }
    BOOST_FOREACH( std::string marker, this->data()->dirichletVelocityMarkers() )
    {
        LOG(INFO) << "[OpusModelFluidOseen::add weakbc boundary terms velocity] boundary " << marker << " id : " << mesh->markerName( marker ) << "\n";
        LOG(INFO) << "[OpusModelFluidOseen::add weakbc boundary terms velocity] " << mesh->markerName( marker )
              << " : nelts: " << std::distance( markedfaces( mesh,marker ).template get<1>(),
                                                markedfaces( mesh,marker ).template get<2>() ) << "\n";
        std::cout << "[OpusModelFluidOseen::add weakbc boundary terms velocity] " << mesh->markerName( marker )
                  << " : nelts: " << std::distance( markedfaces( mesh,marker ).template get<1>(),
                          markedfaces( mesh,marker ).template get<2>() ) << "\n";
        std::cerr << " -- bdy " << marker << " terms 1" << std::endl;
        *M_oplin +=
            integrate( _range=markedfaces( mesh, marker ),
                       _expr=-trans( SigmaNt )*id( v )-trans( SigmaN )*idt( u ) );
        std::cerr << " -- bdy " << marker << " terms 1 done" << std::endl;
        std::cerr << " -- bdy " << marker << " terms 2" << std::endl;
        *M_oplin +=
            integrate( _range=markedfaces( mesh, marker ),
                       _expr=this->data()->gammaBc()*trans( idt( u ) )*id( v )/hFace() );
        std::cerr << " -- bdy " << marker << " terms 2 done" << std::endl;
        LOG(INFO) << "[OpusModelFluidOseen::initLinearOperators] oplin marked faces with marker " << marker << " integration done in " << ti.elapsed() << "s\n";
        ti.restart();
    }

    M_oplin->close();
    LOG(INFO) << "[OpusModelFluidOseen::initLinearOperators] oplin close in " << ti.elapsed() << "s\n";
    ti.restart();

}

template<typename SpaceType>
void
OpusModelFluidOseen<SpaceType>::updateResidual( const vector_ptrtype& X,
        vector_ptrtype& R )
{
    using namespace vf;
    boost::timer ti;
    LOG(INFO) << "[updateResidual] start\n";

    mesh_ptrtype mesh = M_Xh->mesh();

    fluid_element_type U( M_Xh, "U" );
    fluid_element_type V( M_Xh, "V" );
    fluid_element_0_type u = U.template element<0>();
    fluid_element_0_type v = V.template element<0>();
    fluid_element_1_type p = U.template element<1>();
    fluid_element_1_type q = V.template element<1>();

    //fluid_element_0_type un = Un->template element<0>();
    //fluid_element_0_type un1 = Un1->template element<0>();

    U = *X;

    auto deft = sym( gradt( u ) );
    auto def = sym( grad( u ) );
    auto Id = ( mat<Dim,Dim>( cst( 1 ), cst( 0 ), cst( 0 ), cst( 1. ) ) );
    auto SigmaNt = ( -idt( p )*N()+2*this->nu()*deft*N() );
    auto SigmaN = ( -id( p )*N()+2*this->nu()*def*N() );
    auto beta= idv( u );
    // add the right hand side contribution from the non-homogeneous
    // Dirichlet contribution
    *M_residual =
        integrate( _range=elements( mesh ),
                   _expr=this->rho()* ( //trans(idv(u))*id(v)*M_bdf->polyDerivCoefficient(0) +
                             trans( gradv( u )*idv( u ) )*id( v )
                             //- trans(idv( M_bdf->polyDeriv().template element<0>() ) ) *id(v)
                         ) );

    double e_AIR = this->data()->component( "AIR" ).e();
    double e_PCB = this->data()->component( "PCB" ).e();
    double e_IC = this->data()->component( "IC1" ).e();
    double L_IC = this->data()->component( "IC1" ).h();
    LOG(INFO) << "[OpusModelFluidOseen] e_AIR = " << e_AIR << "\n";
    LOG(INFO) << "[OpusModelFluidOseen] e_PCB = " << e_PCB << "\n";
    LOG(INFO) << "[OpusModelFluidOseen] e_IC = " << e_IC << "\n";
    LOG(INFO) << "[OpusModelFluidOseen] L_IC = " << L_IC << "\n";
    LOG(INFO) << "[OpusModelFluidOseen] flow rate = " << M_current_flow_rate << "\n";
    LOG(INFO) << "[OpusModelFluidOseen] gamm bc = " << this->data()->gammaBc() << "\n";

    auto ft = ( constant( 1.0 )-vf::exp( -M_time ) );
    auto vy = ( constant( 3. )/( 2.*( e_AIR ) ) )*M_current_flow_rate*( 1.0-vf::pow( ( Px()-( ( e_AIR )/2.+e_PCB ) )/( ( e_AIR )/2. ),2. ) )*ft;

    *M_residual +=
        integrate( markedfaces( mesh, "Gamma_4_AIR" ),
                   //- this->gammaBc()*max(sqrt(trans(beta)*beta),this->nu()/hFace())*(trans(idf(this->inflow(time)))*N())*(trans(id(v))*N())
                   //-trans(vec( 4*this->Um()*Py()*(this->H()-Py())/math::pow(this->H(),2), constant(0.)))*( -SigmaN+this->gammaBc()*id(v)/hFace() ) );
                   - trans( vec( constant( 0. ),vy ) )*( -SigmaN+this->data()->gammaBc()*id( v )/hFace() ) );
    //trans(vec(constant(0),(constant(3)/(2*(e_AIR)))*M_current_flow_rate*(1-vf::pow((Px()-((e_AIR)/2+e_PCB))/((e_AIR)/2),2))))*( -SigmaN+this->gammaBc()*id(v)/hFace() ) );

#if 0
    *M_residual +=
        integrate( markedfaces( mesh,mesh->markerName( "Gamma_3_AIR4" ) ),
                   trans( 101.56*1e3 * N() ) * id( v ) );
#endif

    FsFunctionalLinear<fluid_functionspace_type> flin( M_Xh, M_backend );
    M_oplin->apply( U, flin );

    M_residual->add( flin );
    M_residual->close();
    *R = M_residual->container();

    if ( this->vm().count( "export-matlab" ) )
        R->printMatlab( "R.m" );

    LOG(INFO) << "[updateResidual] done in " << ti.elapsed() << "s\n";
}

template<typename SpaceType>
void
OpusModelFluidOseen<SpaceType>::updateJacobian( const vector_ptrtype& X,
        sparse_matrix_ptrtype& J )
{
    using namespace vf;
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

    auto convectionTerms = ( trans( val( gradv( u )*this->rho() )*idt( u ) ) + trans( gradt( u )*val( idv( u )*this->rho() ) ) )*id( v );

    if ( is_init == false )
    {
        *M_jac = integrate( elements( mesh ),
                            //this->rho()*trans(idt(u))*id(v)*M_bdf->polyDerivCoefficient( 0 ) +
                            convectionTerms );
        //0*trans(gradt(u)*idv(u))*id(v) );
        *M_jac += integrate( elements( mesh ),
                             +0*div( v )*idt( p )+0*divt( u )*id( q ) +0*idt( p )*id( q ) );
#if 0
        *M_jac += integrate( elements( mesh ),
                             +0*idt( lambda )*id( p )+0*id( lambda )*idt( p ) );
#endif
        //this->rho()*trans(gradt(u)*idv(u))*id(v) );
        is_init = true;
    }

    else
    {
        M_jac->matPtr()->zero();
        *M_jac += integrate( elements( mesh ),
                             //this->rho()*trans(idt(u))*id(v)*M_bdf->polyDerivCoefficient( 0 ) +
                             convectionTerms );
        //0*trans(gradt(u)*idv(u))*id(v) );
    }

#if 0
    AUTO( beta, idv( u ) );
    *M_jac += integrate( boundaryfaces( mesh ),
                         this->gammaBc()*max( sqrt( trans( beta )*beta ),this->nu()/hFace() )*( trans( idt( u ) )*N() )*( trans( id( v ) )*N() ) );
#endif
    M_jac->close();

    M_jac->matPtr()->addMatrix( 1.0, M_oplin->mat() );
    J = M_jac->matPtr();

    LOG(INFO) << "[updateJacobian] done in " << ti.elapsed() << "s\n";

}

template<typename SpaceType>
void
OpusModelFluidOseen<SpaceType>::updateResidualJacobian( const vector_ptrtype& X, vector_ptrtype& R,
        sparse_matrix_ptrtype& J )
{
}

#if 0
template<typename SpaceType>
template< typename MassExpr,
          typename DiffExpr,
          typename ConvExpr,
          typename RhsExpr >
void
OpusModelFluidOseen<SpaceType>::update( MassExpr const& mass_coeff, // rank-0
                                        DiffExpr const& diff_coeff, // rank-0
                                        ConvExpr const& conv_coeff, // rank-1 (column vector)
                                        RhsExpr const& rhs_coeff    // rank-0
                                      )
{
}
#else
template<typename SpaceType>
void
OpusModelFluidOseen<SpaceType>::update( double time )
{
    M_time = time;
}
#endif

template<typename SpaceType>
void
OpusModelFluidOseen<SpaceType>::solve( fluid_element_type& U )
{
    int Ncont=10+std::ceil( std::log10( M_flow_rate ) );

    for ( int i = 0; i < Ncont; ++i )
    {
        int denom = ( Ncont==1 )?1:Ncont-1;
        M_current_flow_rate = math::exp( math::log( 1e-4 )+i*( math::log( M_flow_rate )-math::log( 1e-4 ) )/denom );
        LOG(INFO) << "[oseen] Computing fluid flow for flow rate:" << M_current_flow_rate << " ("  << i << ") / " << M_flow_rate << "("  << Ncont << ")\n";

        auto R = M_backend->newVector( U.functionSpace() );
        auto J = M_backend->newMatrix( U.functionSpace(), U.functionSpace() );
        M_backend->nlSolve( _jacobian=J, _solution=U, _residual=R );
        using namespace Feel::vf;
        double inflow = integrate( markedfaces( M_Xh->mesh(),"Gamma_4_AIR" ),
                                   -trans( idv( U.template element<0>() ) )*N() ).evaluate()( 0,0 );
        double outflow = integrate( markedfaces( M_Xh->mesh(),"Gamma_3_AIR" ),
                                    trans( idv( U.template element<0>() ) )*N() ).evaluate()( 0,0 );
        double wallflow = integrate( markedfaces( M_Xh->mesh(),"Gamma_1" ),
                                     trans( idv( U.template element<0>() ) )*N() ).evaluate()( 0,0 );
        wallflow += integrate( markedfaces( M_Xh->mesh(),"Gamma_2" ),
                               trans( idv( U.template element<0>() ) )*N() ).evaluate()( 0,0 );
        LOG(INFO) << "[oseen] inflow=" << inflow <<  " outflow=" << outflow << "  wallflow ="  << wallflow << " umax = " << U.template element<0>().linftyNorm() << "\n";
    }
}

/** \\@} */

} // Feel

#endif



