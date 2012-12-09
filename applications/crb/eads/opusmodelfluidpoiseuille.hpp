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
#ifndef _OPUSMODELFLUIDPOISEUILLE_HPP_
#define _OPUSMODELFLUIDPOISEUILLE_HPP_ 1

#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/foreach.hpp>

#include <feel/feel.hpp>

//#include <feel/options.hpp>
//#include <feel/feelalg/backend.hpp>
//#include <feel/feeldiscr/functionspace.hpp>
//#include <feel/feeldiscr/region.hpp>
//#include <feel/feeldiscr/operatorlagrangep1.hpp>
//#include <feel/feelpoly/im.hpp>
//#include <feel/feeldiscr/bdf2.hpp>
//#include <feel/feelvf/vf.hpp>

namespace Feel
{
/**
 * \addtogroup Models
 * \\@{
 */
Feel::po::options_description opusModelFluidOptions();
/**
 * \class OpusModelFluidPoiseuille
 * \brief Opus Poiseuille flow model
 *
 * This class is implementing a poiseuille flow for the Opus
 * benchmark.
 *
 * @author Christophe Prud'homme
 */
template<typename SpaceType>
class OpusModelFluidPoiseuille : public OpusModelBase
{
    typedef OpusModelBase super;
public:
#define Entity Simplex
    /**
     * Typedefs  and Constants
     */
    static const uint16_type Dim = SpaceType::nDim;
    static const uint16_type Order = SpaceType::basis_0_type::nOrder;

    static const uint16_type imOrder = 2*Order;
    static const uint16_type GeoOrder = 1;
    typedef OpusModelFluidPoiseuille<SpaceType> self_type;

    typedef double value_type;

    typedef SpaceType functionspace_type;
    typedef boost::shared_ptr<SpaceType> functionspace_ptrtype;

    typedef typename functionspace_type::element_type element_type;
    typedef typename element_type::template sub_element<0>::type element_0_type;
    typedef typename element_type::template sub_element<1>::type element_1_type;

    /**
     * constructor: Xh space and some space functions are initialized
     */
    OpusModelFluidPoiseuille( po::variables_map const& vm, functionspace_ptrtype const& Xh );

    ~OpusModelFluidPoiseuille() {}

    void update( double time );

    void solve( element_type& T );

    void setFluidFlowRate( double r )
    {
        M_flow_rate = r;
    }
private:

    po::variables_map vm;
    value_type M_flow_rate;
    value_type M_time;

    functionspace_ptrtype M_Xh;
}; // OpusModelFluidPoiseuille

template<typename SpaceType>
OpusModelFluidPoiseuille<SpaceType>::OpusModelFluidPoiseuille(  po::variables_map const& _vm,
        functionspace_ptrtype const& Xh )
    :
    super(),
    vm( _vm ),
    M_flow_rate( vm["fluid.flow-rate"].template as<double>() ), // m.s^-1
    M_time( 0. ),
    M_Xh( Xh )
{
    LOG(INFO) << "[OpusModelFluidPoiseuille] constructor starts\n";
    FEELPP_ASSERT( M_Xh != 0 ).error( "[OpusModelFluidPoiseuille] invalid functionspace_ptrtype" );
    LOG(INFO) << "[OpusModelFluidPoiseuille] constructor done\n";
}


template<typename SpaceType>
void
OpusModelFluidPoiseuille<SpaceType>::update( double time )
{
    M_time = time;

}

template<typename SpaceType>
void
OpusModelFluidPoiseuille<SpaceType>::solve( element_type& U )
{
    using namespace vf;
    M_flow_rate = this->data()->component( "AIR" ).flowRate();
    double e_AIR = this->data()->component( "AIR" ).e();
    double e_PCB = this->data()->component( "PCB" ).e();
    double e_IC = this->data()->component( "IC1" ).e();
    //double L_IC = this->data()->component("IC1").h();

    LOG(INFO) << "flow_rate = " << M_flow_rate << "\n";
    LOG(INFO) << "e_AIR = " << e_AIR << "\n";
    LOG(INFO) << "e_PCB = " << e_PCB << "\n";
    LOG(INFO) << "e_IC = " << e_IC << "\n";

    AUTO( chi_AIR, chi( Px() >= e_PCB+e_IC ) );
    AUTO( ft, ( constant( 1.0-math::exp( -M_time/3.0 ) ) ) );
    //AUTO( ft, (constant(1.0)) );
    AUTO( vy, ( constant( 3. )/( 2.*( e_AIR-e_IC ) ) )*M_flow_rate*( 1.-vf::pow( ( Px()-( ( e_AIR+e_IC )/2+e_PCB ) )/( ( e_AIR-e_IC )/2 ),2 ) ) );
    //double x_mid = e_PCB+(e_IC+e_AIR)/2;
    //AUTO( vy, (constant(3)/(2*e_AIR))*M_flow_rate*(1-vf::pow((Px()-(x_mid))/(e_AIR/2),2))*ft*chi_AIR );

    element_0_type u = U.template element<0>();
    u = vf::project( M_Xh->template functionSpace<0>(),
                     markedelements( M_Xh->mesh(),M_Xh->mesh()->markerName( "AIR4" ) ),
                     vec( constant( 0. ), vy*ft*chi_AIR ) );

    double intu = integrate( markedfaces( M_Xh->mesh(),
                                          M_Xh->mesh()->markerName( "Gamma_4_AIR4" ) ),
                             -trans( idv( u ) )*N() ).evaluate()( 0,0 );
    double flowr = intu;
    LOG(INFO) << "[poiseuille] D=" << flowr << " umax = " << u.linftyNorm() << "\n";
}
/** \\@} */
} // Feel

#endif



