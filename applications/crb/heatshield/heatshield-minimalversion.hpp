/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-05-11

  Copyright (C) 2009-2014 Feel++ Consortium

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
   \file heatshield-minimalversion.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-05-11
 */
#ifndef FEELPP_HEATSHIELDMINIMALVERSION_HPP
#define FEELPP_HEATSHIELDMINIMALVERSION_HPP 1

#include <boost/timer.hpp>
#include <feel/options.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/pch.hpp>

#include <feel/feelfilters/gmsh.hpp>

#include <feel/feelcrb/modelcrbbase.hpp>

#include <feel/feelts/bdf.hpp>


namespace Feel
{

po::options_description
makeHeatShieldMinimalVersionOptions()
{
    return  bdf_options( "heatshield" );
}
AboutData
makeHeatShieldMinimalVersionAbout( std::string const& str = "heatShield" )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "heat shield Benchmark (minimal code version)",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2009-2014 Feel++ Consortium");
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Stephane Veys", "developer", "", "" );
    return about;
}

/**
 * \class HeatShieldMinimalVersion
 * \brief brief description
 *
 * @author Christophe Prud'homme
 * @see
 */
template<int Order>
class HeatShieldMinimalVersion : public ModelCrbBase< ParameterSpace<2>, decltype(Pch<Order>(Mesh<Simplex<2>>::New())) , TimeDependent >
{
public:

    typedef ModelCrbBase< ParameterSpace<2>, decltype(Pch<Order>(Mesh<Simplex<2>>::New())) , TimeDependent > super_type;
    typedef typename super_type::element_type element_type;
    typedef typename super_type::parameter_type parameter_type;
    typedef typename super_type::bdf_ptrtype bdf_ptrtype;
    //! initialization of the model
    void initModel();

    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    double output( int output_index, parameter_type const& mu, element_type &T, bool need_to_solve=false, bool export_outputs=false );

    bdf_ptrtype bdfModel(){ return M_bdf; }

private:
    double surface;
    bdf_ptrtype M_bdf;
};


template<int Order>
void HeatShieldMinimalVersion<Order>::initModel()
{

    CHECK( this->is_linear && this->is_time_dependent ) << "Invalid model is_linear:" << this->is_linear << " is_time_dependent:" << this->is_time_dependent << "\n";
    LOG_IF( WARNING, ((this->Options&NonLinear) == NonLinear) ) << "Invalid model is_linear:" << this->is_linear << " is_time_dependent:" << this->is_time_dependent << "\n";

    this->setFunctionSpaces( Pch<Order>( loadMesh( _mesh=new Mesh<Simplex<2>> ) ) );

    auto mesh = this->Xh->mesh();

    if( Environment::worldComm().isMasterRank() )
    {
        std::cout << "Number of local dof " << this->Xh->nLocalDof() << "\n";
        std::cout << "Number of dof " << this->Xh->nDof() << "\n";
    }

    surface = integrate( _range=elements( mesh ), _expr=cst( 1. ) ).evaluate()( 0,0 );
    double holes = integrate( _range=markedfaces( mesh , "gamma_holes" ), _expr=cst( 1. ) ).evaluate()( 0,0 );
    std::cout<<"holes : "<<holes<<" and surface : "<<surface<<std::endl;

    M_bdf = bdf( _space=this->Xh, _vm=Environment::vm(), _name="heatshield" , _prefix="heatshield" );

    auto mu_min = this->Dmu->element();
    mu_min <<  /* Bi_out */ 1e-2 , /*Bi_in*/1e-3;
    this->Dmu->setMin( mu_min );
    auto mu_max = this->Dmu->element();
    mu_max << /* Bi_out*/0.5   ,  /*Bi_in*/0.1;
    this->Dmu->setMax( mu_max );

    auto u = this->Xh->element();
    auto v = this->Xh->element();

    LOG(INFO) << "Number of dof " << this->Xh->nLocalDof() << "\n";

    //lhs
    auto a0 = form2( _trial=this->Xh, _test=this->Xh);
    a0 = integrate( _range= elements( mesh ), _expr= gradt( u )*trans( grad( v ) ) );
    this->addLhs( { a0 , "1" } );
    auto a1 = form2( _trial=this->Xh, _test=this->Xh);
    a1 = integrate( _range= markedfaces( mesh, "left" ), _expr= idt( u )*id( v ) );
    this->addLhs( { a1 , "mu0" } );
    auto a2 = form2( _trial=this->Xh, _test=this->Xh);
    a2 = integrate( _range= markedfaces( mesh, "gamma_holes" ), _expr= idt( u )*id( v ) );
    this->addLhs( { a2 , "mu1" } );
    //mass matrix
    auto mass = form2( _trial=this->Xh, _test=this->Xh);
    mass = integrate ( _range=elements( mesh ), _expr=idt( u )*id( v ) );
    this->addMass( { mass , "1" } );
    //rhs
    auto f0 = form1( _test=this->Xh );
    f0 = integrate( _range=markedfaces( mesh,"left" ), _expr= id( v ) ) ;
    this->addRhs( { f0, "mu0" } );
    //output
    auto out = form1( _test=this->Xh );
    out = integrate( _range=elements( mesh ), _expr= (1./surface)*id( v ) ) ;
    this->addOutput( { out, "1" } );

    auto energy = form2( _trial=this->Xh, _test=this->Xh);
    energy =
        integrate( elements( mesh ), gradt( u )*trans( grad( v ) ) ) +
        integrate(  markedfaces( mesh, "left" ), 0.01 * idt( u )*id( v ) ) +
        integrate(  markedfaces( mesh, "gamma_holes" ), 0.001 * idt( u )*id( v ) )
        ;
    this->addEnergyMatrix( energy );
    this->addMassMatrix(mass);

} // HeatShield::init


template<int Order>
double HeatShieldMinimalVersion<Order>::output( int output_index, parameter_type const& mu, element_type &u, bool need_to_solve , bool export_outputs )
{
    auto mesh = this->Xh->mesh();
    double output=0;
    if ( output_index==0 )
        output = integrate( _range=markedfaces( mesh,"left" ), _expr= mu(0)*idv( u ) ).evaluate()(0,0);
    if( output_index==1 )
        output = integrate( _range=elements( mesh ), _expr= (1./surface)*idv( u ) ).evaluate()(0,0);
    if( output_index > 1 )
        throw std::logic_error( "[HeatShield::output] error with output_index : only 0 or 1 " );

    return output ;
}

}

#endif /* __HeatShield_H */
