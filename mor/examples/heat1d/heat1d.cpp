//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 14 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!
#include <feel/feelmor/crbplugin.hpp>
#include <heat1d.hpp>


namespace Feel {

void
Heat1D::initBetaQ()
{
    this->M_betaAq.resize( 3 );
    this->M_betaFq.resize( 2 );
    this->M_betaFq[0].resize( 2 );
    this->M_betaFq[1].resize( 1 );
}

Heat1D::super_type::betaq_type
Heat1D::computeBetaQ( parameter_type const& mu )
{
    //std::cout << "computeBetaQ start \n" << mu << std::endl;
    if ( this->M_betaAq.empty() )
        this->initBetaQ();

    this->M_betaAq = {1, mu(0), mu(1) };
    this->M_betaFq[0] = { mu( 2 ), mu(3) };
    this->M_betaFq[1] = { 1. };
    //std::cout << "computeBetaQ finish \n";
    return boost::make_tuple( this->M_betaAq, this->M_betaFq );
}

void
Heat1D::initModel()
{

    CHECK( is_linear && !is_time_dependent ) << "Invalid model is_linear:" << is_linear << " is_time_dependent:" << is_time_dependent << "\n";
    LOG_IF( WARNING, ((Options&NonLinear) == NonLinear) ) << "Invalid model is_linear:" << is_linear << " is_time_dependent:" << is_time_dependent << "\n";
    LOG_IF( WARNING, ((Options&TimeDependent) == TimeDependent) ) << "Invalid model is_linear:" << is_linear << " is_time_dependent:" << is_time_dependent << "\n";
    /*
     * First we create the mesh
     */
    this->setFunctionSpaces( Pch<5>( loadMesh( _mesh=new Mesh<Simplex<1>> ) ) );

    Dmu->setDimension( 4 );
    //static const int N = 2;
    auto mu_min = Dmu->element();
    mu_min << 0.2, 0.2, 0.01, 0.1;
    Dmu->setMin( mu_min );
    auto mu_max = Dmu->element();
    mu_max << 50, 50, 5, 5;
    Dmu->setMax( mu_max );

    auto u = Xh->element();
    auto v = Xh->element();
    auto mesh = Xh->mesh();
    //lhs
    auto a0 = form2( _trial=Xh, _test=Xh);
    a0 = integrate( _range=elements( mesh ), _expr=0.1*( gradt( u )*trans( grad( v ) ) ) ) +
        integrate( _range=markedfaces( mesh,"right" ), _expr=idt( u )*id( v ) );
    this->addLhs( { a0 , "1" } );

    auto a1 = form2( _trial=Xh, _test=Xh);
    a1 = integrate( _range=markedelements( mesh,"k1_1" ), _expr=gradt( u )*trans( grad( v ) )  );
    this->addLhs( { a1 , "mu0" } );

    auto a2 = form2( _trial=Xh, _test=Xh);
    a2 = integrate( _range=markedelements( mesh,"k2_1" ), _expr=gradt( u )*trans( grad( v ) )  );
    this->addLhs( { a2 , "mu1" } );

    //rhs
    auto f0 = form1( _test=Xh );
    f0 = integrate( _range=markedfaces( mesh,"left" ), _expr=id( v ) );
    this->addRhs( { f0, "mu2" } );
    auto f1 = form1( _test=Xh );
    f1 =  integrate( _range=elements( mesh ), _expr=id( v ) );
    this->addRhs( { f1, "mu3" } );

    //output
    auto out = form1( _test=Xh );
    out = integrate( _range=markedelements( mesh,"k1_2" ), _expr=id( v )/0.2 ) +
          integrate( _range=markedelements( mesh,"k2_1" ), _expr=id( v )/0.2 );
    this->addOutput( { out, "1" } );

    auto energy = form2( _trial=Xh, _test=Xh);
    energy = integrate( _range=elements( mesh ), _expr=0.1*( gradt( u )*trans( grad( v ) ) ) ) +
             integrate( _range=markedfaces( mesh,"right" ), _expr=idt( u )*id( v ) ) +
             integrate( _range=markedelements( mesh,"k1_1" ), _expr=0.2 * gradt( u )*trans( grad( v ) ) )  +
             integrate( _range=markedelements( mesh,"k2_1" ), _expr=0.2 * gradt( u )*trans( grad( v ) ) )  ;
    this->addEnergyMatrix( energy );
}


double
Heat1D::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve )
{

//CHECK( ! need_to_solve ) << "The model need to have the solution to compute the output\n";

    auto mesh = Xh->mesh();
    double output=0;
    // right hand side (compliant)
    if ( output_index == 0 )
    {
        output  = integrate( _range=markedfaces( mesh,"left" ), _expr=mu(2)*idv( u ) ).evaluate()( 0,0 );
        output += integrate( _range=elements( mesh ), _expr=mu(3)*idv( u ) ).evaluate()( 0,0 );
    }
    // output
    else if ( output_index == 1 )
    {
        output = integrate( _range=elements( mesh ),
                            _expr=chi( ( Px() >= -0.1 ) && ( Px() <= 0.1 ) )*idv( u ) ).evaluate()( 0,0 )/0.2;
    }
    else
        throw std::logic_error( "[Heat1d::output] error with output_index : only 0 or 1 " );
    return output;

}


FEELPP_CRB_PLUGIN( Heat1D, heat1d )
}

