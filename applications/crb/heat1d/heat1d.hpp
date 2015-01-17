/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-11-13

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
   \file heat1d.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-11-13
 */
#ifndef FEELPP_HEAT1D_HPP
#define FEELPP_HEAT1D_HPP 1

#include <boost/timer.hpp>

#include <feel/options.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/gmsh.hpp>

#include <feel/feelcrb/modelcrbbase.hpp>


namespace Feel
{

po::options_description
makeHeat1DOptions()
{
    po::options_description heat1doptions( "Heat1D options" );
    // heat1doptions.add_options()
    // ( "mu1", po::value<double>()->default_value( 0.2 ), "mu1" )
    // ( "mu2", po::value<double>()->default_value( 0.2 ), "mu2" )
    // ( "mu3", po::value<double>()->default_value(-1.0 ), "mu3" )
    // ( "mu4", po::value<double>()->default_value( 0.1 ), "mu4" )
    // ;
    return heat1doptions;
}
AboutData
makeHeat1DAbout( std::string const& str = "heat1d" )
{
    Feel::AboutData about( /*AppName  */ str.c_str(),
                           /*ProgName */ str.c_str(),
                           /*Version  */ "0.1",
                           /*ShortDesc*/ "1D Heat Benchmark",
                           /*Licence  */ Feel::AboutData::License_GPL,
                           /*Copyright*/ "Copyright (c) 2009-2014 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}

/**
 * \class Heat1D
 * \brief brief description
 *
 * @author Christophe Prud'homme
 * @see
 */
class Heat1D : public ModelCrbBase<ParameterSpace<4>, decltype(Pch<5>(Mesh<Simplex<1>>::New()))>
{
public:

    typedef ModelCrbBase<ParameterSpace<4>, decltype(Pch<5>(Mesh<Simplex<1>>::New()))> super_type;

    //! initialisation of the model
    void initModel();

    beta_vector_light_type beta;
    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false);
};


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
    a0 = integrate( elements( mesh ), 0.1*( gradt( u )*trans( grad( v ) ) ) ) +
        integrate( markedfaces( mesh,"right" ), idt( u )*id( v ) );
    this->addLhs( { a0 , "1" } );

    auto a1 = form2( _trial=Xh, _test=Xh);
    a1 = integrate( markedelements( mesh,"k1_1" ),  gradt( u )*trans( grad( v ) )  );
    this->addLhs( { a1 , "mu0" } );

    auto a2 = form2( _trial=Xh, _test=Xh);
    a2 = integrate( markedelements( mesh,"k2_1" ),  gradt( u )*trans( grad( v ) )  );
    this->addLhs( { a2 , "mu1" } );

    //rhs
    auto f0 = form1( _test=Xh );
    f0 = integrate( markedfaces( mesh,"left" ), id( v ) );
    this->addRhs( { f0, "mu2" } );
    auto f1 = form1( _test=Xh );
    f1 =  integrate( elements( mesh ), id( v ) );
    this->addRhs( { f1, "mu3" } );

    //output
    auto out = form1( _test=Xh );
    out = integrate( markedelements( mesh,"k1_2" ), id( v )/0.2 ) +
          integrate( markedelements( mesh,"k2_1" ), id( v )/0.2 );
    this->addOutput( { out, "1" } );

    auto energy = form2( _trial=Xh, _test=Xh);
    energy = integrate( elements( mesh ), 0.1*( gradt( u )*trans( grad( v ) ) ) ) +
             integrate( markedfaces( mesh,"right" ), idt( u )*id( v ) ) +
             integrate( markedelements( mesh,"k1_1" ),  0.2 * gradt( u )*trans( grad( v ) ) )  +
             integrate( markedelements( mesh,"k2_1" ),  0.2 * gradt( u )*trans( grad( v ) ) )  ;
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
        output  = integrate( markedfaces( mesh,"left" ), mu(2)*idv( u ) ).evaluate()( 0,0 );
        output += integrate( elements( mesh ), mu(3)*idv( u ) ).evaluate()( 0,0 );
    }
    // output
    else if ( output_index == 1 )
    {
        output = integrate( elements( mesh ),
                            chi( ( Px() >= -0.1 ) && ( Px() <= 0.1 ) )*idv( u ) ).evaluate()( 0,0 )/0.2;
    }
    else
        throw std::logic_error( "[Heat1d::output] error with output_index : only 0 or 1 " );
    return output;

}

}

#endif /* FEELPP_HEAT1D_HPP */
