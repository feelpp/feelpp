/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Abdoulaye Samake <abdoulaye.samake1@ujf-grenoble.fr>
   Date: 2011-06-02

   Copyright (C) 2009-2014 Feel++ Consortium

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef FEELPP_ThermalBlockMinimalVersion_HPP
#define FEELPP_ThermalBlockMinimalVersion_HPP 1

#include <boost/timer.hpp>
#include <feel/options.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/pch.hpp>

#include <feel/feelfilters/gmsh.hpp>

#include <feel/feelcrb/modelcrbbase.hpp>

namespace Feel
{

/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Feel::Application subclass.
 *
 * \return the list of options
 */

po::options_description
makeThermalBlockMinimalVersionOptions()
{
    po::options_description thermalblockoptionsminimalversion( "ThermalBlock options" );
    thermalblockoptionsminimalversion.add_options()
    ( "gamma_dir", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary Dirichlet formulation" )
    ;
    return thermalblockoptionsminimalversion;
}

/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Feel::Application subclass.
 *
 * \return some data about the application.
 */
AboutData
makeThermalBlockMinimalVersionAbout( std::string const& str = "thermalBlockminimalversion" )
{
    AboutData about( "thermalblockminimalversion" ,
                     "thermalblockminimalversion" ,
                     "0.1",
                     "2D Heterogeneous Thermal Block Problem",
                     Feel::AboutData::License_GPL,
                     "Copyright (C) 2009-2014 Feel++ Consortium");
    about.addAuthor( "Abdoulaye Samake", "main developer", "abdoulaye.samake@ujf-grenoble.fr", "" );
    about.addAuthor( "Christophe Prud'homme", "contributor", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Stephane Veys", "contributor", "", "" );
    return about;

}


/**
 * \class ThermalBlock  Minimal Version
 *
 * ThermalBlock Solver using continuous approximation spaces
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 *
 */
class ThermalBlockMinimalVersion : public ModelCrbBase< ParameterSpace<8> , decltype(Pch<1>(Mesh<Simplex<2>>::New())) >
{
public:

    typedef ModelCrbBase< ParameterSpace<8> , decltype(Pch<1>(Mesh<Simplex<2>>::New())) > super_type;

    static const uint16_type ParameterSpaceDimension = 8;

    typedef typename super_type::element_type element_type;
    typedef typename super_type::parameter_type parameter_type;

    //! initialisation of the model and definition of parameters values
    void initModel();

    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false);


private:
    value_type gamma_dir;
}; // ThermalBlockMinimalVersion





void
ThermalBlockMinimalVersion::initModel()
{

    gamma_dir=doption(_name="gamma_dir");

    this->setFunctionSpaces( Pch<1>( loadMesh( _mesh=new Mesh<Simplex<2>> ) ) );

    auto mesh = Xh->mesh();

    if( Environment::worldComm().isMasterRank() )
    {
        std::cout << "Number of local dof " << Xh->nLocalDof() << "\n";
        std::cout << "Number of dof " << Xh->nDof() << "\n";
    }

    auto mu_min = Dmu->element();
    auto mu_max = Dmu->element();
    mu_min <<  0.1 , 0.1 , 0.1 , 0.1 , 0.1 , 0.1 , 0.1 , 0.1 ;
    mu_max <<  10 , 10 , 10 , 10 , 10 , 10 , 10 , 10 ;
    Dmu->setMin( mu_min );
    Dmu->setMax( mu_max );

    auto u = Xh->element();
    auto v = Xh->element();

    auto M = backend()->newMatrix( Xh , Xh );

    double mu_min_coeff=0.1;
    // on boundary north we have u=0 so term from weak dirichlet condition
    // vanish in the right hand side
    //rhs
    auto f0 = form1( _test=Xh );
    f0 =  integrate( _range=markedfaces( mesh,"south_domain-1" ), _expr= id( v ) )
        + integrate( _range=markedfaces( mesh,"south_domain-2" ), _expr= id( v ) )
        + integrate( _range=markedfaces( mesh,"south_domain-3" ), _expr= id( v ) );
    this->addRhs( { f0, "1" } );

    //lhs
    auto a0 = form2( _trial=Xh, _test=Xh);
    a0 = integrate( markedelements( mesh, "domain-1" ), gradt( u )*trans( grad( v ) ) );
    this->addLhs( { a0 , "1" } );
    auto a1 = form2( _trial=Xh, _test=Xh);
    a1 = integrate( markedelements( mesh, "domain-2" ), gradt( u )*trans( grad( v ) ) );
    this->addLhs( { a1 , "mu0" } );
    auto a2 = form2( _trial=Xh, _test=Xh);
    a2 = integrate( markedelements( mesh, "domain-3" ), gradt( u )*trans( grad( v ) ) );
    this->addLhs( { a2 , "mu1" } );
    auto a3 = form2( _trial=Xh, _test=Xh);
    a3 = integrate( markedelements( mesh, "domain-4" ), gradt( u )*trans( grad( v ) ) );
    this->addLhs( { a3 , "mu2" } );
    auto a4 = form2( _trial=Xh, _test=Xh);
    a4 = integrate( markedelements( mesh, "domain-5" ), gradt( u )*trans( grad( v ) ) );
    this->addLhs( { a4 , "mu3" } );
    auto a5 = form2( _trial=Xh, _test=Xh);
    a5 = integrate( markedelements( mesh, "domain-6" ), gradt( u )*trans( grad( v ) ) );
    this->addLhs( { a5 , "mu4" } );
    auto a6 = form2( _trial=Xh, _test=Xh);
    a6 = integrate( markedelements( mesh, "domain-7" ), gradt( u )*trans( grad( v ) ) )
        +integrate( markedfaces( mesh, "north_domain-7" ),
                   -gradt( u )*vf::N()*id( v )
                   -grad( u )*vf::N()*idt( v )
                   );
    this->addLhs( { a6 , "mu5" } );
    auto a7 = form2( _trial=Xh, _test=Xh);
    a7 = integrate( markedelements( mesh, "domain-8" ), gradt( u )*trans( grad( v ) ) )
        +integrate( markedfaces( mesh, "north_domain-8" ),
                   -gradt( u )*vf::N()*id( v )
                   -grad( u )*vf::N()*idt( v )
                   );
    this->addLhs( { a7 , "mu6" } );
    auto a8 = form2( _trial=Xh, _test=Xh);
    a8 = integrate( markedelements( mesh, "domain-9" ), gradt( u )*trans( grad( v ) ) )
        +integrate( markedfaces( mesh, "north_domain-9" ),
                   -gradt( u )*vf::N()*id( v )
                   -grad( u )*vf::N()*idt( v )
                   );
    this->addLhs( { a8 , "mu7" } );
    auto a9 = form2( _trial=Xh, _test=Xh);
    a9 = integrate( markedfaces( mesh, "north_domain-7" ),gamma_dir*idt( u )*id( v )/h() )
        +integrate( markedfaces( mesh, "north_domain-8" ),gamma_dir*idt( u )*id( v )/h() )
        +integrate( markedfaces( mesh, "north_domain-9" ),gamma_dir*idt( u )*id( v )/h() );
    this->addLhs( { a9 , "1" } );


    form2( _test=Xh, _trial=Xh, _matrix=M ) = integrate( markedelements( mesh, "domain-1" ), gradt( u )*trans( grad( v ) )  );
    form2( _test=Xh, _trial=Xh, _matrix=M ) += integrate( markedelements( mesh, "domain-2" ), gradt( u )*trans( grad( v ) ) * mu_min_coeff );
    form2( _test=Xh, _trial=Xh, _matrix=M ) += integrate( markedelements( mesh, "domain-3" ), gradt( u )*trans( grad( v ) ) * mu_min_coeff );
    form2( _test=Xh, _trial=Xh, _matrix=M ) += integrate( markedelements( mesh, "domain-4" ), gradt( u )*trans( grad( v ) ) * mu_min_coeff );
    form2( _test=Xh, _trial=Xh, _matrix=M ) += integrate( markedelements( mesh, "domain-5" ), gradt( u )*trans( grad( v ) ) * mu_min_coeff );
    form2( _test=Xh, _trial=Xh, _matrix=M ) += integrate( markedelements( mesh, "domain-6" ), gradt( u )*trans( grad( v ) ) * mu_min_coeff );
    form2( _test=Xh, _trial=Xh, _matrix=M ) += integrate( markedelements( mesh, "domain-7" ), gradt( u )*trans( grad( v ) ) * mu_min_coeff );
    form2( _test=Xh, _trial=Xh, _matrix=M ) += integrate( markedelements( mesh, "domain-8" ), gradt( u )*trans( grad( v ) ) * mu_min_coeff );
    form2( _test=Xh, _trial=Xh, _matrix=M ) += integrate( markedelements( mesh, "domain-9" ), gradt( u )*trans( grad( v ) ) * mu_min_coeff );
    form2( _test=Xh, _trial=Xh, _matrix=M ) +=  integrate( markedfaces( mesh, "north_domain-7" ),
                                      -gradt( u )*vf::N()*id( v ) * mu_min_coeff
                                      -grad( u )*vf::N()*idt( v ) * mu_min_coeff
                                      );
    form2( _test=Xh, _trial=Xh, _matrix=M ) +=  integrate( markedfaces( mesh, "north_domain-8" ),
                                      -gradt( u )*vf::N()*id( v ) * mu_min_coeff
                                      -grad( u )*vf::N()*idt( v ) * mu_min_coeff
                                      );
    form2( _test=Xh, _trial=Xh, _matrix=M ) +=  integrate( markedfaces( mesh, "north_domain-9" ),
                                      -gradt( u )*vf::N()*id( v ) * mu_min_coeff
                                      -grad( u )*vf::N()*idt( v ) * mu_min_coeff
                                      );
    form2( _test=Xh, _trial=Xh, _matrix=M ) += integrate( markedfaces( mesh, "north_domain-7" ),gamma_dir*idt( u )*id( v )/h() );
    form2( _test=Xh, _trial=Xh, _matrix=M ) += integrate( markedfaces( mesh, "north_domain-8" ),gamma_dir*idt( u )*id( v )/h() );
    form2( _test=Xh, _trial=Xh, _matrix=M ) += integrate( markedfaces( mesh, "north_domain-9" ),gamma_dir*idt( u )*id( v )/h() );
    this->addEnergyMatrix( M );

}//initModel()

double
ThermalBlockMinimalVersion::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve)
{

    double output=0;
    auto mesh = Xh->mesh();

    if ( output_index==0 )
    {
        output =  integrate( _range=markedfaces( mesh,"south_domain-1" ), _expr= idv( u ) ).evaluate()(0,0)
                + integrate( _range=markedfaces( mesh,"south_domain-2" ), _expr= idv( u ) ).evaluate()(0,0)
                + integrate( _range=markedfaces( mesh,"south_domain-3" ), _expr= idv( u ) ).evaluate()(0,0);
    }
    else
    {
        throw std::logic_error( "[ThermalBlock::output] error with output_index : only 0 " );
    }

    return output;
}//output


} // Feel

#endif /* __ThermalBlock_H */
