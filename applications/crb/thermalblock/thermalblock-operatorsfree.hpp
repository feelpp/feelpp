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
#ifndef FEELPP_ThermalBlockFree_HPP
#define FEELPP_ThermalBlockFree_HPP 1

#include <boost/timer.hpp>
#include <feel/options.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/pch.hpp>

#include <feel/feelfilters/gmsh.hpp>

#include <feel/feelcrb/crbmodelbase.hpp>

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
makeThermalBlockFreeOptions()
{
    po::options_description thermalblockoptionsfree( "ThermalBlock options" );
    thermalblockoptionsfree.add_options()
    ( "gamma_dir", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary Dirichlet formulation" )
    ;
    return thermalblockoptionsfree;
}

/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Feel::Application subclass.
 *
 * \return some data about the application.
 */
AboutData
makeThermalBlockFreeAbout( std::string const& str = "thermalBlockfree" )
{
    AboutData about( "thermalblockfree" ,
                     "thermalblockfree" ,
                     "0.1",
                     "2D Heterogeneous Thermal Block Problem",
                     Feel::AboutData::License_GPL,
                     "Copyright (C) 2009-2014 Feel++ Consortium");
    about.addAuthor( "Abdoulaye Samake", "main developer", "abdoulaye.samake@ujf-grenoble.fr", "" );
    about.addAuthor( "Christophe Prud'homme", "contributor", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Stephane Veys", "contributor", "", "" );
    return about;

}

class FunctionSpaceDefinition
{
public :
    typedef double value_type;

    static const uint16_type Order = 1;

    //! geometry entities type composing the mesh, here Simplex in Dimension 2 of Order 1
    typedef Simplex<2> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;

    //! the basis type of our approximation space
    typedef bases<Lagrange<Order,Scalar> > basis_type;

    //! the approximation function space type
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;

    typedef typename space_type::element_type element_type;
};


/**
 * \class ThermalBlock using operators free
 *
 * ThermalBlock Solver using continuous approximation spaces
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 *
 */
class ThermalBlockFree : public CRBModelBase< ParameterSpace<8> , FunctionSpaceDefinition >
{
public:
    typedef CRBModelBase<ParameterSpace<8>, FunctionSpaceDefinition> super_type;

    //! initialisation of the model and definition of parameters values
    void initModel();

    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false);

}; // ThermalBlockFree


void
ThermalBlockFree::initModel()
{
    double gamma_dir=doption(_name="gamma_dir");

    auto mesh = loadMesh( new mesh_type );

    // Definition of the function space
    auto Xh = functionspace_type::New( mesh );
    this->setFunctionSpaces( Xh );
    if( Environment::worldComm().isMasterRank() )
    {
        std::cout << "Number of local dof " << Xh->nLocalDof() << std::endl;
        std::cout << "Number of dof " << Xh->nDof() << std::endl;
    }
    u = Xh->element();
    v = Xh->element();

    // Seting the parameter space
    auto mu_min = Dmu->element();
    auto mu_max = Dmu->element();
    for ( int i=0; i<8; i++ )
    {
        mu_min[i]=0.1;
        mu_max[i]=10;
    }
    Dmu->setMin( mu_min );
    Dmu->setMax( mu_max );


    // on boundary north we have u=0 so term from weak dirichlet condition
    // vanish in the right hand side
    auto expr_f00 =
          integrate( _range=markedfaces( mesh,"south_domain-1" ), _expr= id( v ) )
        + integrate( _range=markedfaces( mesh,"south_domain-2" ), _expr= id( v ) )
        + integrate( _range=markedfaces( mesh,"south_domain-3" ), _expr= id( v ) );
    addRhs( expr_f00, "1" );


    auto expr_a0 = integrate( markedelements( mesh, "domain-2" ), gradt( u )*trans( grad( v ) ) );
    addLhs( expr_a0, "mu0");

    auto expr_a1 = integrate( markedelements( mesh, "domain-3" ), gradt( u )*trans( grad( v ) ) );
    addLhs( expr_a1, "mu1");

    auto expr_a2 = integrate( markedelements( mesh, "domain-4" ), gradt( u )*trans( grad( v ) ) );
    addLhs( expr_a2, "mu2");

    auto expr_a3 = integrate( markedelements( mesh, "domain-5" ), gradt( u )*trans( grad( v ) ) );
    addLhs( expr_a3, "mu3");

    auto expr_a4 = integrate( markedelements( mesh, "domain-6" ), gradt( u )*trans( grad( v ) ) );
    addLhs( expr_a4, "mu4");

    auto expr_a5 = integrate( markedelements( mesh, "domain-7" ), gradt( u )*trans( grad( v ) ) )
        +integrate( markedfaces( mesh, "north_domain-7" ),
                    -gradt( u )*vf::N()*id( v )
                    -grad( u )*vf::N()*idt( v ) );
    addLhs( expr_a5, "mu5");

    auto expr_a6 = integrate( markedelements( mesh, "domain-8" ), gradt( u )*trans( grad( v ) ) )
        +integrate( markedfaces( mesh, "north_domain-8" ),
                    -gradt( u )*vf::N()*id( v )
                    -grad( u )*vf::N()*idt( v ) );
    addLhs( expr_a6, "mu6");

    auto expr_a7 = integrate( markedelements( mesh, "domain-9" ), gradt( u )*trans( grad( v ) ) )
        +integrate( markedfaces( mesh, "north_domain-9" ),
                    -gradt( u )*vf::N()*id( v )
                    -grad( u )*vf::N()*idt( v )  );
    addLhs( expr_a7, "mu7");

    auto expr_a8 = integrate( markedelements( mesh, "domain-1" ), gradt( u )*trans( grad( v ) ) )
        +integrate( markedfaces( mesh, "north_domain-7" ),gamma_dir*idt( u )*id( v )/h() )
        +integrate( markedfaces( mesh, "north_domain-8" ),gamma_dir*idt( u )*id( v )/h() )
        +integrate( markedfaces( mesh, "north_domain-9" ),gamma_dir*idt( u )*id( v )/h() );
    addLhs( expr_a8, "1");



    auto M = backend()->newMatrix( Xh , Xh );
    double mu_min_coeff=0.1;

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

    addEnergyMatrix( M );

}//initModel()

double
ThermalBlockFree::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve)
{
    auto mesh = Xh->mesh();
    double output=0;

    if ( output_index==0 )
    {
        output =
            integrate( markedfaces( mesh,"south_domain-1" ), idv( u ) ).evaluate() (0,0)
            +integrate( markedfaces( mesh,"south_domain-2" ), idv( u ) ).evaluate() (0,0)
            +integrate( markedfaces( mesh,"south_domain-3" ), idv( u ) ).evaluate() (0,0);
    }
    else
    {
        throw std::logic_error( "[ThermalBlock::output] error with output_index : only 0 " );
    }

    return output;
}//output


} // Feel

#endif /* __ThermalBlock_H */
