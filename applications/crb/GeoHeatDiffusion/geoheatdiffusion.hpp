/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  Author(s): Stephane Veys <stephane.veys@imag.fr>
       Date: 2014-01-19

  Copyright (C) 2011-2014 Feel++ Consortium

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
   \file geoheatdiffusion.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Stephane Veys <stephane.veys@imag.fr>

   date 2014-03-15
 */
#ifndef __GeoHeatDiffusion_H
#define __GeoHeatDiffusion_H 1

#include <boost/timer.hpp>
#include <boost/shared_ptr.hpp>

#include <feel/options.hpp>
#include <feel/feelcore/feel.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelalg/solvereigen.hpp>

#include <feel/feelvf/vf.hpp>
#include <feel/feelcrb/parameterspace.hpp>
#include <feel/feelcrb/modelcrbbase.hpp>

#include <feel/feelmesh/meshmover.hpp>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>

#include <feel/feeldiscr/reducedbasisspace.hpp>



namespace Feel
{

po::options_description
makeGeoHeatDiffusionOptions()
{
    po::options_description options( "geoHeatDiffusion options" );
    return options.add( bdf_options( "geoheatdiffusion" ) );
}
AboutData
makeGeoHeatDiffusionAbout( std::string const& str = "GeoHeatDiffusion" )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "Heat Diffusion with geometrical parameter",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2011-2014 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Stephane Veys", "developer", "stephane.veys@imag.fr", "" );
    return about;
}

class FunctionSpaceDefinition
{
public :
    static const uint16_type Order = 1;
    typedef double value_type;

    /*mesh*/
    typedef Simplex<2,1> entity_type; /*dim,order*/
    typedef Mesh<entity_type> mesh_type;

    /*basis*/
    typedef bases<Lagrange<Order, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;

    /*element*/
    typedef typename space_type::element_type element_type;

    static const bool is_time_dependent = true;
    static const bool is_linear = true;

};


/**
 * \class GeoHeatDiffusion
 * @author Christophe Prud'homme
 * @author Stephane Veys
 * @see
 */
class GeoHeatDiffusion : public ModelCrbBase< ParameterSpace<3>, FunctionSpaceDefinition >
{
public:

    typedef ModelCrbBase<ParameterSpace<3>,FunctionSpaceDefinition> super_type;

    typedef double value_type;

    typedef typename FunctionSpaceDefinition::mesh_type mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef typename FunctionSpaceDefinition::basis_type basis_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef backend_type::vector_ptrtype vector_ptrtype;

    /*space*/
    typedef typename FunctionSpaceDefinition::space_type space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef space_type::element_type element_type;

    /* time discretization */
    typedef Bdf<space_type>  bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    GeoHeatDiffusion();

    //! constructor from command line
    GeoHeatDiffusion( po::variables_map const& vm );


    //! copy constructor
    GeoHeatDiffusion( GeoHeatDiffusion const & );
    //! destructor
    ~GeoHeatDiffusion() {}

    //! initialization of the model
    void initModel();
    //@}

    void assemble();

    bool referenceParametersGivenByUser() { return true; }

    parameter_type refParameter()
    {
        auto muref = Dmu->element();
        muref << 1,1,1;
        return muref;
    }

    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu, element_type &T, bool need_to_solve=false, bool export_outputs=false );

    bdf_ptrtype bdfModel(){ return M_bdf; }

    void adaptMesh( parameter_type const& mu )
    {
        //create a vectorial function space to project geometric transformation
        auto Xhv = Pchv<1>( mesh );
        auto transformation_field = vf::project( _space=Xhv, _expr=oneX()*( Px()*mu(0) - Px() ) );
        meshMove( mesh , transformation_field );
    }

private:
    /* mesh, pointers and spaces */
    mesh_ptrtype mesh;
    double surface;
    bdf_ptrtype M_bdf;
};

GeoHeatDiffusion::GeoHeatDiffusion()
{}

GeoHeatDiffusion::GeoHeatDiffusion( po::variables_map const& vm )
{}

void GeoHeatDiffusion::initModel()
{

    // geometry is a ]0,1[x]0,1[
    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 1,1 );
    GeoTool::Rectangle R( doption( "gmsh.hsize"),"Omega",x1,x2 );
    R.setMarker( _type="line",_name="heat",_marker4=true );
    R.setMarker( _type="line",_name="iso",_marker1=true );
    R.setMarker( _type="line",_name="iso",_marker3=true );
    R.setMarker( _type="line",_name="cool",_marker2=true );
    R.setMarker( _type="surface",_name="Omega",_markerAll=true );
    mesh = R.createMesh( _mesh=new mesh_type, _name="Omega" );

    /*
     * The function space and some associate elements are then defined
     */
    Xh = space_type::New( mesh );
    this->setFunctionSpaces( Xh );
    if( Environment::worldComm().isMasterRank() )
    {
        std::cout << "Number of dof " << Xh->nDof() << "\n";
        std::cout << "Number of local dof " << Xh->nLocalDof() << "\n";
    }

    auto mu_min = Dmu->element();
    mu_min <<  /*length*/0.1, /*heat transfer coefficient*/ 0.1 , 0.1;
    Dmu->setMin( mu_min );
    auto mu_max = Dmu->element();
    mu_max << 5, 5 , 5;
    Dmu->setMax( mu_max );

    M_bdf = bdf( _space=Xh, _name="geoheatdiffusion" , _prefix="geoheatdiffusion" );

    assemble();

} // GeoHeatDiffusion::init


void GeoHeatDiffusion::assemble()
{

    auto u = Xh->element();
    auto v = Xh->element();

    double uair=1;
    double ucool=0;
    double k=1;
    auto a0 = form2( _trial=Xh, _test=Xh);
    a0 = integrate( elements( mesh ),k* dxt(u)*dx(v) ) ;
    this->addLhs( { a0 , "1.0/mu0" } );

    auto a1 = form2( _trial=Xh, _test=Xh);
    a1 = integrate( elements( mesh ), k*dyt(u)*dy(v) );
    this->addLhs( { a1 , "mu0" } );

    auto a2 = form2( _trial=Xh, _test=Xh);
    a2 = integrate( markedfaces( mesh,"heat" ), idt( u )*id( v ) ) ;
    this->addLhs( { a2 , "mu0*mu1" } );

    auto a3 = form2( _trial=Xh, _test=Xh);
    a3 = integrate( markedfaces( mesh,"cool" ), idt( u )*id( v ) ) ;
    this->addLhs( { a3 , "mu0*mu2" } );

    auto m = form2( _trial=Xh, _test=Xh);
    m = integrate ( elements( mesh ), idt( u )*id( v ) );
    this->addMass( { m , "1" } );

    //rhs
    auto f0 = form1( _test=Xh );
    f0 = integrate( markedfaces( mesh,"heat" ), uair*id( v ) );
    this->addRhs( { f0 , "mu0*mu1" } );
    auto f1 = form1( _test=Xh );
    f1 = integrate( markedfaces( mesh,"cool" ), ucool*id( v ) );
    this->addRhs( { f1 , "mu0*mu2" } );

    //output
    auto out = form1( _test=Xh );
    out = integrate( elements( mesh ), id( v ) ) ;
    //surface = mu(0) * 1 = mu(0)
    this->addOutput( { out , "1.0/mu0" } );

    //energy matrix
    auto energy = form2( _trial=Xh, _test=Xh);
    energy=
        integrate( elements( mesh ), 1 * dxt(u)*dx(v) ) +
        integrate( elements( mesh ), 1 * dyt(u)*dy(v) ) +
        integrate( markedfaces( mesh,"heat" ), 1 * idt( u )*id( v ) )+
        integrate( markedfaces( mesh,"cool" ), 1 * idt( u )*id( v ) );
    this->addEnergyMatrix( energy );

    auto mass = form2( _trial=Xh, _test=Xh);
    mass = integrate( _range=elements( mesh ), _expr=idt( u ) * id( v ) ) ;
    this->addMassMatrix( mass );

}


double GeoHeatDiffusion::output( int output_index, parameter_type const& mu, element_type &u, bool need_to_solve , bool export_outputs )
{
    CHECK( ! need_to_solve ) << "The model need to have the solution to compute the output\n";

    auto v = Xh->element();

    double output=0;
    if ( output_index == 0 )
    {
        output  = integrate( markedfaces( mesh,"heat" ), mu(0)*mu(1) * idv( u ) ).evaluate()( 0,0 );
    }
    if ( output_index == 1 )
    {
        surface = mu(0);
        output = integrate( elements( mesh ), idv( u )/surface ).evaluate()( 0,0 );
    }
    if ( output_index>=2 )
    {
        throw std::logic_error( "[GeoHeatDiffusion::output] error with output_index : only 0 or 1 " );
    }

    return output ;
}

}

#endif /* __GeoHeatDiffusion_H */
