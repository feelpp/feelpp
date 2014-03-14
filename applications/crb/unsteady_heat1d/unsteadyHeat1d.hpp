/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-11-13

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
   \file unsteadyHeat1d.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-11-13
 */
#ifndef __UnsteadyHeat1D_H
#define __UnsteadyHeat1D_H 1

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

#include <feel/feelts/bdf.hpp>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>

#include <feel/feelcrb/modelcrbbase.hpp>
#include <feel/feeldiscr/reducedbasisspace.hpp>

namespace Feel
{

po::options_description
makeUnsteadyHeat1DOptions()
{
    po::options_description unsteadyheat1doptions( "UnsteadyHeat1D options" );
    unsteadyheat1doptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.01 ), "mesh size" )
    ( "mu1", po::value<double>()->default_value( 0.2 ), "mu1" )
    ( "mu2", po::value<double>()->default_value( 0.2 ), "mu2" )
    ( "mu3", po::value<double>()->default_value( -1 ), "mu3" )
    ( "mu4", po::value<double>()->default_value( 0.1 ), "mu4" )
    ( "alpha" , po::value<double>()->default_value( 1 ), "temporal coefficient" )
    ;
    return unsteadyheat1doptions.add( Feel::feel_options() ).add( bdf_options( "unsteadyHeat1d" ) );
}
AboutData
makeUnsteadyHeat1DAbout( std::string const& str = "unsteadyHeat1d" )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "unsteady 1D Heat Benchmark",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2010,2011 Universite de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Stephane Veys", "developer", "stephane.veys@imag.fr", "" );
    return about;
}



gmsh_ptrtype
createGeo( double meshSize  )
{
    std::ostringstream ostr;
    ostr << "P=" << meshSize << ";\n"
         << "Point(1) = {-0.5,0,0,P};\n"
         << "Point(2) = {-0.1,0,0,P};\n"
         << "Point(3) = {0,0,0,P};\n"
         << "Point(4) = {0.1,0,0.0,P};\n"
         << "Point(5) = {0.5,0,0.0,P};\n"
         << "\n"

         << "Line(1) = {1, 2};\n"
         << "Line(2) = {2, 3};\n"
         << "Line(3) = {3, 4};\n"
         << "Line(4) = {4, 5};\n"
         << "Line Loop(5) = {1,2,3,4};\n"
         << "Physical Line(\"k1_1\") = {1};\n"
         << "Physical Line(\"k1_2\") = {2};\n"
         << "Physical Line(\"k2_1\") = {3};\n"
         << "Physical Line(\"k2_2\") = {4};\n"
         << "Physical Point(\"left\") = {1};\n"
         << "Physical Point(\"right\") = {5};\n"
         << "\n";

    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setPrefix( "unsteadyheat1d" );
    gmshp->setDescription( ostr.str() );
    return gmshp;
}


class ParameterDefinition
{
public :
    static const uint16_type ParameterSpaceDimension = 4;
    typedef ParameterSpace<ParameterSpaceDimension> parameterspace_type;
};
class FunctionSpaceDefinition
{
public :
    static const uint16_type Order = 1;

    typedef double value_type;

    /*mesh*/
    typedef Simplex<1,1> entity_type;
    typedef Mesh<entity_type> mesh_type;

    /*basis*/
    typedef bases<Lagrange<Order, Scalar> > basis_type;

    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;

    static const bool is_time_dependent = true;
    static const bool is_linear = true;

};


/**
 * \class UnsteadyHeat1D
 * \brief brief description
 *
 * @author Christophe Prud'homme
 * @see
 */
class UnsteadyHeat1D : public ModelCrbBase<ParameterDefinition,FunctionSpaceDefinition>,
                       public boost::enable_shared_from_this< UnsteadyHeat1D >
{
public:


    typedef ModelCrbBase<ParameterDefinition,FunctionSpaceDefinition> super_type;
    typedef typename super_type::funs_type funs_type;
    typedef typename super_type::funsd_type funsd_type;
    typedef typename super_type::beta_vector_light_type beta_vector_light_type;
    typedef typename super_type::affine_decomposition_light_type affine_decomposition_light_type;
    using super_type::computeBetaQm;

    /** @name Typedefs
     */
    //@{

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

    /*reduced basis space*/
    typedef ReducedBasisSpace<super_type, mesh_type, basis_type, value_type> rbfunctionspace_type;
    typedef boost::shared_ptr< rbfunctionspace_type > rbfunctionspace_ptrtype;

    /* parameter space */
    typedef ParameterDefinition::parameterspace_type parameterspace_type;
    typedef boost::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef parameterspace_type::element_type parameter_type;

    /* time discretization */
    typedef Bdf<space_type>  bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    UnsteadyHeat1D();

    //! constructor from command line
    UnsteadyHeat1D( po::variables_map const& vm );


    //! copy constructor
    UnsteadyHeat1D( UnsteadyHeat1D const & );
    //! destructor
    ~UnsteadyHeat1D() {}

    //! initialisation of the model
    void initModel();
    //@}

    /** @name Operator overloads
     */
    //@{

    //@}

    /** @name Accessors
     */
    //@{


    /**
     * \brief Returns the function space
     */
    space_ptrtype functionSpace()
    {
        return Xh;
    }
    /**
     * \brief Returns the reduced basis function space
     */
    rbfunctionspace_ptrtype rBFunctionSpace()
    {
        return RbXh;
    }

    //! return the parameter space
    parameterspace_ptrtype parameterSpace() const
    {
        return M_Dmu;
    }


    //@}


    /** @name  Methods
     */
    //@{


    void assemble();

    /**
     * Matrix associated to scalar product used
     */
    sparse_matrix_ptrtype innerProduct ( void )
    {
        return M;
    }

    /**
     * inner product for mass matrix
     */
    sparse_matrix_ptrtype innerProductForMassMatrix ( void )
    {
        return InnerMassMatrix;
    }

    /**
     * inner product for the POD
     */
    sparse_matrix_ptrtype innerProductForPod ( void )
    {
        return Mpod;
    }

    value_type output( int output_index, parameter_type const& mu, element_type &u, bool need_to_solve=false, bool export_outputs=false );

    bdf_ptrtype bdfModel(){ return M_bdf; }

private:

    po::variables_map M_vm;

    backend_ptrtype backend;

    double alpha;

    double meshSize;

    parameterspace_ptrtype M_Dmu;

    mesh_ptrtype mesh;
    space_ptrtype Xh;
    rbfunctionspace_ptrtype RbXh;
    sparse_matrix_ptrtype M,Mpod,InnerMassMatrix;

    bdf_ptrtype M_bdf;

    int M_nb_terms_in_affine_decomposition_a=3;
    int M_nb_terms_in_affine_decomposition_m=1;
    int M_nb_terms_in_affine_decomposition_rhs=2;
    int M_nb_terms_in_affine_decomposition_output=1;
    int M_nb_outputs=2;

};

UnsteadyHeat1D::UnsteadyHeat1D()
    :
    backend( backend_type::build( BACKEND_PETSC ) ),
    alpha( 1 ),
    meshSize( 0.01 ),
    M_Dmu( new parameterspace_type )
{
}


UnsteadyHeat1D::UnsteadyHeat1D( po::variables_map const& vm )
    :
    M_vm( vm ),
    backend( backend_type::build( vm ) ),
    alpha( vm["alpha"].as<double>() ),
    meshSize( vm["hsize"].as<double>() ),
    M_Dmu( new parameterspace_type )
{
}
void
UnsteadyHeat1D::initModel()
{


    /*
     * First we create the mesh
     */
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=createGeo( meshSize ) );



    /*
     * The function space and some associate elements are then defined
     */
    Xh = space_type::New( mesh );
    RbXh = rbfunctionspace_type::New( _model=this->shared_from_this() , _mesh=mesh );

    if( Environment::worldComm().isMasterRank() )
    {
        std::cout << "Number of local dof " << Xh->nLocalDof() << "\n";
        std::cout << "Number of dof " << Xh->nDof() << "\n";
    }


    M_bdf = bdf( _space=Xh, _vm=M_vm, _name="unsteadyHeat1d" , _prefix="unsteadyHeat1d" );

    Feel::ParameterSpace<4>::Element mu_min( M_Dmu );
    mu_min << 0.2, 0.2, 0.01, 0.1;
    M_Dmu->setMin( mu_min );
    Feel::ParameterSpace<4>::Element mu_max( M_Dmu );
    mu_max << 50, 50, 5, 5;
    M_Dmu->setMax( mu_max );

    Mpod = backend->newMatrix( Xh, Xh );

    InnerMassMatrix = backend->newMatrix( Xh,Xh );

    M = backend->newMatrix( Xh, Xh );

    LOG(INFO) << "Number of dof " << Xh->nLocalDof() << "\n";

    assemble();


} // UnsteadyHeat1d::initModel


void
UnsteadyHeat1D::assemble()
{
    auto u = Xh->element();
    auto v = Xh->element();

    //lhs
    auto a0 = form2( _trial=Xh, _test=Xh);
    a0 = integrate( elements( mesh ), 0.1*( gradt( u )*trans( grad( v ) ) ) ) +
         integrate( markedfaces( mesh,"right" ), idt( u )*id( v ) );
    this->addLhs( boost::make_tuple( a0.matrixPtr() , "1" ) );

    auto a1 = form2( _trial=Xh, _test=Xh);
    a1 = integrate( markedelements( mesh,"k1_1" ),  gradt( u )*trans( grad( v ) )  );
    this->addLhs( boost::make_tuple( a1.matrixPtr() , "mu0" ) );

    auto a2 = form2( _trial=Xh, _test=Xh);
    a2 = integrate( markedelements( mesh,"k2_1" ),  gradt( u )*trans( grad( v ) )  );
    this->addLhs( boost::make_tuple( a2.matrixPtr() , "mu1" ) );

    //mass matrix
    auto m = form2( _trial=Xh, _test=Xh);
    m = integrate ( elements( mesh ), alpha*idt( u )*id( v ) );
    this->addMass( boost::make_tuple( m.matrixPtr() , "1" ) );

    //rhs
    auto f0 = form1( _test=Xh );
    auto f1 = form1( _test=Xh );
    f0 = integrate( markedfaces( mesh,"left" ), id( v ) );
    this->addRhs( boost::make_tuple( f0.vectorPtr() , "mu2" ) );
    f1 =  integrate( elements( mesh ), id( v ) );
    this->addRhs( boost::make_tuple( f1.vectorPtr() , "mu3" ) );

    //output
    auto o = form1( _test=Xh );
    o = integrate( markedelements( mesh,"k1_2" ), id( v )/0.2 ) +
        integrate( markedelements( mesh,"k2_1" ), id( v )/0.2 );
    this->addOutput( boost::make_tuple( o.vectorPtr() , "1" ) );


    form2( Xh, Xh, Mpod ) =
        integrate( elements( mesh ), 0.1*( gradt( u )*trans( grad( v ) ) ) ) +
        integrate( markedfaces( mesh,"right" ), idt( u )*id( v ) ) +
        integrate( markedelements( mesh,"k1_1" ),  0.2 * gradt( u )*trans( grad( v ) ) )  +
        integrate( markedelements( mesh,"k2_1" ),  0.2 * gradt( u )*trans( grad( v ) ) )  ;

    form2( Xh, Xh, M ) =
        integrate( elements( mesh ), 0.1*( gradt( u )*trans( grad( v ) ) ) ) +
        integrate( markedfaces( mesh,"right" ), idt( u )*id( v ) ) +
        integrate( markedelements( mesh,"k1_1" ),  0.2 * gradt( u )*trans( grad( v ) ) )  +
        integrate( markedelements( mesh,"k2_1" ),  0.2 * gradt( u )*trans( grad( v ) ) )  ;

    form2( Xh, Xh, InnerMassMatrix ) =
        integrate( _range=elements( mesh ), _expr=idt( u ) * id( v ) ) ;

}


double
UnsteadyHeat1D::output( int output_index, parameter_type const& mu, element_type& u, bool need_to_solve, bool export_output )
{

    CHECK( ! need_to_solve ) << "The model need to have the solution to compute the output\n";

    auto v = Xh->element();

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
        throw std::logic_error( "[unsteadyHeat1d::output] error with output_index : only 0 or 1 " );

    return output;
}

}

#endif /* __Heat1D_H */


