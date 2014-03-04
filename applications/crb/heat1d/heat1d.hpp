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
   \file heat1d.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-11-13
 */
#ifndef __Heat1D_H
#define __Heat1D_H 1

#include <boost/timer.hpp>
#include <boost/shared_ptr.hpp>

#include <feel/options.hpp>
#include <feel/feelcore/feel.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelalg/solvereigen.hpp>

#include <feel/feelvf/vf.hpp>
#include <feel/feelcrb/parameterspace.hpp>

#include <feel/feelcrb/modelcrbbase.hpp>
#include <feel/feeldiscr/reducedbasisspace.hpp>

namespace Feel
{

po::options_description
makeHeat1DOptions()
{
    po::options_description heat1doptions( "Heat1D options" );
    heat1doptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.01 ), "mesh size" )
    ( "mu1", po::value<double>()->default_value( 0.2 ), "mu1" )
    ( "mu2", po::value<double>()->default_value( 0.2 ), "mu2" )
    ( "mu3", po::value<double>()->default_value( -1 ), "mu3" )
    ( "mu4", po::value<double>()->default_value( 0.1 ), "mu4" )
    ;
    return heat1doptions.add( Feel::feel_options() );
}
AboutData
makeHeat1DAbout( std::string const& str = "heat1d" )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "1D Heat Benchmark",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2010-2012 Universite de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
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
    gmshp->setPrefix( "heat1d" );
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
    static const uint16_type Order = 5;

    typedef double value_type;

    /*mesh*/
    typedef Simplex<1,1> entity_type;
    typedef Mesh<entity_type> mesh_type;

    /*basis*/
    typedef bases<Lagrange<Order, Scalar> > basis_type;

    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;

    static const bool is_time_dependent = false;
    static const bool is_linear = true;
};

/**
 * \class Heat1D
 * \brief brief description
 *
 * @author Christophe Prud'homme
 * @see
 */
class Heat1D : public ModelCrbBase<ParameterDefinition, FunctionSpaceDefinition> ,
               public boost::enable_shared_from_this< Heat1D >
{
public:

    typedef ModelCrbBase<ParameterDefinition,FunctionSpaceDefinition> super_type;
    typedef typename super_type::funs_type funs_type;
    typedef typename super_type::funsd_type funsd_type;
    typedef typename super_type::beta_vector_type beta_vector_type;
    typedef typename super_type::affine_decomposition_type affine_decomposition_type;
    //operator free
    typedef typename super_type::operator_type operator_type;
    typedef typename super_type::operator_ptrtype operator_ptrtype;
    typedef typename super_type::operatorcomposite_type operatorcomposite_type;
    typedef typename super_type::operatorcomposite_ptrtype operatorcomposite_ptrtype;
    typedef typename super_type::functionalcomposite_type functionalcomposite_type;
    typedef typename super_type::functionalcomposite_ptrtype functionalcomposite_ptrtype;
    typedef typename super_type::functional_type functional_type;
    typedef typename super_type::functional_ptrtype functional_ptrtype;
    using super_type::computeBetaQm;
    /** @name Constants
     */
    //@{
    //@}

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

    /* parameter space */
    typedef ParameterDefinition::parameterspace_type parameterspace_type;
    typedef boost::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef parameterspace_type::element_type parameter_type;

    /*reduced basis space*/
    typedef ReducedBasisSpace<self_type, mesh_type, basis_type, value_type> rbfunctionspace_type;
    typedef boost::shared_ptr< rbfunctionspace_type > rbfunctionspace_ptrtype;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    Heat1D();

    //! constructor from command line
    Heat1D( po::variables_map const& vm );


    //! copy constructor
    //Heat1D( Heat1D const & );
    //! destructor
    virtual ~Heat1D() {}

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

    /**
     * \brief compute the theta coefficient for both bilinear and linear form
     * \param mu parameter to evaluate the coefficients
     */
    boost::tuple<beta_vector_type, std::vector<beta_vector_type> >
    computeBetaQm( parameter_type const& mu, double time=0 )
    {
        M_betaAqm.resize( M_nb_terms_in_affine_decomposition_a );
        for(int q=0; q<M_nb_terms_in_affine_decomposition_a; q++)
            M_betaAqm[q].resize(1);

        M_betaAqm[0][0]= 1;
        M_betaAqm[1][0] = mu( 0 ); // k_1
        M_betaAqm[2][0] = mu( 1 ); // k_2

        M_betaFqm.resize( M_nb_outputs );
        M_betaFqm[0].resize( M_nb_terms_in_affine_decomposition_rhs );
        M_betaFqm[1].resize( M_nb_terms_in_affine_decomposition_output );
        M_betaFqm[0][0].resize( 1 );
        M_betaFqm[0][1].resize( 1 );
        M_betaFqm[1][0].resize( 1 );
        M_betaFqm[0][0][0] = mu( 2 ); // delta
        M_betaFqm[0][1][0] = mu( 3 ); // phi
        M_betaFqm[1][0][0] = 1;

        return boost::make_tuple( M_betaAqm, M_betaFqm );
    }


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * H1 scalar product
     */
    sparse_matrix_ptrtype innerProduct ( void )
    {
        return M;
    }

    //@}

    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false);

    /**
     * Important note
     * This model uses operators free so no need
     * to implement computeAffineDecomposition function
     * Instead, need to implement operatorCompositeA and functionalCompositeF
     */
    operatorcomposite_ptrtype operatorCompositeA()
    {
        return M_compositeA;
    }
    std::vector< functionalcomposite_ptrtype > functionalCompositeF()
    {
        return M_compositeF;
    }

private:

    po::variables_map M_vm;
    backend_ptrtype backend;

    double meshSize;

    bool M_use_weak_dirichlet;
    double M_gammabc;

    mesh_ptrtype mesh;
    space_ptrtype Xh;
    rbfunctionspace_ptrtype RbXh;
    sparse_matrix_ptrtype M;

    std::vector< std::vector<operator_ptrtype> > M_Aqm_free;
    std::vector< std::vector<std::vector<functional_ptrtype> > > M_Fqm_free;

    operatorcomposite_ptrtype M_compositeA;
    std::vector< functionalcomposite_ptrtype > M_compositeF;

    beta_vector_type M_betaAqm;
    std::vector<beta_vector_type> M_betaFqm;

    parameterspace_ptrtype M_Dmu;

    element_type u,v;

    int M_nb_terms_in_affine_decomposition_a=3;
    int M_nb_terms_in_affine_decomposition_rhs=2;
    int M_nb_terms_in_affine_decomposition_output=1;
    int M_nb_outputs=2;

};

Heat1D::Heat1D()
    :
    backend( backend_type::build( BACKEND_PETSC ) ),
    meshSize( 0.01 ),
    M_Dmu( new parameterspace_type )
{
}


Heat1D::Heat1D( po::variables_map const& vm )
    :
    M_vm( vm ),
    backend( backend_type::build( vm ) ),
    meshSize( vm["hsize"].as<double>() ),
    M_Dmu( new parameterspace_type )
{
}
void
Heat1D::initModel()
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

    //  initialisation de A1 et A2

    M_Aqm_free.resize( M_nb_terms_in_affine_decomposition_a );
    for(int q=0; q<M_nb_terms_in_affine_decomposition_a; q++)
        M_Aqm_free[q].resize( 1 );

    M_Fqm_free.resize( M_nb_outputs );
    M_Fqm_free[0].resize( M_nb_terms_in_affine_decomposition_rhs );
    for(int q=0; q<M_nb_terms_in_affine_decomposition_rhs ; q++)
        M_Fqm_free[0][q].resize(1);

    M_Fqm_free[1].resize( M_nb_terms_in_affine_decomposition_output );
    for(int q=0; q<M_nb_terms_in_affine_decomposition_output ; q++)
        M_Fqm_free[1][q].resize(1);

    using namespace Feel::vf;
    //static const int N = 2;
    Feel::ParameterSpace<4>::Element mu_min( M_Dmu );
    mu_min << 0.2, 0.2, 0.01, 0.1;
    M_Dmu->setMin( mu_min );
    Feel::ParameterSpace<4>::Element mu_max( M_Dmu );
    mu_max << 50, 50, 5, 5;
    M_Dmu->setMax( mu_max );

    u = Xh->element();
    v = Xh->element();

    LOG(INFO) << "Number of dof " << Xh->nLocalDof() << "\n";

    // right hand side
    auto f0 = integrate( markedfaces( mesh, "left" ), id( v ) );
    auto f1 = integrate( elements( mesh ), id( v ) );
    auto f0free = functionalLinearFree( _space=Xh , _expr=f0 , _backend=backend );
    auto f1free = functionalLinearFree( _space=Xh , _expr=f1 , _backend=backend );
    f0free->setName("f0");
    f1free->setName("f1");
    M_Fqm_free[0][0][0]=f0free;
    M_Fqm_free[0][1][0]=f1free;

    // output
    auto l0 = integrate( markedelements( mesh,"k1_2" ), id( v )/0.2 )
        + integrate( markedelements( mesh,"k2_1" ), id( v )/0.2 );
    auto l0free = functionalLinearFree( _space=Xh , _expr=l0 , _backend=backend );
    l0free->setName("l0");
    M_Fqm_free[1][0][0]=l0free;

    auto a0 = integrate( elements( mesh ), 0.1*( gradt( u )*trans( grad( v ) ) ) )
        + integrate( markedfaces( mesh,"right" ), id( u )*idt( v ) );
    auto a0free = opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=a0 , _backend=backend );
    a0free->setName("A0");
    M_Aqm_free[0][0]=a0free;

    auto a1 = integrate( markedelements( mesh, "k1_1"  ), ( gradt( u )*trans( grad( v ) ) ) )
        + integrate( markedelements( mesh,"k1_2"  ), ( gradt( u )*trans( grad( v ) ) ) );
    auto a1free = opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=a1 , _backend=backend );
    a1free->setName("A1");
    M_Aqm_free[1][0]=a1free;

    auto a2 = integrate( markedelements( mesh, "k2_1"  ), ( gradt( u )*trans( grad( v ) ) ) )
        + integrate( markedelements( mesh, "k2_2"  ), ( gradt( u )*trans( grad( v ) ) ) );
    auto a2free = opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=a2 , _backend=backend );
    a2free->setName("A2");
    M_Aqm_free[2][0]=a2free;

    M_compositeA = opLinearComposite( _domainSpace=Xh , _imageSpace=Xh );
    M_compositeA->addList( M_Aqm_free );
    M_compositeF.resize( M_nb_outputs );
    for(int output=0; output<M_nb_outputs; output++)
    {
        M_compositeF[output]=functionalLinearComposite( _space=Xh );
        M_compositeF[output]->addList( M_Fqm_free[output] );
    }


    M = backend->newMatrix( Xh, Xh );

    form2( _test=Xh, _trial=Xh, _matrix=M, _init=true ) =
        integrate( elements( mesh ), id( u )*idt( v ) + grad( u )*trans( gradt( u ) ) );
    M->close();


}


double
Heat1D::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve )
{

    CHECK( ! need_to_solve ) << "The model need to have the solution to compute the output\n";

    using namespace vf;

    double output=0;

    auto fqm = backend->newVector( Xh );
    // right hand side (compliant)
    if ( output_index == 0 )
    {
        for ( int q=0; q<M_nb_terms_in_affine_decomposition_rhs; q++ )
        {
            M_Fqm_free[output_index][q][0]->containerPtr( fqm );
            output += M_betaFqm[output_index][q][0]*dot( *fqm, u );
        }
    }

    // output
    if ( output_index == 1 )
    {
        //double mean = integrate( elements(mesh),
        //chi( (Px() >= -0.1) && (Px() <= 0.1) )*idv(*pT) ).evaluate()(0,0)/0.2;
        //std::cout<<"output1 c1 = "<<mean<<std::endl;

        output = ( integrate( markedelements( mesh,"k1_2" ),idv( u ) ).evaluate()( 0,0 )+
                   integrate( markedelements( mesh,"k2_1" ),idv( u ) ).evaluate()( 0,0 ) )/0.2;
        //std::cout<<"output1 c2 = "<<meanT<<std::endl;
        //std::cout<<"output1 c3= "<< dot( M_Fq[1][0], U ) <<std::endl;
        //return meanT;
    }

    return output;

}

}

#endif /* __Heat1D_H */


