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

#include <feel/feeldiscr/bdf.hpp>

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
    ( "steady", po::value<bool>()->default_value( false ), "if true then steady else unsteady" )
    ( "no-export", "don't export results" )
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
                           "Copyright (c) 2010,2011 Université de Grenoble 1 (Joseph Fourier)" );

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


    /** @name Constants
     */
    //@{

    static const uint16_type Order = 1;
    static const uint16_type ParameterSpaceDimension = 4;

    //@}

    /** @name Typedefs
     */
    //@{

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef backend_type::sparse_matrix_type sparse_matrix_type;
    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef backend_type::vector_type vector_type;
    typedef backend_type::vector_ptrtype vector_ptrtype;

    typedef Eigen::MatrixXd matrixN_type;


    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigen_matrix_type;
    typedef eigen_matrix_type ematrix_type;
    typedef boost::shared_ptr<eigen_matrix_type> eigen_matrix_ptrtype;

    /*mesh*/
    typedef Simplex<1,1> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef FunctionSpace<mesh_type, bases<Lagrange<0, Scalar> >, Discontinuous> p0_space_type;
    typedef p0_space_type::element_type p0_element_type;

    /*basis*/
    typedef bases<Lagrange<Order, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef space_type functionspace_type;
    typedef space_ptrtype functionspace_ptrtype;
    typedef space_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;

    /*reduced basis space*/
    typedef ReducedBasisSpace<super_type, mesh_type, basis_type, value_type> rbfunctionspace_type;
    typedef boost::shared_ptr< rbfunctionspace_type > rbfunctionspace_ptrtype;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;


    /* parameter space */
    typedef ParameterSpace<ParameterSpaceDimension> parameterspace_type;
    typedef boost::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef parameterspace_type::element_type parameter_type;
    typedef parameterspace_type::element_ptrtype parameter_ptrtype;
    typedef parameterspace_type::sampling_type sampling_type;
    typedef parameterspace_type::sampling_ptrtype sampling_ptrtype;

    /* time discretization */
    typedef Bdf<space_type>  bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;


    typedef std::vector< std::vector< double > > beta_vector_type;

    typedef boost::tuple<
        std::vector< std::vector<sparse_matrix_ptrtype> >,
        std::vector< std::vector<sparse_matrix_ptrtype> >,
        std::vector< std::vector<std::vector<vector_ptrtype> > >
        > affine_decomposition_type;

    typedef OperatorLinear< space_type , space_type > operator_type;
    typedef boost::shared_ptr<operator_type> operator_ptrtype;

    typedef OperatorLinearComposite< space_type , space_type > operatorcomposite_type;
    typedef boost::shared_ptr<operatorcomposite_type> operatorcomposite_ptrtype;

    typedef FsFunctionalLinearComposite< space_type > functionalcomposite_type;
    typedef boost::shared_ptr<functionalcomposite_type> functionalcomposite_ptrtype;

    typedef FsFunctionalLinear< space_type > functional_type;
    typedef boost::shared_ptr<functional_type> functional_ptrtype;

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

    // \return the number of terms in affine decomposition of left hand
    // side bilinear form
    int Qa() const
    {
        return 3;
    }

    // \return the number of terms in affine decomposition of bilinear form
    // associated to mass matrix
    int Qm() const
    {
        return 1;
    }

    // \return the number of terms in affine decomposition of left hand
    // side bilinear form ( adjoint model )
    int Qadu() const
    {
        return 3;
    }

    // \return the number of terms in affine decomposition of right hand
    //  linear form ( adjoint model )
    int Qfdu() const
    {
        return 1;
    }


    /**
     * there is at least one output which is the right hand side of the
     * primal problem
     *
     * \return number of outputs associated to the model
     */
    int Nl() const
    {
        return 2;
    }

    /**
     * \param l the index of output
     * \return number of terms  in affine decomposition of the \p q th output term
     */
    int Ql( int l ) const
    {
        if ( l == 0 ) return 2;

        return 1;
    }

    int mMaxA( int q )
    {
        if ( q < this->Qa() )
            return 1;
        else
            throw std::logic_error( "[Model heat1d] ERROR : try to acces to mMaxA(q) with a bad value of q");
    }

    int mMaxM( int q )
    {
        if ( q < this->Qm() )
            return 1;
        else
            throw std::logic_error( "[Model heat1d] ERROR : try to acces to mMaxM(q) with a bad value of q");
    }

    int mMaxF( int output_index, int q)
    {
        if ( q < this->Ql( output_index ) )
            return 1;
        else
            throw std::logic_error( "[Model heat1d] ERROR : try to acces to mMaxF(output_index,q) with a bad value of q");
    }

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
     * \brief compute the beta coefficient for both bilinear and linear form
     * \param mu parameter to evaluate the coefficients
     */


    boost::tuple<beta_vector_type, beta_vector_type, std::vector<beta_vector_type> >
    computeBetaQm( element_type const& T,parameter_type const& mu , double time=1e30 )
    {
        return computeBetaQm( mu , time );
    }

    boost::tuple<beta_vector_type, beta_vector_type, std::vector<beta_vector_type> >
    computeBetaQm( parameter_type const& mu , double time=0 )
    {
        M_betaAqm.resize( Qa() );
        M_betaAqm[0].resize(1);
        M_betaAqm[1].resize(1);
        M_betaAqm[2].resize(1);
        M_betaAqm[0][0] = 1;
        M_betaAqm[1][0] = mu( 0 ); // k_1
        M_betaAqm[2][0] = mu( 1 ); // k_2

        M_betaFqm.resize( Nl() );
        M_betaFqm[0].resize( Ql( 0 ) );
        M_betaFqm[0][0].resize(1);
        M_betaFqm[0][1].resize(1);
        M_betaFqm[0][0][0] = mu( 2 ); // delta
        M_betaFqm[0][1][0] = mu( 3 ); // phi

        M_betaFqm[1].resize( Ql( 1 ) );
        M_betaFqm[1][0].resize(1);
        M_betaFqm[1][0][0]= 1;

        M_betaMqm.resize( Qm() );
        M_betaMqm[0].resize(1);
        M_betaMqm[0][0] = 1;

        return boost::make_tuple( M_betaMqm, M_betaAqm, M_betaFqm );
    }


    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * set the mesh characteristic length to \p s
     */
    void setMeshSize( double s )
    {
        meshSize = s;
    }


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * run the convergence test
     */

    /**
     * create a new matrix
     * \return the newly created matrix
     */
    sparse_matrix_ptrtype newMatrix() const;
    vector_ptrtype newVector() const;
    /**
     * \brief Returns the affine decomposition
     */
    affine_decomposition_type computeAffineDecomposition();

    std::vector< std::vector< element_ptrtype > > computeInitialGuessAffineDecomposition();
    /**
     * \brief solve the model for parameter \p mu
     * \param mu the model parameter
     * \param T the temperature field
     * \param output_index (optional) useful for the adjoint ( POD )
     */
    void solve( parameter_type const& mu, element_ptrtype& T , int output_index=0 );

    //eigen_matrix_type snapshotsMatrix(){return M_snapshots_matrix;}
    //eigen_matrix_type dualSnapshotsMatrix(){return Mdu_snapshots_matrix;}
    //void fillSnapshotsMatrix (parameter_type const& mu, int output);
    void assemble();
    int computeNumberOfSnapshots();
    double timeStep()
    {
        return  M_bdf->timeStep();
    }
    double timeFinal()
    {
        return M_bdf->timeFinal();
    }
    double timeInitial()
    {
        return M_bdf->timeInitial();
    }
    int timeOrder()
    {
        return M_bdf->timeOrder();
    }
    bool isSteady()
    {
        return M_is_steady;
    }
    void initializationField( element_ptrtype& initial_field,parameter_type const& mu );

    /**
     * solve for a given parameter \p mu
     */
    element_type solve( parameter_type const& mu );

    /**
     * solve \f$ M u = f \f$
     */
    void l2solve( vector_ptrtype& u, vector_ptrtype const& f );


    /**
     * update the PDE system with respect to \param mu
     */
    void update( parameter_type const& mu, double bdf_coeff, element_type const& bdf_poly, int output_index=0 );
    //@}

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( double time, element_type& T );

    /**
     * H1 scalar product
     */
    sparse_matrix_ptrtype innerProduct ( void )
    {
        return M;
    }


    void exportResults1d( double time, element_type& T, double s );

    /**
     * returns the scalar product of the boost::shared_ptr vector x and
     * boost::shared_ptr vector y
     */
    double scalarProduct( vector_ptrtype const& X, vector_ptrtype const& Y );

    /**
     * returns the scalar product of the vector x and vector y
     */
    double scalarProduct( vector_type const& x, vector_type const& y );

    /**
     * returns the scalar product used for POD of the boost::shared_ptr vector x and
     * boost::shared_ptr vector y
     */
    double scalarProductForPod( vector_ptrtype const& X, vector_ptrtype const& Y );

    /**
     * returns the scalar product used for POD of the vector x and vector y
     */
    double scalarProductForPod( vector_type const& x, vector_type const& y );

    /**
     * specific interface for OpenTURNS
     *
     * \param X input vector of size N
     * \param N size of input vector X
     * \param Y input vector of size P
     * \param P size of input vector Y
     */
    void run( const double * X, unsigned long N, double * Y, unsigned long P );

    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false, bool export_output=false );

    static const bool is_time_dependent = true;

    operatorcomposite_ptrtype operatorCompositeA()
    {
        return M_compositeA;
    }
    operatorcomposite_ptrtype operatorCompositeM()
    {
        return M_compositeM;
    }
    std::vector< functionalcomposite_ptrtype > functionalCompositeF()
    {
        return M_compositeF;
    }

    parameter_type refParameter()
    {
        return M_Dmu->min();
    }

private:

    po::variables_map M_vm;

    backend_ptrtype backend;

    bool M_is_steady ;

    double alpha;

    double meshSize;

    parameterspace_ptrtype M_Dmu;

    export_ptrtype exporter;

    mesh_ptrtype mesh;
    space_ptrtype Xh;
    rbfunctionspace_ptrtype RbXh;
    sparse_matrix_ptrtype D,M,Mpod;
    vector_ptrtype F;
    element_ptrtype pT;

    std::vector< std::vector<sparse_matrix_ptrtype> > M_Aqm;
    std::vector< std::vector<sparse_matrix_ptrtype> > M_Mqm;
    std::vector< std::vector<std::vector<vector_ptrtype> > > M_Fqm;

    std::vector< std::vector<operator_ptrtype> > M_Aqm_free;
    std::vector< std::vector<operator_ptrtype> > M_Mqm_free;
    std::vector< std::vector<std::vector<functional_ptrtype> > > M_Fqm_free;
    operatorcomposite_ptrtype M_compositeA;
    operatorcomposite_ptrtype M_compositeM;
    std::vector< functionalcomposite_ptrtype > M_compositeF;

    beta_vector_type M_betaAqm;
    beta_vector_type M_betaMqm;
    std::vector<beta_vector_type> M_betaFqm;

    bdf_ptrtype M_bdf;

    element_type u,v;
};

UnsteadyHeat1D::UnsteadyHeat1D()
    :
    backend( backend_type::build( BACKEND_PETSC ) ),
    M_is_steady( false ),
    alpha( 1 ),
    meshSize( 0.01 ),
    M_Dmu( new parameterspace_type )
{
}


UnsteadyHeat1D::UnsteadyHeat1D( po::variables_map const& vm )
    :
    M_vm( vm ),
    backend( backend_type::build( vm ) ),
    M_is_steady( vm["steady"].as<bool>() ),
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

    // allocate an element of Xh
    pT = element_ptrtype( new element_type( Xh ) );

    //M_bdf = bdf_ptrtype( new bdf_type(M_vm ,Xh, " " ) );
    M_bdf = bdf( _space=Xh, _vm=M_vm, _name="unsteadyHeat1d" , _prefix="unsteadyHeat1d" );

    M_Aqm_free.resize( 3 );
    for(int i=0; i<Qa(); i++) M_Aqm_free[i].resize(1);

    M_Fqm_free.resize( 2 );
    M_Fqm_free[0].resize( 2 );
    M_Fqm_free[1].resize( 1 );

    for(int i=0; i<M_Fqm_free[0].size(); i++) M_Fqm_free[0][i].resize(1);
    for(int i=0; i<M_Fqm_free[1].size(); i++) M_Fqm_free[1][i].resize(1);

    M_Mqm_free.resize( 1 );
    for(int i=0; i<M_Mqm_free.size(); i++) M_Mqm_free[i].resize(1);

    D = backend->newMatrix( Xh, Xh );
    F = backend->newVector( Xh );

    u = Xh->element();
    v = Xh->element();

    using namespace Feel::vf;

    Feel::ParameterSpace<4>::Element mu_min( M_Dmu );
    mu_min << 0.2, 0.2, 0.01, 0.1;
    M_Dmu->setMin( mu_min );
    Feel::ParameterSpace<4>::Element mu_max( M_Dmu );
    mu_max << 50, 50, 5, 5;
    M_Dmu->setMax( mu_max );

    Mpod = backend->newMatrix( Xh, Xh );
    form2( Xh, Xh, Mpod, _init=true ) =
        //integrate( elements(mesh), id(u)*idt(v) );
        integrate( elements( mesh ), id( u )*idt( v ) + grad( v )*trans( gradt( u ) ) );
    Mpod->close();


    M = backend->newMatrix( Xh, Xh );
    form2( Xh, Xh, M, _init=true ) =
        integrate( elements( mesh ), id( u )*idt( v ) + grad( v )*trans( gradt( u ) ) );
    //integrate( elements(mesh), id(u)*idt(v) );
    M->close();

    LOG(INFO) << "Number of dof " << Xh->nLocalDof() << "\n";

    assemble();

} // UnsteadyHeat1d::init


void
UnsteadyHeat1D::initializationField( element_ptrtype& initial_field,parameter_type const& mu )
{
    initial_field->zero();
}

void
UnsteadyHeat1D::assemble()
{

    using namespace Feel::vf;

    u = Xh->element();
    v = Xh->element();


    //mass matrix
    auto expr_mass = integrate ( elements( mesh ), alpha*idt( u )*id( v ) );
    auto operator_mass=opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=expr_mass );
    operator_mass->setName("mass operator");
    M_Mqm_free[0][0]=operator_mass;

    // right hand side
    auto expr_f0 = integrate( markedfaces( mesh,"left" ), id( v ) );
    auto f0 = functionalLinearFree( _space=Xh , _expr=expr_f0  );
    f0->setName("F0");
    M_Fqm_free[0][0][0]=f0;
    auto expr_f1 = integrate( elements( mesh ), id( v ) );
    auto f1 = functionalLinearFree( _space=Xh , _expr=expr_f1  );
    f1->setName("F1");
    M_Fqm_free[0][1][0]=f1;

    // output
    auto expr_l0 = integrate( markedelements( mesh,"k1_2" ), id( v )/0.2 ) + integrate( markedelements( mesh,"k2_1" ), id( v )/0.2 );
    auto l0 = functionalLinearFree( _space=Xh , _expr=expr_l0  );
    l0->setName("L0");
    M_Fqm_free[1][0][0]=l0;

    auto expr_A0 = integrate( elements( mesh ), 0.1*( gradt( u )*trans( grad( v ) ) ) ) + integrate( markedfaces( mesh,"right" ), idt( u )*id( v ) );
    auto operator_A0=opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=expr_A0 );
    operator_A0->setName("A0");
    M_Aqm_free[0][0]=operator_A0;

    auto expr_A1 = integrate( markedelements( mesh,"k1_1" ), ( gradt( u )*trans( grad( v ) ) ) )
        + integrate( markedelements( mesh,"k1_2" ), ( gradt( u )*trans( grad( v ) ) ) );
    auto operator_A1=opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=expr_A1 );
    operator_A1->setName("A1");
    M_Aqm_free[1][0]=operator_A1;

    auto expr_A2 = integrate( markedelements( mesh,"k2_1" ), ( gradt( u )*trans( grad( v ) ) ) )
        + integrate( markedelements( mesh,"k2_2" ), ( gradt( u )*trans( grad( v ) ) ) );
    auto operator_A2=opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=expr_A2 );
    operator_A2->setName("A2");
    M_Aqm_free[2][0]=operator_A2;

    M_compositeA = opLinearComposite( _domainSpace=Xh , _imageSpace=Xh );
    M_compositeM = opLinearComposite( _domainSpace=Xh , _imageSpace=Xh );
    M_compositeA->addList( M_Aqm_free );
    M_compositeM->addList( M_Mqm_free );
    M_compositeF.resize( this->Nl() );
    for(int output=0; output<this->Nl(); output++)
    {
        M_compositeF[output]=functionalLinearComposite( _space=Xh );
        M_compositeF[output]->addList( M_Fqm_free[output] );
    }

}




UnsteadyHeat1D::sparse_matrix_ptrtype
UnsteadyHeat1D::newMatrix() const
{
    return backend->newMatrix( Xh, Xh );
}


typename UnsteadyHeat1D::vector_ptrtype
UnsteadyHeat1D::newVector() const
{
    return backend->newVector( Xh );
}


std::vector< std::vector< UnsteadyHeat1D::element_ptrtype > >
UnsteadyHeat1D::computeInitialGuessAffineDecomposition()
{
    std::vector< std::vector<element_ptrtype> > q;
    q.resize(1);
    q[0].resize(1);
    element_ptrtype elt ( new element_type ( Xh ) );
    q[0][0] = elt;
    return q;
}

UnsteadyHeat1D::affine_decomposition_type
UnsteadyHeat1D::computeAffineDecomposition()
{
    return boost::make_tuple( M_Mqm, M_Aqm, M_Fqm );
}


void
UnsteadyHeat1D::exportResults( double time, element_type& T )
{
#if 0

    if ( M_do_export )
    {
        LOG(INFO) << "exportResults starts\n";

        exporter->step( time )->setMesh( T.functionSpace()->mesh() );

        exporter->step( time )->add( "T", T );

        exporter->save();
    }

#endif
} // Heat1d::export


void
UnsteadyHeat1D::exportResults1d( double time, element_type& T, double output )
{
    std::map<double,double> data;

    for ( auto it = T.functionSpace()->mesh()->beginElement( ),
            en = T.functionSpace()->mesh()->endElement( );
            it!=en; ++it )
    {
        for ( size_type i = 0; i < space_type::basis_type::nLocalDof; ++i )
        {
            value_type a = it->point( 0 ).node()[0];
            value_type b = it->point( 1 ).node()[0];
            value_type x = 0;

            if ( i == 0 )
                x=a;

            else if ( i == 1 )
                x=b;

            else
                x= a + ( i-1 )*( b-a )/( space_type::basis_type::nLocalDof-1 );

            data[x] = T.localToGlobal( it->id(), i, 0 );
        }
    }

    std::ostringstream fname_T;
    std::ofstream ofs3( ( boost::format( "Solution_%.3f" ) % time ).str().c_str() );
    BOOST_FOREACH( auto d, data )
    {
        ofs3 <<std::setprecision( 16 )<< d.first << " " << d.second << "\n";
    }
    ofs3.close();


    std::ofstream output_file;
    output_file.open( "Output.dat",std::ios::out | std::ios::app );
    output_file<<std::setprecision( 16 )<<time<<" "<<output<<"\n";

}





void
UnsteadyHeat1D::update( parameter_type const& mu,double bdf_coeff, element_type const& bdf_poly, int output_index )
{
    if (option(_name="crb.stock-matrices"). as<bool>() )
    {

        D->close();
        D->zero();

        for ( size_type q = 0; q < M_Aqm.size(); ++q )
        {
            for ( size_type m = 0; m < mMaxA(q); ++m )
                D->addMatrix( M_betaAqm[q][m] , M_Aqm[q][m] );
        }

        F->close();
        F->zero();

        for ( size_type q = 0; q < M_Fqm[output_index].size(); ++q )
        {
            for ( size_type m = 0; m < mMaxF(output_index,q); ++m )
            F->add( M_betaFqm[output_index][q][m], M_Fqm[output_index][q][m] );
        }

        auto vec_bdf_poly = backend->newVector( Xh );

        //add contribution from mass matrix
        for ( size_type q = 0; q < Qm(); ++q )
        {
            for ( size_type m = 0; m < mMaxM(q); ++m )
            {
                //left hand side
                D->addMatrix( M_betaMqm[q][m]*bdf_coeff, M_Mqm[q][m] );
                //right hand side
                *vec_bdf_poly = bdf_poly;
                vec_bdf_poly->close();
                vec_bdf_poly->scale( M_betaMqm[q][m] );
                F->addVector( *vec_bdf_poly, *M_Mqm[q][m] );
            }
        }
    }
    else
    {
        D->close();
        D->zero();
        F->close();
        F->zero();

        M_compositeA->setScalars( M_betaAqm );
        M_compositeA->sumAllMatrices( D );

        M_compositeF[output_index]->setScalars( M_betaFqm[output_index] );
        M_compositeF[output_index]->sumAllVectors( F );

        auto vec_bdf_poly = backend->newVector( Xh );

        for ( size_type q = 0; q < Qm(); ++q )
        {
            for ( size_type m = 0; m < mMaxM(q); ++m )
            {
                auto matrix = backend->newMatrix( _test=Xh , _trial=Xh );
                M_compositeM->operatorlinear(q,m)->matPtr( matrix );
                //left hand side
                D->addMatrix( M_betaMqm[q][m]*bdf_coeff, matrix );
                //right hand side
                *vec_bdf_poly = bdf_poly;
                vec_bdf_poly->close();
                vec_bdf_poly->scale( M_betaMqm[q][m] );
                F->addVector( *vec_bdf_poly, *matrix );
            }
        }
    }
}



typename UnsteadyHeat1D::element_type
UnsteadyHeat1D::solve( parameter_type const& mu )
{
    element_ptrtype T( new element_type( Xh ) );
    this->solve( mu, T );
    return *T;
}

void
UnsteadyHeat1D::solve( parameter_type const& mu, element_ptrtype& T , int output_index )
{

    using namespace vf;
    element_type v( Xh, "v" );//test functions

    pT->zero();
    T->zero();

    assemble();

    this->computeBetaQm( mu );

    if ( M_is_steady )
        M_bdf->setSteady();

    for ( M_bdf->start(*T); !M_bdf->isFinished(); M_bdf->next() )
    {
        double bdf_coeff = M_bdf->polyDerivCoefficient( 0 );

        auto bdf_poly = M_bdf->polyDeriv();
        this->update( mu , bdf_coeff, bdf_poly );
        auto ret = backend->solve( _matrix=D,  _solution=T, _rhs=F , _reuse_prec=( M_bdf->iteration() >=2 ) );
        if ( !ret.get<0>() )
            LOG(INFO)<<"WARNING : at time "<<M_bdf->time()<<" we have not converged ( nb_it : "<<ret.get<1>()<<" and residual : "<<ret.get<2>() <<" ) \n";

        M_bdf->shiftRight( *T );
    }

}

void
UnsteadyHeat1D::l2solve( vector_ptrtype& u, vector_ptrtype const& f )
{
    //std::cout << "l2solve(u,f)\n";
    backend->solve( _matrix=M,  _solution=u, _rhs=f );
    //std::cout << "l2solve(u,f) done\n";
}


double
UnsteadyHeat1D::scalarProduct( vector_ptrtype const& x, vector_ptrtype const& y )
{
    return M->energy( x, y );
}
double
UnsteadyHeat1D::scalarProduct( vector_type const& x, vector_type const& y )
{
    return M->energy( x, y );
}


double
UnsteadyHeat1D::scalarProductForPod( vector_ptrtype const& x, vector_ptrtype const& y )
{
    return Mpod->energy( x, y );
}
double
UnsteadyHeat1D::scalarProductForPod( vector_type const& x, vector_type const& y )
{
    return Mpod->energy( x, y );
}



void
UnsteadyHeat1D::run( const double * X, unsigned long N, double * Y, unsigned long P )
{
    using namespace vf;
    Feel::ParameterSpace<4>::Element mu( M_Dmu );
    mu << X[0], X[1], X[2], X[3];
    static int do_init = true;

    if ( do_init )
    {
        meshSize = X[4];
        this->initModel();
        do_init = false;
    }

    this->solve( mu, pT );

    double mean = integrate( elements( mesh ),
                             chi( ( Px() >= -0.1 ) && ( Px() <= 0.1 ) )*idv( *pT ) ).evaluate()( 0,0 )/0.2;
    Y[0]=mean;
}



//if we export output we call it in solve function so we don't need to recall it
double
UnsteadyHeat1D::output( int output_index, parameter_type const& mu, element_type& u, bool need_to_solve, bool export_output )
{

    using namespace vf;

    if ( need_to_solve )
        this->solve( mu, pT );
    else
        *pT = u;

    using namespace vf;

    //vector_ptrtype U( backend->newVector( Xh ) );
    //*U = *pT;

    double output=0;
    // right hand side (compliant)
    if ( output_index == 0 )
    {
        for ( int q=0; q<Ql( output_index ); q++ )
        {
            for ( int m=0; m<mMaxF(output_index,q); m++ )
            {
                //element_ptrtype eltF( new element_type( Xh ) );
                //*eltF = *M_Fqm[output_index][q][m];
                //output += M_betaFqm[output_index][q][m]*dot( *eltF, *pT );
                output += M_betaFqm[output_index][q][m]*dot( *M_Fqm[output_index][q][m], *pT );
            }
        }

    }
    // output
    else if ( output_index == 1 )
    {
        output = integrate( elements( mesh ),
                            chi( ( Px() >= -0.1 ) && ( Px() <= 0.1 ) )*idv( *pT ) ).evaluate()( 0,0 )/0.2;
    }
    else
        throw std::logic_error( "[unsteadyHeat1d::output] error with output_index : only 0 or 1 " );

    return output;
}

}

#endif /* __Heat1D_H */


