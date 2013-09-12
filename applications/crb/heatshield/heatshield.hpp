/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-11-13

  Copyright (C) 2009 Universite Joseph Fourier (Grenoble I)

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
   \file heatshield.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-03-28
 */
#ifndef __HeatShield_H
#define __HeatShield_H 1

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

#include <feel/feeldiscr/bdf2.hpp>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>

#include <feel/feelcore/pslogger.hpp>



namespace Feel
{

po::options_description
makeHeatShieldOptions()
{
    po::options_description heatshieldoptions( "HeatShield options" );
    heatshieldoptions.add_options()
    // mesh parameters
    ( "hsize", Feel::po::value<double>()->default_value( 1e-1 ), "first h value to start convergence" )
    ( "steady", Feel::po::value<bool>()->default_value( true ), "if true : steady else unsteady" )
    ( "mshfile", Feel::po::value<std::string>()->default_value( "" ), "name of the gmsh file input")
    ( "do-export", Feel::po::value<bool>()->default_value( false ), "export results if true" )
    ;
    return heatshieldoptions.add( Feel::feel_options() ).add( bdf_options( "heatshield" ) );
}
AboutData
makeHeatShieldAbout( std::string const& str = "heatShield" )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "heat shield Benchmark",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2010,2011 Universite de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Stephane Veys", "developer", "stephane.veys@imag.fr", "" );
    return about;
}

class ParameterDefinition
{
public :
    static const uint16_type ParameterSpaceDimension = 2;
    typedef ParameterSpace<ParameterSpaceDimension> parameterspace_type;
};
class FunctionSpaceDefinition
{
public :
    static const uint16_type Order = 1;
    typedef double value_type;

    /*mesh*/
    typedef Simplex<2,Order> entity_type;
    typedef Mesh<entity_type> mesh_type;

    /*basis*/
    typedef fusion::vector<Lagrange<Order, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
};


/**
 * \class HeatShield
 * \brief brief description
 *
 * @author Christophe Prud'homme
 * @see
 */
class HeatShield : public ModelCrbBase< ParameterDefinition, FunctionSpaceDefinition >
{
public:

    typedef ModelCrbBase<ParameterDefinition, FunctionSpaceDefinition> super_type;
    typedef typename super_type::funs_type funs_type;
    typedef typename super_type::funsd_type funsd_type;


    /** @name Constants
     */
    //@{

    static const uint16_type Order = 1;
    static const uint16_type ParameterSpaceDimension = 2;
    static const bool is_time_dependent = true;

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

    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigen_matrix_type;
    typedef eigen_matrix_type ematrix_type;
    typedef boost::shared_ptr<eigen_matrix_type> eigen_matrix_ptrtype;

    /*mesh*/
    typedef Simplex<2,Order> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef FunctionSpace<mesh_type, fusion::vector<Lagrange<0, Scalar> >, Discontinuous> p0_space_type;
    typedef typename p0_space_type::element_type p0_element_type;

    /*basis*/
    typedef fusion::vector<Lagrange<Order, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef space_type functionspace_type;
    typedef space_ptrtype functionspace_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;

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
    HeatShield();

    //! constructor from command line
    HeatShield( po::variables_map const& vm );


    //! copy constructor
    HeatShield( HeatShield const & );
    //! destructor
    ~HeatShield() {}

    //! initialization of the model
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

    /**
     * there is at least one output which is the right hand side of the
     * primal problem
     *
     * \return number of outputs associated to the model
     * in our case we have a compliant output and 2 others outputs : average temperature on boundaries
     */
    int Nl() const
    {
        return 2;
    }

    /**
     * \param l the index of output
     * \return number of terms  in affine decomposition of the \p q th output term
     * in our case no outputs depend on parameters
     */
    int Ql( int l ) const
    {
        return 1;
    }

    int mMaxA( int q )
    {
        if ( q < 3 )
            return 1;
        else
            throw std::logic_error( "[Model heatshield] ERROR : try to acces to mMaxA(q) with a bad value of q");
    }

    int mMaxM( int q )
    {
        if ( q < 1 )
            return 1;
        else
            throw std::logic_error( "[Model heatshield] ERROR : try to acces to mMaxM(q) with a bad value of q");
    }

    int mMaxF( int output_index, int q)
    {
        if ( q < 1 )
            return 1;
        else
            throw std::logic_error( "[Model heatshield] ERROR : try to acces to mMaxF(output_index,q) with a bad value of q");
    }


    /**
     * \brief Returns the function space
     */
    space_ptrtype functionSpace()
    {
        return Xh;
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
    boost::tuple<beta_vector_type, beta_vector_type, std::vector<beta_vector_type> >
    computeBetaQm( element_type const& T,parameter_type const& mu , double time=1e30 )
    {
        return computeBetaQm( mu , time );
    }

    boost::tuple<beta_vector_type, beta_vector_type, std::vector<beta_vector_type>  >
    computeBetaQm( parameter_type const& mu , double time=1e30 )
    {
        double biot_out   = mu( 0 );
        double biot_in    = mu( 1 );

        M_betaAqm.resize( Qa() );
        M_betaAqm[0].resize( 1 );
        M_betaAqm[1].resize( 1 );
        M_betaAqm[2].resize( 1 );
        M_betaAqm[0][0] = 1 ;
        M_betaAqm[1][0] = biot_out ;
        M_betaAqm[2][0] = biot_in  ;

        M_betaMqm.resize( Qm() );
        M_betaMqm[0].resize( 1 );
        M_betaMqm[0][0] = 1;

        M_betaFqm.resize( Nl() );
        M_betaFqm[0].resize( Ql(0) );
        M_betaFqm[0][0].resize( 1 );
        M_betaFqm[0][0][0] = biot_out;

        M_betaFqm[1].resize( Ql(1) );
        M_betaFqm[1][0].resize( 1 );
        M_betaFqm[1][0][0] = 1./surface;

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

    /**
     * \return the newly created vector
     */
    vector_ptrtype newVector() const;

    /**
     * \brief Returns the affine decomposition
     */
    affine_decomposition_type computeAffineDecomposition();

    void stockAffineDecomposition();

    std::vector< std::vector< element_ptrtype > > computeInitialGuessAffineDecomposition()
    {
        std::vector< std::vector<element_ptrtype> > q;
        q.resize(1);
        q[0].resize(1);
        element_ptrtype elt ( new element_type ( Xh ) );
        q[0][0] = elt;
        return q;
    }
    /**
     * \brief solve the model for parameter \p mu
     * \param mu the model parameter
     * \param T the temperature field
     */
    void solve( parameter_type const& mu, element_ptrtype& T, int output_index=0 );

    void assemble();
    int computeNumberOfSnapshots();
    double timeFinal()
    {
        return M_bdf->timeFinal();
    }
    double timeStep()
    {
        return  M_bdf->timeStep();
    }
    double timeInitial()
    {
        return M_bdf->timeInitial();
    }
    int timeOrder()
    {
        return M_bdf->timeOrder();
    }
    void initializationField( element_ptrtype& initial_field, parameter_type const& mu );
    bool isSteady()
    {
        return M_is_steady;
    }


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
    void update( parameter_type const& mu,double bdf_coeff, element_type const& bdf_poly, int output_index=0 ) ;
    //@}

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( double time, element_type& T , parameter_type const& mu );

    /**
     * inner product
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

    void solve( sparse_matrix_ptrtype& ,element_type& ,vector_ptrtype&  );

     /**
     * returns the scalar product used fior mass matrix ( to solve eigen values problem )
     * of the boost::shared_ptr vector x and boost::shared_ptr vector
     */
    double scalarProductForMassMatrix( vector_ptrtype const& X, vector_ptrtype const& Y );

    /**
     * returns the scalar product for mass matrix of the vector x and vector y
     */
    double scalarProductForMassMatrix( vector_type const& x, vector_type const& y );

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
     * returns the scalar product used in POD of the boost::shared_ptr vector x and
     * boost::shared_ptr vector y
     */
    double scalarProductForPod( vector_ptrtype const& X, vector_ptrtype const& Y );

    /**
     * returns the scalar product used in POD of the vector x and vector y
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
    value_type output( int output_index, parameter_type const& mu, element_type &T, bool need_to_solve=false, bool export_outputs=false );

    gmsh_ptrtype createGeo( double hsize );

    virtual operatorcomposite_ptrtype operatorCompositeA()
    {
        return M_compositeA;
    }
    virtual operatorcomposite_ptrtype operatorCompositeM()
    {
        return M_compositeM;
    }
    virtual std::vector< functionalcomposite_ptrtype > functionalCompositeF()
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
    backend_ptrtype backendl2;
    bool M_is_steady ;

    /* mesh parameters */
    double meshSize;

    int export_number;

    bool do_export;

    parameterspace_ptrtype M_Dmu;

    double surface;

    /* mesh, pointers and spaces */
    mesh_ptrtype mesh;
    space_ptrtype Xh;

    sparse_matrix_ptrtype D,M,Mpod,InnerMassMatrix;
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
HeatShield::HeatShield()
    :
    backend( backend_type::build( BACKEND_PETSC ) ),
    backendl2( backend_type::build( BACKEND_PETSC ) ),
    M_is_steady( false ),
    meshSize( 2e-1 ),
    export_number( 0 ),
    do_export( false ),
    M_Dmu( new parameterspace_type )
{ }

HeatShield::HeatShield( po::variables_map const& vm )
    :
    M_vm( vm ),
    backend( backend_type::build( vm ) ),
    backendl2( backend_type::build( vm ) ),
    M_is_steady( vm["steady"].as<bool>() ),
    meshSize( vm["hsize"].as<double>() ),
    export_number( 0 ),
    do_export( vm["do-export"].as<bool>() ),
    M_Dmu( new parameterspace_type )
{ }





gmsh_ptrtype
HeatShield::createGeo( double hsize )
{
    gmsh_ptrtype gmshp( new Gmsh );
    std::ostringstream ostr;
    double H = hsize;
    double h = hsize*0.5;
    //double h = hsize*1;
    ostr <<"Point (1) = {0,  0, 0, "<<H<<"};\n"
         <<"Point (2) = {10, 0, 0, "<<H<<"};\n"
         <<"Point (3) = {10, 4, 0, "<<H<<"};\n"
         <<"Point (4) = {0,  4, 0, "<<H<<"};\n"
         <<"Point (10) = {1, 1, 0, "<<h<<"};\n"
         <<"Point (11) = {3, 1, 0, "<<h<<"};\n"
         <<"Point (12) = {3, 3, 0, "<<h<<"};\n"
         <<"Point (13) = {1, 3, 0, "<<h<<"};\n"
         <<"Point (20) = {4, 1, 0, "<<h<<"};\n"
         <<"Point (21) = {6, 1, 0, "<<h<<"};\n"
         <<"Point (22) = {6, 3, 0, "<<h<<"};\n"
         <<"Point (23) = {4, 3, 0, "<<h<<"};\n"
         <<"Point (30) = {7, 1, 0, "<<h<<"};\n"
         <<"Point (31) = {9, 1, 0, "<<h<<"};\n"
         <<"Point (32) = {9, 3, 0, "<<h<<"};\n"
         <<"Point (33) = {7, 3, 0, "<<h<<"};\n"
         <<"Line (101) = {1,2};\n"
         <<"Line (102) = {2,3};\n"
         <<"Line (103) = {3,4};\n"
         <<"Line (104) = {4,1};\n"
         <<"Line (110) = {10,11};\n"
         <<"Line (111) = {11,12};\n"
         <<"Line (112) = {12,13};\n"
         <<"Line (113) = {13,10};\n"
         <<"Line (120) = {20,21};\n"
         <<"Line (121) = {21,22};\n"
         <<"Line (122) = {22,23};\n"
         <<"Line (123) = {23,20};\n"
         <<"Line (130) = {30,31};\n"
         <<"Line (131) = {31,32};\n"
         <<"Line (132) = {32,33};\n"
         <<"Line (133) = {33,30};\n"
         <<"Line Loop (201) = {101, 102, 103, 104};\n"
         <<"Line Loop (210) = {110, 111, 112, 113};\n"
         <<"Line Loop (220) = {120, 121, 122, 123};\n"
         <<"Line Loop (230) = {130, 131, 132, 133};\n"
         <<"Plane Surface (300) = {201,-210,-220,-230};\n"
         <<"Physical Line (\"left\") = {104};\n"
         <<"Physical Line (\"right\") = {102};\n"
         <<"Physical Line (\"bottom\") = {101};\n"
         <<"Physical Line (\"top\") = {103};\n"
         <<"Physical Line (\"gamma_holes\") = {110,111,112,113, 120,121,122,123, 130,131,132,133};\n"
         <<"Physical Surface (\"Omega\") = {300};\n"
         ;
    std::ostringstream nameStr;
    nameStr.precision( 3 );
    nameStr << "heatshield_geo";
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( ostr.str() );
    return gmshp;
}




void HeatShield::initializationField( element_ptrtype& initial_field , parameter_type const& mu )
{
    initial_field->setZero();
}

void HeatShield::initModel()
{


    using namespace Feel::vf;

    std::string mshfile_name = option("mshfile").as<std::string>();

    /*
     * First we create the mesh or load it if already exist
     */

    if( mshfile_name=="" )
    {
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                               _desc = createGeo( meshSize ) );
    }
    else
    {
        mesh = loadGMSHMesh( _mesh=new mesh_type,
                             _filename=option("mshfile").as<std::string>(),
                             _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER );
    }


    /*
     * The function space and some associate elements are then defined
     */
    Xh = space_type::New( mesh );
    std::cout << "Number of dof " << Xh->nLocalDof() << "\n";
    LOG(INFO) << "Number of dof " << Xh->nLocalDof() << "\n";

    // allocate an element of Xh
    pT = element_ptrtype( new element_type( Xh ) );

    surface = integrate( _range=elements( mesh ), _expr=cst( 1. ) ).evaluate()( 0,0 );
    //std::cout<<"surface : "<<surface<<std::endl;

    M_bdf = bdf( _space=Xh, _vm=M_vm, _name="heatshield" , _prefix="heatshield" );


    M_Aqm_free.resize( this->Qa() );
    for(int q=0; q<Qa(); q++)
        M_Aqm_free[q].resize( 1 );

    M_Mqm_free.resize( this->Qm() );
    for(int q=0; q<Qm(); q++)
        M_Mqm_free[q].resize( 1 );

    M_Fqm_free.resize( this->Nl() );
    for(int l=0; l<Nl(); l++)
    {
        M_Fqm_free[l].resize( Ql(l) );
        for(int q=0; q<Ql(l) ; q++)
        {
            M_Fqm_free[l][q].resize(1);
        }
    }

    Feel::ParameterSpace<ParameterSpaceDimension>::Element mu_min( M_Dmu );
    mu_min <<  /* Bi_out */ 1e-2 , /*Bi_in*/1e-3;
    M_Dmu->setMin( mu_min );
    Feel::ParameterSpace<ParameterSpaceDimension>::Element mu_max( M_Dmu );
    mu_max << /* Bi_out*/0.5   ,  /*Bi_in*/0.1;
    M_Dmu->setMax( mu_max );

    LOG(INFO) << "Number of dof " << Xh->nLocalDof() << "\n";

    assemble();
    PsLogger ps("ps-Model");
    ps.log("after assemble");
    if (option(_name="crb.stock-matrices"). as<bool>() )
        stockAffineDecomposition();
    ps.log("after stocking matrices");

} // HeatShield::init



void HeatShield::assemble()
{

    using namespace Feel::vf;

    u = Xh->element();
    v = Xh->element();

    auto expr_a00 = integrate( _range= elements( mesh ), _expr= gradt( u )*trans( grad( v ) ) );
    auto operatorfree00=opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=expr_a00 , _backend=backend );
    operatorfree00->setName("A0");
    M_Aqm_free[0][0]=operatorfree00;

    auto expr_a10  = integrate( _range= markedfaces( mesh, "left" ), _expr= idt( u )*id( v ) );
    auto operatorfree10=opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=expr_a10 , _backend=backend );
    operatorfree10->setName("A1");
    M_Aqm_free[1][0]=operatorfree10;

    auto expr_a20  = integrate( _range= markedfaces( mesh, "gamma_holes" ), _expr= idt( u )*id( v ) );
    auto operatorfree20=opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=expr_a20 , _backend=backend );
    operatorfree20->setName("A2");
    M_Aqm_free[2][0]=operatorfree20;

    auto expr_f000 = integrate( _range=markedfaces( mesh,"left" ), _expr= id( v ) ) ;
    auto functionalfree000 = functionalLinearFree( _space=Xh , _expr=expr_f000 , _backend=backend );
    auto expr_f100 = integrate( _range=elements( mesh ), _expr= id( v ) ) ;
    auto functionalfree100 = functionalLinearFree( _space=Xh , _expr=expr_f100 , _backend=backend );
    M_Fqm_free[0][0][0]=functionalfree000;
    M_Fqm_free[1][0][0]=functionalfree100;

    //mass matrix
    auto expr_m00 = integrate ( _range=elements( mesh ), _expr=idt( u )*id( v ) );
    auto operatorfreeM10=opLinearFree( _domainSpace=Xh , _imageSpace=Xh , _expr=expr_m00 , _backend=backend );
    operatorfreeM10->setName("mass");
    M_Mqm_free[0][0]=operatorfreeM10;

    //for scalarProduct
    M = backend->newMatrix( _test=Xh, _trial=Xh );
    form2( Xh, Xh, M ) =
        integrate( _range=elements( mesh ), _expr=gradt( u )*trans( grad( v ) ) ) +
        integrate( _range= markedfaces( mesh, "left" ), _expr= 0.01 * idt( u )*id( v ) ) +
        integrate( _range= markedfaces( mesh, "gamma_holes" ), _expr= 0.001 * idt( u )*id( v ) )
        ;
    M->close();

    //scalar product used for mass matrix
    InnerMassMatrix = backend->newMatrix( _test=Xh, _trial=Xh );
    form2( Xh, Xh, InnerMassMatrix ) =
        integrate( _range=elements( mesh ), _expr=idt( u ) * id( v ) ) ;

    //scalar product used for the POD
    Mpod = backend->newMatrix( _test=Xh, _trial=Xh );
    form2( Xh, Xh, Mpod ) =
        integrate( _range=elements( mesh ), _expr=gradt( u )*trans( grad( v ) ) ) +
        integrate( _range= markedfaces( mesh, "left" ), _expr= 0.01 * idt( u )*id( v ) ) +
        integrate( _range= markedfaces( mesh, "gamma_holes" ), _expr= 0.001 * idt( u )*id( v ) )
        ;
    Mpod->close();

    D = backend->newMatrix( Xh, Xh );
    F = backend->newVector( Xh );


    M_compositeA = opLinearComposite( _domainSpace=Xh , _imageSpace=Xh );
    M_compositeA->addList( M_Aqm_free );
    M_compositeM = opLinearComposite( _domainSpace=Xh , _imageSpace=Xh );
    M_compositeM->addList( M_Mqm_free );
    M_compositeF.resize( this->Nl() );
    for(int output=0; output<this->Nl(); output++)
    {
        M_compositeF[output]=functionalLinearComposite( _space=Xh );
        M_compositeF[output]->addList( M_Fqm_free[output] );
    }



}


typename HeatShield::sparse_matrix_ptrtype
HeatShield::newMatrix() const
{
    return backend->newMatrix( Xh, Xh );
}

typename HeatShield::vector_ptrtype
HeatShield::newVector() const
{
    return backend->newVector( Xh );
}


typename HeatShield::affine_decomposition_type
HeatShield::computeAffineDecomposition()
{
    return boost::make_tuple( M_Mqm, M_Aqm, M_Fqm );
}



void HeatShield::solve( sparse_matrix_ptrtype& D,
                        element_type& u,
                        vector_ptrtype& F )
{

    vector_ptrtype U( backend->newVector( u.functionSpace() ) );
    backend->solve( D, D, U, F );
    u = *U;
}


void HeatShield::exportResults( double time, element_type& T, parameter_type const& mu )
{
    LOG(INFO) << "exportResults starts\n";
    std::string exp_name = "Model_T" + ( boost::format( "_%1%" ) %time ).str();
    export_ptrtype exporter;
    exporter = export_ptrtype( Exporter<mesh_type>::New( "ensight", exp_name  ) );
    exporter->step( time )->setMesh( T.functionSpace()->mesh() );
    std::string mu_str;

    for ( int i=0; i<mu.size(); i++ )
    {
        mu_str= mu_str + ( boost::format( "_%1%" ) %mu[i] ).str() ;
    }

    std::string name = "T_with_parameters_"+mu_str;
    exporter->step( time )->add( name, T );
    exporter->save();
}


void
HeatShield::stockAffineDecomposition()
{
    auto compositeM = operatorCompositeM();
    int q_max = this->Qm();
    M_Mqm.resize( q_max);
    for(int q=0; q<q_max; q++)
    {
        int m_max = this->mMaxM(q);
        M_Mqm[q].resize(m_max);
        for(int m=0; m<m_max;m++)
        {
            auto operatorfree = compositeM->operatorlinear(q,m);
            size_type pattern = operatorfree->pattern();
            auto trial = operatorfree->domainSpace();
            auto test=operatorfree->dualImageSpace();
            M_Mqm[q][m]= backend->newMatrix( _test=test , _trial=trial , _pattern=pattern );
            operatorfree->matPtr(M_Mqm[q][m]);//fill the matrix
        }//m
    }//q

    auto compositeA = operatorCompositeA();
    q_max = this->Qa();
    M_Aqm.resize( q_max);
    for(int q=0; q<q_max; q++)
    {
        int m_max = this->mMaxA(q);
        M_Aqm[q].resize(m_max);
        for(int m=0; m<m_max;m++)
        {
            auto operatorfree = compositeA->operatorlinear(q,m);
            size_type pattern = operatorfree->pattern();
            auto trial = operatorfree->domainSpace();
            auto test=operatorfree->dualImageSpace();
            M_Aqm[q][m]= backend->newMatrix( _test=test , _trial=trial , _pattern=pattern );
            operatorfree->matPtr(M_Aqm[q][m]);//fill the matrix
        }//m
    }//q

    auto vector_compositeF = functionalCompositeF();
    int number_outputs = vector_compositeF.size();
    M_Fqm.resize(number_outputs);
    for(int output=0; output<number_outputs; output++)
    {
        auto composite_f = vector_compositeF[output];
        int q_max = this->Ql(output);
        M_Fqm[output].resize( q_max);
        for(int q=0; q<q_max; q++)
        {
            int m_max = this->mMaxF(output,q);
            M_Fqm[output][q].resize(m_max);
            for(int m=0; m<m_max;m++)
            {
                auto operatorfree = composite_f->functionallinear(q,m);
                auto space = operatorfree->space();
                M_Fqm[output][q][m]= backend->newVector( space );
                operatorfree->containerPtr(M_Fqm[output][q][m]);//fill the vector
            }//m
        }//q
    }//output

}

void HeatShield::update( parameter_type const& mu,double bdf_coeff, element_type const& bdf_poly, int output_index )
{

    if (option(_name="crb.stock-matrices"). as<bool>() )
    {
        D->close();
        D->zero();

        for ( size_type q = 0; q < Qa(); ++q )
        {
            for ( size_type m = 0; m < mMaxA(q); ++m )
                D->addMatrix( M_betaAqm[q][m] , M_Aqm[q][m] );
        }

        F->close();
        F->zero();

        for ( size_type q = 0; q < Ql(output_index); ++q )
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




typename HeatShield::element_type
HeatShield::solve( parameter_type const& mu )
{
    element_ptrtype T( new element_type( Xh ) );
    this->solve( mu, T );
    return *T;
}


void HeatShield::solve( parameter_type const& mu, element_ptrtype& T, int output_index )
{

    using namespace Feel::vf;

    initializationField( T,mu );
    initializationField( pT,mu );

    //assemble();

    if ( M_is_steady )
        M_bdf->setSteady();


    for ( M_bdf->start(*T); !M_bdf->isFinished() ; M_bdf->next() )
    {
        double bdf_coeff = M_bdf->polyDerivCoefficient( 0 );

        this->computeBetaQm( mu, M_bdf->time() );

        auto bdf_poly = M_bdf->polyDeriv();
        this->update( mu , bdf_coeff, bdf_poly );

        backend->solve( _matrix=D,  _solution=T, _rhs=F  );

        if ( do_export )
        {
            exportResults( M_bdf->time(), *T , mu );
            export_number++;
        }

        M_bdf->shiftRight( *T );

    }

}


int
HeatShield::computeNumberOfSnapshots()
{
    return M_bdf->timeFinal()/M_bdf->timeStep();
}


void HeatShield::l2solve( vector_ptrtype& u, vector_ptrtype const& f )
{
    //std::cout << "l2solve(u,f)\n";
    backendl2->solve( _matrix=M,  _solution=u, _rhs=f );
    //std::cout << "l2solve(u,f) done\n";
}


double HeatShield::scalarProduct( vector_ptrtype const& x, vector_ptrtype const& y )
{
    return M->energy( x, y );
}

double HeatShield::scalarProduct( vector_type const& x, vector_type const& y )
{
    return M->energy( x, y );
}

double HeatShield::scalarProductForMassMatrix( vector_ptrtype const& x, vector_ptrtype const& y )
{
    return InnerMassMatrix->energy( x, y );
}

double HeatShield::scalarProductForMassMatrix( vector_type const& x, vector_type const& y )
{
    return InnerMassMatrix->energy( x, y );
}

double HeatShield::scalarProductForPod( vector_ptrtype const& x, vector_ptrtype const& y )
{
    return Mpod->energy( x, y );
}

double HeatShield::scalarProductForPod( vector_type const& x, vector_type const& y )
{
    return Mpod->energy( x, y );
}


void HeatShield::run( const double * X, unsigned long N, double * Y, unsigned long P )
{
    using namespace vf;
    Feel::ParameterSpace<ParameterSpaceDimension>::Element mu( M_Dmu );
    mu << X[0], X[1], X[2];
    static int do_init = true;

    if ( do_init )
    {
        this->initModel();
        do_init = false;
    }

    this->solve( mu, pT );

    double mean = integrate( elements( mesh ),idv( *pT ) ).evaluate()( 0,0 );
    Y[0]=mean;
}



double HeatShield::output( int output_index, parameter_type const& mu, element_type &T, bool need_to_solve , bool export_outputs )
{
    using namespace vf;


    if ( need_to_solve )
        this->solve( mu, pT );
    else
        *pT = T;

    // vector_ptrtype U( backend->newVector( Xh ) );
    //*U = *pT;
    pT->close();
    double s=0;

    auto fqm = backend->newVector( Xh );
    if ( output_index<2 )
    {
        for ( int q=0; q<Ql( output_index ); q++ )
        {
            for ( int m=0; m<mMaxF(output_index,q); m++ )
            {
                M_Fqm_free[output_index][q][m]->containerPtr( fqm );
                s += M_betaFqm[output_index][q][m]*dot( *fqm , *pT );
                // s += M_betaFqm[output_index][q][m]*dot( *M_Fqm[output_index][q][m], *pT );
            }
        }
    }
    else
    {
        throw std::logic_error( "[HeatShield::output] error with output_index : only 0 or 1 " );
    }

    return s ;
}

}

#endif /* __HeatShield_H */


