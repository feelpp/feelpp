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
   \file heatSink.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-11-13
 */
#ifndef __HeatSink2D_H
#define __HeatSink2D_H 1

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

#include <feel/feeldiscr/bdf2.hpp>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>




namespace Feel
{

po::options_description
makeHeatSink2DOptions()
{
    po::options_description heatsink2doptions( "HeatSink2D options" );
    heatsink2doptions.add_options()
    // mesh parameters
    ( "hsize", Feel::po::value<double>()->default_value( 1e-3 ),
      "first h value to start convergence" )
    ( "Lref", Feel::po::value<double>()->default_value( 0.02 ),
      "dimensional length of the sink (in meters)" )
    ( "width", Feel::po::value<double>()->default_value( 0.0005 ),
      "dimensional width of the fin (in meters)" )
    ( "do-export", Feel::po::value<bool>()->default_value( false ),
      "export results if true" )
    ( "steady", Feel::po::value<bool>()->default_value( true ),
      "if true : steady else unsteady" )
    ( "rho", Feel::po::value<double>()->default_value( 8940 ),
      "density in SI unit kg.m^{-3}" )
    ( "C", Feel::po::value<double>()->default_value( 385 ),
      "heat capacity in SI unit J.kg^{-1}.K^{-1}" )
    ( "k_fin", Feel::po::value<double>()->default_value( 386 ),
      "thermal conductivity of the fin in SI unit W.m^{-1}.K^{-1}" )
    ;
    return heatsink2doptions.add( Feel::feel_options() ).add( bdf_options( "heatSink2d" ) );
}
AboutData
makeHeatSink2DAbout( std::string const& str = "heatSink" )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "heat sink Benchmark",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2010,2011 Université de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Stephane Veys", "developer", "stephane.veys@imag.fr", "" );
    return about;
}




/**
 * \class HeatSink2D
 * \brief brief description
 *
 * @author Christophe Prud'homme
 * @see
 */
class HeatSink2D
{
public:


    /** @name Constants
     */
    //@{

    static const uint16_type Order = 1;
    static const uint16_type ParameterSpaceDimension = 3;
    //static const bool is_time_dependent = false;
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


    typedef Eigen::VectorXd theta_vector_type;

    typedef boost::tuple<std::vector<sparse_matrix_ptrtype>, std::vector<sparse_matrix_ptrtype>, std::vector<std::vector<vector_ptrtype>  > > affine_decomposition_type;
    //typedef boost::tuple<std::vector<sparse_matrix_ptrtype>, std::vector<std::vector<vector_ptrtype>  > > affine_decomposition_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    HeatSink2D();

    //! constructor from command line
    HeatSink2D( po::variables_map const& vm );


    //! copy constructor
    HeatSink2D( HeatSink2D const & );
    //! destructor
    ~HeatSink2D() {}

    //! initialization of the model
    void init();
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
        return 4;
    }

    // \return the number of terms in affine decomposition of bilinear form
    // associated to mass matrix
    int Qm() const
    {
        return 2;
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
        return 1*( l==0 ) + 1*( l>0 );
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
    boost::tuple<theta_vector_type, theta_vector_type, std::vector<theta_vector_type> >
    //boost::tuple<theta_vector_type, std::vector<theta_vector_type> >
    computeThetaq( parameter_type const& mu , double time=1e30 )
    {
        double biot      = mu( 0 );
        double L         = mu( 1 );
        double k         = mu( 2 );
        double detJ = L/Lref;
        //coefficient from J^{-T}J^{-1}[4]
        double JJ4 = ( Lref*Lref )/( L*L );
        double t = M_bdf->time();

        M_thetaMq.resize( Qm() );
        M_thetaMq( 0 )=rho*C/k_fin;
        M_thetaMq( 1 )=rho*C/k_fin * detJ;

        M_thetaAq.resize( Qa() );
        M_thetaAq( 0 ) = k ;
        M_thetaAq( 1 ) = 1 ;
        M_thetaAq( 2 ) = JJ4*detJ;
        M_thetaAq( 3 ) = biot * detJ;

        M_thetaFq.resize( Nl() );

        M_thetaFq[0].resize( Ql( 0 ) );
        M_thetaFq[0]( 0 ) = 1-exp( -t );

        M_thetaFq[1].resize( Ql( 1 ) );
        M_thetaFq[1]( 0 ) = 2;


        return boost::make_tuple( M_thetaMq, M_thetaAq, M_thetaFq );
        //return boost::make_tuple( M_thetaAq, M_thetaFq );
    }

    /**
     * \brief return the coefficient vector
     */
    theta_vector_type const& thetaAq() const
    {
        return M_thetaAq;
    }

    /**
     * \brief return the coefficient vector
     */
    theta_vector_type const& thetaMq() const
    {
        return M_thetaMq;
    }


    /**
     * \brief return the coefficient vector
     */
    std::vector<theta_vector_type> const& thetaFq() const
    {
        return M_thetaFq;
    }

    /**
     * \brief return the coefficient vector \p q component
     *
     */
    value_type thetaAq( int q ) const
    {
        return M_thetaAq( q );
    }

    /**
     * \brief return the coefficient vector \p q component
     *
     */
    value_type thetaMq( int q ) const
    {
        return M_thetaMq( q );
    }


    /**
     * \return the \p q -th term of the \p l -th output
     */
    value_type thetaL( int l, int q ) const
    {
        return M_thetaFq[l]( q );
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
     * \brief Returns the affine decomposition
     */
    affine_decomposition_type computeAffineDecomposition();

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
    void solve( parameter_type const& mu );

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

    void exportOutputs( double time, double output1, double output2 );

    void solve( sparse_matrix_ptrtype& ,element_type& ,vector_ptrtype&  );

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
    value_type output( int output_index, parameter_type const& mu, bool export_outputs=false );

    gmsh_ptrtype createGeo( double hsize, double mu2 );

private:

    po::variables_map M_vm;

    backend_ptrtype backend;
    bool M_is_steady ;

    /* mesh parameters */
    double meshSize;

    double Lref;//reference value for geometric parameter

    int export_number;

    bool do_export;

    double rho;
    double C;
    double k_fin;


    parameterspace_ptrtype M_Dmu;

    /* surfaces*/
    double surface_gamma4, surface_gamma1;

    /* mesh, pointers and spaces */
    mesh_ptrtype mesh;
    space_ptrtype Xh;

    sparse_matrix_ptrtype D,M,Mpod;
    vector_ptrtype F;

    element_ptrtype pT;

    std::vector<sparse_matrix_ptrtype> M_Aq;
    std::vector<sparse_matrix_ptrtype> M_Mq;
    std::vector<std::vector<vector_ptrtype> > M_Fq;

    theta_vector_type M_thetaAq;
    theta_vector_type M_thetaMq;
    std::vector<theta_vector_type> M_thetaFq;

    bdf_ptrtype M_bdf;
    int M_Nsnap;

};
HeatSink2D::HeatSink2D()
    :
    backend( backend_type::build( BACKEND_PETSC ) ),
    M_is_steady( false ),
    meshSize( 1e-3 ),
    Lref( 2 ),
    export_number( 0 ),
    do_export( false ),
    rho( 8940 ),
    C( 385 ),
    k_fin( 386 ),
    M_Dmu( new parameterspace_type )
{
    this->init();
}

HeatSink2D::HeatSink2D( po::variables_map const& vm )
    :
    M_vm( vm ),
    backend( backend_type::build( vm ) ),
    M_is_steady( vm["steady"].as<bool>() ),
    meshSize( vm["hsize"].as<double>() ),
    Lref( vm["Lref"].as<double>() ),
    export_number( 0 ),
    do_export( vm["do-export"].as<bool>() ),
    rho( vm["rho"].as<double>() ),
    C( vm["C"].as<double>() ),
    k_fin( vm["k_fin"].as<double>() ),
    M_Dmu( new parameterspace_type )
{
    this->init();
}





gmsh_ptrtype
HeatSink2D::createGeo( double hsize, double mu2 )
{
    int elements_order = 1;
    int dim = 2;
    gmsh_ptrtype gmshp( new Gmsh ( dim , elements_order ) );
    std::ostringstream ostr;


    ostr << "Point (1) = {0   , 0  , 0, " << hsize << "};\n"
         << "Point (2) = {0.5 , 0  , 0, " << hsize << "};\n"
         << "Point (3) = {0.5 , 0.6, 0, " << hsize << "};\n"
         << "Point (4) = {0.15, 0.6, 0, " << hsize << "};\n"
         << "Point (5) = {0.15, 0.6+"<<mu2/2<<", 0, " << hsize << "};\n"
         << "Point (6) = {0.15, 0.6+"<<mu2<<  ", 0, " << hsize << "};\n"
         << "Point (7) = {0, 0.6+"<<mu2<<  ", 0, " << hsize << "};\n"
         << "Point (8) = {0   , 0.6+"<<mu2/2<<", 0, " << hsize << "};\n"
         << "Point (9) = {0   , 0.6, 0, " << hsize << "};\n"
         << "Line (1)  = {1, 2};\n"
         << "Line (2)  = {2, 3};\n"
         << "Line (3)  = {3, 4};\n"
         << "Line (4)  = {4, 9};\n"
         << "Line (5)  = {4, 5};\n"
         << "Line (6)  = {5, 6};\n"
         << "Line (7)  = {6, 7};\n"
         << "Line (8)  = {7, 8};\n"
         << "Line (9)  = {8, 9};\n"
         << "Line (10) = {9, 1};\n"
         << "Line Loop (20) = {1, 2, 3, 4, 10};\n"
         << "Line Loop (21) = {5, 6, 7, 8, 9, -4};\n"
         << "Plane Surface (20) = {20};\n"
         << "Plane Surface (21) = {21};\n"
         << "Physical Line (\"gamma1\")  = {1};\n"
         << "Physical Line (\"gamma2\")  = {2};\n"
         << "Physical Line (\"gamma3\")  = {3};\n"
         << "Physical Line (\"gamma4\")  = {4};\n"
         << "Physical Line (\"gamma5\")  = {5};\n"
         << "Physical Line (\"gamma6\")  = {6};\n"
         << "Physical Line (\"gamma7\")  = {7};\n"
         << "Physical Line (\"gamma8\")  = {8};\n"
         << "Physical Line (\"gamma9\")  = {9};\n"
         << "Physical Line (\"gamma10\") = {10};\n"
         << "Physical Surface (\"spreader_mesh\") = {20};\n"
         << "Physical Surface (\"fin_mesh\") = {21};\n";

    std::ostringstream nameStr;
    nameStr.precision( 3 );
    nameStr << "fin_sink";
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( ostr.str() );
    return gmshp;
}




void HeatSink2D::initializationField( element_ptrtype& initial_field , parameter_type const& mu )
{

    initial_field->setZero();
    //in the case where initial condition is a parameter we can do :
    //initial_field->setOnes();
    //double Tambiant   = mu( 7 );
    //initial_field->scale(Tambiant);
}

void HeatSink2D::init()
{


    using namespace Feel::vf;

    /*
     * First we create the mesh
     */

    mesh = createGMSHMesh ( _mesh = new mesh_type,
                            _desc = createGeo( meshSize, Lref ),
                            _update=MESH_UPDATE_FACES | MESH_UPDATE_EDGES );

    /*
     * The function space and some associate elements are then defined
     */
    Xh = space_type::New( mesh );
    std::cout << "Number of dof " << Xh->nLocalDof() << "\n";

    // allocate an element of Xh
    pT = element_ptrtype( new element_type( Xh ) );


    std::cout<<"evaluation of surface_gamma1 ..."<<std::endl;
    surface_gamma1 = integrate( _range= markedfaces( mesh,"gamma1" ), _expr=cst( 1. ) ).evaluate()( 0,0 );
    std::cout<<"surface_gamma = "<<surface_gamma1<<std::endl;

    M_bdf = bdf( _space=Xh, _vm=M_vm, _name="heatSink2d" , _prefix="heatSink2d" );

    M_Aq.resize( this->Qa() );
    M_Mq.resize( this->Qm() );

    //three outputs : one "compliant" and two others
    M_Fq.resize( this->Nl() );
    //first : the compliant case.
    M_Fq[0].resize( 1 );
    M_Fq[0][0] = backend->newVector( Xh );
    //second output : non compliant case.
    M_Fq[1].resize( 1 );
    M_Fq[1][0] = backend->newVector( Xh );

    D = backend->newMatrix( Xh, Xh );
    F = backend->newVector( Xh );

    Feel::ParameterSpace<ParameterSpaceDimension>::Element mu_min( M_Dmu );
    //mu_min <<  /* Bi */ 0.4 , /*L*/2, /*k*/1;
    mu_min <<  /* Bi */ 0.1 , /*L*/2, /*k*/1;
    M_Dmu->setMin( mu_min );
    Feel::ParameterSpace<ParameterSpaceDimension>::Element mu_max( M_Dmu );
    mu_max << /* Bi */ 0.5  ,  /*L*/8 , /*k*/10;
    //mu_max <<  /* Bi */ 0.1 , /*L*/2, /*k*/1;
    M_Dmu->setMax( mu_max );

    LOG(INFO) << "Number of dof " << Xh->nLocalDof() << "\n";

    assemble();


} // HeatSink2D::init



void HeatSink2D::assemble()
{

    M_bdf->start();
    using namespace Feel::vf;

    element_type u( Xh, "u" );
    element_type v( Xh, "v" );


    M_Aq[0] = backend->newMatrix( _test=Xh, _trial=Xh );
    form2( _test=Xh, _trial=Xh, _matrix=M_Aq[0] ) = integrate( _range= markedfaces( mesh, "spreader_mesh" ), _expr= gradt( u )*trans( grad( v ) ) );

    M_Aq[1] = backend->newMatrix( _test=Xh, _trial=Xh , _matrix=M_Aq[0] );
    M_Aq[2] = backend->newMatrix( _test=Xh, _trial=Xh , _matrix=M_Aq[0] );
    M_Aq[3] = backend->newMatrix( _test=Xh, _trial=Xh , _matrix=M_Aq[0] );

    form2( _test=Xh, _trial=Xh, _matrix=M_Aq[2] ) = integrate( _range= markedelements( mesh,"fin_mesh" ),_expr=  dy( v )*dyt( u )  );
    form2( _test=Xh, _trial=Xh, _matrix=M_Aq[1] ) = integrate( _range= markedelements( mesh,"fin_mesh" ),_expr=  dx( v )*dxt( u )  );
    form2( _test=Xh, _trial=Xh, _matrix=M_Aq[3] )  = integrate( _range= markedfaces( mesh, "gamma5" ), _expr= idt( u )*id( v ) );
    form2( _test=Xh, _trial=Xh, _matrix=M_Aq[3] ) += integrate( _range= markedfaces( mesh, "gamma6" ), _expr= idt( u )*id( v ) );


    form1( _test=Xh, _vector=M_Fq[0][0] ) = integrate( _range=markedfaces( mesh,"gamma1" ), _expr= id( v ) ) ;
    form1( _test=Xh, _vector=M_Fq[1][0] ) = integrate( _range=markedfaces( mesh,"gamma1" ), _expr= id( v ) ) ;

    M_Aq[0]->close();
    M_Aq[1]->close();
    M_Aq[2]->close();
    M_Aq[3]->close();

    M_Fq[0][0]->close();
    M_Fq[1][0]->close();

    //mass matrix
    M_Mq[0] = backend->newMatrix( _test=Xh, _trial=Xh , _matrix=M_Aq[0] );
    M_Mq[1] = backend->newMatrix( _test=Xh, _trial=Xh , _matrix=M_Aq[0] );
    form2( _test=Xh, _trial=Xh, _matrix=M_Mq[0] ) = integrate ( _range=markedelements( mesh,"spreader_mesh" ), _expr=idt( u )*id( v ) );
    form2( _test=Xh, _trial=Xh, _matrix=M_Mq[1] ) = integrate ( _range=markedelements( mesh,"fin_mesh" ),      _expr=idt( u )*id( v ) );
    M_Mq[0]->close();
    M_Mq[1]->close();

    //for scalarProduct
    M = backend->newMatrix( _test=Xh, _trial=Xh );
    form2( Xh, Xh, M ) =
        integrate( elements( mesh ), id( u )*idt( v ) + grad( u )*trans( gradt( u ) ) );
    M->close();

    Mpod = backend->newMatrix( _test=Xh, _trial=Xh );
    form2( Xh, Xh, Mpod ) =
        integrate( elements( mesh ), id( u )*idt( v ) );
    Mpod->close();


}


typename HeatSink2D::sparse_matrix_ptrtype
HeatSink2D::newMatrix() const
{
    return backend->newMatrix( Xh, Xh , M_Aq[0] );
}


typename HeatSink2D::affine_decomposition_type
HeatSink2D::computeAffineDecomposition()
{
    return boost::make_tuple( M_Mq, M_Aq, M_Fq );
    //return boost::make_tuple( M_Aq, M_Fq );
}


void HeatSink2D::solve( sparse_matrix_ptrtype& D,
                        element_type& u,
                        vector_ptrtype& F )
{

    vector_ptrtype U( backend->newVector( u.functionSpace() ) );
    backend->solve( D, D, U, F );
    u = *U;
}


void HeatSink2D::exportResults( double time, element_type& T, parameter_type const& mu )
{
    std::cout<<"STRT"<<std::endl;
    LOG(INFO) << "exportResults starts\n";
    std::string exp_name = "Model_T" + ( boost::format( "_%1%" ) %time ).str();
    //export_ptrtype exp = export_ptrtype( Exporter<mesh_type>::New( exp_name ) );
    export_ptrtype exporter;
    //exporter = export_ptrtype( Exporter<mesh_type>::New(M_vm, exp_name  ) );
    exporter = export_ptrtype( Exporter<mesh_type>::New( "ensight", exp_name  ) );
    //exporter exp = export_ptrtype( Exporter<mesh_type>::New( exp_name ) );
    exporter->step( time )->setMesh( T.functionSpace()->mesh() );
    std::string mu_str;

    for ( int i=0; i<mu.size(); i++ )
    {
        mu_str= mu_str + ( boost::format( "_%1%" ) %mu[i] ).str() ;
    }

    std::string name = "T_with_parameters_"+mu_str;
    exporter->step( time )->add( name, T );
    exporter->save();
    std::cout<<" ====================== export ok"<<std::endl;
}

void HeatSink2D::exportOutputs( double time, double output1, double output2 )
{
    std::ofstream output_file;
    output_file.open( "Output.dat",std::ios::out | std::ios::app );
    output_file<<std::setprecision( 16 )<<time<<" "<<output1<<" "<<output2<<"\n";

}


void HeatSink2D::update( parameter_type const& mu,double bdf_coeff, element_type const& bdf_poly, int output_index )
{

    D->close();
    D->zero();

    *D = *M_Aq[0];
    D->scale( M_thetaAq[0] );

    for ( size_type q = 1; q < M_Aq.size(); ++q )
    {
        D->addMatrix( M_thetaAq[q] , M_Aq[q] );
    }

    F->close();
    F->zero();

    for ( size_type q = 0; q < M_Fq[output_index].size(); ++q )
    {
        F->add( M_thetaFq[output_index][q], M_Fq[output_index][q] );
    }

    auto vec_bdf_poly = backend->newVector( Xh );

    //add contribution from mass matrix
    for ( size_type q = 0; q < M_Mq.size(); ++q )
    {
        //left hand side
        D->addMatrix( M_thetaMq[q]*bdf_coeff, M_Mq[q] );
        //right hand side
        *vec_bdf_poly = bdf_poly;
        vec_bdf_poly->scale( M_thetaMq[q] );
        F->addVector( *vec_bdf_poly, *M_Mq[q] );
    }

}




void HeatSink2D::solve( parameter_type const& mu )
{
    element_ptrtype T( new element_type( Xh ) );
    this->solve( mu, T );
}


void HeatSink2D::solve( parameter_type const& mu, element_ptrtype& T, int output_index )
{



    using namespace Feel::vf;

    initializationField( T,mu );
    initializationField( pT,mu );

    //*T = vf::project( _space=Xh, _expr=cst(Tamb) );
    //*pT = vf::project( _space=Xh, _expr=cst(Tamb) );

    // pT->zero();
    // T->zero();

    assemble();

    element_type v( Xh, "v" );//test functions

    M_bdf->initialize( *T );

    if ( M_is_steady )
    {
        M_bdf->setSteady();
    }

    double bdf_coeff = M_bdf->polyDerivCoefficient( 0 );


    for ( M_bdf->start(); !M_bdf->isFinished() ; M_bdf->next() )
    {

        this->computeThetaq( mu, M_bdf->time() );
        auto bdf_poly = M_bdf->polyDeriv();
        this->update( mu , bdf_coeff, bdf_poly );

        backend->solve( _matrix=D,  _solution=T, _rhs=F );

        do_export=true;

        if ( do_export )
        {
            exportResults( M_bdf->time(), *T , mu );
            export_number++;
        }



#if 0
        std::ofstream file;
        std::string mu_str;

        for ( int i=0; i<mu.size(); i++ )
        {
            mu_str= mu_str + ( boost::format( "_%1%" ) %mu[i] ).str() ;
        }

        std::string number =  ( boost::format( "Exp_%1%" ) %export_number ).str();
        std::string name = "PFEMsolution" + mu_str + number;
        file.open( name,std::ios::out );

        for ( int i=0; i<T->size(); i++ ) file<<T->operator()( i )<<"\n";

        file.close();

        name = "PFEMrhs" + mu_str + number;
        std::ofstream file_rhs;
        file_rhs.open( name,std::ios::out );
        file_rhs<<*F;
        file_rhs.close();
#endif

#if 0
        std::ofstream file_matrix;
        name = "PFEMmatrix" + mu_str + number;
        file_matrix.open( name,std::ios::out );
        file_matrix<<*D;
        file_matrix.close();
#endif

        M_bdf->shiftRight( *T );

    }


}


int
HeatSink2D::computeNumberOfSnapshots()
{
    M_Nsnap = M_bdf->timeFinal()/M_bdf->timeStep();
    return  M_Nsnap;
}


void HeatSink2D::l2solve( vector_ptrtype& u, vector_ptrtype const& f )
{
    //std::cout << "l2solve(u,f)\n";
    backend->solve( _matrix=M,  _solution=u, _rhs=f );
    //std::cout << "l2solve(u,f) done\n";
}


double HeatSink2D::scalarProduct( vector_ptrtype const& x, vector_ptrtype const& y )
{
    return M->energy( x, y );
}

double HeatSink2D::scalarProduct( vector_type const& x, vector_type const& y )
{
    return M->energy( x, y );
}

double HeatSink2D::scalarProductForPod( vector_ptrtype const& x, vector_ptrtype const& y )
{
    return Mpod->energy( x, y );
}

double HeatSink2D::scalarProductForPod( vector_type const& x, vector_type const& y )
{
    return Mpod->energy( x, y );
}


void HeatSink2D::run( const double * X, unsigned long N, double * Y, unsigned long P )
{
    using namespace vf;
    Feel::ParameterSpace<ParameterSpaceDimension>::Element mu( M_Dmu );
    mu << X[0], X[1], X[2];
    static int do_init = true;

    if ( do_init )
    {
        this->init();
        do_init = false;
    }

    this->solve( mu, pT );

    double mean = integrate( elements( mesh ),idv( *pT ) ).evaluate()( 0,0 );
    Y[0]=mean;
}



double HeatSink2D::output( int output_index, parameter_type const& mu, bool export_outputs )
{
    using namespace vf;


    element_type u( Xh, "u" );
    element_type v( Xh, "v" );

    if ( !export_outputs )
    {
        this->solve( mu, pT );
    }

    vector_ptrtype U( backend->newVector( Xh ) );
    *U = *pT;

    //std::cout<<"*U : "<<*U<<std::endl;

    double s=0;

    if ( output_index==0 )
    {
        for ( int i=0; i<Ql( output_index ); i++ )  s += M_thetaFq[output_index]( i )*dot( M_Fq[output_index][i], U );
    }

    else if ( output_index==1 )
    {
        s = M_thetaFq[output_index]( 0 )*dot( M_Fq[output_index][0], U );
        std::cout<<" s model = "<<s<<std::endl;
    }

    else
    {
        throw std::logic_error( "[HeatSink2D::output] error with output_index : only 0 or 1 " );
    }

    return s ;
}

}

#endif /* __HeatSink2D_H */


