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

#include <feel/feelts/bdf.hpp>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>


#include <feel/feelcrb/modelcrbbase.hpp>
#include <feel/feeldiscr/reducedbasisspace.hpp>

namespace Feel
{

po::options_description
makeHeatSink2DOptions()
{
    po::options_description heatsink2doptions( "HeatSink2D options" );
    heatsink2doptions.add_options()
    // mesh parameters
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
    return heatsink2doptions.add( bdf_options( "heatSink2d" ) );
}
AboutData
makeHeatSink2DAbout( std::string const& str = "heatSink" )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "heat sink Benchmark",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2010,2011 Universite de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Stephane Veys", "developer", "stephane.veys@imag.fr", "" );
    return about;
}


class FunctionSpaceDefinition
{
public :
    typedef double value_type;

    /*mesh*/
    typedef Simplex<2,1> entity_type;
    typedef Mesh<entity_type> mesh_type;

    /*basis*/
    typedef fusion::vector<Lagrange<1, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;

    /*element*/
    typedef typename space_type::element_type element_type;

    static const bool is_time_dependent = true;
    static const bool is_linear = true;

};

/**
 * \class HeatSink2D
 * \brief brief description
 *
 * @author Christophe Prud'homme
 * @see
 */
class HeatSink2D : public ModelCrbBase< ParameterSpace<3>, FunctionSpaceDefinition, 1 >
{
public:

    typedef ModelCrbBase<ParameterSpace<3>,FunctionSpaceDefinition, 1> super_type;
    typedef typename super_type::funs_type funs_type;
    typedef typename super_type::funsd_type funsd_type;


    /** @name Constants
     */
    //@{

    static const uint16_type Order = 1;
    static const uint16_type ParameterSpaceDimension = 3;
    //static const bool is_time_dependent = false;

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

    /* time discretization */
    typedef Bdf<space_type>  bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;


    //typedef std::vector< std::vector< double > > beta_vector_type;

    /*typedef boost::tuple<
        std::vector< std::vector<sparse_matrix_ptrtype> >,
        std::vector< std::vector<sparse_matrix_ptrtype> >,
        std::vector< std::vector<std::vector<vector_ptrtype> > >
     > affine_decomposition_type;*/

    //typedef boost::tuple<std::vector<sparse_matrix_ptrtype>, std::vector<std::vector<vector_ptrtype>  > > affine_decomposition_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    HeatSink2D();

    //! copy constructor
    HeatSink2D( HeatSink2D const & );
    //! destructor
    ~HeatSink2D() {}

    //! initialization of the model
    void initModel();
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

    int mMaxA( int q )
    {
        if( q < this->Qa() )
            return 1;
        else
            throw std::logic_error( "[Model heatsink] ERROR : try to acces to mMaxA(q) with a bad value of q");
    }

    int mMaxM( int q )
    {
        if ( q < this->Qm() )
            return 1;
        else
            throw std::logic_error( "[Model heatsink] ERROR : try to acces to mMaxM(q) with a bad value of q");
    }

    int mMaxF( int output_index, int q)
    {
        if ( q < this->Ql(output_index) )
            return 1;
        else
            throw std::logic_error( "[Model heatsink] ERROR : try to acces to mMaxF(output_index,q) with a bad value of q");
    }

    /**
     * \brief compute the beta coefficient for both bilinear and linear form
     * \param mu parameter to evaluate the coefficients
     */

    betaqm_type
    computeBetaQm( element_type const& T, parameter_type const& mu , double time , bool only_terms_time_dependent=false )
    {
        return computeBetaQm( mu, time, only_terms_time_dependent );
    }

    betaqm_type
    computeBetaQm( parameter_type const& mu , double time , bool only_terms_time_dependent=false )
    {
        double biot      = mu( 0 );
        double L         = mu( 1 );
        double k         = mu( 2 );
        double detJ = L/Lref;
        //coefficient from J^{-T}J^{-1}[4]
        double JJ4 = ( Lref*Lref )/( L*L );
        double t = M_bdf->time();

        M_betaMqm.resize( Qm() );
        M_betaMqm[0].resize(1);
        M_betaMqm[1].resize(1);
        M_betaMqm[0][0]=rho*C/k_fin;
        M_betaMqm[1][0]=rho*C/k_fin * detJ;

        M_betaAqm.resize( Qa() );
        for(int i=0; i<Qa(); i++)
            M_betaAqm[i].resize(1);
        M_betaAqm[0][0] = k ;
        M_betaAqm[1][0] = 1 ;
        M_betaAqm[2][0] = JJ4*detJ;
        M_betaAqm[3][0] = biot * detJ;

        M_betaFqm.resize( Nl() );
        for(int i=0; i<Nl(); i++)
        {
            M_betaFqm[i].resize( Ql(i) );
            for(int j=0;j<Ql(i);j++)
                M_betaFqm[i][j].resize(1);
        }

        M_betaFqm[0][0][0] = 1-math::exp( -t );
        M_betaFqm[1][0][0] = 2;

        return boost::make_tuple( M_betaMqm, M_betaAqm, M_betaFqm );
    }

    //@}

    /** @name  Methods
     */
    //@{

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


    int computeNumberOfSnapshots();

    void assemble();
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
     * H1 scalar product
     */
    sparse_matrix_ptrtype energyMatrix ( void )
    {
        return M;
    }


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
    value_type output( int output_index, parameter_type const& mu, element_type& u, bool need_to_solve=false, bool export_outputs=false );

    gmsh_ptrtype createGeo( double mu2 );

    parameter_type refParameter()
    {
        return Dmu->min();
    }

    bdf_ptrtype bdfModel(){ return M_bdf; }


private:
    bool M_is_steady ;

    double Lref;//reference value for geometric parameter

    int export_number;

    bool do_export;

    double rho;
    double C;
    double k_fin;

    /* surfaces*/
    double surface_gamma4, surface_gamma1;

    /* mesh, pointers and spaces */
    mesh_ptrtype mesh;

    sparse_matrix_ptrtype D,M,Mpod;
    vector_ptrtype F;

    element_ptrtype pT;

    bdf_ptrtype M_bdf;
    int M_Nsnap;

};
HeatSink2D::HeatSink2D()
    :
    M_is_steady( boption("steady") ),
    Lref( doption("Lref") ),
    export_number( 0 ),
    do_export( boption("do-export") ),
    rho( doption("rho") ),
    C( doption("C") ),
    k_fin( doption("k_fin") )
 {}


gmsh_ptrtype
HeatSink2D::createGeo(  double mu2 )
{
    double hsize = doption( "gmsh.hsize");
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

void HeatSink2D::initModel()
{

    using namespace Feel::vf;

    /*
     * First we create the mesh
     */

    mesh = createGMSHMesh ( _mesh = new mesh_type,
                            _desc = createGeo( Lref ),
                            _update=MESH_UPDATE_FACES | MESH_UPDATE_EDGES );

    /*
     * The function space and some associate elements are then defined
     */
    Xh = space_type::New( mesh );
    this->setFunctionSpaces( Xh );
    std::cout << "Number of dof " << Xh->nLocalDof() << "\n";

    // allocate an element of Xh
    pT = element_ptrtype( new element_type( Xh ) );

    surface_gamma1 = integrate( _range= markedfaces( mesh,"gamma1" ), _expr=cst( 1. ) ).evaluate()( 0,0 );

    M_bdf = bdf( _space=Xh, _name="heatSink2d" , _prefix="heatSink2d" );

    M_Aqm.resize( this->Qa() );
    for(int q=0; q<Qa(); q++)
    {
        M_Aqm[q].resize( 1 );
        M_Aqm[q][0] = backend()->newMatrix(Xh,Xh);
    }
    M_Mqm.resize( this->Qm() );
    for(int q=0; q<Qm(); q++)
        M_Mqm[q].resize( 1 );

    //three outputs : one "compliant" and two others
    M_Fqm.resize( this->Nl() );
    //first : the compliant case.
    M_Fqm[0].resize( 1 );
    M_Fqm[0][0].resize(1);
    M_Fqm[0][0][0] = backend()->newVector( Xh );
    //second output : non compliant case.
    M_Fqm[1].resize( 1 );
    M_Fqm[1][0].resize(1);
    M_Fqm[1][0][0] = backend()->newVector( Xh );

    D = backend()->newMatrix( Xh, Xh );
    F = backend()->newVector( Xh );

    auto mu_min = Dmu->element();
    mu_min <<  /* Bi */ 0.1 , /*L*/2, /*k*/1;
    Dmu->setMin( mu_min );
    auto mu_max = Dmu->element();
    mu_max << /* Bi */ 0.5  ,  /*L*/8 , /*k*/10;
    Dmu->setMax( mu_max );

    LOG(INFO) << "Number of dof " << Xh->nLocalDof() << "\n";

    assemble();


} // HeatSink2D::init



void HeatSink2D::assemble()
{

    M_bdf->start();
    using namespace Feel::vf;

    element_type u( Xh, "u" );
    element_type v( Xh, "v" );

    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[0][0] )
        = integrate( _range= markedelements( mesh, "spreader_mesh" ),
                     _expr= gradt( u )*trans( grad( v ) ) );
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[1][0] )
        = integrate( _range= markedelements( mesh,"fin_mesh" ),
                     _expr=  dx( v )*dxt( u )  );
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[2][0] )
        = integrate( _range= markedelements( mesh,"fin_mesh" ),
                     _expr=  dy( v )*dyt( u )  );
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[3][0] )
        = integrate( _range= markedfaces( mesh, "gamma5" ),
                     _expr= idt( u )*id( v ) );
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[3][0] )
        += integrate( _range= markedfaces( mesh, "gamma6" ),
                      _expr= idt( u )*id( v ) );
    M_Aqm[0][0]->close();
    M_Aqm[1][0]->close();
    M_Aqm[2][0]->close();
    M_Aqm[3][0]->close();

    form1( _test=Xh, _vector=M_Fqm[0][0][0] ) =
        integrate( _range=markedfaces( mesh,"gamma1" ), _expr= id( v ) ) ;
    form1( _test=Xh, _vector=M_Fqm[1][0][0] ) =
        integrate( _range=markedfaces( mesh,"gamma1" ), _expr= id( v ) ) ;
    M_Fqm[0][0][0]->close();
    M_Fqm[1][0][0]->close();

    //mass matrix
    M_Mqm[0][0] = backend()->newMatrix( _test=Xh, _trial=Xh );
    M_Mqm[1][0] = backend()->newMatrix( _test=Xh, _trial=Xh );
    form2( _test=Xh, _trial=Xh, _matrix=M_Mqm[0][0] ) =
        integrate ( _range=markedelements( mesh,"spreader_mesh" ),
                    _expr=idt( u )*id( v ) );
    form2( _test=Xh, _trial=Xh, _matrix=M_Mqm[1][0] ) =
        integrate ( _range=markedelements( mesh,"fin_mesh" ),
                    _expr=idt( u )*id( v ) );
    M_Mqm[0][0]->close();
    M_Mqm[1][0]->close();

    //for scalarProduct
    M = backend()->newMatrix( _test=Xh, _trial=Xh );
    form2( _test=Xh, _trial=Xh, _matrix=M ) =
        integrate( elements( mesh ), idt( u )*id( v ) + gradt( u )*trans( grad( v ) ) );
    M->close();
    addEnergyMatrix( M );

    Mpod = backend()->newMatrix( _test=Xh, _trial=Xh );
    form2( _test=Xh, _trial=Xh, _matrix=Mpod ) =
        integrate( elements( mesh ), idt( u )*id( v ) );
    Mpod->close();
    addMassMatrix( Mpod );
}


typename HeatSink2D::sparse_matrix_ptrtype
HeatSink2D::newMatrix() const
{
    return backend()->newMatrix( Xh, Xh );
}

typename HeatSink2D::vector_ptrtype
HeatSink2D::newVector() const
{
    return backend()->newVector( Xh );
}


typename HeatSink2D::affine_decomposition_type
HeatSink2D::computeAffineDecomposition()
{
    return boost::make_tuple( M_Mqm, M_Aqm, M_Fqm );
}


void HeatSink2D::solve( sparse_matrix_ptrtype& D,
                        element_type& u,
                        vector_ptrtype& F )
{

    vector_ptrtype U( backend()->newVector( u.functionSpace() ) );
    backend()->solve( D, D, U, F );
    u = *U;
}


void HeatSink2D::exportResults( double time, element_type& T, parameter_type const& mu )
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

    // *D = *M_Aqm[0][0];
    //D->scale( M_betaAq[0] );

    for ( size_type q = 0; q < M_Aqm.size(); ++q )
    {
        for ( size_type m = 0; m < mMaxA(q); ++m )
        {
            D->addMatrix( M_betaAqm[q][m] , M_Aqm[q][m] );
        }
    }

    F->close();
    F->zero();

    for ( size_type q = 0; q < M_Fqm[output_index].size(); ++q )
    {
        for ( size_type m = 0; m < mMaxF(output_index,q); ++m )
        {
            F->add( M_betaFqm[output_index][q][m], M_Fqm[output_index][q][m] );
        }
    }

    auto vec_bdf_poly = backend()->newVector( Xh );

    //add contribution from mass matrix
    for ( size_type q = 0; q < M_Mqm.size(); ++q )
    {
        for ( size_type m = 0; m < mMaxM(q); ++m )
        {
            //left hand side
            D->addMatrix( M_betaMqm[q][m]*bdf_coeff, M_Mqm[q][m] );
            //right hand side
            *vec_bdf_poly = bdf_poly;
            vec_bdf_poly->scale( M_betaMqm[q][m] );
            F->addVector( *vec_bdf_poly, *M_Mqm[q][m] );
        }
    }

}



typename HeatSink2D::element_type
HeatSink2D::solve( parameter_type const& mu )
{
    element_ptrtype T( new element_type( Xh ) );
    this->solve( mu, T );
    return *T;
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

        this->computeBetaQm( mu, M_bdf->time() );
        auto bdf_poly = M_bdf->polyDeriv();
        this->update( mu , bdf_coeff, bdf_poly );

        backend()->solve( _matrix=D,  _solution=T, _rhs=F );

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
    return M_bdf->timeFinal()/M_bdf->timeStep();
}

void HeatSink2D::l2solve( vector_ptrtype& u, vector_ptrtype const& f )
{
    backend()->solve( _matrix=M,  _solution=u, _rhs=f );
}

void HeatSink2D::run( const double * X, unsigned long N, double * Y, unsigned long P )
{
    using namespace vf;
    auto mu = Dmu->element();
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



double HeatSink2D::output( int output_index, parameter_type const& mu, element_type& unknown, bool need_to_solve, bool export_outputs )
{
    using namespace vf;


    element_type u( Xh, "u" );
    element_type v( Xh, "v" );

    if ( need_to_solve )
        this->solve( mu, pT );
    else
        *pT=unknown;

    vector_ptrtype U( backend()->newVector( Xh ) );
    *U = *pT;

    //std::cout<<"*U : "<<*U<<std::endl;

    double s=0;

    if( output_index < this->Nl() )
    {
        for ( int q=0; q<Ql( output_index ); q++ )
        {
            for ( int m=0; m<mMaxF(output_index,q); m++ )
            {
                s += M_betaFqm[output_index][q][m]*dot( M_Fqm[output_index][q][m], U );
            }
        }
    }
    else
        throw std::logic_error( "[HeatSink2d::output] error with output_index : only 0 or 1 " );


    return s ;
}

}

#endif /* __HeatSink2D_H */
