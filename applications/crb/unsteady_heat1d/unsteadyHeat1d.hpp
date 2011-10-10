/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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

#include <feel/feeldiscr/bdf2.hpp>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>



namespace Feel
{

po::options_description
makeUnsteadyHeat1DOptions()
{
    po::options_description unsteadyheat1doptions("UnsteadyHeat1D options");
    unsteadyheat1doptions.add_options()
        ("hsize", po::value<double>()->default_value( 0.01 ), "mesh size")
        ("mu1", po::value<double>()->default_value( 0.2 ), "mu1")
        ("mu2", po::value<double>()->default_value( 0.2 ), "mu2")
        ("mu3", po::value<double>()->default_value( -1 ), "mu3")
        ("mu4", po::value<double>()->default_value( 0.1 ), "mu4")
        ("alpha" , po::value<double>()->default_value( 1 ), "temporal coefficient")
        ("steady", po::value<bool>()->default_value(false), "if true then steady else unsteady")
        ("no-export", "don't export results")
        ;
    return unsteadyheat1doptions.add( Feel::feel_options() ).add(bdf_options("unsteadyHeat1d"));
}
AboutData
makeUnsteadyHeat1DAbout( std::string const& str = "unsteadyHeat1d" )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "unsteady 1D Heat Benchmark",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2010,2011 Université de Grenoble 1 (Joseph Fourier)");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    about.addAuthor("Stephane Veys", "developer", "stephane.veys@imag.fr", "");
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

/**
 * \class UnsteadyHeat1D
 * \brief brief description
 *
 * @author Christophe Prud'homme
 * @see
 */
class UnsteadyHeat1D
{
public:


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

    typedef FunctionSpace<mesh_type, fusion::vector<Lagrange<0, Scalar> >, Discontinuous> p0_space_type;
    typedef p0_space_type::element_type p0_element_type;

    /*basis*/
    typedef fusion::vector<Lagrange<Order, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef space_type functionspace_type;
    typedef space_ptrtype functionspace_ptrtype;
    typedef space_type::element_type element_type;
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
    ~UnsteadyHeat1D(){}

    //! initialisation of the model
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
    int Qa() const { return 3; }

    // \return the number of terms in affine decomposition of bilinear form
    // associated to mass matrix
    int Qm() const { return 1; }

    // \return the number of terms in affine decomposition of left hand
    // side bilinear form ( adjoint model )
    int Qadu() const { return 3; }

    // \return the number of terms in affine decomposition of right hand
    //  linear form ( adjoint model )
    int Qfdu() const { return 1; }


    /**
     * there is at least one output which is the right hand side of the
     * primal problem
     *
     * \return number of outputs associated to the model
     */
    int Nl() const { return 2; }

    /**
     * \param l the index of output
     * \return number of terms  in affine decomposition of the \p q th output term
     */
    int Ql( int l ) const { if ( l == 0 ) return 2; return 1; }

    /**
     * \brief Returns the function space
     */
    space_ptrtype functionSpace() { return Xh; }

    //! return the parameter space
    parameterspace_ptrtype parameterSpace() const { return M_Dmu;}

    /**
     * \brief compute the theta coefficient for both bilinear and linear form
     * \param mu parameter to evaluate the coefficients
     */
    boost::tuple<theta_vector_type, theta_vector_type, std::vector<theta_vector_type> >
    computeThetaq( parameter_type const& mu )
        {
            M_thetaAq.resize( Qa() );
            M_thetaAq( 0 ) = 1;
            M_thetaAq( 1 ) = mu( 0 ); // k_1
            M_thetaAq( 2 ) = mu( 1 ); // k_2

            M_thetaFq.resize( Nl() );
            M_thetaFq[0].resize( Ql(0) );
            M_thetaFq[0]( 0 ) = mu(2); // delta
            M_thetaFq[0]( 1 ) = mu(3); // phi

            M_thetaFq[1].resize( Ql(1) );
            M_thetaFq[1]( 0 ) = 1;

            M_thetaMq.resize( Qm() );
            M_thetaMq( 0 ) = 1;


            M_thetaAq_du.resize( Qadu() );
            M_thetaAq_du( 0 ) = 1;
            M_thetaAq_du( 1 ) = mu( 0 ); // k_1
            M_thetaAq_du( 2 ) = mu( 1 ); // k_2

            M_thetaFq_du.resize( Qfdu() );
            M_thetaFq_du( 0 ) = 1;

            return boost::make_tuple( M_thetaMq, M_thetaAq, M_thetaFq );
        }

    /**
     * \brief return the coefficient vector
     */
    theta_vector_type const& thetaAq() const { return M_thetaAq; }

    /**
     * \brief return the coefficient vector
     */
    theta_vector_type const& thetaMq() const { return M_thetaMq; }


    /**
     * \brief return the coefficient vector
     */
    std::vector<theta_vector_type> const& thetaFq() const { return M_thetaFq; }

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
    void setMeshSize( double s ) { meshSize = s; }


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
     * \param output_index (optional) useful for the adjoint ( POD )
     */
    void solve( parameter_type const& mu, element_ptrtype& T , int output_index=0);

    eigen_matrix_type snapshotsMatrix(){return M_snapshots_matrix;}
    eigen_matrix_type dualSnapshotsMatrix(){return Mdu_snapshots_matrix;}
    void fillSnapshotsMatrix (parameter_type const& mu, int output);
    void assemble();
    int computeNumberOfSnapshots();
    double timeStep() { return  M_bdf->timeStep(); }
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
    void update( parameter_type const& mu, double bdf_coeff, element_type const& bdf_poly, int output_index=0 );
    //@}

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults(double time, element_type& T );


    void exportResults1d( double time, element_type& T, double s);

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
    value_type output( int output_index, parameter_type const& mu , bool export_output=false);

    static const bool is_time_dependent = true;
private:

    bool steady ;
    bool pfem;
    double alpha;

    po::variables_map M_vm;
    backend_ptrtype backend;

    double meshSize;

    bool M_use_weak_dirichlet;
    double M_gammabc;

    bool M_do_export;
    export_ptrtype exporter;

    mesh_ptrtype mesh;
    space_ptrtype Xh;
    sparse_matrix_ptrtype D,M,A_du;
    vector_ptrtype F, F_du;
    element_ptrtype pT;

    std::vector<sparse_matrix_ptrtype> M_Aq;
    std::vector<sparse_matrix_ptrtype> M_Mq;
    std::vector<sparse_matrix_ptrtype> M_Aq_du;
    std::vector<std::vector<vector_ptrtype> > M_Fq;
    std::vector<vector_ptrtype> M_Fq_du;


    parameterspace_ptrtype M_Dmu;
    theta_vector_type M_thetaAq;
    theta_vector_type M_thetaMq;
    std::vector<theta_vector_type> M_thetaFq;
    theta_vector_type M_thetaAq_du;
    theta_vector_type M_thetaFq_du;


	bdf_ptrtype M_bdf;
    bool M_fill_snapshots_matrix;
    int M_Nsnap;
    matrixN_type M_snapshots_matrix;
    matrixN_type Mdu_snapshots_matrix;

};

UnsteadyHeat1D::UnsteadyHeat1D()
    :
    backend( backend_type::build( BACKEND_PETSC ) ),
    steady( false ),
    alpha( 1 ),
    meshSize( 0.01 ),
    M_do_export( true ),
    exporter( Exporter<mesh_type>::New( "ensight" ) ),
    M_Dmu( new parameterspace_type )
{
  this->init();
}


UnsteadyHeat1D::UnsteadyHeat1D( po::variables_map const& vm )
    :
    M_vm( vm ),
    backend( backend_type::build( vm ) ),
    steady( vm["steady"].as<bool>() ),
    alpha( vm["alpha"].as<double>() ),
    meshSize( vm["hsize"].as<double>() ),
    M_do_export( !vm.count( "no-export" ) ),
    exporter( Exporter<mesh_type>::New( vm, "heat1d" ) ),
    M_Dmu( new parameterspace_type )
{
  this->init();
}
void
UnsteadyHeat1D::init()
{

    M_fill_snapshots_matrix = false;
    /*
     * First we create the mesh
     */
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=createGeo(meshSize) );


    /*
     * The function space and some associate elements are then defined
     */
    Xh = space_type::New( mesh );
    // allocate an element of Xh
    pT = element_ptrtype( new element_type( Xh ) );

    //M_bdf = bdf_ptrtype( new bdf_type(M_vm ,Xh, " " ) );
    M_bdf = bdf( _space=Xh, _vm=M_vm, _name="unsteadyHeat1d" , _prefix="unsteadyHeat1d");

    M_Aq.resize( 3 );
    M_Aq[0] = backend->newMatrix( Xh, Xh );
    M_Aq[1] = backend->newMatrix( Xh, Xh );
    M_Aq[2] = backend->newMatrix( Xh, Xh );

    M_Fq.resize( 2 );
    M_Fq[0].resize( 2 );
    M_Fq[0][0] = backend->newVector( Xh );
    M_Fq[0][1] = backend->newVector( Xh );

    M_Fq[1].resize( 1 );
    M_Fq[1][0] = backend->newVector( Xh );

    M_Mq.resize( 1 );
    M_Mq[0] = backend->newMatrix( Xh, Xh );

    M_Aq_du.resize( 3 );
    M_Aq_du[0] = backend->newMatrix( Xh, Xh );
    M_Aq_du[1] = backend->newMatrix( Xh, Xh );
    M_Aq_du[2] = backend->newMatrix( Xh, Xh );

    M_Fq_du.resize(1);
    M_Fq_du[0] = backend->newVector( Xh );

    D = backend->newMatrix( Xh, Xh );
    F = backend->newVector( Xh );
    A_du = backend->newMatrix( Xh, Xh );
    F_du = backend->newVector( Xh );

    using namespace Feel::vf;
    static const int N = 2;
    Feel::ParameterSpace<4>::Element mu_min( M_Dmu );
    mu_min << 0.2, 0.2, 0.01, 0.1;
    M_Dmu->setMin( mu_min );
    Feel::ParameterSpace<4>::Element mu_max( M_Dmu );
    mu_max << 50, 50, 5, 5;
    M_Dmu->setMax( mu_max );


    /*
    Feel::ParameterSpace<4>::Element mu_min( M_Dmu );
    mu_min << 3.5857, 0.2531,  1.4058, 0.53378;
    M_Dmu->setMin( mu_min );
    Feel::ParameterSpace<4>::Element mu_max( M_Dmu );
    mu_max << 3.5857, 0.2531,  1.4058, 0.53378;
    M_Dmu->setMax( mu_max );
    */

    element_type u( Xh, "u" );
    element_type v( Xh, "v" );

    M = backend->newMatrix( Xh, Xh );
    form2( Xh, Xh, M, _init=true ) =
        integrate( elements(mesh), id(u)*idt(v) );
    M->close();


    /*
    M = backend->newMatrix( Xh, Xh );
    form2( Xh, Xh, M, _init=true ) =
        integrate( elements(mesh), id(u)*idt(v) + grad(v)*trans(gradt(u)) );
    M->close();
    */

    Log() << "Number of dof " << Xh->nLocalDof() << "\n";

    assemble();


} // UnsteadyHeat1d::init





int
UnsteadyHeat1D::computeNumberOfSnapshots()
{
    M_Nsnap = M_bdf->timeFinal()/M_bdf->timeStep();
    return  M_Nsnap;
}




void
UnsteadyHeat1D::assemble()
{

    M_bdf->start();

    using namespace Feel::vf;

    element_type u( Xh, "u" );
    element_type v( Xh, "v" );


    //mass matrix
    form2( Xh, Xh, M_Mq[0], _init=true ) = integrate ( elements(mesh), alpha*idt(u)*id(v) );
    M_Mq[0]->close();

    // right hand side
    form1( Xh, M_Fq[0][0], _init=true ) = integrate( markedfaces(mesh,"left"), id(v) );
    form1( _test=Xh, _vector=M_Fq[0][1], _init=true ) = integrate( elements(mesh), id(v) );
    M_Fq[0][0]->close();
    M_Fq[0][1]->close();

    // output
    form1( Xh, M_Fq[1][0], _init=true ) = integrate( markedelements(mesh,"k1_2"), id(v)/0.2 );
    form1( Xh, M_Fq[1][0] ) += integrate( markedelements(mesh,"k2_1"), id(v)/0.2 );
    M_Fq[1][0]->close();

    form2( Xh, Xh, M_Aq[0], _init=true ) = integrate( elements(mesh), 0.1*(gradt(u)*trans(grad(v)) ) );
    form2( Xh, Xh, M_Aq[0] ) += integrate( markedfaces(mesh,"right"), idt(u)*id(v) );

    //if(!steady)
    //{
    //    form2( Xh, Xh, M_Aq[0] ) += integrate( elements(mesh),alpha*idt(u)*id(v)*M_bdf->polyDerivCoefficient(0) );
    //}
    M_Aq[0]->close();

    form2( Xh, Xh, M_Aq[1], _init=true ) = integrate( markedelements(mesh,"k1_1"), (gradt(u)*trans(grad(v)) ) );
    form2( Xh, Xh, M_Aq[1] ) += integrate( markedelements(mesh,"k1_2"), (gradt(u)*trans(grad(v)) ) );
    M_Aq[1]->close();

    form2( Xh, Xh, M_Aq[2], _init=true ) = integrate( markedelements(mesh,"k2_1"), (gradt(u)*trans(grad(v)) ) );
    form2( Xh, Xh, M_Aq[2] ) += integrate( markedelements(mesh,"k2_2"), (gradt(u)*trans(grad(v)) ) );
    M_Aq[2]->close();

    //and now the adjoint
    form2( Xh, Xh, M_Aq_du[0], _init=true) = integrate( elements(mesh), 0.1*(gradt(u)*trans(grad(v)) ) );
    form2( Xh, Xh, M_Aq_du[0] ) += integrate( elements(mesh),-alpha*idt(u)*id(v)*M_bdf->polyDerivCoefficient(0) );
    form2( Xh, Xh, M_Aq_du[0] ) += integrate( markedfaces(mesh,"right"), idt(u)*id(v)*(alpha-1) );
    M_Aq_du[0]->close();

    form2( Xh, Xh, M_Aq_du[1], _init=true ) = integrate( markedelements(mesh,"k1_1"), (gradt(u)*trans(grad(v)) ) );
    form2( Xh, Xh, M_Aq_du[1] ) += integrate( markedelements(mesh,"k1_2"), (gradt(u)*trans(grad(v)) ) );
    M_Aq_du[1]->close();

    form2( Xh, Xh, M_Aq_du[2], _init=true ) = integrate( markedelements(mesh,"k2_1"), (gradt(u)*trans(grad(v)) ) );
    form2( Xh, Xh, M_Aq_du[2] ) += integrate( markedelements(mesh,"k2_2"), (gradt(u)*trans(grad(v)) ) );
    M_Aq_du[2]->close();
}




UnsteadyHeat1D::sparse_matrix_ptrtype
UnsteadyHeat1D::newMatrix() const
{
    return backend->newMatrix( Xh, Xh );
}

UnsteadyHeat1D::affine_decomposition_type
UnsteadyHeat1D::computeAffineDecomposition()
{
    return boost::make_tuple( M_Mq, M_Aq, M_Fq );
}


void
UnsteadyHeat1D::exportResults(double time, element_type& T )
{
    if ( M_do_export )
    {
        Log() << "exportResults starts\n";

        exporter->step(time)->setMesh( T.functionSpace()->mesh() );

        exporter->step(time)->add( "T", T );

        exporter->save();
    }
} // Heat1d::export


void
UnsteadyHeat1D::exportResults1d( double time, element_type& T, double output)
{
    std::map<double,double> data;
    for(auto it = T.functionSpace()->mesh()->beginElement( ),
            en = T.functionSpace()->mesh()->endElement( );
        it!=en; ++it )
    {
        for( size_type i = 0; i < space_type::basis_type::nLocalDof; ++i )
        {
            value_type a = it->point(0).node()[0];
            value_type b = it->point(1).node()[0];
            value_type x = 0;
            if ( i == 0 )
                x=a;
            else if ( i == 1 )
                x=b;
            else
                x= a + (i-1)*(b-a)/(space_type::basis_type::nLocalDof-1);

            data[x] = T.localToGlobal( it->id(), i, 0);
        }
    }

    std::ostringstream fname_T;
    std::ofstream ofs3( (boost::format("Solution_%.3f") % time).str().c_str() );
    BOOST_FOREACH( auto d, data )
    {
        ofs3 <<std::setprecision(16)<< d.first << " " << d.second << "\n";
    }
    ofs3.close();


    std::ofstream output_file;
    output_file.open("Output.dat",std::ios::out | std::ios::app);
    output_file<<std::setprecision(16)<<time<<" "<<output<<"\n";

}





void
UnsteadyHeat1D::update( parameter_type const& mu,double bdf_coeff, element_type const& bdf_poly, int output_index)
{
    //first direct model
    D->close();
    D->zero();
    for( size_type q = 0;q < M_Aq.size(); ++q )
    {
        D->addMatrix( M_thetaAq[q], M_Aq[q] );
    }

    F->close();
    F->zero();

    for( size_type q = 0;q < M_Fq[output_index].size(); ++q )
    {
        F->add( M_thetaFq[output_index][q], M_Fq[output_index][q] );
    }

    auto vec_bdf_poly = backend->newVector(Xh);
    //add contribution from mass matrix
    for( size_type q = 0;q < M_Mq.size(); ++q )
    {
        //left hand side
        D->addMatrix( M_thetaMq[q]*bdf_coeff, M_Mq[q] );

        *vec_bdf_poly = bdf_poly;
        vec_bdf_poly->scale( M_thetaMq [q]);

        F->addVector( *vec_bdf_poly, *M_Mq[q]);
    }

    //now adjoint model
    A_du->close();
    A_du->zero();
    for( size_type q = 0;q < M_Aq_du.size(); ++q )
    {
        A_du->addMatrix( M_thetaAq[q], M_Aq_du[q] );
    }
    F_du->close();
    F_du->zero();

    for( size_type q = 0;q < M_Fq[output_index].size(); ++q )
    {
        //right hand side
        F_du->add( -M_thetaFq[output_index][q], M_Fq[output_index][q] );
    }
    for( size_type q = 0;q < M_Fq_du.size(); ++q )
    {
        //right hand side
        //bondaries conditions of adjoint model
        F_du->add( M_thetaFq_du[q], M_Fq_du[q] );
    }

    //add contribution from mass matrix
    for( size_type q = 0;q < M_Mq.size(); ++q )
    {
        //left hand side
        A_du->addMatrix( M_thetaMq[q]*(-bdf_coeff), M_Mq[q] );
        //right hand side
        *vec_bdf_poly = bdf_poly;
        vec_bdf_poly->scale(M_thetaMq[q]);
        F_du->addVector( *vec_bdf_poly, *M_Mq[q]);
    }



}



void
UnsteadyHeat1D::fillSnapshotsMatrix (parameter_type const& mu, int output_index)
{


    M_fill_snapshots_matrix=true;

    const int Ndof = Xh->nDof();//number of dofs used
    M_snapshots_matrix.resize(Ndof,M_Nsnap);
    Mdu_snapshots_matrix.resize(Ndof,M_Nsnap);

    //note : we need output_index for the adjoint
    solve(mu , pT, output_index);
}

void
UnsteadyHeat1D::solve( parameter_type const& mu )
{
    element_ptrtype T( new element_type( Xh ) );
    this->solve( mu, T );
}

void
UnsteadyHeat1D::solve( parameter_type const& mu, element_ptrtype& T , int output_index )
{

    using namespace vf;
    element_type v( Xh, "v" );//test functions

    pT->zero();
    assemble();

    this->computeThetaq( mu );

    M_bdf->initialize(*T);

    if(steady)
    {
        M_bdf->setSteady();
    }

    double bdf_coeff = M_bdf->polyDerivCoefficient(0);
    int column_index=0;//usefull to fill the snapshots matrix
    std::cout<<"  -- solving primal "<<std::endl;
    for ( M_bdf->start(); !M_bdf->isFinished();M_bdf->next() )
    {
        auto bdf_poly = M_bdf->polyDeriv();
        this->update( mu , bdf_coeff, bdf_poly );
        auto ret = backend->solve( _matrix=D,  _solution=T, _rhs=F, _prec=D, _reuse_prec=(M_bdf->iteration() >=2));
        if ( !ret.get<0>() )
        {
            Log()<<"WARNING : at time "<<M_bdf->time()<<" we have not converged ( nb_it : "<<ret.get<1>()<<" and residual : "<<ret.get<2>() <<" ) \n";
        }

#if( 0 )
        //this->exportResults(M_bdf->time(), *T );
        bool export_output=true;
        double s = output(output_index,mu,export_output);
        if( s != 0 && output_index==1) this->exportResults1d(M_bdf->time(),*T,s);
#endif

        if( M_fill_snapshots_matrix )
        {
            element_type solution=*T;
            for(int i=0;i<solution.size();i++)
            {
                M_snapshots_matrix(i,column_index)=solution[i];
            }
        }
        column_index++;
        M_bdf->shiftRight( *T );
    }

#if (0)
    if( M_fill_snapshots_matrix )
    {
        //and now we compute dual ( z ) solution
        // - alpha dz/dt - k /Delta z = - output
        element_ptrtype Tdu( new element_type( Xh ) );
        column_index=0;

        double dt=-M_bdf->timeStep();
        M_bdf->setTimeStep(dt);
        M_bdf->setTimeInitial(M_bdf->timeFinal());
        M_bdf->setTimeFinal(0);
        std::cout<<"  -- solving dual "<<std::endl;
        for ( M_bdf->start(); M_bdf->time()>=M_bdf->timeFinal(); M_bdf->next() )
        {
            element_type v( Xh, "v" );//test functions
            auto bdf_poly = M_bdf->polyDeriv();

            this->update( mu, bdf_coeff, bdf_poly, output_index );

            auto ret = backend->solve( _matrix=A_du,  _solution=Tdu, _rhs=F_du, _prec=A_du , _reuse_prec=0);
            if ( !ret.get<0>() )
            {
                Log()<<"ADJOINT MODEL WARNING : at time "<<M_bdf->time()<<" we have not converged ( nb_it : "<<ret.get<1>()<<" and residual : "<<ret.get<2>() <<" ) \n";
            }

            //fill the matrix
            element_type adjoint=*Tdu;
            for(int i=0;i<adjoint.size();i++)
            {
                Mdu_snapshots_matrix(i,column_index)=adjoint[i];
            }
            column_index++;
            M_bdf->shiftRight( *Tdu );
        }
        dt=-dt;
        M_bdf->setTimeStep(dt);
        M_bdf->setTimeFinal(M_bdf->timeInitial());
        M_bdf->setTimeInitial(0);
    }
#endif
}

void
UnsteadyHeat1D::l2solve( vector_ptrtype& u, vector_ptrtype const& f )
{
    //std::cout << "l2solve(u,f)\n";
    backend->solve( _matrix=M,  _solution=u, _rhs=f, _prec=M );
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
        this->init();
        do_init = false;
    }
    this->solve( mu, pT );

    double mean = integrate( elements(mesh),
                             chi( (Px() >= -0.1) && (Px() <= 0.1) )*idv(*pT) ).evaluate()(0,0)/0.2;
    Y[0]=mean;
}



//if we export output we call it in solve function so we don't need to recall it
double
UnsteadyHeat1D::output( int output_index, parameter_type const& mu, bool export_output )
{

    using namespace vf;
    if(!export_output)
    {
        this->solve( mu, pT );
    }

    using namespace vf;


    vector_ptrtype U( backend->newVector( Xh ) );
    *U = *pT;

    // right hand side (compliant)
    if( output_index == 0 )
    {
        double s1 = M_thetaFq[0](0)*dot( M_Fq[0][0], U )+M_thetaFq[0](1)*dot( M_Fq[0][1], U );
        //std::cout << "output0 c1 = " << s1 <<"\n";
        double s2 = ( M_thetaFq[0](0)*integrate( markedfaces(mesh,mesh->markerName( "left" )), idv(*pT) ).evaluate()(0,0) +
                      M_thetaFq[0](1)*integrate( elements(mesh), idv(*pT) ).evaluate()(0,0) );
        //std::cout << "output0 c2 = " << s2 <<"\n";
        return s1;
    }
    // output
    if ( output_index == 1 )
    {
        double mean = integrate( elements(mesh),
                                 chi( (Px() >= -0.1) && (Px() <= 0.1) )*idv(*pT) ).evaluate()(0,0)/0.2;
        //std::cout<<"output1 c1 = "<<mean<<std::endl;

        double meanT = ( integrate( markedelements(mesh,"k1_2"),idv(*pT) ).evaluate()(0,0)+
                         integrate( markedelements(mesh,"k2_1"),idv(*pT) ).evaluate()(0,0) )/0.2;
        //std::cout<<"output1 c2 = "<<meanT<<std::endl;
        //std::cout<<"output1 c3= "<< dot( M_Fq[1][0], U ) <<std::endl;
        return meanT;
    }

}

}

#endif /* __Heat1D_H */


