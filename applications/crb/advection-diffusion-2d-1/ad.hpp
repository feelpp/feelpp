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
   \file AdvectionDiffusion.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-11-13
 */
#ifndef __AdvectionDiffusion_H
#define __AdvectionDiffusion_H 1

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
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelcrb/parameterspace.hpp>


namespace Feel
{

po::options_description
makeAdvectionDiffusionOptions()
{
    po::options_description AdvectionDiffusionoptions("AdvectionDiffusion options");
    AdvectionDiffusionoptions.add_options()
        ("hsize", po::value<double>()->default_value( 0.1 ), "mesh size")
        ("mu1", po::value<double>()->default_value( 1 ), "lenght of the channel in [1;10]")
        ("mu2", po::value<double>()->default_value( 0.1 ), "Peclet number in [0.1;100]")
        ("no-export", "don't export results")
        ;
    return AdvectionDiffusionoptions.add( Feel::feel_options() );
}
AboutData
makeAdvectionDiffusionAbout( std::string const& str = "AdvectionDiffusion" )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "2D steady Advection-Diffusion",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2011 Université de Grenoble 1 (Joseph Fourier)");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;
}

/**
 * \class AdvectionDiffusion
 * \brief 2D Advection-Diffusion problem
 *
 * We consider an advection-diffusion example in a rectangular domain \f$\Omega
 * (\mu) = ]0, L[ × ]0, 1[\f$ representing a chan- nel. The governing equation
 * for the passive-scalar field (say, temperature) is the advection-diffusion
 * equation with im- posed Couette velocity \f$(y,0)\f$.
 *
 * Neumann (flux) boundary conditions are imposed on the bottom wall
 * \f$\Gamma_{bot}\f$ ; homogeneous Dirichlet conditions are imposed on the top
 * wall \f$\Gamma_{top}\f$ and on the “inflow” left boundary \f$\Gamma_{in}\f$ ;
 * and homogeneous (zero flux) Neumann conditions are imposed on the “outflow”
 * right boundary \f$\Gamma_{out}\f$ .
 *
 * The output of interest is the integral of the temperature over the heated
 * (bottom) surface \f$\Gamma_{bot}\f$ . This example is a simplified version of
 * a Couette- Graetz problem [1].  We consider two parameters: the length of the
 * channel, \f$L\f$, and the Peclet number, \f$Pe\f$ [2]. Hence \f$P = 2\f$ and
 * \f$\mu = (\mu_1 , \mu_2 )\f$: \f$\mu_1\f$ is the channel length \f$L\f$, and
 * \f$\mu_2\f$ is the Peclet number \f$Pe\f$; the parameter domain is given by
 * \f$D = [1, 10] \times [0.1, 100]\f$. We now choose \f$\mu_{ref} = (1, 1)\f$,
 * which in turn defines the reference domain \f$\Omega = \Omega (\mu_{ref}
 * )\f$.
 *
 *
 * 1.Schiesser WE, Silebi CA (1997) Computational transport phenomena:
 * numerical methods for the solution of transport problems. Cambridge
 * University Press, Cambridge
 *
 * 2. Gunzburger MD (1989) Finite element methods for viscous incompressible
 * flows. Academic Press, San Diego
 *
 *
 * @author Christophe Prud'homme
 * @see
 */
class AdvectionDiffusion
{
public:


    /** @name Constants
     */
    //@{

    static const uint16_type Order = 5;
    static const uint16_type ParameterSpaceDimension = 2;
    static const bool is_time_dependent = false;
    //@}

    /** @name Typedefs
     */
    //@{

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    typedef backend_type::sparse_matrix_type sparse_matrix_type;
    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef backend_type::vector_type vector_type;
    typedef backend_type::vector_ptrtype vector_ptrtype;


    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigen_matrix_type;
    typedef eigen_matrix_type ematrix_type;
    typedef boost::shared_ptr<eigen_matrix_type> eigen_matrix_ptrtype;


    /*mesh*/
    typedef Simplex<2,1> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

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

    typedef Eigen::VectorXd theta_vector_type;


    typedef boost::tuple<std::vector<sparse_matrix_ptrtype>, std::vector<std::vector<vector_ptrtype>  > > affine_decomposition_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    AdvectionDiffusion();

    //! constructor from command line
    AdvectionDiffusion( po::variables_map const& vm );


    //! copy constructor
    //AdvectionDiffusion( AdvectionDiffusion const & );
    //! destructor
    ~AdvectionDiffusion(){}

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
    int Qa() const { return 5; }

    /**
     * there is at least one output which is the right hand side of the
     * primal problem
     *
     * \return number of outputs associated to the model
     */
    int Nl() const { return 1; }

    /**
     * \param l the index of output
     * \return number of terms  in affine decomposition of the \p q th output term
     */
    int Ql( int l ) const { return 1; }

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
    boost::tuple<theta_vector_type, std::vector<theta_vector_type> >
    computeThetaq( parameter_type const& mu , double time=0)
        {
            M_thetaAq.resize( Qa() );
            M_thetaAq( 0 ) = 1;
            M_thetaAq( 1 ) = 1./(mu(0)*mu(1));
            M_thetaAq( 2 ) = mu(0)/mu(1);
            M_thetaAq( 3 ) = 1./(mu(0)*mu(0)*mu(1));
            M_thetaAq( 4 ) = 1./mu(1);

            M_thetaFq.resize( Nl() );
            M_thetaFq[0].resize( Ql(0) );
            M_thetaFq[0]( 0 ) = mu(0);

            return boost::make_tuple( M_thetaAq, M_thetaFq );
        }

    /**
     * \brief return the coefficient vector
     */
    theta_vector_type const& thetaAq() const { return M_thetaAq; }

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
     */
    void solve( parameter_type const& mu, element_ptrtype& T );

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
    void update( parameter_type const& mu );
    //@}

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( element_type& u );

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
    value_type output( int output_index, parameter_type const& mu );

private:

    po::variables_map M_vm;
    backend_ptrtype backend;

    double meshSize;

    bool M_use_weak_dirichlet;
    double M_gammabc;

    bool M_do_export;
    export_ptrtype exporter;

    mesh_ptrtype mesh;
    space_ptrtype Xh;
    sparse_matrix_ptrtype D,M;
    vector_ptrtype F;
    element_ptrtype pT;

    std::vector<sparse_matrix_ptrtype> M_Aq;
    std::vector<std::vector<vector_ptrtype> > M_Fq;

    parameterspace_ptrtype M_Dmu;
    theta_vector_type M_thetaAq;
    std::vector<theta_vector_type> M_thetaFq;
};

AdvectionDiffusion::AdvectionDiffusion()
    :
    backend( backend_type::build( BACKEND_PETSC ) ),
    meshSize( 0.01 ),
    M_do_export( true ),
    exporter( Exporter<mesh_type>::New( "ensight" ) ),
    M_Dmu( new parameterspace_type )
{
  this->init();
}


AdvectionDiffusion::AdvectionDiffusion( po::variables_map const& vm )
    :
    M_vm( vm ),
    backend( backend_type::build( vm ) ),
    meshSize( vm["hsize"].as<double>() ),
    M_do_export( !vm.count( "no-export" ) ),
    exporter( Exporter<mesh_type>::New( vm, "AdvectionDiffusion" ) ),
    M_Dmu( new parameterspace_type )
{
  this->init();
}
void
AdvectionDiffusion::init()
{
    // geometry is a ]0,1[x]0,1[
    GeoTool::Node x1(0,0);
    GeoTool::Node x2(1,1);
    GeoTool::Rectangle R( meshSize,"Omega",x1,x2);
    R.setMarker(_type="line",_name="Inflow",_marker4=true);
    R.setMarker(_type="line",_name="Bottom",_marker1=true);
    R.setMarker(_type="line",_name="Top",_marker3=true);
    R.setMarker(_type="line",_name="Outflow",_marker2=true);
    R.setMarker(_type="surface",_name="Omega",_markerAll=true);
    mesh = R.createMesh<mesh_type>("Omega");


    /*
     * The function space and some associate elements are then defined
     */
    Xh = space_type::New( mesh );
    // allocate an element of Xh
    pT = element_ptrtype( new element_type( Xh ) );

    //  initialisation de A1 et A2
    M_Aq.resize( Qa() );
    M_Aq[0] = backend->newMatrix( Xh, Xh );
    M_Aq[1] = backend->newMatrix( Xh, Xh );
    M_Aq[2] = backend->newMatrix( Xh, Xh );
    M_Aq[3] = backend->newMatrix( Xh, Xh );
    M_Aq[4] = backend->newMatrix( Xh, Xh );


    M_Fq.resize( Nl() );
    M_Fq[0].resize( Ql(0) );
    M_Fq[0][0] = backend->newVector( Xh );

    D = backend->newMatrix( Xh, Xh );
    F = backend->newVector( Xh );

    using namespace Feel::vf;
    static const int N = 2;
    Feel::ParameterSpace<2>::Element mu_min( M_Dmu );
    mu_min << 1, 0.1;
    M_Dmu->setMin( mu_min );
    Feel::ParameterSpace<2>::Element mu_max( M_Dmu );
    mu_max << 10,100;
    M_Dmu->setMax( mu_max );

    element_type u( Xh, "u" );
    element_type v( Xh, "v" );

    Log() << "Number of dof " << Xh->nLocalDof() << "\n";

    // right hand side
    form1( Xh, M_Fq[0][0], _init=true ) = integrate( markedfaces(mesh, "Bottom"), id(v) );
    M_Fq[0][0]->close();

    form2( Xh, Xh, M_Aq[0], _init=true ) = integrate( elements(mesh), Py()*dxt(u)*id(v));
    M_Aq[0]->close();

    form2( Xh, Xh, M_Aq[1], _init=true ) = integrate( elements(mesh), dxt(u)*dx(v));
    M_Aq[1]->close();

    form2( Xh, Xh, M_Aq[2], _init=true ) = integrate( elements(mesh), dyt(u)*dy(v));
    form2( Xh, Xh, M_Aq[2] ) += integrate( markedfaces(mesh,"Top"),
                                           - dyt(u)*Ny()*id(v) - dy(u)*Ny()*idt(v) + 20*idt(u)*id(v)/hFace() );
    M_Aq[2]->close();

    form2( Xh, Xh, M_Aq[3], _init=true ) = integrate( markedfaces(mesh,"Inflow"),
                                                      - dxt(u)*Nx()*id(v) - dx(u)*Nx()*idt(v) );
    M_Aq[3]->close();
    form2( Xh, Xh, M_Aq[4], _init=true ) = integrate( markedfaces(mesh,"Inflow"),
                                                      20*idt(u)*id(v)/hFace() );
    M_Aq[4]->close();
    M = backend->newMatrix( Xh, Xh );

    form2( Xh, Xh, M, _init=true ) =
        integrate( elements(mesh), id(u)*idt(v) + grad(u)*trans(gradt(u)) );
    M->close();

} // AdvectionDiffusion::run

AdvectionDiffusion::sparse_matrix_ptrtype
AdvectionDiffusion::newMatrix() const
{
    return backend->newMatrix( Xh, Xh );
}

AdvectionDiffusion::affine_decomposition_type
AdvectionDiffusion::computeAffineDecomposition()
{
    return boost::make_tuple( M_Aq, M_Fq );
}


void
AdvectionDiffusion::solve( sparse_matrix_ptrtype& D,
                           element_type& u,
                           vector_ptrtype& F )
{
    backend->solve( _matrix=D, _solution=u, _rhs=F );
} // AdvectionDiffusion::solve


void
AdvectionDiffusion::exportResults( element_type& U )
{
    if ( M_do_export )
    {
        Log() << "exportResults starts\n";

        exporter->step(0)->setMesh( U.functionSpace()->mesh() );

        exporter->step(0)->add( "u", U );

        exporter->save();
    }
} // AdvectionDiffusion::export

void
AdvectionDiffusion::update( parameter_type const& mu )
{
    *D = *M_Aq[0];
    for( size_type q = 1;q < M_Aq.size(); ++q )
    {
        //std::cout << "[affine decomp] scale q=" << q << " with " << M_thetaAq[q] << "\n";
        D->addMatrix( M_thetaAq[q], M_Aq[q] );
    }
    F->close();
    F->zero();
    for( size_type q = 0;q < M_Fq[0].size(); ++q )
    {
        //std::cout << "[affine decomp] scale q=" << q << " with " << M_thetaFq[0][q] << "\n";
        F->add( M_thetaFq[0][q], M_Fq[0][q] );
    }
}
void
AdvectionDiffusion::solve( parameter_type const& mu )
{
    //std::cout << "solve(mu) for parameter " << mu << "\n";

    element_ptrtype T( new element_type( Xh ) );
    this->solve( mu, T );
    this->exportResults( *T );

}

void
AdvectionDiffusion::solve( parameter_type const& mu, element_ptrtype& T )
{
    this->computeThetaq( mu );
    this->update( mu );
    backend->solve( _matrix=D,  _solution=T, _rhs=F, _prec=D );
}

void
AdvectionDiffusion::l2solve( vector_ptrtype& u, vector_ptrtype const& f )
{
    //std::cout << "l2solve(u,f)\n";
    backend->solve( _matrix=M,  _solution=u, _rhs=f, _prec=M );
    //std::cout << "l2solve(u,f) done\n";
}

double
AdvectionDiffusion::scalarProduct( vector_ptrtype const& x, vector_ptrtype const& y )
{
    return M->energy( x, y );
}
double
AdvectionDiffusion::scalarProduct( vector_type const& x, vector_type const& y )
{
    return M->energy( x, y );
}

void
AdvectionDiffusion::run( const double * X, unsigned long N, double * Y, unsigned long P )
{
    //std::cout<<"RUN ::::::::::: "<<std::endl;
    using namespace vf;
    Feel::ParameterSpace<2>::Element mu( M_Dmu );
    mu << X[0], X[1];
    static int do_init = true;
    if ( do_init )
    {
        meshSize = X[2];
        this->init();
        do_init = false;
    }
    this->solve( mu, pT );


    Y[0]=M_thetaFq[0](0)*integrate( markedfaces(mesh, "Bottom"), idv(*pT) ).evaluate()(0,0);
}



double
AdvectionDiffusion::output( int output_index, parameter_type const& mu )
{
    using namespace vf;
    this->solve( mu, pT );
    vector_ptrtype U( backend->newVector( Xh ) );
    *U = *pT;

    // right hand side (compliant)
    if( output_index == 0 )
    {
        double s1 = M_thetaFq[0](0)*dot( M_Fq[0][0], U );
        return s1;
    }

}

}

#endif /* __AdvectionDiffusion_H */


