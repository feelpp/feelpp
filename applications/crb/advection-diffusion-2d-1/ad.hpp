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
   \file AdvectionDiffusion.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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

#include <feel/feelcrb/modelcrbbase.hpp>
#include <feel/feeldiscr/reducedbasisspace.hpp>

namespace Feel
{

po::options_description
makeAdvectionDiffusionOptions()
{
    po::options_description AdvectionDiffusionoptions( "AdvectionDiffusion options" );
    AdvectionDiffusionoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
    ( "mu1", po::value<double>()->default_value( 1 ), "lenght of the channel in [1;10]" )
    ( "mu2", po::value<double>()->default_value( 0.1 ), "Peclet number in [0.1;100]" )
    ( "no-export", "don't export results" )
    ;
    return AdvectionDiffusionoptions;
}
AboutData
makeAdvectionDiffusionAbout( std::string const& str = "AdvectionDiffusion" )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "2D steady Advection-Diffusion",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2011 Université de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}

class FunctionSpaceDefinition
{
public :

    typedef double value_type;
    static const uint16_type Order = 5;

    /*mesh*/
    typedef Simplex<2,1> entity_type;
    typedef Mesh<entity_type> mesh_type;

    /*basis*/
    typedef bases<Lagrange<Order, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;

    /*elements*/
    typedef typename space_type::element_type element_type;

    static const bool is_time_dependent = false;
    static const bool is_linear = true;

};
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
class AdvectionDiffusion : public ModelCrbBase<ParameterSpace<2>,FunctionSpaceDefinition>
{
public:

    typedef ModelCrbBase<ParameterSpace<2>,FunctionSpaceDefinition> super_type;
    typedef typename super_type::funs_type funs_type;
    typedef typename super_type::funsd_type funsd_type;

    /** @name Constants
     */
    //@{

    static const uint16_type Order = 5;
    static const uint16_type ParameterSpaceDimension = 2;
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
    typedef bases<Lagrange<Order, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef space_type functionspace_type;
    typedef space_ptrtype functionspace_ptrtype;
    typedef space_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;

    //typedef ReducedBasisSpace<super_type, mesh_type, basis_type, value_type> rbfunctionspace_type;
    //typedef boost::shared_ptr< rbfunctionspace_type > rbfunctionspace_ptrtype;

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


    typedef std::vector< std::vector< double > > beta_vector_type;

    typedef boost::tuple<
        std::vector< std::vector<sparse_matrix_ptrtype> >,
        std::vector< std::vector<std::vector<vector_ptrtype> > >
        > affine_decomposition_type;

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
    ~AdvectionDiffusion() {}

    //! initialisation of the model
    void initModel();
    //@}

    /** @name Accessors
     */
    //@{

    // \return the number of terms in affine decomposition of left hand
    // side bilinear form
    int Qa() const
    {
        return 5;
    }

    /**
     * there is at least one output which is the right hand side of the
     * primal problem
     *
     * \return number of outputs associated to the model
     */
    int Nl() const
    {
        return 1;
    }

    /**
     * \param l the index of output
     * \return number of terms  in affine decomposition of the \p q th output term
     */
    int Ql( int l ) const
    {
        return 1;
    }

    int mMaxA( int q )
    {
        if ( q < Qa() )
            return 1;
        else
            throw std::logic_error( "[Model] ERROR : try to acces to mMaxA(q) with a bad value of q");
    }

    int mMaxF( int output_index, int q)
    {
        if ( q < 1 )
            return 1;
        else
            throw std::logic_error( "[Model] ERROR : try to acces to mMaxF(output_index,q) with a bad value of q");
    }

    /**
     * \brief compute the beta coefficient for both bilinear and linear form
     * \param mu parameter to evaluate the coefficients
     */

    betaqm_type
    computeBetaQm(element_type const& T,  parameter_type const& mu )
    {
        return computeBetaQm( mu );
    }

    betaqm_type
    computeBetaQm( parameter_type const& mu )
    {
        M_betaAqm.resize( Qa() );
        for(int i=0; i<Qa(); i++)
            M_betaAqm[i].resize(1);
        M_betaAqm[0][0] = 1;
        M_betaAqm[1][0] = 1./( mu( 0 )*mu( 1 ) );
        M_betaAqm[2][0] = mu( 0 )/mu( 1 );
        M_betaAqm[3][0] = 1./( mu( 0 )*mu( 0 )*mu( 1 ) );
        M_betaAqm[4][0] = 1./mu( 1 );

        M_betaFqm.resize( Nl() );
        M_betaFqm[0].resize( Ql( 0 ) );
        M_betaFqm[0][0].resize(1);
        M_betaFqm[0][0][0] = mu( 0 );

        return boost::make_tuple( M_betaAqm, M_betaFqm );
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
     * create a new matrix
     * \return the newly created matrix
     */
    sparse_matrix_ptrtype newMatrix() const;

    /**
     * create a new vector
     * \return the newly created vector
     */
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
    void solve( parameter_type const& mu, element_ptrtype& T );

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
    void update( parameter_type const& mu );
    //@}

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( element_type& u , parameter_type const& mu );

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
    value_type output( int output_index, parameter_type const& mu , element_type &u, bool need_to_solve=false);


    parameter_type refParameter()
    {
        return Dmu->min();
    }

private:
    po::variables_map M_vm;

    double meshSize;

    bool M_do_export;

    int export_number;

    bool M_use_weak_dirichlet;

    double M_gammabc;

    mesh_ptrtype mesh;
    sparse_matrix_ptrtype D,M;
    vector_ptrtype F;
    element_ptrtype pT;

    std::vector< std::vector<sparse_matrix_ptrtype> > M_Aqm;
    std::vector< std::vector<std::vector<vector_ptrtype> > > M_Fqm;

    beta_vector_type M_betaAqm;
    std::vector<beta_vector_type> M_betaFqm;

};

AdvectionDiffusion::AdvectionDiffusion()
    :
    meshSize( 0.01 ),
    M_do_export( true ),
    export_number( 0 )
{}


AdvectionDiffusion::AdvectionDiffusion( po::variables_map const& vm )
    :
    M_vm( vm ),
    meshSize( vm["hsize"].as<double>() ),
    M_do_export( !vm.count( "no-export" ) ),
    export_number( 0 )
{}
void
AdvectionDiffusion::initModel()
{
    // geometry is a ]0,1[x]0,1[
    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 1,1 );
    GeoTool::Rectangle R( meshSize,"Omega",x1,x2 );
    R.setMarker( _type="line",_name="Inflow",_marker4=true );
    R.setMarker( _type="line",_name="Bottom",_marker1=true );
    R.setMarker( _type="line",_name="Top",_marker3=true );
    R.setMarker( _type="line",_name="Outflow",_marker2=true );
    R.setMarker( _type="surface",_name="Omega",_markerAll=true );
    mesh = R.createMesh( _mesh=new mesh_type, _name="Omega" );

    /*
     * The function space and some associate elements are then defined
     */
    auto Xh = space_type::New( mesh );
    this->setFunctionSpaces( Xh );
    if (Environment::isMasterRank() )
        std::cout << "Number of dof : "<< Xh->nDof() << std::endl;
    // allocate an element of Xh
    pT = element_ptrtype( new element_type( Xh ) );

    //  initialisation de A1 et A2
    M_Aqm.resize( Qa() );
    for(int q=0; q<Qa(); q++)
    {
        M_Aqm[q].resize( 1 );
        M_Aqm[q][0] = backend()->newMatrix( Xh, Xh );
    }

    M_Fqm.resize( this->Nl() );
    for(int l=0; l<Nl(); l++)
    {
        M_Fqm[l].resize( Ql(l) );
        for(int q=0; q<Ql(l) ; q++)
        {
            M_Fqm[l][q].resize(1);
            M_Fqm[l][q][0] = backend()->newVector( Xh );
        }
    }

    D = backend()->newMatrix( Xh, Xh );
    F = backend()->newVector( Xh );

    auto mu_min = Dmu->element();
    mu_min << 1, 0.1;
    Dmu->setMin( mu_min );
    auto mu_max = Dmu->element();
    mu_max << 10,100;
    Dmu->setMax( mu_max );

    element_type u( Xh, "u" );
    element_type v( Xh, "v" );

    std::cout << "Number of dof " << Xh->nLocalDof() << "\n";

    // right hand side
    form1( _test=Xh, _vector=M_Fqm[0][0][0] )
        = integrate( markedfaces( mesh, "Bottom" ), id( v ) );
    M_Fqm[0][0][0]->close();

    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[0][0] )
        = integrate( elements( mesh ), Py()*dxt( u )*id( v ) );
    M_Aqm[0][0]->close();

    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[1][0] )
        = integrate( elements( mesh ), dxt( u )*dx( v ) );
    M_Aqm[1][0]->close();

    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[2][0] )
        = integrate( elements( mesh ), dyt( u )*dy( v ) );
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[2][0] )
        += integrate( markedfaces( mesh,"Top" ),
                      - dyt( u )*Ny()*id( v ) - dy( u )*Ny()*idt( v ) + 20*idt( u )*id( v )/hFace() );
    M_Aqm[2][0]->close();

    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[3][0] )
        = integrate( markedfaces( mesh,"Inflow" ),
                     - dxt( u )*Nx()*id( v ) - dx( u )*Nx()*idt( v ) );
    M_Aqm[3][0]->close();

    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[4][0] )
        = integrate( markedfaces( mesh,"Inflow" ),
                     20*idt( u )*id( v )/hFace() );
    M_Aqm[4][0]->close();

    M = backend()->newMatrix( Xh, Xh );

    form2( _test=Xh, _trial=Xh, _matrix=M )
        = integrate( elements( mesh ), id( u )*idt( v ) + grad( u )*trans( gradt( u ) ) );
    M->close();

} // AdvectionDiffusion::run

AdvectionDiffusion::sparse_matrix_ptrtype
AdvectionDiffusion::newMatrix() const
{
    return backend()->newMatrix( Xh, Xh );
}

AdvectionDiffusion::vector_ptrtype
AdvectionDiffusion::newVector() const
{
    return backend()->newVector( Xh );
}

AdvectionDiffusion::affine_decomposition_type
AdvectionDiffusion::computeAffineDecomposition()
{
    return boost::make_tuple( M_Aqm, M_Fqm  );
}


void
AdvectionDiffusion::solve( sparse_matrix_ptrtype& D,
                           element_type& u,
                           vector_ptrtype& F )
{
    backend()->solve( _matrix=D, _solution=u, _rhs=F );
} // AdvectionDiffusion::solve


void
AdvectionDiffusion::exportResults( element_type& U , parameter_type const& mu )
{

    if ( M_do_export )
    {
        LOG(INFO) << "exportResults starts\n";

        std::string exp_name;
        export_ptrtype exporter;
        std::string mu_str;

        for ( int i=0; i<mu.size(); i++ )
        {
            mu_str= mu_str + ( boost::format( "_%1%" ) %mu[i] ).str() ;
        }

        exp_name = "solution_with_parameters_" + mu_str;

        exporter = export_ptrtype( Exporter<mesh_type>::New( "ensight", exp_name  ) );
        exporter->step( 0 )->setMesh( U.functionSpace()->mesh() );
        exporter->step( 0 )->add( "u", U );
        exporter->save();
    }
} // AdvectionDiffusion::export

void
AdvectionDiffusion::update( parameter_type const& mu )
{
    *D = *M_Aqm[0][0];

    for ( size_type q = 1; q < M_Aqm.size(); ++q )
    {
        for ( size_type m = 0; m < mMaxA(q); ++m )
        {
            D->addMatrix( M_betaAqm[q][m] , M_Aqm[q][m] );
            //std::cout << "[affine decomp] scale q=" << q << " with " << M_betaAqm[q][0] << "\n";
        }
    }

    F->close();
    F->zero();

    for ( size_type q = 0; q < M_Fqm[0].size(); ++q )
    {
        for ( size_type m = 0; m < mMaxF(0,q); ++m )
        {
            //std::cout << "[affine decomp] scale q=" << q << " with " << M_betaFqm[0][q][0] << "\n";
            F->add( M_betaFqm[0][q][m], M_Fqm[0][q][m] );
        }
    }
}

AdvectionDiffusion::element_type
AdvectionDiffusion::solve( parameter_type const& mu )
{
    //std::cout << "solve(mu) for parameter " << mu << "\n";

    element_ptrtype T( new element_type( Xh ) );
    this->solve( mu, T );
    return *T;
    //this->exportResults( *T );

}

void
AdvectionDiffusion::solve( parameter_type const& mu, element_ptrtype& T )
{
    this->computeBetaQm( mu );
    this->update( mu );
    backend()->solve( _matrix=D,  _solution=T, _rhs=F );
    export_number++;
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


    std::cout<<"pfem solution ok"<<std::endl;
    std::ofstream file_matrix;
    name = "PFEMmatrix" + mu_str + number;
    file_matrix.open( name,std::ios::out );
    file_matrix<<*D;
    file_matrix.close();

    std::cout<<"pfem matrix ok"<<std::endl;
    name = "PFEMrhs" + mu_str + number;
    std::ofstream file_rhs;
    file_rhs.open( name,std::ios::out );
    file_rhs<<*F;
    file_rhs.close();
    std::cout<<"pfem rhs ok"<<std::endl;
#endif

}

void
AdvectionDiffusion::l2solve( vector_ptrtype& u, vector_ptrtype const& f )
{
    //std::cout << "l2solve(u,f)\n";
    backend()->solve( _matrix=M,  _solution=u, _rhs=f );
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
    auto mu = Dmu->element();
    mu << X[0], X[1];
    static int do_init = true;

    if ( do_init )
    {
        meshSize = X[2];
        this->initModel();
        do_init = false;
    }

    this->solve( mu, pT );


    Y[0]=M_betaFqm[0][0][0]*integrate( markedfaces( mesh, "Bottom" ), idv( *pT ) ).evaluate()( 0,0 );
}



double
AdvectionDiffusion::output( int output_index, parameter_type const& mu, element_type &u, bool need_to_solve )
{
    using namespace vf;
    if( need_to_solve )
        this->solve( mu, pT );
    else
        *pT = u;

    double output=0;

    // right hand side (compliant)
    if ( output_index == 0 )
    {
        //output = M_betaFqm[0][0][0]*dot( M_Fqm[0][0][0], U );
        for ( int q=0; q<Ql( output_index ); q++ )
        {
            for ( int m=0; m<mMaxF(output_index,q); m++ )
            {
                //element_ptrtype eltF( new element_type( Xh ) );
                //*eltF = *M_Fqm[output_index][q][m];
                //output += M_betaFqm[output_index][q][m]*dot( *eltF, *pT );
                output += M_betaFqm[output_index][q][m]*dot( *M_Fqm[output_index][q][m] , *pT );
            }
        }
    }

    else
    {
        throw std::logic_error( "[AdvectionDiffusion::output] error with output_index : only 0 " );
    }

    return output;
}

}

#endif /* __AdvectionDiffusion_H */
