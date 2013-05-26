/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-06-18

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file mic.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-06-18
 */
#ifndef __Microphone_H
#define __Microphone_H 1

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

namespace Feel
{

po::options_description
makeMicrophoneOptions()
{
    po::options_description Microphoneoptions( "helmholtz Microphoneoptions" );
    Microphoneoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
    ( "mu1", po::value<double>()->default_value( 1 ), "Height of the microphon" )
    ( "mu2", po::value<double>()->default_value( 0.1 ), "wave speed [10;50]" )
    ( "gamma-dir", po::value<double>()->default_value( 50 ), "penalisation coefficient for weak Dirichlet condition" )
    ;
    return Microphoneoptions;
}
AboutData
makeMicrophoneAbout( std::string const& str = "mic" )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "Helmholtz microphone problem (2D)",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2011 Université de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}

class ParameterDefinition
{
public :
    static const uint16_type ParameterSpaceDimension = 2;
    typedef ParameterSpace<ParameterSpaceDimension> parameterspace_type;
};


/**
 * \class Microphone
 * @author Christophe Prud'homme
 * @see
 */
class Microphone : public ModelCrbBase< ParameterDefinition >
{
public:

    typedef ModelCrbBase<ParameterDefinition> super_type;
    typedef typename super_type::funs_type funs_type;
    typedef typename super_type::funsd_type funsd_type;


    /** @name Constants
     */
    //@{

    static const uint16_type Order = 3;
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

    /*mesh*/
    typedef Simplex<2,1> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigen_matrix_type;
    typedef eigen_matrix_type ematrix_type;
    typedef boost::shared_ptr<eigen_matrix_type> eigen_matrix_ptrtype;

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

    typedef Eigen::VectorXd beta_vector_type;

    typedef boost::tuple<
        std::vector< std::vector<sparse_matrix_ptrtype> >,
        std::vector< std::vector<std::vector<vector_ptrtype> > > ,
        std::vector< std::vector< element_ptrtype > > > affine_decomposition_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    Microphone();

    //! constructor from command line
    Microphone( po::variables_map const& vm );


    //! copy constructor
    //Microphone( Microphone const & );
    //! destructor
    ~Microphone() {}

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
     * \brief compute the beta coefficient for both bilinear and linear form
     * \param mu parameter to evaluate the coefficients
     */
    boost::tuple<beta_vector_type, std::vector<beta_vector_type> >
    computeBetaqm( parameter_type const& mu , double time=0 )
    {
        M_betaAqm.resize( Qa() );
        M_betaAqm( 0 ) = 1;
        M_betaAqm( 1 ) = -mu( 1 );
        M_betaAqm( 2 ) =  mu( 0 );
        M_betaAqm( 3 ) =1./mu( 0 );
        M_betaAqm( 4 ) = -mu( 0 )*mu( 1 );

        M_betaFqm.resize( Nl() );
        M_betaFqm[0].resize( Ql( 0 ) );
        M_betaFqm[0]( 0 ) = 1;
        M_betaFqm[0]( 1 ) = mu( 0 );

        return boost::make_tuple( M_betaAqm, M_betaFqm );
    }

    /**
     * \brief return the coefficient vector
     */
    beta_vector_type const& betaAqm() const
    {
        return M_betaAqm;
    }

    /**
     * \brief return the coefficient vector
     */
    std::vector<beta_vector_type> const& betaFqm() const
    {
        return M_betaFqm;
    }

    /**
     * \brief return the coefficient vector \p q component
     *
     */
    value_type betaAqm( int q ) const
    {
        return M_betaAqm( q );
    }

    /**
     * \return the \p q -th term of the \p l -th output
     */
    value_type betaL( int l, int q ) const
    {
        return M_betaFqm[l]( q );
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
    bool M_is_initialized;

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

    std::vector< std::vector<sparse_matrix_ptrtype> >M_Aqm;
    std::vector< std::vector<std::vector<vector_ptrtype> > > M_Fqm;
    std::vector< std::vector< element_ptrtype> > M_InitialGuessQm;

    parameterspace_ptrtype M_Dmu;
    beta_vector_type M_betaAqm;
    std::vector<beta_vector_type> M_betaFqm;
};

Microphone::Microphone()
    :
    M_is_initialized( false ),
    backend( backend_type::build( BACKEND_PETSC ) ),
    meshSize( 0.01 ),
    M_gammabc( 50 ),
    M_do_export( true ),
    exporter( Exporter<mesh_type>::New( "ensight" ) ),
    M_Dmu( new parameterspace_type )
{
    this->init();
}


Microphone::Microphone( po::variables_map const& vm )
    :
    M_is_initialized( false ),
    M_vm( vm ),
    backend( backend_type::build( vm ) ),
    meshSize( vm["hsize"].as<double>() ),
    M_gammabc( vm["gamma-dir"].as<double>() ),
    M_do_export( !vm.count( "no-export" ) ),
    exporter( Exporter<mesh_type>::New( vm, "Microphone" ) ),
    M_Dmu( new parameterspace_type )
{
    this->init();
}
void
Microphone::init()
{
    if ( M_is_initialized  )
        return;

    M_is_initialized = true;
    LOG(INFO) << "[MIC::init] starting...\n";
    // geometry is a ]0,1[x]0,1[
    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 1.5,0.25 );
    GeoTool::Node x3( 0.5,0.25 );
    GeoTool::Node x4( 1.5,0.65 );
    GeoTool::Rectangle R1( meshSize,"R1",x1,x2 );
    R1.setMarker( _type="line",_name="In",_marker4=true );
    R1.setMarker( _type="surface",_name="Bottom",_markerAll=true );
    GeoTool::Rectangle R2( meshSize,"R2",x3,x4 );
    R2.setMarker( _type="surface",_name="Top",_markerAll=true );
    mesh = ( R1+R2 ).createMesh<mesh_type>( "Omega" );

#if 0
    /*
     * The function space and some associate elements are then defined
     */
    Xh = space_type::New( mesh );
    // allocate an element of Xh
    pT = element_ptrtype( new element_type( Xh ) );

    //  initialisation de A1 et A2
    M_Aqm.resize( Qa() );
    for(int q=0; q<Qa(); q++)
        M_Aqm[q].resize( 1 );

    M_Aqm[0][0] = backend->newMatrix( Xh, Xh );
    M_Aqm[1][0] = backend->newMatrix( Xh, Xh );
    M_Aqm[2][0] = backend->newMatrix( Xh, Xh );


    M_Fqm.resize( Nl() );
    M_Fqm[0].resize( Ql( 0 ) );
    M_Fqm[0][0].resize(1);
    M_Fqm[0][0][0] = backend->newVector( Xh );

    D = backend->newMatrix( Xh, Xh );
    F = backend->newVector( Xh );

    using namespace Feel::vf;
    static const int N = 2;
    Feel::ParameterSpace<2>::Element mu_min( M_Dmu );
    mu_min << 0.8, 10;
    M_Dmu->setMin( mu_min );
    Feel::ParameterSpace<2>::Element mu_max( M_Dmu );
    mu_max << 1.2,50;
    M_Dmu->setMax( mu_max );

    element_type u( Xh, "u" );
    element_type v( Xh, "v" );

    LOG(INFO) << "Number of dof " << Xh->nLocalDof() << "\n";

    // right hand side
    form1( Xh, M_Fqm[0][0][0], _init=true ) = integrate( elements( mesh ), id( v ) );
    M_Fqm[0][0[0]->close();

    form2( Xh, Xh, M_Aqm[0][0], _init=true ) = integrate( elements( mesh ), dxt( u )*dx( v ) );
    M_Aqm[0][0]->close();

    form2( Xh, Xh, M_Aqm[1][0], _init=true ) = integrate( elements( mesh ), dyt( u )*dy( v ) );
    // Dirichlet condition apply to Bottom, only y-dir terms non zero because of normal being N()=(Nx(),Ny()) = (0,-1)
    // thus the simplification below with the signs which should -Ny() = +1
    form2( Xh, Xh, M_Aqm[1][0] ) += integrate( markedfaces( mesh,"Bottom" ), dyt( u )*id( v )+dy( u )*idt( v )+M_gammabc*idt( u )*id( v )/hFace() );
    M_Aqm[1][0]->close();

    form2( Xh, Xh, M_Aqm[2][0], _init=true ) = integrate( elements( mesh ), idt( u )*id( v ) );
    M_Aqm[2][0]->close();

    M = backend->newMatrix( Xh, Xh );

    form2( Xh, Xh, M, _init=true ) =
        integrate( elements( mesh ), id( u )*idt( v ) + grad( u )*trans( gradt( u ) ) );
    M->close();

#endif
    LOG(INFO) << "[MIC::init] done\n";

} // Microphone::run

Microphone::sparse_matrix_ptrtype
Microphone::newMatrix() const
{
    return backend->newMatrix( Xh, Xh );
}

Microphone::affine_decomposition_type
Microphone::computeAffineDecomposition()
{
    return boost::make_tuple( M_Aqm, M_Fqm , M_InitialGuessQm );
}


void
Microphone::solve( sparse_matrix_ptrtype& D,
                   element_type& u,
                   vector_ptrtype& F )
{
    backend->solve( _matrix=D, _solution=u, _rhs=F );
} // Microphone::solve


void
Microphone::exportResults( element_type& U )
{
    if ( M_do_export )
    {
        LOG(INFO) << "exportResults starts\n";

        exporter->step( 0 )->setMesh( U.functionSpace()->mesh() );

        exporter->step( 0 )->add( "u", U );

        exporter->save();
    }
} // Microphone::export

void
Microphone::update( parameter_type const& mu )
{
    *D = *M_Aqm[0];

    for ( size_type q = 1; q < M_Aqm.size(); ++q )
    {
        //std::cout << "[affine decomp] scale q=" << q << " with " << M_betaAqm[q] << "\n";
        D->addMatrix( M_betaAqm[q], M_Aqm[q] );
    }

    F->close();
    F->zero();

    for ( size_type q = 0; q < M_Fqm[0].size(); ++q )
    {
        //std::cout << "[affine decomp] scale q=" << q << " with " << M_betaFqm[0][q] << "\n";
        F->add( M_betaFqm[0][q], M_Fqm[0][q] );
    }
}
void
Microphone::solve( parameter_type const& mu )
{
    //std::cout << "solve(mu) for parameter " << mu << "\n";

    element_ptrtype T( new element_type( Xh ) );
    this->solve( mu, T );
    this->exportResults( *T );

}

void
Microphone::solve( parameter_type const& mu, element_ptrtype& T )
{
    this->computeBetaq( mu );
    this->update( mu );
    backend->solve( _matrix=D,  _solution=T, _rhs=F );
}

void
Microphone::l2solve( vector_ptrtype& u, vector_ptrtype const& f )
{
    //std::cout << "l2solve(u,f)\n";
    backend->solve( _matrix=M,  _solution=u, _rhs=f );
    //std::cout << "l2solve(u,f) done\n";
}

double
Microphone::scalarProduct( vector_ptrtype const& x, vector_ptrtype const& y )
{
    return M->energy( x, y );
}
double
Microphone::scalarProduct( vector_type const& x, vector_type const& y )
{
    return M->energy( x, y );
}

void
Microphone::run( const double * X, unsigned long N, double * Y, unsigned long P )
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


    Y[0]=M_betaFqm[0]( 0 )*integrate( elements( mesh ), idv( *pT ) ).evaluate()( 0,0 );
}



double
Microphone::output( int output_index, parameter_type const& mu )
{
    using namespace vf;
    this->solve( mu, pT );
    vector_ptrtype U( backend->newVector( Xh ) );
    *U = *pT;

    // right hand side (compliant)
    if ( output_index == 0 )
    {
        double s1 = M_betaFqm[0]( 0 )*dot( M_Fqm[0][0], U );
        return s1;
    }

}

}

#endif /* __Microphone_H */



