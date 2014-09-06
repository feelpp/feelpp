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
   \file aw.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-06-18
 */
#ifndef __AnisotropicWavespeed_H
#define __AnisotropicWavespeed_H 1

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
makeAnisotropicWavespeedOptions()
{
    po::options_description AnisotropicWavespeedoptions( "helmholtz Anisotropic Wavespeed options" );
    AnisotropicWavespeedoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
    ( "mu1", po::value<double>()->default_value( 1 ), "diffusion in y direction [0.8;1.2]" )
    ( "mu2", po::value<double>()->default_value( 0.1 ), "wave speed [10;50]" )
    ( "gamma-dir", po::value<double>()->default_value( 50 ), "penalisation coefficient for weak Dirichlet condition" )
    ;
    return AnisotropicWavespeedoptions;
}
AboutData
makeAnisotropicWavespeedAbout( std::string const& str = "aw" )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "Helmholtz anisotropic wavespeed (2D)",
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

class FunctionSpaceDefinition
{
public :
    static const uint16_type Order = 3;

    typedef double value_type;

    /*mesh*/
    typedef Simplex<2,1> entity_type;
    typedef Mesh<entity_type> mesh_type;

    /*basis*/
    typedef fusion::vector<Lagrange<Order, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;

    static const bool is_time_dependent = false;
    static const bool is_linear = true;

};

/**
 * \class AnisotropicWavespeed
 * @author Christophe Prud'homme
 * @see
 */
class AnisotropicWavespeed : public ModelCrbBase< ParameterDefinition , FunctionSpaceDefinition > ,
                             public boost::enable_shared_from_this< AnisotropicWavespeed >
{
public:

    typedef ModelCrbBase<ParameterDefinition, FunctionSpaceDefinition > super_type;
    typedef typename super_type::funs_type funs_type;
    typedef typename super_type::funsd_type funsd_type;


    /** @name Constants
     */
    //@{

    static const uint16_type Order = 3;
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
    typedef fusion::vector<Lagrange<Order, Scalar> > basis_type;

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
    AnisotropicWavespeed();

    //! constructor from command line
    AnisotropicWavespeed( po::variables_map const& vm );


    //! copy constructor
    //AnisotropicWavespeed( AnisotropicWavespeed const & );
    //! destructor
    ~AnisotropicWavespeed() {}

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
    boost::tuple<beta_vector_type, std::vector<beta_vector_type> >
    computeBetaQm( element_type const& T,parameter_type const& mu , double time=1e30 )
    {
        return computeBetaQm( mu , time );
    }

    boost::tuple<beta_vector_type, std::vector<beta_vector_type> >
    computeBetaQm( parameter_type const& mu , double time=0 )
    {
        M_betaAqm.resize( Qa() );
        M_betaAqm[0].resize( 1 );
        M_betaAqm[1].resize( 1 );
        M_betaAqm[2].resize( 1 );
        M_betaAqm[0][0] = 1;
        M_betaAqm[1][0] = mu( 0 );
        M_betaAqm[2][0] = -mu( 1 );

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
     * run the convergence test
     */

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
    value_type output( int output_index, parameter_type const& mu ,  element_type &T, bool need_to_solve=false);


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
    rbfunctionspace_ptrtype RbXh;
    sparse_matrix_ptrtype D,M;
    vector_ptrtype F;
    element_ptrtype pT;

    std::vector< std::vector<sparse_matrix_ptrtype> > M_Aqm;
    std::vector< std::vector<std::vector<vector_ptrtype> > > M_Fqm;

    parameterspace_ptrtype M_Dmu;
    beta_vector_type M_betaAqm;
    std::vector<beta_vector_type> M_betaFqm;

};

AnisotropicWavespeed::AnisotropicWavespeed()
    :
    M_is_initialized( false ),
    backend( backend_type::build( BACKEND_PETSC ) ),
    meshSize( 0.01 ),
    M_gammabc( 50 ),
    M_do_export( true ),
    exporter( Exporter<mesh_type>::New( "ensight" ) ),
    M_Dmu( new parameterspace_type )
{
}


AnisotropicWavespeed::AnisotropicWavespeed( po::variables_map const& vm )
    :
    M_is_initialized( false ),
    M_vm( vm ),
    backend( backend_type::build( vm ) ),
    meshSize( vm["hsize"].as<double>() ),
    M_gammabc( vm["gamma-dir"].as<double>() ),
    M_do_export( !vm.count( "no-export" ) ),
    exporter( Exporter<mesh_type>::New( vm, "AnisotropicWavespeed" ) ),
    M_Dmu( new parameterspace_type )
{
}
void
AnisotropicWavespeed::initModel()
{
    if ( M_is_initialized  )
        return;

    M_is_initialized = true;
    LOG(INFO) << "[AW::init] starting...\n";
    // geometry is a ]0,1[x]0,1[
    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 1,1 );
    GeoTool::Rectangle R( meshSize,"Omega",x1,x2 );
    R.setMarker( _type="line",_name="Top",_marker4=true,_marker2=true,_marker3=true );
    R.setMarker( _type="line",_name="Bottom",_marker1=true );
    R.setMarker( _type="surface",_name="Omega",_markerAll=true );
    mesh = R.createMesh( _mesh=mesh, _name="Omega" );

    /*
     * The function space and some associate elements are then defined
     */
    Xh = space_type::New( mesh );
    RbXh = rbfunctionspace_type::New( _model=this->shared_from_this() , _mesh=mesh );
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
    M_Fqm[0][0][0]->close();

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

    LOG(INFO) << "[AW::init] done\n";

} // AnisotropicWavespeed::run

AnisotropicWavespeed::sparse_matrix_ptrtype
AnisotropicWavespeed::newMatrix() const
{
    return backend->newMatrix( Xh, Xh );
}

AnisotropicWavespeed::vector_ptrtype
AnisotropicWavespeed::newVector() const
{
    return backend->newVector( Xh );
}


AnisotropicWavespeed::affine_decomposition_type
AnisotropicWavespeed::computeAffineDecomposition()
{
    return boost::make_tuple( M_Aqm, M_Fqm  );
}


void
AnisotropicWavespeed::solve( sparse_matrix_ptrtype& D,
                             element_type& u,
                             vector_ptrtype& F )
{
    backend->solve( _matrix=D, _solution=u, _rhs=F );
} // AnisotropicWavespeed::solve


void
AnisotropicWavespeed::exportResults( element_type& U )
{
    if ( M_do_export )
    {
        LOG(INFO) << "exportResults starts\n";

        exporter->step( 0 )->setMesh( U.functionSpace()->mesh() );

        exporter->step( 0 )->add( "u", U );

        exporter->save();
    }
} // AnisotropicWavespeed::export

void
AnisotropicWavespeed::update( parameter_type const& mu )
{
    *D = *M_Aqm[0][0];

    for ( size_type q = 1; q < M_Aqm.size(); ++q )
    {
        for ( size_type m = 0; m < mMaxA(q); ++m )
        {
            D->addMatrix( M_betaAqm[q][m] , M_Aqm[q][m] );
            //std::cout << "[affine decomp] scale q=" << q << " with " << M_betaAqm[q] << "\n";
        }
    }

    F->close();
    F->zero();

    for ( size_type q = 0; q < M_Fqm[0].size(); ++q )
    {
        for ( size_type m = 0; m < mMaxF(0,q); ++m )
        {
            //std::cout << "[affine decomp] scale q=" << q << " with " << M_betaFqm[0][q] << "\n";
            F->add( M_betaFqm[0][q][m], M_Fqm[0][q][m] );
        }
    }
}

typename AnisotropicWavespeed::element_type
AnisotropicWavespeed::solve( parameter_type const& mu )
{
    //std::cout << "solve(mu) for parameter " << mu << "\n";

    element_ptrtype T( new element_type( Xh ) );
    this->solve( mu, T );
    this->exportResults( *T );
    return *T;
}

void
AnisotropicWavespeed::solve( parameter_type const& mu, element_ptrtype& T )
{
    this->computeBetaQm( mu );
    this->update( mu );
    backend->solve( _matrix=D,  _solution=T, _rhs=F );
}

void
AnisotropicWavespeed::l2solve( vector_ptrtype& u, vector_ptrtype const& f )
{
    //std::cout << "l2solve(u,f)\n";
    backend->solve( _matrix=M,  _solution=u, _rhs=f );
    //std::cout << "l2solve(u,f) done\n";
}

double
AnisotropicWavespeed::scalarProduct( vector_ptrtype const& x, vector_ptrtype const& y )
{
    return M->energy( x, y );
}
double
AnisotropicWavespeed::scalarProduct( vector_type const& x, vector_type const& y )
{
    return M->energy( x, y );
}

void
AnisotropicWavespeed::run( const double * X, unsigned long N, double * Y, unsigned long P )
{
    //std::cout<<"RUN ::::::::::: "<<std::endl;
    using namespace vf;
    Feel::ParameterSpace<2>::Element mu( M_Dmu );
    mu << X[0], X[1];
    static int do_init = true;

    if ( do_init )
    {
        meshSize = X[2];
        this->initModel();
        do_init = false;
    }

    this->solve( mu, pT );


    Y[0]=M_betaFqm[0][0][0]*integrate( elements( mesh ), idv( *pT ) ).evaluate()( 0,0 );
}



double
AnisotropicWavespeed::output( int output_index, parameter_type const& mu ,  element_type &T, bool need_to_solve )
{
    using namespace vf;
    if( need_to_solve )
        this->solve( mu, pT );
    else
        *pT = T;
    //vector_ptrtype U( backend->newVector( Xh ) );
    //*U = *pT;
}

}




#endif /* __AnisotropicWavespeed_H */



