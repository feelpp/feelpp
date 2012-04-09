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
   \file eeg.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-11-13
 */
#ifndef __EEG_H
#define __EEG_H 1

#include <feel/options.hpp>
#include <feel/feelcore/application.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feeldiscr/operatorlinear.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelpoly/polynomialset.hpp>

#include <feel/feelvf/vf.hpp>

#include <feel/feelcrb/parameterspace.hpp>

namespace Feel
{
using namespace vf;

po::options_description
makeEEGOptions()
{
    po::options_description eegoptions( "EEG options" );
    eegoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.01 ), "mesh size" )
    /*   ("white", po::value<double>()->default_value( 0.14 ), "white")
        ("whiterad", po::value<double>()->default_value( 0.65 ), "whiterad")
        ("whitetang", po::value<double>()->default_value( 0.065 ), "whitetang")
        ("gray", po::value<double>()->default_value( 0.33 ), "gray")
        ("csf", po::value<double>()->default_value( 1.79 ), "csf")
        ("skullrad", po::value<double>()->default_value( 0.0042 ), "skullrad")
        ("skulltang", po::value<double>()->default_value( 0.042 ), "skulltang")
        ("scalp", po::value<double>()->default_value( 0.33 ), "scalp")
        ("no-export", "don't export results")*/
    ;
    return eegoptions.add( Feel::feel_options() );
}
AboutData
makeEEGAbout( std::string const& str = "eeg" )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "EEG model",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2010,2011 Université de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Sylvain Vallaghé", "developer", "sylvain.vallaghe@gmail.com", "" );
    return about;
}




/**
 * \class EEG
 * \brief brief description
 *
 * @author Sylvain Vallaghé
 * @see
 */
class EEG
{
public:


    /** @name Constants
     */
    //@{

    static const uint16_type Order = 1;
    static const uint16_type Dim = 3;
    static const uint16_type ParameterSpaceDimension = 8;
    static const uint16_type nbtissue = 8;
    static const bool is_time_dependent = false;
    //@}

    /** @name Typedefs
     */
    //@{

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;


    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigen_matrix_type;
    typedef eigen_matrix_type ematrix_type;
    typedef boost::shared_ptr<eigen_matrix_type> eigen_matrix_ptrtype;


    /*mesh*/
    typedef Simplex<Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef fusion::vector<Lagrange<Order, Scalar>,
            Lagrange<0, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> functionspace_type;
    typedef functionspace_type space_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef functionspace_ptrtype space_ptrtype;

    typedef typename functionspace_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;
    typedef typename element_type::sub_element<0>::type element_0_type;
    typedef typename element_type::sub_element<1>::type element_1_type;

    //typedef FunctionSpace<mesh_type, fusion::vector<Lagrange<0, Scalar> >, Discontinuous>
    typedef FunctionSpace<mesh_type, fusion::vector<Lagrange<0, Scalar, Discontinuous> > > p0_space_type;
    typedef boost::shared_ptr<p0_space_type> p0_space_ptrtype;
    typedef typename p0_space_type::element_type p0_element_type;


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
    EEG();

    //! constructor from command line
    EEG( po::variables_map const& vm );


    //! copy constructor
    //EEG( EEG const & );
    //! destructor
    ~EEG() {}

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
        return nbtissue+1;
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
        if ( l==0 ) return nbtissue+2;

        return 0;
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

    // return matrix of the functionspace scalar product
    sparse_matrix_ptrtype ScalarProductMatrix() const
    {
        return M;
    }

    /**
     * \brief compute the theta coefficient for both bilinear and linear form
     * \param mu parameter to evaluate the coefficients
     */
    boost::tuple<theta_vector_type, std::vector<theta_vector_type> >
    computeThetaq( parameter_type const& mu , double time=0 )
    {
        std::cout << "compute thetaq for mu " << mu << "\n" ;
        M_thetaAq.resize( Qa() );

        for ( int i=0; i<nbtissue; i++ ) M_thetaAq( i ) = mu ( i );

        M_thetaAq( nbtissue ) = 1.0;

        M_thetaFq.resize( Nl() );
        M_thetaFq[0].resize( Ql( 0 ) );

        for ( int i=0; i<nbtissue; i++ ) M_thetaFq[0]( i ) = 1.0-mu( i )/mu( grey );

        M_thetaFq[0]( nbtissue ) = 1.0 ;
        M_thetaFq[0]( nbtissue+1 ) = 1.0/mu( grey ) ;

        return boost::make_tuple( M_thetaAq, M_thetaFq );
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
    void update( parameter_type const& mu , double time=0 );
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
    backend_ptrtype M_backend;

    bool M_do_export;
    export_ptrtype exporter;

    mesh_ptrtype mesh;
    space_ptrtype Xh;
    sparse_matrix_ptrtype D,M;
    vector_ptrtype F;
    element_ptrtype pT;
    element_ptrtype ginf ;

    std::vector<sparse_matrix_ptrtype> Aq;
    std::vector<std::vector<vector_ptrtype> > Lq;

    parameterspace_ptrtype M_Dmu;
    theta_vector_type M_thetaAq;
    std::vector<theta_vector_type> M_thetaFq;

    enum Tissue {white, whiterad, whitetang, grey , csf, skullrad, skulltang, scalp } ;
    unsigned int domain[nbtissue];

};

EEG::EEG()
    :
    M_backend( backend_type::build( BACKEND_PETSC ) ),
    M_do_export( true ),
    exporter( Exporter<mesh_type>::New( "ensight" ) ),
    mesh( new mesh_type ),
    M_Dmu( new parameterspace_type )
{
    this->init();
}


EEG::EEG( po::variables_map const& vm )
    :
    M_vm( vm ),
    M_backend( backend_type::build( vm ) ),
    M_do_export( !vm.count( "no-export" ) ),
    exporter( Exporter<mesh_type>::New( vm, "eeg" ) ),
    mesh( new mesh_type ),
    M_Dmu( new parameterspace_type )
{
    this->init();
}
void
EEG::init()
{

    // Different head tissue domains

    domain[scalp]=7;
    domain[skullrad]=8;
    domain[skulltang]=8;
    domain[csf]=9;
    domain[grey]=10;
    domain[whiterad]=12;
    domain[whitetang]=12;
    domain[white]=11;

    /*
     * Loading mesh
     */
    ImporterGmsh<mesh_type> import( "real.msh" );
    mesh->accept( import );
    mesh->setComponents( MESH_CHECK | MESH_UPDATE_FACES );
    mesh->updateForUse();

    // Loading anisotropy tensor

    p0_space_ptrtype P0h = p0_space_type::New( mesh );
    p0_element_type k1( P0h, "k1" );
    p0_element_type k2( P0h, "k2" );
    p0_element_type k3( P0h, "k3" );
    Log() << "print space P0 info\n";
    P0h->printInfo();

    Log() << "Taille K " << k1.size() << "\n" ;
    std::cout << "Lecture fichier anisotropie\n" ;
    std::ifstream inFile( "diraniso" );

    if ( !inFile )
    {
        std::cout << std::endl << "Failed to open normals file " << std::endl ;
    }

    else
        for ( size_type i=0; i<k1.size(); ++i )
        {
            inFile >> k1( i ) ;
            inFile >> k2( i ) ;
            inFile >> k3( i ) ;
        }

    auto K = vec( idv( k1 ),idv( k2 ),idv( k3 ) )*trans( vec( idv( k1 ),idv( k2 ),idv( k3 ) ) ) ;
    //auto K = vec(constant(0),constant(0),constant(1))*trans(vec(constant(0),constant(0),constant(1))) ;


    // Defining homogeneous potential

    const double qx=0.0,qy=0.0,qz=-1.0,px=124.496, py=170.731, pz=77.292;
    auto r0 = vec( constant( px ),constant( py ),constant( pz ) );
    auto q = vec( constant( qx ),constant( qy ),constant( qz ) );
    auto posrel = P()-r0 ;
    auto posrelsc = trans( posrel )*posrel ;
    Log() << "posrelsc" << integrate ( elements( mesh ) , posrelsc ).evaluate()( 0, 0 ) << "\n" ;
    auto momentposrelsc = trans( q )*posrel ;
    Log() << "momentposrelsc" << integrate ( elements( mesh ) , momentposrelsc ).evaluate()( 0, 0 ) << "\n" ;
    auto g = 1.0/( 4*M_PI )*( momentposrelsc )/vf::pow( posrelsc,1.5 ) ;
    //auto g = constant(0.0) ;
    auto grad_g = 1.0/( 4*M_PI*vf::pow( posrelsc,1.5 ) )*q-( 3.0*momentposrelsc )/( 4*M_PI*vf::pow( posrelsc,2.5 ) )*posrel ;
    //auto grad_g = vec(constant(0.0),constant(0.0),constant(0.0)) ;



    /*
     * The function space and some associate elements are then defined
     */
    Xh = space_type::New( mesh );
    // allocate an element of Xh
    pT = element_ptrtype( new element_type( Xh ) );
    ginf = element_ptrtype( new element_type( Xh ) );
    ginf->element<0>() = project( Xh->functionSpace<0>(), elements( mesh ), g );

    D = M_backend->newMatrix( Xh, Xh );
    F = M_backend->newVector( Xh );

    using namespace Feel::vf;
    Feel::ParameterSpace<8>::Element mu_min( M_Dmu );
    mu_min << 0.07, 0.325, 0.0325, 0.165, 0.895, 0.0021, 0.021, 0.165;
    M_Dmu->setMin( mu_min );
    Feel::ParameterSpace<8>::Element mu_max( M_Dmu );
    mu_max << 0.28, 1.3, 0.13, 0.66, 3.58, 0.0084, 0.084, 0.66;
    M_Dmu->setMax( mu_max );

    element_type U( Xh, "u" );
    element_type V( Xh, "v" );
    // trial functions
    element_0_type u = U.element<0>() ;
    element_1_type lambda = U.element<1>() ;
    // test functions
    element_0_type v = V.element<0>() ;
    element_1_type nu = V.element<1>() ;

    Log() << "Number of dof " << Xh->nLocalDof() << "\n";

    //form2( Xh, Xh, D,_init=true )=integrate(elements(mesh), 0*idt(v)*id(v),_Q<0>());
    //D->close();
    form1( Xh, F,_init=true )=integrate( elements( mesh ), 0*id( v ),_Q<0>() );
    F->close();


    Aq.resize( nbtissue+1 );

    for ( int i=0; i<nbtissue; i++ )
    {
        Log() << "assembling Aq[" << i << "]\n" ;
        Aq[i]=M_backend->newMatrix( Xh, Xh );

        if ( i==skullrad || i==whiterad ) form2( Xh, Xh, Aq[i], _init=true ) = integrate( markedelements( mesh, domain[i] ), gradt( u )*K*trans( grad( v ) ) ) ;

        else if ( i==skulltang || i==whitetang ) form2( Xh, Xh, Aq[i], _init=true ) = integrate( markedelements( mesh, domain[i] ), gradt( u )*trans( grad( v ) )-gradt( u )*K*trans( grad( v ) ) ) ;

        else form2( Xh, Xh, Aq[i], _init=true ) = integrate( markedelements( mesh, domain[i] ), gradt( u )*trans( grad( v ) ) ) ;

        Aq[i]->close();
    }

    Aq[nbtissue]=M_backend->newMatrix( Xh, Xh );
    form2( Xh, Xh, Aq[nbtissue], _init=true ) = integrate( elements( mesh ),
            id( v )*idt( lambda ) + idt( u )*id( nu ) + 0*idt( lambda )*id( nu ) ) ;
    Aq[nbtissue]->close();

    double area = integrate( elements( mesh ), constant( 1.0 ) ).evaluate()( 0, 0 );
    double mean = integrate( elements( mesh ),-g, _Q<5>() ).evaluate()( 0, 0 )/area;
    //const double meangana = integrate( elements(mesh), gana, _Q<5>() ).evaluate()( 0, 0)/area;
    //auto ganabase = gana - constant(meangana) ;
    //Log() << " int gana " << meangana << "\n" ;
    //Log() << "int g  = " << mean << "\n";


    Lq.resize( 1 );
    Lq[0].resize( nbtissue+2 );
    vector_ptrtype F( M_backend->newVector( Xh ) );

    for ( int i=0; i<nbtissue; i++ )
    {
        Log() << "assembling Lq[" << i << "]\n" ;
        Lq[0][i]=M_backend->newVector( Xh );

        if ( i==skullrad || i==whiterad ) form1( Xh, Lq[0][i], _init=true ) = integrate( markedelements( mesh, domain[i] ),  grad( v )*K*grad_g, _Q<5>() ) ;

        else if ( i==skulltang || i==whitetang ) form1( Xh, Lq[0][i], _init=true ) = integrate( markedelements( mesh, domain[i] ), grad( v )*grad_g-grad( v )*K*grad_g, _Q<5>() ) ;

        else form1( Xh, Lq[0][i], _init=true ) = integrate( markedelements( mesh, domain[i] ),  grad( v )*grad_g, _Q<5>() ) ;

        Lq[0][i]->close();
    }

    Lq[0][nbtissue]=M_backend->newVector( Xh );
    form1( Xh, Lq[0][nbtissue], _init=true ) = integrate( boundaryfaces( mesh ), -trans( grad_g )*N()*id( v ),_Q<5>() ) ;
    Lq[0][nbtissue]->close();
    Lq[0][nbtissue+1]=M_backend->newVector( Xh );
    form1( Xh, Lq[0][nbtissue+1], _init=true ) = integrate( elements( mesh ),  mean*id( nu ) ) ;
    Lq[0][nbtissue+1]->close();

    M = M_backend->newMatrix( Xh, Xh );
    form2( Xh, Xh, M, _init=true ) =
        integrate( elements( mesh ), id( u )*idt( v ) + grad( u )*trans( gradt( u ) ) + idt( lambda )*id( nu ) );
    M->close();


    Log() << "Assembly done\n" ;

}

EEG::sparse_matrix_ptrtype
EEG::newMatrix() const
{
    return M_backend->newMatrix( Xh, Xh );
}

EEG::affine_decomposition_type
EEG::computeAffineDecomposition()
{
    return boost::make_tuple( Aq, Lq );
}


void
EEG::solve( sparse_matrix_ptrtype& D,
            element_type& u,
            vector_ptrtype& F )
{

    vector_ptrtype U( M_backend->newVector( u.functionSpace() ) );
    M_backend->solve( D, D, U, F );
    u = *U;
} // Heat1d::solve


void
EEG::exportResults( element_type& U )
{
    if ( M_do_export )
    {
        Log() << "exportResults starts\n";
        exporter->step( 0 )->setMesh( U.functionSpace()->mesh() );

        exporter->step( 0 )->add( "u", U.element<0>() );

        exporter->save();
    }
} // Heat1d::export

void
EEG::update( parameter_type const& mu )
{

    std::cout << "update for parameter " << mu << "\n" ;
    *D = *Aq[nbtissue];
    std::cout << "Merging Aq\n" ;

    for ( size_type q = 0; q < Aq.size()-1; ++q )
    {
        std::cout << "[affine decomp] scale q=" << q << " with " << M_thetaAq[q] << "\n";
        D->addMatrix( M_thetaAq[q], Aq[q] );
    }

    std::cout << "Merging Aq done\n" ;
    F->zero();
    std::cout << "Merging Lq\n" ;

    for ( size_type q = 0; q < Lq[0].size(); ++q )
    {
        std::cout << "[affine decomp] scale q=" << q << " with " << M_thetaFq[0][q] << "\n";
        F->add( M_thetaFq[0][q], Lq[0][q] );
    }

    std::cout << "Merging Lq done\n" ;

}
void
EEG::solve( parameter_type const& mu )
{
    std::cout << "solve(mu) for parameter " << mu << "\n";

    element_ptrtype T( new element_type( Xh ) );
    this->solve( mu, T );
    element_0_type u = T->element<0>() ;
    element_0_type g = ginf->element<0>() ;
    element_type E( Xh, "e" );
    element_0_type e = E.element<0>();
    e = vf::project( Xh->functionSpace<0>(), elements( mesh ), idv( u )+1.0/mu( grey )*idv( g ) );
    this->exportResults( E );

}

void
EEG::solve( parameter_type const& mu, element_ptrtype& T )
{
    this->computeThetaq( mu );
    this->update( mu );
    Log() << "Solving system starts\n" ;
    M_backend->solve( _matrix=D,  _solution=T, _rhs=F );
    Log() << "Solving system done\n" ;
}

void
EEG::l2solve( vector_ptrtype& u, vector_ptrtype const& f )
{
    //std::cout << "l2solve(u,f)\n";
    M_backend->solve( _matrix=M,  _solution=u, _rhs=f );
    //std::cout << "l2solve(u,f) done\n";
}

double
EEG::scalarProduct( vector_ptrtype const& x, vector_ptrtype const& y )
{
    return M->energy( x, y );
}
double
EEG::scalarProduct( vector_type const& x, vector_type const& y )
{
    return M->energy( x, y );
}

void
EEG::run( const double * X, unsigned long N, double * Y, unsigned long P )
{
    //std::cout<<"RUN ::::::::::: "<<std::endl;
    using namespace vf;
    Feel::ParameterSpace<8>::Element mu( M_Dmu );
    mu << X[0], X[1], X[2], X[3], X[4], X[5], X[6], X[7];
    static int do_init = true;

    if ( do_init )
    {
        this->init();
        do_init = false;
    }

    this->solve( mu, pT );

    double mean = 0.0 ;
    Y[0]=mean;
}

double
EEG::output( int output_index, parameter_type const& mu )
{
    using namespace vf;
    this->solve( mu, pT );
    vector_ptrtype U( M_backend->newVector( Xh ) );
    *U = *pT;

    return 0;
}





}

#endif /* __EEG_H */



