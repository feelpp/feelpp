/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  Author(s): Stephane Veys <stephane.veys@imag.fr>
       Date: 2014-01-19

  Copyright (C) 2011-2014 Feel++ Consortium

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
   \file benchmarkgrepl.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Stephane Veys <stephane.veys@imag.fr>

   date 2014-01-19
 */
#ifndef __BenchmarkGrepl_H
#define __BenchmarkGrepl_H 1

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

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>

#include <feel/feeldiscr/reducedbasisspace.hpp>



namespace Feel
{

po::options_description
makeBenchmarkGreplOptions()
{
    po::options_description bgoptions( "BenchmarkGrepl options" );
    bgoptions.add_options()
        ( "mshfile", Feel::po::value<std::string>()->default_value( "" ), "name of the gmsh file input")
        ( "do-export", Feel::po::value<bool>()->default_value( false ), "export results if true" )
    ;
    return bgoptions.add( Feel::feel_options() ).add( backend_options("backendl2") );
}
AboutData
makeBenchmarkGreplAbout( std::string const& str = "benchmarkGrepl" )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "Benchmark Grepl",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2011-2014 Feel++ Consortium" );

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

template<int Order>
class FunctionSpaceDefinition
{
public :
    //static const uint16_type Order = 1;
    typedef double value_type;

    /*mesh*/
    typedef Simplex<2,1> entity_type; /*dim,order*/
    typedef Mesh<entity_type> mesh_type;

    /*basis*/
    typedef bases<Lagrange<Order, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
};

template <typename ParameterDefinition, typename FunctionSpaceDefinition >
class EimDefinition
{
public :
    typedef typename ParameterDefinition::parameterspace_type parameterspace_type;
    typedef typename FunctionSpaceDefinition::space_type space_type;


    /* EIM */
    // Scalar continuous //
    typedef EIMFunctionBase<space_type, space_type , parameterspace_type> fun_type;
    // Scalar Discontinuous //
    typedef EIMFunctionBase<space_type, space_type , parameterspace_type> fund_type;

};

/**
 * \class BenchmarkGrepl
 * \brief brief description
 *
 * This is from the paper
 * EFFICIENT REDUCED-BASIS TREATMENT OF NONAFFINE
 * AND NONLINEAR PARTIAL DIFFERENTIAL EQUATIONS
 *  authors :
 * Martin A. Grepl, Yvon Maday, Ngoc C. Nguyen and Anthony T. Patera
 * ESAIM: Mathematical Modelling and Numerical Analysis
 * --
 * 3. Nonaffine linear coercive elliptic equations
 *
 * @author Christophe Prud'homme
 * @author Stephane Veys
 * @see
 */
template<int Order>
class BenchmarkGrepl : public ModelCrbBase< ParameterDefinition, FunctionSpaceDefinition<Order> , EimDefinition<ParameterDefinition, FunctionSpaceDefinition<Order> > >,
                       public boost::enable_shared_from_this< BenchmarkGrepl<Order> >
{
public:

    typedef ModelCrbBase<ParameterDefinition, FunctionSpaceDefinition<Order>, EimDefinition<ParameterDefinition,FunctionSpaceDefinition<Order> > > super_type;
    typedef typename super_type::funs_type funs_type;
    typedef typename super_type::funsd_type funsd_type;


    /** @name Constants
     */
    //@{

    //static const uint16_type Order = 1;
    static const uint16_type ParameterSpaceDimension = 2;
    static const bool is_time_dependent = false;

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
    typedef Simplex<2,1> entity_type; /*dim,order*/
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/
    typedef bases<Lagrange<Order, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef space_type functionspace_type;
    typedef space_ptrtype functionspace_ptrtype;
    typedef typename space_type::element_type element_type;
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
    typedef typename parameterspace_type::element_type parameter_type;
    typedef typename parameterspace_type::element_ptrtype parameter_ptrtype;
    typedef typename parameterspace_type::sampling_type sampling_type;
    typedef typename parameterspace_type::sampling_ptrtype sampling_ptrtype;


    typedef Eigen::VectorXd vectorN_type;
    typedef std::vector< std::vector< double > > beta_vector_type;

    typedef boost::tuple<
        std::vector< std::vector<sparse_matrix_ptrtype> >,
        std::vector< std::vector<std::vector<vector_ptrtype> > >
        > affine_decomposition_type;

    typedef Preconditioner<double> preconditioner_type;
    typedef boost::shared_ptr<preconditioner_type> preconditioner_ptrtype;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    BenchmarkGrepl();

    //! constructor from command line
    BenchmarkGrepl( po::variables_map const& vm );


    //! copy constructor
    BenchmarkGrepl( BenchmarkGrepl const & );
    //! destructor
    ~BenchmarkGrepl() {}

    //! initialization of the model
    void initModel();
    //@}

    std::string modelName()
    {
        std::ostringstream ostr;
        ostr << "BenchMarkGrepl" <<  Order;
        return ostr.str();
    }

    /** @name Operator overloads
     */
    //@{

    //@}

    /** @name Accessors
     */
    //@{

    //\return the list of EIM objects
    virtual funs_type scalarContinuousEim() const
    {
        return M_funs;
    }

    // \return the number of terms in affine decomposition of left hand
    // side bilinear form
    int Qa() const
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
        return 1;
    }

    int mMaxA( int q )
    {
        if ( q==0 )
            return 1;
        if( q==1 )
        {
            auto eim_g = M_funs[0];
            return eim_g->mMax();
        }
        else
            throw std::logic_error( "[Model Benchmark Grepl] ERROR : try to acces to mMaxA(q) with a bad value of q");
    }

    int mMaxF( int output_index, int q)
    {
        if ( q == 0 )
            return 1;
        if( q == 1 )
        {
            auto eim_g = M_funs[0];
            return eim_g->mMax();
        }
        else
            throw std::logic_error( "[Model Benchmark Grepl] ERROR : try to acces to mMaxF(output_index,q) with a bad value of q");
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
     * \brief compute the theta coefficient for both bilinear and linear form
     * \param mu parameter to evaluate the coefficients
     */
    boost::tuple<beta_vector_type,  std::vector<beta_vector_type> >
    computeBetaQm( element_type const& T,parameter_type const& mu , double time=1e30 )
    {
        return computeBetaQm( mu , time );
    }

    boost::tuple<beta_vector_type,  std::vector<beta_vector_type>  >
    computeBetaQm( parameter_type const& mu , double time=1e30 )
    {
        double mu0   = mu( 0 );
        double mu1    = mu( 1 );
        M_betaAqm.resize( Qa() );
        M_betaAqm[0].resize( 1 );
        M_betaAqm[0][0]=1;

        auto eim_g = M_funs[0];
        int M_g = eim_g->mMax();
        vectorN_type beta_g = eim_g->beta( mu );
        M_betaAqm[1].resize( M_g );
        for(int m=0; m<M_g; m++)
        {
            M_betaAqm[1][m] = beta_g(m);
        }

        M_betaFqm.resize( Nl() );
        M_betaFqm[0].resize( Ql(0) );
        M_betaFqm[0][0].resize( M_g );
        for(int m=0; m<M_g; m++)
        {
            M_betaFqm[0][0][m] = beta_g(m);
        }

        M_betaFqm[1].resize( Ql(1) );
        M_betaFqm[1][0].resize( 1 );
        M_betaFqm[1][0][0] = 1;

        return boost::make_tuple(  M_betaAqm, M_betaFqm );
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

    std::vector< std::vector< element_ptrtype > > computeInitialGuessAffineDecomposition()
    {
        std::vector< std::vector<element_ptrtype> > q;
        q.resize(1);
        q[0].resize(1);
        element_ptrtype elt ( new element_type ( Xh ) );
        q[0][0] = elt;
        return q;
    }


    void assemble();

    /**
     * solve for a given parameter \p mu
     */
    element_type solve( parameter_type const& mu );

    /**
     * solve \f$ M u = f \f$
     */
    void l2solve( vector_ptrtype& u, vector_ptrtype const& f );


    //@}

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( element_type& T , parameter_type const& mu );

    /**
     * inner product
     */
    sparse_matrix_ptrtype innerProduct ( void )
    {
        return M;
    }

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
    value_type output( int output_index, parameter_type const& mu, element_type &T, bool need_to_solve=false, bool export_outputs=false );

    parameter_type refParameter()
    {
        return M_Dmu->min();
    }

private:

    po::variables_map M_vm;

    backend_ptrtype M_backend;
    backend_ptrtype M_backendl2;

    preconditioner_ptrtype M_preconditionerl2;

    bool do_export;

    parameterspace_ptrtype M_Dmu;

    /* mesh, pointers and spaces */
    mesh_ptrtype mesh;
    space_ptrtype Xh;
    rbfunctionspace_ptrtype RbXh;

    sparse_matrix_ptrtype D,M;
    vector_ptrtype F;

    element_ptrtype pT;

    std::vector< std::vector<sparse_matrix_ptrtype> > M_Aqm;
    std::vector< std::vector<std::vector<vector_ptrtype> > > M_Fqm;

    beta_vector_type M_betaAqm;
    std::vector<beta_vector_type> M_betaFqm;

    element_type u,v;

    parameter_type M_mu;

    funs_type M_funs;

};

template<int Order>
BenchmarkGrepl<Order>::BenchmarkGrepl()
    :
    M_backend( backend_type::build( BACKEND_PETSC ) ),
    M_backendl2( backend_type::build( BACKEND_PETSC ) ),
    do_export( false ),
    M_Dmu( new parameterspace_type ),
    M_mu( M_Dmu->element() )
{ }

template<int Order>
BenchmarkGrepl<Order>::BenchmarkGrepl( po::variables_map const& vm )
    :
    M_vm( vm ),
    M_backend( backend_type::build( vm ) ),
    M_backendl2( backend_type::build( vm , "backendl2" ) ),
    do_export( option(_name="do-export").template as<bool>() ),
    M_Dmu( new parameterspace_type ),
    M_mu( M_Dmu->element() )
{
        M_preconditionerl2 = preconditioner(_pc=(PreconditionerType) M_backendl2->pcEnumType(), // by default : lu in seq or wirh mumps, else gasm in parallel
                                            _backend= M_backendl2,
                                            _pcfactormatsolverpackage=(MatSolverPackageType) M_backendl2->matSolverPackageEnumType(),// mumps if is installed ( by defaut )
                                            _worldcomm=M_backendl2->comm(),
                                            _prefix=M_backendl2->prefix() ,
                                            _rebuild=true);
}


template<int Order>
void BenchmarkGrepl<Order>::initModel()
{

    using namespace Feel::vf;

    std::string mshfile_name = option("mshfile").as<std::string>();

    /*
     * First we create the mesh or load it if already exist
     */

    if( mshfile_name=="" )
    {
        mesh = unitSquare();
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
    RbXh = rbfunctionspace_type::New( _model=this->shared_from_this() , _mesh=mesh );
    if( Environment::worldComm().isMasterRank() )
    {
        std::cout << "Number of dof " << Xh->nLocalDof() << "\n";
    }

    // allocate an element of Xh
    pT = element_ptrtype( new element_type( Xh ) );

    typename Feel::ParameterSpace<ParameterSpaceDimension>::Element mu_min( M_Dmu );
    mu_min <<  -1, -1;
    M_Dmu->setMin( mu_min );
    typename Feel::ParameterSpace<ParameterSpaceDimension>::Element mu_max( M_Dmu );
    mu_max << -0.01, -0.01;
    M_Dmu->setMax( mu_max );

    u = Xh->element();
    v = Xh->element();

    M_mu = M_Dmu->element();

    auto Pset = M_Dmu->sampling();
    //specify how many elements we take in each direction
    std::vector<int> N(2);
    //40 elements in each direction
    N[0]=40; N[1]=40;
    Pset->equidistributeProduct( N );

    auto eim_g = eim( _model=eim_no_solve(this->shared_from_this()),
                      _element=*pT,
                      _space=Xh,
                      _parameter=M_mu,
                      _expr=1./sqrt( (Px()-cst_ref(M_mu(0)))*(Px()-cst_ref(M_mu(0))) + (Py()-cst_ref(M_mu(1)))*(Py()-cst_ref(M_mu(1))) ),
                      _sampling=Pset,
                      _name="eim_g" );
    M_funs.push_back( eim_g );

    assemble();

} // BenchmarkGrepl::init


template<int Order>
void BenchmarkGrepl<Order>::assemble()
{
    using namespace Feel::vf;

    auto eim_g = M_funs[0];

    M_Aqm.resize( Qa() );
    M_Aqm[0].resize( 1 );
    M_Aqm[0][0] = M_backend->newMatrix( Xh, Xh );
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[0][0] ) =
        integrate( elements( mesh ), gradt( u )*trans( grad( v ) ) );

    int M_g = eim_g->mMax();
    M_Aqm[1].resize( M_g );
    for(int m=0; m<M_g; m++)
    {
        M_Aqm[1][m] = M_backend->newMatrix( Xh, Xh );
        form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[1][m] ) =
            integrate( elements( mesh ), idt( u )* id( v ) * idv( eim_g->q(m) ) );
    }

    M_Fqm.resize( Nl() );
    M_Fqm[0].resize( Ql(0) );
    M_Fqm[1].resize( Ql(1) );

    M_Fqm[0][0].resize(M_g);
    for(int m=0; m<M_g; m++)
    {
        M_Fqm[0][0][m] = M_backend->newVector( Xh );
        form1( Xh, M_Fqm[0][0][m] ) = integrate( elements( mesh ), id( v ) * idv( eim_g->q(m) ) );
    }
    M_Fqm[1][0].resize(1);
    M_Fqm[1][0][0] = M_backend->newVector( Xh );
    form1( Xh, M_Fqm[1][0][0] ) = integrate( elements( mesh ), id( v ) );



    //for scalarProduct
    auto mu = refParameter();
    vectorN_type beta_g = eim_g->beta( mu );

    M = M_backend->newMatrix( _test=Xh, _trial=Xh );
    form2( Xh, Xh, M ) = integrate( _range=elements( mesh ), _expr=gradt( u )*trans( grad( v ) ) );
    for(int m=0; m<M_g; m++)
    {
        auto q = eim_g->q(m);
        q.scale( beta_g(m) );
        form2( Xh, Xh, M ) +=  integrate( _range=elements( mesh ), _expr= idt( u )*id( v ) * idv( q ) );
    }

    M_preconditionerl2->setMatrix( M );

}

template<int Order>
typename BenchmarkGrepl<Order>::sparse_matrix_ptrtype
BenchmarkGrepl<Order>::newMatrix() const
{
    return M_backend->newMatrix( Xh, Xh );
}

template<int Order>
typename BenchmarkGrepl<Order>::vector_ptrtype
BenchmarkGrepl<Order>::newVector() const
{
    return M_backend->newVector( Xh );
}

template<int Order>
typename BenchmarkGrepl<Order>::affine_decomposition_type
BenchmarkGrepl<Order>::computeAffineDecomposition()
{
    return boost::make_tuple( M_Aqm, M_Fqm );
}



template<int Order>
void BenchmarkGrepl<Order>::exportResults( element_type& solution, parameter_type const& mu )
{
    LOG(INFO) << "exportResults starts\n";
    std::string exp_name = "Model_solution";
    export_ptrtype exporter;
    exporter = export_ptrtype( Exporter<mesh_type>::New( "ensight", exp_name  ) );
    exporter->step( 0 )->setMesh( solution.functionSpace()->mesh() );
    std::string mu_str;

    for ( int i=0; i<mu.size(); i++ )
    {
        mu_str= mu_str + ( boost::format( "_%1%" ) %mu[i] ).str() ;
    }

    std::string name = "solution_with_parameters_"+mu_str;
    exporter->step( 0 )->add( name, solution );
    exporter->save();
}



template<int Order>
typename BenchmarkGrepl<Order>::element_type
BenchmarkGrepl<Order>::solve( parameter_type const& mu )
{
    std::cout<<"[BenchMark solve]"<<std::endl;
    auto solution = Xh->element();
    auto exprg = 1./sqrt( (Px()-mu(0))*(Px()-mu(0)) + (Py()-mu(1))*(Py()-mu(1)) );
    auto A = M_backend->newMatrix( Xh, Xh );
    auto F = M_backend->newVector( Xh );
    form2( Xh, Xh, A ) = integrate( _range=elements( mesh ), _expr=gradt( u )*trans( grad( v ) ) + idt( u )*id( v )*exprg );
    form1( Xh, F ) = integrate( _range=elements(mesh) , _expr=id( v ) * exprg );
    M_backend->solve( _matrix=A, _solution=solution, _rhs=F );
    std::cout<<"[BenchMark solve] finished"<<std::endl;
    return solution;
}



template<int Order>
void BenchmarkGrepl<Order>::l2solve( vector_ptrtype& u, vector_ptrtype const& f )
{
    M_backendl2->solve( _matrix=M,  _solution=u, _rhs=f , _prec=M_preconditionerl2 );
}

template<int Order>
double BenchmarkGrepl<Order>::scalarProduct( vector_ptrtype const& x, vector_ptrtype const& y )
{
    return M->energy( x, y );
}

template<int Order>
double BenchmarkGrepl<Order>::scalarProduct( vector_type const& x, vector_type const& y )
{
    return M->energy( x, y );
}


template<int Order>
void BenchmarkGrepl<Order>::run( const double * X, unsigned long N, double * Y, unsigned long P )
{
    using namespace vf;
    typename Feel::ParameterSpace<ParameterSpaceDimension>::Element mu( M_Dmu );
    mu << X[0], X[1], X[2];
    static int do_init = true;

    if ( do_init )
    {
        this->initModel();
        do_init = false;
    }

    *pT = this->solve( mu );

    double mean = integrate( elements( mesh ),idv( *pT ) ).evaluate()( 0,0 );
    Y[0]=mean;
}


template<int Order>
double BenchmarkGrepl<Order>::output( int output_index, parameter_type const& mu, element_type &solution, bool need_to_solve , bool export_outputs )
{
    using namespace vf;


    if ( need_to_solve )
        *pT = this->solve( mu );
    else
        *pT = solution;

    pT->close();
    double s=0;

    if ( output_index<2 )
    {
        for ( int q=0; q<Ql( output_index ); q++ )
        {
            for ( int m=0; m<mMaxF(output_index,q); m++ )
            {
                s += M_betaFqm[output_index][q][m]*dot( *M_Fqm[output_index][q][m] , *pT );
            }
        }
    }
    else
    {
        throw std::logic_error( "[BenchmarkGrepl::output] error with output_index : only 0 or 1 " );
    }

    return s ;
}

}

#endif /* __BenchmarkGrepl_H */


