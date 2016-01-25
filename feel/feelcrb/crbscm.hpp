/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-08-07

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
   \file opusscm.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-08-07
 */
#ifndef __CRBSCM_H
#define __CRBSCM_H 1

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/timer.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/format.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/foreach.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/parameter.hpp>
#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelcrb/parameterspace.hpp>
#include <feel/feelcrb/crbdb.hpp>
#include <feel/feelcore/serialization.hpp>
#if defined(FEELPP_HAS_GLPK_H)
#include <glpk.h>
#endif /* FEELPP_HAS_GLPK_H */


namespace Feel
{
/**
 * \class CRBSCM
 * \brief SCM algorithm
 * \anchor crbscm
 *
 * \tparam TruthModelType the type of the class that describes the truth (e.g. fem) model problem
 *
 * TruthModelType must support a certain interface: affine decomposition, solve
 * for fm problem and problem associated with with affine decomposition. As to
 * the SCM, it should problem
 *  - eigensolves for the full Truth problem
 *  - eigensolves associated with the affine decomposition
 *
 * @author Christophe Prud'homme
 * @see page \ref scm
 */
template<typename TruthModelType>
class CRBSCM : public CRBDB
{
    typedef CRBDB super;
public:


    /** @name Constants
     */
    //@{


    //@}

    /** @name Typedefs
     */
    //@{

    typedef TruthModelType truth_model_type;
    typedef truth_model_type model_type;
    typedef boost::shared_ptr<truth_model_type> truth_model_ptrtype;

    typedef double value_type;
    typedef boost::tuple<double,double> bounds_type;

    typedef ParameterSpace<TruthModelType::ParameterSpaceDimension> parameterspace_type;
    typedef boost::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef typename parameterspace_type::element_type parameter_type;
    typedef typename parameterspace_type::element_ptrtype parameter_ptrtype;
    typedef typename parameterspace_type::sampling_type sampling_type;
    typedef typename parameterspace_type::sampling_ptrtype sampling_ptrtype;

    typedef boost::tuple<double, parameter_type, size_type, double, double> relative_error_type;

    typedef typename model_type::backend_type backend_type;
    typedef typename model_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename model_type::vector_ptrtype vector_ptrtype;
    typedef typename model_type::beta_vector_type beta_vector_type;


    typedef Eigen::VectorXd y_type;
    typedef Eigen::VectorXd vector_type;
    typedef std::vector< std::vector<y_type> > y_set_type;
    typedef std::vector< std::vector<boost::tuple<double,double> > > y_bounds_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    CRBSCM()
        :
        super(),
        M_is_initialized( false ),
        M_model(),
        M_tolerance( 1e-2 ),
        M_iter_max( 3 ),
        M_Mplus( 10 ),
        M_Malpha( 4 ),
        M_Dmu( new parameterspace_type ),
        M_Xi( M_Dmu ),
        M_C( M_Dmu, 1, M_Xi ),
        M_C_complement( M_Dmu, 1, M_Xi ),
        M_scm_for_mass_matrix( false ),
        M_mu_ref( M_Dmu->element() ),
        M_use_scm( true )
    {
    }

    //! constructor from command line options
    CRBSCM( std::string const& name,
            po::variables_map const& vm )
        :
        super( "scm",
               ( boost::format( "%1%" ) % name ).str(),
               ( boost::format( "%1%" ) % name ).str(),
               vm ),
        M_is_initialized( false ),
        M_model(),
        M_tolerance( vm["crb.scm.tol"].template as<double>() ),
        M_iter_max( vm["crb.scm.iter-max"].template as<int>() ),
        M_Mplus( vm["crb.scm.Mplus"].template as<int>() ),
        M_Malpha( vm["crb.scm.Malpha"].template as<int>()  ),
        M_Dmu( new parameterspace_type ),
        M_Xi( new sampling_type( M_Dmu ) ),
        M_C( new sampling_type( M_Dmu, 1, M_Xi ) ),
        M_C_complement( new sampling_type( M_Dmu, 1, M_Xi ) ),
        M_vm( vm ),
        M_scm_for_mass_matrix( false ),
        M_mu_ref( M_Dmu->element() ),
        M_use_scm( vm["crb.scm.use-scm"].template as<bool>() )
    {
        if ( this->loadDB() )
            LOG( INFO ) << "Database " << this->lookForDB() << " available and loaded";
    }


    //! constructor from command line options
    CRBSCM( std::string const& name,
            po::variables_map const& vm ,
            truth_model_ptrtype const & model )
        :
        super( "scm",
               ( boost::format( "%1%" ) % name ).str(),
               ( boost::format( "%1%" ) % name ).str(),
               vm ),
        M_is_initialized( false ),
        M_model(),
        M_tolerance( vm["crb.scm.tol"].template as<double>() ),
        M_iter_max( vm["crb.scm.iter-max"].template as<int>() ),
        M_Mplus( vm["crb.scm.Mplus"].template as<int>() ),
        M_Malpha( vm["crb.scm.Malpha"].template as<int>()  ),
        M_Dmu( new parameterspace_type ),
        M_Xi( new sampling_type( M_Dmu ) ),
        M_C( new sampling_type( M_Dmu, 1, M_Xi ) ),
        M_C_complement( new sampling_type( M_Dmu, 1, M_Xi ) ),
        M_vm( vm ),
        M_scm_for_mass_matrix( false ),
        M_mu_ref( M_Dmu->element() ),
        M_use_scm( vm["crb.scm.use-scm"].template as<bool>() )
    {
        this->setTruthModel( model );
        if ( this->loadDB() )
            LOG( INFO ) << "Database " << this->lookForDB() << " available and loaded";
    }

    //! constructor from command line options
    CRBSCM( std::string const& name,
            po::variables_map const& vm ,
            truth_model_ptrtype const & model ,
            bool scm_for_mass_matrix )
        :
        super( "scm",
               ( boost::format( "%1%" ) % name ).str(),
               ( boost::format( "%1%" ) % name ).str(),
               vm ),
        M_is_initialized( false ),
        M_model(),
        M_tolerance( vm["crb.scm.tol"].template as<double>() ),
        M_iter_max( vm["crb.scm.iter-max"].template as<int>() ),
        M_Mplus( vm["crb.scm.Mplus"].template as<int>() ),
        M_Malpha( vm["crb.scm.Malpha"].template as<int>()  ),
        M_Dmu( new parameterspace_type ),
        M_Xi( new sampling_type( M_Dmu ) ),
        M_C( new sampling_type( M_Dmu, 1, M_Xi ) ),
        M_C_complement( new sampling_type( M_Dmu, 1, M_Xi ) ),
        M_vm( vm ),
        M_scm_for_mass_matrix( scm_for_mass_matrix ),
        M_mu_ref( M_Dmu->element() ),
        M_use_scm( vm["crb.scm.use-scm"].template as<bool>() )
    {
        this->setTruthModel( model );
        if ( this->loadDB() )
            LOG( INFO ) << "Database " << this->lookForDB() << " available and loaded";
    }

    //! copy constructor
    CRBSCM( CRBSCM const & o )
        :
        super( o ),
        M_is_initialized( o.M_is_initialized ),
        M_tolerance( o.M_tolerance ),
        M_iter_max( o.M_iter_max ),
        M_Mplus( o.M_Mplus ),
        M_Malpha( o.M_Malpha ),
        M_Dmu( o.M_Dmu ),
        M_Xi( o.M_Xi ),
        M_C( o.M_C ),
        M_C_complement( o.M_C_complement ),
        M_vm( o.M_vm ),
        M_scm_for_mass_matrix( o.M_scm_for_mass_matrix ),
        M_mu_ref( o.M_mu_ref ),
        M_use_scm( o.M_use_scm )
    {
    }

    //! destructor
    ~CRBSCM()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    CRBSCM& operator=( CRBSCM const & o )
    {
        if ( this != &o )
        {
        }

        return *this;
    }
    //@}

    /** @name Accessors
     */
    //@{

    //! Returns maximum value of K
    size_type KMax() const
    {
        return M_C->size();
    }

    //! return max iterations
    int maxIter() const
    {
        return M_iter_max;
    }

    //! return the parameter space
    parameterspace_ptrtype Dmu() const
    {
        return M_Dmu;
    }

    //! get M+
    int Mplus() const
    {
        return M_Mplus;
    }

    //! get Malpha
    int Malpha()
    {
        return M_Malpha;
    }

    //@}

    /** @name  Mutators
     */
    //@{

    //! set offline tolerance
    void setTolerance( double tolerance )
    {
        M_tolerance = tolerance;
    }

    //! set M+
    void setMplus( int Mplus )
    {
        M_Mplus = Mplus;
    }

    //! set Malpha
    void setMalpha( int Malpha )
    {
        M_Malpha = Malpha;
    }

    //! set the truth offline model
    void setTruthModel( truth_model_ptrtype const& model )
    {
        M_model = model;
        M_Dmu = M_model->parameterSpace();

        if ( ! loadDB() )
        {
            M_Xi = sampling_ptrtype( new sampling_type( M_Dmu ) );
            M_C = sampling_ptrtype( new sampling_type( M_Dmu ) );
            M_C_complement = sampling_ptrtype( new sampling_type( M_Dmu ) );
        }
    }

    //! set Kmax
    void setMaxIter( int K )
    {
        M_iter_max = K;
    }

    //! set bool ( if we do scm for mass matrix or not )
    void setScmForMassMatrix( bool b )
    {
        M_scm_for_mass_matrix = b;
    }

    int mMax( int q ) const
    {
        if ( M_scm_for_mass_matrix )
            return M_model->mMaxM( q );
        else
            return M_model->mMaxA( q );
    }

    int nb_decomposition_terms_q ( void ) const
    {
        if ( M_scm_for_mass_matrix )
            return M_model->Qm();
        else
            return M_model->Qa();
    }

    //total number of decomposition terms
    int nb_decomposition_terms_qm ( void ) const
    {
        int count=0;
        for(int q=0;q<nb_decomposition_terms_q();q++)
        {
            for(int m=0;m<mMax(q);m++)
                count++;
        }
        return count;
    }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * generate offline the space \f$C_K\f$
     */
    void generate();

    /**
     * check \f$C_K\f$ properties
     * \param K dimension of \f$C_K\f$
     */
    void checkC( size_type K ) const;

    /**
     * check if eigen value and eigen vector given are solution
     * of the generalized eigenvalue problem
     * A w = \lambda B w
     * where w is an eigenvector and \lambda is an eigenvalue
     * \param A : matrix
     * \param B : matrix
     * \param eigenvector : eigenvector (like w)
     * \param eigenvalue : eigenvalue (like \lambda)
     */
    void checkEigenVectorEigenValue( sparse_matrix_ptrtype const& A, sparse_matrix_ptrtype const& B, vector_ptrtype const& eigenvector, double eigenvalue ) const;

    boost::tuple<value_type,value_type> ex( parameter_type const& mu ) const;

    /**
     * Returns the lower bound of the coercive constant given a parameter \p
     * \f$\mu\f$
     *
     * \param mu \f$ \mu\f$ the parameter at which to evaluate the coercive constant
     * \param K the dimension of \f$Y_{\mathrm{UB}}\f$
     *
     *\return compute online the lower bounmd
     */
    boost::tuple<value_type,value_type> lb( parameter_type const& mu, size_type K = invalid_size_type_value, int indexmu = -1 ) const;
    boost::tuple<value_type,value_type> lbSCM( parameter_type const& mu, size_type K = invalid_size_type_value, int indexmu = -1 ) const;
    boost::tuple<value_type,value_type> lbNoSCM( parameter_type const& mu, size_type K = invalid_size_type_value, int indexmu = -1 ) const;
    /**
     * Returns the lower bound of the coercive constant given a parameter \p
     * \f$\mu\f$
     *
     * \param mu \f$ \mu\f$ the parameter at which to evaluate the coercive constant
     * \param K the dimension of \f$Y_{\mathrm{UB}}\f$
     *
     *\return compute online the lower bounmd
     */
    boost::tuple<value_type,value_type> lb( parameter_ptrtype const& mu, size_type K = invalid_size_type_value, int indexmu = -1 ) const
    {
        return lb( *mu, K, indexmu );
    }


    /**
     * Returns the upper bound of the coercive constant given a parameter \p
     * \f$\mu\f$
     *
     * \param mu \f$ \mu\f$ the parameter at which to evaluate the coercive constant
     * \param K the dimension of \f$Y_{\mathrm{UB}}\f$
     *
     *\return compute online the lower bounmd
     */
    boost::tuple<value_type,value_type> ub( parameter_type const& mu, size_type K = invalid_size_type_value ) const;

    /**
     * Returns the upper bound of the coercive constant given a parameter \p
     * \f$\mu\f$
     *
     * \param mu \f$ \mu\f$ the parameter at which to evaluate the coercive constant
     * \param K the dimension of \f$Y_{\mathrm{UB}}\f$
     *
     *\return compute online the lower bounmd
     */
    boost::tuple<value_type,value_type> ub( parameter_ptrtype const& mu, size_type K = invalid_size_type_value ) const
    {
        return ub( *mu, K );
    }


    /**
     * Offline computation
     */
    std::vector<boost::tuple<double,double,double> > offline();
    std::vector<boost::tuple<double,double,double> > offlineSCM();
    std::vector<boost::tuple<double,double,double> > offlineNoSCM();

    /**
     * Online computation
     */
    bounds_type
    online() const
    {
        return bounds_type();
    }

    /**
     * \brief Retuns maximum value of the relative error
     * \param Xi fine sampling of the parameter space
     */
    relative_error_type maxRelativeError( size_type K ) const;

    /**
     * \brief compute the bounds for \f$y \in \mathbb{R}^{Q_a}\f$
     */
    void computeYBounds();

    /**
     * \param K size of C_K
     * \param mu parameter set
     */
    std::vector<double> run( parameter_type const& mu, int K );

    /**
     * run scm for a set of parameter
     * \param X - input parameter
     * \param N - number of input parameter
     * \param Y - outputs
     * \param P - number of outputs
     *
     * the output is the lower/upper bounds, hence P=2
     */
    void run( const double * X, unsigned long N, double * Y, unsigned long P );

    /**
     * save the CRB database
     */
    void saveDB();

    /**
     * load the CRB database
     */
    bool loadDB();

    bool rebuildDB();

    bool doScmForMassMatrix();
    //@}


    sampling_ptrtype c() const
    {
        return M_C;
    }

protected:

private:
private:
    friend class boost::serialization::access;

    template<class Archive>
    void save( Archive & ar, const unsigned int version ) const;

    template<class Archive>
    void load( Archive & ar, const unsigned int version ) ;

    BOOST_SERIALIZATION_SPLIT_MEMBER()

private:

    bool M_is_initialized;


    truth_model_ptrtype M_model;

    double M_tolerance;
    int M_iter_max;
    size_type M_Mplus;
    size_type M_Malpha;

    // parameter space
    parameterspace_ptrtype M_Dmu;

    // sampling of parameter space
    sampling_ptrtype M_Xi;

    // sampling of parameter space (subset of M_Xi later on)
    sampling_ptrtype M_C, M_C_complement;

    //! eigenvalues associated with the parameters used to build M_C
    std::map<size_type,value_type> M_C_eigenvalues;
    std::vector<double> M_eig;
    mutable std::map<size_type,std::map<size_type,value_type> > M_C_alpha_lb;

    //! y_set_type \f$\mathcal{U}_{\mathrm{UB}}\f$
    y_set_type M_Y_ub;

    //! bounds for y
    y_bounds_type M_y_bounds;
    std::vector< std::vector<double> > M_y_bounds_0;
    std::vector< std::vector<double> > M_y_bounds_1;

    po::variables_map M_vm;

    bool M_scm_for_mass_matrix;
    bool M_print_matrix;

    parameter_type M_mu_ref;
    bool M_use_scm;
};

template<typename TruthModelType>
std::vector<boost::tuple<double,double,double> >
CRBSCM<TruthModelType>::offline()
{
    if( M_use_scm )
        return this->offlineSCM();
    else
        return this->offlineNoSCM();
}

template<typename TruthModelType>
std::vector<boost::tuple<double,double,double> >
CRBSCM<TruthModelType>::offlineNoSCM()
{
    sparse_matrix_ptrtype inner_prod,sym,Matrix;
    M_mu_ref = M_model->refParameter();
    M_model->computeAffineDecomposition();
    M_model->countAffineDecompositionTerms();
    if ( M_scm_for_mass_matrix )
    {
        inner_prod = M_model->massMatrix();
        boost::tie( Matrix, boost::tuples::ignore, boost::tuples::ignore ) = M_model->update( M_mu_ref );
    }
    else
    {
        inner_prod = M_model->energyMatrix();
        boost::tie( boost::tuples::ignore, Matrix, boost::tuples::ignore ) = M_model->update( M_mu_ref );
    }
    sym = M_model->newMatrix();sym->close();
    Matrix->symmetricPart( sym );
    // solve  for eigenvalue problem at mu_ref
    SolverEigen<double>::eigenmodes_type modes;

    modes=
        eigs( _matrixA=sym,
              _matrixB=inner_prod,
              _solver=( EigenSolverType )M_vm["crb.scm.solvereigen-solver-type"].template as<int>(),
              _spectrum=SMALLEST_REAL,
              //_spectrum=LARGEST_MAGNITUDE,
              _transform=SINVERT,
              _ncv=M_vm["crb.scm.solvereigen-ncv"].template as<int>(),
              _nev=M_vm["crb.scm.solvereigen-nev"].template as<int>(),
              _tolerance=M_vm["crb.scm.solvereigen-tol"].template as<double>(),
              _maxit=M_vm["crb.scm.solvereigen-maxiter"].template as<int>()
              );
    double eigen_value = modes.begin()->second.template get<0>();
#if 0
    if( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
    {
        std::cout<<"-------------------------------------------"<<std::endl;
        std::cout<<"eigenvalue ( min ) for mu_ref : "<<eigen_value<<std::endl;
        std::cout<<"-------------------------------------------"<<std::endl;
    }
#endif
    LOG( INFO )<<"eigenvalue ( min ) for mu_ref : "<<eigen_value;


    if( boption(_name="crb.scm.check-eigenvector") )
    {
        auto eigen_vector = modes.begin()->second.template get<2>();
        checkEigenVectorEigenValue( sym, inner_prod, eigen_vector, eigen_value );
    }

    //store the eigen value in M_C_eigenvalues
    M_C_eigenvalues[0] = modes.begin()->second.template get<0>();

    saveDB();

    //only to have the same signature that offlineSCM
    std::vector<boost::tuple<double,double,double> > ckconv;
    return ckconv;
}

template<typename TruthModelType>
std::vector<boost::tuple<double,double,double> >
CRBSCM<TruthModelType>::offlineSCM()
{
    std::ofstream os_y( "y.m" );
    std::ofstream os_C( "C.m" );


    std::vector<boost::tuple<double,double,double> > ckconv;
    // do the affine decomposition
    M_model->computeAffineDecomposition();
    M_model->countAffineDecompositionTerms();
    // random sampling
    bool all_procs_have_same_sampling=true;
    M_Xi->randomize( M_vm["crb.scm.sampling-size"].template as<int>() , all_procs_have_same_sampling );
    //M_Xi->logEquidistribute( M_vm["crb.scm.sampling-size"].template as<int>() );
    M_C->setSuperSampling( M_Xi );
    parameter_type mu( M_Dmu );


    M_print_matrix = M_vm["crb.scm.print-matrix"].template as<bool>() ;

#if 0

    //Test coercicivty
    //std::cout << "Testing coercivity at sample points\n" ;
    for ( int i=0; i<M_Xi->size(); i++ )
    {
        double testalpha = ex( M_Xi->at( i ) );
        //std::cout << "mu[" << i+1 << "]\n" << M_Xi->at(i) << "\nalpha = " << testalpha << "\n" ;
    }

#endif
    double relative_error = 1e30;


    // compute the bounds for y in R^Q
    this->computeYBounds();

    // empty sets
    M_C->clear();
    M_Y_ub.clear();
    M_C_alpha_lb.clear();

    size_type index;

    bool use_predefined_C = boption(_name="crb.scm.use-predefined-C");
    int N_log_equi = this->vm()["crb.scm.use-logEquidistributed-C"].template as<int>() ;
    int N_equi = this->vm()["crb.scm.use-equidistributed-C"].template as<int>() ;
    std::vector<int> index_vector;

    if( N_log_equi > 0 || N_equi > 0 )
        use_predefined_C = true;

    if( use_predefined_C )
    {
        std::string file_name = ( boost::format("SamplingC") ).str();
        std::ifstream file ( file_name );
        if( ! file )
            throw std::logic_error( "[CRBSCM::offline] ERROR the file SamplingC doesn't exist so it's impossible to known which parameters you want to use to build the database" );
        else
        {
            M_C->clear();
            index_vector = M_C->closestSamplingFromFile(file_name);
            int sampling_size = index_vector.size();
            M_iter_max = sampling_size;
        }
        mu = M_C->at( 0 ); // first element
        index = index_vector[0];

        if( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
            std::cout<<"[CRBSCM::offline] read sampling C ( sampling size : "<<M_iter_max<<" )"<<std::endl;
    }
    else
    {
        // start with M_C = { arg min mu, mu \in Xi }
        boost::tie( mu, index ) = M_Xi->min();
        M_C->push_back( mu, index );
    }

    M_C_complement = M_C->complement();
    //std::cout << " -- start with mu = " << mu << "\n";
    //std::cout << " -- C size :  " << M_C->size() << "\n";
    //std::cout << " -- C complement size :  " << M_C_complement->size() << "\n";

    sparse_matrix_ptrtype B,symmMatrix,Matrix;
    std::vector<vector_ptrtype> F;

    // dimension of Y_UB and U_LB
    int K = 1;


    std::vector< std::vector<sparse_matrix_ptrtype> > Matrixq;

    if ( M_scm_for_mass_matrix )
        boost::tie( Matrixq, boost::tuples::ignore, boost::tuples::ignore ) = M_model->computeAffineDecomposition();
    else
        boost::tie( boost::tuples::ignore, Matrixq, boost::tuples::ignore ) = M_model->computeAffineDecomposition();

    int Qmax = nb_decomposition_terms_q();

    while ( relative_error > M_tolerance && K <= M_iter_max )
    {
        if( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
        {
            std::cout << "============================================================\n";
            std::cout << "K=" << K << "\n";
        }

        os_C << M_C->at( K-1 ) << "\n";

        // resize y_ub
        M_Y_ub.resize( K );
        M_Y_ub[K-1].resize( Qmax );
        for(int q=0; q<Qmax; q++)
            M_Y_ub[K-1][q].resize( mMax(q) );

        M_model->solve( mu );


        // for a given parameter \p mu assemble the left and right hand side
        if ( M_scm_for_mass_matrix )
        {
            B = M_model->massMatrix();
            boost::tie( Matrix, boost::tuples::ignore, F ) = M_model->update( mu );
        }
        else
        {
            B = M_model->energyMatrix();
            boost::tie( boost::tuples::ignore, Matrix, F ) = M_model->update( mu );
        }

        std::string mu_str;
        for(int i=0;i<mu.size();i++)
            mu_str= mu_str + (boost::format("_%1%") %mu[i]).str() ;

        std::string file_name;

        symmMatrix = M_model->newMatrix();symmMatrix->close();
        Matrix->symmetricPart( symmMatrix );

        if( M_print_matrix && ( Environment::worldComm().globalSize() == 1 ) )
        {
            if( M_scm_for_mass_matrix)
            {
                file_name = "offline_Matrix_M" + mu_str +".m";
                Matrix->printMatlab( file_name );
                file_name = "offline_symmMatrix_M" + mu_str +".m";
                symmMatrix->printMatlab( file_name );
            }
            else
            {
                file_name = "offline_Matrix_A" + mu_str +".m";
                Matrix->printMatlab( file_name );
                file_name = "offline_symmMatrix_A" + mu_str +".m";
                symmMatrix->printMatlab( file_name );
            }
            file_name = "offline_B" + mu_str +".m";
            B->printMatlab( file_name );
        }

        //
        // Build Y_UB
        //
        vector_ptrtype eigenvector;
        // solve  for eigenvalue problem at \p mu
        SolverEigen<double>::eigenmodes_type modes;

        modes=
            eigs( _matrixA=symmMatrix,
                  _matrixB=B,
                  _solver=( EigenSolverType )M_vm["crb.scm.solvereigen-solver-type"].template as<int>(),
                  _spectrum=SMALLEST_REAL,
                  //_spectrum=LARGEST_MAGNITUDE,
                  _transform=SINVERT,
                  _ncv=M_vm["crb.scm.solvereigen-ncv"].template as<int>(),
                  _nev=M_vm["crb.scm.solvereigen-nev"].template as<int>(),
                  _tolerance=M_vm["crb.scm.solvereigen-tol"].template as<double>(),
                  _maxit=M_vm["crb.scm.solvereigen-maxiter"].template as<int>()
                );

        if ( modes.empty()  )
        {
            LOG(INFO) << "eigs failed to converge\n";
            return ckconv;
        }

        LOG( INFO ) << "[fe eig] mu=" << std::setprecision( 4 ) << mu ;
        LOG( INFO ) << "[fe eig] eigmin : " << std::setprecision( 16 ) << modes.begin()->second.template get<0>() ;
        LOG( INFO ) << "[fe eig] ndof:" << M_model->functionSpace()->nDof() ;

        // extract the eigenvector associated with the smallest eigenvalue
        eigenvector = modes.begin()->second.template get<2>();

        // store the  eigenvalue associated with \p mu
        M_C_eigenvalues[index] = modes.begin()->second.template get<0>();
        typedef std::pair<size_type,value_type> key_t;

        if( boption(_name="crb.scm.check-eigenvector") )
        {
            auto eigenvalue = modes.begin()->second.template get<0>();
            checkEigenVectorEigenValue( symmMatrix, B, eigenvector, eigenvalue );
        }

        //BOOST_FOREACH( key_t eig, M_C_eigenvalues )
	    //std::cout << "[fe eig] stored/map eig=" << eig.second <<" ( " << eig.first << " ) " << "\n";

        M_eig.push_back( M_C_eigenvalues[index] );
        //BOOST_FOREACH( value_type eig, M_eig )
	    //std::cout << "[fe eig] stored/vec eig=" << eig << "\n";


        /*
         * now apply eigenvector to the Aq to compute
         * y( eigenvector ) = ( a_q( eigenvector, eigenvector )/ ||eigenvector||^2
         */
        for ( size_type q = 0; q < nb_decomposition_terms_q() ; ++q )
        {
            for ( size_type m_eim = 0; m_eim < mMax(q) ; ++m_eim )
            {
                value_type aqmw = Matrixq[q][m_eim]->energy( eigenvector, eigenvector );
                value_type bw = B->energy( eigenvector, eigenvector );
                LOG( INFO ) << "[scm_offline] q=" << q << " aqmw = " << aqmw << ", bw = " << bw ;
                M_Y_ub[K-1][ q ](m_eim) = aqmw/bw;
            }

            //checkC( K );
            //std::cout << "[scm_offline] M_Y_ub[" << K-1 << "]=" << M_Y_ub[K-1] << "\n";

            // save Y
            os_y << M_Y_ub[K-1][q] << "\n";
        }//q


	if( M_scm_for_mass_matrix )
	    LOG( INFO ) <<"scm is done for mass matrix";
	else
	    LOG( INFO )<<"scm is done for a( . , . ; mu )";

        double minerr, meanerr;
        if( use_predefined_C )
        {
            relative_error = M_tolerance+10;
            minerr = relative_error;
            meanerr = relative_error;
            if( K < M_iter_max )
            {
                mu = M_C->at( K );
                index = index_vector[ K ];
            }
        }
        else
            boost::tie( relative_error, mu, index, minerr, meanerr ) = maxRelativeError( K );
#if 0
        std::cout << " -- max relative error = " << relative_error
                  << " at mu = " << mu
                  << " at index = " << index << "\n";
#endif
        //ofs << K << " " << std::setprecision(16) << relative_error << std::endl;
        ckconv.push_back( boost::make_tuple( relative_error, minerr, meanerr ) );

        // could be that the max relative error is smaller than the tolerance if
        // the coercivity constant is independant of the parameter set
        if ( relative_error > M_tolerance && K < M_iter_max  && ! use_predefined_C )
        {
            if( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
                std::cout << " -- inserting mu - index : "<<index<<" -  in C (" << M_C->size() << ")\n";
            M_C->push_back( mu, index );

            //for ( size_type _i =0; _i < M_C->size(); ++_i )
	        //std::cout << " -- mu [" << _i << "]=" << M_C->at( _i ) << std::endl;

            M_C_complement = M_C->complement();

            //for ( size_type _i =0; _i < M_C_complement->size(); ++_i )
            //std::cout << " -- mu complement [" << _i << "]=" << M_C_complement->at( _i ) << std::endl;
        }

        ++K;
        if( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
            std::cout << "============================================================\n";
    }

    //before call saveDB we have to split the vector of tuple M_y_bounds
    M_y_bounds_0.resize( Qmax );
    M_y_bounds_1.resize( Qmax );
    for ( int q=0; q<Qmax; q++ )
    {
        for(int m=0; m<mMax(q); m++)
        {
            M_y_bounds_0[q].push_back( M_y_bounds[q][m].template get<0>() );
            M_y_bounds_1[q].push_back( M_y_bounds[q][m].template get<1>() );
        }
    }

    saveDB();
    return ckconv;
}

template<typename TruthModelType>
void
CRBSCM<TruthModelType>::checkEigenVectorEigenValue( sparse_matrix_ptrtype const& A, sparse_matrix_ptrtype const& B, vector_ptrtype const& eigenvector, double eigenvalue ) const
{
    //auto backend = backend_type::build( BACKEND_PETSC, Environment::worldComm() );
    auto Xh = M_model->functionSpace() ;
    auto Aw =  backend()->newVector( Xh );
    auto Bw =  backend()->newVector( Xh );
    A->multVector( eigenvector, Aw );
    B->multVector( eigenvector, Bw );
    //we should have Aw = eigen_value Bw
    Bw->scale( eigenvalue );
    double energyAwAw = A->energy( Aw , Aw );
    double energyAwBw = A->energy( Aw , Bw );
    double energyBwBw = A->energy( Bw , Bw );
    double tol = doption(_name="crb.scm.check-eigenvector-tol");
    CHECK( math::abs(energyAwAw - energyAwBw) <  tol )<<"eigen vector and/or eigen value not satisfy generalized eigenvalue problem : math::abs(energyAwAw - energyAwBw) = "<<math::abs(energyAwAw - energyAwBw)<<std::endl;
    CHECK( math::abs(energyAwAw - energyAwBw) <  tol )<<"eigen vector and/or eigen value not satisfy generalized eigenvalue problem : math::abs(energyAwAw - energyAwBw) = "<<math::abs(energyAwAw - energyAwBw)<<std::endl;
}

template<typename TruthModelType>
void
CRBSCM<TruthModelType>::checkC( size_type K ) const
{
    y_type err( M_Xi->size() );

    // check that for each mu in C_K the lb and ub are very close
    for ( int k = 0; k < K; ++k )
    {
        std::cout << "[checkC] k= " << k << " K=" << K << " index in superSampling: " << M_C->indexInSuperSampling( k ) <<"\n";
        parameter_type const& mu = M_C->at( k );
        size_type index = M_C->indexInSuperSampling( k );
        double _lb = lb( mu, K, index  ).template get<0>();
        double _lblb = lb( mu, std::max( K-1,size_type( 0 ) ), index ).template get<0>();

        if ( _lblb - _lb > 1e-10 )
        {
            LOG(INFO) << "[lberror] the lower bound is decreasing\n"
                  << "[lberror] _lblb = " << _lblb << "\n"
                  << "[lberror] _lb = " << _lb << "\n";
        }

        double _ub = ub( mu, K ).template get<0>();
        double _ex = M_C_eigenvalues.find( index )->second;
        std::cout << "_lb = " << _lb << " _lblb=" << _lblb << " ub = " << _ub << " _ex = " << _ex << "\n";

        err( k ) = 1. - _lb/_ub;

        if ( std::abs( err( k ) ) > 1e-6 )
        {
            std::ostringstream ostr;
            ostr << "error is not small as it should be " << k <<  " " << err( k ) << "\n";
            ostr << "[checkC] relative error : " << mu << "\n"
                 << "[checkC] lb=" << std::setprecision( 16 ) << _lb << "\n"
                 << "[checkC] ub=" << std::setprecision( 16 ) << _ub << "\n"
                 << "[checkC] relerr(" << k << ")=" << std::setprecision( 16 ) << err( k ) << "\n"
                 << "[checkC] exact(" << k << ")=" << std::setprecision( 16 ) << _ex << "\n";
            throw std::logic_error( ostr.str() );
        }

#if 1
        std::cout << "[checkC] relative error : " << mu << "\n"
                  << "[checkC] lb=" << std::setprecision( 16 ) << _lb << "\n"
                  << "[checkC] ub=" << std::setprecision( 16 ) << _ub << "\n"
                  << "[checkC] relerr(" << k << ")=" << std::setprecision( 16 ) << err( k ) << "\n"
                  << "[checkC] exerr(" << k << ")=" << std::setprecision( 16 ) << _ex << "\n";
#endif
    }
}


template<typename TruthModelType>
boost::tuple<typename CRBSCM<TruthModelType>::value_type,
      typename CRBSCM<TruthModelType>::value_type>
      CRBSCM<TruthModelType>::ex( parameter_type const& mu ) const
{
    boost::timer ti;
    sparse_matrix_ptrtype Matrix,symmMatrix;
    sparse_matrix_ptrtype M ;
    std::vector<vector_ptrtype> F;

    if ( M_scm_for_mass_matrix )
    {
        M = M_model->massMatrix();
        boost::tie( Matrix, boost::tuples::ignore, F ) = M_model->update( mu );
    }
    else
    {
        M = M_model->energyMatrix();
        boost::tie( boost::tuples::ignore, Matrix, F ) = M_model->update( mu );
    }

    symmMatrix = M_model->newMatrix();
    symmMatrix->close();
    Matrix->symmetricPart( symmMatrix );

    SolverEigen<double>::eigenmodes_type modesmin=
        eigs( _matrixA=Matrix,
              _matrixB=M,
              _solver=( EigenSolverType )M_vm["crb.scm.solvereigen-solver-type"].template as<int>(),
              //_spectrum=LARGEST_MAGNITUDE,
              //_spectrum=LARGEST_REAL,
              _spectrum=SMALLEST_REAL,
              //_spectrum=SMALLEST_MAGNITUDE,
              _transform=SINVERT,
              _ncv=M_vm["crb.scm.solvereigen-ncv"].template as<int>(),
              _nev=M_vm["crb.scm.solvereigen-nev"].template as<int>(),
              _tolerance=M_vm["crb.scm.solvereigen-tol"].template as<double>(),
              _maxit=M_vm["crb.scm.solvereigen-maxiter"].template as<int>()
            );

    if ( modesmin.empty() )
    {
        std::cout << "no eigenmode converged: increase --solvereigen-ncv\n";
        return 0.;
    }

    double eigmin = modesmin.begin()->second.template get<0>();
    //std::cout << std::setprecision(4) << mu << " eigmin = "
    //          << std::setprecision(16) << eigmin << "\n";

    return boost::make_tuple( eigmin, ti.elapsed() );
#if 0
    SolverEigen<double>::eigenmodes_type modesmax=
        eigs( _matrixA=Matrix,
              _matrixB=M,
              _solver=( EigenSolverType )M_vm["crb.scm.solvereigen-solver-type"].as<int>(),
              _spectrum=LARGEST_MAGNITUDE,
              _ncv=M_vm["crb.scm.solvereigen-ncv"].as<int>(),
              _nev=M_vm["crb.scm.solvereigen-nev"].as<int>(),
              _tolerance=M_vm["crb.scm.solvereigen-tol"].as<double>(),
              _maxit=M_vm["crb.scm.solvereigen-maxiter"].as<int>()
            );

    if ( modesmax.empty() )
    {
        std::cout << "no eigenmode converged: increase --solvereigen-ncv\n";
        return;
    }

    double eigmax = modesmax.rbegin()->second.template get<0>();

    //std::cout << "[crbmodel] " << std::setprecision(4) << mu << " "
    //          << std::setprecision(16) << eigmin << " " << eigmax
    //          << " " << Xh->nDof() << "\n";

#endif

}



template<typename TruthModelType>
boost::tuple<typename CRBSCM<TruthModelType>::value_type, double>
CRBSCM<TruthModelType>::lb( parameter_type const& mu ,size_type K ,int indexmu ) const
{
    if( M_use_scm )
        return this->lbSCM( mu,  K , indexmu );
    else
        return this->lbNoSCM( mu,  K , indexmu  );
}

template<typename TruthModelType>
boost::tuple<typename CRBSCM<TruthModelType>::value_type, double>
CRBSCM<TruthModelType>::lbNoSCM( parameter_type const& mu ,size_type K ,int indexmu ) const
{

    boost::mpi::timer ti;

    vector_type vec_min_coeff;
    auto all_beta = M_model->computeBetaQm( mu );
    auto all_beta_ref = M_model->computeBetaQm( M_mu_ref );
    beta_vector_type beta_mu;
    beta_vector_type beta_mu_ref;
    if ( M_scm_for_mass_matrix )
    {
        beta_mu = all_beta.template get<0>();
        beta_mu_ref = all_beta_ref.template get<0>();
    }
    else
    {
        beta_mu = all_beta.template get<1>();
        beta_mu_ref = all_beta_ref.template get<1>();
    }

    int Q = this->nb_decomposition_terms_q();
    vec_min_coeff.resize(Q);
    for(int q=0; q<Q; q++)
    {
        //for each q, search the min coeff
        vector_type vec_local_min_coeff;
        int M = beta_mu[q].size();
        vec_local_min_coeff.resize( M );
        for(int m=0; m<M; m++)
        {
            double __beta = beta_mu[q][m]/beta_mu_ref[q][m];
            vec_local_min_coeff(m) = __beta ;
        }
        double local_min = vec_local_min_coeff.minCoeff();
        vec_min_coeff(q) =  local_min ;
    }
    double min = vec_min_coeff.minCoeff();
    double lower_bound = min * M_C_eigenvalues.find(0)->second;

    return boost::make_tuple( lower_bound, ti.elapsed() );
}

template<typename TruthModelType>
boost::tuple<typename CRBSCM<TruthModelType>::value_type, double>
CRBSCM<TruthModelType>::lbSCM( parameter_type const& mu ,size_type K ,int indexmu ) const
{

#if defined(FEELPP_HAS_GLPK_H)
    if ( K == invalid_size_type_value ) K = this->KMax();

    if ( K > this->KMax() ) K = this->KMax();

    boost::mpi::timer ti;

    // value if K==0
    if ( K <= 0 ) return 0.0;

    if ( indexmu >= 0 && ( M_C_alpha_lb[indexmu].find( K ) !=  M_C_alpha_lb[indexmu].end() ) )
        return  M_C_alpha_lb[indexmu][K] ;

    //if ( K == std::max(size_type(0),M_C->size()-(size_type)M_vm["crb.scm.level"].template as<int>() ) ) return 0.0;
    //int level = M_vm["crb.scm.level"].template as<int>();

    //std::cout << "[CRBSCM::lb] Alphalb size " << M_C_alpha_lb.size() << "\n";

    beta_vector_type beta_qm;
    glp_prob *lp;


    lp = glp_create_prob();
    glp_set_prob_name( lp, ( boost::format( "scm_%1%" ) % K ).str().c_str() );
    glp_set_obj_dir( lp, GLP_MIN );

    int Malpha = std::min( M_Malpha,std::min( K,M_Xi->size() ) );

    int Mplus = std::min( M_Mplus,M_Xi->size()-K );


    if ( M_vm["crb.scm.strategy"].template as<int>()==2 )
        Mplus = std::min( M_Mplus,std::min( K, M_Xi->size()-K ) );

    // we have exactely Qa*(M+ + Malpha) entries in the matrix


    int nnz = nb_decomposition_terms_qm()*( Malpha+Mplus );

    int ia[1+10000], ja[1+10000];
    double ar[1+10000];
    int nnz_index = 1;

    // set the auxiliary variables: we have first Malpha of them from C_K and
    // Mplus from its complement
    glp_add_rows( lp,Malpha+Mplus );
    // search the the M_Malpha closest points in C_K, M_Malpha must be < K
    // index_vector will contain index of neighbors
    std::vector<int> index_vector;
    sampling_ptrtype C_neighbors =  M_C->searchNearestNeighbors( mu, Malpha , index_vector);

    //std::cout << "[CRBSCM::lb] C_neighbors size = " << C_neighbors->size() << "\n";

    //std::cout << "[CRBSCM::lb] add rows associated with C_K\n";
    // first the constraints associated with C_K: make sure that Malpha is not
    // greater by taking the min of M_Malpha and K

    for ( int m=0; m < Malpha; ++m )
    {

        //std::cout << "[CRBSCM::lb] add row " << m << " from C_K\n";
        parameter_type mup = C_neighbors->at( m );

        // update the theta_q associated with mup
        if ( M_scm_for_mass_matrix )
        {
            boost::tie( beta_qm , boost::tuples::ignore , boost::tuples::ignore ) = M_model->computeBetaQm( mup );
        }
        else
        {
            boost::tie( boost::tuples::ignore, beta_qm, boost::tuples::ignore ) = M_model->computeBetaQm( mup );
        }

        //std::cout << "[CRBSCM::lb] thetaq = " << theta_q << "\n";

        //std::cout << "[CRBSCM::lb] row name : " << (boost::format( "c_%1%_%2%" ) % K % m).str() << "\n";
        glp_set_row_name( lp, m+1, ( boost::format( "c_%1%_%2%" ) % K % m ).str().c_str() );

        //std::cout << "[CRBSCM::lb] row bounds\n";
        //std::cout << "[CRBSCM::lb] index in super sampling: " << C_neighbors->indexInSuperSampling( m ) << "\n";
        //std::cout << "[CRBSCM::lb] eig value: " << M_C_eigenvalues.find( C_neighbors->indexInSuperSampling( m ) )->second << "\n";
        glp_set_row_bnds( lp,
                          m+1,
                          GLP_LO,
                          M_C_eigenvalues.find( C_neighbors->indexInSuperSampling( m ) )->second,
                          0.0 );


        //std::cout << "[CRBSCM::lb] constraints matrix\n";
        int count=1;
        for ( int q = 0; q < nb_decomposition_terms_q(); ++q )
        {
            for ( int m_eim = 0; m_eim < mMax( q ); ++m_eim, ++nnz_index, ++count )
            {
                //std::cout << "[CRBSCM::lb] constraints matrix q,m = " << q << ","<< m <<"\n";
                ia[nnz_index]=m+1;
                ja[nnz_index]=count;
                ar[nnz_index]=beta_qm[q][m_eim];
            }
        }

    }
    //std::cout << "[CRBSCM::lb] add rows associated with C_K done. nnz=" << nnz_index << "\n";

    // search the the Mplus closest points in Xi\C_K
    sampling_ptrtype Xi_C_neighbors =  M_C_complement->searchNearestNeighbors( mu, Mplus , index_vector);

    //std::cout << "[CRBSCM::lb] C_complement size = " << Xi_C_neighbors->size() << "\n";

    //std::cout << "[CRBSCM::lb] add rows associated with Xi\\C_K\n";
    //std::cout << "[CRBSCM::lb] Mplus =" << Mplus << " , neighbors= " << Xi_C_neighbors->size() << "\n";
    // second the monotonicity constraints associated with the complement of C_K
    for ( int m=0; m < Mplus; ++m )
        //for( int m=0;m < Xi_C_neighbors->size(); ++m )
    {

        //std::cout << "[CRBSCM::lb] add row " << m << " from Xi\\C_K\n";
        parameter_type mup = Xi_C_neighbors->at( m );

        double _lb;

#if 1
        //std::cout << "[CRBSCM::lb] Entering find alphamap\n" << Xi_C_neighbors->indexInSuperSampling( m ) << std::endl ;
        //std::cout << "[CRBSCM::lb] Size of alphalb " << M_C_alpha_lb.find( Xi_C_neighbors->indexInSuperSampling( m ) )->second.size() << "\n" ;

        if ( M_C_alpha_lb[Xi_C_neighbors->indexInSuperSampling( m ) ].find( K-1 ) !=  M_C_alpha_lb[Xi_C_neighbors->indexInSuperSampling( m ) ].end() )
        {
            //std::cout << "[CRBSCM::lb] lb déjà calculée\n" ;
        }
        else
        {
            //std::cout << "[CRBSCM::lb] Nouvelle lb\n" ;
            M_C_alpha_lb[  Xi_C_neighbors->indexInSuperSampling( m )][K-1] = lb( mup, K-1,  Xi_C_neighbors->indexInSuperSampling( m ) ).template get<0>();
        }

        _lb = M_C_alpha_lb[ Xi_C_neighbors->indexInSuperSampling( m )][K-1];
#endif

        // update the theta_q associated with mup
        if ( M_scm_for_mass_matrix )
        {
            boost::tie( beta_qm, boost::tuples::ignore, boost::tuples::ignore ) = M_model->computeBetaQm( mup );
        }

        else
        {
            boost::tie( boost::tuples::ignore, beta_qm, boost::tuples::ignore ) = M_model->computeBetaQm( mup );
        }

        glp_set_row_name( lp, Malpha+m+1, ( boost::format( "xi_c_%1%_%2%" ) % K % m ).str().c_str() );

        switch ( M_vm["crb.scm.strategy"].template as<int>() )
        {
            // Patera
        case 0:
        {
            //std::cout << "Patera\n" ;
            glp_set_row_bnds( lp, Malpha+m+1, GLP_LO, 0.0, 0.0 );
        }
        break;
        // Maday

        case 1:

            // Prud'homme
        case 2:
        {
            //std::cout << "Prud'homme\n" ;
            //    std::cout << "Lb = " << _lb << std::endl ;
            glp_set_row_bnds( lp, Malpha+m+1, GLP_LO, _lb, 0.0 );
        }
        break;
        }

        int count=1;
        int nb_already_done = nb_decomposition_terms_qm();
        for ( int q = 0; q < nb_decomposition_terms_q(); ++q )
        {
            for ( int m_eim = 0; m_eim < this->mMax(q); ++m_eim, ++nnz_index,++count )
            {
                //std::cout << "[CRBSCM::lb] constraints matrix q = " << q <<" , "<<m<< "\n";
                ia[nnz_index]=Malpha+m+1;
                ja[nnz_index]=count;
                ar[nnz_index]=beta_qm[q][m_eim];
            }
        }
    }

    //std::cout << "[CRBSCM::lb] add rows associated with C_K complement done. nnz=" << nnz_index << "\n";

    // set the structural variables, we have M_model->Qa() of them
    if ( M_scm_for_mass_matrix )
    {
        boost::tie( beta_qm,boost::tuples::ignore, boost::tuples::ignore ) = M_model->computeBetaQm( mu );
    }
    else
    {
        boost::tie( boost::tuples::ignore, beta_qm, boost::tuples::ignore ) = M_model->computeBetaQm( mu );
    }

    //nb_columns
#if 0
    int nb_col = 0;
    for(int q=0; q<nb_decomposition_terms_q(); q++)
    {
        for(int m_eim=0; m_eim<this->mMax(q); m_eim++)
            nb_col++;
    }
#endif
    // glp_add_cols( lp, nb_col );
    glp_add_cols( lp, nb_decomposition_terms_qm() );
    int count=0;
    for ( int q = 0; q < nb_decomposition_terms_q(); ++q )
    {
        for(int m_eim=0; m_eim<this->mMax(q); ++m_eim,++count)
        {
            glp_set_col_name( lp, count+1, ( boost::format( "y_%1%-%2%" ) % q %m_eim ).str().c_str() );
            glp_set_col_bnds( lp, count+1, GLP_DB,
                              M_y_bounds[q][m_eim].template get<0>(),
                              M_y_bounds[q][m_eim].template get<1>() );
            glp_set_obj_coef( lp, count+1, beta_qm[q][m_eim] );
        }
    }

    // now we need to build the matrix
    glp_load_matrix( lp, nnz, ia, ja, ar );
    glp_smcp parm;
    glp_init_smcp( &parm );
    parm.msg_lev = GLP_MSG_ERR;

    // use the simplex method and solve the LP
    glp_simplex( lp, &parm );

    // retrieve the minimum
    double Jobj = glp_get_obj_val( lp );
    //std::cout << "Jobj = " << Jobj << "\n";
#if 0

    for ( size_type q = 0; q < M_model->Qa(); ++q )
    {
        double y = glp_get_col_prim( lp,q+1 );
        //std::cout << "y" << q << " = " << y << "\n";
    }

#endif
    glp_delete_prob( lp );

    if ( indexmu >= 0 )
    {

        M_C_alpha_lb[ indexmu ][K] = Jobj;

    }
#else // FEELPP_HAS_GLPK_H

    double Jobj = 0;
    boost::timer ti;
#endif /* FEELPP_HAS_GLPK_H */

    return boost::make_tuple( Jobj, ti.elapsed() );

}


template<typename TruthModelType>
boost::tuple<typename CRBSCM<TruthModelType>::value_type, double>
CRBSCM<TruthModelType>::ub( parameter_type const& mu ,size_type K ) const
{
    if ( K == invalid_size_type_value ) K = this->KMax();

    if ( K > this->KMax() ) K = this->KMax();

    boost::timer ti;
    beta_vector_type beta_qm;

    if ( M_scm_for_mass_matrix )
    {
        boost::tie( beta_qm ,boost::tuples::ignore, boost::tuples::ignore ) = M_model->computeBetaQm( mu );
    }

    else
    {
        boost::tie( boost::tuples::ignore, beta_qm, boost::tuples::ignore ) = M_model->computeBetaQm( mu );
    }

    //std::cout << "[CRBSCM<TruthModelType>::ub] theta_q = " << theta_q << "\n";
    y_type y( K );

    //for ( size_type k = 0; k < K; ++k )
    //{
    //    y( k ) = ( beta_qm.array()*M_Y_ub[k].array() ).sum();
    //}
    for ( size_type k = 0; k < K; ++k )
    {
        y( k ) = 0;
        for(int q=0; q<nb_decomposition_terms_q(); q++)
        {
            //y( k ) = ( beta_qm[q].array()*M_Y_ub[k][q].array() ).sum();
            for(int m=0; m<mMax(q); m++)
                y( k ) +=  beta_qm[q][m]*M_Y_ub[k][q][m];
        }
    }

    //std::cout << "[CRBSCM<TruthModelType>::ub] y = " << y << "\n";
    return boost::make_tuple( y.minCoeff(), ti.elapsed() );

}

template<typename TruthModelType>
typename CRBSCM<TruthModelType>::relative_error_type
CRBSCM<TruthModelType>::maxRelativeError( size_type K ) const
{
    //std::cout << "==================================================\n";
#if 0
    y_type err( M_Xi->size() );

    for ( size_type k = 0; k < M_Xi->size(); ++k )
    {
        parameter_type const& mu = M_Xi->at( k );
#else
    y_type err( M_C_complement->size() );

    for ( size_type k = 0; k < M_C_complement->size(); ++k )
    {
        parameter_type const& mu = M_C_complement->at( k );
#endif

        //std::cout << "[maxRelativeError] Calcul de lb pour mu[" <<  k << "]\n" ;
        double _lb = lb( mu, K, k ).template get<0>();
        //std::cout << "[maxRelativeError] Calcul de lblb pour mu[" << k << "]\n" ;
        double _lblb = lb( mu, std::max( K-1,size_type( 0 ) ), k ).template get<0>();

        if ( _lblb - _lb > 1e-10 )
        {
            //LOG(INFO) << "[lberror] the lower bound is decreasing\n"
            //	<< "[lberror] _lblb = " << _lblb << "\n"
            //	<< "[lberror] _lb = " << _lb << "\n";
        }

        double _ub = ub( mu, K ).template get<0>();

        err( k ) = 1. - _lb/_ub;
#if 0
        //std::cout << "[maxRelativeError] k = " << k << "\n";
        //std::cout << "[maxRelativeError] parameter : " << mu << "\n";
        //std::cout << "[maxRelativeError] lb = " << std::setprecision(16) << _lb << "\n";
        //std::cout << "[maxRelativeError] ub = " << std::setprecision(16) << _ub << "\n";
        //std::cout << "[maxRelativeError] rel err(" << k << ")=" << std::setprecision(16) << err(k) << "\n";
#endif
    }


    Eigen::MatrixXf::Index index;
    double maxerr = err.array().abs().maxCoeff( &index );
    //std::cout << "[maxRelativeError] K=" << K << " max Error = " << maxerr << " at index = " << index << "\n";
    Eigen::MatrixXf::Index indexmin;
    double minerr = err.array().abs().minCoeff( &indexmin );
    //std::cout << "[maxRelativeError] K=" << K << " min Error = " << minerr << " at index = " << indexmin << "\n";
    double meanerr = err.array().abs().sum()/err.size();
    //std::cout << "[maxRelativeError] K=" << K << " mean Error = " << meanerr << "\n";
    //std::cout << "==================================================\n";
    //return boost::make_tuple( maxerr, M_Xi->at( index ), index, minerr, meanerr );

#if 0
    parameter_type const& mu = M_C_complement->at( index );
    double _lb = lb( mu, K, index ).template get<0>();
    double _ub = ub( mu, K ).template get<0>();
    std::cout<<"_lb = "<<_lb<<" and _ub = "<<_ub<<std::endl;
#endif

    return boost::make_tuple( maxerr, M_C_complement->at( index ), M_C_complement->indexInSuperSampling( index ), minerr, meanerr );
}

template<typename TruthModelType>
void
CRBSCM<TruthModelType>::computeYBounds()
{
    //std::cout << "************************************************************\n";
    LOG(INFO) << "[CRBSCM<TruthModelType>::computeYBounds()] start...\n";
    int Qmax = nb_decomposition_terms_q();
    M_y_bounds.resize(Qmax);
    sparse_matrix_ptrtype Matrix, symmMatrix=M_model->newMatrix(), B;


    if ( M_scm_for_mass_matrix )
        B=M_model->massMatrix();
    else
        B=M_model->energyMatrix();

    B->close();


    // solve 2 * Q_a eigenproblems
    for ( int q = 0; q < nb_decomposition_terms_q() ; ++q )
    {

        for ( int m = 0; m < mMax(q) ; ++m )
        {
            //std::cout << "================================================================================\n";
            //std::cout << "[ComputeYBounds] = q = " << q << " / " << M_model->Qa() << "\n";
            LOG(INFO) << "[CRBSCM<TruthModelType>::computeYBounds()] q = " << q << "/" << M_model->Qa() << "\n";
            // for a given parameter \p mu assemble the left and right hand side
            std::ostringstream os;

            if ( M_scm_for_mass_matrix )
                Matrix = M_model->Mqm( q , m );
            else
                Matrix = M_model->Aqm( q , m );

            Matrix->close();
            Matrix->symmetricPart( symmMatrix );
            if( Environment::worldComm().globalSize() == 1 )
            {
                os << "yb_Matrix" << q << " - "<< m << ".m";
                Matrix->printMatlab( os.str() );
                os.str( "" );
                os << "yb_symmMatrix" << q << " - "<< m <<".m";
                symmMatrix->printMatlab( os.str() );
            }
            if ( symmMatrix->l1Norm()==0.0 )
            {
                std::cout << "matrix is null\n" ;
                M_y_bounds[q].push_back( boost::make_tuple( 0.0, 1e-10 ) );
            }

            else
            {

#if 0
                // solve  for eigenvalue problem at \p mu
                boost::tie( nconv, eigenvalue_lb, boost::tuples::ignore, boost::tuples::ignore ) =
                    eigs( _matrixA=symmA,
                          _matrixB=B,
                          _spectrum=SMALLEST_MAGNITUDE,
                          _ncv=15 );

                //std::cout << " -- lower bounds q=" << q
                //          << " nconv=" << nconv
                //          << " eigenvalue_lb = " << eigenvalue_lb << "\n";
                double eigmin=eigenvalue_lb;

                boost::tie( nconv, eigenvalue_ub, boost::tuples::ignore, boost::tuples::ignore ) =
                    eigs( _matrixA=symmA,
                          _matrixB=B,
                          _spectrum=LARGEST_MAGNITUDE );

                //std::cout << " -- upper bounds q=" << q
                //          << " nconv=" << nconv
                //          << " eigenvalue_ub = " << eigenvalue_ub << "\n";
                double eigmax=eigenvalue_ub;
#else

                SolverEigen<double>::eigenmodes_type modes;
#if 1
                // solve  for eigenvalue problem at \p mu
                modes =
                    eigs( _matrixA=symmMatrix,
                          _matrixB=B,
                          //_problem=(EigenProblemType)PGNHEP,
                          _problem=( EigenProblemType )GHEP,
                          _solver=( EigenSolverType )M_vm["crb.scm.solvereigen-solver-type"].template as<int>(),
                          //_spectrum=SMALLEST_REAL,
                          _spectrum=SMALLEST_MAGNITUDE,
                          _ncv=M_vm["crb.scm.solvereigen-ncv"].template as<int>(),
                          _nev=M_vm["crb.scm.solvereigen-nev"].template as<int>(),
                          _tolerance=M_vm["crb.scm.solvereigen-tol"].template as<double>(),
                          _maxit=M_vm["crb.scm.solvereigen-maxiter"].template as<int>()
                          );

#endif

                if ( modes.empty() )
                {
                    LOG(INFO) << "[Computeybounds] eigmin did not converge for q=" << q << " (set to 0)\n";
                    if( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
                        std::cout << "[Computeybounds] eigmin did not converge for q=" << q << " (set to 0)"<<std::endl;
                }

                double eigmin = modes.empty()?0:modes.begin()->second.template get<0>();
#if 1
                modes=
                    eigs( _matrixA=symmMatrix,
                          _matrixB=B,
                          //_problem=(EigenProblemType)PGNHEP,
                          _problem=( EigenProblemType )GHEP,
                          _solver=( EigenSolverType )M_vm["crb.scm.solvereigen-solver-type"].template as<int>(),
                          _spectrum=LARGEST_REAL,
                          //_spectrum=LARGEST_MAGNITUDE,
                          _ncv=M_vm["crb.scm.solvereigen-ncv"].template as<int>(),
                          //_ncv=20,
                          _nev=M_vm["crb.scm.solvereigen-nev"].template as<int>(),
                          _tolerance=M_vm["crb.scm.solvereigen-tol"].template as<double>(),
                          //_tolerance=1e-7,
                          _maxit=M_vm["crb.scm.solvereigen-maxiter"].template as<int>()
                          );

#endif

                if ( modes.empty() )
                {
                    LOG(INFO) << "[Computeybounds] eigmax did not converge for q=" << q << " (set to 0)\n";
                    if( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
                        std::cout << "[Computeybounds] eigmax did not converge for q=" << q << " (set to 0)"<<std::endl;
                }

                double eigmax = modes.empty()?0:modes.rbegin()->second.template get<0>();
                if( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
                    std::cout<<"[computeYBounds] bounds for (q,m) = ("<<q<<","<<m<<") [ "<<eigmin<<" ; "<<eigmax<<"]"<<std::endl;
                LOG(INFO)<<"[computeYBounds] bounds for (q,m) = ("<<q<<","<<m<<") [ "<<eigmin<<" ; "<<eigmax<<"]\n";
                //std::cout << "[Computeybounds] q= " << q << " eigmin=" << std::setprecision(16) << eigmin << " eigmax=" << std::setprecision(16) << eigmax << "\n";
                //std::cout << std::setprecision(16) << eigmin << " " << eigmax << "\n";

#endif

                if ( eigmin==0.0 && eigmax==0.0 ) throw std::logic_error( "eigs null\n" );
                M_y_bounds[q].push_back( boost::make_tuple( eigmin, eigmax ) );
            }//m
        }//q
    }

    LOG(INFO) << "[CRBSCM<TruthModelType>::computeYBounds()] stop.\n";
    //std::cout << "************************************************************\n";
}

template<typename TruthModelType>
std::vector<double>
CRBSCM<TruthModelType>::run( parameter_type const& mu, int K )
{
    std::cout << "------------------------------------------------------------\n";
    double alpha_lb,alpha_lbti;
    boost::tie( alpha_lb, alpha_lbti ) = this->lb( mu, K );
    double alpha_ub,alpha_ubti;
    if( M_use_scm )
        boost::tie( alpha_ub, alpha_ubti ) = this->ub( mu, K );
    else
    {
        alpha_ub = alpha_lb;
        alpha_ubti = alpha_lbti;
    }
    double alpha_ex, alpha_exti;
    boost::tie( alpha_ex, alpha_exti ) = this->ex( mu );
    LOG( INFO ) << "alpha_lb=" << alpha_lb << " alpha_ub=" << alpha_ub << " alpha_ex=" << alpha_ex << "\n";
    LOG( INFO ) << ( alpha_ex-alpha_lb )/( alpha_ub-alpha_lb ) << "\n";
    LOG( INFO ) << K << " "
                << std::setprecision( 16 ) << alpha_lb << " "
                << std::setprecision( 3 ) << alpha_lbti << " "
                << std::setprecision( 16 ) << alpha_ub << " "
                << std::setprecision( 3 ) << alpha_ubti << " "
                << std::setprecision( 16 ) << alpha_ex << " "
                << std::setprecision( 16 ) << alpha_exti << " "
                << std::setprecision( 16 ) << ( alpha_ub-alpha_lb )/( alpha_ub ) << " "
                << std::setprecision( 16 ) << ( alpha_ex-alpha_lb )/( alpha_ex ) << " "
                << std::setprecision( 16 ) << ( alpha_ub-alpha_ex )/( alpha_ex ) << " "
        ;
    LOG( INFO ) << "------------------------------------------------------------\n";
    double rel_diff = (alpha_ex - alpha_lb)/alpha_ex;
    return boost::assign::list_of( alpha_lb )( alpha_lbti )( alpha_ub )( alpha_ubti )( alpha_ex )( alpha_exti )( rel_diff );
}

template<typename TruthModelType>
void
CRBSCM<TruthModelType>::run( const double * X, unsigned long N, double * Y, unsigned long P )
{
    parameter_type mu( M_Dmu );

    // the last parameter is the max error
    for ( unsigned long p= 0; p < N-3; ++p )
        mu( p ) = X[p];

    for ( unsigned long i=0; i<N; i++ ) std::cout<<"X["<<i<<"] = "<<X[i]<<std::endl;

    double meshSize  = X[N-3];
    //M_model->setMeshSize( meshSize );

    size_type K = this->KMax();
    double alpha_lb,lbti;
    boost::tie( alpha_lb, lbti ) = this->lb( mu, K );
    double alpha_ub,ubti;
    boost::tie( alpha_ub, ubti ) = this->ub( mu, K );
    double alpha_ex, alpha_exti;
    boost::tie( alpha_ex, alpha_exti ) = this->ex( mu );
    std::cout << "lb=" << alpha_lb << " ub=" << alpha_ub << " ex=" << alpha_ex << "\n";
    std::cout << ( alpha_ex-alpha_lb )/( alpha_ub-alpha_lb ) << "\n";
    Y[0] = alpha_lb;
    Y[1] = alpha_ub;

}

template<typename TruthModelType>
template<class Archive>
void
CRBSCM<TruthModelType>::save( Archive & ar, const unsigned int version ) const
{
    if( M_use_scm )
    {
        ar & M_use_scm;
        ar & M_Malpha;
        ar & M_Mplus;
        ar & M_C_alpha_lb;
        ar & M_C;
        ar & M_C_complement;
        ar & M_C_eigenvalues;
        ar & M_y_bounds_0;
        ar & M_y_bounds_1;
        ar & M_Y_ub;
        ar & M_Xi;
    }
    else
    {
        ar & M_use_scm;
        ar & M_mu_ref;
        ar & M_C_eigenvalues;
    }
}

template<typename TruthModelType>
template<class Archive>
void
CRBSCM<TruthModelType>::load( Archive & ar, const unsigned int version )
{

    ar & M_use_scm ;
    bool use_scm =  boption(_name="crb.scm.use-scm") ;
    bool rebuild =  boption(_name="crb.scm.rebuild-database") ;
    if( M_use_scm != use_scm && rebuild==false)
    {
        if( use_scm )
            throw std::logic_error( "[CRBSCM::load] ERROR the database created is not appropriate to use SCM. Use option crb.scm.rebuild-database=true");
        else
            throw std::logic_error( "[CRBSCM::load] ERROR the database was created to use SCM. Use option crb.scm.rebuild-database=true");
    }

    if( M_use_scm )
    {
        ar & M_Malpha;
        ar & M_Mplus;
        ar & M_C_alpha_lb;
        ar & M_C;
        ar & M_C_complement;
        ar & M_C_eigenvalues;
        ar & M_y_bounds_0;
        ar & M_y_bounds_1;
        ar & M_Y_ub;
        ar & M_Xi;
        int Qmax = this->nb_decomposition_terms_q();
        M_y_bounds.resize( Qmax );
        for ( int q=0; q<Qmax; q++ )
        {
            for(int m=0; m<mMax(q); m++)
                M_y_bounds[q].push_back( boost::make_tuple( M_y_bounds_0[q][m] , M_y_bounds_1[q][m] ) );
        }
    }
    else
    {
        ar & M_mu_ref;
        ar & M_C_eigenvalues;
    }

}

template<typename TruthModelType>
bool
CRBSCM<TruthModelType>::doScmForMassMatrix()
{
    bool b = this->vm()["crb.scm.do-scm-for-mass-matrix"].template as<bool>();
    return b;
}

template<typename TruthModelType>
bool
CRBSCM<TruthModelType>::rebuildDB()
{
    bool rebuild = this->vm()["crb.scm.rebuild-database"].template as<bool>();
    return rebuild;
}

template<typename TruthModelType>
void
CRBSCM<TruthModelType>::saveDB()
{
    fs::ofstream ofs( this->dbLocalPath() / this->dbFilename() );

    if ( ofs )
    {
        boost::archive::text_oarchive oa( ofs );
        // write class instance to archive
        oa << *this;
        // archive and stream closed when destructors are called
    }
}
template<typename TruthModelType>
bool
CRBSCM<TruthModelType>::loadDB()
{
    if ( this->rebuildDB() )
        return false;

    fs::path db = this->lookForDB();

    if ( db.empty() )
        return false;

    if ( !fs::exists( db ) )
        return false;

    LOG( INFO ) << "Loading " << db << "...";
    fs::ifstream ifs( db );
    if ( ifs )
    {
        boost::archive::text_iarchive ia( ifs );
        // write class instance to archive
        ia >> *this;
        LOG( INFO ) << "Loading " << db << " done...";
        this->setIsLoaded( true );
        // archive and stream closed when destructors are called
        return true;
    }

    return false;
}

} // Feel

namespace boost
{
namespace serialization
{
template< typename T>
struct version< Feel::CRBSCM<T> >
{
    // at the moment the version of the CRBSCM DB is 0. if any changes is done
    // to the format it is mandatory to increase the version number below
    // and use the new version number of identify the new entries in the DB
    typedef mpl::int_<0> type;
    typedef mpl::integral_c_tag tag;
    static const unsigned int value = version::type::value;
};
template<typename T> const unsigned int version<Feel::CRBSCM<T> >::value;
}
}
#endif /* __CRBSCM_H */
