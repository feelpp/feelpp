/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-11-24

  Copyright (C) 2009-2012 Université Joseph Fourier (Grenoble I)

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
   \file crb_trilinear.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Elisa Schenone
   \author Stephane Veys
   \date 2012-11-10
 */
#ifndef __CRBTrilinear_H
#define __CRBTrilinear_H 1

#include <boost/multi_array.hpp>
#include <boost/tuple/tuple.hpp>
#include "boost/tuple/tuple_io.hpp"
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/bimap.hpp>
#include <boost/bimap/support/lambda.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <fstream>


#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_member.hpp>

#include <vector>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>

#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/parameter.hpp>
#include <feel/feelcrb/parameterspace.hpp>
#include <feel/feelcrb/crbdb.hpp>
#include <feel/feelcrb/crbscm.hpp>
#include <feel/feelcore/serialization.hpp>
#include <feel/feelfilters/exporter.hpp>

namespace Feel
{

/**
 * \class CRBTrilinear
 * \brief Certifed Reduced Basis for Trilinear forms class
 *
 * Implements the certified reduced basis method for treat trilinear forms
 *
 *
 * @author Elisa Schenone, Stephane Veys
 * @see
 */
template<typename TruthModelType>
class CRBTrilinear : public CRBDB
{
    typedef  CRBDB super;
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

    typedef ParameterSpace<TruthModelType::ParameterSpaceDimension> parameterspace_type;
    typedef boost::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef typename parameterspace_type::element_type parameter_type;
    typedef typename parameterspace_type::element_ptrtype parameter_ptrtype;
    typedef typename parameterspace_type::sampling_type sampling_type;
    typedef typename parameterspace_type::sampling_ptrtype sampling_ptrtype;

    typedef boost::bimap< int, boost::tuple<double,double,double> > convergence_type;

    typedef double value_type;

    typedef typename convergence_type::value_type convergence;

    //! function space type
    typedef typename model_type::functionspace_type functionspace_type;
    typedef typename model_type::functionspace_ptrtype functionspace_ptrtype;

    //! element of the functionspace type
    typedef typename model_type::element_type element_type;
    typedef typename model_type::element_ptrtype element_ptrtype;

    typedef typename model_type::backend_type backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename model_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename model_type::vector_ptrtype vector_ptrtype;
    typedef typename model_type::beta_vector_type beta_vector_type;

    typedef std::vector<element_type> wn_type;
    typedef boost::tuple< std::vector<wn_type> , std::vector<std::string> > export_vector_wn_type;

    typedef Eigen::VectorXd vectorN_type;
    typedef Eigen::MatrixXd matrixN_type;

    //! mesh type
    typedef typename model_type::mesh_type mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //! space type
    typedef typename model_type::space_type space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;

    typedef Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > map_dense_matrix_type;
    typedef Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > map_dense_vector_type;

    // ! export
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    typedef CRBTrilinear self_type;


    //! scm
    typedef CRBSCM<truth_model_type> scm_type;
    typedef boost::shared_ptr<scm_type> scm_ptrtype;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    CRBTrilinear()
        :
        super(),
        M_nlsolver( SolverNonLinear<double>::build( SOLVERS_PETSC, Environment::worldComm() ) ),
        M_model(),
        M_output_index( 0 ),
        M_tolerance( 1e-2 ),
        M_iter_max( 3 ),
        M_error_type( CRB_NO_RESIDUAL ),
        M_Dmu( new parameterspace_type ),
        M_Xi( new sampling_type( M_Dmu ) ),
        M_WNmu( new sampling_type( M_Dmu, 1, M_Xi ) ),
        M_WNmu_complement(),
        M_scm( new scm_type( name, vm ) ),
        exporter( Exporter<mesh_type>::New( "ensight" ) )
    {

    }

    //! constructor from command line options
    CRBTrilinear( std::string  name, po::variables_map const& vm )
        :
        super( ( boost::format( "%1%" ) % vm["crb.error-type"].template as<int>() ).str(),
               name,
               ( boost::format( "%1%-%2%-%3%" ) % name % vm["crb.output-index"].template as<int>() % vm["crb.error-type"].template as<int>() ).str(),
               vm ),
        M_nlsolver( SolverNonLinear<double>::build( SOLVERS_PETSC, Environment::worldComm() ) ),
        M_model(),
        M_backend( backend_type::build( vm ) ),
        M_output_index( vm["crb.output-index"].template as<int>() ),
        M_tolerance( vm["crb.error-max"].template as<double>() ),
        M_iter_max( vm["crb.dimension-max"].template as<int>() ),
        M_error_type( CRBErrorType( vm["crb.error-type"].template as<int>() ) ),
        M_Dmu( new parameterspace_type ),
        M_Xi( new sampling_type( M_Dmu ) ),
        M_WNmu( new sampling_type( M_Dmu, 1, M_Xi ) ),
        M_WNmu_complement(),
        exporter( Exporter<mesh_type>::New( vm, "BasisFunction" ) )
    {
    }

    //! constructor from command line options
    CRBTrilinear( std::string  name,
         po::variables_map const& vm,
         truth_model_ptrtype const & model )
        :
        super( ( boost::format( "%1%" ) % vm["crb.error-type"].template as<int>() ).str(),
               name,
               ( boost::format( "%1%-%2%-%3%" ) % name % vm["crb.output-index"].template as<int>() % vm["crb.error-type"].template as<int>() ).str(),
               vm ),
        M_nlsolver( SolverNonLinear<double>::build( SOLVERS_PETSC, Environment::worldComm() ) ),
        M_model(),
        M_backend( backend_type::build( vm ) ),
        M_output_index( vm["crb.output-index"].template as<int>() ),
        M_tolerance( vm["crb.error-max"].template as<double>() ),
        M_iter_max( vm["crb.dimension-max"].template as<int>() ),
        M_error_type( CRBErrorType( vm["crb.error-type"].template as<int>() ) ),
        M_Dmu( new parameterspace_type ),
        M_Xi( new sampling_type( M_Dmu ) ),
        M_WNmu( new sampling_type( M_Dmu, 1, M_Xi ) ),
        M_WNmu_complement(),
        M_scm( new scm_type( name, vm ) ),
        exporter( Exporter<mesh_type>::New( vm, "BasisFunction" ) )
    {
        this->setTruthModel( model );
        if ( this->loadDB() )
            LOG(INFO) << "Database " << this->lookForDB() << " available and loaded\n";
    }


    //! copy constructor
    CRBTrilinear( CRBTrilinear const & o )
        :
        super( o ),
        M_output_index( o.M_output_index ),
        M_tolerance( o.M_tolerance ),
        M_iter_max( o.M_iter_max ),
        M_error_type( o.M_error_type ),
        M_maxerror( o.M_maxerror ),
        M_Dmu( o.M_Dmu ),
        M_Xi( o.M_Xi ),
        M_WNmu( o.M_WNmu ),
        M_WNmu_complement( o.M_WNmu_complement )
    {}

    //! destructor
    ~CRBTrilinear()
    {}

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * find the nearest neighbor of mu in the sampling WNmu
     * return th neighbor and the index of this neighbor in WNmu
     */
    void findNearestNeighborInWNmu( parameter_type const& mu, parameter_type & neighbor, int & index ) const;

    /**
     * Returns the lower bound of the output
     *
     * \param mu \f$ \mu\f$ the parameter at which to evaluate the output
     * \param N the size of the reduced basis space to use
     * \param uN primal solution
     *
     *\return compute online the lower bound
     *\and also condition number of matrix A
     */

    boost::tuple<double,double> lb( size_type N, parameter_type const& mu, vectorN_type& uN ) const;

    /**
     * Offline computation
     *
     * \return the convergence history (max error)
     */
    convergence_type offline();

    /**
     * \param mu : parameters
     * \param N : dimension of the reduced basis used
     */
    void updateLinearTerms( parameter_type const& mu , int N ) const;

    /**
     * Update the Jacobian Matrix for Newton Solver
     *
     */
    void updateJacobian( const map_dense_vector_type& X, map_dense_matrix_type& J , parameter_type const& mu , int N ) const;

    /**
     * Update the Residual of the Newton Solver
     *
     */
    void updateResidual( const map_dense_vector_type& X, map_dense_vector_type& R , parameter_type const& mu , int N ) const;


    void displayVector(const map_dense_vector_type& V ) const ;
    void displayVector(const vectorN_type& V ) const ;
    void displayMatrix(const matrixN_type& M ) const ;
    // -------- from crb

    void orthonormalize( size_type N, wn_type& wn, int Nm = 1 );
    void checkOrthonormality( int N, const wn_type& wn ) const;

    /**
     * if true, print the max error (absolute) during the offline stage
     */
    bool printErrorDuringOfflineStep()
    {
        bool print = this->vm()["crb.print-error-during-rb-construction"].template as<bool>();
        return print;
    }

    CRBErrorType errorType() const
    {
        return CRB_NO_RESIDUAL;
    }

    //! \return the scm object (only to compile)
    scm_ptrtype scm() const
    {
        return M_scm;
    }

    sampling_ptrtype wnmu ( ) const
    {
        return M_WNmu;
    }

    bool useWNmu()
    {
        bool use = this->vm()["crb.run-on-WNmu"].template as<bool>();
        return use;
    }



    /**
     * print max errors (total error and also primal and dual contributions)
     * during offline stage
     * this function is here only to compile !!!
     */
    void printErrorsDuringRbConstruction( void ){};

    /**
     * return the crb expansion at parameter \p \mu, ie \f$\sum_{i=0}^N u^N_i
     * \phi_i\f$ where $\phi_i, i=1...N$ are the basis function of the reduced
     * basis space
     * if N>0 take the N^th first elements, else take all elements
     */
    element_type expansion( parameter_type const& mu , int N=-1);

    /**
     * return the crb expansion at parameter \p \mu, ie \f$\sum_{i=0}^N u^N_i
     * \phi_i\f$ where $\phi_i, i=1...N$ are the basis function of the reduced
     * basis space
     */
    element_type expansion( vectorN_type const& u , int const N) const;


    /**
     * export basis functions to visualize it
     * \param wn : tuple composed of a vector of wn_type and a vector of string (used to name basis)
     */
    void exportBasisFunctions( const export_vector_wn_type& wn )const ;

    boost::tuple<double,double,double,double, vectorN_type > run( parameter_type const& mu, double eps = 1e-6, int N = -1 );
    void run( const double * X, unsigned long N, double * Y, unsigned long P ){};

    //! set the truth offline model
    void setTruthModel( truth_model_ptrtype const& model )
    {
        M_model = model;
        M_Dmu = M_model->parameterSpace();
        M_Xi = sampling_ptrtype( new sampling_type( M_Dmu ) );

        if ( ! loadDB() )
            M_WNmu = sampling_ptrtype( new sampling_type( M_Dmu ) );
        else
        {
            LOG(INFO) << "Database " << this->lookForDB() << " available and loaded\n";
        }
        M_scm->setTruthModel( M_model );
    }


    /**
     * if true, show the mu selected during the offline stage
     */
    bool showMuSelection()
    {
        bool show = this->vm()["crb.show-mu-selection"].template as<bool>();
        return show;
    }

    /**
     * print parameters set mu selected during the offline stage
     */
    void printMuSelection( void );

    //! \return the dimension of the reduced basis space
    int dimension() const
    {
        return M_N;
    }

    bool rebuildDB()
    {
        bool rebuild = this->vm()["crb.rebuild-database"].template as<bool>();
        return rebuild;
    }

    /**
     * save the CRB database
     */
    void saveDB();

    /**
     * load the CRB database
     */
    bool loadDB();

    WorldComm const& worldComm() const { return Environment::worldComm() ; }
    //@}


private:
    std::vector < std::vector < matrixN_type> >  M_Aqm_tril_pr;
    mutable matrixN_type M_bilinear_terms;
    mutable vectorN_type M_linear_terms;
    boost::shared_ptr<SolverNonLinear<double> > M_nlsolver;


    // ------ from crb

    truth_model_ptrtype M_model;
    backend_ptrtype M_backend;
    int M_output_index;
    double M_tolerance;
    size_type M_iter_max;
    CRBErrorType M_error_type;
    double M_maxerror;
    // parameter space
    parameterspace_ptrtype M_Dmu;
    // fine sampling of the parameter space
    sampling_ptrtype M_Xi;
    // sampling of parameter space to build WN
    sampling_ptrtype M_WNmu;
    sampling_ptrtype M_WNmu_complement;

    scm_ptrtype M_scm;

    //export
    export_ptrtype exporter;


    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void save( Archive & ar, const unsigned int version ) const;

    template<class Archive>
    void load( Archive & ar, const unsigned int version ) ;

    BOOST_SERIALIZATION_SPLIT_MEMBER()

    // reduced basis space
    wn_type M_WN;

    size_type M_N;

    bool orthonormalize_primal;
    bool orthonormalize_dual;
    bool solve_dual_problem;

    convergence_type M_rbconv;

    // left hand side ( bilinear )
    std::vector < std::vector<matrixN_type> > M_Aqm_pr;

    // right hand side
    std::vector < std::vector<vectorN_type> > M_Fqm_pr;
    // output
    std::vector < std::vector<vectorN_type> > M_Lqm_pr;

    std::vector<int> M_index;
    int M_mode_number;

    parameter_type M_current_mu;
    int M_no_residual_index;


};

template<typename TruthModelType>
typename CRBTrilinear<TruthModelType>::convergence_type
CRBTrilinear<TruthModelType>::offline()
{

    int proc_number = this->worldComm().globalRank();

    bool rebuild_database = this->vm()["crb.rebuild-database"].template as<bool>() ;
    orthonormalize_primal = this->vm()["crb.orthonormalize-primal"].template as<bool>() ;

    boost::timer ti;
    if( proc_number == 0 ) std::cout << "Offline CRB starts, this may take a while until Database is computed...\n";
    LOG(INFO) << "[CRB::offline] Starting offline for output " << M_output_index << "\n";
    LOG(INFO) << "[CRB::offline] initialize underlying finite element model\n";
    M_model->init();
    if( proc_number == 0 ) std::cout << " -- model init done in " << ti.elapsed() << "s\n";

    parameter_type mu( M_Dmu );

    double delta_pr;
    double delta_du;
    size_type index;
    //if M_N == 0 then there is not an already existing database
    if ( rebuild_database || M_N == 0)
    {

        ti.restart();

        LOG(INFO) << "[CRB::offline] compute random sampling\n";

        int sampling_size = this->vm()["crb.sampling-size"].template as<int>();
        std::string file_name = ( boost::format("M_Xi_%1%") % sampling_size ).str();
        std::ifstream file ( file_name );
        if( ! file )
        {
            // random sampling
            M_Xi->randomize( sampling_size );
            //M_Xi->equidistribute( this->vm()["crb.sampling-size"].template as<int>() );
            M_Xi->writeOnFile(file_name);
        }
        else
        {
            M_Xi->clear();
            M_Xi->readFromFile(file_name);
        }

        M_WNmu->setSuperSampling( M_Xi );

        if( proc_number == 0 ) std::cout<<"[CRB offline] M_error_type = "<<M_error_type<<std::endl;

        std::cout << " -- sampling init done in " << ti.elapsed() << "s\n";
        ti.restart();

        // empty sets
        M_WNmu->clear();
        if( M_error_type == CRB_NO_RESIDUAL )
            mu = M_Dmu->element();
        else
        {
            // start with M_C = { arg min mu, mu \in Xi }
            boost::tie( mu, index ) = M_Xi->min();
        }


        int size = mu.size();
        if( proc_number == 0 )
        {
            std::cout << "  -- start with mu = [ ";
            for ( int i=0; i<size-1; i++ ) std::cout<<mu( i )<<" ";
            std::cout<<mu( size-1 )<<" ]"<<std::endl;
        }
        //std::cout << " -- WN size :  " << M_WNmu->size() << "\n";

        // dimension of reduced basis space
        M_N = 0;

        M_maxerror = 1e10;
        delta_pr = 0;
        delta_du = 0;
        //boost::tie( M_maxerror, mu, index ) = maxErrorBounds( N );

        LOG(INFO) << "[CRB::offline] allocate reduced basis data structures\n";
        M_Aqm_pr.resize( M_model->Qa() );
        for(int q=0; q<M_model->Qa(); q++)
        {
            M_Aqm_pr[q].resize( 1 );
        }

        M_Aqm_tril_pr.resize( M_model->QaTri() );
        for(int q=0; q<M_model->QaTri(); q++)

        M_Fqm_pr.resize( M_model->Ql( 0 ) );

        for(int q=0; q<M_model->Ql( 0 ); q++)
        {
            M_Fqm_pr[q].resize( 1 );
        }

        M_Lqm_pr.resize( M_model->Ql( M_output_index ) );
        for(int q=0; q<M_model->Ql( M_output_index ); q++)
            M_Lqm_pr[q].resize( 1 );

    }//end of if( rebuild_database )
#if 1
    else
    {
        mu = M_current_mu;
        if( proc_number == 0 )
        {
            std::cout<<"we are going to enrich the reduced basis"<<std::endl;
            std::cout<<"there are "<<M_N<<" elements in the database"<<std::endl;
        }
        LOG(INFO) <<"we are going to enrich the reduced basis"<<std::endl;
        LOG(INFO) <<"there are "<<M_N<<" elements in the database"<<std::endl;
    }//end of else associated to if ( rebuild_databse )
#endif

    LOG(INFO) << "[CRBTrilinear::offline] compute affine decomposition\n";
    std::vector< std::vector<sparse_matrix_ptrtype> > Aqm;
    std::vector< std::vector<sparse_matrix_ptrtype> > Aqm_tril;
    std::vector< std::vector<std::vector<vector_ptrtype> > > Fqm;

    boost::tie( boost::tuples::ignore, Aqm, Fqm , boost::tuples::ignore ) = M_model->computeAffineDecomposition();

    element_ptrtype u( new element_type( M_model->functionSpace() ) );

    LOG(INFO) << "[CRBTrilinear::offline] starting offline adaptive loop\n";

    bool reuse_prec = this->vm()["crb.reuse-prec"].template as<bool>() ;

    bool use_predefined_WNmu = this->vm()["crb.use-predefined-WNmu"].template as<bool>() ;
    int N_log_equi = this->vm()["crb.use-logEquidistributed-WNmu"].template as<int>() ;
    int N_equi = this->vm()["crb.use-equidistributed-WNmu"].template as<int>() ;

    if( N_log_equi > 0 || N_equi > 0 )
        use_predefined_WNmu = true;

    if ( use_predefined_WNmu )
    {
        std::string file_name = ( boost::format("SamplingWNmu") ).str();
        std::ifstream file ( file_name );
        if( ! file )
        {
            throw std::logic_error( "[CRBTrilinear::offline] ERROR the file SamplingWNmu doesn't exist so it's impossible to known which parameters you want to use to build the database" );
        }
        else
        {
            M_WNmu->clear();
            int sampling_size = M_WNmu->readFromFile(file_name);
            M_iter_max = sampling_size;
        }
        mu = M_WNmu->at( M_N ); // first element
        //std::cout<<" [use_predefined_WNmu] mu = \n"<<mu<<std::endl;

        if( proc_number == this->worldComm().masterRank() )
            std::cout<<"[CRB::offline] read WNmu ( sampling size : "<<M_iter_max<<" )"<<std::endl;

    }

    if( M_error_type == CRB_NO_RESIDUAL || use_predefined_WNmu )
    {
        //in this case it makes no sens to check the estimated error
        M_maxerror = 1e10;
    }


    LOG(INFO) << "[CRBTrilinear::offline] strategy "<< M_error_type <<"\n";
    if( proc_number == 0 ) std::cout << "[CRBTrilinear::offline] strategy "<< M_error_type <<"\n";

    while ( M_maxerror > M_tolerance && M_N < M_iter_max )
    {

        boost::timer timer, timer2;
        LOG(INFO) <<"========================================"<<"\n";

        if( proc_number == 0 )
        {
            std::cout << "============================================================"<<std::endl;
            std::cout << "N=" << M_N << "/"  << M_iter_max <<" ( nb proc : "<<worldComm().globalSize()<<")"<<std::endl;
        }
        LOG(INFO) << "N=" << M_N << "/"  << M_iter_max << "\n";

        // for a given parameter \p mu assemble the left and right hand side
        u->setName( ( boost::format( "fem-primal-N%1%-proc%2%" ) % (M_N)  % proc_number ).str() );

        mu.check();
        u->zero();

        timer2.restart();

        LOG(INFO) << "[CRB::offline] solving primal" << "\n";
        *u = M_model->solve( mu );

        if( proc_number == this->worldComm().masterRank() ) std::cout << "  -- primal problem solved in " << timer2.elapsed() << "s\n";
        timer2.restart();


        if( ! use_predefined_WNmu )
            M_WNmu->push_back( mu, index );

        M_WNmu_complement = M_WNmu->complement();

        M_WN.push_back( *u );
	    int number_of_added_elements=1;
        M_N+=number_of_added_elements;

        if ( orthonormalize_primal )
        {
            orthonormalize( M_N, M_WN, number_of_added_elements );
            orthonormalize( M_N, M_WN, number_of_added_elements );
            orthonormalize( M_N, M_WN, number_of_added_elements );
        }

        LOG(INFO) << "[CRB::offline] compute Aq_pr, Aq_du, Aq_pr_du" << "\n";
        for  (size_type q = 0; q < M_model->Qa(); ++q )
        {
            M_Aqm_pr[q][0].conservativeResize( M_N, M_N );

            // only compute the last line and last column of reduced matrices
            for ( size_type i = M_N-number_of_added_elements; i < M_N; i++ )
            {
                for ( size_type j = 0; j < M_N; ++j )
                {
                    M_Aqm_pr[q][0]( i, j ) = Aqm[q][0]->energy( M_WN[i], M_WN[j] );
                }
            }

            for ( size_type j=M_N-number_of_added_elements; j < M_N; j++ )
            {
                for ( size_type i = 0; i < M_N; ++i )
                {
                    M_Aqm_pr[q][0]( i, j ) = Aqm[q][0]->energy( M_WN[i], M_WN[j] );
                }
            }
        }//loop over q


        LOG(INFO) << "[CRBTrilinear::offline] compute Fq_pr" << "\n";

        for ( size_type q = 0; q < M_model->Ql( 0 ); ++q )
        {
            M_Fqm_pr[q][0].conservativeResize( M_N );

            for ( size_type l = 1; l <= number_of_added_elements; ++l )
            {
                int index = M_N-l;
                M_Fqm_pr[q][0]( index ) = M_model->Fqm( 0, q, 0, M_WN[index] );
            }
        }//loop over q

        LOG(INFO) << "[CRB::offline] compute Lq_pr" << "\n";

        for ( size_type q = 0; q < M_model->Ql( M_output_index ); ++q )
        {
            M_Lqm_pr[q][0].conservativeResize( M_N );

            for ( size_type l = 1; l <= number_of_added_elements; ++l )
            {
                int index = M_N-l;
                M_Lqm_pr[q][0]( index ) = M_model->Fqm( M_output_index, q, 0, M_WN[index] );
            }
        }//loop over q

        sparse_matrix_ptrtype trilinear_form;
        for  (size_type q = 0; q < M_model->QaTri(); ++q )
        {
            M_Aqm_tril_pr[q].resize( M_N );
            for (int k=0 ; k<M_N; k++)
            {
                //bring back the matrix associated to the trilinear form for a given basis function
                //we do this here to use only one matrix
                trilinear_form  = M_model->computeTrilinearForm( M_WN[k] );

                M_Aqm_tril_pr[q][k].conservativeResize( M_N, M_N );
                for ( int i = 0; i < M_N; ++i )
                {
                    for ( int j = 0; j < M_N; ++j )
                    {
                        M_Aqm_tril_pr[q][k]( i, j ) = trilinear_form->energy( M_WN[j], M_WN[i] );
                    }//j
                }//i
            }//k
        }// q

        timer2.restart();

        if ( ! use_predefined_WNmu )
        {
            bool already_exist;
            do
            {
                //initialization
                already_exist=false;
                //pick randomly an element
                mu = M_Dmu->element();
                //make sure that the new mu is not already is M_WNmu
                BOOST_FOREACH( auto _mu, *M_WNmu )
                {
                    if( mu == _mu )
                        already_exist=true;
                }
            }
            while( already_exist );
            M_current_mu = mu;
        }
        else
        {
            //remmber that in this case M_iter_max = sampling size
            if( M_N < M_iter_max )
            {
                mu = M_WNmu->at( M_N );
                M_current_mu = mu;
            }
        }

        M_rbconv.insert( convergence( M_N, boost::make_tuple(M_maxerror,delta_pr,delta_du) ) );

        timer2.restart();
        LOG(INFO) << "time: " << timer.elapsed() << "\n";
        if( proc_number == 0 ) std::cout << "============================================================\n";
        LOG(INFO) <<"========================================"<<"\n";

        //save DB after adding an element
        this->saveDB();
    }


    if( proc_number == 0 )
        std::cout<<"number of elements in the reduced basis : "<<M_N<<" ( nb proc : "<<worldComm().globalSize()<<")"<<std::endl;
    LOG(INFO) << " index choosen : ";
    BOOST_FOREACH( auto id, M_index )
    LOG(INFO)<<id<<" ";
    LOG(INFO)<<"\n";
    bool visualize_basis = this->vm()["crb.visualize-basis"].template as<bool>() ;

    if ( visualize_basis )
    {
        std::vector<wn_type> wn;
        std::vector<std::string> names;
        wn.push_back( M_WN );
        names.push_back( "primal" );
        exportBasisFunctions( boost::make_tuple( wn ,names ) );

        if ( orthonormalize_primal )
        {
            std::cout<<"[CRB::offline] Basis functions have been exported but warning elements have been orthonormalized"<<std::endl;
        }
    }

    if( proc_number == 0 ) std::cout << "Offline CRB is done\n";

    return M_rbconv;

}


template<typename TruthModelType>
void
CRBTrilinear<TruthModelType>::findNearestNeighborInWNmu( parameter_type const& mu, parameter_type & neighbor, int & index ) const
{
    std::vector<int> index_vector;
    sampling_ptrtype S =  M_WNmu->searchNearestNeighbors( mu, 1 , index_vector);
    neighbor = S->at( 0 );
    index = index_vector[0];
    //std::cout<<"[CRBTrilinear::findNearestNeighborInWNmu] for Gr = "<<mu(0)<<" th nearest neighbor in WNmu is "<<neighbor(0)<<" at index "<<index<<std::endl;
}

template<typename TruthModelType>
boost::tuple<double,double>
CRBTrilinear<TruthModelType>::lb( size_type N, parameter_type const& mu, vectorN_type & uN ) const
{

    std::cout<<"\nCRBTrilinear::lb starts with N = "<<N<<std::endl;
    google::FlushLogFiles(google::GLOG_INFO);

    if ( N > M_N ) N = M_N;
    beta_vector_type betaAqm;
    std::vector<beta_vector_type> betaFqm, betaLqm;

    //-- end of initialization step

    //vector containing outputs from time=time_step until time=time_for_output
    double output;

    boost::tie( boost::tuples::ignore,  betaAqm, betaFqm  , boost::tuples::ignore ) = M_model->computeBetaQm( mu );

    google::FlushLogFiles(google::GLOG_INFO);

    /*
     ------> Here add Newton Method and continuous Gr,Pr
     */
    using namespace vf;
    //Feel::ParameterSpace<2>::Element current_mu( mu );
    parameter_type current_mu ;

    backend_ptrtype backend_primal_problem = backend_type::build( BACKEND_PETSC );

    vectorN_type R( (int) N );
    matrixN_type J( (int) N , (int) N );

    //initialization of uN
    //we look for the nearest neighbor of mu in the sampling WNmu
    //let i the index of this neighbor in WNmu, we will set zeros in uN except at the i^th component where we will set 1
    uN.setZero( (int) N );
    parameter_type neighbor( M_Dmu );
    int index;
    findNearestNeighborInWNmu(  mu,  neighbor, index );
    if( this->vm()["crb.cvg-study"].template as<bool>() == true )
    {
        //in this case, index may be smaller than uN.size
        //so we do nothing
    }
    else
        uN( index ) = 1;

    double *r_data = R.data();
    double *j_data = J.data();
    double *uN_data = uN.data();

    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic , 1> > map_R ( r_data, N );
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic , 1> > map_uN ( uN_data, N );
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic , Eigen::Dynamic> > map_J ( j_data, N , N );

    double gr = mu( 0 );
    double pr = mu( 1 );

    bool use_continuity = this->vm()["crb.use-continuity"].template as<bool>();
    int Nmax=1;

    if( use_continuity )
        Nmax=std::max( 1.0,std::max( std::ceil( std::log( gr ) ),std::ceil( std::log( pr )-std::log( 1.e-2 ) ) ) );

    for ( int i = 0; i < Nmax; ++i )
    {

        double current_Grashofs;
        double current_Prandtl;
        if( use_continuity )
        {
            int denom = ( Nmax==1 )?1:Nmax-1;
            current_Grashofs = math::exp( math::log( 1. )+i*( math::log( gr )-math::log( 1. ) )/denom );
            current_Prandtl = math::exp( math::log( 1.e-2 )+i*( math::log( pr )-math::log( 1.e-2 ) )/denom );
            std::cout << "[CRBTrilinear::lb] i/N = " << i+1 << "/" << Nmax <<std::endl;
            std::cout << "[CRBTrilinear::lb] intermediary Grashof = " << current_Grashofs<<std::endl;
            std::cout << "[CRBTrilinear::lb] and Prandtl = " << current_Prandtl <<" \n" <<std::endl;
        }
        else
        {
            current_Grashofs = gr;
            current_Prandtl = pr;
        }

        current_mu << current_Grashofs, current_Prandtl;

        this->updateLinearTerms( current_mu , N );

        //M_nlsolver->setRelativeResidualTol( 1e-12 );
        M_nlsolver->map_dense_jacobian = boost::bind( &self_type::updateJacobian, boost::ref( *this ), _1, _2  , current_mu , N );
        M_nlsolver->map_dense_residual = boost::bind( &self_type::updateResidual, boost::ref( *this ), _1, _2  , current_mu , N );
        M_nlsolver->setType( TRUST_REGION );
        M_nlsolver->solve( map_J , map_uN , map_R, 1e-12, 100);
    }

    LOG(INFO) << "[CRBTrilinear::lb] solve with Newton done";

    //compute conditioning of matrix J
    Eigen::SelfAdjointEigenSolver< matrixN_type > eigen_solver;
    eigen_solver.compute( map_J );
    int number_of_eigenvalues =  eigen_solver.eigenvalues().size();
    //we copy eigenvalues in a std::vector beacause it's easier to manipulate it
    std::vector<double> eigen_values( number_of_eigenvalues );
    for ( int i=0; i<number_of_eigenvalues; i++ )
    {
        if ( imag( eigen_solver.eigenvalues()[i] )>1e-12 )
        {
            throw std::logic_error( "[CRBTrilinear::lb] ERROR : complex eigenvalues were found" );
        }

        eigen_values[i]=real( eigen_solver.eigenvalues()[i] );
    }
    int position_of_largest_eigenvalue=number_of_eigenvalues-1;
    int position_of_smallest_eigenvalue=0;
    double eig_max = eigen_values[position_of_largest_eigenvalue];
    double eig_min = eigen_values[position_of_smallest_eigenvalue];
    double condition_number = eig_max / eig_min;
    //end of computation of conditionning

    LOG( INFO ) <<"[CRBTrilinear::lb] compute condition number of jacobian done\n";
    vectorN_type L ( ( int )N );
    L.setZero( N );

    for ( size_type q = 0; q < M_model->Ql( M_output_index ); ++q )
    {
        L += betaFqm[M_output_index][q][0]*M_Lqm_pr[q][0].head( N );
    }
    output = L.dot( uN );
    LOG(INFO) << "[CRBTrilinear::lb] computation of the output done";

    google::FlushLogFiles(google::GLOG_INFO);

    //std::cout<<"[CRBTrilinear uN] : \n"<<uN<<std::endl;

    return boost::make_tuple( output, condition_number );

}
template<typename TruthModelType>
void
CRBTrilinear<TruthModelType>::updateLinearTerms( parameter_type const& mu , int N ) const
{

    LOG(INFO) << "update linear terms \n";

    beta_vector_type betaAqm;
    std::vector<beta_vector_type> betaFqm, betaLqm;

    boost::tie( boost::tuples::ignore,  betaAqm, betaFqm , boost::tuples::ignore  ) = M_model->computeBetaQm( mu );

    google::FlushLogFiles(google::GLOG_INFO);

    M_bilinear_terms.setZero( N , N );

    for ( size_type q = 0; q < M_model->Qa(); ++q )
    {
        M_bilinear_terms += betaAqm[q][0]*M_Aqm_pr[q][0].block( 0, 0, N, N );
    }

    M_linear_terms.setZero( N );

    for ( size_type q = 0; q < M_model->Ql( 0 ); ++q )
    {
        M_linear_terms += betaFqm[0][q][0]*M_Fqm_pr[q][0].head( N );
    }

}
template<typename TruthModelType>
void
CRBTrilinear<TruthModelType>::updateJacobian( const map_dense_vector_type& map_X, map_dense_matrix_type& map_J , const parameter_type & mu , int N) const
{
    LOG(INFO) << "updateJacobian \n";
    map_J = M_bilinear_terms;

    bool enable = this->vm()["crb.enable-convection-terms"].template as<bool>();
    if( enable )
    {
        for ( size_type q = 0; q < M_model->QaTri(); ++q )
        {
            for (int k = 0 ; k < N; ++k)
            {
                for ( int i = 0; i < N; ++i )
                {
                    map_J( i, k ) += ( M_Aqm_tril_pr[q][k].row( i ).head( N ) ).dot(map_X);
                    map_J( i, k ) += ( M_Aqm_tril_pr[q][k].col( i ).head( N ) ).dot(map_X);
                }
            }
        }
    }

    if ( this->vm()["crb.compute-error-on-reduced-residual-jacobian"].template as<bool>() )
    {
        //bring the jacobian matrix from the model and then project it into the reduced basis
        auto expansionX = expansion( map_X , N);
        auto J = M_model->jacobian( expansionX );
        matrixN_type model_reduced_jacobian( N , N );
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
                model_reduced_jacobian(i,j) = J->energy( M_WN[i], M_WN[j]);
        }
        //compute difference
        matrixN_type diff = map_J - model_reduced_jacobian;
        double max = diff.maxCoeff();
        std::cout<<std::setprecision(14)<<"[CRB::updateJacobian] with X : "; this->displayVector(map_X); std::cout<<" the max coeff of the difference jacobian matrix : "<<max<<std::endl;
        if( math::abs(max) > 1e-14  )
        {
            std::cout<<"here is the jacobian matrix containing difference between reduced jacobian from CRB and the model";
            this->displayMatrix( diff );
        }
        std::cout<<std::endl;
    }
}



template<typename TruthModelType>
void
CRBTrilinear<TruthModelType>::updateResidual( const map_dense_vector_type& map_X, map_dense_vector_type& map_R , const parameter_type & mu, int N ) const
{
    LOG(INFO) << " updateResidual \n";

    matrixN_type temp ( N , N );
    temp.setZero( N , N );

    map_R = M_linear_terms;
    map_R += M_bilinear_terms * map_X ;

    bool enable = this->vm()["crb.enable-convection-terms"].template as<bool>();
    if( enable )
    {
        for ( size_type q = 0; q < M_model->QaTri(); ++q )
        {
            for (int k = 0 ; k < N; ++k)
            {
                for ( int i = 0; i < N; ++i )
                {
                    temp( k, i ) = map_X.dot( M_Aqm_tril_pr[q][k].row( i ).head( N ) ) ;
                }
            }
            map_R += temp * map_X ;
        }
    }

    if ( this->vm()["crb.compute-error-on-reduced-residual-jacobian"].template as<bool>() )
    {
        //bring the residual matrix from the model and then project it into the reduced basis
        auto expansionX = expansion( map_X , N);
        auto R = M_model->residual( expansionX );
        vectorN_type model_reduced_residual( N );
        element_ptrtype eltR( new element_type( M_model->functionSpace() ) );
        for(int i=0; i<eltR->localSize();i++)
            eltR->operator()(i)=R->operator()(i);
        for(int i=0; i<N; i++)
            model_reduced_residual(i) = inner_product( *eltR , M_WN[i] );
        //compute difference
        vectorN_type diff = map_R - model_reduced_residual;
        double max = diff.maxCoeff();
        std::cout<<std::setprecision(14)<<"[CRB::updateResidual] with X :"; this->displayVector(map_X); std::cout<<" Residual :";this->displayVector(map_R);
        std::cout<<" the max error  : "<<max<<std::endl;
        if( math::abs(max) > 1e-14  )
        {
            std::cout<<"here is the residual containing difference between reduced residual from CRB and the model";
            this->displayVector( diff );
        }
        std::cout<<std::endl;
    }
}

template<typename TruthModelType>
void
CRBTrilinear<TruthModelType>::displayVector( const map_dense_vector_type& V ) const
{
    int size=V.size();
    std::cout<<std::setprecision(14)<<" ( ";
    for(int i=0; i<size-1;i++)
        std::cout<<V(i)<<" , ";
    std::cout<<V(size-1)<<" ) ";
}

template<typename TruthModelType>
void
CRBTrilinear<TruthModelType>::displayVector( const vectorN_type& V ) const
{
    int size=V.size();
    std::cout<<std::setprecision(14)<<" ( ";
    for(int i=0; i<size-1;i++)
        std::cout<<V(i)<<" , ";
    std::cout<<V(size-1)<<" ) ";
}

template<typename TruthModelType>
void
CRBTrilinear<TruthModelType>::displayMatrix( const matrixN_type& M ) const
{
    std::cout<<std::setprecision(14)<<"\n[ ";
    int cols = M.cols();
    int rows = M.rows();
    for(int i=0; i<rows-1;i++)
    {
        for(int j=0; j<cols-1;j++)
            std::cout<<M(i,j)<<" , ";
        std::cout<<M(i,cols-1)<<" ]"<<std::endl;std::cout<<"[ ";
    }
    for(int j=0; j<cols-1;j++)
        std::cout<<M(rows-1,j)<<" , ";
    std::cout<<M(rows-1,cols-1)<<" ]"<<std::endl;std::cout<<" ";


}


//---------------- from crb

template<typename TruthModelType>
void
CRBTrilinear<TruthModelType>::orthonormalize( size_type N, wn_type& wn, int Nm )
{
    int proc_number = this->worldComm().globalRank();
    if( proc_number == 0 ) std::cout << "  -- orthonormalization (Gram-Schmidt)\n";
    DVLOG(2) << "[CRB::orthonormalize] orthonormalize basis for N=" << N << "\n";
    DVLOG(2) << "[CRB::orthonormalize] orthonormalize basis for WN="
                    << wn.size() << "\n";
    DVLOG(2) << "[CRB::orthonormalize] starting ...\n";

    for ( size_type i = 0; i < N; ++i )
    {
        for ( size_type j = std::max( i+1,N-Nm ); j < N; ++j )
        {
            value_type __rij_pr = M_model->scalarProduct(  wn[i], wn[ j ] );
            wn[j].add( -__rij_pr, wn[i] );
        }
    }

    // normalize
    for ( size_type i =N-Nm; i < N; ++i )
    {
        value_type __rii_pr = math::sqrt( M_model->scalarProduct(  wn[i], wn[i] ) );
        wn[i].scale( 1./__rii_pr );
    }

    DVLOG(2) << "[CRB::orthonormalize] finished ...\n";
    DVLOG(2) << "[CRB::orthonormalize] copying back results in basis\n";

    if ( this->vm()["crb.check.gs"].template as<int>() )
        checkOrthonormality( N , wn );

}

template <typename TruthModelType>
void
CRBTrilinear<TruthModelType>::checkOrthonormality ( int N, const wn_type& wn ) const
{

    if ( wn.size()==0 )
    {
        throw std::logic_error( "[CRB::checkOrthonormality] ERROR : size of wn is zero" );
    }

    if ( orthonormalize_primal*orthonormalize_dual==0 && this->worldComm().globalRank() == this->worldComm().masterRank() )
    {
        std::cout<<"Warning : calling checkOrthonormality is called but ";
        std::cout<<" orthonormalize_dual = "<<orthonormalize_dual;
        std::cout<<" and orthonormalize_primal = "<<orthonormalize_primal<<std::endl;
    }

    matrixN_type A, I;
    A.setZero( N, N );
    I.setIdentity( N, N );

    for ( int i = 0; i < N; ++i )
    {
        for ( int j = 0; j < N; ++j )
        {
            A( i, j ) = M_model->scalarProduct(  wn[i], wn[j] );
        }
    }

    A -= I;
    DVLOG(2) << "orthonormalization: " << A.norm() << "\n";
    if( this->worldComm().globalRank() == this->worldComm().masterRank() )
        std::cout << "    o check : " << A.norm() << " (should be 0)\n";
    //FEELPP_ASSERT( A.norm() < 1e-14 )( A.norm() ).error( "orthonormalization failed.");
}

template <typename TruthModelType>
void
CRBTrilinear<TruthModelType>::exportBasisFunctions( const export_vector_wn_type& export_vector_wn )const
{


    std::vector<wn_type> vect_wn=export_vector_wn.template get<0>();
    std::vector<std::string> vect_names=export_vector_wn.template get<1>();

    if ( vect_wn.size()==0 )
    {
        throw std::logic_error( "[CRB::exportBasisFunctions] ERROR : there are no wn_type to export" );
    }


    auto first_wn = vect_wn[0];
    auto first_element = first_wn[0];

    exporter->step( 0 )->setMesh( first_element.functionSpace()->mesh() );
    int basis_number=0;
    BOOST_FOREACH( auto wn , vect_wn )
    {

        if ( wn.size()==0 )
        {
            throw std::logic_error( "[CRB::exportBasisFunctions] ERROR : there are no element to export" );
        }

        int element_number=0;
        parameter_type mu;

        BOOST_FOREACH( auto element, wn )
        {

            std::string basis_name = vect_names[basis_number];
            std::string number = ( boost::format( "%1%_with_parameters" ) %element_number ).str() ;
            mu = M_WNmu->at( element_number );
            std::string mu_str;

            for ( int i=0; i<mu.size(); i++ )
            {
                mu_str= mu_str + ( boost::format( "_%1%" ) %mu[i] ).str() ;
            }

            std::string name =   basis_name + number + mu_str;
            exporter->step( 0 )->add( name, element );
            element_number++;
        }
        basis_number++;
    }

    exporter->save();

}


template<typename TruthModelType>
void
CRBTrilinear<TruthModelType>::printMuSelection( void )
{

    LOG(INFO)<<" List of parameter selectionned during the offline algorithm \n";
    for(int k=0;k<M_WNmu->size();k++)
    {
        LOG(INFO)<<" mu "<<k<<" = [ ";
        parameter_type const& _mu = M_WNmu->at( k );
        for( int i=0; i<_mu.size()-1; i++ )
            LOG(INFO)<<_mu(i)<<" , ";
        LOG(INFO)<<_mu( _mu.size()-1 )<<" ] \n";
    }

    if( this->worldComm().globalRank() == this->worldComm().masterRank() )
    {
        std::cout<<" List of parameter selectionned during the offline algorithm"<<std::endl;
        for(int k=0;k<M_WNmu->size();k++)
        {
            std::cout<<" mu "<<k<<" = [ ";
            parameter_type const& _mu = M_WNmu->at( k );
            for( int i=0; i<_mu.size()-1; i++ )
                std::cout<<_mu(i)<<" , ";
            std::cout<<_mu( _mu.size()-1 )<<" ] "<<std::endl;
        }
    }

}

template<typename TruthModelType>
boost::tuple<double,double,double,double,typename CRBTrilinear<TruthModelType>::vectorN_type >
CRBTrilinear<TruthModelType>::run( parameter_type const& mu, double eps , int N)
{

    //int Nwn = M_N;
    int Nwn_max = vm()["crb.dimension-max"].template as<int>();

    vectorN_type uN;

    int Nwn;
    if( N > 0 )
    {
        //in this case we want to fix Nwn
        Nwn = N;
    }
    else
    {
        //Nwn = Nwn_max;
        //M_N may be different of dimension-max
        Nwn = M_N;
    }
    auto o = lb( Nwn, mu, uN );
    double output = o.template get<0>();

    double e = 0;
    double condition_number = o.template get<1>();

    return boost::make_tuple( output , e, Nwn , condition_number, uN);
}


template<typename TruthModelType>
typename CRBTrilinear<TruthModelType>::element_type
CRBTrilinear<TruthModelType>::expansion( parameter_type const& mu , int N)
{
    int Nwn;

    if( N > 0 )
        Nwn = N;
    else
        Nwn = M_N;

    vectorN_type uN;

    auto o = lb( Nwn, mu, uN );
    int size = uN.size();
    return Feel::expansion( M_WN, uN , Nwn);
}


template<typename TruthModelType>
typename CRBTrilinear<TruthModelType>::element_type
CRBTrilinear<TruthModelType>::expansion( vectorN_type const& u , int const N) const
{
    int Nwn;

    if( N > 0 )
        Nwn = N;
    else
        Nwn = M_N;

    //FEELPP_ASSERT( N == u.size() )( N )( u.size() ).error( "invalid expansion size");
    FEELPP_ASSERT( Nwn == u.size() )( Nwn )( u.size() ).error( "invalid expansion size");
    return Feel::expansion( M_WN, u, N );
}



template<typename TruthModelType>
template<class Archive>
void
CRBTrilinear<TruthModelType>::save( Archive & ar, const unsigned int version ) const
{
    int proc_number = this->worldComm().globalRank();

    LOG(INFO) <<"[CRBTrilinear::save] version : "<<version<<std::endl;

    auto mesh = mesh_type::New();
    auto is_mesh_loaded = mesh->load( _name="mymesh",_path=this->dbLocalPath(),_type="binary" );

    if ( ! is_mesh_loaded )
    {
        auto first_element = M_WN[0];
        mesh = first_element.functionSpace()->mesh() ;
        mesh->save( _name="mymesh",_path=this->dbLocalPath(),_type="binary" );
    }

    auto Xh = space_type::New( mesh );

    ar & boost::serialization::base_object<super>( *this );
    ar & BOOST_SERIALIZATION_NVP( M_output_index );
    ar & BOOST_SERIALIZATION_NVP( M_N );
    ar & BOOST_SERIALIZATION_NVP( M_rbconv );
    ar & BOOST_SERIALIZATION_NVP( M_error_type );
    ar & BOOST_SERIALIZATION_NVP( M_Xi );
    ar & BOOST_SERIALIZATION_NVP( M_WNmu );
    ar & BOOST_SERIALIZATION_NVP( M_Aqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Aqm_tril_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Fqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Lqm_pr );

    ar & BOOST_SERIALIZATION_NVP( M_current_mu );
    ar & BOOST_SERIALIZATION_NVP( M_no_residual_index );

    for(int i=0; i<M_N; i++)
        ar & BOOST_SERIALIZATION_NVP( M_WN[i] );

    ar & BOOST_SERIALIZATION_NVP( M_maxerror );
}

template<typename TruthModelType>
template<class Archive>
void
CRBTrilinear<TruthModelType>::load( Archive & ar, const unsigned int version )
{

    int proc_number = this->worldComm().globalRank();

    LOG(INFO) <<"[CRBTrilinear::load] version"<< version <<std::endl;

    mesh_ptrtype mesh;
    space_ptrtype Xh;

    if ( !M_model )
    {
        LOG(INFO) << "[load] model not initialized, loading fdb files...\n";
        mesh = mesh_type::New();

        bool is_mesh_loaded = mesh->load( _name="mymesh",_path=this->dbLocalPath(),_type="binary" );

        Xh = space_type::New( mesh );
        LOG(INFO) << "[load] loading fdb files done.\n";
    }
    else
    {
        LOG(INFO) << "[load] get mesh/Xh from model...\n";
        mesh = M_model->functionSpace()->mesh();
        Xh = M_model->functionSpace();
        LOG(INFO) << "[load] get mesh/Xh from model done.\n";
    }

    ar & boost::serialization::base_object<super>( *this );
    ar & BOOST_SERIALIZATION_NVP( M_output_index );
    ar & BOOST_SERIALIZATION_NVP( M_N );

	ar & BOOST_SERIALIZATION_NVP( M_rbconv );

    ar & BOOST_SERIALIZATION_NVP( M_error_type );
    ar & BOOST_SERIALIZATION_NVP( M_Xi );
    ar & BOOST_SERIALIZATION_NVP( M_WNmu );
    ar & BOOST_SERIALIZATION_NVP( M_Aqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Aqm_tril_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Fqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Lqm_pr );


    ar & BOOST_SERIALIZATION_NVP( M_current_mu );
    ar & BOOST_SERIALIZATION_NVP( M_no_residual_index );

    element_type temp = Xh->element();

    M_WN.resize( M_N );

    for( int i = 0 ; i < M_N ; i++ )
    {
        temp.setName( (boost::format( "fem-primal-%1%" ) % ( i ) ).str() );
        ar & BOOST_SERIALIZATION_NVP( temp );
        M_WN[i] = temp;
    }

    ar & BOOST_SERIALIZATION_NVP( M_maxerror );

    LOG(INFO) << "[CRBTrilinear::load] end of load function" << std::endl;
}


template<typename TruthModelType>
void
CRBTrilinear<TruthModelType>::saveDB()
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
CRBTrilinear<TruthModelType>::loadDB()
{
    if ( this->rebuildDB() )
        return false;

    fs::path db = this->lookForDB();

    if ( db.empty() )
        return false;

    if ( !fs::exists( db ) )
        return false;

    //std::cout << "Loading " << db << "...\n";
    fs::ifstream ifs( db );

    if ( ifs )
    {
        boost::archive::text_iarchive ia( ifs );
        // write class instance to archive
        ia >> *this;
        //std::cout << "Loading " << db << " done...\n";
        this->setIsLoaded( true );
        // archive and stream closed when destructors are called
        return true;
    }

    return false;
}


} // Feel

#endif /* __CRBTrilinear_H */
