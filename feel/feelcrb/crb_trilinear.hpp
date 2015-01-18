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
//class CRBTrilinear : public CRBDB,
//                     public CRB<TruthModelType>
template<typename TruthModelType>
class CRBTrilinear : public CRB<TruthModelType>
{
    //typedef  CRBDB super_crbdb;
    typedef  CRB<TruthModelType> super_crb;
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

    typedef boost::tuple< std::vector<vectorN_type> , std::vector<vectorN_type> , std::vector<vectorN_type>, std::vector<vectorN_type> > solutions_tuple;
    typedef boost::tuple< double,double,double , std::vector< std::vector< double > > , std::vector< std::vector< double > > > upper_bounds_tuple;
    typedef boost::tuple< double,double > matrix_info_tuple; //conditioning, determinant

    // ! export
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    typedef CRBTrilinear self_type;

    //! scm
    typedef CRBSCM<truth_model_type> scm_type;
    typedef boost::shared_ptr<scm_type> scm_ptrtype;

    //! elements database
    typedef CRBElementsDB<truth_model_type> crb_elements_db_type;
    typedef boost::shared_ptr<crb_elements_db_type> crb_elements_db_ptrtype;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    CRBTrilinear()
        :
        super_crb(),
        M_crbdb()
    {

    }

    //! constructor from command line options
    CRBTrilinear( std::string  name,
         po::variables_map const& vm,
         truth_model_ptrtype const & model )
        :
        super_crb( name, vm , model ),
        M_crbdb( ( boost::format( "%1%" ) % ioption(_name="crb.error-type") ).str(),
                 name,
                 ( boost::format( "%1%-%2%-%3%-trilinear" ) % name % iioption(_name="crb.output-index") % option(_name="crb.error-type") ).str(),
                 vm )
    {
        this->setTruthModel( model );
        if ( M_crbdb.loadDB() )
            LOG(INFO) << "Database " << M_crbdb.lookForDB() << " available and loaded\n";

        //this will be in the offline step (it's only when we enrich or create the database that we want to have access to elements of the RB)
        this->M_elements_database.setMN( this->M_N );
        if( this->M_elements_database.loadDB() )
        {
            LOG(INFO) << "database for basis functions " << this->M_elements_database.lookForDB() << " available and loaded\n";
            auto basis_functions = this->M_elements_database.wn();
            this->M_model->rBFunctionSpace()->setBasis( basis_functions );
        }
        else
        {
            LOG( INFO ) <<"no database for basis functions loaded. Start from the begining";
        }

    }


    //! copy constructor
    CRBTrilinear( CRBTrilinear const & o )
        :
        super_crb( o ),
        M_crbdb( o )
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


    WorldComm const& worldComm() const { return Environment::worldComm() ; }

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
    boost::tuple<std::vector<double>,matrix_info_tuple> lb( size_type N, parameter_type const& mu, std::vector< vectorN_type >& uN, std::vector< vectorN_type >& uNdu ,
                                               std::vector<vectorN_type> & uNold, std::vector<vectorN_type> & uNduold, bool print_rb_matrix=false, int K=0 ) const;

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

    /**
     * save the CRB database
     */
    void saveDB();

    /**
     * load the CRB database
     */
    bool loadDB();

    //@}


private:
    CRBDB M_crbdb;

    //crb_elements_db_type M_elements_database;

    std::vector < std::vector < matrixN_type> >  M_Aqm_tril_pr;
    mutable matrixN_type M_bilinear_terms;
    mutable vectorN_type M_linear_terms;
    //boost::shared_ptr<SolverNonLinear<double> > M_nlsolver;


    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void save( Archive & ar, const unsigned int version ) const;

    template<class Archive>
    void load( Archive & ar, const unsigned int version ) ;

    BOOST_SERIALIZATION_SPLIT_MEMBER()


};

template<typename TruthModelType>
typename CRBTrilinear<TruthModelType>::convergence_type
CRBTrilinear<TruthModelType>::offline()
{

    int proc_number = this->worldComm().globalRank();

    bool rebuild_database = boption(_name="crb.rebuild-database") ;
    bool orthonormalize_primal = boption(_name="crb.orthonormalize-primal") ;

    boost::timer ti;
    if( this->worldComm().isMasterRank() )
        std::cout << "Offline CRBTrilinear starts, this may take a while until Database is computed..."<<std::endl;
    LOG(INFO) << "[CRBTrilinear::offline] Starting offline for output " << this->M_output_index << "\n";
    LOG(INFO) << "[CRBTrilinear::offline] initialize underlying finite element model\n";
    //M_model->initModel();
    LOG( INFO )<< " -- model init done in " << ti.elapsed() << "s";

    parameter_type mu( this->M_Dmu );

    double delta_pr;
    double delta_du;
    size_type index;
    //if M_N == 0 then there is not an already existing database
    if ( rebuild_database || this->M_N == 0)
    {

        ti.restart();

        LOG(INFO) << "[CRBTrilinear::offline] compute random sampling\n";

        int total_proc = this->worldComm().globalSize();
        std::string sampling_mode = soption("crb.sampling-mode");
        bool all_proc_same_sampling=boption("crb.all-procs-have-same-sampling");
        int sampling_size = ioption("crb.sampling-size");
        std::string file_name = ( boost::format("M_Xi_%1%_"+sampling_mode+"-proc%2%on%3%") % sampling_size %proc_number %total_proc ).str();
        if( all_proc_same_sampling )
            file_name+="-all-proc-have-same-sampling";

        std::ifstream file ( file_name );

        if( ! file )
        {

            // random sampling
            std::string supersamplingname =(boost::format("Dmu-%1%-generated-by-master-proc") %sampling_size ).str();

            if( sampling_mode == "log-random" )
                this->M_Xi->randomize( sampling_size , all_proc_same_sampling , supersamplingname );
            else if( sampling_mode == "log-equidistribute" )
                this->M_Xi->logEquidistribute( sampling_size , all_proc_same_sampling , supersamplingname );
            else if( sampling_mode == "equidistribute" )
                this->M_Xi->equidistribute( sampling_size , all_proc_same_sampling , supersamplingname );
            else
                throw std::logic_error( "[CRBTrilinear::offline] ERROR invalid option crb.sampling-mode, please select between log-random, log-equidistribute or equidistribute" );
            //M_Xi->equidistribute( this->vm()["crb.sampling-size"].template as<int>() );

            this->M_Xi->writeOnFile(file_name);

        }
        else
        {
            this->M_Xi->clear();
            this->M_Xi->readFromFile(file_name);
        }

        this->M_WNmu->setSuperSampling( this->M_Xi );


        LOG( INFO )<<"[CRBTrilinear offline] M_error_type = "<<this->M_error_type<<std::endl;

        LOG(INFO) << " -- sampling init done in " << ti.elapsed() << "s";
        ti.restart();

        // empty sets
        this->M_WNmu->clear();
        if( this->M_error_type == CRB_NO_RESIDUAL )
            mu = this->M_Dmu->element();
        else
        {
            // start with M_C = { arg min mu, mu \in Xi }
            boost::tie( mu, index ) = this->M_Xi->min();
        }


        int size = mu.size();
        //std::cout << " -- WN size :  " << M_WNmu->size() << "\n";

        // dimension of reduced basis space
        this->M_N = 0;

        this->M_maxerror = 1e10;
        delta_pr = 0;
        delta_du = 0;
        //boost::tie( M_maxerror, mu, index ) = maxErrorBounds( N );

        LOG(INFO) << "[CRBTrilinear::offline] allocate reduced basis data structures\n";

        this->M_Aqm_pr.resize( this->M_model->Qa() );
        for(int q=0; q<this->M_model->Qa(); q++)
        {
            this->M_Aqm_pr[q].resize( 1 );
        }

        M_Aqm_tril_pr.resize( this->M_model->QaTri() );

        //for(int q=0; q<this->M_model->QaTri(); q++)
        this->M_Fqm_pr.resize( this->M_model->Ql( 0 ) );

        for(int q=0; q<this->M_model->Ql( 0 ); q++)
        {
            this->M_Fqm_pr[q].resize( 1 );
        }

        this->M_Lqm_pr.resize( this->M_model->Ql( this->M_output_index ) );
        for(int q=0; q<this->M_model->Ql( this->M_output_index ); q++)
            this->M_Lqm_pr[q].resize( 1 );

    }//end of if( rebuild_database )
#if 1
    else
    {
        mu = this->M_current_mu;
        if( proc_number == 0 )
        {
            std::cout<<"we are going to enrich the reduced basis"<<std::endl;
            std::cout<<"there are "<<this->M_N<<" elements in the database"<<std::endl;
        }
        LOG(INFO) <<"we are going to enrich the reduced basis"<<std::endl;
        LOG(INFO) <<"there are "<<this->M_N<<" elements in the database"<<std::endl;
    }//end of else associated to if ( rebuild_databse )
#endif

    LOG(INFO) << "[CRBTrilinear::offline] compute affine decomposition\n";

    std::vector< std::vector<sparse_matrix_ptrtype> > Aqm;
    std::vector< std::vector<sparse_matrix_ptrtype> > Aqm_tril;
    std::vector< std::vector<std::vector<vector_ptrtype> > > Fqm;

    boost::tie( boost::tuples::ignore, Aqm, Fqm ) = this->M_model->computeAffineDecomposition();

    element_ptrtype u( new element_type( this->M_model->functionSpace() ) );

    LOG(INFO) << "[CRBTrilinear::offline] starting offline adaptive loop\n";

    bool reuse_prec = this->vm()["crb.reuse-prec"].template as<bool>() ;

    bool use_predefined_WNmu = this->vm()["crb.use-predefined-WNmu"].template as<bool>() ;
    int N_log_equi = this->vm()["crb.use-logEquidistributed-WNmu"].template as<int>() ;
    int N_equi = this->vm()["crb.use-equidistributed-WNmu"].template as<int>() ;
    int N_random = ioption( "crb.use-random-WNmu" );

    /*    if( N_log_equi > 0 || N_equi > 0 )
     use_predefined_WNmu = true;*/

    // file where the sampling is savec
    std::string file_name = ( boost::format("SamplingWNmu") ).str();
    std::ifstream file ( file_name );

    this->M_WNmu->clear();

    if ( use_predefined_WNmu ) // In this case we want to read the sampling
    {
        if( ! file ) // The user forgot to give the sampling file
            throw std::logic_error( "[CRB::offline] ERROR the file SamplingWNmu doesn't exist so it's impossible to known which parameters you want to use to build the database" );
        else
        {
            int sampling_size = this->M_WNmu->readFromFile(file_name);
            if( Environment::isMasterRank() )
                std::cout<<"[CRB::offline] Read WNmu ( sampling size : "
                         << sampling_size <<" )"<<std::endl;
            LOG( INFO )<<"[CRB::offline] Read WNmu ( sampling size : "
                       << sampling_size <<" )";
        }
    }
    else // We generate the sampling with choosen strategy
    {
        if ( N_log_equi>0 )
        {
            this->M_WNmu->logEquidistribute( N_log_equi , true );
            if( Environment::isMasterRank() )
                std::cout<<"[CRB::offline] Log-Equidistribute WNmu ( sampling size : "
                         <<N_log_equi<<" )"<<std::endl;
            LOG( INFO )<<"[CRB::offline] Log-Equidistribute WNmu ( sampling size : "
                       <<N_log_equi<<" )";
        }
        else if ( N_equi>0 )
        {
            this->M_WNmu->equidistribute( N_equi , true );
            if( Environment::isMasterRank() )
                std::cout<<"[CRB::offline] Equidistribute WNmu ( sampling size : "
                         <<N_equi<<" )"<<std::endl;
            LOG( INFO )<<"[CRB::offline] Equidistribute WNmu ( sampling size : "
                       <<N_equi<<" )";
        }
        else if ( N_random>0 )
        {
            this->M_WNmu->randomize( N_random , true );
            if( Environment::isMasterRank() )
                std::cout<<"[CRB::offline] Randomize WNmu ( sampling size : "
                         <<N_random<<" )"<<std::endl;
            LOG( INFO )<<"[CRB::offline] Randomize WNmu ( sampling size : "
                       <<N_random<<" )";
        }
        else // In this case we don't know what sampling to use
            throw std::logic_error( "[CRB::offline] ERROR : You have to choose an appropriate strategy for the offline sampling : random, equi, logequi or predefined" );

        this->M_WNmu->writeOnFile(file_name);

        /*        if( ! file )
        {
            this->M_WNmu->clear();
            std::vector< parameter_type > V;
            parameter_type __mu;
            __mu = this->M_Dmu->element();
            __mu(0)= 1      ; __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)= 111112 ; __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)= 222223 ; __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)= 333334 ; __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)= 444445 , __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)= 555556 ; __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)= 666667 ; __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)= 777778 ; __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)= 888889 ; __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)= 1e+06  ; __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)= 8123   ; __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)= 9123   ; __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)=1.123e4 ; __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)=2.123e4 ; __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)=4.123e4 ; __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)=912     ; __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)=1.123e3 ; __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)=4.123e3 ; __mu(1)= 1  ; V.push_back( __mu );
         __mu(0)=7.123e4 ; __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)=2123    ; __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)=6.123e3 ; __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)=3.123e3 ; __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)=3.123e4 ; __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)=5.123e4 ; __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)=9.123e4 ; __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)=812     ; __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)=5.111e3 ; __mu(1)= 1  ; V.push_back( __mu );
            __mu(0)=5.124e2 ; __mu(1)= 1  ; V.push_back( __mu );
            this->M_WNmu->setElements( V );
            this->M_iter_max = this->M_WNmu->size();
            this->M_WNmu->writeOnFile(file_name);
         }*/
        use_predefined_WNmu=true;
    } //build sampling

    this->M_iter_max = this->M_WNmu->size();
    mu = this->M_WNmu->at( this->M_N ); // first element


    if( this->M_error_type == CRB_NO_RESIDUAL || use_predefined_WNmu )
    {
        //in this case it makes no sens to check the estimated error
        this->M_maxerror = 1e10;
    }


    LOG(INFO) << "[CRBTrilinear::offline] strategy "<< this->M_error_type <<"\n";

    while ( this->M_maxerror > this->M_tolerance && this->M_N < this->M_iter_max )
    {

        boost::timer timer, timer2;
        LOG(INFO) <<"========================================"<<"\n";
        if( proc_number == this->worldComm().masterRank() )
            std::cout<<"construction of "<<this->M_N<<"/"<<this->M_iter_max<<" basis "<<std::endl;
        LOG(INFO) << "N=" << this->M_N << "/"  << this->M_iter_max << "( nb proc : "<<worldComm().globalSize()<<")";

        // for a given parameter \p mu assemble the left and right hand side
        u->setName( ( boost::format( "fem-primal-N%1%-proc%2%" ) % (this->M_N)  % proc_number ).str() );

        mu.check();
        u->zero();

        timer2.restart();

        LOG(INFO) << "[CRB::offline] solving primal" << "\n";
        *u = this->M_model->solve( mu );

        //if( proc_number == this->worldComm().masterRank() ) std::cout << "  -- primal problem solved in " << timer2.elapsed() << "s\n";
        timer2.restart();


        if( ! use_predefined_WNmu )
            this->M_WNmu->push_back( mu, index );

        this->M_WNmu_complement = this->M_WNmu->complement();

        this->M_model->rBFunctionSpace()->addPrimalBasisElement( *u );
        //WARNING : the dual element is not the real dual solution !
        //no dual problem was solved
        this->M_model->rBFunctionSpace()->addDualBasisElement( *u );

	    int number_of_added_elements=1;
        this->M_N+=number_of_added_elements;

        if ( orthonormalize_primal )
        {
            this->orthonormalize( this->M_N, this->M_model->rBFunctionSpace()->primalRB() , number_of_added_elements );
            this->orthonormalize( this->M_N, this->M_model->rBFunctionSpace()->primalRB() , number_of_added_elements );
            this->orthonormalize( this->M_N, this->M_model->rBFunctionSpace()->primalRB() , number_of_added_elements );
        }

        LOG(INFO) << "[CRB::offline] compute Aq_pr, Aq_du, Aq_pr_du" << "\n";
        for  (size_type q = 0; q < this->M_model->Qa(); ++q )
        {
            this->M_Aqm_pr[q][0].conservativeResize( this->M_N, this->M_N );

            // only compute the last line and last column of reduced matrices
            for ( size_type i = this->M_N-number_of_added_elements; i < this->M_N; i++ )
            {
                for ( size_type j = 0; j < this->M_N; ++j )
                {
                    this->M_Aqm_pr[q][0]( i, j ) = Aqm[q][0]->energy( this->M_model->rBFunctionSpace()->primalBasisElement(i) , this->M_model->rBFunctionSpace()->primalBasisElement(j) );
                }
            }

            for ( size_type j=this->M_N-number_of_added_elements; j < this->M_N; j++ )
            {
                for ( size_type i = 0; i < this->M_N; ++i )
                {
                    this->M_Aqm_pr[q][0]( i, j ) = Aqm[q][0]->energy( this->M_model->rBFunctionSpace()->primalBasisElement(i), this->M_model->rBFunctionSpace()->primalBasisElement(j) );
                }
            }
        }//loop over q


        LOG(INFO) << "[CRBTrilinear::offline] compute Fq_pr" << "\n";

        for ( size_type q = 0; q < this->M_model->Ql( 0 ); ++q )
        {
            this->M_Fqm_pr[q][0].conservativeResize( this->M_N );

            for ( size_type l = 1; l <= number_of_added_elements; ++l )
            {
                int index = this->M_N-l;
                this->M_Fqm_pr[q][0]( index ) = this->M_model->Fqm( 0, q, 0, this->M_model->rBFunctionSpace()->primalBasisElement(index) );
            }
        }//loop over q

        LOG(INFO) << "[CRB::offline] compute Lq_pr" << "\n";

        for ( size_type q = 0; q < this->M_model->Ql( this->M_output_index ); ++q )
        {
            this->M_Lqm_pr[q][0].conservativeResize( this->M_N );

            for ( size_type l = 1; l <= number_of_added_elements; ++l )
            {
                int index = this->M_N-l;
                this->M_Lqm_pr[q][0]( index ) = this->M_model->Fqm( this->M_output_index, q, 0, this->M_model->rBFunctionSpace()->primalBasisElement(index) );
            }
        }//loop over q

        sparse_matrix_ptrtype trilinear_form;
        for  (size_type q = 0; q < this->M_model->QaTri(); ++q )
        {
            M_Aqm_tril_pr[q].resize( this->M_N );
            for (int k=0 ; k<this->M_N; k++)
            {
                //bring back the matrix associated to the trilinear form for a given basis function
                //we do this here to use only one matrix
                trilinear_form  = this->M_model->computeTrilinearForm( this->M_model->rBFunctionSpace()->primalBasisElement(k) );

                M_Aqm_tril_pr[q][k].conservativeResize( this->M_N, this->M_N );
                for ( int i = 0; i < this->M_N; ++i )
                {
                    for ( int j = 0; j < this->M_N; ++j )
                    {
                        M_Aqm_tril_pr[q][k]( i, j ) = trilinear_form->energy( this->M_model->rBFunctionSpace()->primalBasisElement(j), this->M_model->rBFunctionSpace()->primalBasisElement(i) );
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
                mu = this->M_Dmu->element();
                //make sure that the new mu is not already is M_WNmu
                BOOST_FOREACH( auto _mu, *this->M_WNmu )
                {
                    if( mu == _mu )
                        already_exist=true;
                }
            }
            while( already_exist );
            this->M_current_mu = mu;
        }
        else
        {
            //remmber that in this case M_iter_max = sampling size
            if( this->M_N < this->M_iter_max )
            {
                mu = this->M_WNmu->at( this->M_N );
                this->M_current_mu = mu;
            }
        }

        this->M_rbconv.insert( convergence( this->M_N, boost::make_tuple(this->M_maxerror,delta_pr,delta_du) ) );

        timer2.restart();
        LOG(INFO) << "time: " << timer.elapsed() << "\n";
        LOG(INFO) <<"========================================"<<"\n";

        //save DB after adding an element
        saveDB();
        this->M_elements_database.setWn( boost::make_tuple( this->M_model->rBFunctionSpace()->primalRB() , this->M_model->rBFunctionSpace()->dualRB() ) );
        this->M_elements_database.saveDB();
    }


    LOG( INFO )<<"number of elements in the reduced basis : "<<this->M_N<<" ( nb proc : "<<worldComm().globalSize()<<")";
    bool visualize_basis = this->vm()["crb.visualize-basis"].template as<bool>() ;

    if ( visualize_basis )
    {
        std::vector<wn_type> wn;
        std::vector<std::string> names;
        wn.push_back( this->M_model->rBFunctionSpace()->primalRB() );
        names.push_back( "primal" );
        this->exportBasisFunctions( boost::make_tuple( wn ,names ) );

        if ( orthonormalize_primal )
            LOG(INFO)<<"[CRB::offline] Basis functions have been exported but warning elements have been orthonormalized";
    }

    LOG( INFO ) << "Offline CRB is done";

    return this->M_rbconv;

}


template<typename TruthModelType>
void
CRBTrilinear<TruthModelType>::findNearestNeighborInWNmu( parameter_type const& mu,
                                                         parameter_type & neighbor,
                                                         int & index ) const
{
    std::vector<int> index_vector;
    sampling_ptrtype S =  this->M_WNmu->searchNearestNeighbors( mu, 1 , index_vector);
    neighbor = S->at( 0 );
    index = index_vector[0];
    //std::cout<<"[CRBTrilinear::findNearestNeighborInWNmu] for Gr = "<<mu(0)<<" th nearest neighbor in WNmu is "<<neighbor(0)<<" at index "<<index<<std::endl;
}

template<typename TruthModelType>
typename boost::tuple<std::vector<double>,typename CRBTrilinear<TruthModelType>::matrix_info_tuple >
CRBTrilinear<TruthModelType>::lb( size_type N, parameter_type const& mu,
                                  std::vector< vectorN_type >& uN,
                                  std::vector< vectorN_type >& uNdu,
                                  std::vector<vectorN_type> & uNold,
                                  std::vector<vectorN_type> & uNduold,
                                  bool print_rb_matrix, int K ) const
{
    uN.resize(1);
    if ( N > this->M_N ) N =this->M_N;
    beta_vector_type betaAqm;
    std::vector<beta_vector_type> betaFqm, betaLqm;

    //-- end of initialization step

    std::vector<double>output_vector(1);

    boost::tie( boost::tuples::ignore,  betaAqm, betaFqm  ) = this->M_model->computeBetaQm( mu );

    /*
     ------> Here add Newton Method and continuous Gr,Pr
     */
    using namespace vf;
    //Feel::ParameterSpace<2>::Element current_mu( mu );
    parameter_type current_mu ;

    //backend_ptrtype backend_primal_problem = backend_type::build( BACKEND_PETSC );

    vectorN_type R( (int) N );
    matrixN_type J( (int) N , (int) N );

    //initialization of uN
    //we look for the nearest neighbor of mu in the sampling WNmu
    //let i the index of this neighbor in WNmu, we will set zeros in uN except at the i^th component where we will set 1
    uN[0].setZero( (int) N );
    int number_of_neighbors = this->M_N - N + 1;
    std::vector<int> index_vector;
    sampling_ptrtype S = this->M_WNmu->searchNearestNeighbors( mu, number_of_neighbors, index_vector);
    int n_index=0;
    int index;
    //with this loop we check that the index of the nearest neighbor is not out
    //of range. This case only happens when you do not want to use all the reduced basis :
    // online Wn size < offline Wn size
    do
    {
        index=index_vector[n_index];
        n_index++;
    }while( index>=N );

    //parameter_type neighbor( this->M_Dmu );
    //int index;
    //    findNearestNeighborInWNmu(  mu,  neighbor, index );

     if( this->vm()["crb.cvg-study"].template as<bool>() == true )
    {
        //in this case, index may be smaller than uN.size
        //so we do nothing
    }
    else
        uN[0]( index ) = 1;
    double *r_data = R.data();
    double *j_data = J.data();
    double *uN_data = uN[0].data();

    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic , 1> > map_R ( r_data, N );
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic , 1> > map_uN ( uN_data, N );
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic , Eigen::Dynamic> > map_J ( j_data, N , N );

    /*double gr = mu( 0 );
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
     //LOG( INFO ) << "[CRBTrilinear::lb] i/N = " << i+1 << "/" << Nmax ;
     //LOG( INFO ) << "[CRBTrilinear::lb] intermediary Grashof = " << current_Grashofs;
     //LOG( INFO ) << "[CRBTrilinear::lb] and Prandtl = " << current_Prandtl ;
     }
     else
     {
     current_Grashofs = gr;
     current_Prandtl = pr;
     }*/

    //current_mu << current_Grashofs, current_Prandtl;
    current_mu = mu;

    this->updateLinearTerms( current_mu , N );

    //this->M_nlsolver->setRelativeResidualTol( 1e-12 );
    this->M_nlsolver->map_dense_jacobian = boost::bind( &self_type::updateJacobian, boost::ref( *this ), _1, _2  , current_mu , N );
    this->M_nlsolver->map_dense_residual = boost::bind( &self_type::updateResidual, boost::ref( *this ), _1, _2  , current_mu , N );
    this->M_nlsolver->setType( TRUST_REGION );
    this->M_nlsolver->solve( map_J , map_uN , map_R, 1e-12, 100);
    //}

    LOG(INFO) << "[CRBTrilinear::lb] solve with Newton done";

    double condition_number = 0;
    double determinant=0;
    LOG( INFO ) <<"[CRBTrilinear::lb] The condition number of jacobian done is not computed ";
    vectorN_type L ( ( int )N );
    L.setZero( N );

    int output_index=this->M_output_index;
    int qoutput = this->M_model->Ql( output_index );

    for ( size_type q = 0; q < qoutput; ++q )
    {
        L += betaFqm[output_index][q][0] * this->M_Lqm_pr[q][0].head( N );
    }
    output_vector[0] = L.dot( uN[0] );
    LOG(INFO) << "[CRBTrilinear::lb] computation of the output done";

    //std::cout<<"[CRBTrilinear uN] : \n"<<uN<<std::endl;

    auto matrix_info = boost::make_tuple( condition_number, determinant );
    return boost::make_tuple( output_vector, matrix_info);

}
template<typename TruthModelType>
void
CRBTrilinear<TruthModelType>::updateLinearTerms( parameter_type const& mu , int N ) const
{

    LOG(INFO) << "update linear terms \n";

    beta_vector_type betaAqm;
    std::vector<beta_vector_type> betaFqm, betaLqm;

    boost::tie( boost::tuples::ignore,  betaAqm, betaFqm ) = this->M_model->computeBetaQm( mu );

    M_bilinear_terms.setZero( N , N );

    for ( size_type q = 0; q < this->M_model->Qa(); ++q )
    {
        M_bilinear_terms += betaAqm[q][0]*this->M_Aqm_pr[q][0].block( 0, 0, N, N );
    }

    M_linear_terms.setZero( N );

    for ( size_type q = 0; q < this->M_model->Ql( 0 ); ++q )
    {
        M_linear_terms += betaFqm[0][q][0]*this->M_Fqm_pr[q][0].head( N );
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
        for ( size_type q = 0; q < this->M_model->QaTri(); ++q )
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
        auto expansionX = this->expansion( map_X , N , this->M_model->rBFunctionSpace()->primalRB() );
        auto J = this->M_model->jacobian( expansionX );
        matrixN_type model_reduced_jacobian( N , N );
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
                model_reduced_jacobian(i,j) = J->energy( this->M_model->rBFunctionSpace()->primalBasisElement(i), this->M_model->rBFunctionSpace()->primalBasisElement(j) );
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

    map_R = -M_linear_terms;
    map_R += M_bilinear_terms * map_X ;

    bool enable = this->vm()["crb.enable-convection-terms"].template as<bool>();
    if( enable )
    {
        int qatri = this->M_model->QaTri();
        for ( size_type q = 0; q < qatri; ++q )
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
        auto expansionX = this->expansion( map_X , N , this->M_model->rBFunctionSpace()->primalRB() );
        auto R = this->M_model->residual( expansionX );
        vectorN_type model_reduced_residual( N );
        element_ptrtype eltR( new element_type( this->M_model->functionSpace() ) );
        for(int i=0; i<eltR->localSize();i++)
            eltR->operator()(i)=R->operator()(i);
        for(int i=0; i<N; i++)
            model_reduced_residual(i) = inner_product( *eltR , this->M_model->rBFunctionSpace()->primalBasisElement(i) );
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



template<typename TruthModelType>
template<class Archive>
void
CRBTrilinear<TruthModelType>::save( Archive & ar, const unsigned int version ) const
{
    int proc_number = this->worldComm().globalRank();

    LOG(INFO) <<"[CRBTrilinear::save] version : "<<version<<std::endl;

    ar & BOOST_SERIALIZATION_NVP( this->M_output_index );
    ar & BOOST_SERIALIZATION_NVP( this->M_N );
    ar & BOOST_SERIALIZATION_NVP( this->M_rbconv );
    ar & BOOST_SERIALIZATION_NVP( this->M_error_type );
    ar & BOOST_SERIALIZATION_NVP( this->M_Xi );
    ar & BOOST_SERIALIZATION_NVP( this->M_WNmu );
    ar & BOOST_SERIALIZATION_NVP( this->M_Aqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Aqm_tril_pr );
    ar & BOOST_SERIALIZATION_NVP( this->M_Fqm_pr );
    ar & BOOST_SERIALIZATION_NVP( this->M_Lqm_pr );

    ar & BOOST_SERIALIZATION_NVP( this->M_current_mu );
    ar & BOOST_SERIALIZATION_NVP( this->M_no_residual_index );

    ar & BOOST_SERIALIZATION_NVP( this->M_maxerror );

}

template<typename TruthModelType>
template<class Archive>
void
CRBTrilinear<TruthModelType>::load( Archive & ar, const unsigned int version )
{

    int proc_number = this->worldComm().globalRank();

    LOG(INFO) <<"[CRBTrilinear::load] version"<< version <<std::endl;

    ar & BOOST_SERIALIZATION_NVP( this->M_output_index );
    ar & BOOST_SERIALIZATION_NVP( this->M_N );

	ar & BOOST_SERIALIZATION_NVP( this->M_rbconv );

    ar & BOOST_SERIALIZATION_NVP( this->M_error_type );
    ar & BOOST_SERIALIZATION_NVP( this->M_Xi );
    ar & BOOST_SERIALIZATION_NVP( this->M_WNmu );
    ar & BOOST_SERIALIZATION_NVP( this->M_Aqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Aqm_tril_pr );
    ar & BOOST_SERIALIZATION_NVP( this->M_Fqm_pr );
    ar & BOOST_SERIALIZATION_NVP( this->M_Lqm_pr );


    ar & BOOST_SERIALIZATION_NVP( this->M_current_mu );
    ar & BOOST_SERIALIZATION_NVP( this->M_no_residual_index );

    ar & BOOST_SERIALIZATION_NVP( this->M_maxerror );
}


template<typename TruthModelType>
void
CRBTrilinear<TruthModelType>::saveDB()
{
    super_crb::saveDB();

    fs::ofstream ofs( M_crbdb.dbLocalPath() / M_crbdb.dbFilename() );

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

    super_crb::loadDB();

    fs::path db = M_crbdb.lookForDB();

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
struct version< Feel::CRBTrilinear<T> >
{
    // at the moment the version of the CRBTrilinear DB is 0. if any changes is done
    // to the format it is mandatory to increase the version number below
    // and use the new version number of identify the new entries in the DB
    typedef mpl::int_<0> type;
    typedef mpl::integral_c_tag tag;
    static const unsigned int value = version::type::value;
};
template<typename T> const unsigned int version<Feel::CRBTrilinear<T> >::value;
}
}

#endif /* __CRBTrilinear_H */
