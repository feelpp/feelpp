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
        super(),
        M_elements_database(),
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
#if 0
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
#endif
    //! constructor from command line options
    CRBTrilinear( std::string  name,
         po::variables_map const& vm,
         truth_model_ptrtype const & model )
        :
        super( ( boost::format( "%1%" ) % vm["crb.error-type"].template as<int>() ).str(),
               name,
               ( boost::format( "%1%-%2%-%3%" ) % name % vm["crb.output-index"].template as<int>() % vm["crb.error-type"].template as<int>() ).str(),
               vm ),
        M_elements_database(
                            ( boost::format( "%1%" ) % vm["crb.error-type"].template as<int>() ).str(),
                            name,
                            ( boost::format( "%1%-%2%-%3%-elements" ) % name % vm["crb.output-index"].template as<int>() % vm["crb.error-type"].template as<int>() ).str(),
                            vm ,
                            model ),
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

        //this will be in the offline step (it's only when we enrich or create the database that we want to have access to elements of the RB)
        M_elements_database.setMN( M_N );
        if( M_elements_database.loadDB() )
        {
            LOG(INFO) << "database for basis functions " << M_elements_database.lookForDB() << " available and loaded\n";
            auto basis_functions = M_elements_database.wn();
            M_model->rBFunctionSpace()->setBasis( basis_functions );
        }
        else
        {
            LOG( INFO ) <<"no database for basis functions loaded. Start from the begining";
        }

    }


    //! copy constructor
    CRBTrilinear( CRBTrilinear const & o )
        :
        super( o ),
        M_elements_database( o.M_elements_database ),
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

    void computationalTimeStatistics( std::string );

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

    void checkResidual( parameter_type const& mu, std::vector< std::vector<double> > const& primal_residual_coeffs,
                        std::vector< std::vector<double> > const& dual_residual_coeffs , element_type & u, element_type & udu ) const
    {
        LOG( INFO ) << "CRBTrilinear::checkResidual do nothing";
    }
    void compareResidualsForTransientProblems(int N, parameter_type const& mu, std::vector<element_type> const & Un,
                                              std::vector<element_type> const & Unold, std::vector<element_type> const& Undu,
                                              std::vector<element_type> const & Unduold,
                                              std::vector< std::vector<double> > const& primal_residual_coeffs,
                                              std::vector< std::vector<double> > const& dual_residual_coeffs  ) const
    {
        LOG( INFO ) << "CRBTrilinear::compareResidualsForTransientProblems do nothing";
    }


    void loadSCMDB()
    {
        LOG( INFO ) << "CRBTrilinear::loadSCMDB do nothing";
    }

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


    //only to compile
    wn_type wndu() const
    {
        return M_model->rBFunctionSpace()->dualRB();
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
    element_type expansion( vectorN_type const& u , int const N, wn_type const & WN ) const;


    /**
     * export basis functions to visualize it
     * \param wn : tuple composed of a vector of wn_type and a vector of string (used to name basis)
     */
    void exportBasisFunctions( const export_vector_wn_type& wn )const ;

    boost::tuple<double,double, solutions_tuple, matrix_info_tuple, double, double, upper_bounds_tuple > run( parameter_type const& mu, double eps = 1e-6, int N = -1, bool print_rb_matrix=false );

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

    //! set boolean indicates if we are in offline_step or not
    void setOfflineStep( bool b )
    {
        M_offline_step = b;
    }

    wn_type wn() const
    {
        return M_model->rBFunctionSpace()->primalRB();
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
    crb_elements_db_type M_elements_database;

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

    bool M_offline_step;

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
    //M_model->initModel();
    LOG( INFO )<< " -- model init done in " << ti.elapsed() << "s";

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

        LOG( INFO )<<"[CRB offline] M_error_type = "<<M_error_type<<std::endl;

        LOG(INFO) << " -- sampling init done in " << ti.elapsed() << "s";
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

    boost::tie( boost::tuples::ignore, Aqm, Fqm ) = M_model->computeAffineDecomposition();

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
            M_WNmu->clear();
            std::vector< parameter_type > V;
            parameter_type __mu;
            __mu = M_Dmu->element();
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
            M_WNmu->setElements( V );
            M_iter_max = M_WNmu->size();
            M_WNmu->writeOnFile(file_name);
        }
        else
        {
            M_WNmu->clear();
            int sampling_size = M_WNmu->readFromFile(file_name);
            M_iter_max = sampling_size;
        }
        mu = M_WNmu->at( M_N ); // first element
        //std::cout<<" [use_predefined_WNmu] mu = \n"<<mu<<std::endl;

        LOG( INFO )<<"[CRB::offline] read WNmu ( sampling size : "<<M_iter_max<<" )";

    }

    if( M_error_type == CRB_NO_RESIDUAL || use_predefined_WNmu )
    {
        //in this case it makes no sens to check the estimated error
        M_maxerror = 1e10;
    }


    LOG(INFO) << "[CRBTrilinear::offline] strategy "<< M_error_type <<"\n";

    while ( M_maxerror > M_tolerance && M_N < M_iter_max )
    {

        boost::timer timer, timer2;
        LOG(INFO) <<"========================================"<<"\n";
        if( proc_number == this->worldComm().masterRank() )
            std::cout<<"construction of "<<M_N<<"/"<<M_iter_max<<" basis "<<std::endl;
        LOG(INFO) << "N=" << M_N << "/"  << M_iter_max << "( nb proc : "<<worldComm().globalSize()<<")";

        // for a given parameter \p mu assemble the left and right hand side
        u->setName( ( boost::format( "fem-primal-N%1%-proc%2%" ) % (M_N)  % proc_number ).str() );

        mu.check();
        u->zero();

        timer2.restart();

        LOG(INFO) << "[CRB::offline] solving primal" << "\n";
        *u = M_model->solve( mu );

        //if( proc_number == this->worldComm().masterRank() ) std::cout << "  -- primal problem solved in " << timer2.elapsed() << "s\n";
        timer2.restart();


        if( ! use_predefined_WNmu )
            M_WNmu->push_back( mu, index );

        M_WNmu_complement = M_WNmu->complement();

        M_model->rBFunctionSpace()->addPrimalBasisElement( *u );
        //WARNING : the dual element is not the real dual solution !
        //no dual problem was solved
        M_model->rBFunctionSpace()->addDualBasisElement( *u );

	    int number_of_added_elements=1;
        M_N+=number_of_added_elements;

        if ( orthonormalize_primal )
        {
            orthonormalize( M_N, M_model->rBFunctionSpace()->primalRB() , number_of_added_elements );
            orthonormalize( M_N, M_model->rBFunctionSpace()->primalRB() , number_of_added_elements );
            orthonormalize( M_N, M_model->rBFunctionSpace()->primalRB() , number_of_added_elements );
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
                    M_Aqm_pr[q][0]( i, j ) = Aqm[q][0]->energy( M_model->rBFunctionSpace()->primalBasisElement(i) , M_model->rBFunctionSpace()->primalBasisElement(j) );
                }
            }

            for ( size_type j=M_N-number_of_added_elements; j < M_N; j++ )
            {
                for ( size_type i = 0; i < M_N; ++i )
                {
                    M_Aqm_pr[q][0]( i, j ) = Aqm[q][0]->energy( M_model->rBFunctionSpace()->primalBasisElement(i), M_model->rBFunctionSpace()->primalBasisElement(j) );
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
                M_Fqm_pr[q][0]( index ) = M_model->Fqm( 0, q, 0, M_model->rBFunctionSpace()->primalBasisElement(index) );
            }
        }//loop over q

        LOG(INFO) << "[CRB::offline] compute Lq_pr" << "\n";

        for ( size_type q = 0; q < M_model->Ql( M_output_index ); ++q )
        {
            M_Lqm_pr[q][0].conservativeResize( M_N );

            for ( size_type l = 1; l <= number_of_added_elements; ++l )
            {
                int index = M_N-l;
                M_Lqm_pr[q][0]( index ) = M_model->Fqm( M_output_index, q, 0, M_model->rBFunctionSpace()->primalBasisElement(index) );
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
                trilinear_form  = M_model->computeTrilinearForm( M_model->rBFunctionSpace()->primalBasisElement(k) );

                M_Aqm_tril_pr[q][k].conservativeResize( M_N, M_N );
                for ( int i = 0; i < M_N; ++i )
                {
                    for ( int j = 0; j < M_N; ++j )
                    {
                        M_Aqm_tril_pr[q][k]( i, j ) = trilinear_form->energy( M_model->rBFunctionSpace()->primalBasisElement(j), M_model->rBFunctionSpace()->primalBasisElement(i) );
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
        LOG(INFO) <<"========================================"<<"\n";

        //save DB after adding an element
        this->saveDB();
        M_elements_database.setWn( boost::make_tuple( M_model->rBFunctionSpace()->primalRB() , M_model->rBFunctionSpace()->dualRB() ) );
        M_elements_database.saveDB();
    }


    LOG( INFO )<<"number of elements in the reduced basis : "<<M_N<<" ( nb proc : "<<worldComm().globalSize()<<")";
    LOG(INFO) << " index choosen : ";
    BOOST_FOREACH( auto id, M_index )
    LOG(INFO)<<id<<" ";
    LOG(INFO)<<"\n";
    bool visualize_basis = this->vm()["crb.visualize-basis"].template as<bool>() ;

    if ( visualize_basis )
    {
        std::vector<wn_type> wn;
        std::vector<std::string> names;
        wn.push_back( M_model->rBFunctionSpace()->primalRB() );
        names.push_back( "primal" );
        exportBasisFunctions( boost::make_tuple( wn ,names ) );

        if ( orthonormalize_primal )
            LOG(INFO)<<"[CRB::offline] Basis functions have been exported but warning elements have been orthonormalized";
    }

    LOG( INFO ) << "Offline CRB is done";

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

    if ( N > M_N ) N = M_N;
    beta_vector_type betaAqm;
    std::vector<beta_vector_type> betaFqm, betaLqm;

    //-- end of initialization step

    //vector containing outputs from time=time_step until time=time_for_output
    double output;

    boost::tie( boost::tuples::ignore,  betaAqm, betaFqm  ) = M_model->computeBetaQm( mu );

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
            LOG( INFO ) << "[CRBTrilinear::lb] i/N = " << i+1 << "/" << Nmax ;
            LOG( INFO ) << "[CRBTrilinear::lb] intermediary Grashof = " << current_Grashofs;
            LOG( INFO ) << "[CRBTrilinear::lb] and Prandtl = " << current_Prandtl ;
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

    double condition_number = 0;

    LOG( INFO ) <<"[CRBTrilinear::lb] The condition number of jacobian done is not computed ";
    vectorN_type L ( ( int )N );
    L.setZero( N );

    for ( size_type q = 0; q < M_model->Ql( M_output_index ); ++q )
    {
        L += betaFqm[M_output_index][q][0]*M_Lqm_pr[q][0].head( N );
    }
    output = L.dot( uN );
    LOG(INFO) << "[CRBTrilinear::lb] computation of the output done";

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

    boost::tie( boost::tuples::ignore,  betaAqm, betaFqm ) = M_model->computeBetaQm( mu );

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
        auto expansionX = expansion( map_X , N , M_model->rBFunctionSpace()->primalRB() );
        auto J = M_model->jacobian( expansionX );
        matrixN_type model_reduced_jacobian( N , N );
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
                model_reduced_jacobian(i,j) = J->energy( M_model->rBFunctionSpace()->primalBasisElement(i), M_model->rBFunctionSpace()->primalBasisElement(j) );
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
        auto expansionX = expansion( map_X , N , M_model->rBFunctionSpace()->primalRB() );
        auto R = M_model->residual( expansionX );
        vectorN_type model_reduced_residual( N );
        element_ptrtype eltR( new element_type( M_model->functionSpace() ) );
        for(int i=0; i<eltR->localSize();i++)
            eltR->operator()(i)=R->operator()(i);
        for(int i=0; i<N; i++)
            model_reduced_residual(i) = inner_product( *eltR , M_model->rBFunctionSpace()->primalBasisElement(i) );
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
}

template<typename TruthModelType>
typename boost::tuple<double,double, typename CRBTrilinear<TruthModelType>::solutions_tuple, typename CRBTrilinear<TruthModelType>::matrix_info_tuple,
                      double, double, typename CRBTrilinear<TruthModelType>::upper_bounds_tuple >
CRBTrilinear<TruthModelType>::run( parameter_type const& mu, double eps , int N, bool print_rb_matrix )
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

    std::vector< std::vector<double> > coefficients;
    auto upper_bounds = boost::make_tuple(-1 , -1 , -1 , coefficients , coefficients );
    std::vector< vectorN_type > vector_uN( 1 );
    vector_uN[0] = uN ;
    auto solutions = boost::make_tuple( vector_uN, vector_uN, vector_uN, vector_uN );
    auto matrix_info = boost::make_tuple(-1,-1);
    return boost::make_tuple( output , Nwn , solutions, matrix_info , -1 , -1 , upper_bounds );

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
    return Feel::expansion( M_model->rBFunctionSpace()->primalRB(), uN , Nwn);
}


template<typename TruthModelType>
typename CRBTrilinear<TruthModelType>::element_type
CRBTrilinear<TruthModelType>::expansion( vectorN_type const& u , int const N,  wn_type const & WN) const
{
    int Nwn;

    if( N > 0 )
        Nwn = N;
    else
        Nwn = M_N;

    //FEELPP_ASSERT( N == u.size() )( N )( u.size() ).error( "invalid expansion size");
    //FEELPP_ASSERT( Nwn == u.size() )( Nwn )( u.size() ).error( "invalid expansion size");
    //int size = uN.size();
    FEELPP_ASSERT( Nwn <= WN.size() )( Nwn )( WN.size() ).error( "invalid expansion size ( N and M_WN ) ");
    return Feel::expansion( WN , u, Nwn );
}


template<typename TruthModelType>
void
CRBTrilinear<TruthModelType>::computationalTimeStatistics(std::string appname)
{

    double min=0,max=0,mean=0,standard_deviation=0;

    int n_eval = option(_name="crb.computational-time-neval").template as<int>();

    Eigen::Matrix<double, Eigen::Dynamic, 1> time_crb;
    time_crb.resize( n_eval );

    sampling_ptrtype Sampling( new sampling_type( M_Dmu ) );
    Sampling->logEquidistribute( n_eval  );

    bool cvg = option(_name="crb.cvg-study").template as<bool>();
    int dimension = this->dimension();
    double tol = option(_name="crb.online-tolerance").template as<double>();

    int N=dimension;//by default we perform only one time statistics

    if( cvg ) //if we want to compute time statistics for every crb basis then we start at 1
        N=1;

    int proc_number =  Environment::worldComm().globalRank();
    int master =  Environment::worldComm().masterRank();

    //write on a file
    std::string file_name = "cvg-timing-crb.dat";

    std::ofstream conv;
    if( proc_number == master )
    {
        conv.open(file_name, std::ios::app);
        conv << "NbBasis" << "\t" << "min" <<"\t"<< "max" <<"\t"<< "mean"<<"\t"<<"standard_deviation" << "\n";
    }
    bool print_rb_matrix=false;
    //loop over basis functions (if cvg option)
    for(; N<=dimension; N++)
    {

        int mu_number = 0;
        BOOST_FOREACH( auto mu, *Sampling )
        {
            boost::mpi::timer tcrb;
            auto o = this->run( mu, tol , N, print_rb_matrix);
            time_crb( mu_number ) = tcrb.elapsed() ;
            mu_number++;
        }

        auto stat = M_model->computeStatistics( time_crb , appname );

        min=stat(0);
        max=stat(1);
        mean=stat(2);
        standard_deviation=stat(3);

        if( proc_number == master )
            conv << N << "\t" << min << "\t" << max<< "\t"<< mean<< "\t"<< standard_deviation<<"\n";
    }//loop over basis functions
    conv.close();
}

template<typename TruthModelType>
template<class Archive>
void
CRBTrilinear<TruthModelType>::save( Archive & ar, const unsigned int version ) const
{
    int proc_number = this->worldComm().globalRank();

    LOG(INFO) <<"[CRBTrilinear::save] version : "<<version<<std::endl;

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

    ar & BOOST_SERIALIZATION_NVP( M_maxerror );

}

template<typename TruthModelType>
template<class Archive>
void
CRBTrilinear<TruthModelType>::load( Archive & ar, const unsigned int version )
{

    int proc_number = this->worldComm().globalRank();

    LOG(INFO) <<"[CRBTrilinear::load] version"<< version <<std::endl;

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
