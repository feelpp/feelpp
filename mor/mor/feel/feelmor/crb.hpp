/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-11-24

  Copyright (C) 2009-2012 Universite Joseph Fourier (Grenoble I)

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
   \file crb.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Alexandre Ancel <alexandre.ancel@cemosis.fr>
   \author Cecile Daversin <daversin@math.unistra.fr>
   \author Stephane Veys
   \date 2009-11-24
 */
#ifndef __CRB_H
#define __CRB_H 1


#if 0
#include <boost/multi_array.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/bimap.hpp>
#include <boost/bimap/support/lambda.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/algorithm/string.hpp>

#include <fstream>


#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_member.hpp>
#endif

#include <vector>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>

#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/solvereigen.hpp>
//#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>
//#include <feel/feelcore/parameter.hpp>
#include <feel/feelcore/serialization.hpp>

#include <feel/feeldiscr/functionspace.hpp>
//#include <feel/feeldiscr/subelements.hpp>
#include <feel/feeldiscr/expansion.hpp>
#include <feel/feelts/bdf.hpp>

#include <feel/feelmor/crbenums.hpp>
#include <feel/feelmor/crbdata.hpp>
#include <feel/feelmor/options.hpp>
#include <feel/feelmor/parameterspace.hpp>
#include <feel/feelmor/crbdb.hpp>
#include <feel/feelmor/crbscm.hpp>
#include <feel/feelmor/crbelementsdb.hpp>
#include <feel/feelmor/pod.hpp>

//#include <feel/feelfilters/exporter.hpp>

#include <feel/feelcore/pslogger.hpp>

#if defined(FEELPP_HAS_HARTS)

#include "hartsconfig.h"
#include "HARTS.h"
#if defined(HARTS_HAS_OPENCL)

#define __CL_ENABLE_EXCEPTIONS

#ifdef __APPLE__
// cl.hpp is not included on OS X, have to rely on a custom cl.hpp file
// provided by viennacl
//#include <OpenCL/cl.hpp>
#include "cl.hpp"
#else
#include <CL/cl.hpp>
#endif

#include <feel/feelmor/crbclcontext.hpp>
//#include "feel/feelmor/crb.cl.hpp"

// declare that we want to use a custom context
#define VIENNACL_WITH_OPENCL

// ViennaCL includes
//
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/coordinate_matrix.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/ilu.hpp"
#include "viennacl/linalg/jacobi_precond.hpp"
#include "viennacl/linalg/cg.hpp"
#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/linalg/gmres.hpp"
#include "viennacl/io/matrix_market.hpp"
#include "viennacl/ocl/backend.hpp"

#endif
#endif


namespace Feel
{

template <typename SpaceType, typename ParameterSpaceType>
class CRBBase
{
public :
    typedef SpaceType space_type;
    typedef ParameterSpaceType parameterspace_type;

    typedef typename space_type::element_type element_type;
    typedef typename parameterspace_type::element_type parameter_type;

    typedef typename space_type::value_type value_type;
    typedef Backend<value_type> backend_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    using vectorN_type = Feel::vectorN_type;
    typedef boost::tuple< double,double > matrix_info_tuple;


    virtual boost::tuple<std::vector<double>,matrix_info_tuple> lb( size_type N,
                                                                    parameter_type const& mu,
                                                                    std::vector< vectorN_type >& uN,
                                                                    std::vector< vectorN_type >& uNdu ,
                                                                    std::vector<vectorN_type> & uNold,
                                                                    std::vector<vectorN_type> & uNduold,
                                                                    bool print_rb_matrix=false, int K=0,
                                                                    bool computeOutput = true ) const =0;


    virtual double computeRieszResidualNorm( parameter_type const& mu ) const=0;
    virtual double computeRieszResidualNorm( parameter_type const& mu, std::vector<vectorN_type> const& uN ) const=0;
    virtual element_type expansion( vectorN_type const& u , int N = -1, bool dual=false ) const=0;
    virtual int dimension() const=0;
    virtual std::pair<int,double> online_iterations()=0;
    virtual int WNmuSize() const=0;
    virtual element_type solveFemModelUsingAffineDecomposition( parameter_type const& mu ) const=0;

};

/**
 * \brief Certifed Reduced Basis class
 *
 * Implements the certified reduced basis method
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename TruthModelType>
class CRB : public CRBDB,
            public CRBBase<typename TruthModelType::functionspace_type, typename TruthModelType::parameterspace_type>
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
    typedef std::shared_ptr<truth_model_type> truth_model_ptrtype;

    typedef double value_type;
    typedef boost::tuple<double,double> bounds_type;

    typedef typename TruthModelType::parameterspace_type parameterspace_type;
    typedef typename TruthModelType::parameterspace_ptrtype parameterspace_ptrtype;
    typedef typename parameterspace_type::element_type parameter_type;
    typedef typename parameterspace_type::element_ptrtype parameter_ptrtype;
    typedef typename parameterspace_type::sampling_type sampling_type;
    typedef typename parameterspace_type::sampling_ptrtype sampling_ptrtype;

    //typedef boost::tuple<double, parameter_type, size_type, double, double> relative_error_type;
    typedef boost::tuple<double, parameter_type, double, double> relative_error_type;
    typedef relative_error_type max_error_type;

    typedef boost::tuple<std::vector<double>, std::vector< std::vector<double> > , std::vector< std::vector<double> >, double, double > error_estimation_type;
    typedef boost::tuple<double, std::vector<double> > residual_error_type;

    typedef boost::bimap< int, boost::tuple<double,double,double> > convergence_type;

    typedef typename convergence_type::value_type convergence;

    typedef CRB self_type;
    using self_ptrtype = std::shared_ptr<self_type>;

    //! scm
    typedef CRBSCM<truth_model_type> scm_type;
    typedef std::shared_ptr<scm_type> scm_ptrtype;

    //! elements database
    typedef CRBElementsDB<truth_model_type> crb_elements_db_type;
    typedef std::shared_ptr<crb_elements_db_type> crb_elements_db_ptrtype;

    //! POD
    typedef POD<truth_model_type> pod_type;
    typedef std::shared_ptr<pod_type> pod_ptrtype;
    typedef typename pod_type::mode_set_type mode_set_type;

    //! function space type
    typedef typename model_type::functionspace_type functionspace_type;
    typedef typename model_type::functionspace_ptrtype functionspace_ptrtype;


    //! element of the functionspace type
    typedef typename model_type::element_type element_type;
    typedef typename model_type::element_ptrtype element_ptrtype;

    typedef typename model_type::backend_type backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;
    typedef typename model_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename model_type::vector_ptrtype vector_ptrtype;
    typedef typename model_type::beta_vector_type beta_vector_type;


    typedef Eigen::VectorXd y_type;
    typedef std::vector<y_type> y_set_type;
    typedef std::vector<boost::tuple<double,double> > y_bounds_type;

    typedef typename model_type::rbfunctionspace_type::rb_basis_type wn_type;
    typedef boost::tuple< std::vector<wn_type> , std::vector<std::string> > export_vector_wn_type;

    typedef std::vector<double> vector_double_type;
    typedef std::shared_ptr<vector_double_type> vector_double_ptrtype;

    using vectorN_type = Feel::vectorN_type;
    using matrixN_type = Feel::matrixN_type;

    typedef Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > map_dense_matrix_type;
    typedef Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > map_dense_vector_type;

    typedef boost::tuple< std::vector<vectorN_type> , std::vector<vectorN_type> , std::vector<vectorN_type>, std::vector<vectorN_type> > solutions_tuple;
    typedef boost::tuple< std::vector<double>,double,double , std::vector< std::vector< double > > , std::vector< std::vector< double > > > upper_bounds_tuple;
    typedef boost::tuple< double,double > matrix_info_tuple; //conditioning, determinant

    typedef boost::multi_array<value_type, 2> array_2_type;
    typedef boost::multi_array<vectorN_type, 2> array_3_type;
    typedef boost::multi_array<matrixN_type, 2> array_4_type;

    //! mesh type
    typedef typename model_type::mesh_type mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    //! space type
    typedef typename model_type::space_type space_type;
    typedef std::shared_ptr<space_type> space_ptrtype;

    //! time discretization
    typedef Bdf<space_type>  bdf_type;
    typedef std::shared_ptr<bdf_type> bdf_ptrtype;

    // ! export
    typedef Exporter<mesh_type> export_type;
    typedef std::shared_ptr<export_type> export_ptrtype;

    typedef Preconditioner<double> preconditioner_type;
    typedef std::shared_ptr<preconditioner_type> preconditioner_ptrtype;

    //here a fusion vector containing sequence 0 ... nb_spaces
    //useful to acces to a component of a composite space in ComputeIntegrals
    static const int nb_spaces = functionspace_type::nSpaces;
    typedef typename mpl::if_< boost::is_same< mpl::int_<nb_spaces> , mpl::int_<2> > , fusion::vector< mpl::int_<0>, mpl::int_<1> >  ,
                       typename mpl::if_ < boost::is_same< mpl::int_<nb_spaces> , mpl::int_<3> > , fusion::vector < mpl::int_<0> , mpl::int_<1> , mpl::int_<2> > ,
                                  typename mpl::if_< boost::is_same< mpl::int_<nb_spaces> , mpl::int_<4> >, fusion::vector< mpl::int_<0>, mpl::int_<1>, mpl::int_<2>, mpl::int_<3> >,
                                                     fusion::vector< mpl::int_<0>, mpl::int_<1>, mpl::int_<2>, mpl::int_<3>, mpl::int_<4> >
                                                     >::type >::type >::type index_vector_type;

#if defined(FEELPP_HAS_HARTS) && defined(HARTS_HAS_OPENCL)
    mutable crbCLContext clContext_;
#endif

    //@}

    /** @name Constructors, destructor
     */
    //@{

    static self_ptrtype New( crb::stage stage = crb::stage::online )
        {
            return New( "noname", std::make_shared<truth_model_type>(stage), stage );
        }

    static self_ptrtype New( std::string const& name, crb::stage stage = crb::stage::online )
        {
            return New( name, std::make_shared<truth_model_type>(stage), stage );
        }

    static self_ptrtype New( std::string const& name,
                             truth_model_ptrtype const & model,
                             crb::stage stage = crb::stage::online,
                             std::string const& prefixExt = "" )
        {
            auto crb = std::make_shared<self_type>( name, model, stage, prefixExt );
            if ( stage == crb::stage::offline )
                crb->init();
            return crb;
        }

    //! constructor from command line options
    CRB( std::string const& name,
         truth_model_ptrtype const & model,
         crb::stage stage = crb::stage::online,
         std::string const& prefixExt = "" )
        :
        super( name, prefixvm(prefixExt, "crb"), model->worldCommPtr()),
        M_prefix(model->prefix()),
        M_elements_database( name,
                             prefixvm(prefixExt,"elements"),
                             this->worldCommPtr() ),
        M_nlsolver( SolverNonLinear<double>::build( "petsc", "", Environment::worldCommSeqPtr() ) ),// this->worldComm() ) ),
        M_model( model ),
        M_output_index( ioption(_prefix=M_prefix,_name="crb.output-index") ),
        M_tolerance( doption(_prefix=M_prefix,_name="crb.error-max") ),
        M_iter_max( ioption(_prefix=M_prefix,_name="crb.dimension-max") ),
        M_factor( ioption(_prefix=M_prefix,_name="crb.factor") ),
        M_error_type( CRBErrorType( ioption(_prefix=M_prefix,_name="crb.error-type" ) ) ),
        // M_Dmu( new parameterspace_type ),
        // M_Xi( new sampling_type( M_Dmu ) ),
        // M_WNmu( new sampling_type( M_Dmu, 0, M_Xi ) ),
        // M_WNmu_complement(),
        // M_primal_apee_mu( new sampling_type( M_Dmu, 0, M_Xi ) ),
        // M_dual_apee_mu( new sampling_type( M_Dmu, 0, M_Xi ) ),
        M_scmA( new scm_type( name, prefixvm(prefixExt,"scma"), false /*not scm for mass mastrix*/, this->worldCommPtr() )  ),
        M_scmM( new scm_type( name, prefixvm(prefixExt,"scmm"), true /*scm for mass matrix*/, this->worldCommPtr() ) ),
        M_N( 0 ),
        M_solve_dual_problem( boption(_prefix=M_prefix,_name="crb.solve-dual-problem") ),
        M_orthonormalize_primal( boption(_prefix=M_prefix,_name="crb.orthonormalize-primal") ),
        M_orthonormalize_dual( M_solve_dual_problem && boption(_prefix=M_prefix,_name="crb.orthonormalize-dual") ),
        M_compute_variance( boption(_prefix=M_prefix,_name="crb.compute-variance") ),
        M_computeMatrixInfo( boption(_prefix=M_prefix,_name="crb.compute-matrix-information") ),
        M_database_contains_variance_info( boption(_prefix=M_prefix,_name="crb.save-information-for-variance") ),
        M_use_newton( boption(_prefix=M_prefix,_name="crb.use-newton") ),
        M_offline_step( false ),
        M_rebuild( boption(_prefix=M_prefix,_name="crb.rebuild-database") || ioption(_prefix=M_prefix,_name="crb.restart-from-N") == 0 ),
        M_SER_adapt( false ),
        M_SER_maxerr( 0 ),
        M_SER_errorEstimation( boption(_prefix=M_prefix,_name="ser.error-estimation") ),
        M_SER_useGreedyInRb( boption(_prefix=M_prefix,_name="ser.use-greedy-in-rb") ),
        M_model_executed_in_steady_mode( boption(_prefix=M_prefix,_name="crb.is-model-executed-in-steady-mode") ),
        M_save_output_behavior( boption(_prefix=M_prefix,_name="crb.save-output-behavior") ),
        M_fixedpointMaxIterations( ioption(_prefix=M_prefix,_name="crb.fixedpoint.maxit") ),
        M_fixedpointIncrementTol( doption(_prefix=M_prefix,_name="crb.fixedpoint.increment-tol") ),
        M_fixedpointVerbose( boption(_prefix=M_prefix,_name="crb.fixedpoint.verbose") ),
        M_fixedpointCriticalValue( doption(_prefix=M_prefix,_name="crb.fixedpoint.critical-value") ),
        M_fixedpointUseAitken( boption(_prefix=M_prefix,_name="crb.fixedpoint.aitken") ),
        M_loadElementsDb( boption(_prefix=M_prefix,_name="crb.load-elements-database") ),
        M_useRbSpaceContextEim( boption(_prefix=M_prefix,_name="crb.use-fast-eim") ),
        M_hasRbSpaceContextEim( false ),
        M_useAccurateApee( boption(_prefix=M_prefix,_name="crb.use-accurate-apee") ),
        M_computeApeeForEachTimeStep( boption(_prefix=M_prefix,_name="crb.compute-apee-for-each-time-step") ),
        M_seekMuInComplement( boption(_prefix=M_prefix,_name="crb.seek-mu-in-complement") ),
        M_showResidual( boption(_prefix=M_prefix,_name="crb.show-residual") ),
        M_check_cvg( false )
        {
            this->setTruthModel( model );
            if( stage == crb::stage::offline )
            {
                // define offline backend and preconditioner
                M_backend =  backend();
                M_backend_primal = backend(_name="backend-primal");
                M_backend_dual = backend(_name="backend-dual");

                if( boption(_prefix=M_prefix,_name="crb.use-primal-pc") )
                {
                    M_preconditioner_primal = preconditioner(_pc=(PreconditionerType) M_backend_primal->pcEnumType(), // by default : lu in seq or wirh mumps, else gasm in parallel
                                                             _backend= M_backend_primal,
                                                             _pcfactormatsolverpackage=(MatSolverPackageType) M_backend_primal->matSolverPackageEnumType(),// mumps if is installed ( by defaut )
                                                             _worldcomm=M_backend_primal->comm(),
                                                             _prefix=M_backend_primal->prefix() ,
                                                             _rebuild=M_model->useSER());
                }

                M_preconditioner_dual = preconditioner(_pc=(PreconditionerType) M_backend_dual->pcEnumType(), // by default : lu in seq or wirh mumps, else gasm in parallel
                                                       _backend= M_backend_dual,
                                                       _pcfactormatsolverpackage=(MatSolverPackageType) M_backend_dual->matSolverPackageEnumType(),// mumps if is installed ( by defaut )
                                                       _worldcomm=M_backend_dual->comm(),
                                                       _prefix=M_backend_dual->prefix() ,
                                                       _rebuild=M_model->useSER());
            }
        }



    //! Copy constructor
    CRB( CRB const & o )
        :
        super( o ),
        M_elements_database( o.M_elements_database ),
        M_nlsolver( o.M_nlsolver ),
        M_output_index( o.M_output_index ),
        M_tolerance( o.M_tolerance ),
        M_iter_max( o.M_iter_max ),
        M_factor( o.M_factor ),
        M_error_type( o.M_error_type ),
        M_maxerror( o.M_maxerror ),
        M_Dmu( o.M_Dmu ),
        M_Xi( o.M_Xi ),
        M_WNmu( o.M_WNmu ),
        M_WNmu_complement( o.M_WNmu_complement ),
        M_primal_apee_mu( o.M_primal_apee_mu ),
        M_dual_apee_mu( o.M_dual_apee_mu ),
        M_C0_pr( o.M_C0_pr ),
        M_C0_du( o.M_C0_du ),
        M_Lambda_pr( o.M_Lambda_pr ),
        M_Lambda_du( o.M_Lambda_du ),
        M_Gamma_pr( o.M_Gamma_pr ),
        M_Gamma_du( o.M_Gamma_du ),
        M_C0_pr_eim( o.M_C0_pr ),
        M_C0_du_eim( o.M_C0_du ),
        M_Lambda_pr_eim( o.M_Lambda_pr ),
        M_Lambda_du_eim( o.M_Lambda_du ),
        M_Gamma_pr_eim( o.M_Gamma_pr ),
        M_Gamma_du_eim( o.M_Gamma_du ),
        M_Cmf_pr( o.M_Cmf_pr ),
        M_Cmf_du( o.M_Cmf_du ),
        M_Cmf_du_ini( o.M_Cmf_du_ini ),
        M_Cma_pr( o.M_Cma_pr ),
        M_Cma_du( o.M_Cma_du ),
        M_Cmm_pr( o.M_Cmm_pr ),
        M_Cmm_du( o.M_Cmm_du ),
        M_Cmf_pr_eim( o.M_Cmf_pr ),
        M_Cmf_du_eim( o.M_Cmf_du ),
        M_Cmf_du_ini_eim( o.M_Cmf_du_ini ),
        M_Cma_pr_eim( o.M_Cma_pr ),
        M_Cma_du_eim( o.M_Cma_du ),
        M_Cmm_pr_eim( o.M_Cmm_pr ),
        M_Cmm_du_eim( o.M_Cmm_du ),
        M_coeff_pr_ini_online( o.M_coeff_pr_ini_online ),
        M_coeff_du_ini_online( o.M_coeff_du_ini_online ),
        M_rebuild( o.M_rebuild ),
        M_SER_adapt( o.M_SER_adapt ),
        M_SER_maxerr( o.M_SER_maxerr )
        {}

public:
    //! destructor
    ~CRB() override
        {}

    void init()
        {
            if ( !M_rebuild && !fs::exists(this->dbLocalPath()/fs::path(this->jsonFilename()) ) )
            {
                M_rebuild = true;
                this->worldComm().barrier();
                LOG(INFO) << fmt::format("CRB::init() : rebuild database");
            }
            if ( !M_rebuild )
            {
                LOG(INFO) << fmt::format("CRB::init() : load database from {}", (this->dbLocalPath()/fs::path(this->jsonFilename())).string() );
                this->setupOfflineFromDB();
            }
            else
            {
                M_scmM->setId( this->id() );
                M_scmA->setId( this->id() );
                M_elements_database.setId( this->id() );
                LOG(INFO) << fmt::format("CRB::init() : build database, initialized SCM_M and SCM_A");
            }
            if( this->worldComm().isMasterRank() )
                std::cout << "CRB::init() : Use DB id " << this->id() << std::endl;
            LOG(INFO) << fmt::format("CRB::init() : Use DB id {}", this->id() );

            if ( M_N == 0 )
            {
                if( this->worldComm().isMasterRank() )
                    std::cout << "CRB::init() : Databases does not exist or incomplete -> Start from the begining\n";
                LOG( INFO ) <<"CRB::init() : Databases does not exist or incomplete -> Start from the begining";
            }

            // fe vector is requiert in online : must not be TODO
            if ( M_use_newton && M_loadElementsDb && M_Rqm.empty() )
                boost::tie( boost::tuples::ignore, boost::tuples::ignore/*M_Jqm*/, M_Rqm ) = M_model->computeAffineDecomposition();


        }

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    CRB& operator=( CRB const & o )
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

    truth_model_ptrtype & model() { return M_model; }
    truth_model_ptrtype const& model() const { return M_model; }

    std::string const& prefix() const { return M_prefix; }

    //! return factor
    int factor() const
        {
            return factor;
        }
    //! \return max iterations
    int maxIter() const
        {
            return M_iter_max;
        }

    //! \return the parameter space
    parameterspace_ptrtype const& Dmu() const
        {
            return M_Dmu;
        }

    //! \return the output index
    int outputIndex() const
        {
            return M_output_index;
        }

    //! \return the output name
    std::string const& outputName() const
    {
        return M_output_name;
    }

    //! \return the dimension of the reduced basis space
    int dimension() const override
        {
            return M_N;
        }

    //! \return the current size of WNmu:
    //! may differ from the dimension of XN for multiphysic problems
    int WNmuSize() const override
        {
            return M_N;
        }

    //! \return the train sampling used to generate the reduced basis space
    sampling_ptrtype const& trainSampling() const
        {
            return M_Xi;
        }

    //! \return the error type
    CRBErrorType errorType() const
        {
            return M_error_type;
        }

    //! \return true if solveDualProblem
    bool solveDualProblem() const
        {
            return M_solve_dual_problem;
        }

    //! \return the scm object
    scm_ptrtype const& scm() const
        {
            return M_scmA;
        }

    void loadSCMDB();

    //@}

    /** @name  Mutators
     */
    //@{

    //! set the output index
    void setOutputIndex( uint16_type oindex )
        {
            int M_prev_o = M_output_index;
            M_output_index = oindex;

            if ( ( size_type )M_output_index >= M_model->Nl()  )
                M_output_index = M_prev_o;

            //std::cout << " -- crb set output index to " << M_output_index << " (max output = " << M_model->Nl() << ")\n";

            if ( M_output_index != M_prev_o )
                this->loadDB();

            //std::cout << "Database " << this->lookForDB() << " available and loaded\n";

        }
    /**
     * @brief Set the Output Name 
     * 
     * @param oname name of the outout
     */
    void setOutputName( std::string const& oname )
    {
        M_output_name = oname;
    }
    //! set the crb error type
    void setCRBErrorType( CRBErrorType error )
        {
            M_error_type = error;
        }

    //! set offline tolerance
    void setTolerance( double tolerance )
        {
            M_tolerance = tolerance;
        }

    //! set solve dual problem
    void setSolveDualProblem( bool b )
        {
            M_solve_dual_problem = b;
        }

    //! set the truth offline model
    void setTruthModel( truth_model_ptrtype const& model )
        {
            if ( !model )
                return;
            M_model = model;
            this->setDBDirectory( M_model->uuid() );
            M_Dmu = M_model->parameterSpace();
            M_Xi = std::make_shared<sampling_type>( M_Dmu );
            M_WNmu = std::make_shared<sampling_type>( M_Dmu, 0, M_Xi );
            //M_WNmu_complement(),
            M_primal_apee_mu = std::make_shared<sampling_type>( M_Dmu, 0, M_Xi );
            M_dual_apee_mu = std::make_shared<sampling_type>( M_Dmu, 0, M_Xi );

            M_elements_database.setModel( model );
#if 0

            if ( ! loadDB() )
                M_WNmu = sampling_ptrtype( new sampling_type( M_Dmu ) );
            else
            {
                LOG(INFO) << "Database " << this->lookForDB() << " available and loaded\n";
            }
#endif
            if ( M_scmA )
                M_scmA->setTruthModel( M_model );
            if ( M_scmM )
                M_scmM->setTruthModel( M_model );
        }

    //!
    //! return the SER Level, default is 0
    //!
    int level() const { return M_model->level(); }

    //! set max iteration number
    void setMaxIter( int K )
        {
            M_iter_max = K;
        }

    //! set factor
    void setFactor( int Factor )
        {
            M_factor = Factor;
        }

    //!
    //! set load finite element basis functions from DB  to \p r
    //!
    void setLoadBasisFromDB( bool r )
        {
            M_loadElementsDb = r;
        }

    //!
    //! @return boolean to load element basis functions from DB
    //!
    bool loadBasisFromDB() const
    {
        return M_loadElementsDb;
    }

    //! set boolean indicates if we are in offline_step or not
    void setOfflineStep( bool b )
        {
            M_offline_step = b;
        }
    //! get  boolean indicates if we are in offline_step or not
    bool offlineStep() const
        {
            return M_offline_step;
        }

    virtual std::vector<bool> hasZeroMean()
        { return {false}; }

    void setCheckCvg( bool check )
        { M_check_cvg = check; }

    bool onlineHasConverged()
        {return M_last_online_converged; }

    virtual std::map<std::string,double> timerMap() const
        {
            std::map<std::string,double> m;
            return m;
        }

    struct ComputePhi
    {

        ComputePhi( wn_type & v)
            :
            M_curs( 0 ),
            M_vect( v )
            {}

        template< typename T >
        void
        operator()( const T& t ) const
            {
                mesh_ptrtype mesh = t.functionSpace()->mesh();
                double surface = integrate( _range=elements(mesh), _expr=vf::cst(1.) ).evaluate()(0,0);
                double mean = integrate( _range=elements(mesh), _expr=vf::idv( t ) ).evaluate()(0,0);
                mean /= surface;
                auto element_mean = vf::project(t.functionSpace(), elements(mesh), vf::cst(mean) );
                auto element_t = vf::project(t.functionSpace(), elements(mesh), vf::idv(t) );
                int first_dof = t.start();
                int nb_dof = t.functionSpace()->nLocalDof();
                for(int dof=0 ; dof<nb_dof; dof++)
                    M_vect[M_curs][first_dof + dof] = element_t[dof] - element_mean[dof];

                ++M_curs;
            }

    private :
        mutable int M_curs;
        wn_type & M_vect;
    };


    struct ComputeIntegrals
    {

        ComputeIntegrals( element_type const& composite_e1 ,
                          element_type const& composite_e2,
                          std::vector<double> & vect )
            :
            M_composite_e1 ( composite_e1 ),
            M_composite_e2 ( composite_e2 ),
            M_vect( vect )
            {}

        template< typename T >
        void
        operator()( const T& t ) const
            {
                auto e1 = M_composite_e1.template element< T::value >();
                auto e2 = M_composite_e2.template element< T::value >();
                mesh_ptrtype mesh = e1.functionSpace()->mesh();
                double integral = integrate( _range=elements(mesh) , _expr=trans( vf::idv( e1 ) ) * vf::idv( e2 ) ).evaluate()(0,0);

                M_vect.push_back( integral );
            }

    private :
        std::vector<double> & M_vect;
        element_type const& M_composite_e1;
        element_type const& M_composite_e2;
    };

    struct ComputeIntegralsSquare
    {

        ComputeIntegralsSquare( element_type const& composite_e1 ,
                                element_type const& composite_e2,
                                vectorN_type & error )
            :
            M_composite_e1 ( composite_e1 ),
            M_composite_e2 ( composite_e2 ),
            M_error( error )
            {}

        template< typename T >
        void
        operator()( const T& t ) const
            {
                int i = T::value;
                if( i == 0 )
                    M_error.resize( 1 );
                else
                    M_error.conservativeResize( i+1 );

                auto e1 = M_composite_e1.template element< T::value >();
                auto e2 = M_composite_e2.template element< T::value >();

                auto Xh = M_composite_e1.functionSpace();
                mesh_ptrtype mesh = Xh->mesh();
                double integral = integrate( _range=elements(mesh) ,
                                             _expr=trans( vf::idv( e1 ) - vf::idv( e2 ) )
                                             *     ( vf::idv( e1 ) - vf::idv( e2 ) )
                                             ).evaluate()(0,0);
                M_error(i) = math::sqrt( integral ) ;
            }

    private:
        vectorN_type & M_error;
        element_type const& M_composite_e1;
        element_type const& M_composite_e2;
    };

    //@}

    /** @name  Methods
     */
    //@{
    /**
     * orthonormalize the basis
     * return the norm of the matrix A(i,j)=M_model->scalarProduct( WN[j], WN[i] ), should be 0
     */
    double orthonormalize( size_type N, wn_type& wn, int Nm = 1 );

    void checkResidual( parameter_type const& mu, std::vector< std::vector<double> > const& primal_residual_coeffs,
                        std::vector< std::vector<double> > const& dual_residual_coeffs , element_type & u, element_type & udu ) const;

    virtual void testResidual() {}

    void compareResidualsForTransientProblems(int N, parameter_type const& mu, std::vector<element_type> const & Un,
                                              std::vector<element_type> const & Unold, std::vector<element_type> const& Undu,
                                              std::vector<element_type> const & Unduold,
                                              std::vector< std::vector<double> > const& primal_residual_coeffs,
                                              std::vector< std::vector<double> > const& dual_residual_coeffs  ) const ;


    void buildFunctionFromRbCoefficients(int N, std::vector< vectorN_type > const & RBcoeff, wn_type const & WN, std::vector<element_ptrtype> & FEMsolutions );

    /*
     * check orthonormality
     * return the norm of the matrix A(i,j)=M_model->scalarProduct( WN[j], WN[i] ), should be 0
     */
    double checkOrthonormality( int N, const wn_type& wn ) const;

    /**
     * check the reduced basis space invariant properties
     * \param N dimension of \f$W_N\f$
     */
    void check( size_type N )  const;


    /**
     * export basis functions to visualize it
     * \param wn : tuple composed of a vector of wn_type and a vector of string (used to name basis)
     */
    void exportBasisFunctions( const export_vector_wn_type& wn )const ;
    virtual void exportBasisFunctions()
    {
        wn_type WN = this->wn();
        std::vector<wn_type> WN_vec = std::vector<wn_type>();
        WN_vec.push_back(WN);
        std::vector<std::string> name_vec = std::vector<std::string>(1, "primal");
        export_vector_wn_type exportWn = boost::make_tuple(WN_vec, name_vec);
        this->exportBasisFunctions(exportWn);
    }

    /**
     * Returns the lower bound of the output
     *
     * \param mu \f$ \mu\f$ the parameter at which to evaluate the output
     * \param N the size of the reduced basis space to use
     * \param uN primal solution
     * \param uNdu dual solution
     * \param K : index of time ( time = K*dt) at which we want to evaluate the output
     * \param computeOutput : if true compute quantity of interest (output)
     * Note : K as a default value for non time-dependent problems
     *
     *\return compute online the lower bound
     *\and also condition number of matrix A
     */

    //    boost::tuple<double,double> lb( size_type N, parameter_type const& mu, std::vector< vectorN_type >& uN, std::vector< vectorN_type >& uNdu , std::vector<vectorN_type> & uNold=std::vector<vectorN_type>(), std::vector<vectorN_type> & uNduold=std::vector<vectorN_type>(), int K=0) const;
    boost::tuple<std::vector<double>,matrix_info_tuple> lb( size_type N, parameter_type const& mu, std::vector< vectorN_type >& uN, std::vector< vectorN_type >& uNdu ,
                                                                    std::vector<vectorN_type> & uNold, std::vector<vectorN_type> & uNduold, bool print_rb_matrix=false, int K=0,
                                                                    bool computeOutput = true ) const override;


    /*
     * update the jacobian ( offline step )
     * \param X : solution
     * \param J : jacobian
     * \param mu : current parameter
     */
    void offlineUpdateJacobian( const vector_ptrtype & X, sparse_matrix_ptrtype & J , const parameter_type & mu) ;

    /*
     * update the jacobian (online step )
     * \param map_X : solution
     * \param map_J : jacobian
     * \param mu : current parameter
     * \param N : dimension of the reduced basis
     */
    virtual void updateJacobian( const map_dense_vector_type& map_X, map_dense_matrix_type& map_J , const parameter_type& mu, int N) const ;

    /*
     * update the residual ( offline step )
     * \param X : solution
     * \param R : residual
     * \param mu : current parameter
     */
    void offlineUpdateResidual( const vector_ptrtype & X,  vector_ptrtype& R , const parameter_type & mu) ;

    /*
     * update the residual
     * \param map_X : solution
     * \param map_R : residual
     * \param mu : current parameter
     * \param N : dimension of the reduced basis
     */
    virtual void updateResidual( const map_dense_vector_type& map_X, map_dense_vector_type& map_R , const parameter_type& mu, int N ) const ;

    /*
     * compute the projection of the initial guess
     * \param mu : current parameter
     * \param N : dimension of the reduced basis
     * \param initial guess
     */
    void computeProjectionInitialGuess( const parameter_type & mu, int N , vectorN_type& initial_guess ) const ;

    element_type femSolve( parameter_type& mu );

    virtual void offlineSolve( element_type& u, element_type& udu, parameter_type& mu, element_ptrtype & dual_initial_field );

    /*
     * newton for primal problem ( offline step )
     * \param mu : current parameter
     */
    element_type offlineNewtonPrimal( parameter_type const& mu ) ;

    /*
     * newton algorithm to solve non linear problem
     * \param N : dimension of the reduced basis
     * \param mu : current parameter
     * \param uN : reduced solution ( vectorN_type )
     * \param condition number of the reduced jacobian
     */
    matrix_info_tuple newton(  size_type N, parameter_type const& mu , vectorN_type & uN, double& output) const ;

    /*
     * fixed point for primal problem ( offline step )
     * \param mu : current parameter
     */
    element_type offlineFixedPointPrimal( parameter_type const& mu );//, sparse_matrix_ptrtype & A ) ;

    /*
     * \param mu : current parameter
     * \param dual_initial_field : to be filled
     */
    element_type offlineFixedPointDual( parameter_type const& mu , element_ptrtype & dual_initial_field );//, const sparse_matrix_ptrtype & A, const element_type & u ) ;

    /*
     * fixed point ( primal problem ) - ONLINE step
     * \param N : dimension of the reduced basis
     * \param mu :current parameter
     * \param uN : dual reduced solution ( vectorN_type )
     * \param uNold : dual old reduced solution ( vectorN_type )
     * \param condition number of the reduced matrix A
     * \param output : vector of outpus at each time step
     * \param K : number of time step ( default value, must be >0 if used )
     * \param computeOutput : if true compute quantity of interest (output)
     */
    virtual matrix_info_tuple fixedPointPrimal( size_type N, parameter_type const& mu, std::vector< vectorN_type > & uN,  std::vector<vectorN_type> & uNold,
                                                std::vector< double > & output_vector, int K=0, bool print_rb_matrix=false, bool computeOutput=true ) const;


    /**
     * solve fem model using the affine decomposition
     * \param mu :current parameter
     */
    element_type solveFemModelUsingAffineDecomposition( parameter_type const& mu ) const override
        {
            if ( M_use_newton )
                return M_model->solveFemUsingAffineDecompositionNewton( mu );
            else
                return M_model->solveFemUsingAffineDecompositionFixedPoint( mu );
        }


    /*
     * Dump data array into a file
     * \param out : Name of the output file
     * \param prefix : Prefix used for current data
     * \param array : Array to dump
     * \param nbelem : Number of elements to dump
     */
    void dumpData(std::string out, std::string prefix, const double * array, int nbelem) const ;

#if defined(FEELPP_HAS_HARTS) && defined(HARTS_HAS_OPENCL)
    /*
     * fixed point ( primal problem ) - ONLINE step with OpenCL
     * \param N : dimension of the reduced basis
     * \param mu :current parameter
     * \param uN : dual reduced solution ( vectorN_type )
     * \param uNold : dual old reduced solution ( vectorN_type )
     * \param condition number of the reduced matrix A
     * \param output : vector of outpus at each time step
     * \param K : number of time step ( default value, must be >0 if used )
     */
    matrix_info_tuple fixedPointPrimalCL( size_type N, parameter_type const& mu, std::vector< vectorN_type > & uN,  std::vector<vectorN_type> & uNold,
                                          std::vector< double > & output_vector, int K=0, bool print_rb_matrix=false) const;
#endif

    /*
     * fixed point ( dual problem ) - ONLINE step
     * \param N : dimension of the reduced basis
     * \param mu :current parameter
     * \param uN : primal reduced solution ( vectorN_type )
     * \param uNdu : dual reduced solution ( vectorN_type )
     * \param uNduold : dual old reduced solution ( vectorN_type )
     * \param output : vector of outpus at each time step
     * \param K : number of time step ( default value, must be >0 if used )
     */
    virtual void fixedPointDual(  size_type N, parameter_type const& mu, std::vector< vectorN_type > const& uN,
                          std::vector< vectorN_type > & uNdu,  std::vector<vectorN_type> & uNduold, std::vector< double > & output_vector, int K=0) const;

    /**
     * fixed point ( main ) - ONLINE step
     * \param N : dimension of the reduced basis
     * \param mu :current parameter
     * \param uN : dual reduced solution ( vectorN_type )
     * \param uNdu : dual reduced solution ( vectorN_type )
     * \param uNold : dual old reduced solution ( vectorN_type )
     * \param uNduold : dual old reduced solution ( vectorN_type )
     * \param output : vector of outpus at each time step
     * \param K : number of time step ( default value, must be >0 if used )
     * \param computeOutput : if true compute quantity of interest (output)
     */
    matrix_info_tuple fixedPoint(  size_type N, parameter_type const& mu, std::vector< vectorN_type > & uN, std::vector< vectorN_type > & uNdu,
                                   std::vector<vectorN_type> & uNold, std::vector<vectorN_type> & uNduold,
                                   std::vector< double > & output_vector, int K=0, bool print_rb_matrix=false, bool computeOutput = true ) const;

    virtual matrix_info_tuple onlineSolve(  size_type N, parameter_type const& mu, std::vector< vectorN_type > & uN, std::vector< vectorN_type > & uNdu, std::vector<vectorN_type> & uNold, std::vector<vectorN_type> & uNduold, std::vector< double > & output_vector, int K=0, bool print_rb_matrix=false, bool computeOutput = true ) const;

    /**
     * computation of the conditioning number of a given matrix
     * \param A : reduced matrix
     */
    double computeConditioning( matrixN_type & A ) const ;


    /**
     * Returns the lower bound of the output
     *
     * \param mu \f$ \mu\f$ the parameter at which to evaluate the output
     * \param N the size of the reduced basis space to use
     * \param uN primal solution
     * \param uNdu dual solution
     *
     *\return compute online the lower bound and condition number of matrix A
     */
    boost::tuple<std::vector<double>,double> lb( parameter_ptrtype const& mu, size_type N, std::vector< vectorN_type >& uN, std::vector< vectorN_type >& uNdu ) const
        {
            return lb( N, *mu, uN, uNdu );
        }

    /**
     * Returns the upper bound of the output associed to \f$\mu\f$
     *
     * \param mu \f$ \mu\f$ the parameter at which to evaluate the output
     * \param N the dimension of \f$W_N\f$
     * \param uN primal solution
     * \param uNdu dual solution
     *
     *\return compute online the lower bound
     */


    value_type ub( size_type N, parameter_type const& mu, std::vector< vectorN_type >& uN, std::vector< vectorN_type >& uNdu ) const
        {
            auto o = lb( N, mu, uN, uNdu );
            auto e = delta( N, mu, uN, uNdu );
            auto output_vector=o.template get<0>();
            double output_vector_size=output_vector.size();
            double output = output_vector[output_vector_size-1];
            return output + e.template get<0>();
        }

    /**
     * Returns the error bound on the output
     *
     * \param mu \f$ \mu\f$ the parameter at which to evaluate the output
     * \param N the size of the reduced basis space to use
     * \param uN primal solution
     * \param uNdu dual solution
     *
     *\return compute online the lower bound
     */
    value_type delta( size_type N, parameter_ptrtype const& mu, std::vector< vectorN_type > const& uN, std::vector< vectorN_type > const& uNdu, std::vector<vectorN_type> const& uNold,  std::vector<vectorN_type> const& uNduold , int k=0 ) const
        {
            auto e = delta( N, mu, uN, uNdu );
            return e.template get<0>();
        }

    /**
     * Returns the error bound on the output associed to \f$\mu\f$
     *
     * \param mu \f$ \mu\f$ the parameter at which to evaluate the output
     * \param N the dimension of \f$W_N\f$
     * \param uN primal solution
     * \param uNdu dual solution
     *
     *\return compute online the lower bound
     */
    virtual error_estimation_type delta( size_type N, parameter_type const& mu, std::vector< vectorN_type > const& uN, std::vector< vectorN_type > const& uNdu, std::vector<vectorN_type> const& uNold, std::vector<vectorN_type> const& uNduold, int k=0 ) const;

    /**
     * Returns the upper bound of the output associed to \f$\mu\f$
     *
     * \param mu \f$ \mu\f$ the parameter at which to evaluate the output
     * \param N the dimension of \f$W_N\f$
     *
     *\return compute online the lower bound
     */
    value_type ub( size_type K, parameter_ptrtype const& mu, std::vector< vectorN_type >& uN, std::vector< vectorN_type >& uNdu ) const
        {
            return ub( K, *mu, uN, uNdu );
        }


    vector_ptrtype computeRieszResidual( parameter_type const& mu, std::vector<vectorN_type> const& uN ) const;
    double computeRieszResidualNorm( parameter_type const& mu ) const override;
    double computeRieszResidualNorm( parameter_type const& mu, std::vector<vectorN_type> const& uN ) const override;

    /**
     * Offline computation
     *
     * \return the convergence history (max error)
     */
    virtual convergence_type offline();

    /**
     * \brief update size of Affine Decomposition (needed by co-build since size of eim decomposition can change)
     */
    virtual void updateAffineDecompositionSize();

    /**
     * \brief Retuns maximum value of the relative error
     * \param N number of elements in the reduced basis <=> M_N
     */
    virtual max_error_type maxErrorBounds( size_type N ) const;

    /**
     * evaluate online the residual
     */

    residual_error_type transientPrimalResidual( int Ncur, parameter_type const& mu,  vectorN_type const& Un, vectorN_type const& Unold=vectorN_type(), double time_step=1, double time=1e30 ) const;
    residual_error_type transientPrimalResidualEim( int Ncur, parameter_type const& mu,  vectorN_type const& Un, vectorN_type const& Unold=vectorN_type(), double time_step=1, double time=1e30 ) const;
    residual_error_type steadyPrimalResidual( int Ncur, parameter_type const& mu,  vectorN_type const& Un, double time=0 ) const;
    residual_error_type steadyPrimalResidualEim( int Ncur, parameter_type const& mu,  vectorN_type const& Un, double time=0 ) const;

    residual_error_type transientDualResidual( int Ncur, parameter_type const& mu,  vectorN_type const& Un, vectorN_type const& Unold=vectorN_type(), double time_step=1, double time=1e30 ) const;
    residual_error_type transientDualResidualEim( int Ncur, parameter_type const& mu,  vectorN_type const& Un, vectorN_type const& Unold=vectorN_type(), double time_step=1, double time=1e30 ) const;
    residual_error_type steadyDualResidual( int Ncur, parameter_type const& mu,  vectorN_type const& Un, double time=0 ) const;
    residual_error_type steadyDualResidualEim( int Ncur, parameter_type const& mu,  vectorN_type const& Un, double time=0 ) const;


    value_type initialDualResidual( int Ncur, parameter_type const& mu, vectorN_type const& Uduini, double time_step ) const ;


    /**
     * generate offline the residual
     */
    virtual void offlineResidual( int Ncur , int number_of_added_elements=1 );
    void offlineResidual( int Ncur, mpl::bool_<true> ,int number_of_added_elements=1 );
    void offlineResidual( int Ncur, mpl::bool_<false> , int number_of_added_elements=1 );
    void offlineResidualV0( int Ncur, mpl::bool_<false> , int number_of_added_elements=1 );
    void offlineResidualV1( int Ncur, mpl::bool_<false> , int number_of_added_elements=1 );

    void offlineResidualEim( int Ncur , int number_of_added_elements=1 );
    void offlineResidualEim( int Ncur, mpl::bool_<true> ,int number_of_added_elements=1 );
    void offlineResidualEim( int Ncur, mpl::bool_<false> , int number_of_added_elements=1 );

    /*
     * compute empirical error estimation, ie : |S_n - S{n-1}|
     * \param Ncur : number of elements in the reduced basis <=> M_N
     * \param mu : parameters value (one value per parameter)
     * \output : empirical error estimation
     */
    value_type empiricalErrorEstimation ( int Ncur, parameter_type const& mu, int k ) const ;


    /**
     * return the crb expansion at parameter \p \mu, ie \f$\sum_{i=0}^N u^N_i
     * \phi_i\f$ where $\phi_i, i=1...N$ are the basis function of the reduced
     * basis space
     * if N>0 take the N^th first elements, else take all elements
     */
    virtual element_type runWithExpansion( parameter_type const& mu , int N=-1, int time_index=-1);

    /**
     * return the crb expansion at parameter \p \mu, ie \f$\sum_{i=0}^N u^N_i
     * \phi_i\f$ where $\phi_i, i=1...N$ are the basis function of the reduced
     * basis space
     */
    element_type expansion( vectorN_type const& u , int N = -1, bool dual=false ) const override;

    // Summary of number of iterations (at the current step)
    std::pair<int,double> online_iterations() override {return online_iterations_summary;}
    std::pair<int,double> offline_iterations(){return offline_iterations_summary;}

    void checkInitialGuess( const element_type expansion_uN , parameter_type const& mu, vectorN_type & error ) const ;
    void checkInitialGuess( const element_type expansion_uN , parameter_type const& mu, vectorN_type & error , mpl::bool_<true> ) const ;
    void checkInitialGuess( const element_type expansion_uN , parameter_type const& mu, vectorN_type & error , mpl::bool_<false> ) const ;

    /**
     * run the certified reduced basis with P parameters and returns 1 output
     */
    //boost::tuple<double,double,double> run( parameter_type const& mu, double eps = 1e-6 );
    //boost::tuple<double,double,double,double> run( parameter_type const& mu, double eps = 1e-6 );
    //by default N=-1 so we take dimension-max but if N>0 then we take N basis functions toperform online step
    //boost::tuple<std::vector<double>,double, solutions_tuple, matrix_info_tuple, double, double, upper_bounds_tuple >
    CRBResults
    run( parameter_type const& mu,
         vectorN_type & time,
         double eps = 1e-6,
         int N = -1,
         bool print_rb_matrix=false );

    CRBResults
    run( parameter_type const& mu,
         double eps = 1e-6,
         int N = -1,
         bool print_rb_matrix=false )
        {
            vectorN_type times;
            return this->run( mu, times, eps, N, print_rb_matrix );
        }
    std::vector<CRBResults>
    run( std::vector<parameter_type> const& S,
         double eps = 1e-6,
         int N = -1,
         bool print_rb_matrix=false )
        {
            std::vector<CRBResults> res;
            res.reserve( S.size() );
            for( auto const& mu : S )
            {
                auto r = this->run( mu, eps, N, print_rb_matrix );
                res.push_back( r );
            }
            return res;
        }

    /**
     * run the certified reduced basis with P parameters and returns 1 output
     */
    void run( const double * X, unsigned long N, double * Y, unsigned long P );
    void run( const double * X, unsigned long N, double * Y, unsigned long P, mpl::bool_<true> );
    void run( const double * X, unsigned long N, double * Y, unsigned long P, mpl::bool_<false> );

    /**
     * \return a random sampling
     */
    sampling_type randomSampling( int N )
        {
            bool same_sampling_on_all_proc=false;
            M_Xi->randomize( N , same_sampling_on_all_proc );
            return *M_Xi;
        }

    /**
     * \return a equidistributed sampling
     */
    sampling_type equidistributedSampling( int N )
        {
            bool same_sampling_on_all_proc=false;
            M_Xi->equidistribute( N , same_sampling_on_all_proc );
            return *M_Xi;
        }

    sampling_ptrtype const& wnmu ( ) const
        {
            return M_WNmu;
        }

    wn_type const& wn() const
        {
            //return M_WN;
            return M_model->rBFunctionSpace()->primalRB();
        }

    wn_type const& wndu() const
        {
            //return M_WNdu;
            return M_model->rBFunctionSpace()->dualRB();
        }

    //!
    //! load RB database given crb.json file \p filename and load strategy \p l
    //! @param filename filename (with crb.json extension) including path
    //! @param l loading strategy (reduced basis, finite element, all)
    //!
    void loadDB( std::string const& filename, crb::load l ) override;


    /**
     * \brief load Crb from json
     * \param filename : input json filename
     * \param loadingContext : 0 minimal crb online run, 1 with fe reduced basis
     */
    void loadJson( std::string const& filename, size_type loadingContext = 0 );

    void saveJson();
    /**
     * \brief setup Crb from property_tree::ptree
     * \param ptree : input property_tree::ptree
     * \param loadingContext : 0 minimal crb online run, 1 with fe reduced basis
     */
    void setup( boost::property_tree::ptree const& ptree, size_type loadingContext, std::string const& dbDir );

    /**
     * save the CRB database
     */
    void saveDB() override;

    /**
     * load the CRB database
     */
    bool loadDB() override;

    /**
     *  do the projection on the POD space of u (for transient problems)
     *  \param u : the solution to project (input parameter)
     *  \param projection : the projection (output parameter)
     *  \param name_of_space : primal or dual
     */
    void projectionOnPodSpace( const element_type & u , element_ptrtype& projection ,const std::string& name_of_space="primal" );


    bool useWNmu() const
        {
            bool use = boption(_prefix=M_prefix,_name="crb.run-on-WNmu");
            return use;
        }


    /**
     * if true, rebuild the database (if already exist)
     */
    bool rebuild() const { return M_rebuild; }

    void setRebuild( bool b ) { M_rebuild = b; }

    void setAdaptationSER( bool b ) const { M_SER_adapt = b; }
    bool adaptationSER() const { return M_SER_adapt; }
    /**
     * if true, show the mu selected during the offline stage
     */
    bool showMuSelection() ;

    /**
     * if true, print the max error (absolute) during the offline stage
     */
    bool printErrorDuringOfflineStep();

    /**
     * print parameters set mu selected during the offline stage
     */
    void printMuSelection( void );

    /**
     * print informations on offline Picard iterations
     */
    void printRbPicardIterations();
    /**
     * print informations on online Picard iterations
     */
    void printOnlineRbPicardIterations(parameter_type const& mu) const;

    /**
     * print max errors (total error and also primal and dual contributions)
     * during offline stage
     */
    void printErrorsDuringRbConstruction( void );
    /*
     * compute correction terms for output
     * \param mu \f$ \mu\f$ the parameter at which to evaluate the output
     * \param uN : primal solution
     * \param uNdu : dual solution
     * \pram uNold : old primal solution
     * \param K : time index
     */
    virtual double correctionTerms(parameter_type const& mu, std::vector< vectorN_type > const & uN, std::vector< vectorN_type > const & uNdu , std::vector<vectorN_type> const & uNold,  int const K=0) const;

    /*
     * build matrix to store functions used to compute the variance output
     * \param N : number of elements in the reduced basis
     */
    void buildVarianceMatrixPhi(int const N);
    void buildVarianceMatrixPhi(int const N , mpl::bool_<true> );
    void buildVarianceMatrixPhi(int const N , mpl::bool_<false> );

    void printRBMatrix(matrixN_type const& A, parameter_type const& mu) const
        {
            int rows=A.rows();
            int cols=A.cols();
            int N = rows;
            std::string file_name = fmt::format("RBmatrix-N{}-mu-{}", N, mu[0]);
            if( this->worldComm().isMasterRank() )
            {
                std::ofstream file;
                file.open( file_name,std::ios::out );
                file<<"# name: A\n";
                file<<"# type: matrix\n";
                file<<"# rows: "<<rows<<"\n";
                file<<"# columns: "<<cols<<"\n";
                for(int i=0; i<rows; i++)
                {
                    for(int j=0; j<cols; j++)
                        file<<std::setprecision(15)<<A(i,j)<<" ";
                    file<<"\n";
                }
            }
        }

    /**
     * evaluate online time via the option crb.computational-time-neval
     */
    void computationalTimeStatistics( std::string appname );

    /**
     * update basis built to have error estimators from F.Casenave's paper
     * \param N : current number of elements in the RB
     * \param new parameters that we must select from N-1 to N
     * Be carreful, M_primal_apee_basis and M_dual_apee_basis don't have the same size because affine decomposition can be different
     */
    void updateApeeBasisConstruction( int N , std::vector< parameter_type > const & new_primal_parameters, std::vector< parameter_type > const & new_dual_parameters) const;
    void updateApeeOfflineResidualComputation( int N , std::vector< parameter_type > const & new_primal_parameters, std::vector< parameter_type > const & new_dual_parameters ) ;
    void computePrimalApeeBasis( parameter_type const& mu)  const ;
    void computeDualApeeBasis( parameter_type const& mu)  const ;
    double computeSquareDualNormOfPrimalResidual( parameter_type const& mu, element_type const & u);
    double computeSquareDualNormOfDualResidual( parameter_type const& mu, element_type const & udu);
    //compute X-vector in F.Casenave's paper
    vectorN_type computeOnlinePrimalApeeVector( parameter_type const& mu , vectorN_type const & uN, vectorN_type const & uNold=vectorN_type(), double dt=1e30, double time=1e30 ) const;
    vectorN_type computeOnlineDualApeeVector( parameter_type const& mu , vectorN_type const & uNdu, vectorN_type const & uNduold=vectorN_type(), double dt=1e30, double time=1e30 ) const ;

    void assemblePrimalAndDualApeeBasisMatrices( int N ) const ;

    /**
     * Greedy algorithm for apee basis
     * \param N : current number of elements in the RB
     * \param new_parameters : vector to be enrich by new parameters associated to N
     */
    void selectPrimalApeeParameters( int N, std::vector< parameter_type > & new_parameters );
    void selectDualApeeParameters( int N, std::vector< parameter_type > & new_parameters );

    double computeOnlinePrimalApee( int N , parameter_type const& mu , vectorN_type const & uN, vectorN_type const & uNold=vectorN_type(), double dt=1e30, double time=1e30 ) const ;
    double computeOnlineDualApee( int N , parameter_type const& mu , vectorN_type const & uNdu, vectorN_type const & uNduold=vectorN_type(), double dt=1e30, double time=1e30 ) const ;
    //@}


    virtual void updateRbInDeim()
        {
            M_model->updateRbInDeim();
        }

protected:
    /**
     * generate the super sampling M_Xi depending of the options
     **/
    void generateSuperSampling();
    bool buildSampling();
    virtual void addBasis( element_type& u, element_type& udu, parameter_type& mu );
    virtual void orthonormalizeBasis( int number_of_added_elements );
    virtual void buildRbMatrix( int number_of_added_elements, parameter_type& mu, element_ptrtype dual_initial_field );
    virtual void buildRbMatrixTrilinear( int number_of_added_elements, parameter_type& mu )
        {}

    virtual void saveRB();

    std::string M_prefix;

    crb_elements_db_type M_elements_database;

    std::shared_ptr<SolverNonLinear<double> > M_nlsolver;

    truth_model_ptrtype M_model;

    backend_ptrtype M_backend;
    backend_ptrtype M_backend_primal;
    backend_ptrtype M_backend_dual;

    int M_output_index;
    std::string M_output_name;

    double M_tolerance;

    size_type M_iter_max;

    int M_factor;

    CRBErrorType M_error_type;

    double M_maxerror;

    // parameter space
    parameterspace_ptrtype M_Dmu;

    // fine sampling of the parameter space
    sampling_ptrtype M_Xi;

    // sampling of parameter space to build WN
    sampling_ptrtype M_WNmu;
    sampling_ptrtype M_WNmu_complement;

    //sampling of parameter space to build a posteriori error estimators (paper of F.Casenave)
    //Accurate a posteriori error evaluation in the reduced basis method - 2012
    sampling_ptrtype M_primal_apee_mu;
    sampling_ptrtype M_dual_apee_mu;

    //scm
    scm_ptrtype M_scmA;
    scm_ptrtype M_scmM;

    //export
    mutable export_ptrtype M_exporter;

#if 0
    array_2_type M_C0_pr;
    array_2_type M_C0_du;
    array_3_type M_Lambda_pr;
    array_3_type M_Lambda_du;
    array_4_type M_Gamma_pr;
    array_4_type M_Gamma_du;
    array_3_type M_Cmf_pr;
    array_3_type M_Cmf_du;
    array_4_type M_Cma_pr;
    array_4_type M_Cma_du;
    array_4_type M_Cmm_pr;
    array_4_type M_Cmm_du;
#endif
    std::vector< std::vector< std::vector< std::vector< double > > > > M_C0_pr;
    std::vector< std::vector< std::vector< std::vector< double > > > > M_C0_du;
    std::vector< std::vector< std::vector< std::vector< vectorN_type > > > > M_Lambda_pr;
    std::vector< std::vector< std::vector< std::vector< vectorN_type > > > > M_Lambda_du;
    std::vector< std::vector< std::vector< std::vector< matrixN_type > > > > M_Gamma_pr;
    std::vector< std::vector< std::vector< std::vector< matrixN_type > > > > M_Gamma_du;
    std::vector< std::vector< std::vector< vector_ptrtype > > > M_precomputeResidualPrA;
    std::vector< std::vector< std::vector< vector_ptrtype > > > M_precomputeResidualDuA;
    std::vector< std::vector< vector_ptrtype > > M_precomputeResidualPrF;
    std::vector< std::vector< vector_ptrtype > > M_precomputeResidualDuF;

    std::vector< std::vector< double > > M_C0_pr_eim;
    std::vector< std::vector< double > > M_C0_du_eim;
    std::vector< std::vector< vectorN_type > > M_Lambda_pr_eim;
    std::vector< std::vector< vectorN_type > > M_Lambda_du_eim;
    std::vector< std::vector< matrixN_type > > M_Gamma_pr_eim;
    std::vector< std::vector< matrixN_type > > M_Gamma_du_eim;
    std::vector< std::vector< std::vector< std::vector< vectorN_type > > > > M_Cmf_pr;
    std::vector< std::vector< std::vector< std::vector< vectorN_type > > > > M_Cmf_du;
    std::vector< std::vector< std::vector< std::vector< vectorN_type > > > > M_Cmf_du_ini;
    std::vector< std::vector< std::vector< std::vector< matrixN_type > > > > M_Cma_pr;
    std::vector< std::vector< std::vector< std::vector< matrixN_type > > > > M_Cma_du;
    std::vector< std::vector< std::vector< std::vector< matrixN_type > > > > M_Cmm_pr;
    std::vector< std::vector< std::vector< std::vector< matrixN_type > > > > M_Cmm_du;
    std::vector< std::vector< vectorN_type > > M_Cmf_pr_eim;
    std::vector< std::vector< vectorN_type > > M_Cmf_du_eim;
    std::vector< std::vector< vectorN_type > > M_Cmf_du_ini_eim;
    std::vector< std::vector< matrixN_type > > M_Cma_pr_eim;
    std::vector< std::vector< matrixN_type > > M_Cma_du_eim;
    std::vector< std::vector< matrixN_type > > M_Cmm_pr_eim;
    std::vector< std::vector< matrixN_type > > M_Cmm_du_eim;

    //X( \mu_r ) in F.casnave's paper
    mutable std::vector< vectorN_type > M_primal_apee_basis;
    mutable std::vector< vectorN_type > M_dual_apee_basis;
    //V_r in F.casnave's paper ( contains square of dual norm of residuals )
    vectorN_type M_primal_V;
    vectorN_type M_dual_V;
    //matrices who have X( \mu_r ) as columns
    mutable matrixN_type M_primal_T;
    mutable matrixN_type M_dual_T;

    vectorN_type M_coeff_pr_ini_online;
    vectorN_type M_coeff_du_ini_online;

    // inner_product operator
    matrixN_type M_algebraicInnerProductPrimal;


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
    //wn_type M_WN;
    //wn_type M_WNdu;

    size_type M_N;
    size_type M_Nm;
    size_type number_of_added_elements;

    bool M_solve_dual_problem;
    bool M_orthonormalize_primal;
    bool M_orthonormalize_dual;

    convergence_type M_rbconv;

    //time
    bdf_ptrtype M_bdf_primal;
    bdf_ptrtype M_bdf_primal_save;
    bdf_ptrtype M_bdf_dual;
    bdf_ptrtype M_bdf_dual_save;

    // left hand side
    std::vector < std::vector<matrixN_type> > M_Aqm_pr;
    std::vector < std::vector< std::vector<vector_ptrtype> > > M_hAqm;
    std::vector < std::vector<matrixN_type> > M_Aqm_du;
    std::vector < std::vector<matrixN_type> > M_Aqm_pr_du;

    //jacobian
    std::vector < std::vector<matrixN_type> > M_Jqm_pr;
    //residual
    std::vector < std::vector<vectorN_type> > M_Rqm_pr;

    //mass matrix
    std::vector < std::vector<matrixN_type> > M_Mqm_pr;
    std::vector < std::vector<matrixN_type> > M_Mqm_du;
    std::vector < std::vector<matrixN_type> > M_Mqm_pr_du;

    // right hand side
    std::vector < std::vector<vectorN_type> > M_Fqm_pr;
    std::vector < std::vector<vector_ptrtype> > M_hFqm;
    std::vector < std::vector<vectorN_type> > M_Fqm_du;
    // output
    std::vector < std::vector<vectorN_type> > M_Lqm_pr;
    std::vector < std::vector<vectorN_type> > M_Lqm_du;

    //initial guess
    std::vector < std::vector<vectorN_type> > M_InitialGuessV_pr;

    //std::vector<int> M_index;
    int M_mode_number;

    std::vector < matrixN_type > M_variance_matrix_phi;

    bool M_compute_variance;
    bool M_computeMatrixInfo;
    //bool M_rbconv_contains_primal_and_dual_contributions;
    parameter_type M_current_mu;
    int M_no_residual_index;

    bool M_database_contains_variance_info;

    bool M_use_newton;

    bool M_offline_step;

    bool M_rebuild;
    //int M_SER_groupsize;
    mutable bool M_SER_adapt;
    mutable double M_SER_maxerr;
    bool M_SER_errorEstimation;
    bool M_SER_useGreedyInRb;

    preconditioner_ptrtype M_preconditioner;
    preconditioner_ptrtype M_preconditioner_primal;
    preconditioner_ptrtype M_preconditioner_dual;

    //true if the model is executed in steady mode
    bool M_model_executed_in_steady_mode;

    std::vector< std::vector<element_ptrtype> > M_InitialGuessV;
    std::vector< std::vector<sparse_matrix_ptrtype> > M_Aqm;
    std::vector< std::vector<std::vector<vector_ptrtype> > > M_Fqm;
    std::vector< std::vector<sparse_matrix_ptrtype> > M_Jqm;
    std::vector< std::vector< std::vector<vector_ptrtype> > > M_Rqm;
    std::vector< std::vector<std::vector<vector_ptrtype> > > M_Lqm;
    std::vector< std::vector<sparse_matrix_ptrtype> > M_Mqm;

    bool M_save_output_behavior;

    int M_fixedpointMaxIterations;
    double M_fixedpointIncrementTol;
    bool M_fixedpointVerbose;
    double M_fixedpointCriticalValue;
    bool M_fixedpointUseAitken;

    bool M_loadElementsDb;
    bool M_useRbSpaceContextEim;
    bool M_hasRbSpaceContextEim;

    bool M_useAccurateApee;
    bool M_computeApeeForEachTimeStep;
    bool M_seekMuInComplement;
    bool M_showResidual;

    bool M_check_cvg;
    mutable bool M_last_online_converged;

    mutable std::pair<int,double> offline_iterations_summary;
    mutable std::pair<int,double> online_iterations_summary;
private:
    //!
    //! setup CRB from DB
    //!
    void setupOfflineFromDB();
};



template<typename TruthModelType>
typename CRB<TruthModelType>::element_type
CRB<TruthModelType>::offlineFixedPointPrimal(parameter_type const& mu )//, sparse_matrix_ptrtype & A )
{
    if ( M_model->isSteady() )
        return M_model->solveFemUsingAffineDecompositionFixedPoint( mu );

    auto u = M_model->functionSpace()->element();

    sparse_matrix_ptrtype M = M_model->newMatrix();
    int nl = M_model->Nl();  //number of outputs
    std::vector< vector_ptrtype > F( nl );
    for(int l=0; l<nl; l++)
        F[l]=M_model->newVector();

    //M_backend_primal = backend_type::build( BACKEND_PETSC );
    //bool reuse_prec = boption(_prefix=M_prefix,_name="crb.reuse-prec") ;

    M_bdf_primal = bdf( _space=M_model->functionSpace(), _vm=Environment::vm() , _name="bdf_primal" );
    M_bdf_primal_save = bdf( _space=M_model->functionSpace(), _vm=Environment::vm() , _name="bdf_primal_save" );

    //set parameters for time discretization
    M_bdf_primal->setTimeInitial( M_model->timeInitial() );
    M_bdf_primal->setTimeStep( M_model->timeStep() );
    M_bdf_primal->setTimeFinal( M_model->timeFinal() );
    M_bdf_primal->setOrder( M_model->timeOrder() );

    M_bdf_primal_save->setTimeInitial( M_model->timeInitial() );
    M_bdf_primal_save->setTimeStep( M_model->timeStep() );
    M_bdf_primal_save->setTimeFinal( M_model->timeFinal() );
    M_bdf_primal_save->setOrder( M_model->timeOrder() );

    M_bdf_primal_save->setRankProcInNameOfFiles( true );
    M_bdf_primal->setRankProcInNameOfFiles( true );

    //initialization of unknown
    auto elementptr = M_model->functionSpace()->elementPtr();
    M_model->initializationField( elementptr, mu );

    u = *elementptr;

    auto Apr = M_model->newMatrix();

    int max_fixedpoint_iterations  = ioption(_prefix=M_prefix,_name="crb.fixedpoint.maxit");
    double increment_fixedpoint_tol  = doption(_prefix=M_prefix,_name="crb.fixedpoint.increment-tol");
    double fixedpoint_critical_value  = doption(_prefix=M_prefix,_name="crb.fixedpoint.critical-value");
    int iteration=0;
    double increment_norm=1e3;
    bool is_linear=M_model->isLinear();

    if( is_linear )
        increment_norm = 0;

    double bdf_coeff;

    auto vec_bdf_poly = M_backend_primal->newVector( M_model->functionSpace() );

    //assemble the initial guess for the given mu
    if ( M_model->isSteady() && ! is_linear )
    {
        elementptr = M_model->assembleInitialGuess( mu ) ;
        u = *elementptr ;
    }
    auto uold = M_model->functionSpace()->element();
    auto bdf_poly = M_model->functionSpace()->element();


    vector_ptrtype Rhs( M_backend_primal->newVector( M_model->functionSpace() ) );

    bool POD_WN = boption(_prefix=M_prefix,_name="crb.apply-POD-to-WN") ;

    for ( M_bdf_primal->start(u),M_bdf_primal_save->start(u);
          !M_bdf_primal->isFinished() && !M_bdf_primal_save->isFinished();
          M_bdf_primal->next(u) , M_bdf_primal_save->next() )
    {
        int bdf_iter = M_bdf_primal->iteration();

        if ( ! M_model->isSteady() )
        {
            bdf_coeff = M_bdf_primal->polyDerivCoefficient( 0 );
            bdf_poly = M_bdf_primal->polyDeriv();
        }

        do
        {
            if( is_linear )
            {
                bool compute_only_terms_time_dependent=false;
                if ( bdf_iter == 1 )
                {
                    boost::tie( M, Apr, F) = M_model->update( mu , M_bdf_primal->time() , compute_only_terms_time_dependent );
                }
                else
                {
                    compute_only_terms_time_dependent=true;
                    boost::tie( boost::tuples::ignore , boost::tuples::ignore , F) = M_model->update( mu , M_bdf_primal->time() , compute_only_terms_time_dependent );
                }
            }
            else
            {
                bool compute_only_terms_time_dependent=false;
                if ( bdf_iter == 1 )
                {
                    boost::tie( M, Apr, F) = M_model->update( mu , u, M_bdf_primal->time() , compute_only_terms_time_dependent );
                }
                else
                {
                    compute_only_terms_time_dependent=true;
                    boost::tie( boost::tuples::ignore, boost::tuples::ignore, F) = M_model->update( mu , u, M_bdf_primal->time() , compute_only_terms_time_dependent );
                }
            }

            if ( ! M_model->isSteady() )
            {
                bdf_coeff = M_bdf_primal->polyDerivCoefficient( 0 );
                if ( bdf_iter == 1 )
                {
                    Apr->addMatrix( bdf_coeff, M );
                }
                //auto bdf_poly = M_bdf_primal->polyDeriv();
                *Rhs = *F[0];
                *vec_bdf_poly = bdf_poly;
                Rhs->addVector( *vec_bdf_poly, *M );
            }
            else
            {
                *Rhs = *F[0];
            }
            //backup for non linear problems
            uold = u;

            //solve
            if( boption(_prefix=M_prefix,_name="crb.use-primal-pc") )
            {
                M_preconditioner_primal->setMatrix( Apr );
                auto ret = M_backend_primal->solve( _matrix=Apr, _solution=u, _rhs=Rhs,  _prec=M_preconditioner_primal, _reuse_prec=( bdf_iter >= 2 ) );
                if  ( !ret.template get<0>() )
                    LOG(INFO)<<"[CRB] WARNING : at time "<<M_bdf_primal->time()<<" we have not converged ( nb_it : "<<ret.template get<1>()<<" and residual : "<<ret.template get<2>() <<" ) \n";
            }
            else
            {
                auto ret = M_backend_primal->solve( _matrix=Apr, _solution=u, _rhs=Rhs, _rebuild=true);
                if  ( !ret.template get<0>() )
                    LOG(INFO)<<"[CRB] WARNING : at time "<<M_bdf_primal->time()<<" we have not converged ( nb_it : "<<ret.template get<1>()<<" and residual : "<<ret.template get<2>() <<" ) \n";
            }

            //on each subspace the norme of the increment is computed and then we perform the sum
            if( is_linear )
                increment_norm = 0;
            else
                increment_norm = M_model->computeNormL2( u , uold );

            LOG(INFO) << "iteration " << iteration << ", increment_norm = " <<  increment_norm << "\n";
            this->offline_iterations_summary.first = iteration;
            this->offline_iterations_summary.second = increment_norm;
            iteration++;

        }while( increment_norm > increment_fixedpoint_tol && iteration < max_fixedpoint_iterations );

        if ( ! M_model->isSteady() )
        {
            u.close();
            if(POD_WN )
            {
                M_bdf_primal_save->shiftRight( u );
            }
            else
            {
                element_ptrtype projection ( new element_type ( M_model->functionSpace() ) );
                projectionOnPodSpace ( u , projection, "primal" );
                auto e = u;
                e-=*projection;
                M_bdf_primal_save->shiftRight( e );
            }
        }

        if( increment_norm > fixedpoint_critical_value )
            throw std::logic_error( fmt::format("[CRB::offlineFixedPointPrimal]  at time {} ERROR : increment > critical value ", M_bdf_primal->time() ) );

        for ( size_type l = 0; l < M_model->Nl(); ++l )
        {
            F[l]->close();
            element_ptrtype eltF( new element_type( M_model->functionSpace() ) );
            *eltF = *F[l];
            LOG(INFO) << "u^T F[" << l << "]= " << inner_product( u, *eltF ) << " at time : "<<M_bdf_primal->time()<<"\n";
        }
        LOG(INFO) << "[CRB::offlineWithErrorEstimation] energy = " << Apr->energy( u, u ) << "\n";

    }//end of loop over time

    return u;
#if 0
    //*u = M_model->solve( mu );
    do
    {
        //merge all matrices/vectors contributions from affine decomposition
        //result : a tuple : M , A and F
        auto merge = M_model->update( mu );
        A = merge.template get<1>();
        F = merge.template get<2>();
        //backup
        uold = un;

        //solve
        //M_backend->solve( _matrix=A , _solution=un, _rhs=F[0]);
        M_preconditioner_primal->setMatrix( A );
        M_backend_primal->solve( _matrix=A , _solution=un, _rhs=F[0] , _prec=M_preconditioner );

        //on each subspace the norme of the increment is computed and then we perform the sum
        increment_norm = M_model->computeNormL2( un , uold );
        iteration++;

    } while( increment_norm > increment_fixedpoint_tol && iteration < max_fixedpoint_iterations );

    element_ptrtype eltF( new element_type( M_model->functionSpace() ) );
    *eltF = *F[0];
    LOG(INFO) << "[CRB::offlineFixedPoint] u^T F = " << inner_product( *u, *eltF ) ;
    LOG(INFO) << "[CRB::offlineFixedPoint] energy = " << A->energy( *u, *u ) ;
#endif
}//offline fixed point


template<typename TruthModelType>
typename CRB<TruthModelType>::element_type
CRB<TruthModelType>::offlineFixedPointDual(parameter_type const& mu, element_ptrtype & dual_initial_field) //const sparse_matrix_ptrtype & A, const element_type & u )
{
    if ( M_model->isSteady() )
        return M_model->solveFemDualUsingAffineDecompositionFixedPoint( mu );

    //M_backend_dual = backend_type::build( BACKEND_PETSC );
    bool reuse_prec = boption(_prefix=M_prefix,_name="crb.reuse-prec") ;

    auto udu = M_model->functionSpace()->element();

    sparse_matrix_ptrtype M = M_model->newMatrix();
    int nl = M_model->Nl();  //number of outputs
    std::vector< vector_ptrtype > F( nl );
    for(int l=0; l<nl; l++)
        F[l]=M_model->newVector();


    M_bdf_dual = bdf( _space=M_model->functionSpace(), _vm=this->vm() , _name="bdf_dual" );
    M_bdf_dual_save = bdf( _space=M_model->functionSpace(), _vm=this->vm() , _name="bdf_dual_save" );

    M_bdf_dual->setTimeInitial( M_model->timeFinal()+M_model->timeStep() );

    M_bdf_dual->setTimeStep( -M_model->timeStep() );
    M_bdf_dual->setTimeFinal( M_model->timeInitial()+M_model->timeStep() );
    M_bdf_dual->setOrder( M_model->timeOrder() );

    M_bdf_dual_save->setTimeInitial( M_model->timeFinal()+M_model->timeStep() );
    M_bdf_dual_save->setTimeStep( -M_model->timeStep() );
    M_bdf_dual_save->setTimeFinal( M_model->timeInitial()+M_model->timeStep() );
    M_bdf_dual_save->setOrder( M_model->timeOrder() );

    M_bdf_dual_save->setRankProcInNameOfFiles( true );
    M_bdf_dual->setRankProcInNameOfFiles( true );

    auto Adu = M_model->newMatrix();
    auto Apr = M_model->newMatrix();

    double dt = M_model->timeStep();

    int max_fixedpoint_iterations  = ioption(_prefix=M_prefix,_name="crb.fixedpoint.maxit");
    double increment_fixedpoint_tol  = doption(_prefix=M_prefix,_name="crb.fixedpoint.increment-tol");
    double fixedpoint_critical_value  = doption(_prefix=M_prefix,_name="crb.fixedpoint.critical-value");
    int iteration=0;
    double increment_norm=1e3;

    vector_ptrtype Rhs( M_backend_dual->newVector( M_model->functionSpace() ) );
    bool is_linear = M_model->isLinear();
    if( is_linear )
        increment_norm = 0;

    double bdf_coeff;

    auto vec_bdf_poly = M_backend_dual->newVector( M_model->functionSpace() );

    bool POD_WN = boption(_prefix=M_prefix,_name="crb.apply-POD-to-WN") ;

    if ( M_model->isSteady() )
        udu.zero() ;
    else
    {
        boost::tie( M, Apr, F) = M_model->update( mu , M_bdf_dual->timeInitial() );

#if 0
        Apr->addMatrix( 1./dt, M );
        Apr->transpose( Adu );
        *Rhs=*F[M_output_index];
        Rhs->scale( 1./dt );
        // M->scale(1./dt);

        M_preconditioner_dual->setMatrix( Adu );
        M_backend_dual->solve( _matrix=Adu, _solution=dual_initial_field, _rhs=Rhs, _prec=M_preconditioner_dual );
#else
        *Rhs=*F[M_output_index];
        //Rhs->scale( 1./dt );
        //M->scale(1./dt);
        M_preconditioner_dual->setMatrix( M );
        M_backend_dual->solve( _matrix=M, _solution=dual_initial_field, _rhs=Rhs, _prec=M_preconditioner_dual );
#endif
        udu=*dual_initial_field;
    }

    auto uold = M_model->functionSpace()->element();
    auto bdf_poly = M_model->functionSpace()->element();


    for ( M_bdf_dual->start(udu),M_bdf_dual_save->start(udu);
          !M_bdf_dual->isFinished() && !M_bdf_dual_save->isFinished();
          M_bdf_dual->next(udu) , M_bdf_dual_save->next() )
    {

        int bdf_iter = M_bdf_dual->iteration();

        if ( ! M_model->isSteady() )
        {
            bdf_coeff = M_bdf_dual->polyDerivCoefficient( 0 );
            bdf_poly = M_bdf_dual->polyDeriv();
        }

        do
        {
            if( is_linear )
            {
                bool compute_only_terms_time_dependent=false;
                if ( bdf_iter == 1 )
                {
                    boost::tie( M, Apr, F) = M_model->update( mu , M_bdf_dual->time() , compute_only_terms_time_dependent );
                }
                else
                {
                    compute_only_terms_time_dependent=true;
                    boost::tie( boost::tuples::ignore, boost::tuples::ignore, F) = M_model->update( mu , M_bdf_dual->time() , compute_only_terms_time_dependent );
                }
            }
            else
            {
                bool compute_only_terms_time_dependent=false;
                if ( bdf_iter == 1 )
                {
                    boost::tie( M, Apr, F) = M_model->update( mu , udu, M_bdf_dual->time() , compute_only_terms_time_dependent );
                }
                else
                {
                    compute_only_terms_time_dependent=true;
                    boost::tie( boost::tuples::ignore, boost::tuples::ignore, F) = M_model->update( mu , udu, M_bdf_dual->time() , compute_only_terms_time_dependent );
                }
            }

            if( ! M_model->isSteady() )
            {
                if ( bdf_iter == 1 )
                {
                    Apr->addMatrix( bdf_coeff, M );
                }
                Rhs->zero();
                *vec_bdf_poly = bdf_poly;
                Rhs->addVector( *vec_bdf_poly, *M );
            }
            else
            {
                *Rhs = *F[M_output_index];
                Rhs->close();
                Rhs->scale( -1 );
            }


            if( boption(_prefix=M_prefix,_name="crb.use-symmetric-matrix") )
                Adu = Apr;
            else
                Apr->transpose( Adu );
            //Adu->close();
            //Rhs->close();

            //backup for non linear problems
            uold = udu;

            //solve
            M_preconditioner_dual->setMatrix( Adu );
            auto ret = M_backend_dual->solve( _matrix=Adu, _solution=udu, _rhs=Rhs,  _prec=M_preconditioner_dual, _reuse_prec=( bdf_iter >=2 ) );
            if  ( !ret.template get<0>() )
                LOG(INFO)<<"[CRB] WARNING : at time "<<M_bdf_dual->time()<<" we have not converged ( nb_it : "<<ret.template get<1>()<<" and residual : "<<ret.template get<2>() <<" ) \n";

            //on each subspace the norme of the increment is computed and then we perform the sum
            if( is_linear )
                increment_norm = 0;
            else
                increment_norm = M_model->computeNormL2( udu , uold );

            iteration++;

        }while( increment_norm > increment_fixedpoint_tol && iteration < max_fixedpoint_iterations );

        //check dual property
        //double term1 = A->energy( udu, u );
        //double term2 = Adu->energy( u, udu );
        //double diff = math::abs( term1-term2 );
        //LOG(INFO) << "< A u , udu > - < u , A* udu > = "<<diff<<"\n";

        if ( ! M_model->isSteady() )
        {
            udu.close();
            if( POD_WN )
            {
                M_bdf_dual_save->shiftRight( udu );
            }
            else
            {
                element_ptrtype projection ( new element_type ( M_model->functionSpace() ) );
                projectionOnPodSpace ( udu , projection, "dual" );
                auto e=udu;
                e-=*projection;
                M_bdf_dual_save->shiftRight( e );
            }
        }

        if( increment_norm > fixedpoint_critical_value )
            throw std::logic_error( fmt::format("[CRB::offlineFixedPointDual]  at time {} ERROR : increment > critical value ", M_bdf_dual->time() ) );

    }//end of loop over time

    return udu;
}//offline fixed point

template<typename TruthModelType>
typename CRB<TruthModelType>::element_type
CRB<TruthModelType>::femSolve( parameter_type& mu )
{
    if ( M_model->hasEim() )
        return M_model->solve(mu);

    auto u = M_model->functionSpace()->element();
    auto udu = M_model->functionSpace()->element();
    element_ptrtype dual_initial_field( new element_type( M_model->functionSpace() ) );
    this->offlineSolve( u, udu, mu, dual_initial_field );

    return u;
}


template<typename TruthModelType>
void
CRB<TruthModelType>::offlineSolve( element_type& u, element_type& udu, parameter_type& mu, element_ptrtype & dual_initial_field )
{
    if ( M_model->isSteady()  )
    {
        tic();
        if ( boption(_prefix=M_prefix,_name="crb.solve-fem-monolithic") )
            u = M_model->solve(mu);
        else if( M_use_newton )
            u = offlineNewtonPrimal( mu );
        else
            u = offlineFixedPointPrimal( mu );
        toc("Solve primal problem");

        if( M_solve_dual_problem )
        {
            if ( M_use_newton )
                return;
            tic();
            udu = offlineFixedPointDual( mu , dual_initial_field );
            toc("Solve dual problem");
        }

#if 0
        // TODO VINCENT : fix use of dual problem
        if ( M_model->hasEim() && boption(_prefix=M_prefix,_name="crb.solve-fem-monolithic") )
        {
            u = M_model->solve(mu);
        }
        else if( ! M_use_newton )
        {
            u = offlineFixedPointPrimal( mu );//, A  );
            toc("Solve primal problem");

            if( M_solve_dual_problem )
            {
                tic();
                udu = offlineFixedPointDual( mu , dual_initial_field );//,  A , u );
                toc("Solve dual problem");
            }
        }
        else
        {
            tic();
            u = offlineNewtonPrimal( mu );
            toc("Solve primal problem");
        }
#endif
    }//steady
    else
    {
        tic();
        u = offlineFixedPointPrimal( mu  );
        toc("Solve primal problem");

        if ( M_solve_dual_problem || M_error_type==CRB_RESIDUAL || M_error_type == CRB_RESIDUAL_SCM )
        {
            tic();
            udu = offlineFixedPointDual( mu , dual_initial_field );
            toc("Solve dual problem");
        }
    }//transient
}


template<typename TruthModelType>
typename CRB<TruthModelType>::element_type
CRB<TruthModelType>::offlineNewtonPrimal( parameter_type const& mu )
{
    return M_model->solveFemUsingAffineDecompositionNewton( mu );
#if 0
    sparse_matrix_ptrtype J = M_model->newMatrix();
    vector_ptrtype R = M_model->newVector();

    auto initialguess = M_model->functionSpace()->elementPtr();
    initialguess = M_model->assembleInitialGuess( mu ) ;

    beta_vector_type betaJqm;
    std::vector< beta_vector_type > betaRqm;

    boost::tie( boost::tuples::ignore, M_Jqm, M_Rqm ) = M_model->computeAffineDecomposition();

    M_backend_primal->nlSolver()->jacobian = std::bind( &self_type::offlineUpdateJacobian,
                                                          std::ref( *this ), std::placeholders::_1, std::placeholders::_2, mu );
    M_backend_primal->nlSolver()->residual = std::bind( &self_type::offlineUpdateResidual,
                                                          std::ref( *this ), std::placeholders::_1, std::placeholders::_2, mu );
#if 0
    M_backend_primal->nlSolver()->setType( TRUST_REGION );
#endif
    auto solution = M_model->functionSpace()->element();
    M_backend_primal->nlSolve(_jacobian=J, _solution=solution, _residual=R);

    M_model->solve( mu );// supp

    return solution;
#endif
}

template<typename TruthModelType>
void
CRB<TruthModelType>::offlineUpdateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype & J , const parameter_type & mu)
{
    J->zero();

    beta_vector_type betaJqm;
    boost::tie( boost::tuples::ignore, betaJqm, boost::tuples::ignore ) = M_model->computeBetaQm( X , mu , 0 );
    M_model->updateJacobian(X,M_Jqm);

    for ( size_type q = 0; q < M_model->Qa(); ++q )
    {
        for(int m=0; m<M_model->mMaxA(q); m++)
        {
            J->addMatrix( betaJqm[q][m], M_Jqm[q][m] );
        }
    }
}

template<typename TruthModelType>
void
CRB<TruthModelType>::offlineUpdateResidual( const vector_ptrtype& X, vector_ptrtype& R , const parameter_type & mu)
{
    R->zero();
    std::vector< beta_vector_type > betaRqm;
    boost::tie( boost::tuples::ignore, boost::tuples::ignore, betaRqm ) = M_model->computeBetaQm( X , mu , 0 );
    M_model->updateResidual(X,M_Rqm);
    for ( size_type q = 0; q < M_model->Ql( 0 ); ++q )
    {
        for(int m=0; m<M_model->mMaxF(0,q); m++)
            R->add( betaRqm[0][q][m] , *M_Rqm[0][q][m] );
    }
}

template<typename TruthModelType>
typename CRB<TruthModelType>::vector_ptrtype
CRB<TruthModelType>::computeRieszResidual( parameter_type const& mu, std::vector<vectorN_type> const& uN ) const
{
    auto all_beta = M_model->computePicardBetaQm( this->expansion( uN[0] ), mu );
    auto beta_A = all_beta.template get<1>();
    auto beta_F = all_beta.template get<2>();

    vector_ptrtype YR( backend()->newVector(M_model->functionSpace()) );
    YR->zero();

    int Qa = beta_A.size();
    for ( size_type q = 0; q < Qa; ++q )
    {
        int mMaxA = M_hAqm[q].size();
        for(size_type m = 0; m < mMaxA; ++m )
        {
            for(size_type n=0; n<M_N; ++n)
            {
                YR->add( beta_A[q][m]*uN[0](n), M_hAqm[q][m][n] );
            }
        }
    }
    YR->scale( -1.0 );

    int Ql = beta_F[0].size();
    for ( size_type q = 0; q < Ql; ++q )
    {
        int mMaxF = M_hFqm[q].size();
        for ( size_type m = 0; m < mMaxF; ++m )
        {
            YR->add( beta_F[0][q][m] , M_hFqm[q][m] );
        }
    }

    YR->close();
    return YR;
}

template<typename TruthModelType>
double
CRB<TruthModelType>::computeRieszResidualNorm( parameter_type const& mu ) const
{
    std::vector<vectorN_type> uN;
    std::vector<vectorN_type> uNdu;
    std::vector<vectorN_type> uNold;
    std::vector<vectorN_type> uNduold;
    auto o = lb( this->dimension(), mu, uN, uNdu, uNold, uNduold  );
    return this->computeRieszResidualNorm( mu, uN );
}

template<typename TruthModelType>
double
CRB<TruthModelType>::computeRieszResidualNorm( parameter_type const& mu, std::vector<vectorN_type> const& uN ) const
{
    vector_ptrtype YR( this->computeRieszResidual( mu, uN ) );
    double dual_norm = math::sqrt( M_model->scalarProduct(YR,YR) );
    //Feel::cout << "[computeRieszResidualNorm] norm = " << dual_norm << std::endl;
    return dual_norm;
}

template<typename TruthModelType>
void
CRB<TruthModelType>::loadSCMDB()
{
    if ( M_error_type == CRB_RESIDUAL_SCM )
    {
        M_scmA->setScmForMassMatrix( false );
        if( ! M_scmA->loadDB() )
            std::vector<boost::tuple<double,double,double> > M_rbconv2 = M_scmA->offline();

        if ( ! M_model->isSteady() )
        {
            M_scmM->setScmForMassMatrix( true );
            if( ! M_scmM->loadDB() )
                std::vector<boost::tuple<double,double,double> > M_rbconv3 = M_scmM->offline();
        }
    }
}

template<typename TruthModelType>
typename CRB<TruthModelType>::convergence_type
CRB<TruthModelType>::offline()
{
    this->setOfflineStep( true );

    if ( this->worldComm().isMasterRank() )
    {
        if ( !fs::exists( this->dbLocalPath() ) )
            fs::create_directories( this->dbLocalPath() );
    }
    this->worldComm().barrier();

    int proc_number = this->worldComm().globalRank();
    int master_proc = this->worldComm().masterRank();
    //M_rbconv_contains_primal_and_dual_contributions = true;

    M_Nm = ioption(_prefix=M_prefix,_name="crb.Nm") ;
    number_of_added_elements = 1;
    bool seek_mu_in_complement = M_seekMuInComplement;

    Feel::Timer ti;
    LOG(INFO)<< "[CRB::offline] Offline CRB starts, this may take a while until Database is computed...\n";
    //LOG(INFO) << "[CRB::offline] Starting offline for output " << M_output_index << "\n";
    //LOG(INFO) << "[CRB::offline] initialize underlying finite element model\n";
    M_model->init();
    //LOG(INFO) << " -- model init done in " << ti.elapsed() << "s\n";

    parameter_type mu( M_Dmu );

    double delta_pr;
    double delta_du;
    size_type index = 0;
    bool use_predefined_WNmu = false;
    int Nrestart = ioption(_prefix=M_prefix,_name="crb.restart-from-N");
    int Frestart = ioption(_prefix=M_prefix,_name="ser.rb-rebuild-freq");


#if 0
    //we do affine decomposition here to then have access to Qa, mMax ect...
    //and then resize data structure if necessary
    std::vector< std::vector<sparse_matrix_ptrtype> > Jqm;
    std::vector< std::vector<std::vector<vector_ptrtype> > > Rqm;
    LOG(INFO) << "[CRB::offline] compute affine decomposition\n";
    std::vector< std::vector<sparse_matrix_ptrtype> > Aqm;
    std::vector< std::vector<sparse_matrix_ptrtype> > Mqm;
    //to project the initial guess on the reduced basis we solve A u = F
    //with A = \int_\Omega u v ( mass matrix )
    //F = \int_\Omega initial_guess v
    //so InitialGuessV is the mu-independant part of the vector F
    std::vector< std::vector<element_ptrtype> > InitialGuessV;
    std::vector< std::vector<std::vector<vector_ptrtype> > > Fqm,Lqm;
#endif

    bool model_is_linear = M_model->isLinear();
    if( ! model_is_linear )
        M_InitialGuessV = M_model->computeInitialGuessVAffineDecomposition();

    if( M_use_newton )
    {
        //boost::tie( Mqm , Jqm, Rqm ) = M_model->computeAffineDecomposition();
        boost::tie( M_Mqm , M_Jqm, M_Rqm ) = M_model->computeAffineDecomposition();
    }
    else
    {
        if( boption(_prefix=M_prefix,_name="crb.stock-matrices") )
            //boost::tie( Mqm, Aqm, Fqm ) = M_model->computeAffineDecomposition();
            boost::tie( M_Mqm, M_Aqm, M_Fqm ) = M_model->computeAffineDecomposition();
    }

    M_model->countAffineDecompositionTerms();
    this->updateAffineDecompositionSize();

    int Nold = M_N; // store current M_N in case of restart
    bool freq_restart = ( Frestart > 0 && (M_N-1) % Frestart == 0 );
    if( Nrestart == 0 || freq_restart )
        this->setRebuild( true );
    if( Nrestart > 0 )
        this->setRebuild( false );

    //if M_N == 0 then there is not an already existing database
    if( M_rebuild || M_N == 0 )
    {
        ti.start();
        //scm_ptrtype M_scm = scm_ptrtype( new scm_type( M_vm ) );
        //M_scm->setTruthModel( M_model );
        //    std::vector<boost::tuple<double,double,double> > M_rbconv2 = M_scm->offline();

        M_coeff_pr_ini_online.resize(0);
        M_coeff_du_ini_online.resize(0);

        this->generateSuperSampling();

        if( this->worldComm().isMasterRank() )
            std::cout<<"[CRB offline] M_error_type = "<<M_error_type<<std::endl;

        if( this->worldComm().isMasterRank() )
            std::cout << " -- sampling init done in " << ti.elapsed() << "s\n";
        ti.start();

        if( this->worldComm().isMasterRank() )
            std::cout << " -- residual data init done in " << ti.elapsed() << std::endl;
        ti.start();

        // empty sets
        M_WNmu->clear();
        M_primal_apee_mu->clear();
        M_dual_apee_mu->clear();

        use_predefined_WNmu = buildSampling();
        if( M_error_type == CRB_NO_RESIDUAL )
        {
            if( this->worldComm().isMasterRank() )
            {
                //we don't want broadcast the element
                bool broadcast=false;
                mu = M_Dmu->element( broadcast );
            }
            //every proc must have the same parameter mu
            boost::mpi::broadcast( this->worldComm() , mu , master_proc );
        }
        else
        {
            // start with M_C = { arg min mu, mu \in Xi }
            //the min is a global min so we can check
            bool check=true;
            boost::tie( mu, index ) = M_Xi->min( check );
        }
        M_current_mu = mu;

        int size = mu.size();
        if( this->worldComm().isMasterRank() )
        {
            std::cout << "  -- start with mu = [ ";
            for ( int i=0; i<size-1; i++ ) std::cout<<mu( i )<<" ";
            std::cout<<mu( size-1 )<<" ]"<<std::endl;
        }
        //std::cout << " -- WN size :  " << M_WNmu->size() << "\n";

        // dimension of reduced basis space
        M_N = 0;

        // scm offline stage: build C_K
        this->loadSCMDB();

        M_maxerror = 1e10;
        delta_pr = 0;
        delta_du = 0;
        //boost::tie( M_maxerror, mu, index ) = maxErrorBounds( N );

        // save mesh of rbspace
        if ( true )
        {
            std::string meshFilenameBase = fmt::format("{}_mesh_p{}.json",this->name(),this->worldComm().size());
            std::string meshFilename = (M_elements_database.dbLocalPath() / fs::path(meshFilenameBase)).string();
            if( this->worldComm().isMasterRank() )
                std::cout << "save Mesh : " << meshFilename << std::endl;
            M_model->rBFunctionSpace()->saveMesh( meshFilename );
        }
        // save the rb function space
        if ( true )
        {
            std::string functionSpaceFilenameBase = fmt::format("{}_functionSpace_p{}.json",this->name(),this->worldComm().size());
            std::string functionSpaceFilename = (M_elements_database.dbLocalPath() / fs::path(functionSpaceFilenameBase)).string();
            if( this->worldComm().isMasterRank() )
                std::cout << "save Function Space : " << functionSpaceFilename << std::endl;
            M_model->rBFunctionSpace()->functionSpace()->save( functionSpaceFilename );
        }
        M_model->copyAdditionalModelFiles( this->dbDirectory() );



    }//end of if( rebuild_database )
    else
    {
        // scm offline stage: build C_K
        this->loadSCMDB();

        bool restart = false;
        int wnmu_nb_elements=M_WNmu->nbElements();
        restart = ( Nrestart > 1 && Nrestart < wnmu_nb_elements );

        if( restart )
        {
            int number_of_elem_to_remove = wnmu_nb_elements - Nrestart ;

            for(int i=0; i<number_of_elem_to_remove; i++)
            {
                if( i==number_of_elem_to_remove-1 )
                    mu = M_WNmu->lastElement();
                M_WNmu->pop_back();
                M_rbconv.left.erase( M_N-i );
            }

            M_model->rBFunctionSpace()->deleteLastPrimalBasisElements( number_of_elem_to_remove );
            if( M_solve_dual_problem )
                M_model->rBFunctionSpace()->deleteLastDualBasisElements( number_of_elem_to_remove );

            M_N = Nrestart;
            if( this->worldComm().isMasterRank() )
                std::cout<<"Restart the RB construction at N = "<<Nrestart<<std::endl;
            LOG( INFO ) << "Restart the RB construction at N = "<<Nrestart;
        }
        else
        {
            mu = M_current_mu;
            if( this->worldComm().isMasterRank() && M_N < ioption(_prefix=M_prefix,_name="crb.dimension-max") - 1 )
            {
                std::cout<<"We are going to enrich the reduced basis"<<std::endl;
                std::cout<<"There are "<<M_N<<" elements in the database"<<std::endl;
            }
            LOG(INFO) <<"we are going to enrich the reduced basis"<<std::endl;
            LOG(INFO) <<"there are "<<M_N<<" elements in the database"<<std::endl;
        }
    }//end of else associated to if ( rebuild_database )

    // Need to update affine decomposition size if it has changed
    this->updateAffineDecompositionSize();

#if 0
    sparse_matrix_ptrtype Aq_transpose = M_model->newMatrix();
    sparse_matrix_ptrtype A = M_model->newMatrix();
    int nl = M_model->Nl();
    std::vector< vector_ptrtype > F( nl );
    for(int l=0; l<nl; l++)
        F[l]=M_model->newVector();
#endif
    element_ptrtype dual_initial_field( new element_type( M_model->functionSpace() ) );
    auto u = M_model->functionSpace()->element();
    auto udu = M_model->functionSpace()->element();

    LOG(INFO) << "[CRB::offline] starting offline adaptive loop\n";

    //bool reuse_prec = boption(_prefix=M_prefix,_name="crb.reuse-prec") ;



    LOG(INFO) << "[CRB::offline] strategy "<< M_error_type <<"\n";
    if( this->worldComm().isMasterRank() )
        std::cout << "[CRB::offline] strategy "<< M_error_type <<std::endl;

    use_predefined_WNmu = this->M_error_type==CRB_NO_RESIDUAL || boption(_prefix=M_prefix,_name="crb.use-predefined-WNmu");
    if( M_error_type == CRB_NO_RESIDUAL || use_predefined_WNmu )
    {
        //in this case it makes no sens to check the estimated error
        M_maxerror = 1e10;
    }

    bool all_procs = boption(_prefix=M_prefix,_name="crb.system-memory-evolution-on-all-procs") ;
    PsLogger ps ("PsLogCrbOffline" , this->worldComm() , "rss pmem pcpu" , all_procs );

    bool only_master=boption(_prefix=M_prefix,_name="crb.system-memory-evolution");
    bool only_one_proc= only_master * ( this->worldComm().isMasterRank() );
    bool write_memory_evolution = all_procs || only_one_proc ;

    int user_max = ioption(_prefix=M_prefix,_name="crb.dimension-max");
    int rb_frequency = ioption(_prefix=M_prefix,_name="ser.rb-frequency");
    if( rb_frequency != 0 ) // SER
    {
        if( M_N == 0 || (this->adaptationSER() && Nold < user_max) )
            M_iter_max = Nold + 1; //Initialization step || add a suppl. basis is needed
        else if( Nold + rb_frequency <= user_max )
            M_iter_max = Nold + rb_frequency;
        else
            M_iter_max = user_max;
        // if( proc_number == this->worldComm().masterRank() )
        //     std::cout << "[crb - SER] M_N = " << M_N << ", N_old = " << Nold << ", M_iter_max = " << M_iter_max << std::endl;
    }

    if ( use_predefined_WNmu )
    {
        if ( M_iter_max>this->M_WNmu->size() )
            M_iter_max = this->M_WNmu->size();
        mu = this->M_WNmu->at( std::min( M_N, this->M_iter_max-1 ) );
        M_current_mu =mu;
    }

    while ( M_maxerror > M_tolerance && M_N < M_iter_max  )
    {
        tic();
        M_mode_number=1;

        std::string pslogname = fmt::format("N-{}",M_N );
        if( write_memory_evolution )
            ps.log(pslogname);

        if( this->worldComm().isMasterRank() )
        {
            std::cout << "N = " << M_N+1 << "/"  << M_iter_max << "( max = " << user_max << ")\n";
            int size = mu.size();
            std::cout << "  -- mu = [ ";
            for ( int i=0; i< size-1; i++ )
                std::cout<< M_current_mu( i ) <<" ";
            std::cout << M_current_mu( size-1 ) << " ]" <<std::endl;
        }

        Feel::Timer timer2, timer3;
        LOG(INFO) <<"========================================"<<"\n";

        if ( M_error_type == CRB_NO_RESIDUAL )
            LOG(INFO) << "N=" << M_N << "/"  << M_iter_max << " ( nb proc : "<<worldComm().globalSize()<<")";
        else
            LOG(INFO) << "N=" << M_N << "/"  << M_iter_max << " maxerror=" << M_maxerror << " / "  << M_tolerance << "( nb proc : "<<worldComm().globalSize()<<")";

        // for a given parameter \p mu assemble the left and right hand side
        u.setName( fmt::format( "fem-primal-N{}-proc{}", M_N, proc_number ) );
        udu.setName( fmt::format( "fem-dual-N{}-proc{}", M_N, proc_number) );

#if !defined(NDEBUG)
        mu.check();
#endif

        double tpr=0,tdu=0;

        this->offlineSolve( u, udu, mu, dual_initial_field );

        if( ! use_predefined_WNmu )
            M_WNmu->push_back( mu, index );
        M_WNmu_complement = M_WNmu->complement();

        bool norm_zero = false;
        timer2.start();
        timer3.start();
        if ( M_model->isSteady() )
        {
            this->addBasis( u, udu, mu );
        }//end of steady case

        else
        {
            if ( M_N == 0 && M_orthonormalize_primal==false )
            {
                //add initial solution as element in the reduced basis if it's not zero
                //else there will be problems during orthonormalization step
                //note : at this step *u contains solution of primal probleme at final time
                //so it's the initial solution of the dual problem
                //in the case where orthonormalize_primal==true, even with initial solution in
                //the reduced basis we will commit an approximation in the online part
                //when initializing uN.
                //if initial solution is zero then the system the matrix used in online will be singular
                element_ptrtype primal_initial_field ( new element_type ( M_model->functionSpace() ) );

                M_model->initializationField( primal_initial_field, mu ); //fill initial_field
                int norm_ini = M_model->scalarProduct( *primal_initial_field, *primal_initial_field );

                if ( norm_ini != 0 )
                {
                    M_model->rBFunctionSpace()->addPrimalBasisElement( *primal_initial_field );
                    M_model->rBFunctionSpace()->addDualBasisElement( *dual_initial_field );
                    //M_WN.push_back( *primal_initial_field );
                    //M_WNdu.push_back( *dual_initial_field );
                }

                else
                {
                    norm_zero=true;
                }
            } // end of if ( M_N == 0 )

            //POD in time
            LOG(INFO)<<"[CRB::offline] start of POD \n";

            pod_ptrtype POD = pod_ptrtype( new pod_type(  ) );
            bool POD_WN = boption(_prefix=M_prefix,_name="crb.apply-POD-to-WN") ;

            if ( seek_mu_in_complement ) // M_mode_number == 1 )
            {
                if( POD_WN )
                {
                    POD->setNm( -1 );
                }
                else
                {
                    //in this case, it's the first time that we add mu
                    POD->setNm( M_Nm );
                }
            }
            else
            {
                if( POD_WN )
                {
                    POD->setNm( -1 );
                }
                else
                {
                    //in this case we have to count mu occurrences in WMmu (mode_number)
                    //to add the mode_number^th mode in the basis
                    M_mode_number=1;
                    for( auto _mu : *M_WNmu )
                    {
                        if( mu == _mu )
                            M_mode_number++;
                    }

                    POD->setNm( M_mode_number*M_Nm );
                    int size = mu.size();
                    LOG(INFO)<<"... CRB M_mode_number = "<<M_mode_number<<"\n";
                    LOG(INFO)<<"for mu = [ ";

                    for ( int i=0; i<size-1; i++ ) LOG(INFO)<<mu[i]<<" , ";

                    LOG(INFO)<<mu[ size-1 ];
                    LOG(INFO)<<" ]\n";

                    double Tf = M_model->timeFinal();
                    double dt = M_model->timeStep();
                    int nb_mode_max = Tf/dt;

                    if ( M_mode_number>=nb_mode_max-1 )
                    {
                        if( this->worldComm().isMasterRank() )
                        {
                            std::cout<<"Error : we access to "<<M_mode_number<<"^th mode"<<std::endl;
                            std::cout<<"parameter choosen : [ ";
                            for ( int i=0; i<size-1; i++ ) std::cout<<mu[i]<<" , ";
                            std::cout<<mu[ size-1 ]<<" ] "<<std::endl;
                        }
                        throw std::logic_error( "[CRB::offline] ERROR during the construction of the reduced basis, one parameter has been choosen too many times" );
                    }
                }
            }

            POD->setBdf( M_bdf_primal_save );
            POD->setModel( M_model );
            POD->setTimeInitial( M_model->timeInitial() );
            mode_set_type ModeSet;

            timer2.start();
            size_type number_max_of_mode = POD->pod( ModeSet,true );
            tpr=timer2.elapsed();

            if ( number_max_of_mode < M_Nm )
            {
                LOG( INFO )<<"With crb.Nm = "<<M_Nm<<" there was too much too small eigenvalues so the value"
                           <<" of crb.Nm as been changed to "<<number_max_of_mode;
                M_Nm = number_max_of_mode;
            }

            //now : loop over number modes per mu

            if ( !seek_mu_in_complement )
            {

                if( POD_WN )
                {
                    int imax = ModeSet.size();
                    for ( size_type i=0; i<imax; i++ )
                    {
                        M_model->rBFunctionSpace()->addPrimalBasisElement( ModeSet[i] );
                    }
                }
                else
                {
                    for ( size_type i=0; i<M_Nm; i++ )
                    {
                        M_model->rBFunctionSpace()->addPrimalBasisElement( ModeSet[M_mode_number*M_Nm-1+i] );
                        //M_WN.push_back( ModeSet[M_mode_number*M_Nm-1+i] ) ;
                        //M_WN.push_back( ModeSet[M_mode_number-1] ) ;
                    }
                }
            }//! seek in complement
            else
            {
                int imax = ModeSet.size();
                for ( size_type i=0; i<imax; i++ )
                {
                    M_model->rBFunctionSpace()->addPrimalBasisElement( ModeSet[i] );
                    //M_WN.push_back( ModeSet[i] ) ;
                }
            }

            //and now the dual
            if ( M_solve_dual_problem || M_error_type==CRB_RESIDUAL || M_error_type == CRB_RESIDUAL_SCM )
            {
                POD->setBdf( M_bdf_dual );
                POD->setTimeInitial( M_model->timeFinal()+M_model->timeStep() );
                mode_set_type ModeSetdu;

                timer2.start();
                POD->pod( ModeSetdu,false );
                tdu=timer2.elapsed();

                if ( !seek_mu_in_complement )
                {
                    if( POD_WN )
                    {
                        int imax = ModeSetdu.size();
                        for ( size_type i=0; i<imax; i++ )
                        {
                            M_model->rBFunctionSpace()->addDualBasisElement( ModeSet[i] );
                        }
                    }
                    else
                    {
                        for ( size_type i=0; i<M_Nm; i++ )
                        {
                            M_model->rBFunctionSpace()->addDualBasisElement( ModeSetdu[M_mode_number*M_Nm-1+i] );
                            //M_WNdu.push_back( ModeSetdu[M_mode_number*M_Nm-1+i] ) ;
                        }
                    }
                }//! seek mu in complement
                else
                {
                    int imax = ModeSetdu.size();
                    for ( size_type i=0; i<imax; i++ )
                    {
                        M_model->rBFunctionSpace()->addDualBasisElement( ModeSetdu[i] );
                        //M_WNdu.push_back( ModeSetdu[i] ) ;
                    }
                }
            }

            else
            {
                element_ptrtype element_zero ( new element_type ( M_model->functionSpace() ) );

                for ( size_type i=0; i<M_Nm; i++ )
                {
                    M_model->rBFunctionSpace()->addDualBasisElement( *element_zero );
                    //M_WNdu.push_back( *element_zero ) ;
                }
            }


        }//end of transient case



        //in the case of transient problem, we can add severals modes for a same mu
        //Moreover, if the case where the initial condition is not zero and we don't orthonormalize elements in the basis,
        //we add the initial condition in the basis (so one more element)
        number_of_added_elements = M_Nm + ( M_N==0 && M_orthonormalize_primal==false && norm_zero==false && !M_model->isSteady() );

        //in the case of steady problems, we add only one element
        if( M_model->isSteady() )
            number_of_added_elements=1;

        M_N+=number_of_added_elements;

        this->orthonormalizeBasis( number_of_added_elements );

        this->buildRbMatrix( number_of_added_elements, mu, dual_initial_field );

        if( M_useAccurateApee )
        {
            if( M_solve_dual_problem )
            {
                std::vector< parameter_type > new_primal_parameters;
                selectPrimalApeeParameters( M_N,  new_primal_parameters );
                std::vector< parameter_type > new_dual_parameters;
                selectDualApeeParameters( M_N,  new_dual_parameters );
                updateApeeBasisConstruction( M_N , new_primal_parameters , new_dual_parameters );
                updateApeeOfflineResidualComputation( M_N , new_primal_parameters , new_dual_parameters );
                assemblePrimalAndDualApeeBasisMatrices( M_N );
            }
            else
                throw std::logic_error( "[CRB::offline] if you want accurate a posteriori error estimation make sure that crb.solve-dual-problem=true " );
        }



        if ( M_database_contains_variance_info )
            throw std::logic_error( "[CRB::offline] ERROR : build variance is not actived" );
        //buildVarianceMatrixPhi( M_N );

        if ( M_useRbSpaceContextEim )
        {
            M_model->updateRbSpaceContextEim();
            M_hasRbSpaceContextEim = true;
        }
        this->updateRbInDeim();

        if ( M_error_type==CRB_RESIDUAL || M_error_type == CRB_RESIDUAL_SCM )
        {
            if( this->worldComm().isMasterRank() )
                std::cout << "  -- offlineResidual update starts\n";
            timer2.start();
            offlineResidual( M_N, number_of_added_elements );
            LOG(INFO)<<"[CRB::offline] end of call offlineResidual and M_N = "<< M_N <<"\n";
            if( this->worldComm().isMasterRank() )
                std::cout << "  -- offlineResidual updated in " << timer2.elapsed() << "s\n";
            bool model_has_eim_error = M_model->hasEimError();
            if( model_has_eim_error )
                offlineResidualEim( M_N, number_of_added_elements );
            timer2.start();
        }

        bool ser_greedy = M_SER_useGreedyInRb;
        if ( M_error_type == CRB_NO_RESIDUAL && ! use_predefined_WNmu && !ser_greedy )
        {
            //M_maxerror=M_iter_max-M_N;

            if( this->worldComm().isMasterRank() )
            {
                bool already_exist;
                do
                {
                    //initialization
                    already_exist=false;
                    //pick randomly an element
                    bool broadcast=false;
                    mu = M_Dmu->element( broadcast );
                    //make sure that the new mu is not already is M_WNmu
                    for( auto _mu : *M_WNmu )
                    {
                        if( mu == _mu )
                            already_exist=true;
                    }
                }
                while( already_exist );

            }

            boost::mpi::broadcast( this->worldComm() , mu , master_proc );
            M_current_mu = mu;

        }
        else if ( use_predefined_WNmu )
        {
            //remmber that in this case M_iter_max = sampling size
            if( M_N < M_WNmu->size() )
            {
                mu = M_WNmu->at( M_N );
                M_current_mu = mu;
            }
        }
        else
        {
            timer2.start();
            boost::tie( M_maxerror, mu , delta_pr , delta_du ) = maxErrorBounds( M_N );
            auto time=timer2.elapsed();
            M_current_mu = mu;

            if( this->worldComm().isMasterRank() )
            {
                int size = mu.size();
                std::cout << "  -- max error bound (" << M_maxerror << ") computed in " << time << "s"
                          << "  -- mu = [ ";
                for ( int i=0; i< size-1; i++ )
                    std::cout<< M_current_mu( i ) <<" ";
                std::cout << M_current_mu( size-1 ) << " ]" <<std::endl;
            }
        }


        M_rbconv.insert( convergence( M_N, boost::make_tuple(M_maxerror,delta_pr,delta_du) ) );
        //mu = M_Xi->at( M_N );//M_WNmu_complement->min().template get<0>();

        if ( ioption(_prefix=M_prefix,_name="crb.check.rb") == 1 )
        {
            timer2.start();
            check( M_WNmu->size() );
            std::cout << "  -- check reduced basis done in " << timer2.elapsed() << "s\n";
        }

        //save DB after adding an element
        tic();
        this->saveDB();
        this->saveRB();
        toc("Saving the Database");

        toc("Total Time");
        Feel::cout << "============================================================\n";
        LOG(INFO) <<"========================================"<<"\n";
    }

    if( this->worldComm().isMasterRank() )
        std::cout<<"number of elements in the reduced basis : "<<M_N<<" ( nb proc : "<<worldComm().globalSize()<<")"<<std::endl;


    if ( boption(_prefix=M_prefix,_name="crb.check.residual") )
        this->testResidual();

    if( M_maxerror <= M_tolerance || M_N >= user_max  )
    {
        this->setOfflineStep( false );
        if (boption(_prefix=M_prefix,_name="crb.visualize-basis"))
            this->exportBasisFunctions();
    }
    return M_rbconv;
}

template<typename TruthModelType>
void
CRB<TruthModelType>::updateAffineDecompositionSize()
{
    if ( M_error_type == CRB_RESIDUAL || M_error_type == CRB_RESIDUAL_SCM )
    {

        int __QLhs = M_model->Qa();
        int __QRhs = M_model->Ql( 0 );
        int __QOutput = M_model->Ql( M_output_index );
        int __Qm = M_model->Qm();

        //typename array_2_type::extent_gen extents2;
        //M_C0_pr.resize( extents2[__QRhs][__QRhs] );
        //M_C0_du.resize( extents2[__QOutput][__QOutput] );
        M_C0_pr.resize( __QRhs );
        for( int __q1=0; __q1< __QRhs; __q1++)
        {
            int __mMaxQ1=M_model->mMaxF(0,__q1);
            M_C0_pr[__q1].resize( __mMaxQ1 );
            for( int __m1=0; __m1< __mMaxQ1; __m1++)
            {
                M_C0_pr[__q1][__m1].resize(  __QRhs );
                for( int __q2=0; __q2< __QRhs; __q2++)
                {
                    int __mMaxQ2=M_model->mMaxF(0,__q2);
                    M_C0_pr[__q1][__m1][__q2].resize( __mMaxQ2 );
                }
            }
        }

        M_C0_du.resize( __QOutput );
        for( int __q1=0; __q1< __QOutput; __q1++)
        {
            int __mMaxQ1=M_model->mMaxF(M_output_index,__q1);
            M_C0_du[__q1].resize( __mMaxQ1 );
            for( int __m1=0; __m1< __mMaxQ1; __m1++)
            {
                M_C0_du[__q1][__m1].resize(  __QOutput );
                for( int __q2=0; __q2< __QOutput; __q2++)
                {
                    int __mMaxQ2=M_model->mMaxF(M_output_index,__q2);
                    M_C0_du[__q1][__m1][__q2].resize( __mMaxQ2 );
                }
            }
        }

        //typename array_3_type::extent_gen extents3;
        //M_Lambda_pr.resize( extents3[__QLhs][__QRhs] );
        //M_Lambda_du.resize( extents3[__QLhs][__QOutput] );
        M_Lambda_pr.resize( __QLhs );
        for( int __q1=0; __q1< __QLhs; __q1++)
        {
            int __mMaxQ1=M_model->mMaxA(__q1);
            M_Lambda_pr[__q1].resize( __mMaxQ1 );
            for( int __m1=0; __m1< __mMaxQ1; __m1++)
            {
                M_Lambda_pr[__q1][__m1].resize(  __QRhs );
                for( int __q2=0; __q2< __QRhs; __q2++)
                {
                    int __mMaxQ2=M_model->mMaxF(0,__q2);
                    M_Lambda_pr[__q1][__m1][__q2].resize( __mMaxQ2 );
                }
            }
        }

        M_Lambda_du.resize( __QLhs );
        for( int __q1=0; __q1< __QLhs; __q1++)
        {
            int __mMaxQ1=M_model->mMaxA(__q1);
            M_Lambda_du[__q1].resize( __mMaxQ1 );
            for( int __m1=0; __m1< __mMaxQ1; __m1++)
            {
                M_Lambda_du[__q1][__m1].resize(  __QOutput );
                for( int __q2=0; __q2< __QOutput; __q2++)
                {
                    int __mMaxQ2=M_model->mMaxF(M_output_index,__q2);
                    M_Lambda_du[__q1][__m1][__q2].resize( __mMaxQ2 );
                }
            }
        }

        //typename array_4_type::extent_gen extents4;
        //M_Gamma_pr.resize( extents4[__QLhs][__QLhs] );
        //M_Gamma_du.resize( extents4[__QLhs][__QLhs] );
        M_Gamma_pr.resize( __QLhs );
        for( int __q1=0; __q1< __QLhs; __q1++)
        {
            int __mMaxQ1=M_model->mMaxA(__q1);
            M_Gamma_pr[__q1].resize( __mMaxQ1 );
            for( int __m1=0; __m1< __mMaxQ1; __m1++)
            {
                M_Gamma_pr[__q1][__m1].resize(  __QLhs );
                for( int __q2=0; __q2< __QLhs; __q2++)
                {
                    int __mMaxQ2=M_model->mMaxA(__q2);
                    M_Gamma_pr[__q1][__m1][__q2].resize( __mMaxQ2 );
                }
            }
        }

        M_Gamma_du.resize( __QLhs );
        for( int __q1=0; __q1< __QLhs; __q1++)
        {
            int __mMaxQ1=M_model->mMaxA(__q1);
            M_Gamma_du[__q1].resize( __mMaxQ1 );
            for( int __m1=0; __m1< __mMaxQ1; __m1++)
            {
                M_Gamma_du[__q1][__m1].resize(  __QLhs );
                for( int __q2=0; __q2< __QLhs; __q2++)
                {
                    int __mMaxQ2=M_model->mMaxA(__q2);
                    M_Gamma_du[__q1][__m1][__q2].resize( __mMaxQ2 );
                }
            }
        }


        M_C0_pr_eim.resize( __QRhs );
        for( int __q1=0; __q1< __QRhs; __q1++)
        {
            M_C0_pr_eim[__q1].resize(  __QRhs );
        }

        M_C0_du_eim.resize( __QOutput );
        for( int __q1=0; __q1< __QOutput; __q1++)
        {
            M_C0_du_eim[__q1].resize(  __QOutput );
        }

        M_Lambda_pr_eim.resize( __QLhs );
        for( int __q1=0; __q1< __QLhs; __q1++)
        {
            M_Lambda_pr_eim[__q1].resize(  __QRhs );
        }

        M_Lambda_du_eim.resize( __QLhs );
        for( int __q1=0; __q1< __QLhs; __q1++)
        {
            M_Lambda_du_eim[__q1].resize(  __QOutput );
        }

        M_Gamma_pr_eim.resize( __QLhs );
        for( int __q1=0; __q1< __QLhs; __q1++)
        {
            M_Gamma_pr_eim[__q1].resize( __QLhs );
        }

        M_Gamma_du_eim.resize( __QLhs );
        for( int __q1=0; __q1< __QLhs; __q1++)
        {
            M_Gamma_du_eim[__q1].resize( __QLhs );
        }


        if ( model_type::is_time_dependent )
        {
            // M_Cmf_pr.resize( extents3[__Qm][__QRhs] );
            // M_Cmf_du.resize( extents3[__Qm][__QRhs] );
            M_Cmf_pr.resize( __Qm );
            for( int __q1=0; __q1< __Qm; __q1++)
            {
                int __mMaxQ1=M_model->mMaxM(__q1);
                M_Cmf_pr[__q1].resize( __mMaxQ1 );
                for( int __m1=0; __m1< __mMaxQ1; __m1++)
                {
                    M_Cmf_pr[__q1][__m1].resize(  __QRhs );
                    for( int __q2=0; __q2< __QRhs; __q2++)
                    {
                        int __mMaxQ2=M_model->mMaxF(0,__q2);
                        M_Cmf_pr[__q1][__m1][__q2].resize( __mMaxQ2 );
                    }
                }
            }

            M_Cmf_du.resize( __Qm );
            for( int __q1=0; __q1< __Qm; __q1++)
            {
                int __mMaxQ1=M_model->mMaxM(__q1);
                M_Cmf_du[__q1].resize( __mMaxQ1 );
                for( int __m1=0; __m1< __mMaxQ1; __m1++)
                {
                    M_Cmf_du[__q1][__m1].resize( __QOutput );
                    for( int __q2=0; __q2< __QOutput; __q2++)
                    {
                        int __mMaxQ2=M_model->mMaxF(M_output_index,__q2);
                        M_Cmf_du[__q1][__m1][__q2].resize( __mMaxQ2 );
                    }
                }
            }
            M_Cmf_du_ini.resize( __Qm );
            for( int __q1=0; __q1< __Qm; __q1++)
            {
                int __mMaxQ1=M_model->mMaxM(__q1);
                M_Cmf_du_ini[__q1].resize( __mMaxQ1 );
                for( int __m1=0; __m1< __mMaxQ1; __m1++)
                {
                    M_Cmf_du_ini[__q1][__m1].resize( __QOutput );
                    for( int __q2=0; __q2< __QOutput; __q2++)
                    {
                        int __mMaxQ2=M_model->mMaxF(M_output_index,__q2);
                        M_Cmf_du_ini[__q1][__m1][__q2].resize( __mMaxQ2 );
                    }
                }
            }

            // M_Cma_pr.resize( extents4[__Qm][__QLhs] );
            // M_Cma_du.resize( extents4[__Qm][__QLhs] );
            M_Cma_pr.resize( __Qm );
            for( int __q1=0; __q1< __Qm; __q1++)
            {
                int __mMaxQ1=M_model->mMaxM(__q1);
                M_Cma_pr[__q1].resize( __mMaxQ1 );
                for( int __m1=0; __m1< __mMaxQ1; __m1++)
                {
                    M_Cma_pr[__q1][__m1].resize(  __QLhs );
                    for( int __q2=0; __q2< __QLhs; __q2++)
                    {
                        int __mMaxQ2=M_model->mMaxA(__q2);
                        M_Cma_pr[__q1][__m1][__q2].resize( __mMaxQ2 );
                    }
                }
            }

            M_Cma_du.resize( __Qm );
            for( int __q1=0; __q1< __Qm; __q1++)
            {
                int __mMaxQ1=M_model->mMaxM(__q1);
                M_Cma_du[__q1].resize( __mMaxQ1 );
                for( int __m1=0; __m1< __mMaxQ1; __m1++)
                {
                    M_Cma_du[__q1][__m1].resize( __QLhs );
                    for( int __q2=0; __q2< __QLhs; __q2++)
                    {
                        int __mMaxQ2=M_model->mMaxA(__q2);
                        M_Cma_du[__q1][__m1][__q2].resize( __mMaxQ2 );
                    }
                }
            }

            // M_Cmm_pr.resize( extents4[__Qm][__Qm] );
            // M_Cmm_du.resize( extents4[__Qm][__Qm] );
            M_Cmm_pr.resize( __Qm );
            for( int __q1=0; __q1< __Qm; __q1++)
            {
                int __mMaxQ1=M_model->mMaxM(__q1);
                M_Cmm_pr[__q1].resize( __mMaxQ1 );
                for( int __m1=0; __m1< __mMaxQ1; __m1++)
                {
                    M_Cmm_pr[__q1][__m1].resize(  __Qm );
                    for( int __q2=0; __q2< __Qm; __q2++)
                    {
                        int __mMaxQ2=M_model->mMaxM(__q2);
                        M_Cmm_pr[__q1][__m1][__q2].resize( __mMaxQ2 );
                    }
                }
            }

            M_Cmm_du.resize( __Qm );
            for( int __q1=0; __q1< __Qm; __q1++)
            {
                int __mMaxQ1=M_model->mMaxM(__q1);
                M_Cmm_du[__q1].resize( __mMaxQ1 );
                for( int __m1=0; __m1< __mMaxQ1; __m1++)
                {
                    M_Cmm_du[__q1][__m1].resize( __Qm );
                    for( int __q2=0; __q2< __Qm; __q2++)
                    {
                        int __mMaxQ2=M_model->mMaxM(__q2);
                        M_Cmm_du[__q1][__m1][__q2].resize( __mMaxQ2 );
                    }
                }
            }

        }//end of if ( model_type::is_time_dependent )
    }//end of if ( M_error_type == CRB_RESIDUAL || M_error_type == CRB_RESIDUAL_SCM )

    LOG(INFO) << "[CRB::offline] allocate reduced basis data structures\n";
    M_Aqm_pr.resize( M_model->Qa() );
    M_Aqm_du.resize( M_model->Qa() );
    M_Aqm_pr_du.resize( M_model->Qa() );
    if( M_use_newton )
        M_Jqm_pr.resize( M_model->Qa() );

    for(int q=0; q<M_model->Qa(); q++)
    {
        M_Aqm_pr[q].resize( M_model->mMaxA(q) );
        M_Aqm_du[q].resize( M_model->mMaxA(q) );
        M_Aqm_pr_du[q].resize( M_model->mMaxA(q) );
        if(M_use_newton)
            M_Jqm_pr[q].resize( M_model->mMaxA(q) );
    }

    M_Mqm_pr.resize( M_model->Qm() );
    M_Mqm_du.resize( M_model->Qm() );
    M_Mqm_pr_du.resize( M_model->Qm() );
    for(int q=0; q<M_model->Qm(); q++)
    {
        M_Mqm_pr[q].resize( M_model->mMaxM(q) );
        M_Mqm_du[q].resize( M_model->mMaxM(q) );
        M_Mqm_pr_du[q].resize( M_model->mMaxM(q) );
    }

    int QInitialGuessV = M_model->QInitialGuess();
    M_InitialGuessV_pr.resize( QInitialGuessV );
    for(int q=0; q<QInitialGuessV; q++)
    {
        M_InitialGuessV_pr[q].resize( M_model->mMaxInitialGuess(q) );
    }

    M_Fqm_pr.resize( M_model->Ql( 0 ) );
    M_Fqm_du.resize( M_model->Ql( 0 ) );
    M_Lqm_pr.resize( M_model->Ql( M_output_index ) );
    M_Lqm_du.resize( M_model->Ql( M_output_index ) );
    if(M_use_newton)
        M_Rqm_pr.resize( M_model->Ql( 0 ) );

    for(int q=0; q<M_model->Ql( 0 ); q++)
    {
        M_Fqm_pr[q].resize( M_model->mMaxF( 0 , q) );
        M_Fqm_du[q].resize( M_model->mMaxF( 0 , q) );
        if(M_use_newton)
            M_Rqm_pr[q].resize( M_model->mMaxF( 0 , q) );
    }
    for(int q=0; q<M_model->Ql( M_output_index ); q++)
    {
        M_Lqm_pr[q].resize( M_model->mMaxF( M_output_index , q) );
        M_Lqm_du[q].resize( M_model->mMaxF( M_output_index , q) );
    }
}


template<typename TruthModelType>
void
CRB<TruthModelType>::checkInitialGuess( const element_type expansion_uN , parameter_type const& mu, vectorN_type & error) const
{
    static const bool is_composite = functionspace_type::is_composite;
    checkInitialGuess( expansion_uN, mu, error , mpl::bool_< is_composite >() );
}

template<typename TruthModelType>
void
CRB<TruthModelType>::checkInitialGuess( const element_type expansion_uN , parameter_type const& mu, vectorN_type & error, mpl::bool_<false>) const
{
    error.resize(1);
    const element_ptrtype initial_guess = M_model->assembleInitialGuess( mu );
    auto Xh = expansion_uN.functionSpace();
    auto mesh = Xh->mesh();
    error(0) = math::sqrt(
        integrate( _range=elements(mesh) ,
                   _expr=vf::idv( initial_guess ) * vf::idv( expansion_uN )
                   * vf::idv( initial_guess ) * vf::idv( expansion_uN )
                   ).evaluate()(0,0)
                          );


}

template<typename TruthModelType>
void
CRB<TruthModelType>::checkInitialGuess( const element_type expansion_uN , parameter_type const& mu, vectorN_type & error, mpl::bool_<true>) const
{
    //using namespace Feel::vf;
    index_vector_type index_vector;
    const element_ptrtype initial_guess = M_model->assembleInitialGuess( mu );
    fusion::for_each( index_vector , ComputeIntegralsSquare( *initial_guess , expansion_uN, error ) );
}


template<typename TruthModelType>
void
CRB<TruthModelType>::buildVarianceMatrixPhi( int const N )
{
    static const bool is_composite = functionspace_type::is_composite;
    return buildVarianceMatrixPhi( N , mpl::bool_< is_composite >() );
}
template<typename TruthModelType>
void
CRB<TruthModelType>::buildVarianceMatrixPhi( int const N , mpl::bool_<true> )
{

    std::vector<std::string> s;
    static const int nb_spaces = functionspace_type::nSpaces;

    if( N == 1 )
        M_variance_matrix_phi.resize( nb_spaces );


    for(int i=0;i<nb_spaces;i++)
        M_variance_matrix_phi[i].conservativeResize( M_N , M_N );

    //introduction of the new variable phi
    //int size_WN = M_WN.size();
    int size_WN = M_model->rBFunctionSpace()->size();
    wn_type phi;

    for(int i = 0; i < size_WN; ++i)
    {
        //auto sub = subelements( M_WN[i] , s );
        auto sub = subelements( M_model->rBFunctionSpace()->primalBasisElement(i) , s );
        element_type global_element = M_model->functionSpace()->element();
        wn_type vect( nb_spaces , global_element);
        fusion::for_each( sub , ComputePhi( vect ) );

        //now we want to have only one element_type (global_element)
        //which sum the contribution of each space
        global_element.zero();
        for( auto element : vect )
            global_element += element;
        phi.push_back( global_element );
    }


    index_vector_type index_vector;
    for(int space=0; space<nb_spaces; space++)
        M_variance_matrix_phi[space].conservativeResize( N , N );


    std::vector<double> vectIntegrals;
    for(int i = 0; i < M_N; ++i)
    {
        for(int j = i+1; j < M_N; ++j)
        {
            vectIntegrals.clear();
            fusion::for_each( index_vector , ComputeIntegrals( phi[i] , phi[j], vectIntegrals ) );
            for(int space=0; space<nb_spaces; space++)
                M_variance_matrix_phi[space](i,j)=vectIntegrals[space];
        } //j
        vectIntegrals.clear();
        fusion::for_each( index_vector, ComputeIntegrals(  phi[i] , phi[i], vectIntegrals ) );
        for(int space=0; space<nb_spaces; space++)
            M_variance_matrix_phi[space](i,i)=vectIntegrals[space];
    }// i

}


template<typename TruthModelType>
void
CRB<TruthModelType>::buildVarianceMatrixPhi( int const N , mpl::bool_<false> )
{
    std::vector<std::string> s;
    int nb_spaces = functionspace_type::nSpaces;

    if( N == 1 )
        M_variance_matrix_phi.resize( nb_spaces );

    M_variance_matrix_phi[0].conservativeResize( M_N , M_N );

    //introduction of the new variable phi
    int size_WN = M_model->rBFunctionSpace()->size();
    //wn_type phi ( size_WN );
    wn_type phi;

    mesh_ptrtype mesh = M_model->rBFunctionSpace()->primalBasisElement(0).functionSpace()->mesh();

    for(int i = 0; i < size_WN; ++i)
    {
        double surface = integrate( _range=elements(mesh), _expr=vf::cst(1.) ).evaluate()(0,0);
        double mean =  integrate( _range=elements(mesh), _expr=vf::idv( M_model->rBFunctionSpace()->primalBasisElement(i) ) ).evaluate()(0,0);
        mean /= surface;
        auto element_mean = vf::project(M_model->rBFunctionSpace()->primalBasisElement(i).functionSpace(), elements(mesh), vf::cst(mean) );
        phi.push_back( M_model->rBFunctionSpace()->primalBasisElement(i) - element_mean );
    }

    for(int space=0; space<nb_spaces; space++)
        M_variance_matrix_phi[space].conservativeResize( N , N );

    for(int i = 0; i < M_N; ++i)
    {
        for(int j = i+1; j < M_N; ++j)
        {
            M_variance_matrix_phi[0]( i , j ) = integrate( _range=elements(mesh) , _expr=vf::idv( phi[i] ) * vf::idv( phi[j] ) ).evaluate()(0,0);
            M_variance_matrix_phi[0]( j , i ) =  M_variance_matrix_phi[0]( i , j );
        } //j

        M_variance_matrix_phi[0]( i , i ) = integrate( _range=elements(mesh) , _expr=vf::idv( phi[i] ) * vf::idv( phi[i] ) ).evaluate()(0,0);
    }// i

}



template<typename TruthModelType>
void
CRB<TruthModelType>::buildFunctionFromRbCoefficients(int N, std::vector< vectorN_type > const & RBcoeff, wn_type const & WN, std::vector<element_ptrtype> & FEMsolutions )
{

    if( WN.size() == 0 )
        throw std::logic_error( "[CRB::buildFunctionFromRbCoefficients] ERROR : reduced basis space is empty" );


    int nb_solutions = RBcoeff.size();

    for( int i = 0; i < nb_solutions; i++ )
    {
        element_ptrtype FEMelement ( new element_type( M_model->functionSpace() ) );
        FEMelement->setZero();
        for( int j = 0; j < N; j++ )
            FEMelement->add( RBcoeff[i](j) , WN[j] );
        FEMsolutions.push_back( FEMelement );
    }
}

template<typename TruthModelType>
void
CRB<TruthModelType>::compareResidualsForTransientProblems( int N, parameter_type const& mu, std::vector<element_type> const & Un, std::vector<element_type> const & Unold,
                                                           std::vector<element_type> const& Undu, std::vector<element_type> const & Unduold,
                                                           std::vector< std::vector<double> > const& primal_residual_coeffs,  std::vector < std::vector<double> > const& dual_residual_coeffs ) const
{

    LOG( INFO ) <<"\n compareResidualsForTransientProblems \n";
    //backend_ptrtype backend = backend_type::build( BACKEND_PETSC, this->worldComm() ) ;

    if ( M_model->isSteady() )
    {
        throw std::logic_error( "[CRB::compareResidualsForTransientProblems] ERROR : to check residual in a steady case, use checkResidual and not compareResidualsForTransientProblems" );
    }

    sparse_matrix_ptrtype A,AM,M,Adu;
    //vector_ptrtype MF;
    std::vector<vector_ptrtype> F,L;

    vector_ptrtype Rhs( backend()->newVector( M_model->functionSpace() ) );
    vector_ptrtype Aun( backend()->newVector( M_model->functionSpace() ) );
    vector_ptrtype AduUn( backend()->newVector( M_model->functionSpace() ) );

    vector_ptrtype Mun( backend()->newVector( M_model->functionSpace() ) );
    vector_ptrtype Munold( backend()->newVector( M_model->functionSpace() ) );
    vector_ptrtype Frhs( backend()->newVector( M_model->functionSpace() ) );
    vector_ptrtype un( backend()->newVector( M_model->functionSpace() ) );
    vector_ptrtype unold( backend()->newVector( M_model->functionSpace() ) );
    vector_ptrtype undu( backend()->newVector( M_model->functionSpace() ) );
    vector_ptrtype unduold( backend()->newVector( M_model->functionSpace() ) );

    //set parameters for time discretization
    auto bdf_primal = bdf( _space=M_model->functionSpace(), _vm=this->vm() , _name="bdf_primal_check_residual_transient" );
    bdf_primal->setTimeInitial( M_model->timeInitial() );
    bdf_primal->setTimeStep( M_model->timeStep() );
    bdf_primal->setTimeFinal( M_model->timeFinal() );
    bdf_primal->setOrder( M_model->timeOrder() );


    double time_step = bdf_primal->timeStep();
    int time_index=0;
    for ( bdf_primal->start() ; !bdf_primal->isFinished(); bdf_primal->next() )
    {

        auto bdf_poly = bdf_primal->polyDeriv();
        boost::tie( M, A, F) = M_model->update( mu , bdf_primal->time() );

        A->close();

        *un = Un[time_index];
        *unold = Unold[time_index];

        A->multVector( un, Aun );
        M->multVector( un, Mun );

        M->multVector( unold, Munold );
        Aun->scale( -1 );
        Mun->scale( -1 );
        *Frhs = *F[0];

        vector_ptrtype __ef_pr(  backend()->newVector( M_model->functionSpace() ) );
        vector_ptrtype __ea_pr(  backend()->newVector( M_model->functionSpace() ) );
        vector_ptrtype __emu_pr(  backend()->newVector( M_model->functionSpace() ) );
        vector_ptrtype __emuold_pr(  backend()->newVector( M_model->functionSpace() ) );
        M_model->l2solve( __ef_pr, Frhs );
        M_model->l2solve( __ea_pr, Aun );
        M_model->l2solve( __emu_pr, Mun );
        M_model->l2solve( __emuold_pr, Munold );

        double check_Cff_pr = M_model->scalarProduct( __ef_pr,__ef_pr );
        double check_Caf_pr = 2*M_model->scalarProduct( __ef_pr,__ea_pr );
        double check_Caa_pr = M_model->scalarProduct( __ea_pr,__ea_pr );
        double check_Cmf_pr = 2./time_step*( M_model->scalarProduct( __emu_pr , __ef_pr )+M_model->scalarProduct( __emuold_pr , __ef_pr ) );
        double check_Cma_pr = 2./time_step*( M_model->scalarProduct( __emu_pr , __ea_pr )+M_model->scalarProduct( __emuold_pr , __ea_pr ) );
        double check_Cmm_pr = 1./(time_step*time_step)*( M_model->scalarProduct( __emu_pr , __emu_pr ) + 2*M_model->scalarProduct( __emu_pr , __emuold_pr ) + M_model->scalarProduct( __emuold_pr , __emuold_pr ) );

        double Cff_pr = primal_residual_coeffs[time_index][0];
        double Caf_pr = primal_residual_coeffs[time_index][1];
        double Caa_pr = primal_residual_coeffs[time_index][2];
        double Cmf_pr = primal_residual_coeffs[time_index][3];
        double Cma_pr = primal_residual_coeffs[time_index][4];
        double Cmm_pr = primal_residual_coeffs[time_index][5];

        LOG(INFO)<<" --- time : "<<bdf_primal->time();
        LOG(INFO)<<"Cff : "<< check_Cff_pr <<"  -  "<<Cff_pr<<"  =>  "<<check_Cff_pr-Cff_pr;
        LOG(INFO)<<"Caf : "<< check_Caf_pr <<"  -  "<<Caf_pr<<"  =>  "<<check_Caf_pr-Caf_pr;
        LOG(INFO)<<"Caa : "<< check_Caa_pr <<"  -  "<<Caa_pr<<"  =>  "<<check_Caa_pr-Caa_pr;
        LOG(INFO)<<"Cmf : "<< check_Cmf_pr <<"  -  "<<Cmf_pr<<"  =>  "<<check_Cmf_pr-Cmf_pr;
        LOG(INFO)<<"Cma : "<< check_Cma_pr <<"  -  "<<Cma_pr<<"  =>  "<<check_Cma_pr-Cma_pr;
        LOG(INFO)<<"Cmm : "<< check_Cmm_pr <<"  -  "<<Cmm_pr<<"  =>  "<<check_Cmm_pr-Cmm_pr;
        time_index++;
    }

    time_index--;

    //bool solve_dual_problem = boption(_prefix=M_prefix,_name="crb.solve-dual-problem");
    //if( this->worldComm().globalSize() > 1 )
    //    solve_dual_problem=false;

    double sum=0;
    if( M_solve_dual_problem )
    {

        Adu = M_model->newMatrix();

        auto bdf_dual = bdf( _space=M_model->functionSpace(), _vm=this->vm() , _name="bdf_dual_check_residual_transient" );

        bdf_dual->setTimeInitial( M_model->timeFinal()+M_model->timeStep() );
        bdf_dual->setTimeStep( -M_model->timeStep() );
        bdf_dual->setTimeFinal( M_model->timeInitial()+M_model->timeStep() );
        bdf_dual->setOrder( M_model->timeOrder() );
        bdf_dual->start();

        //element_ptrtype dual_initial_field( new element_type( M_model->functionSpace() ) );
        LOG(INFO)<<"**********dual problem************* "<<std::endl;

        boost::tie( M, A, F) = M_model->update( mu , bdf_dual->timeInitial() );

        vectorN_type dual_initial ( N );
        for(int i=0; i<N; i++)
            dual_initial(i) = M_coeff_du_ini_online(i);
        auto dual_initial_field = this->expansion( dual_initial , N , true );

        *undu = dual_initial_field;
        M->multVector( undu, Mun );
        *Frhs = *F[0];


        auto R = backend()->newVector( M_model->functionSpace() );
        R = Frhs;
        R->add( -1 , *Mun ); //R -= Mun;
        //std::cout<<"[COMPARE] R->l2Norm() : "<<R->l2Norm()<<std::endl;

        vector_ptrtype __ef_du(  backend()->newVector( M_model->functionSpace() ) );
        vector_ptrtype __emu_du(  backend()->newVector( M_model->functionSpace() ) );
        M_model->l2solve( __ef_du, Frhs );
        M_model->l2solve( __emu_du, Mun );
        double check_Cff_du = M_model->scalarProduct( __ef_du,__ef_du );
        double check_Cmf_du = 2*M_model->scalarProduct( __ef_du,__emu_du );
        double check_Cmm_du = M_model->scalarProduct( __emu_du,__emu_du );
        double residual_final_condition = math::abs( check_Cff_du + check_Cmf_du + check_Cmm_du );
        //std::cout<<"residual on final condition : "<<residual_final_condition<<std::endl;
        //initialization
        time_step = bdf_dual->timeStep();


        for ( bdf_dual->start(); !bdf_dual->isFinished() ; bdf_dual->next() )
        {
            auto bdf_poly = bdf_dual->polyDeriv();

            boost::tie( M, A, F ) = M_model->update( mu , bdf_dual->time() );
            if( boption(_prefix=M_prefix,_name="crb.use-symmetric-matrix") )
                Adu = A;
            else
                A->transpose( Adu );
            *undu = Undu[time_index];
            *unduold = Unduold[time_index];
            Adu->multVector( undu, AduUn );
            M->multVector( undu, Mun );
            M->multVector( unduold, Munold );
            AduUn->scale( -1 );
            Munold->scale( -1 );
            *Frhs = *F[0];

            vector_ptrtype __ea_du(  backend()->newVector( M_model->functionSpace() ) );
            vector_ptrtype __emu_du(  backend()->newVector( M_model->functionSpace() ) );
            vector_ptrtype __emuold_du(  backend()->newVector( M_model->functionSpace() ) );
            M_model->l2solve( __ea_du, AduUn );
            M_model->l2solve( __emu_du, Mun );
            M_model->l2solve( __emuold_du, Munold );
            double check_Caa_du = M_model->scalarProduct( __ea_du,__ea_du );
            double check_Cma_du = 2./time_step*( M_model->scalarProduct( __emu_du , __ea_du )+M_model->scalarProduct( __emuold_du , __ea_du ) );
            double check_Cmm_du = 1./(time_step*time_step)*( M_model->scalarProduct( __emu_du , __emu_du ) + 2*M_model->scalarProduct( __emu_du , __emuold_du ) + M_model->scalarProduct( __emuold_du , __emuold_du ) );


            double Cff_du =  dual_residual_coeffs[time_index][0];
            double Caf_du =  dual_residual_coeffs[time_index][1];
            double Caa_du =  dual_residual_coeffs[time_index][2];
            double Cmf_du =  dual_residual_coeffs[time_index][3];
            double Cma_du =  dual_residual_coeffs[time_index][4];
            double Cmm_du =  dual_residual_coeffs[time_index][5];
            LOG(INFO)<<" --- time : "<<bdf_dual->time()<<std::endl;
            LOG(INFO)<<"Caa : "<< check_Caa_du <<"  -  "<<Caa_du<<"  =>  "<<check_Caa_du-Caa_du<<std::endl;
            LOG(INFO)<<"Cma : "<< check_Cma_du <<"  -  "<<Cma_du<<"  =>  "<<check_Cma_du-Cma_du<<std::endl;
            LOG(INFO)<<"Cmm : "<< check_Cmm_du <<"  -  "<<Cmm_du<<"  =>  "<<check_Cmm_du-Cmm_du<<std::endl;
            time_index--;
            //std::cout<<"[CHECK] ------ time "<<bdf_dual->time()<<std::endl;
            //std::cout<<"[CHECK] Caa_du : "<<check_Caa_du<<std::endl;
            //std::cout<<"[CHECK] Cma_du : "<<check_Cma_du<<std::endl;
            //std::cout<<"[CHECK] Cmm_du : "<<check_Cmm_du<<std::endl;
            sum += math::abs( check_Caa_du + check_Cma_du + check_Cmm_du );
        }
        //std::cout<<"[CHECK] dual_sum : "<<sum<<std::endl;
    }//solve-dual-problem
}


template<typename TruthModelType>
void
CRB<TruthModelType>::checkResidual( parameter_type const& mu, std::vector< std::vector<double> > const& primal_residual_coeffs,
                                    std::vector< std::vector<double> > const& dual_residual_coeffs , element_type & u, element_type & udu ) const
{


    if ( 0 )// orthonormalize_primal || orthonormalize_dual )
    {
        throw std::logic_error( "[CRB::checkResidual] ERROR : to check residual don't use orthonormalization" );
    }

    if ( !M_model->isSteady() )
    {
        throw std::logic_error( "[CRB::checkResidual] ERROR : to check residual select steady state" );
    }

    if ( M_error_type==CRB_NO_RESIDUAL || M_error_type==CRB_EMPIRICAL )
    {
        throw std::logic_error( "[CRB::checkResidual] ERROR : to check residual set option crb.error-type to 0 or 1" );
    }

    int size = mu.size();

    if ( 0 )
    {
        std::cout<<"[CRB::checkResidual] use mu = [";
        for ( int i=0; i<size-1; i++ ) std::cout<< mu[i] <<" , ";
        std::cout<< mu[size-1]<<" ]"<<std::endl;
    }

    LOG( INFO ) <<"[CRB::checkResidual] use mu = \n"<<mu;
    sparse_matrix_ptrtype A,At,M;
    std::vector<vector_ptrtype> F,L;

    //backend_ptrtype backendA = backend_type::build( BACKEND_PETSC );
    //backend_ptrtype backendAt = backend_type::build( BACKEND_PETSC );

    //element_ptrtype u( new element_type( M_model->functionSpace() ) );
    //element_ptrtype udu( new element_type( M_model->functionSpace() ) );
    vector_ptrtype U( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype Rhs( M_backend->newVector( M_model->functionSpace() ) );

    Feel::Timer timer, timer2;

    boost::tie( boost::tuples::ignore, A, F ) = M_model->update( mu );

    LOG(INFO) << "  -- updated model for parameter in " << timer2.elapsed() << "s\n";
    timer2.start();

    LOG(INFO) << "[CRB::checkResidual] transpose primal matrix" << "\n";
    At = M_model->newMatrix();
    if( boption(_prefix=M_prefix,_name="crb.use-symmetric-matrix") )
        At = A;
    else
        A->transpose( At );
    //u->setName( ( boost::format( "fem-primal-%1%" ) % ( M_N ) ).str() );
    //udu->setName( ( boost::format( "fem-dual-%1%" ) % ( M_N ) ).str() );

    //LOG(INFO) << "[CRB::checkResidual] solving primal" << "\n";
    //backendA->solve( _matrix=A,  _solution=u, _rhs=F[0] );
    //LOG(INFO) << "  -- primal problem solved in " << timer2.elapsed() << "s\n";
    timer2.start();
    *Rhs = *F[M_output_index];
    Rhs->close();
    Rhs->scale( -1 );
    //backendAt->solve( _matrix=At,  _solution=udu, _rhs=Rhs );
    //LOG(INFO) << "  -- dual problem solved in " << timer2.elapsed() << "s\n";
    timer2.start();

    vector_ptrtype Aun( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype Atun( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype Un( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype Undu( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype Frhs( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype Lrhs( M_backend->newVector( M_model->functionSpace() ) );
    *Un = u;
    *Undu = udu;
    A->multVector( Un, Aun );
    At->multVector( Undu, Atun );
    Aun->close();
    Atun->close();
    Aun->scale( -1 );
    Atun->scale( -1 );
    *Frhs = *F[0];
    *Lrhs = *F[M_output_index];

#if 0
    LOG(INFO) << "[CRB::checkResidual] residual (f,f) " << M_N-1 << ":=" << M_model->scalarProduct( Frhs, Frhs ) << "\n";
    LOG(INFO) << "[CRB::checkResidual] residual (f,A) " << M_N-1 << ":=" << 2*M_model->scalarProduct( Frhs, Aun ) << "\n";
    LOG(INFO) << "[CRB::checkResidual] residual (A,A) " << M_N-1 << ":=" << M_model->scalarProduct( Aun, Aun ) << "\n";

    LOG(INFO) << "[CRB::checkResidual] residual (l,l) " << M_N-1 << ":=" << M_model->scalarProduct( Lrhs, Lrhs ) << "\n";
    LOG(INFO) << "[CRB::checkResidual] residual (l,At) " << M_N-1 << ":=" << 2*M_model->scalarProduct( Lrhs, Atun ) << "\n";
    LOG(INFO) << "[CRB::checkResidual] residual (At,At) " << M_N-1 << ":=" << M_model->scalarProduct( Atun, Atun ) << "\n";
#endif

    Lrhs->close();
    Lrhs->scale( -1 );

    vector_ptrtype __ef_pr(  M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __ea_pr(  M_backend->newVector( M_model->functionSpace() ) );
    M_model->l2solve( __ef_pr, Frhs );
    M_model->l2solve( __ea_pr, Aun );
    double check_C0_pr = M_model->scalarProduct( __ef_pr,__ef_pr );
    double check_Lambda_pr = 2*M_model->scalarProduct( __ef_pr,__ea_pr );
    double check_Gamma_pr = M_model->scalarProduct( __ea_pr,__ea_pr );

    vector_ptrtype __ef_du(  M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __ea_du(  M_backend->newVector( M_model->functionSpace() ) );
    M_model->l2solve( __ef_du, Lrhs );
    M_model->l2solve( __ea_du, Atun );
    double check_C0_du = M_model->scalarProduct( __ef_du,__ef_du );
    double check_Lambda_du = 2*M_model->scalarProduct( __ef_du,__ea_du );
    double check_Gamma_du = M_model->scalarProduct( __ea_du,__ea_du );

    double primal_sum =  check_C0_pr + check_Lambda_pr + check_Gamma_pr ;
    double dual_sum   =  check_C0_du + check_Lambda_du + check_Gamma_du ;

    //LOG(INFO)<<"[CRB::checkResidual] primal_sum = "<<check_C0_pr<<" + "<<check_Lambda_pr<<" + "<<check_Gamma_pr<<" = "<<primal_sum<<"\n";
    //LOG(INFO)<<"[CRB::checkResidual] dual_sum = "<<check_C0_du<<" + "<<check_Lambda_du<<" + "<<check_Gamma_du<<" = "<<dual_sum<<"\n";


    int time_index=0;
    double C0_pr = primal_residual_coeffs[time_index][0];
    double Lambda_pr = primal_residual_coeffs[time_index][1];
    double Gamma_pr = primal_residual_coeffs[time_index][2];

    double C0_du,Lambda_du,Gamma_du;
    C0_du = dual_residual_coeffs[time_index][0];
    Lambda_du = dual_residual_coeffs[time_index][1];
    Gamma_du = dual_residual_coeffs[time_index][2];

    double err_C0_pr = math::abs( C0_pr - check_C0_pr ) ;
    double err_Lambda_pr = math::abs( Lambda_pr - check_Lambda_pr ) ;
    double err_Gamma_pr = math::abs( Gamma_pr - check_Gamma_pr ) ;
    double err_C0_du = math::abs( C0_du - check_C0_du ) ;
    double err_Lambda_du = math::abs( Lambda_du - check_Lambda_du ) ;
    double err_Gamma_du = math::abs( Gamma_du - check_Gamma_du ) ;

    int start_dual_index = 6;

    if( math::abs(check_C0_pr - C0_pr)>1e-14 )
    {
        LOG( INFO )<<std::setprecision( 15 )<<" C0_pr without decomposition : "<<check_C0_pr;
        LOG( INFO )<<std::setprecision( 15 )<<" C0_pr with decomposition    : "<<C0_pr;
    }
    if( math::abs(check_C0_pr - C0_pr)>1e-14 )
    {
        LOG( INFO )<<std::setprecision( 15 )<<" Lambda_pr without decomposition : "<<check_Lambda_pr;
        LOG( INFO )<<std::setprecision( 15 )<<" Lamnda_pr with decomposition    : "<<Lambda_pr;
    }
    if( math::abs(check_Gamma_pr - Gamma_pr)>1e-14 )
    {
        LOG( INFO )<<std::setprecision( 15 )<<" Gamma_pr without decomposition : "<<check_Gamma_pr;
        LOG( INFO )<<std::setprecision( 15 )<<" Gamma_pr with decomposition    : "<<Gamma_pr;
    }
    if( math::abs(check_C0_du - C0_du)>1e-14 )
    {
        LOG( INFO )<<std::setprecision( 15 )<<" C0_du without decomposition : "<<check_C0_du;
        LOG( INFO )<<std::setprecision( 15 )<<" C0_du with decomposition    : "<<C0_du;
    }
    if( math::abs(check_C0_pr - C0_pr)>1e-14 )
    {
        LOG( INFO )<<std::setprecision( 15 )<<" Lambda_du without decomposition : "<<check_Lambda_du;
        LOG( INFO )<<std::setprecision( 15 )<<" Lamnda_du with decomposition    : "<<Lambda_du;
    }
    if( math::abs(check_Gamma_pr - Gamma_pr)>1e-14 )
    {
        LOG( INFO )<<std::setprecision( 15 )<<" Gamma_du without decomposition : "<<check_Gamma_du;
        LOG( INFO )<<std::setprecision( 15 )<<" Gamma_du with decomposition    : "<<Gamma_du;
    }

#if 0
    LOG(INFO)<<"[CRB::checkResidual]";
    LOG(INFO)<<std::setprecision( 16 )<<"prima sum ( without decomposition ) = "<<check_C0_pr+check_Lambda_pr+check_Gamma_pr;
    LOG(INFO)<<std::setprecision( 16 )<<"primal sum with decomposition       = "<<C0_pr+Lambda_pr+Gamma_pr;
    LOG(INFO)<<std::setprecision( 16 )<<"dual sum ( without decomposition )  = "<<check_C0_du+check_Lambda_du+check_Gamma_du;
    LOG(INFO)<<std::setprecision( 16 )<<"dual sum with decomposition         = "<<C0_du+Lambda_du+Gamma_du;

    LOG(INFO)<<"errors committed on coefficients ";
    LOG(INFO)<<std::setprecision( 16 )<<"C0_pr : "<<err_C0_pr<<"\tLambda_pr : "<<err_Lambda_pr<<"\tGamma_pr : "<<err_Gamma_pr;
    LOG(INFO)<<std::setprecision( 16 )<<"C0_du : "<<err_C0_du<<"\tLambda_du : "<<err_Lambda_du<<"\tGamma_du : "<<err_Gamma_du;
    LOG(INFO)<<"and now relative error : ";
    double errC0pr = err_C0_pr/check_C0_pr;
    double errLambdapr = err_Lambda_pr/check_Lambda_pr;
    double errGammapr = err_Gamma_pr/check_Gamma_pr;
    double errC0du = err_C0_du/check_C0_du;
    double errLambdadu = err_Lambda_pr/check_Lambda_pr;
    double errGammadu = err_Gamma_pr/check_Gamma_pr;
    LOG(INFO)<<std::setprecision( 16 )<<errC0pr<<"\t"<<errLambdapr<<"\t"<<errGammapr;
    LOG(INFO)<<std::setprecision( 16 )<<errC0du<<"\t"<<errLambdadu<<"\t"<<errGammadu;
#endif

#if 1
    auto primal_residual = Frhs ;
    primal_residual->add( *Aun );
    auto dual_residual = Lrhs ;
    dual_residual->add( *Atun );

    vector_ptrtype __e_pr(  M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __e_du(  M_backend->newVector( M_model->functionSpace() ) );
    M_model->l2solve( __e_pr, primal_residual );
    M_model->l2solve( __e_du, dual_residual );

    double check_dual_norm_pr = math::sqrt( M_model->scalarProduct( __e_pr,__e_pr ) );
    double dual_norm_pr = math::sqrt( math::abs( C0_pr + Lambda_pr + Gamma_pr) );
    double check_dual_norm_du = math::sqrt( M_model->scalarProduct( __e_du,__e_du ) );
    double dual_norm_du = math::sqrt( math::abs(C0_du + Lambda_du + Gamma_du) );

    double sum_pr =   C0_pr + Lambda_pr + Gamma_pr ;
    double sum_du =   C0_du + Lambda_du + Gamma_du  ;
    double check_sum_pr =   check_C0_pr + check_Lambda_pr + check_Gamma_pr ;
    double check_sum_du =   check_C0_du + check_Lambda_du + check_Gamma_du ;

    double alpha=1;
    double t;
    if ( M_error_type == CRB_RESIDUAL_SCM )
    {
        M_scmA->setScmForMassMatrix( false );
        boost::tie( alpha, t ) = M_scmA->lb( mu );
    }
    LOG( INFO )<< "primal sum without affine decomposition (using ef_pr and ea_pr) : "<<std::setprecision( 15 )<<check_sum_pr;
    LOG( INFO )<< "primal sum without affine decomposition (using only e_pr)       : "<<std::setprecision( 15 )<<M_model->scalarProduct( __e_pr,__e_pr );
    LOG( INFO )<< "primal sum with affine decomposition  (from CRB code)           : "<<std::setprecision( 15 )<<sum_pr;
    LOG( INFO )<< "dual sum without affine decomposition (using ef_du and ea_du)   : "<<std::setprecision( 15 )<<check_sum_du;
    LOG( INFO )<< "dual sum without affine decomposition (using only e_du)         : "<<std::setprecision( 15 )<<M_model->scalarProduct( __e_du,__e_du );
    LOG( INFO )<< "dual sum with affine decomposition  (from CRB code)             : "<<std::setprecision( 15 )<<sum_du;
    LOG( INFO )<< " ----------------------- ";
    LOG( INFO )<< "dual norm of primal residual without affine decomposition (e_pr)            : "<<std::setprecision( 15 )<<check_dual_norm_pr;
    LOG( INFO )<< "dual norm of primal residual without affine decomposition (ef_pr and ea_pr) : "<<std::setprecision( 15 )<<math::sqrt(check_sum_pr);
    LOG( INFO )<< "dual norm of primal residual with affine decomposition                      : "<<std::setprecision( 15 )<<dual_norm_pr;
    LOG( INFO )<< "dual norm of dual residual without affine decomposition (e_du)            : "<<std::setprecision( 15 )<<check_dual_norm_du;
    LOG( INFO )<< "dual norm of dual residual without affine decomposition (ef_du and ea_du) : "<<std::setprecision( 15 )<<math::sqrt(check_sum_du);
    LOG( INFO )<< "dual norm of dual residual with affine decomposition                      : "<<std::setprecision( 15 )<<dual_norm_du;
    LOG( INFO )<< "Primal Error Estimator without affine decomposition (e_pr)            : "<<std::setprecision( 15 )<<check_dual_norm_pr/alpha;
    LOG( INFO )<< "Primal Error Estimator without affine decomposition (ef_pr and ea_pr) : "<<std::setprecision( 15 )<<math::sqrt(check_sum_pr)/alpha;
    LOG( INFO )<< "Primal Error Estimator with affine decomposition                      : "<<std::setprecision( 15 )<<dual_norm_pr/alpha;
    LOG( INFO )<< "Dual Error Estimator without affine decomposition (e_du)            : "<<std::setprecision( 15 )<<check_dual_norm_du/alpha;
    LOG( INFO )<< "Dual Error Estimator without affine decomposition (ef_du and ea_du) : "<<std::setprecision( 15 )<<math::sqrt(check_sum_du)/alpha;
    LOG( INFO )<< "Dual Error Estimator with affine decomposition                      : "<<std::setprecision( 15 )<<dual_norm_du/alpha;
    double epsilon =2.22e-16; //machine precision
    LOG( INFO )<<"Precision bound for primal error estimator : "<< C0_pr/alpha*math::sqrt( epsilon);
    LOG( INFO )<<"Precision bound for dual error estimator : "<< C0_du/alpha*math::sqrt( epsilon);
#endif
#if 0

    //residual r(v)
    Aun->add( *Frhs );
    //Lrhs->scale( -1 );
    Atun->add( *Lrhs );
    //to have dual norm of residual we need to solve ( e , v ) = r(v) and then it's given by ||e||
    vector_ptrtype __e_pr(  M_backend->newVector( M_model->functionSpace() ) );
    M_model->l2solve( __e_pr, Aun );
    double dual_norm_pr = math::sqrt ( M_model->scalarProduct( __e_pr,__e_pr ) );
    LOG(INFO)<<"[CRB::checkResidual] dual norm of primal residual without isolate terms (c0,lambda,gamma) = "<<dual_norm_pr<<"\n";
    //idem for the dual equation
    vector_ptrtype __e_du(  M_backend->newVector( M_model->functionSpace() ) );
    M_model->l2solve( __e_du, Atun );
    double dual_norm_du = math::sqrt ( M_model->scalarProduct( __e_du,__e_du ) );
    LOG(INFO) <<"[CRB::checkResidual] dual norm of dual residual without isolate terms = "<<dual_norm_du<<"\n";

    double err_primal = math::sqrt ( M_model->scalarProduct( Aun, Aun ) );
    double err_dual = math::sqrt ( M_model->scalarProduct( Atun, Atun ) );
    LOG(INFO) << "[CRB::checkResidual] true primal residual for reduced basis function " << M_N-1 << ":=" << err_primal << "\n";
    LOG(INFO) << "[CRB::checkResidual] true dual residual for reduced basis function " << M_N-1 << ":=" << err_dual << "\n";
#endif
}


template<typename TruthModelType>
void
CRB<TruthModelType>::check( size_type N ) const
{

    if ( ioption(_prefix=M_prefix,_name="crb.check.rb") == 0)
        return;

    std::cout << "  -- check reduced basis\n";


    LOG(INFO) << "----------------------------------------------------------------------\n";

    // check that for each mu associated to a basis function of \f$W_N\f$
    //for( int k = std::max(0,(int)N-2); k < N; ++k )
    for ( size_type k = 0; k < N; ++k )
    {
        LOG(INFO) << "**********************************************************************\n";
        parameter_type const& mu = M_WNmu->at( k );
        std::vector< vectorN_type > uN; //uN.resize( N );
        std::vector< vectorN_type > uNdu; //( N );
        std::vector< vectorN_type > uNold;
        std::vector< vectorN_type > uNduold;
        auto tuple = lb( N, mu, uN, uNdu, uNold, uNduold );
        auto output_vector=tuple.template get<0>();
        double output_vector_size=output_vector.size();
        double s = output_vector[output_vector_size-1];
        auto error_estimation = delta( N, mu, uN, uNdu , uNold, uNduold );
        auto vector_err = error_estimation.template get<0>();
        int size=vector_err.size();
        double err = vector_err[size-1];

        for (int i=0;i< uN[uN.size()-1].size();++i)
            std::cout << "uN(" << i << ") =" << uN[uN.size()-1](i) << "\n";

#if 0
        //if (  err > 1e-5 )
        // {
        std::cout << "[check] error bounds are not < 1e-10\n";
        std::cout << "[check] k = " << k << "\n";
        std::cout << "[check] mu = " << mu << "\n";
        std::cout << "[check] delta = " <<  err << "\n";
        std::cout << "[check] uN( " << k << " ) = " << uN( k ) << "\n";
#endif
        // }
        element_type u_fem; bool need_to_solve=false;
        u_fem = M_model->solveFemUsingOfflineEim ( mu );
        double sfem = M_model->output( M_output_index, mu , u_fem , need_to_solve );
        size = mu.size();
        std::cout<<"    o mu = [ ";

        for ( int i=0; i<size-1; i++ ) std::cout<< mu[i] <<" , ";

        std::cout<< mu[size-1]<<" ]"<<std::endl;

        LOG(INFO) << "[check] s= " << s << " +- " << err  << " | sfem= " << sfem << " | abs(sfem-srb) =" << math::abs( sfem - s ) << "\n";
        std::cout <<"[check] s = " << s << " +- " << err  << " | sfem= " << sfem << " | abs(sfem-srb) =" << math::abs( sfem - s )<< "\n";

    }

    LOG(INFO) << "----------------------------------------------------------------------\n";

}


template< typename TruthModelType>
double
CRB<TruthModelType>::correctionTerms(parameter_type const& mu, std::vector< vectorN_type > const & uN, std::vector< vectorN_type > const & uNdu,  std::vector<vectorN_type> const & /*uNold*/, int const k ) const
{
    int N = uN[0].size();

    matrixN_type Aprdu ( (int)N, (int)N ) ;
    vectorN_type Fdu ( (int)N );
    vectorN_type du ( (int)N );
    vectorN_type pr ( (int)N );

    beta_vector_type betaAqm;
    std::vector<beta_vector_type> betaFqm;

    bool is_linear = M_model->isLinear();

    double correction=0;

    if( M_model->isSteady() )
    {
        Aprdu.setZero( N , N );
        Fdu.setZero( N );

#if 0
        boost::tie( betaMqm, betaAqm, betaFqm ) = M_model->computeBetaQm( mu ,time);
#else
        if ( is_linear )
            boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) = M_model->computeBetaQm( mu/* ,time*/);
        else
        {
            if ( M_useRbSpaceContextEim && M_hasRbSpaceContextEim )
                boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) = M_model->computeBetaQm( uN[0], mu/*, N*/ );
            else
                boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) =
                    M_model->computeBetaQm( this->expansion( uN[0], N ), mu );
        }
#endif

        for(size_type q = 0;q < M_model->Ql(0); ++q)
        {
            for(int m=0; m < M_model->mMaxF(0,q); m++)
                Fdu += betaFqm[0][q][m]*M_Fqm_du[q][m].head(N);
        }
        for(size_type q = 0;q < M_model->Qa(); ++q)
        {
            for(int m=0; m < M_model->mMaxA(q); m++)
                Aprdu += betaAqm[q][m]*M_Aqm_pr_du[q][m].block(0,0,N,N);
        }

        du = uNdu[0];
        pr = uN[0];
        correction = -( Fdu.dot( du ) - du.dot( Aprdu*pr )  );
    }
    else
    {
        matrixN_type Mprdu ( (int)N, (int)N ) ;
        vectorN_type oldpr ( (int)N );
        beta_vector_type betaMqm;

        double time = 1e30;
        double dt = M_model->timeStep();
        double Tf = M_model->timeFinal();
        int K = Tf/dt;
        int primal_time_index;
        int dual_time_index;
        int kpd=0;//kp dual
        for( int kp=1; kp<=k; kp++)
        {

            Aprdu.setZero( N , N );
            Mprdu.setZero( N , N );
            Fdu.setZero( N );

            kpd=K-kp+1;
            dual_time_index = kpd;
            primal_time_index = kp;
            time = primal_time_index*dt;

            boost::tie( betaMqm, betaAqm, betaFqm) = M_model->computeBetaQm( mu ,time);

            //time_index--;


            for(size_type q = 0;q < M_model->Ql(0); ++q)
            {
                for(int m=0; m < M_model->mMaxF(0,q); m++)
                    Fdu += betaFqm[0][q][m]*M_Fqm_du[q][m].head(N);
            }


            for(size_type q = 0;q < M_model->Qa(); ++q)
            {
                for(int m=0;  m < M_model->mMaxA(q); m++)
                    Aprdu += betaAqm[q][m]*M_Aqm_pr_du[q][m].block(0,0,N,N);
            }


            for(size_type q = 0;q < M_model->Qm(); ++q)
            {
                for(int m=0; m<M_model->mMaxM(q); m++)
                    Mprdu += betaMqm[q][m]*M_Mqm_pr_du[q][m].block(0,0,N,N);
            }

            //du = uNdu[K-dual_time_index];
            du = uNdu[dual_time_index];
            pr = uN[primal_time_index];
            oldpr = uN[primal_time_index-1];//uNold[primal_time_index];
            correction += dt*( Fdu.dot( du ) - du.dot( Aprdu*pr ) ) - du.dot(Mprdu*pr) + du.dot(Mprdu*oldpr) ;
        }
    }

    return correction;

}


template<typename TruthModelType>
void
CRB<TruthModelType>::computeProjectionInitialGuess( const parameter_type & mu, int N , vectorN_type& initial_guess ) const
{
    VLOG(2) <<"Compute projection of initial guess\n";
    beta_vector_type betaMqm;
    beta_vector_type beta_initial_guess;

    matrixN_type Mass ( ( int )N, ( int )N ) ;
    vectorN_type F ( ( int )N );

    //beta coefficients of the initial guess ( mu-dependant part)
    beta_initial_guess = M_model->computeBetaInitialGuess( mu );

    //in steady case
    //for the mass matrix, the beta coefficient is 1
    //and the mu-indepenant part is assembled in crbmodel
    //WARNING : for unsteady case, don't call computeBetaQm to have beta coefficients
    //of the mass matrix, because computeBetaQm will compute all beta coefficients ( M A and F )
    //and beta coeff need to solve a model if we don't give an approximation of the unknown
    Mass.setZero( N,N );
    int q_max = M_Mqm_pr.size();
    for ( size_type q = 0; q < q_max; ++q )
    {
        int m_max = M_Mqm_pr[q].size();
        for(int m=0; m<m_max; m++)
        {
            Mass += 1*M_Mqm_pr[q][m].block( 0,0,N,N );
        }
    }

    F.setZero( N );
    q_max = M_InitialGuessV_pr.size();
    for ( size_type q = 0; q < q_max; ++q )
    {
        int m_max = M_InitialGuessV_pr[q].size();
        for(int m=0; m<m_max; m++)
        {
            F += beta_initial_guess[q][m]*M_InitialGuessV_pr[q][m].head( N );
        }
    }

    initial_guess = Mass.lu().solve( F );

#if 0
    export_ptrtype exporter;
    exporter = export_ptrtype( Exporter<mesh_type>::New( "FE_initial_guess" ) );

    auto FE_initial_guess_after = this->expansion( initial_guess , N );
    //auto min = FE_initial_guess_after.min();
    auto max = FE_initial_guess_after.max();

    if( this->worldComm().globalRank() == 0 )
        std::cout << "FE initial guess maxT = " << max << std::endl;

    exporter->step( 0 )->setMesh( FE_initial_guess_after.functionSpace()->mesh() );
    exporter->step( 0 )->add( "InitialGuess",  FE_initial_guess_after );
    exporter->save();
#endif

}


template<typename TruthModelType>
void
CRB<TruthModelType>::updateJacobian( const map_dense_vector_type& map_X, map_dense_matrix_type& map_J , const parameter_type& mu , int N) const
{
    //map_J.setZero( N , N );
    map_J.setZero( );
    beta_vector_type betaJqm;
    typedef std::vector< std::vector<sparse_matrix_ptrtype> > Jqm_type;
    typedef std::vector < std::vector<matrixN_type> > Jqm_pr_type;
    bool up = false;

    if( M_loadElementsDb )
    {
        boost::tie( boost::tuples::ignore, betaJqm, boost::tuples::ignore ) = M_model->computeBetaQm( this->expansion( map_X , N  ), mu , 0 );
        up = M_model->updateJacobian( this->expansion( map_X , N  ), const_cast<Jqm_type&>(this->M_Jqm) );
    }
    else
        boost::tie( boost::tuples::ignore,betaJqm, boost::tuples::ignore ) = M_model->computeBetaQm( mu , 0 );

    if( up ) // If affine decomposition of J needed bo te updated
    {
        for  (size_type q = 0; q < M_model->Qa(); ++q )
        {
            for( size_type m = 0; m < M_model->mMaxA(q); ++m )
            {
                for ( size_type i = 0; i < N; i++ )
                {
                    for ( size_type j = 0; j < N; ++j )
                    {
                        const_cast<Jqm_pr_type&>(M_Jqm_pr)[q][m]( i, j ) = M_Jqm[q][m]->energy( M_model->rBFunctionSpace()->primalBasisElement(i),
                                                                                                M_model->rBFunctionSpace()->primalBasisElement(j) );
                    }
                }
            }
        }
    }

    for ( size_type q = 0; q < M_model->Qa(); ++q )
    {
        for(int m=0; m<M_model->mMaxA(q); m++)
            map_J += betaJqm[q][m]*M_Jqm_pr[q][m].block( 0,0,N,N );
    }

}

template<typename TruthModelType>
void
CRB<TruthModelType>::updateResidual( const map_dense_vector_type& map_X, map_dense_vector_type& map_R , const parameter_type& mu, int N ) const
{
    map_R.setZero( );
    std::vector<beta_vector_type> betaRqm;

    typedef std::vector< std::vector<std::vector<vector_ptrtype> > > Rqm_type;
    typedef std::vector < std::vector<vectorN_type> > Rqm_pr_type;

    bool up = false;

    if( M_loadElementsDb )
    {
        boost::tie( boost::tuples::ignore, boost::tuples::ignore, betaRqm ) = M_model->computeBetaQm( this->expansion( map_X , N ), mu , 0 );
        up = M_model->updateResidual( this->expansion( map_X , N  ), const_cast<Rqm_type&>(this->M_Rqm) );
    }
    else
        boost::tie( boost::tuples::ignore, boost::tuples::ignore, betaRqm ) = M_model->computeBetaQm( mu , 0 );

    if( up ) // If affine decomposition of R needed bo te updated
    {
        for ( size_type q = 0; q < M_model->Ql( 0 ); ++q )
        {
            for( size_type m = 0; m < M_model->mMaxF( 0, q ); ++m )
            {
                for( int basis = 0; basis<N; basis++ )
                    const_cast<Rqm_pr_type&>(M_Rqm_pr)[q][m](basis) = inner_product( *M_Rqm[0][q][m],
                                                                                     M_model->rBFunctionSpace()->primalBasisElement(basis) );

            }
        }
    }

    for ( size_type q = 0; q < M_model->Ql( 0 ); ++q )
    {
        for(int m=0; m<M_model->mMaxF(0,q); m++)
        {
            map_R += betaRqm[0][q][m]*M_Rqm_pr[q][m].head( N );
            double norm=M_Rqm_pr[q][m].norm();
        }
    }
}

template<typename TruthModelType>
typename CRB<TruthModelType>::matrix_info_tuple
CRB<TruthModelType>::newton(  size_type N, parameter_type const& mu , vectorN_type & uN, double& output) const
{

    matrixN_type J ( ( int )N, ( int )N ) ;
    vectorN_type R ( ( int )N );

    double *r_data = R.data();
    double *j_data = J.data();
    double *uN_data = uN.data();

    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic , 1> > map_R ( r_data, N );
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic , 1> > map_uN ( uN_data, N );
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic , Eigen::Dynamic> > map_J ( j_data, N , N );

    typedef std::vector< std::vector<std::vector<vector_ptrtype> > > Rqm_type;
    typedef std::vector < std::vector<vectorN_type> > Lqm_type;
    bool up = false;

    computeProjectionInitialGuess( mu , N , uN );

    M_nlsolver->map_dense_jacobian = std::bind( &self_type::updateJacobian, std::ref( *this ), std::placeholders::_1, std::placeholders::_2  , mu , N );
    M_nlsolver->map_dense_residual = std::bind( &self_type::updateResidual, std::ref( *this ), std::placeholders::_1, std::placeholders::_2  , mu , N );
    M_nlsolver->setType( TRUST_REGION );
    M_nlsolver->solve( map_J , map_uN , map_R, 1e-12, 100);

    double conditioning=0;
    double determinant=0;
    if( M_computeMatrixInfo )
    {
        conditioning = computeConditioning( J );
        determinant = J.determinant();
    }

    std::vector<beta_vector_type> betaRqm;
    if( M_loadElementsDb )
    {
        boost::tie( boost::tuples::ignore, boost::tuples::ignore, betaRqm ) = M_model->computeBetaQm( this->expansion( uN , N  ), mu , 0 );
        up = M_model->updateResidual( this->expansion( uN , N ), const_cast<Rqm_type&>(this->M_Rqm) );
    }
    else
        boost::tie( boost::tuples::ignore, boost::tuples::ignore, betaRqm ) = M_model->computeBetaQm( mu , 0 );

    //update Lqm_pr if needed
    if( up )
    {
        for ( size_type q = 0; q < M_model->Ql( M_output_index ); ++q )
        {
            for( size_type m = 0; m < M_model->mMaxF( M_output_index, q ); ++m )
            {
                for( int basis = 0; basis<N; basis++ )
                {
                    const_cast<Lqm_type&>(M_Lqm_pr)[q][m]( basis ) = inner_product( *M_Rqm[M_output_index][q][m],
                                                                                    M_model->rBFunctionSpace()->primalBasisElement(basis) );
                    const_cast<Lqm_type&>(M_Lqm_du)[q][m]( basis ) = inner_product( *M_Rqm[M_output_index][q][m],
                                                                                    M_model->rBFunctionSpace()->dualBasisElement(basis) );

                }

            }//loop over m
        }//loop over q
    }

    //compute output
    vectorN_type L ( ( int )N );
    L.setZero( N );
    for ( size_type q = 0; q < M_model->Ql( M_output_index ); ++q )
    {
        for(int m=0; m < M_model->mMaxF(M_output_index,q); m++)
            L += betaRqm[M_output_index][q][m]*M_Lqm_pr[q][m].head( N );
    }

    output = L.dot( uN );

    auto matrix_info = boost::make_tuple(conditioning,determinant);
    return matrix_info;

}

template<typename TruthModelType>
double
CRB<TruthModelType>::computeConditioning( matrixN_type & A ) const
{
    Eigen::SelfAdjointEigenSolver< matrixN_type > eigen_solver;
    eigen_solver.compute( A );
    int number_of_eigenvalues =  eigen_solver.eigenvalues().size();
    //we copy eigenvalues in a std::vector beacause it's easier to manipulate it
    std::vector<double> eigen_values( number_of_eigenvalues );


    for ( int i=0; i<number_of_eigenvalues; i++ )
    {
        if ( imag( eigen_solver.eigenvalues()[i] )>1e-12 )
        {
            throw std::logic_error( "[CRB::lb] ERROR : complex eigenvalues were found" );
        }

        eigen_values[i]=real( eigen_solver.eigenvalues()[i] );
    }

    //LOG( INFO ) << " °°°°°°°°°°°°°°°°° EIGENVALUES ";
    //for(int i=0; i<number_of_eigenvalues; i++)
    //    LOG( INFO ) << "eig "<<i<<" : "<<eigen_values[i];
    int position_of_largest_eigenvalue=number_of_eigenvalues-1;
    int position_of_smallest_eigenvalue=0;
    double eig_max = eigen_values[position_of_largest_eigenvalue];
    double eig_min = eigen_values[position_of_smallest_eigenvalue];

    return eig_max / eig_min;
}



template<typename TruthModelType>
void
CRB<TruthModelType>::fixedPointDual(  size_type N, parameter_type const& mu, std::vector< vectorN_type > const& uN,
                                      std::vector< vectorN_type > & uNdu,  std::vector<vectorN_type> & uNduold, std::vector< double > & output_vector, int K ) const
{
    beta_vector_type betaAqm;
    beta_vector_type betaMqm;
    beta_vector_type betaMFqm;
    std::vector<beta_vector_type> betaFqm, betaLqm;

    matrixN_type Adu ( ( int )N, ( int )N ) ;
    matrixN_type Mdu ( ( int )N, ( int )N ) ;
    vectorN_type Fdu ( ( int )N );
    vectorN_type Ldu ( ( int )N );

    matrixN_type Aprdu( ( int )N, ( int )N );
    matrixN_type Mprdu( ( int )N, ( int )N );

    int Qa=M_model->Qa();
    int Ql=M_model->Ql(M_output_index);
    std::vector<int> mMaxA(Qa);
    std::vector<int> mMaxF( Ql );
    for ( size_type q = 0; q < Qa; ++q )
    {
        mMaxA[q]=M_model->mMaxA(q);
    }
    for ( size_type q = 0; q < Ql; ++q )
    {
        mMaxF[q]=M_model->mMaxF(M_output_index,q);
    }

    double increment = M_fixedpointIncrementTol;
    bool is_linear = M_model->isLinear();

    if ( M_model->isSteady() )
    {
        if ( is_linear )
        {
            boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) = M_model->computeBetaQm( mu );
            Adu.setZero( N,N );
            for ( size_type q = 0; q < Qa; ++q )
            {
                for(int m=0; m < mMaxA[q]; m++)
                    Adu += betaAqm[q][m]*M_Aqm_du[q][m].block( 0,0,N,N );
            }
            Fdu.setZero( N );
            for ( size_type q = 0; q < Ql ; ++q )
            {
                for(int m=0; m < mMaxF[q]; m++)
                    Fdu -= betaFqm[M_output_index][q][m]*M_Lqm_du[q][m].head( N );
            }
            uNdu[0] = Adu.lu().solve( Fdu );
        }
        else
        {
            vectorN_type next_uNdu( M_N );
            uNdu[0].setZero( N );
            int fi=0;
            do
            {
                // backup uNdu
                next_uNdu = uNdu[0];
                // update coefficients of affine decomposition
                if ( M_useRbSpaceContextEim && M_hasRbSpaceContextEim )
                    boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) = M_model->computeBetaQm( uN[0], mu/*, N*/ );
                else
                    boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) =
                        M_model->computeBetaQm( this->expansion( uN[0], N )/*dualRB*/, mu );
                // assemble rb matrix
                Adu.setZero( N,N );
                for ( size_type q = 0; q < Qa; ++q )
                {
                    for(int m=0; m < mMaxA[q]; m++)
                        Adu += betaAqm[q][m]*M_Aqm_du[q][m].block( 0,0,N,N );
                }
                // assemble rb rhs
                Fdu.setZero( N );
                for ( size_type q = 0; q < Ql ; ++q )
                {
                    for(int m=0; m < mMaxF[q]; m++)
                        Fdu -= betaFqm[M_output_index][q][m]*M_Lqm_du[q][m].head( N );
                }
                // solve rb system
                uNdu[0] = Adu.lu().solve( Fdu );

                increment = (uNdu[0]-next_uNdu).norm();
                if( M_fixedpointVerbose  && this->worldComm().isMasterRank() )
                {
                    VLOG(2)<<"[CRB::fixedPointDual] fixedpoint iteration " << fi << " increment error: " << increment;
                    std::cout<<"[CRB::fixedPointDual] fixedpoint iteration " << fi << " increment error: " << increment << "\n";
                }
                fi++;

            } while ( increment > M_fixedpointIncrementTol && fi<M_fixedpointMaxIterations );

        }
        return;
    }


    double time_for_output = 1e30;
    int number_of_time_step = M_model->numberOfTimeStep();
    double time_step = M_model->timeStep();
    //time_final = M_model->timeFinal();
    if ( K > 0 )
    {
        time_for_output = K * time_step;
        number_of_time_step = K;
    }
    else
    {
        //number_of_time_step = (time_final / time_step)+1;
        time_for_output = (number_of_time_step-1) * time_step;
    }
    int time_index = number_of_time_step-1;
    double time = time_for_output;

    size_type Qm = M_model->Qm();
    std::vector<int> mMaxM(Qm);
    for ( size_type q = 0; q < Qm; ++q )
    {
        mMaxM[q]=M_model->mMaxM(q);
    }


    for ( size_type n=0; n<N; n++ )
    {
        uNduold[time_index]( n ) = M_coeff_du_ini_online(n);
    }

    //uNdu[0] = Adu.lu().solve( -Ldu );
    //Adu.setZero( N,N );

    int time_iter=0;
    double tini = M_model->timeInitial();
    for ( time=time_for_output; math::abs(time - tini) > 1e-9; time-=time_step )
    {
        time_iter++;
        int fi=0;
        vectorN_type next_uNdu( M_N );

        do
        {
            if( time_iter == 1 )
            {
                bool only_terms_time_dependent=false;
                boost::tie( betaMqm, betaAqm, betaFqm ) = M_model->computeBetaQm( mu ,time , only_terms_time_dependent );

                Adu.setZero( N,N );
                for ( size_type q = 0; q < Qa; ++q )
                {
                    for(int m=0; m < mMaxA[q]; m++)
                        Adu += betaAqm[q][m]*M_Aqm_du[q][m].block( 0,0,N,N );
                }
                for ( size_type q = 0; q < Qm; ++q )
                {
                    for(int m=0; m < mMaxM[q]; m++)
                        Adu += betaMqm[q][m]*M_Mqm_du[q][m].block( 0,0,N,N )/time_step;
                }
            }
            else
            {
                bool only_terms_time_dependent=true;
                boost::tie( boost::tuples::ignore, boost::tuples::ignore, betaFqm ) = M_model->computeBetaQm( mu ,time , only_terms_time_dependent );
            }

            Fdu.setZero( N );
            //Ldu.setZero( N );
            //No Rhs for adjoint problem except mass contribution
            for ( size_type q = 0; q < Qm; ++q )
            {
                for(int m=0; m < mMaxM[q]; m++)
                {
                    Fdu += betaMqm[q][m]*M_Mqm_du[q][m].block( 0,0,N,N )*uNduold[time_index]/time_step;
                }
            }

            // backup uNdu
            next_uNdu = uNdu[time_index];

            uNdu[time_index] = Adu.lu().solve( Fdu );

            fi++;

            if( is_linear )
                break;

            //next_uNdu=uNdu[time_index];

            increment = (uNdu[time_index]-next_uNdu).norm();

            if( M_fixedpointVerbose  && this->worldComm().isMasterRank() )
            {
                VLOG(2)<<"[CRB::fixedPointDual] fixedpoint iteration " << fi << " increment error: " << increment;
                std::cout<<"[CRB::fixedPointDual] fixedpoint iteration " << fi << " increment error: " << increment << "\n";
            }

        }while ( increment > M_fixedpointIncrementTol && fi<M_fixedpointMaxIterations );

        if( increment > M_fixedpointIncrementTol )
            DVLOG(2)<<"[CRB::fixedPointDual] fixed point dual, proc "<<this->worldComm().globalRank()
                    <<" fixed point has no converged : norm of increment = "<<increment
                    <<" and tolerance : "<<M_fixedpointIncrementTol<<" so "<<M_fixedpointMaxIterations<<" iterations were done";

        if( increment > M_fixedpointCriticalValue )
            throw std::logic_error( "[CRB::fixedPointDual] fixed point ERROR : increment > critical value " );

        if( time_index > 0 )
            uNduold[time_index-1] = uNdu[time_index];

        time_index--;

    }//end of non steady case

#if 0
    double initial_dual_time = time_for_output+time_step;
    //std::cout<<"initial_dual_time = "<<initial_dual_time<<std::endl;
    boost::tie( betaMqm, betaAqm, betaFqm, betaMFqm ) = M_model->computeBetaQm( mu ,initial_dual_time );
    Mdu.setZero( N,N );

    for ( size_type q = 0; q < M_model->Qm(); ++q )
    {
        for ( size_type m = 0; m < M_model->mMaxM( q ); ++m )
        {
            for(int m=0; m < M_model->mMaxM(q); m++)
                Mdu += betaMqm[q][m]*M_Mqm_du[q][m].block( 0,0,N,N );
        }
    }

    Ldu.setZero( N );

    for ( size_type q = 0; q < M_model->Ql( M_output_index ); ++q )
    {
        for ( size_type m = 0; m < M_model->mMaxF(  M_output_index , q ); ++m )
        {
            for(int m=0; m < M_model->mMaxF(M_output_index,q); m++)
                Ldu += betaFqm[M_output_index][q][m]*M_Lqm_du[q][m].head( N );
        }
    }

    /*
     vectorN_type coeff(N);
     for(int i=0; i<N; i++) coeff(i) = M_coeff_du_ini_online[i];
     vectorN_type diff2 = uNduold[time_index] - coeff;
     std::cout<<"et maintenant le deuxieme diff = \n"<<diff2<<"\n";
     */
#endif

}

template<typename TruthModelType>
typename CRB<TruthModelType>::matrix_info_tuple
CRB<TruthModelType>::fixedPointPrimal(  size_type N, parameter_type const& mu, std::vector< vectorN_type > & uN,  std::vector<vectorN_type> & /*uNold*/,
                                        std::vector< double > & output_vector , int K, bool print_rb_matrix, bool computeOutput ) const
{
    //size_type Qm = 0;
    double output = 0;


    beta_vector_type betaAqm;
    beta_vector_type betaMqm;
    beta_vector_type betaMFqm;
    beta_vector_type beta_initial_guess;

    std::vector<beta_vector_type> betaFqm, betaLqm;

    int Qa=M_model->Qa();
    int Ql=M_model->Ql(M_output_index);
    int Qf=M_model->Ql(0);
    int Qm=M_model->Qm();
    std::vector<int> mMaxA(Qa);
    std::vector<int> mMaxM(Qm);
    std::vector<int> mMaxL( Ql );
    std::vector<int> mMaxF( Qf );
    for ( size_type q = 0; q < Qa; ++q )
    {
        mMaxA[q]=M_model->mMaxA(q);
    }
    for ( size_type q = 0; q < Qm; ++q )
    {
        mMaxM[q]=M_model->mMaxM(q);
    }
    for ( size_type q = 0; q < Qf; ++q )
    {
        mMaxF[q]=M_model->mMaxF(0,q);
    }
    for ( size_type q = 0; q < Ql; ++q )
    {
        mMaxL[q]=M_model->mMaxF(M_output_index,q);
    }

    matrixN_type A ( ( int )N, ( int )N ) ;
    vectorN_type F ( ( int )N );
    vectorN_type L ( ( int )N );

    double increment = M_fixedpointIncrementTol;
    bool is_linear = M_model->isLinear();

    if ( M_model->isSteady() )
    {
        if ( is_linear )
        {
            boost::tie( betaMqm, betaAqm, betaFqm ) = M_model->computeBetaQm( mu );
            A.setZero( N,N );
            for ( size_type q = 0; q < Qa; ++q )
            {
                for(int m=0; m<mMaxA[q]; m++)
                {
                    A += betaAqm[q][m]*M_Aqm_pr[q][m].block( 0,0,N,N );
                }
            }
            F.setZero( N );
            for ( size_type q = 0; q < Qf; ++q )
            {
                for(int m=0; m<mMaxF[q]; m++)
                {
                    F += betaFqm[0][q][m]*M_Fqm_pr[q][m].head( N );
                }
            }
            uN[0] = A.lu().solve( F );
        }
        else // nonlinear
        {
            vectorN_type previous_uN( N );
            //uN[0].setZero( N );
            computeProjectionInitialGuess( mu , N , uN[0] );
            int fi=0;

#if 1 // AITKEN
            bool useAitkenRelaxation = M_fixedpointUseAitken;
            vectorN_type rbaitkenAux( ( int )N );
            vectorN_type rbPreviousResidual( ( int )N );
            double rbAitkenTheta = 1.;//aitkenRelax.theta();
#endif
            bool fixPointIsFinished = false;
            do
            {
                // backup uNdu
                previous_uN = uN[0];
                // update coefficients of affine decomposition
                if ( M_useRbSpaceContextEim && M_hasRbSpaceContextEim )
                    boost::tie( betaMqm, betaAqm, betaFqm ) = M_model->computeBetaQm( uN[0], mu/*, N*/ );
                else
                    boost::tie( betaMqm, betaAqm, betaFqm ) =
                        M_model->computeBetaQm( this->expansion( uN[0], N ), mu );
                // assemble rb matrix
                A.setZero( N,N );
                for ( size_type q = 0; q < Qa; ++q )
                {
                    for(int m=0; m<mMaxA[q]; m++)
                    {
                        A += betaAqm[q][m]*M_Aqm_pr[q][m].block( 0,0,N,N );
                    }
                }
                // assemble rb rhs
                F.setZero( N );
                for ( size_type q = 0; q < Qf; ++q )
                {
                    for(int m=0; m<mMaxF[q]; m++)
                    {
                        F += betaFqm[0][q][m]*M_Fqm_pr[q][m].head( N );
                    }
                }
                // solve rb system
                if ( true )
                    uN[0] = A.fullPivLu().solve( F );
                else
                    uN[0] = A.lu().solve( F );

                if ( useAitkenRelaxation )
                {
                    auto rbresidual = uN[0]-previous_uN;

                    double residualConvergence = 1;
                    if ( fi >= 1 )
                    {
                        double oldEltL2Norm = previous_uN.norm();
                        if ( oldEltL2Norm > 1e-13 )
                            residualConvergence = rbresidual.norm()/oldEltL2Norm;
                        else
                            residualConvergence = rbresidual.norm();

                        if ( residualConvergence <  M_fixedpointIncrementTol || fi>=M_fixedpointMaxIterations )
                            fixPointIsFinished=true;
                        else
                        {
                            rbaitkenAux = rbresidual- rbPreviousResidual;
                            double scalar = rbaitkenAux.dot( M_algebraicInnerProductPrimal.block( 0,0,N,N )*rbaitkenAux );
                            rbaitkenAux *= ( 1.0/scalar );
                            scalar =  -rbAitkenTheta*rbPreviousResidual.dot(M_algebraicInnerProductPrimal.block( 0,0,N,N )*rbaitkenAux);
                            if ( scalar > 1 || scalar < 1e-4 )
                                scalar = 1.;

                            rbAitkenTheta = scalar;
                        }
                    }
                    if( M_fixedpointVerbose  && this->worldComm().isMasterRank() )
                        std::cout << "iteration=" << fi << " theta="<< rbAitkenTheta << " residualConv=" << residualConvergence <<"\n";

                    rbPreviousResidual = rbresidual;
                    // apply relaxation
                    uN[0] = previous_uN + rbAitkenTheta*rbresidual;
                }
                else
                {
                    double oldEltL2Norm = previous_uN.norm();
                    increment = (uN[0]-previous_uN).norm();
                    if ( oldEltL2Norm > 1e-13 )
                        increment /= oldEltL2Norm;

                    fixPointIsFinished = increment < M_fixedpointIncrementTol || fi>=M_fixedpointMaxIterations;

                    this->online_iterations_summary.first = fi;
                    this->online_iterations_summary.second = increment;

                    if( M_fixedpointVerbose  && this->worldComm().isMasterRank() )
                    {
                        DVLOG(2) << "iteration " << fi << " increment error: " << increment << "\n";
                        VLOG(2)<<"[CRB::fixedPointPrimal] fixedpoint iteration " << fi << " increment : " << increment <<std::endl;
                        double residual_norm = (A * uN[0] - F).norm() ;
                        VLOG(2) << " residual_norm :  "<<residual_norm;
                        std::cout << "[CRB::fixedPointPrimal] fixedpoint iteration " << fi << " increment : " << increment << " residual_norm :  "<<residual_norm <<std::endl;
                    }
                }
                ++fi;
            }
            while ( !fixPointIsFinished );// increment > M_fixedpointIncrementTol && fi<M_fixedpointMaxIterations );

        } // nonlinear

        if ( computeOutput )
        {
            L.setZero( N );
            for ( size_type q = 0; q < Ql; ++q )
            {
                for(int m=0; m < mMaxL[q]; m++)
                {
                    L += betaFqm[M_output_index][q][m]*M_Lqm_pr[q][m].head( N );
                }
            }
            output = L.dot( uN[0] );
            output_vector[0] = output;
        }

    }
    else // unsteady case
    {
        double time_for_output = 1e30;
        double time_step = 1e30;
        //double time_final = 1e30;
        int time_index = 0;
        int number_of_time_step = M_model->numberOfTimeStep();

        time_step = M_model->timeStep();
        //time_final = M_model->timeFinal();

        if ( K > 0 )
            time_for_output = K * time_step;
        else
        {
            //number_of_time_step = (time_final / time_step)+1;
            time_for_output = (number_of_time_step-1) * time_step;
        }

        for ( size_type n=0; n<N; n++ )
        {
            //uNold[1]( n ) = M_coeff_pr_ini_online(n);
            uN[0]( n ) = M_coeff_pr_ini_online(n);
        }

        if ( time_index<number_of_time_step-1 )
            time_index++;


        int time_iter=0;

        // if( !is_linear )
        // {
        //     computeProjectionInitialGuess( mu , N , uN[0] );
        // }
#if 1 // AITKEN
        bool useAitkenRelaxation = M_fixedpointUseAitken;
        vectorN_type rbaitkenAux( ( int )N );
        vectorN_type rbPreviousResidual( ( int )N );
        double rbAitkenTheta = 1.;//aitkenRelax.theta();
#endif
        //for ( double time=time_step; time<time_for_output+time_step; time+=time_step )
        for ( double time=time_step; math::abs(time - time_for_output - time_step) > 1e-9; time+=time_step )
        {

            time_iter++;
            //computeProjectionInitialGuess( mu , N , uN[time_index] );

            //vectorN_type error;
            //const element_type expansion_uN = this->expansion( uN[time_index] , N , M_WN);
            //checkInitialGuess( expansion_uN , mu , error);
            //std::cout<<"***************************************************************error.sum : "<<error.sum()<<std::endl;

            VLOG(2) << "lb: start fix point\n";

            vectorN_type previous_uN( M_N );
            int fi=0;

            bool fixPointIsFinished = false;
            do
            {
                if( is_linear )
                {
                    if( time_iter==1 )
                    {
                        bool only_terms_time_dependent=false;
                        boost::tie( betaMqm, betaAqm, betaFqm ) = M_model->computeBetaQm( mu ,time , only_terms_time_dependent );
                        A.setZero( N,N );
                        for ( size_type q = 0; q < Qa; ++q )
                        {
                            for(int m=0; m<mMaxA[q]; m++)
                            {
                                A += betaAqm[q][m]*M_Aqm_pr[q][m].block( 0,0,N,N );
                            }
                        }
                        for ( size_type q = 0; q < Qm; ++q )
                        {
                            for(int m=0; m<mMaxM[q]; m++)
                            {
                                A += betaMqm[q][m]*M_Mqm_pr[q][m].block( 0,0,N,N )/time_step;
                            }
                        }
                    }
                    else
                    {
                        bool only_terms_time_dependent=true;
                        boost::tie( boost::tuples::ignore, boost::tuples::ignore, betaFqm ) = M_model->computeBetaQm( mu ,time , only_terms_time_dependent );
                    }
                }
                else
                {
                    //important note :
                    //when lambda expressions will be totally operational
                    //we will call computeBetaQm( uN, mu, tim )
                    //and the test if( load_elements_db ) will disappear
                    if( M_loadElementsDb )
                    {
                        if( time_iter==1 )
                        {
                            bool only_terms_time_dependent=false;
                            boost::tie( betaMqm, betaAqm, betaFqm ) = M_model->computeBetaQm( this->expansion( uN[time_index], N ),
                                                                                              mu , time, only_terms_time_dependent );

                            A.setZero( N,N );
                            for ( size_type q = 0; q < Qa; ++q )
                            {
                                for(int m=0; m<mMaxA[q]; m++)
                                {
                                    A += betaAqm[q][m]*M_Aqm_pr[q][m].block( 0,0,N,N );
                                }
                            }
                            for ( size_type q = 0; q < Qm; ++q )
                            {
                                for(int m=0; m<mMaxM[q]; m++)
                                {
                                    A += betaMqm[q][m]*M_Mqm_pr[q][m].block( 0,0,N,N )/time_step;
                                }
                            }
                        }
                        else
                        {
                            bool only_terms_time_dependent=true;
                            boost::tie( boost::tuples::ignore, boost::tuples::ignore, betaFqm ) =
                                M_model->computeBetaQm( this->expansion( uN[time_index] , N  ),
                                                        mu ,time , only_terms_time_dependent );
                        }
                    }
                    else
                    {
                        bool only_terms_time_dependent=false;
                        boost::tie( betaMqm, betaAqm, betaFqm ) = M_model->computeBetaQm( mu ,time , only_terms_time_dependent );

                        A.setZero( N,N );
                        for ( size_type q = 0; q < Qa; ++q )
                        {
                            for(int m=0; m< mMaxA[q]; m++)
                            {
                                A += betaAqm[q][m]*M_Aqm_pr[q][m].block( 0,0,N,N );
                            }
                        }
                        for ( size_type q = 0; q < Qm; ++q )
                        {
                            for(int m=0; m< mMaxM[q]; m++)
                            {
                                A += betaMqm[q][m]*M_Mqm_pr[q][m].block( 0,0,N,N )/time_step;
                            }
                        }
                    }
                }

                F.setZero( N );
                for ( size_type q = 0; q < Qf; ++q )
                {
                    for(int m=0; m<mMaxF[q]; m++)
                    {
                        F += betaFqm[0][q][m]*M_Fqm_pr[q][m].head( N );
                    }
                }

                for ( size_type q = 0; q < Qm; ++q )
                {
                    for(int m=0; m<mMaxM[q]; m++)
                    {
                        //F += betaMqm[q][m]*M_Mqm_pr[q][m].block( 0,0,N,N )*uNold[time_index]/time_step;
                        F += betaMqm[q][m]*M_Mqm_pr[q][m].block( 0,0,N,N )*uN[time_index-1]/time_step;
                    }
                }

#if defined(FEELPP_HAS_HARTS) && defined(HARTS_HAS_OPENCL)
                if(ioption(_name="parallel.debug"))
                {
                    int nbelemAq = N * N * M_model->Qa();
                    double * arrayAq = new double[nbelemAq];
                    for( size_type q = 0; q < nbelemAq; ++q )
                        memcpy(arrayAq + q * N * N, M_Aqm_pr[q][0].block( 0,0,N,N ).data(), N * N * sizeof(double));
                    this->dumpData("./out.cpu.dump", "[CPU] Aq: ", arrayAq, nbelemAq);
                    delete[] arrayAq;

                    int nbelemFq = N * M_model->Ql(0);
                    double * arrayFq = new double[nbelemFq];
                    for( size_type q = 0; q < nbelemFq; ++q )
                        memcpy(arrayFq + q * N, M_Fqm_pr[q][0].head(N).data(), N * sizeof(double));
                    this->dumpData("./out.cpu.dump", "[CPU] Aq: ", arrayFq, nbelemFq);
                    delete[] arrayFq;

                    int nbelemBetaAq = M_model->Qa();
                    double * arrayBetaAq = new double[nbelemBetaAq];
                    for( size_type q = 0; q < nbelemBetaAq; ++q )
                        memcpy(arrayBetaAq + q, &(betaAqm[q][0]), sizeof(double));
                    this->dumpData("./out.gpu.dump", "[CPU] betaAq: ", arrayBetaAq, nbelemBetaAq);
                    delete[] arrayBetaAq;

                    int nbelemBetaFq = M_model->Ql(0);
                    double * arrayBetaFq = new double[nbelemBetaFq];
                    for( size_type q = 0; q < nbelemBetaFq; ++q )
                        memcpy(arrayBetaFq + q, &(betaFqm[0][q][0]), sizeof(double));
                    this->dumpData("./out.gpu.dump", "[CPU] betaFq: ", arrayBetaFq, nbelemBetaFq);
                    delete[] arrayBetaFq;

                    this->dumpData("./out.cpu.dump", "[CPU] A: ", A.data(), N * N);
                    this->dumpData("./out.cpu.dump", "[CPU] F: ", F.data(), N);
                }
#endif

                // backup uN
                previous_uN = uN[time_index];

                // solve for new fix point iteration
                uN[time_index] = A.lu().solve( F );

                //vectorN_type full_lu; full_lu.resize(2);
                //full_lu=A.fullPivLu().solve( F );
                //LOG( INFO ) << " oooooooooooooooooooooooooo mu = \n"<<mu;
                //LOG( INFO )<<std::setprecision(14)<<" norm of full LU : "<<full_lu.norm();
                //LOG( INFO )<<std::setprecision(14)<<" norm of  LU     : "<<uN[time_index].norm();

                //if ( time_index<number_of_time_step-1 )
                //    uNold[time_index+1] = uN[time_index];

                if( is_linear )
                    break;
                //previous_uN=uN[time_index];

                if ( useAitkenRelaxation )
                {
#if 1 // AITKEN
                    auto rbresidual = uN[time_index]-previous_uN;

                    double residualConvergence = 1;
                    if ( fi >= 1 )
                    {
                        double oldEltL2Norm = previous_uN.norm();
                        if ( oldEltL2Norm > 1e-13 )
                            residualConvergence = rbresidual.norm()/oldEltL2Norm;
                        else
                            residualConvergence = rbresidual.norm();

                        if ( residualConvergence <  M_fixedpointIncrementTol || fi>=M_fixedpointMaxIterations )
                            fixPointIsFinished=true;
                        else
                        {
                            rbaitkenAux = rbresidual- rbPreviousResidual;
                            double scalar = rbaitkenAux.dot( M_algebraicInnerProductPrimal.block( 0,0,N,N )*rbaitkenAux );
                            rbaitkenAux *= ( 1.0/scalar );
                            scalar =  -rbAitkenTheta*rbPreviousResidual.dot(M_algebraicInnerProductPrimal*rbaitkenAux);
                            if ( scalar > 1 || scalar < 1e-4 )
                                scalar = 1.;

                            rbAitkenTheta = scalar;
                        }
                    }
                    if( M_fixedpointVerbose  && this->worldComm().isMasterRank() )
                        std::cout << "iteration=" << fi << " theta="<< rbAitkenTheta << " residualConv=" << residualConvergence <<"\n";

                    rbPreviousResidual = rbresidual;
                    // apply relaxation
                    uN[time_index] = previous_uN + rbAitkenTheta*rbresidual;
#endif
                }
                else
                {
                    increment = (uN[time_index]-previous_uN).norm();
                    auto increment_abs = (uN[time_index]-previous_uN).array().abs();
                    fixPointIsFinished = increment < M_fixedpointIncrementTol || fi>=M_fixedpointMaxIterations;
                    //auto increment_abs2 = increment_abs1.array();
                    //auto increment_abs = increment_abs2.norm();

                    this->online_iterations_summary.first = fi;
                    this->online_iterations_summary.second = increment;

                    if( M_fixedpointVerbose  && this->worldComm().isMasterRank() )
                    {
                        DVLOG(2) << "iteration " << fi << " increment error: " << increment << "\n";
                        VLOG(2)<<"[CRB::fixedPointPrimal] fixedpoint iteration " << fi << " increment : " << increment <<std::endl;
                        double residual_norm = (A * uN[time_index] - F).norm() ;
                        VLOG(2) << " residual_norm :  "<<residual_norm;
                    }
                }
                ++fi;
            }
            while ( !fixPointIsFinished );// increment > M_fixedpointIncrementTol && fi<M_fixedpointMaxIterations );
            //while ( math::abs(output - old_output) >  output_fixedpoint_tol && fi < max_fixedpoint_iterations );

            if( increment > M_fixedpointIncrementTol )
                DVLOG(2)<<"[CRB::fixedPointPrimal] fixed point, proc "<<this->worldComm().globalRank()
                        <<" fixed point has no converged : increment = "<<increment
                        <<" and tolerance : "<<M_fixedpointIncrementTol<<" so "<<M_fixedpointMaxIterations<<" iterations were done"<<std::endl;

            if( increment > M_fixedpointCriticalValue )
                throw std::logic_error( "[CRB::fixedPointPrimal] fixed point ERROR : increment > critical value " );

            if ( computeOutput )
            {
                L.setZero( N );
                for ( size_type q = 0; q < Ql; ++q )
                {
                    for(int m=0; m < mMaxL[q]; m++)
                    {
                        L += betaFqm[M_output_index][q][m]*M_Lqm_pr[q][m].head( N );
                    }
                }
                output = L.dot( uN[time_index] );
                output_vector[time_index] = output;
            }

            if ( time_index<number_of_time_step-1 )
                time_index++;

        }
    } // unsteady case
    double condition_number = 0;
    double determinant = 0;
    if( M_computeMatrixInfo )
    {
        condition_number = computeConditioning( A );
        determinant = A.determinant();
    }

    auto matrix_info = boost::make_tuple(condition_number,determinant);

    if( print_rb_matrix && !M_offline_step )
        this->printRBMatrix( A,mu );

    return matrix_info;
}

template<typename TruthModelType>
void CRB<TruthModelType>::dumpData(std::string out, std::string prefix, const double * array, int nbelem) const
{
    std::ofstream ofs(out, std::ofstream::out | std::ofstream::app);

    if(nbelem)
    {
        ofs << prefix << array[0];
        for(int i = 1; i < nbelem; i++)
        {
            ofs << ";" << array[i];
        }
    }
    else
    {
        std::cout << ";;";
    }
    ofs << std::endl;

    ofs.close();
}

#if defined(FEELPP_HAS_HARTS) && defined(HARTS_HAS_OPENCL)
template<typename TruthModelType>
typename CRB<TruthModelType>::matrix_info_tuple
CRB<TruthModelType>::fixedPointPrimalCL(  size_type N, parameter_type const& mu, std::vector< vectorN_type > & uN,  std::vector<vectorN_type> & uNold,
                                          std::vector< double > & output_vector , int K, bool print_rb_matrix) const
{
    int i;
    int devID;
    int nDevOK;
    int nbelem;
    int nDevNeeded = 1;
    size_t devPWSM;
    size_t devLMS;

    double * dbuf;

    cl_device_fp_config fpConfig;
    cl_int err;
    cl_double dzero = 0.0;

    LOG( INFO ) << "[CRB::fixedPointPrimalCL] Checking for OpenCL support\n";

    /* initialize context and create queues for needed devices */
    if(soption(_name="parallel.opencl.device") == "gpu")
    {
        nDevOK = clContext_.init(CL_DEVICE_TYPE_GPU, nDevNeeded);
    }
    else
    {
        nDevOK = clContext_.init(CL_DEVICE_TYPE_CPU, nDevNeeded);
    }

    /* revert back to classical implementation */
    /* if no GPU is available */
    if(nDevOK != nDevNeeded)
    {
        LOG( INFO ) << "[CRB::fixedPointPrimalCL] Reverting to classic implementation\n";
        return fixedPointPrimal(N, mu, uN, uNold, output_vector, K, print_rb_matrix);
    }

    /* get back the device ID */
    //devID = clContext_.getDeviceID();
    cl::Context * context = clContext_.getContext();
    cl::CommandQueue * queue = clContext_.getCommandQueue(0);
    cl::Device * device = clContext_.getDevice(0);

    // TODO Typechecking of matrices
    // TODO Check whether its better to pass the matrices row-major or column-major

    /* Get beta factors */
    /* beta_vector_type is vector<vector<double> > */
    int time_index=0;
    double time = 1e30;
    beta_vector_type betaAqm;
    beta_vector_type betaMqm;
    std::vector<beta_vector_type> betaFqm;

    bool is_linear=M_model->isLinear();
    if( is_linear )
        boost::tie( betaMqm, betaAqm, betaFqm ) = M_model->computeBetaQm( mu ,time );
    else
        boost::tie( betaMqm, betaAqm, betaFqm ) = M_model->computeBetaQm( this->expansion( uN[time_index] , N ), mu ,time );

    cl::Event event;

    /*
     std::cout << "Params: N=" << N << "; mu=" << mu << "; uN.size()=" << uN.size() << "; uNold.size()=" << uNold.size()
     << "; output_vector.size()=" << output_vector.size() << "; K=" << K << "; print_rb_matrix=" << print_rb_matrix << std::endl;
     std::cout << "M_model->Qa(): " << M_model->Qa() << std::endl;
     std::cout << "M_model->Ql(0): " << M_model->Ql(0) << std::endl;
     */

    //gpuList[devID].getInfo(CL_DEVICE_LOCAL_MEM_SIZE, &devLMS);
    //std::cout << "Local Mem Size: " << devLMS << std::endl;

    /*
     gpuList[devID].getInfo(CL_DEVICE_DOUBLE_FP_CONFIG, &fpConfig);
     std::cout << "Double support: "
     << (fpConfig >= (CL_FP_FMA | CL_FP_ROUND_TO_NEAREST | CL_FP_ROUND_TO_ZERO | CL_FP_ROUND_TO_INF | CL_FP_INF_NAN | CL_FP_DENORM) ? "OK" : "KO") << std::endl;
     */

    /* create buffers on the GPU */
    /* we add one more matrix to store results */
    nbelem = N * N * M_model->Qa();
    dbuf = new double[nbelem];
    for( size_type q = 0; q < M_model->Qa(); ++q )
    {
        memcpy(dbuf + q * N * N, M_Aqm_pr[q][0].block( 0,0,N,N ).data(), N * N * sizeof(double));
    }
#if defined(FEELPP_HAS_HARTS) && defined(HARTS_HAS_OPENCL)
    if(ioption(_name="parallel.debug"))
    {
        this->dumpData("./out.gpu.dump", "[CPU] Aq: ", dbuf, nbelem);
    }
#endif
    cl::Buffer * Aq = clContext_.getBuffer("Aq", CL_MEM_READ_WRITE, nbelem * sizeof(double), NULL, &err);
    OPENCL_CHECK_ERR(err, "Could not allocate buffer");
    OPENCL_CHECK_ERR(queue->enqueueWriteBuffer(*Aq, CL_TRUE,
                                               0,
                                               nbelem * sizeof(double),
                                               dbuf,
                                               NULL, NULL), "Failed to write to buffer");

    delete[] dbuf;


    nbelem = N * M_model->Ql(0);
    dbuf = new double[nbelem];
    for( size_type q = 0; q < M_model->Ql(0); ++q )
    {
        memcpy(dbuf + q * N, M_Fqm_pr[q][0].head(N).data(), N * sizeof(double));
    }
#if defined(FEELPP_HAS_HARTS) && defined(HARTS_HAS_OPENCL)
    if(ioption(_name="parallel.debug"))
    {
        this->dumpData("./out.gpu.dump", "[CPU] Fq: ", dbuf, nbelem);
    }
#endif
    /* we add one more vector to store results */
    cl::Buffer * Fq = clContext_.getBuffer("Fq", CL_MEM_READ_WRITE, nbelem * sizeof(double), NULL, &err);
    OPENCL_CHECK_ERR(err, "Could not allocate buffer");
    OPENCL_CHECK_ERR(queue->enqueueWriteBuffer(*Fq, CL_TRUE,
                                               0,
                                               nbelem * sizeof(double),
                                               dbuf,
                                               NULL, NULL), "Failed to write to buffer");
    delete[] dbuf;


    nbelem = M_model->Qa();
    dbuf = new double[nbelem];
    for( size_type q = 0; q < M_model->Qa(); ++q )
    {
        dbuf[q] = betaAqm[q][0];
    }
#if defined(FEELPP_HAS_HARTS) && defined(HARTS_HAS_OPENCL)
    if(ioption(_name="parallel.debug"))
    {
        this->dumpData("./out.gpu.dump", "[CPU] betaAq: ", dbuf, nbelem);
        //this->dumpData("./out.gpu.dump", "[CPU] betaAq: ", betaAqm[0].data(), M_model->Qa());
    }
#endif
    cl::Buffer * betaAq = clContext_.getBuffer("betaAq", CL_MEM_READ_WRITE, nbelem * sizeof(double), NULL, &err);
    OPENCL_CHECK_ERR(err, "Could not allocate buffer");
    queue->enqueueWriteBuffer(*betaAq, CL_TRUE,
                              0, nbelem * sizeof(double),
                              dbuf,
                              NULL, NULL);
    delete[] dbuf;


    nbelem = M_model->Ql(0);
    dbuf = new double[nbelem];
    for( size_type q = 0; q < M_model->Ql(0); ++q )
    {
        dbuf[q] = betaFqm[0][q][0];
    }
#if defined(FEELPP_HAS_HARTS) && defined(HARTS_HAS_OPENCL)
    if(ioption(_name="parallel.debug"))
    {
        this->dumpData("./out.gpu.dump", "[CPU] betaFq: ", dbuf, nbelem);
        //this->dumpData("./out.gpu.dump", "[CPU] betaFq: ", betaFqm[0][0].data(), M_model->Ql(0));
    }
#endif
    cl::Buffer * betaFq = clContext_.getBuffer("betaFq", CL_MEM_READ_WRITE, nbelem * sizeof(double), NULL, &err);
    OPENCL_CHECK_ERR(err, "Could not allocate buffer");
    queue->enqueueWriteBuffer(*betaFq, CL_TRUE,
                              0, nbelem * sizeof(double),
                              dbuf,
                              NULL, NULL);
    delete[] dbuf;

    dbuf = new double[N * N];
    for(int i = 0; i < N * N; i++)
    { dbuf[i] = 0.0; }

    cl::Buffer * A = clContext_.getBuffer("A", CL_MEM_READ_WRITE, N * N * sizeof(double), NULL, &err);
    OPENCL_CHECK_ERR(err, "Could not allocate buffer");
    //queue.enqueueFillBuffer<double>(A, dzero, 0, N * N * sizeof(double), NULL, NULL);
    queue->enqueueWriteBuffer(*A, CL_TRUE,
                              0, N * N * sizeof(double),
                              dbuf,
                              NULL, NULL);

    cl::Buffer * F = clContext_.getBuffer("F", CL_MEM_READ_WRITE, N * sizeof(double), NULL, &err);
    OPENCL_CHECK_ERR(err, "Could not allocate buffer");
    //queue.enqueueFillBuffer<double>(F, dzero, 0, N * sizeof(double), NULL, NULL);
    queue->enqueueWriteBuffer(*F, CL_TRUE,
                              0, N * sizeof(double),
                              dbuf,
                              NULL, NULL);
    delete[] dbuf;

    cl::Program * program = clContext_.getProgram("CRB");
    if(!program)
    {
#if 0
        cl::Program::Sources source(1, std::make_pair(crb_kernels, strlen(crb_kernels)+1));
#endif
        std::string clsrc = Environment::findFile("crb.cl");
        /* If we found the cl source file, */
        /* we build the cl code */
        if(clsrc != "")
        {
            std::ifstream file(clsrc);
            err = file.is_open() ? CL_SUCCESS : -1;
            OPENCL_CHECK_ERR(err, "Could not open .cl file");

            std::string prog(std::istreambuf_iterator<char>(file), (std::istreambuf_iterator<char>()));

            clContext_.addProgram("CRB", prog.c_str(), prog.length());
            program = clContext_.getProgram("CRB");
        }
        /* otherwise we revert to the standard implementation */
        else
        {
            LOG( INFO ) << "[CRB::fixedPointPrimalCL] Could not find .cl source file. Reverting to classic implementation\n";
            /* Clean all previously allocated data */
            clContext_.clean();
            return fixedPointPrimal(N, mu, uN, uNold, output_vector, K, print_rb_matrix);
        }
    }

    /* Scalar * Matrices */
#if defined(FEELPP_HAS_HARTS) && defined(HARTS_HAS_OPENCL)
    if(ioption(_name="parallel.debug"))
    {
        int nbelem = N * N * M_model->Qa();
        double * array = new double[nbelem];
        OPENCL_CHECK_ERR(queue->enqueueReadBuffer(*Aq, CL_TRUE, 0, nbelem * sizeof(double), array, NULL, NULL), "Cannor read back data");

        this->dumpData("./out.gpu.dump", "[GPU] Aq: ", array, nbelem);

        delete[] array;
    }
#endif

    int nM = M_model->Qa();
    int nV = M_model->Ql(0);
    int NN = N * N;
    cl::Kernel smk(*program, "SVProd");

    smk.getWorkGroupInfo(*device, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, &devPWSM);
    //std::cout << "Preferred work group size: " << devPWSM << std::endl;

    OPENCL_CHECK_ERR(smk.setArg(0, *betaAq), "Could not add argument: betaAq");
    OPENCL_CHECK_ERR(smk.setArg(1, *Aq), "Could not add argument: Aq");
    OPENCL_CHECK_ERR(smk.setArg(2, sizeof(int), (void *)(&NN)), "Could not add argument: NN");

    OPENCL_CHECK_ERR(queue->enqueueNDRangeKernel(
                         smk,
                         cl::NullRange,
                         cl::NDRange(nM * devPWSM),
                         cl::NDRange(devPWSM),
                         NULL,
                         NULL), "Could not launch kernel");
    //event.wait();

#if defined(FEELPP_HAS_HARTS) && defined(HARTS_HAS_OPENCL)
    if(ioption(_name="parallel.debug"))
    {
        int nbelem = N * N * M_model->Qa();
        double * array = new double[nbelem];
        err = queue->enqueueReadBuffer(*Aq, CL_TRUE, 0, nbelem * sizeof(double), array, NULL, NULL);

        this->dumpData("./out.gpu.dump", "[GPU] Aq: ", array, nbelem);

        delete[] array;
    }
#endif

    /* Matrix sum */
    cl::Kernel msum(*program, "VSum");

    OPENCL_CHECK_ERR(msum.setArg(0, *A), "Could not add argument: A");
    OPENCL_CHECK_ERR(msum.setArg(1, *Aq), "Could not add argument: Aq");
    OPENCL_CHECK_ERR(msum.setArg(2, sizeof(int), (void *)(&NN)), "Could not add argument: NN");
    OPENCL_CHECK_ERR(msum.setArg(3, sizeof(int), (void *)(&nM)), "Could not add argument: nM");

    err = queue->enqueueNDRangeKernel(
        msum,
        cl::NullRange,
        cl::NDRange(devPWSM),
        cl::NDRange(devPWSM),
        NULL,
        NULL);
    OPENCL_CHECK_ERR(err, "Could not launch kernel");
    //event.wait();

#if defined(FEELPP_HAS_HARTS) && defined(HARTS_HAS_OPENCL)
    if(ioption(_name="parallel.debug"))
    {
        int nbelem = N * N;
        double * array = new double[nbelem];
        err = queue->enqueueReadBuffer(*A, CL_TRUE, 0, nbelem * sizeof(double), array, NULL, NULL);

        this->dumpData("./out.gpu.dump", "[GPU] A: ", array, nbelem);

        delete[] array;
    }
#endif

    /* Scalar * Vector */
#if defined(FEELPP_HAS_HARTS) && defined(HARTS_HAS_OPENCL)
    if(ioption(_name="parallel.debug"))
    {
        int nbelem = N * M_model->Ql( 0 );
        double * array = new double[nbelem];
        err = queue->enqueueReadBuffer(*Fq, CL_TRUE, 0, nbelem * sizeof(double), array, NULL, NULL);

        this->dumpData("./out.gpu.dump", "[GPU] Fq: ", array, nbelem);

        delete[] array;
    }
#endif
    cl::Kernel svk(*program, "SVProd");

    OPENCL_CHECK_ERR(svk.setArg(0, *betaFq), "Could not add argument: betaFq");
    OPENCL_CHECK_ERR(svk.setArg(1, *Fq), "Could not add argument: Fq");
    OPENCL_CHECK_ERR(svk.setArg(2, sizeof(int), (void *)(&N)), "Could not add argument: N");

    err = queue->enqueueNDRangeKernel(
        svk,
        cl::NullRange,
        cl::NDRange(nV * devPWSM),
        cl::NDRange(devPWSM),
        NULL,
        NULL);
    OPENCL_CHECK_ERR(err, "Could not launch kernel");
    //event.wait();

#if defined(FEELPP_HAS_HARTS) && defined(HARTS_HAS_OPENCL)
    if(ioption(_name="parallel.debug"))
    {
        int nbelem = N * M_model->Ql( 0 );
        double * array = new double[nbelem];
        err = queue->enqueueReadBuffer(*Fq, CL_TRUE, 0, nbelem * sizeof(double), array, NULL, NULL);

        this->dumpData("./out.gpu.dump", "[GPU] Fq: ", array, nbelem);

        delete[] array;
    }
#endif

    /* Vector Sum */
    cl::Kernel vsum(*program, "VSum");

    OPENCL_CHECK_ERR(vsum.setArg(0, *F), "Could not add argument: F");
    OPENCL_CHECK_ERR(vsum.setArg(1, *Fq), "Could not add argument: Fq");
    OPENCL_CHECK_ERR(vsum.setArg(2, sizeof(int), (void *)(&N)), "Could not add argument: N");
    OPENCL_CHECK_ERR(vsum.setArg(3, sizeof(int), (void *)(&nV)), "Could not add argument: nV");

    err = queue->enqueueNDRangeKernel(
        vsum,
        cl::NullRange,
        cl::NDRange(devPWSM),
        cl::NDRange(devPWSM),
        NULL,
        NULL);
    OPENCL_CHECK_ERR(err, "Could not launch kernel");
    //event.wait();

#if defined(FEELPP_HAS_HARTS) && defined(HARTS_HAS_OPENCL)
    if(ioption(_name="parallel.debug"))
    {
        int nbelem = N;
        double * array = new double[nbelem];
        err = queue->enqueueReadBuffer(*F, CL_TRUE, 0, nbelem * sizeof(double), array, NULL, NULL);

        this->dumpData("./out.gpu.dump", "[GPU] F: ", array, nbelem);

        delete[] array;
    }
#endif

    /* setup ViennaCL with current context */
    if(!(clContext_.isLinkedToVCL()))
    {
        viennacl::ocl::setup_context(0, (*context)(), (*device)(), (*queue)());
        clContext_.setLinkedToVCL();
    }

    /* wrap existing data */
    viennacl::vector<double> vclF((*F)(), N);
    viennacl::matrix<double> vclA((*A)(), N, N);
    viennacl::vector<double> vcl_result;
    std::vector<double> cpures;

    vcl_result = viennacl::linalg::solve(vclA, vclF, viennacl::linalg::cg_tag());

    cpures.reserve(vcl_result.size());
    viennacl::copy(vcl_result, cpures);

    // uN is a vector of vectorN_type, aka eigen's VectorXd
    for(int i = 0; i < vcl_result.size(); i++)
    {
        uN[time_index][i] = vcl_result[i];
    }

    // backup uN
    //previous_uN = uN[time_index];

    // fixedpointprimal code to merge
#if 0
    // solve for new fix point iteration
    uN[time_index] = A.lu().solve( F );

    if ( time_index<number_of_time_step-1 )
        uNold[time_index+1] = uN[time_index];

    L.setZero( N );
    for ( size_type q = 0; q < M_model->Ql( M_output_index ); ++q )
    {
        for(int m=0; m < M_model->mMaxF(M_output_index,q); m++)
        {
            L += betaFqm[M_output_index][q][m]*M_Lqm_pr[q][m].head( N );
        }
    }
    old_output = output;
    output = L.dot( uN[time_index] );

    if( is_linear )
        previous_uN=uN[time_index];

    increment = (uN[time_index]-previous_uN).norm();

    //output_vector.push_back( output );
    output_vector[time_index] = output;
    DVLOG(2) << "iteration " << fi << " increment error: " << increment << "\n";
    fi++;

    if( fixedpoint_verbose  && this->worldComm().globalRank()==this->worldComm().masterRank() )
        VLOG(2)<<"[CRB::fixedPointPrimal] fixedpoint iteration " << fi << " increment : " << increment <<std::endl;

    double residual_norm = (A * uN[time_index] - F).norm() ;
    VLOG(2) << " residual_norm :  "<<residual_norm;
#endif

    double condition_number = 0;
    double determinant = 0;

    auto matrix_info = boost::make_tuple(condition_number,determinant);

    return matrix_info;
}
#endif

template<typename TruthModelType>
typename CRB<TruthModelType>::matrix_info_tuple
CRB<TruthModelType>::onlineSolve(  size_type N, parameter_type const& mu, std::vector< vectorN_type > & uN, std::vector< vectorN_type > & uNdu, std::vector<vectorN_type> & uNold, std::vector<vectorN_type> & uNduold, std::vector< double > & output_vector, int K, bool print_rb_matrix, bool computeOutput ) const
{
    if( M_use_newton )
        return newton( N , mu , uN[0], output_vector[0] );
    else
        return fixedPoint( N ,  mu , uN , uNdu , uNold , uNduold , output_vector , K , print_rb_matrix, computeOutput );
}


template<typename TruthModelType>
typename CRB<TruthModelType>::matrix_info_tuple
CRB<TruthModelType>::fixedPoint(  size_type N, parameter_type const& mu, std::vector< vectorN_type > & uN, std::vector< vectorN_type > & uNdu,
                                  std::vector<vectorN_type> & uNold, std::vector<vectorN_type> & uNduold,
                                  std::vector< double > & output_vector, int K, bool print_rb_matrix, bool computeOutput ) const
{
    matrix_info_tuple matrix_info;

#if defined(FEELPP_HAS_HARTS) && defined(HARTS_HAS_OPENCL)
    if(boption(_name="parallel.opencl.enable"))
    {
        matrix_info = fixedPointPrimalCL( N, mu , uN , uNold, output_vector, K , print_rb_matrix) ;

        if(ioption(_name="parallel.debug"))
        {
            this->dumpData("./out.gpu.dump", "[CPU] uN: ", uN[0].data(), N);
        }
    }
    else
    {
#endif

        matrix_info = fixedPointPrimal( N, mu , uN , uNold, output_vector, K , print_rb_matrix, computeOutput ) ;

#if defined(FEELPP_HAS_HARTS) && defined(HARTS_HAS_OPENCL)
        if(ioption(_name="parallel.debug"))
        {
            this->dumpData("./out.cpu.dump", "[CPU] uN: ", uN[0].data(), N);
        }
    }
#endif


    if( M_solve_dual_problem )
    {
        fixedPointDual( N, mu , uN, uNdu , uNduold , output_vector , K ) ;

        if ( computeOutput )
        {
            // apply correction terms
            if ( M_model->isSteady() )
            {
                output_vector[0]+=correctionTerms(mu, uN , uNdu, uNold, 0 );
            }
            else
            {
                int number_of_time_step = M_model->numberOfTimeStep();
                double time_step = M_model->timeStep();
                double time_for_output = (number_of_time_step-1)*time_step;
                if ( K > 0 )
                    time_for_output = K*time_step;

                int time_index=0;
                for ( double time=0; math::abs(time-time_for_output-time_step)>1e-9; time+=time_step )
                {
                    output_vector[time_index]+=correctionTerms(mu, uN , uNdu, uNold, time_index );
                    ++time_index;
                }
            }
        }
    }

    return matrix_info;
}

template<typename TruthModelType>
typename boost::tuple<std::vector<double>,typename CRB<TruthModelType>::matrix_info_tuple>
CRB<TruthModelType>::lb( size_type N, parameter_type const& mu, std::vector< vectorN_type > & uN, std::vector< vectorN_type > & uNdu,
                         std::vector<vectorN_type> & uNold, std::vector<vectorN_type> & uNduold, bool print_rb_matrix, int K,
                         bool computeOutput ) const
{
    VLOG(2) << fmt::format("CRB::lb() : Start lb functions with mu={}, N={}, print_rb_matrix={}, K={}, computeOutput={}",mu.toString(),N,print_rb_matrix,K,computeOutput);
    if ( N > M_N ) N = M_N;

    int number_of_time_step = M_model->numberOfTimeStep();
    uN.resize( number_of_time_step );
    uNdu.resize( number_of_time_step );
    uNold.resize( number_of_time_step );
    uNduold.resize( number_of_time_step );

    for( int index=0;index<number_of_time_step;++index)
    {
        uN[index].resize( N );
        uNdu[index].resize( N );
        uNold[index].resize( N );
        uNduold[index].resize( N );
    }
    //vector containing outputs from time=time_step until time=time_for_output
    std::vector<double> output_vector;
    output_vector.resize( number_of_time_step, 0. );

    double conditioning=0;
    double determinant=0;
    // init by 1, the model could provide better init
    uN[0].setOnes(N);

    auto matrix_info = onlineSolve( N ,  mu , uN , uNdu , uNold , uNduold , output_vector , K , print_rb_matrix, computeOutput );

    if ( M_compute_variance || M_save_output_behavior )
    {
        double time_step = M_model->timeStep();
        double time_final = M_model->timeFinal();
        //if K>0 then the time at which we want to evaluate output is defined by
        //time_for_output = K * time_step
        //else it's the default value and in this case we take final time
        double time_for_output = time_step;
        if ( !M_model->isSteady() )
        {
            time_for_output = (number_of_time_step-1)*time_step;
            if ( K > 0 )
                time_for_output = K*time_step;
        }

        if( M_compute_variance )
        {
            if( ! M_database_contains_variance_info )
                throw std::logic_error( "[CRB::offline] ERROR there are no information available in the DataBase for variance computing, please set option crb.save-information-for-variance=true and rebuild the database" );

            int nb_spaces = functionspace_type::nSpaces;

            int space=0;

            int time_index=0;
            for ( double time=time_step; time<=time_for_output; time+=time_step )
            {
                vectorN_type uNsquare = uN[time_index].array().pow(2);
                double first = uNsquare.dot( M_variance_matrix_phi[space].block(0,0,N,N).diagonal() );

                double second = 0;
                for(int k = 1; k <= N-1; ++k)
                {
                    for(int j = 1; j <= N-k; ++j)
                    {
                        second += 2 * uN[time_index](k-1) * uN[time_index](k+j-1) * M_variance_matrix_phi[space](k-1 , j+k-1) ;
                    }
                }
                output_vector[time_index] = first + second;

                ++time_index;
            }
        }

        if ( M_save_output_behavior && this->worldComm().isMasterRank() )
        {
            int time_index=0;
            std::ofstream file_output;
            std::string mu_str;

            for ( int i=0; i<mu.size(); i++ )
            {
                mu_str= mu_str +  fmt::format( "_{}" ,mu[i] ) ;
            }

            std::string name = "output-evolution" + mu_str;
            file_output.open( name.c_str(),std::ios::out );
            file_output<<"time \t outputRB\n";
            for ( double time=0; math::abs(time-time_for_output-time_step)>1e-9; time+=time_step )
            {
                file_output<<time<<"\t"<<output_vector[time_index]<<"\n";
                ++time_index;
            }
            file_output.close();
        }
    }
    VLOG(2) <<"CRB::lb() : end lb function";
    return boost::make_tuple( output_vector, matrix_info);

}

template<typename TruthModelType>
typename CRB<TruthModelType>::error_estimation_type
CRB<TruthModelType>::delta( size_type N,
                            parameter_type const& mu,
                            std::vector<vectorN_type> const& uN,
                            std::vector<vectorN_type> const& uNdu,
                            std::vector<vectorN_type> const& /*uNold*/,
                            std::vector<vectorN_type> const& uNduold,
                            int k ) const
{
    VLOG(2) << "CRB::delta : Start compute delta with mu="<<mu.toString()<<", N="<<N;
    std::vector< std::vector<double> > primal_residual_coeffs;
    std::vector< std::vector<double> > dual_residual_coeffs;
    std::vector<double> output_upper_bound;
    double delta_pr=0;
    double delta_du=0;

    if( M_SER_errorEstimation && M_SER_useGreedyInRb ) // If SER is used : use Riesz residual norm as error indicator
    {
        VLOG(2) << "CRB::delta : Use SER error estimation \n";
        output_upper_bound.resize(1);
        output_upper_bound[0] = computeRieszResidualNorm( mu, uN );
        return boost::make_tuple( output_upper_bound , primal_residual_coeffs, dual_residual_coeffs , delta_pr, delta_du);
    }

    if ( M_error_type == CRB_NO_RESIDUAL )
    {
        output_upper_bound.resize(1);
        output_upper_bound[0]=-1;
        return boost::make_tuple( output_upper_bound ,primal_residual_coeffs,dual_residual_coeffs,delta_pr,delta_du );
    }

    else if ( M_error_type == CRB_EMPIRICAL )
    {
        output_upper_bound.resize(1);
        output_upper_bound[0]= empiricalErrorEstimation ( N, mu , k );
        return boost::make_tuple( output_upper_bound , primal_residual_coeffs, dual_residual_coeffs , delta_pr, delta_du);
    }

    else
    {
        double Tf = M_model->timeFinal();
        double Ti = M_model->timeInitial();
        double dt = M_model->timeStep();

        int time_index=1;
        int shift=1;
        int restart_time_index=1;
        if( M_model->isSteady() )
        {
            time_index=0;
            shift=0;
            restart_time_index=0;
        }

        double primal_sum=0;
        double dual_sum=0;
        double primal_sum_eim=0;
        double dual_sum_eim=0;

        //vectors to store residual coefficients
        int K = Tf/dt;
        primal_residual_coeffs.resize( K );
        dual_residual_coeffs.resize( K );

        double alphaA=1,alphaM=1;

        if ( M_error_type == CRB_RESIDUAL_SCM )
        {
            double alphaA_up, lbti;
            boost::tie( alphaA, lbti ) = M_scmA->lb( mu );
            if( M_scmA->useScm() )
                boost::tie( alphaA_up, lbti ) = M_scmA->ub( mu );
            //LOG( INFO ) << "alphaA_lo = " << alphaA << " alphaA_hi = " << alphaA_up ;

            if ( ! M_model->isSteady() )
            {
                double alphaM_up, lbti;
                boost::tie( alphaM, lbti ) = M_scmM->lb( mu );
                if( M_scmM->useScm() )
                    boost::tie( alphaM_up, lbti ) = M_scmM->ub( mu );
                //LOG( INFO ) << "alphaM_lo = " << alphaM << " alphaM_hi = " << alphaM_up ;
            }
        }

        //index associated to the output time
        int global_time_index=0;
        output_upper_bound.resize(K+1);

        bool compute_error_for_each_time_step = M_computeApeeForEachTimeStep;

        bool model_has_eim_error = M_model->hasEimError();

        double start_time=0;
        if( ! compute_error_for_each_time_step )
        {
            start_time=Tf;
            global_time_index=K;
        }

        for(double output_time=start_time; math::abs(output_time-Tf-dt)>1e-9; output_time+=dt)
        {
            time_index=restart_time_index;
            if( M_useAccurateApee )
            {
                //in this case, we use a different way to compute primal_sum (i.e. square of dual norm of primal residual)
                for( double time=dt; math::abs(time-output_time-dt)>1e-9; time+=dt )
                {
                    primal_sum += computeOnlinePrimalApee( N , mu , uN[time_index], uN[time_index]/*uNold[time_index]*/, dt, time );
                }
            }
            else
            {
                for ( double time=dt; math::abs(time-output_time-dt)>1e-9; time+=dt )
                {
                    int time_old_index = ( M_model->isSteady() )? 0 : time_index-1;
                    auto pr = transientPrimalResidual( N, mu, uN[time_index], uN[time_old_index]/*uNold[time_index]*/, dt, time );
                    primal_sum += pr.template get<0>();
                    if( global_time_index==K )
                    {
                        primal_residual_coeffs[time_index-shift].resize( pr.template get<1>().size() );
                        primal_residual_coeffs[time_index-shift] = pr.template get<1>() ;
                    }

                    if( model_has_eim_error )
                    {
                        auto preim = transientPrimalResidualEim( N, mu, uN[time_index], uN[time_old_index]/*uNold[time_index]*/, dt, time );
                        primal_sum_eim += preim.template get<0>();
                    }
                    time_index++;
                }//end of time loop for primal problem
            }

            time_index--;

            double dual_residual=0;

            if( M_solve_dual_problem )
            {
                if ( !M_model->isSteady() )
                    dual_residual = initialDualResidual( N,mu,uNduold[time_index],dt );

                if( M_useAccurateApee )
                {
                    //in this case, we use a different way to compute primal_sum (i.e. square of dual norm of primal residual)
                    for ( double time=output_time; math::abs(time - dt + dt) > 1e-9 ; time-=dt )
                    {
                        dual_sum += computeOnlineDualApee( N , mu , uNdu[time_index], uNduold[time_index], dt, time );
                    }
                }
                else
                {
                    for ( double time=output_time; math::abs(time - dt + dt) > 1e-9 ; time-=dt )
                    {
                        auto du = transientDualResidual( N, mu, uNdu[time_index], uNduold[time_index], dt, time );
                        dual_sum += du.template get<0>();
                        if( global_time_index == K )
                        {
                            dual_residual_coeffs[time_index-shift].resize( du.template get<1>().size() );
                            dual_residual_coeffs[time_index-shift] = du.template get<1>();
                        }

                        if( model_has_eim_error )
                        {
                            int time_old_pr_index = ( M_model->isSteady() )? 0 : time_index-1;
                            auto dueim = transientDualResidualEim( N, mu, uN[time_index], uN[time_old_pr_index]/*uNold[time_index]*/, dt, time );
                            dual_sum_eim += dueim.template get<0>();
                        }

                        time_index--;
                    }//end of time loop for dual problem
                }//not with accurate apee
            }//solve dual problem

            double solution_upper_bound;
            double solution_dual_upper_bound;
            //alphaA=1;
            //dual_residual=0;
            if ( M_model->isSteady() )
            {
                if( model_has_eim_error )
                {
                    double r = math::sqrt( primal_sum );
                    double reim = math::sqrt( primal_sum_eim );
                    delta_pr =  ( r + reim ) /  alphaA ;
                    if( M_solve_dual_problem )
                    {
                        double rdu = math::sqrt( dual_sum );
                        double rdueim = math::sqrt( dual_sum_eim );
                        delta_du =  ( rdu + rdueim ) / alphaA;
                    }
                    else
                        delta_du = 1;
                    output_upper_bound[global_time_index] = alphaA * delta_pr * delta_du;
                }
                else
                {
                    delta_pr = math::sqrt( primal_sum ) /  alphaA ;
                    if( M_solve_dual_problem )
                        delta_du = math::sqrt( dual_sum ) / alphaA;
                    else
                        delta_du = 1;
                    output_upper_bound[global_time_index] = alphaA * delta_pr * delta_du;
                }
                //solution_upper_bound =  delta_pr;
                //solution_dual_upper_bound =  delta_du;
            }
            else
            {
                delta_pr = math::sqrt( dt/alphaA * primal_sum );
                if( M_solve_dual_problem )
                    delta_du = math::sqrt( dt/alphaA * dual_sum + dual_residual/alphaM );
                else
                    delta_du = 1;
                output_upper_bound[global_time_index] = delta_pr * delta_du;
                //solution_upper_bound = delta_pr;
                //solution_dual_upper_bound =  delta_du;
            }

            global_time_index++;
        }//end of loop over output time

        if( !M_offline_step && M_showResidual )
        {
            double sum=0;
            time_index=1;
            if( M_model->isSteady() )
                time_index=0;
            LOG( INFO ) <<" =========== Residual with "<<N<<" basis functions - seek mu in complement of WNmu : "<< M_seekMuInComplement<<"============ \n";
            for ( double time=dt; time<=Tf; time+=dt )
            {
                int time_old_index = ( M_model->isSteady() )? 0 : time_index-1;
                auto pr = transientPrimalResidual( N, mu, uN[time_index], uN[time_old_index]/*uNold[time_index]*/, dt, time );
                //LOG(INFO) << "primal residual at time "<<time<<" : "<<pr.template get<0>()<<"\n";
                sum+=pr.template get<0>();
                time_index++;
            }
            LOG(INFO) << "sum of primal residuals  "<<sum<<std::endl;

            time_index--;
            sum=0;

            if (M_solve_dual_problem )
            {
                for ( double time=Tf; time>=dt; time-=dt )
                {
                    auto du = transientDualResidual( N, mu, uNdu[time_index], uNduold[time_index], dt, time );
                    //LOG(INFO) << "dual residual at time "<<time<<" : "<<du.template get<0>()<<"\n";
                    sum += du.template get<0>();
                    time_index--;
                }
            }
            LOG(INFO) << "sum of dual residuals  "<<sum<<std::endl;
            LOG( INFO ) <<" ================================= \n";
            //std::cout<<"[REAL ] duam_sum : "<<sum<<std::endl;
        }//if show_residual_convergence

        //return boost::make_tuple( output_upper_bound, primal_residual_coeffs, dual_residual_coeffs , delta_pr, delta_du , solution_upper_bound, solution_dual_upper_bound);
        return boost::make_tuple( output_upper_bound, primal_residual_coeffs, dual_residual_coeffs , delta_pr, delta_du );

    }//end of else
}

template<typename TruthModelType>
typename CRB<TruthModelType>::max_error_type
CRB<TruthModelType>::maxErrorBounds( size_type N ) const
{
    int proc = this->worldComm().globalRank();
    int master_proc = this->worldComm().masterRank();

    std::vector< vectorN_type > uN;
    std::vector< vectorN_type > uNdu;
    std::vector< vectorN_type > uNold;
    std::vector< vectorN_type > uNduold;

    y_type err( M_Xi->size() );
    y_type vect_delta_pr( M_Xi->size() );
    y_type vect_delta_du( M_Xi->size() );
    std::vector<double> check_err( M_Xi->size() );

    double delta_pr=0;
    double delta_du=0;

    if ( M_error_type == CRB_EMPIRICAL )
    {
        if ( M_WNmu->size()==1 )
        {
            parameter_type mu ( M_Dmu );
            size_type id;
            boost::tie ( mu, id ) = M_Xi->max();

            return boost::make_tuple( 1e5, mu , delta_pr, delta_du);
        }

        else
        {
            err.resize( M_WNmu_complement->size() );
            check_err.resize( M_WNmu_complement->size() );

            for ( size_type k = 0; k < M_WNmu_complement->size(); ++k )
            {
                parameter_type const& mu = M_WNmu_complement->at( k );
                lb( N, mu, uN, uNdu , uNold ,uNduold );
                auto error_estimation = delta( N, mu, uN, uNdu, uNold, uNduold, k );
                auto vector_err = error_estimation.template get<0>();
                int size=vector_err.size();
                double _err = vector_err[size-1];
                vect_delta_pr( k ) = error_estimation.template get<3>();
                vect_delta_du( k ) = error_estimation.template get<4>();
                err( k ) = _err;
                check_err[k] = _err;
            }
        }

    }//end of if ( M_error_type == CRB_EMPIRICAL)

    else
    {
        if ( M_seekMuInComplement )
        {
            err.resize( M_WNmu_complement->size() );
            check_err.resize( M_WNmu_complement->size() );
            for ( size_type k = 0; k < M_WNmu_complement->size(); ++k )
            {
                parameter_type const& mu = M_WNmu_complement->at( k );
                lb( N, mu, uN, uNdu , uNold ,uNduold, false, 0, false );
                auto error_estimation = delta( N, mu, uN, uNdu, uNold, uNduold, k );
                auto const& vector_err = error_estimation.template get<0>();
                //int size=vector_err.size();
                //double _err = vector_err[size-1];
                double _err = vector_err.back();
                vect_delta_pr( k ) = error_estimation.template get<3>();
                vect_delta_du( k ) = error_estimation.template get<4>();
                err( k ) = _err;
                check_err[k] = _err;
            }

        }//end of seeking mu in the complement

        else
        {
            for ( size_type k = 0; k < M_Xi->size(); ++k )
            {
                //std::cout << "--------------------------------------------------\n";
                parameter_type const& mu = M_Xi->at( k );
                //std::cout << "[maxErrorBounds] mu=" << mu << "\n";
                lb( N, mu, uN, uNdu , uNold ,uNduold, false, 0, false );
                //auto tuple = lb( N, mu, uN, uNdu , uNold ,uNduold );
                //double o = tuple.template get<0>();
                //std::cout << "[maxErrorBounds] output=" << o << "\n";
                auto error_estimation = delta( N, mu, uN, uNdu, uNold, uNduold, k );
                auto vector_err = error_estimation.template get<0>();
                int size=vector_err.size();
                double _err = vector_err[size-1];
                vect_delta_pr( k ) = error_estimation.template get<3>();
                vect_delta_du( k ) = error_estimation.template get<4>();
                //std::cout << "[maxErrorBounds] error=" << _err << "\n";
                err( k ) = _err;
                check_err[k] = _err;
            }
        }//else ( seek_mu_in_complement )
    }//else

    Eigen::MatrixXf::Index index;
    double maxerr = err.array().abs().maxCoeff( &index );
    delta_pr = vect_delta_pr( index );
    delta_du = vect_delta_du( index );
    parameter_type mu;

#if 0
    std::vector<double>::iterator it = std::max_element( check_err.begin(), check_err.end() );
    int check_index = it - check_err.begin() ;
    double check_maxerr = *it;

    if ( index != check_index || maxerr != check_maxerr )
    {
        std::cout<<"[CRB::maxErrorBounds] index = "<<index<<" / check_index = "<<check_index<<"   and   maxerr = "<<maxerr<<" / "<<check_maxerr<<std::endl;
        throw std::logic_error( "[CRB::maxErrorBounds] index and check_index have different values" );
    }
#endif
    int _index=0;

    if ( M_seekMuInComplement )
    {
        LOG(INFO) << "[maxErrorBounds] WNmu_complement N=" << N << " max Error = " << maxerr << " at index = " << index << "\n";
        mu = M_WNmu_complement->at( index );
        _index = M_WNmu_complement->indexInSuperSampling( index );
    }
    else
    {
        LOG(INFO) << "[maxErrorBounds] N=" << N << " max Error = " << maxerr << " at index = " << index << "\n";
        mu = M_Xi->at( index );
        _index = index;
    }

    //do communications to have global max
    int world_size = this->worldComm().globalSize();
    std::vector<double> max_world( world_size );
    mpi::all_gather( this->worldComm().globalComm(),
                     maxerr,
                     max_world );
    auto it_max = std::max_element( max_world.begin() , max_world.end() );
    int proc_having_good_mu = it_max - max_world.begin();

    auto tuple = boost::make_tuple( mu , delta_pr , delta_du );

    boost::mpi::broadcast( this->worldComm() , tuple , proc_having_good_mu );

    mu = tuple.template get<0>();
    delta_pr = tuple.template get<1>();
    delta_du = tuple.template get<2>();
    mu.check();
    //now all proc have same maxerr
    maxerr = *it_max;

    // SER : criterion for RB r-adaptation
    if( M_SER_errorEstimation )
    {
        double increment = math::abs( maxerr - M_SER_maxerr );
        double inc_relative = increment/math::abs( M_SER_maxerr );
        // Increase size of next RB group if necessary
        if( increment > 1e-10 && inc_relative > 0 && inc_relative < doption(_prefix=M_prefix,_name="ser.radapt-rb-rtol") )
            this->setAdaptationSER( true );

        if( this->worldComm().isMasterRank() )
            std::cout << "[RB] SER adaptation : " << this->adaptationSER()  << std::endl;
        if( !adaptationSER() )
            M_SER_maxerr = maxerr;
    }

    if( this->worldComm().isMasterRank() )
        std::cout<< std::setprecision(15)<<"[CRB maxerror] proc "<< proc
                 <<" delta_pr : "<<delta_pr<<" -- delta_du : "<<delta_du
                 <<" -- output error : "<<maxerr<<std::endl;
    //lb( N, mu, uN, uNdu , uNold ,uNduold );

    if( proc == master_proc )
    {
        int size = mu.size();
        std::cout << "  [maxErrorBounds] -- max error bound (" << maxerr << ")"
                  << "  -- mu = [ ";
        for ( int i=0; i< size-1; i++ )
            std::cout<< mu( i ) <<" ";
        std::cout << mu( size-1 ) << " ]" <<std::endl;
    }

    return boost::make_tuple( maxerr, mu , delta_pr, delta_du);


}

template<typename TruthModelType>
double
CRB<TruthModelType>::orthonormalize( size_type N, wn_type& wn, int Nm )
{
    int proc_number = this->worldComm().globalRank();
    if( proc_number == 0 ) std::cout << "  -- orthonormalization (Gram-Schmidt)\n";
    DVLOG(2) << "[CRB::orthonormalize] orthonormalize basis for N = " << N << "\n";
    DVLOG(2) << "[CRB::orthonormalize] orthonormalize basis for WN = " << wn.size() << "\n";
    DVLOG(2) << "[CRB::orthonormalize] starting ...\n";

#if 0
    for ( size_type i = 0; i < N; ++i )
    {
        for( size_type j=0; j<i; ++j )
        {
            value_type __rij_pr = M_model->scalarProduct(  wn[i], wn[j] );
            wn[i].add( -__rij_pr, wn[j] );
        }
    }
#else
    for ( size_type i =N-Nm; i < N; ++i )
    {
        auto & wni = unwrap_ptr( wn[i] );
        for ( size_type j = 0; j < i; ++j )
        {
            auto const& wnj = unwrap_ptr( wn[j] );
            value_type __rij_pr = M_model->scalarProduct( wni, wnj );
            wni.add( -__rij_pr, wnj );
        }
    }
#endif
    // normalize
    for ( size_type i =N-Nm; i < N; ++i )
    {
        auto & wni = unwrap_ptr( wn[i] );
        value_type __rii_pr = math::sqrt( M_model->scalarProduct( wni, wni ) );
        wni.scale( 1./__rii_pr );
    }

    DVLOG(2) << "[CRB::orthonormalize] finished ...\n";
    DVLOG(2) << "[CRB::orthonormalize] copying back results in basis\n";

    return checkOrthonormality( N , wn );
}

template <typename TruthModelType>
double
CRB<TruthModelType>::checkOrthonormality ( int N, const wn_type& wn ) const
{

    if ( wn.size()==0 )
    {
        throw std::logic_error( "[CRB::checkOrthonormality] ERROR : size of wn is zero" );
    }

    if ( M_orthonormalize_primal*M_orthonormalize_dual==0 )
    {
        LOG(INFO)<<"Warning : calling checkOrthonormality is called but ";
        LOG(INFO)<<" orthonormalize_dual = "<<M_orthonormalize_dual;
        LOG(INFO)<<" and orthonormalize_primal = "<<M_orthonormalize_primal;
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
    if ( this->worldComm().isMasterRank() )
    {
        LOG( INFO ) << "    o check : " << A.norm() << " (should be 0)";
    }
    //FEELPP_ASSERT( A.norm() < 1e-14 )( A.norm() ).error( "orthonormalization failed.");

    return A.norm();
}


template <typename TruthModelType>
void
CRB<TruthModelType>::exportBasisFunctions( const export_vector_wn_type& export_vector_wn )const
{


    std::vector<wn_type> vect_wn=export_vector_wn.template get<0>();
    std::vector<std::string> vect_names=export_vector_wn.template get<1>();

    if ( vect_wn.size()==0 )
    {
        throw std::logic_error( "[CRB::exportBasisFunctions] ERROR : there are no wn_type to export" );
    }

    auto first_wn = vect_wn[0];
    auto first_element = first_wn[0];
    if ( !M_exporter )
        M_exporter = exporter( _mesh=unwrap_ptr(first_element).functionSpace()->mesh(), _name="BasisFunctions");

    int basis_number=0;
    for( auto const& wn : vect_wn )
    {

        if ( wn.size()==0 )
        {
            throw std::logic_error( "[CRB::exportBasisFunctions] ERROR : there are no element to export" );
        }

        int element_number=0;
        parameter_type mu;

        for( auto const& element : wn )
        {

            std::string basis_name = vect_names[basis_number];
            std::string number = fmt::format( "{}_with_parameters",element_number );
            mu = M_WNmu->at( element_number );
            std::string mu_str;

            for ( int i=0; i<mu.size(); i++ )
            {
                mu_str= mu_str + fmt::format( "_{}",mu[i] ) ;
            }

            std::string name =   basis_name + number + mu_str;
            M_exporter->step( 0 )->add( name, unwrap_ptr( element ) );
            element_number++;
        }
        basis_number++;
    }


    M_exporter->save();

}


//error estimation only to build reduced basis
template <typename TruthModelType>
typename CRB<TruthModelType>::value_type
CRB<TruthModelType>::empiricalErrorEstimation ( int Nwn, parameter_type const& mu , int second_index ) const
{
    std::vector<vectorN_type> Un;
    std::vector<vectorN_type> Undu;
    std::vector<vectorN_type> Unold;
    std::vector<vectorN_type> Unduold;
    auto o = lb( Nwn, mu, Un, Undu, Unold, Unduold  );
    auto output_vector=o.template get<0>();
    double output_vector_size=output_vector.size();
    double sn = output_vector[output_vector_size-1];//output at last time

    int nb_element =Nwn/M_factor*( M_factor>0 ) + ( Nwn+M_factor )*( M_factor<0 && ( int )Nwn>( -M_factor ) ) + 1*( M_factor<0 && ( int )Nwn<=( -M_factor ) )  ;

    std::vector<vectorN_type> Un2;
    std::vector<vectorN_type> Undu2;//( nb_element );

    //double output_smaller_basis = lb(nb_element, mu, Un2, Undu2, Unold, Unduold);
    auto tuple = lb( nb_element, mu, Un2, Undu2, Unold, Unduold );
    output_vector=tuple.template get<0>();
    output_vector_size=output_vector.size();
    double output_smaller_basis = output_vector[output_vector_size-1];

    double error_estimation = math::abs( sn-output_smaller_basis );

    return error_estimation;

}




template<typename TruthModelType>
typename CRB<TruthModelType>::value_type
CRB<TruthModelType>::initialDualResidual( int Ncur, parameter_type const& mu, vectorN_type const& Unduini, double time_step ) const
{
    beta_vector_type betaAqm;
    beta_vector_type betaMqm;
    beta_vector_type betaMFqm;
    std::vector<beta_vector_type> betaFqm;
    double time = M_model->timeFinal();
    boost::tie( betaMqm, betaAqm, betaFqm) = M_model->computeBetaQm( mu, time );

    int __QLhs = M_model->Qa();
    int __QOutput = M_model->Ql( M_output_index );
    int __Qm = M_model->Qm();
    int __N = Ncur;


    value_type __c0_du = 0.0;

    for ( int __q1 = 0; __q1 < __QOutput; ++__q1 )
    {
        for ( int __m1 = 0; __m1 < M_model->mMaxF(M_output_index,__q1); ++__m1 )
        {
            value_type fq1 = betaFqm[M_output_index][__q1][__m1];

            for ( int __q2 = 0; __q2 < __QOutput; ++__q2 )
            {
                for ( int __m2 = 0; __m2 < M_model->mMaxF(M_output_index,__q2); ++__m2 )
                {
                    value_type fq2 = betaFqm[M_output_index][__q2][__m2];
                    //__c0_du += 1./(time_step*time_step) * M_C0_du[__q1][__q2]*fq1*fq2;
                    //__c0_du += 1./(time_step*time_step) * M_C0_du[__q1][__q2]*fq1*fq2;
                    __c0_du +=  M_C0_du[__q1][__m1][__q2][__m2]*fq1*fq2;
                }//end of loop __m2
            }//end of loop __q2
        }//end of loop __m1
    }//end of loop __q1

#if 0
    value_type __Caf_du = 0.0;
    value_type __Caa_du = 0.0;

    for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
    {
        for ( int __m1 = 0; __m1 < M_model->mMaxA(__q1) ; ++__m1 )
        {

            value_type a_q1 = betaAqm[__q1][__m1];

            for ( int __q2 = 0; __q2 < __QOutput; ++__q2 )
            {

                for ( int __m2 = 0; __m2 < M_model->mMaxF(M_output_index,__q2); ++__m2 )
                {
                    value_type f_q2 = betaFqm[M_output_index][__q2][__m2];
                    __Caf_du += 1./time_step * a_q1*f_q2*M_Lambda_du[__q1][__m1][__q2][__m2].dot( Unduini );
                }
            }


            for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
            {
                for ( int __m2 = 0; __m2 < M_model->mMaxA(__q2); ++__m2 )
                {
                    value_type a_q2 = betaAqm[__q2][__m2];
                    auto m = M_Gamma_du[__q1][__m1][__q2][__m2].block( 0,0,__N,__N )*Unduini;
                    __Caa_du += a_q1 * a_q2 * Unduini.dot( m );
                }
            }
        }
    }


    for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
    {
        for ( int __m2=0 ; __m2< M_model->mMaxA(__q2); ++__m2 )
        {
            value_type a_q2 = betaAqm[__q2][__m2];
            auto m = M_Cma_du[__q1][__m1][__q2][__m2].block( 0,0,__N,__N )*Unduini;
            //__Cma_du += 1./time_step * m_q1 * a_q2 * Unduini.dot(m);
            __Cma_du +=  m_q1 * a_q2 * Unduini.dot( m );
        }//m2
    }//q2

#endif

    value_type __Cmf_du=0;
    value_type __Cmm_du=0;
    value_type __Cma_du=0;

    for ( int __q1=0 ; __q1<__Qm; ++__q1 )
    {
        for ( int __m1=0 ; __m1< M_model->mMaxM(__q1); ++__m1 )
        {

            value_type m_q1 = betaMqm[__q1][__m1];

            for ( int __q2=0 ; __q2<__QOutput; ++__q2 )
            {
                for ( int __m2=0 ; __m2< M_model->mMaxF(M_output_index,__q2); ++__m2 )
                {
                    value_type f_q2 = betaFqm[M_output_index][__q2][__m2];
                    //__Cmf_du +=  1./(time_step*time_step) * m_q1 * f_q2 * M_Cmf_du[__q1][__m1][__q2][__m2].head(__N).dot( Unduini );
                    __Cmf_du +=   m_q1 * f_q2 * M_Cmf_du_ini[__q1][__m1][__q2][__m2].head( __N ).dot( Unduini );
                }//m2
            }//q2

            for ( int __q2 = 0; __q2 < __Qm; ++__q2 )
            {
                for ( int __m2=0 ; __m2< M_model->mMaxM(__q2); ++__m2 )
                {
                    value_type m_q2 = betaMqm[__q2][__m2];
                    auto m = M_Cmm_du[__q1][__m1][__q2][__m2].block( 0,0,__N,__N )*Unduini;
                    //__Cmm_du += 1./(time_step*time_step) * m_q1 * m_q2 * Unduini.dot(m);
                    __Cmm_du += m_q1 * m_q2 * Unduini.dot( m );
                }//m2
            }// q2

        }//m1
    }//q1

    //return  math::abs(__c0_du+__Cmf_du+__Caf_du+__Cmm_du+__Cma_du+__Caa_du) ;
    return  math::abs( __c0_du+__Cmf_du+__Cmm_du ) ;

}


template<typename TruthModelType>
typename CRB<TruthModelType>::residual_error_type
CRB<TruthModelType>::transientPrimalResidualEim( int Ncur,parameter_type const& mu,  vectorN_type const& Un ,vectorN_type const& Unold , double time_step, double time ) const
{
    /*
     * transient part needs to be implemented !
     */
    residual_error_type steady_residual_contribution = steadyPrimalResidualEim( Ncur, mu, Un, time );
    std::vector<double> steady_coeff_vector = steady_residual_contribution.template get<1>();
    double delta_pr = steady_residual_contribution.template get<0>();
    value_type __c0_pr     = steady_coeff_vector[0];
    value_type __lambda_pr = steady_coeff_vector[1];
    value_type __gamma_pr  = steady_coeff_vector[2];

    std::vector<double> transient_coeffs_vector;
    transient_coeffs_vector.push_back( __c0_pr );
    transient_coeffs_vector.push_back( __lambda_pr );
    transient_coeffs_vector.push_back( __gamma_pr );

    return boost::make_tuple( delta_pr , transient_coeffs_vector ) ;
}


template<typename TruthModelType>
typename CRB<TruthModelType>::residual_error_type
CRB<TruthModelType>::transientPrimalResidual( int Ncur,parameter_type const& mu,  vectorN_type const& Un ,vectorN_type const& Unold , double time_step, double time ) const
{
    VLOG(2) << fmt::format("CRB : Start transientPrimalresidual with mu={}, Ncur={}, Un.size={}, Unold.size={}, time_step={}, time={}",mu.toString(), Ncur, Un.size(), Unold.size(), time_step, time);
    beta_vector_type betaAqm;
    beta_vector_type betaMqm;
    std::vector<beta_vector_type> betaFqm;
    if( M_model->isLinear() )
    {
        boost::tie( betaMqm, betaAqm, betaFqm ) = M_model->computeBetaQm( mu ,time );
    }
    else
    {
        if ( M_useRbSpaceContextEim && M_hasRbSpaceContextEim )
            boost::tie( boost::tuples::ignore, betaAqm, betaFqm )
                = M_model->computeBetaQm( Un, mu, time );
        else
            boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) =
                M_model->computeBetaQm( this->expansion( Un, Ncur ), mu, time );
    }


    residual_error_type steady_residual_contribution = steadyPrimalResidual( Ncur, mu, Un, time );
    std::vector<double> steady_coeff_vector = steady_residual_contribution.template get<1>();
    value_type __c0_pr     = steady_coeff_vector[0];
    value_type __lambda_pr = steady_coeff_vector[1];
    value_type __gamma_pr  = steady_coeff_vector[2];

    value_type __Cmf_pr=0;
    value_type __Cma_pr=0;
    value_type __Cmm_pr=0;

    //the model can be time-dependant and be executed in steady mode
    //so in that case, we don't need to compute this.
    if ( ! M_model->isSteady() && ! M_model_executed_in_steady_mode )
    {
        size_type __QLhs = M_model->Qa();
        size_type __QRhs = M_model->Ql( 0 );
        size_type __Qm = M_model->Qm();
        size_type __N = Ncur;

        const vectorN_type UnTimeDeriv  = (Un - Unold)/time_step;
        for ( size_type __q1=0 ; __q1<__Qm; ++__q1 )
        {
            for ( size_type __m1=0 ; __m1< M_model->mMaxM(__q1); ++__m1 )
            {
                const value_type m_q1 = betaMqm[__q1][__m1];

                for ( size_type __q2=0 ; __q2<__QRhs; ++__q2 )
                {
                    for ( size_type __m2=0 ; __m2< M_model->mMaxF(0,__q2); ++__m2 )
                    {
                        const value_type f_q2 = betaFqm[0][__q2][__m2];
#if 0
                        __Cmf_pr +=  1./time_step * m_q1 * f_q2  * M_Cmf_pr[__q1][__m1][__q2][__m2].head( __N ).dot( Un );
                        __Cmf_pr -=  1./time_step * m_q1 * f_q2  * M_Cmf_pr[__q1][__m1][__q2][__m2].head( __N ).dot( Unold );
#else
                        __Cmf_pr +=  m_q1 * f_q2  * M_Cmf_pr[__q1][__m1][__q2][__m2].head( __N ).dot( UnTimeDeriv );
#endif
                    }//m2
                }//q2

                for ( size_type __q2 = 0; __q2 < __QLhs; ++__q2 )
                {
                    for ( size_type __m2 = 0; __m2 < M_model->mMaxA(__q2); ++__m2 )
                    {
                        const value_type a_q2 = betaAqm[__q2][__m2];
                        auto const m = M_Cma_pr[__q1][__m1][__q2][__m2].block( 0,0,__N,__N )*Un;
#if 0
                        __Cma_pr += 1./time_step * m_q1 * a_q2 * Un.dot( m );
                        __Cma_pr -= 1./time_step * m_q1 * a_q2 * Unold.dot( m );
#else
                        __Cma_pr += m_q1 * a_q2 * UnTimeDeriv.dot( m );
#endif
                    }//m2
                }//q2

                for ( size_type __q2 = 0; __q2 < __Qm; ++__q2 )
                {
                    for ( size_type __m2 = 0; __m2 < M_model->mMaxM(__q2); ++__m2 )
                    {
                        const value_type m_q2 = betaMqm[__q2][__m2];
#if 0
                        auto m1 = M_Cmm_pr[__q1][__m1][__q2][__m2].block( 0,0,__N,__N )*Un;
                        auto m2 = M_Cmm_pr[__q1][__m1][__q2][__m2].block( 0,0,__N,__N )*Unold;
                        __Cmm_pr += 1./( time_step*time_step ) * m_q1 * m_q2 * Un.dot( m1 );
                        __Cmm_pr -= 1./( time_step*time_step ) * m_q1 * m_q2 * Un.dot( m2 );
                        __Cmm_pr -= 1./( time_step*time_step ) * m_q1 * m_q2 * Unold.dot( m1 );
                        __Cmm_pr += 1./( time_step*time_step ) * m_q1 * m_q2 * Unold.dot( m2 );
#else
                        auto m = M_Cmm_pr[__q1][__m1][__q2][__m2].block( 0,0,__N,__N )*UnTimeDeriv;
                        __Cmm_pr += m_q1 * m_q2 * UnTimeDeriv.dot( m );
#endif
                    }//m2
                }//q2
            }//m1
        }//q1
    }//end of if(! M_model->isSteady() && ! M_model_executed_in_steady_mode )

#if 0
    std::cout<<"[transientResidual]  time "<<time<<std::endl;
    std::cout<<"Un = \n"<<Un<<std::endl;
    std::cout<<"Unold = \n"<<Unold<<std::endl;
    std::cout<<"__C0_pr = "<<__c0_pr<<std::endl;
    std::cout<<"__lambda_pr = "<<__lambda_pr<<std::endl;
    std::cout<<"__gamma_pr = "<<__gamma_pr<<std::endl;
    std::cout<<"__Cmf_pr = "<<__Cmf_pr<<std::endl;
    std::cout<<"__Cma_pr = "<<__Cma_pr<<std::endl;
    std::cout<<"__Cmm_pr = "<<__Cmm_pr<<std::endl;
    std::cout<<"__N = "<<__N<<std::endl;
#endif



    value_type delta_pr =  math::abs( __c0_pr+__lambda_pr+__gamma_pr+__Cmf_pr+__Cma_pr+__Cmm_pr ) ;
    std::vector<double> transient_coeffs_vector;
    transient_coeffs_vector.push_back( __c0_pr );
    transient_coeffs_vector.push_back( __lambda_pr );
    transient_coeffs_vector.push_back( __gamma_pr );
    transient_coeffs_vector.push_back( __Cmf_pr );
    transient_coeffs_vector.push_back( __Cma_pr );
    transient_coeffs_vector.push_back( __Cmm_pr );

    return boost::make_tuple( delta_pr , transient_coeffs_vector ) ;
}



template<typename TruthModelType>
typename CRB<TruthModelType>::residual_error_type
CRB<TruthModelType>::steadyPrimalResidualEim( int Ncur,parameter_type const& mu, vectorN_type const& Un, double time) const
{
    int __QLhs = M_model->Qa();
    int __QRhs = M_model->Ql( 0 );
    int __N = Ncur;

    auto all_eim_interpolation_errors = M_model->eimInterpolationErrorEstimation( mu , Un );

    auto eim_interpolation_errors_A = all_eim_interpolation_errors.template get<1>() ;
    auto eim_interpolation_errors_F = all_eim_interpolation_errors.template get<2>() ;

    std::map<int,double>::iterator it;
    auto endA = eim_interpolation_errors_A.end();
    auto endF = eim_interpolation_errors_F[0].end();

    value_type __c0_pr_eim = 0.0;
    value_type __lambda_pr_eim = 0.0;
    value_type __gamma_pr_eim = 0.0;

    value_type interpolation_error_q1=0;
    value_type interpolation_error_q2=0;

    for ( int __q1 = 0; __q1 < __QRhs; ++__q1 )
    {
        auto itq1 = eim_interpolation_errors_F[0].find(__q1);
        if( itq1 != endF )
            interpolation_error_q1 = itq1->second;
        else
            interpolation_error_q1 = 0;

        for ( int __q2 = 0; __q2 < __QRhs; ++__q2 )
        {
            auto itq2 = eim_interpolation_errors_F[0].find(__q2);
            if( itq2 != endF )
                interpolation_error_q2 = itq2->second;
            else
                interpolation_error_q2 = 0;

            __c0_pr_eim += interpolation_error_q1*interpolation_error_q2 * M_C0_pr_eim[__q1][__q2];
        }//q2
    }//q1

    for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
    {
        auto itq1 = eim_interpolation_errors_A.find(__q1);
        if( itq1 != endA )
            interpolation_error_q1 = itq1->second;
        else
            interpolation_error_q1 = 0;

        for ( int __q2 = 0; __q2 < __QRhs; ++__q2 )
        {
            auto itq2 = eim_interpolation_errors_F[0].find(__q2);
            if( itq2 != endF )
                interpolation_error_q2 = itq2->second;
            else
                interpolation_error_q2 = 0;

            __lambda_pr_eim += interpolation_error_q1*interpolation_error_q2 * M_Lambda_pr_eim[__q1][__q2].head( __N ).dot( Un );

        }//q2
        for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
        {
            auto itq2 = eim_interpolation_errors_A.find(__q2);
            if( itq2 != endA )
                interpolation_error_q2 = itq2->second;
            else
                interpolation_error_q2 = 0;

            auto m = M_Gamma_pr_eim[__q1][__q2].block( 0,0,__N,__N )*Un;
            __gamma_pr_eim += interpolation_error_q1*interpolation_error_q2 * Un.dot( m );

        }//q2
    }//q1

    value_type delta_pr_eim = math::abs( __c0_pr_eim+__lambda_pr_eim+__gamma_pr_eim ) ;
    std::vector<double> coeffs_vector;
    coeffs_vector.push_back( __c0_pr_eim );
    coeffs_vector.push_back( __lambda_pr_eim );
    coeffs_vector.push_back( __gamma_pr_eim );

    return boost::make_tuple( delta_pr_eim,coeffs_vector );

}

template<typename TruthModelType>
typename CRB<TruthModelType>::residual_error_type
CRB<TruthModelType>::steadyPrimalResidual( int Ncur,parameter_type const& mu, vectorN_type const& Un, double time ) const
{

    beta_vector_type betaAqm;
    std::vector<beta_vector_type> betaFqm, betaLqm;
    if( M_model->isLinear() )
    {
        boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) = M_model->computeBetaQm( mu ,time );
    }
    else
    {
        if ( M_useRbSpaceContextEim && M_hasRbSpaceContextEim )
            boost::tie( boost::tuples::ignore, betaAqm, betaFqm )
                = M_model->computeBetaQm( Un, mu, time );
        else
            boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) =
                M_model->computeBetaQm( this->expansion( Un, Ncur ), mu, time );
    }

    int __QLhs = M_model->Qa();
    int __QRhs = M_model->Ql( 0 );
    int __N = Ncur;

    // primal terms
    value_type __c0_pr = 0.0;

    for ( int __q1 = 0; __q1 < __QRhs; ++__q1 )
    {
        for ( int __m1 = 0; __m1 < M_model->mMaxF(0,__q1); ++__m1 )
        {
            value_type fq1 = betaFqm[0][__q1][__m1];
            for ( int __q2 = 0; __q2 < __QRhs; ++__q2 )
            {
                for ( int __m2 = 0; __m2 < M_model->mMaxF(0,__q2); ++__m2 )
                {
                    value_type fq2 = betaFqm[0][__q2][__m2];
                    __c0_pr += M_C0_pr[__q1][__m1][__q2][__m2]*fq1*fq2;
                }//m2
            }//q2
        }//m1
    }//q1

    value_type __lambda_pr = 0.0;
    value_type __gamma_pr = 0.0;


    for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
    {
        for ( int __m1 = 0; __m1 < M_model->mMaxA(__q1); ++__m1 )
        {
            value_type a_q1 = betaAqm[__q1][__m1];

            for ( int __q2 = 0; __q2 < __QRhs; ++__q2 )
            {
                for ( int __m2 = 0; __m2 < M_model->mMaxF(0,__q2); ++__m2 )
                {
                    value_type f_q2 = betaFqm[0][__q2][__m2];
                    __lambda_pr += a_q1*f_q2*M_Lambda_pr[__q1][__m1][__q2][__m2].head( __N ).dot( Un );
                }//m2
            }//q2

            for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
            {
                for ( int __m2 = 0; __m2 < M_model->mMaxA(__q2); ++__m2 )
                {
                    value_type a_q2 = betaAqm[__q2][__m2];
                    auto m = M_Gamma_pr[__q1][__m1][__q2][__m2].block( 0,0,__N,__N )*Un;
                    __gamma_pr += a_q1 * a_q2 * Un.dot( m );
                }//m2
            }//q2
        }//m1
    }//q1

    //value_type delta_pr = math::sqrt( math::abs(__c0_pr+__lambda_pr+__gamma_pr) );
    value_type delta_pr = math::abs( __c0_pr+__lambda_pr+__gamma_pr ) ;
#if 0

    if ( !boost::math::isfinite( delta_pr ) )
    {
        std::cout << "delta_pr is not finite:\n"
                  << "  - c0_pr = " << __c0_pr << "\n"
                  << "  - lambda_pr = " << __lambda_pr << "\n"
                  << "  - gamma_pr = " << __gamma_pr << "\n";
        std::cout << " - betaAqm = " << betaAqm << "\n"
                  << " - betaFqm = " << betaFqm << "\n";
        std::cout << " - Un : " << Un << "\n";
    }

    //std::cout << "delta_pr=" << delta_pr << std::endl;
#endif

    std::vector<double> coeffs_vector;
    coeffs_vector.push_back( __c0_pr );
    coeffs_vector.push_back( __lambda_pr );
    coeffs_vector.push_back( __gamma_pr );

    return boost::make_tuple( delta_pr,coeffs_vector );
}


template<typename TruthModelType>
typename CRB<TruthModelType>::residual_error_type
CRB<TruthModelType>::transientDualResidualEim( int Ncur,parameter_type const& mu,  vectorN_type const& Un ,vectorN_type const& Unold , double time_step, double time ) const
{
    /*
     * transient part needs to be implemented !
     */
    residual_error_type steady_residual_contribution = steadyDualResidualEim( Ncur, mu, Un, time );
    std::vector<double> steady_coeff_vector = steady_residual_contribution.template get<1>();
    double delta_du = steady_residual_contribution.template get<0>();
    value_type __c0_du     = steady_coeff_vector[0];
    value_type __lambda_du = steady_coeff_vector[1];
    value_type __gamma_du  = steady_coeff_vector[2];

    std::vector<double> transient_coeffs_vector;
    transient_coeffs_vector.push_back( __c0_du );
    transient_coeffs_vector.push_back( __lambda_du );
    transient_coeffs_vector.push_back( __gamma_du );

    return boost::make_tuple( delta_du , transient_coeffs_vector ) ;
}

template<typename TruthModelType>
typename CRB<TruthModelType>::residual_error_type
CRB<TruthModelType>::transientDualResidual( int Ncur,parameter_type const& mu,  vectorN_type const& Undu ,vectorN_type const& Unduold , double time_step, double time ) const
{

    if ( M_model->isSteady() )
    {
        residual_error_type steady_residual_contribution = steadyDualResidual( Ncur, mu, Undu, time );
        std::vector<double> steady_coeff_vector = steady_residual_contribution.template get<1>();
        steady_coeff_vector.push_back(0.);
        steady_coeff_vector.push_back(0.);
        steady_coeff_vector.push_back(0.);
        return boost::make_tuple( steady_residual_contribution.template get<0>(), steady_coeff_vector );
    }

    beta_vector_type betaAqm;
    beta_vector_type betaMqm;
    std::vector<beta_vector_type> betaFqm;
    boost::tie( betaMqm, betaAqm, betaFqm) = M_model->computeBetaQm( mu, time );

    value_type __gamma_du = 0.0;
    value_type __Cmf_du=0;
    value_type __Cma_du=0;
    value_type __Cmm_du=0;

    int __QLhs = M_model->Qa();
    int __QOutput = M_model->Ql( M_output_index );
    int __Qm = M_model->Qm();
    int __N = Ncur;

    for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
    {
        for ( int __m1 = 0; __m1 < M_model->mMaxA(__q1); ++__m1 )
        {
            value_type a_q1 = betaAqm[__q1][__m1];

            for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
            {
                for ( int __m2 = 0; __m2 < M_model->mMaxA(__q2); ++__m2 )
                {
                    value_type a_q2 = betaAqm[__q2][__m2]*a_q1;
                    auto m = M_Gamma_du[ __q1][ __m1][ __q2][ __m2].block( 0,0,__N,__N )*Undu;
                    __gamma_du += a_q2*Undu.dot( m );
                }//m2
            }//q2
        }//m1
    }//q1

    const vectorN_type UnduTimeDeriv  = (Undu - Unduold)/time_step;
    for ( int __q1=0 ; __q1<__Qm; ++__q1 )
    {
        for ( int __m1=0 ; __m1< M_model->mMaxM(__q1); ++__m1 )
        {
            const value_type m_q1 = betaMqm[__q1][__m1];

            for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
            {
                for ( int __m2 = 0; __m2 < M_model->mMaxA(__q2); ++__m2 )
                {
                    const value_type a_q2 = betaAqm[__q2][__m2];
                    auto const m = M_Cma_du[__q1][__m1][__q2][__m2].block( 0,0,__N,__N )*Undu;
#if 0
                    __Cma_du += 1./time_step * m_q1 * a_q2 * Undu.dot( m );
                    __Cma_du -= 1./time_step * m_q1 * a_q2 * Unduold.dot( m );
#else
                    __Cma_du += m_q1 * a_q2 * UnduTimeDeriv.dot( m );
#endif
                }//m2
            }//q2

            for ( int __q2 = 0; __q2 < __Qm; ++__q2 )
            {
                for ( int __m2 = 0; __m2 < M_model->mMaxM(__q2); ++__m2 )
                {
                    const value_type m_q2 = betaMqm[__q2][__m2];
#if 0
                    auto m1 = M_Cmm_du[__q1][__m1][__q2][__m2].block( 0,0,__N,__N )*Undu;
                    auto m2 = M_Cmm_du[__q1][__m1][__q2][__m2].block( 0,0,__N,__N )*Unduold;
                    __Cmm_du += 1./( time_step*time_step ) * m_q1 * m_q2 * Undu.dot( m1 );
                    __Cmm_du -= 1./( time_step*time_step ) * m_q1 * m_q2 * Undu.dot( m2 );
                    __Cmm_du -= 1./( time_step*time_step ) * m_q1 * m_q2 * Unduold.dot( m1 );
                    __Cmm_du += 1./( time_step*time_step ) * m_q1 * m_q2 * Unduold.dot( m2 );
#else
                    auto const m = M_Cmm_du[__q1][__m1][__q2][__m2].block( 0,0,__N,__N )*UnduTimeDeriv;
                    __Cmm_du += m_q1 * m_q2 * UnduTimeDeriv.dot( m );
#endif
                } //m2
            } //q2

        } //m1
    }//end of loop over q1


    value_type delta_du =  math::abs( __gamma_du+__Cma_du+__Cmm_du ) ;

#if 0
    double UnduNorm = Undu.norm();
    double UnduoldNorm = Unduold.norm();
    bool print=false;
    if ( !boost::math::isfinite( delta_du ) )
        print=true;
    if( delta_du > 100 || UnduNorm == 0 )
        print=true;
    if( 0 )
    {
        if( this->worldComm().isMasterRank() )
        {
            std::cout<<"[transientDualResidual] time "<<time<<std::endl;
            std::cout<<"Undu = "<<UnduNorm<<std::endl;
            std::cout<<"Unduold = "<<UnduoldNorm<<std::endl;
            std::cout<<"__c0_du = "<<__c0_du<<std::endl;
            std::cout<<"__lambda_du = "<<__lambda_du<<std::endl;
            std::cout<<"__gamma_du = "<<__gamma_du<<std::endl;
            std::cout<<"__Cma_du = "<<__Cma_du<<std::endl;
            std::cout<<"__Cmm_du = "<<__Cmm_du<<std::endl;
            std::cout<<"delta du : "<<delta_du<<std::endl;
            std::cout<<"time step : "<<time_step<<" donc 1/dt*dt = "<<1./( time_step*time_step )<<std::endl;
        }
    }
#endif


    std::vector<double> coeffs_vector;
    coeffs_vector.push_back( 0./*__c0_du*/ );
    coeffs_vector.push_back( 0./*__lambda_du*/ );
    coeffs_vector.push_back( __gamma_du );
    coeffs_vector.push_back( __Cmf_du );
    coeffs_vector.push_back( __Cma_du );
    coeffs_vector.push_back( __Cmm_du );

    return boost::make_tuple( delta_du,coeffs_vector );
}



template<typename TruthModelType>
typename CRB<TruthModelType>::residual_error_type
CRB<TruthModelType>::steadyDualResidualEim( int Ncur,parameter_type const& mu, vectorN_type const& Un, double time) const
{

    int __QLhs = M_model->Qa();
    int __QOutput = M_model->Ql( M_output_index );
    int __N = Ncur;

    auto all_eim_interpolation_errors = M_model->eimInterpolationErrorEstimation( mu , Un );
    auto eim_interpolation_errors_A = all_eim_interpolation_errors.template get<1>() ;
    auto eim_interpolation_errors_F = all_eim_interpolation_errors.template get<2>() ;

    std::map<int,double>::iterator it;
    auto endA = eim_interpolation_errors_A.end();
    auto endF = eim_interpolation_errors_F[M_output_index].end();

    value_type __c0_du_eim = 0.0;
    value_type __lambda_du_eim = 0.0;
    value_type __gamma_du_eim = 0.0;

    value_type interpolation_error_q1=0;
    value_type interpolation_error_q2=0;

    for ( int __q1 = 0; __q1 < __QOutput; ++__q1 )
    {
        auto itq1 = eim_interpolation_errors_F[M_output_index].find(__q1);
        if( itq1 != endF )
            interpolation_error_q1 = itq1->second;
        else
            interpolation_error_q1 = 0;

        for ( int __q2 = 0; __q2 < __QOutput; ++__q2 )
        {
            auto itq2 = eim_interpolation_errors_F[M_output_index].find(__q2);
            if( itq2 != endF )
                interpolation_error_q2 = itq2->second;
            else
                interpolation_error_q2 = 0;

            __c0_du_eim += interpolation_error_q1*interpolation_error_q2 * M_C0_du_eim[__q1][__q2];
        }//q2
    }//q1

    for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
    {
        auto itq1 = eim_interpolation_errors_A.find(__q1);
        if( itq1 != endA )
            interpolation_error_q1 = itq1->second;
        else
            interpolation_error_q1 = 0;

        for ( int __q2 = 0; __q2 < __QOutput; ++__q2 )
        {
            auto itq2 = eim_interpolation_errors_F[M_output_index].find(__q2);
            if( itq2 != endF )
                interpolation_error_q2 = itq2->second;
            else
                interpolation_error_q2 = 0;

            __lambda_du_eim += interpolation_error_q1*interpolation_error_q2 * M_Lambda_du_eim[__q1][__q2].head( __N ).dot( Un );
        }//q2
        for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
        {
            auto itq2 = eim_interpolation_errors_A.find(__q2);
            if( itq2 != endA )
                interpolation_error_q2 = itq2->second;
            else
                interpolation_error_q2 = 0;

            auto m = M_Gamma_du_eim[__q1][__q2].block( 0,0,__N,__N )*Un;
            __gamma_du_eim += interpolation_error_q1*interpolation_error_q2 * Un.dot( m );
        }//q2
    }//q1

    value_type delta_du_eim = math::abs( __c0_du_eim+__lambda_du_eim+__gamma_du_eim ) ;
    std::vector<double> coeffs_vector;
    coeffs_vector.push_back( __c0_du_eim );
    coeffs_vector.push_back( __lambda_du_eim );
    coeffs_vector.push_back( __gamma_du_eim );

    return boost::make_tuple( delta_du_eim,coeffs_vector );

}

template<typename TruthModelType>
typename CRB<TruthModelType>::residual_error_type
CRB<TruthModelType>::steadyDualResidual( int Ncur,parameter_type const& mu, vectorN_type const& Undu, double time ) const
{

    int __QLhs = M_model->Qa();
    int __QOutput = M_model->Ql( M_output_index );
    int __N = Ncur;

    beta_vector_type betaAqm;
    beta_vector_type betaMqm;
    std::vector<beta_vector_type> betaFqm, betaLqm;
    boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) = M_model->computeBetaQm( mu , time );

    value_type __c0_du = 0.0;

    for ( int __q1 = 0; __q1 < __QOutput; ++__q1 )
    {
        for ( int __m1 = 0; __m1 < M_model->mMaxF(M_output_index,__q1); ++__m1 )
        {
            value_type fq1 = betaFqm[M_output_index][__q1][__m1];
            for ( int __q2 = 0; __q2 < __QOutput; ++__q2 )
            {
                for ( int __m2 = 0; __m2 < M_model->mMaxF(M_output_index,__q2); ++__m2 )
                {
                    value_type fq2 = betaFqm[M_output_index][__q2][__m2];
                    __c0_du += M_C0_du[__q1][__m1][__q2][__m2]*fq1*fq2;
                }//m2
            }//q2
        } //m1
    } //q1

    value_type __lambda_du = 0.0;
    value_type __gamma_du = 0.0;

    for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
    {
        for ( int __m1 = 0; __m1 < M_model->mMaxA(__q1); ++__m1 )
        {
            value_type a_q1 = betaAqm[__q1][__m1];
            for ( int __q2 = 0; __q2 < __QOutput; ++__q2 )
            {
                for ( int __m2 = 0; __m2 < M_model->mMaxF(M_output_index,__q2); ++__m2 )
                {
                    value_type f_q2 = betaFqm[M_output_index][__q2][__m2]*a_q1;
                    __lambda_du += f_q2 * M_Lambda_du[__q1][__m1][__q2][__m2].head( __N ).dot( Undu );
                }//m2
            }//q2

            for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
            {
                for ( int __m2 = 0; __m2 < M_model->mMaxA(__q2); ++__m2 )
                {
                    value_type a_q2 = betaAqm[__q2][__m2]*a_q1;
                    auto m = M_Gamma_du[ __q1][ __m1][ __q2][ __m2].block( 0,0,__N,__N )*Undu;
                    __gamma_du += a_q2*Undu.dot( m );
                }//m2
            }//q2
        }//m1
    }//q1

    value_type delta_du = math::abs( __c0_du+__lambda_du+__gamma_du );

#if(0)

    if ( !boost::math::isfinite( delta_du ) )
    {
        std::cout << "delta_du is not finite:\n"
                  << "  - c0_du = " << __c0_du << "\n"
                  << "  - lambda_du = " << __lambda_du << "\n"
                  << "  - gamma_du = " << __gamma_du << "\n";
        std::cout << " - betaAqm = " << betaAqm << "\n"
                  << " - betaFqm = " << betaFqm << "\n";
        std::cout << " - Undu : " << Undu << "\n";
    }

#endif

    std::vector<double> coeffs_vector;
    coeffs_vector.push_back( __c0_du );
    coeffs_vector.push_back( __lambda_du );
    coeffs_vector.push_back( __gamma_du );
    return boost::make_tuple( delta_du,coeffs_vector );
}

template<typename TruthModelType>
void
CRB<TruthModelType>::offlineResidual( int Ncur ,int number_of_added_elements )
{
    return offlineResidual( Ncur, mpl::bool_<model_type::is_time_dependent>(), number_of_added_elements );
}

template<typename TruthModelType>
void
CRB<TruthModelType>::offlineResidualEim( int Ncur ,int number_of_added_elements )
{
    return offlineResidualEim( Ncur, mpl::bool_<model_type::is_time_dependent>(), number_of_added_elements );
}

template<typename TruthModelType>
void
CRB<TruthModelType>::offlineResidual( int Ncur, mpl::bool_<true>, int number_of_added_elements )
{

    Feel:Timer ti;
    int __QLhs = M_model->Qa();
    int __QRhs = M_model->Ql( 0 );
    int __QOutput = M_model->Ql( M_output_index );
    int __Qm = M_model->Qm();
    int __N = Ncur;

    if( this->worldComm().isMasterRank() )
        std::cout << "     o N=" << Ncur << " QLhs=" << __QLhs
                  << " QRhs=" << __QRhs << " Qoutput=" << __QOutput
                  << " Qm=" << __Qm << "\n";

    vector_ptrtype __X( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __Y( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __Fdu( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __Z1(  M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __Z2(  M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __W(  M_backend->newVector( M_model->functionSpace() ) );
    namespace ublas = boost::numeric::ublas;

#if 0
    std::vector< std::vector<sparse_matrix_ptrtype> > Aqm;
    std::vector< std::vector<sparse_matrix_ptrtype> > Mqm;
    std::vector< std::vector<std::vector<vector_ptrtype> > > Fqm;
    std::vector<std::vector<vector_ptrtype> > MFqm;
#endif
    //boost::tie( Mqm, Aqm, Fqm ) = M_model->computeAffineDecomposition();
    boost::tie( M_Mqm, M_Aqm, M_Fqm ) = M_model->computeAffineDecomposition();
    __X->zero();
    __X->add( 1.0 );
    //std::cout << "measure of domain= " << M_model->scalarProduct( __X, __X ) << "\n";
#if 0
    ublas::vector<value_type> mu( P );

    for ( int q = 0; q < P; ++q )
    {
        mu[q] = mut( 0.0 );
    }

#endif

    offlineResidual( Ncur, mpl::bool_<false>(), number_of_added_elements );

    bool optimize = boption(_prefix=M_prefix,_name="crb.optimize-offline-residual") ;

    //the model can be time-dependant and be executed in steady mode
    //so in that case, we don't need to compute this.
    ti.start();
    if( ! M_model_executed_in_steady_mode )
    {

        for ( int elem=__N-number_of_added_elements; elem<__N; elem++ )
        {
            *__X=M_model->rBFunctionSpace()->primalBasisElement(elem);
            for ( int __q1 = 0; __q1 < __Qm; ++__q1 )
            {
                for ( int __m1 = 0; __m1 < M_model->mMaxM(__q1); ++__m1 )
                {
                    M_Mqm[__q1][__m1]->multVector(  __X, __W );
                    __W->scale( -1. );
                    M_model->l2solve( __Z1, __W );

                    for ( int __q2 = 0; __q2 < __QRhs; ++__q2 )
                    {
                        for ( int __m2 = 0; __m2 < M_model->mMaxF(0,__q2); ++__m2 )
                        {

                            M_Cmf_pr[__q1][__m1][__q2][__m2].conservativeResize( __N );
                            M_model->l2solve( __Z2, M_Fqm[0][__q2][__m2] );
                            M_Cmf_pr[ __q1][ __m1][ __q2][ __m2]( elem ) = 2.0*M_model->scalarProduct( __Z1, __Z2 );
                        }// elem
                    }//m2
                }//q2
            }//m1
        }//q1

        if( this->worldComm().isMasterRank() )
            std::cout << "     o M_Cmf_pr updated in " << ti.elapsed() << "s\n";
        ti.start();

        for ( int elem=__N-number_of_added_elements; elem<__N; elem++ )
        {
            *__X=M_model->rBFunctionSpace()->primalBasisElement(elem);

            for ( int __q1 = 0; __q1 < __Qm; ++__q1 )
            {
                for ( int __m1 = 0; __m1 < M_model->mMaxM(__q1); ++__m1 )
                {
                    M_Mqm[__q1][__m1]->multVector(  __X, __W );
                    __W->scale( -1. );
                    M_model->l2solve( __Z1, __W );

                    for ( int __l = 0; __l < ( int )__N; ++__l )
                    {

                        *__Y=M_model->rBFunctionSpace()->primalBasisElement(__l);

                        for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
                        {
                            for ( int __m2 = 0; __m2 < M_model->mMaxA(__q2); ++__m2 )
                            {
                                M_Cma_pr[__q1][__m1][__q2][__m2].conservativeResize( __N, __N );
                                M_Aqm[__q2][__m2]->multVector(  __Y, __W );
                                __W->scale( -1. );
                                M_model->l2solve( __Z2, __W );
                                M_Cma_pr[ __q1][ __m1][ __q2][ __m2]( elem,__l ) = 2.0*M_model->scalarProduct( __Z1, __Z2 );
                            }// m2
                        }// q2
                    }//__l
                }//m1
            }// q1
        }//elem
        for ( int __j = 0; __j < ( int )__N; ++__j )
        {
            *__X=M_model->rBFunctionSpace()->primalBasisElement(__j);

            for ( int __q1 = 0; __q1 < __Qm; ++__q1 )
            {
                for ( int __m1 = 0; __m1 < M_model->mMaxM(__q1); ++__m1 )
                {
                    M_Mqm[__q1][__m1]->multVector(  __X, __W );
                    __W->scale( -1. );
                    M_model->l2solve( __Z1, __W );

                    for ( int elem=__N-number_of_added_elements; elem<__N; elem++ )
                    {

                        *__Y=M_model->rBFunctionSpace()->primalBasisElement(elem);

                        for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
                        {
                            for ( int __m2 = 0; __m2 < M_model->mMaxA(__q2); ++__m2 )
                            {

                                M_Aqm[__q2][ __m2]->multVector(  __Y, __W );
                                __W->scale( -1. );
                                M_model->l2solve( __Z2, __W );
                                M_Cma_pr[ __q1][ __m1][ __q2][ __m2]( __j,elem ) = 2.0*M_model->scalarProduct( __Z1, __Z2 );
                            }//m2
                        }// q2

                    }// elem
                }// m1
            }// q1
        } // __j

        if( this->worldComm().isMasterRank() )
            std::cout << "     o M_Cma_pr updated in " << ti.elapsed() << "s\n";
        ti.start();

        for ( int elem=__N-number_of_added_elements; elem<__N; elem++ )
        {

            for ( int __q1 = 0; __q1 < __Qm; ++__q1 )
            {
                *__X=M_model->rBFunctionSpace()->primalBasisElement(elem);

                for ( int __m1 = 0; __m1 < M_model->mMaxM(__q1); ++__m1 )
                {
                    M_Mqm[__q1][__m1]->multVector(  __X, __W );
                    __W->scale( -1. );
                    M_model->l2solve( __Z1, __W );

                    for ( int __l = 0; __l < ( int )__N; ++__l )
                    {
                        *__Y=M_model->rBFunctionSpace()->primalBasisElement(__l);

                        if( optimize )
                        {
                            //case we don't use EIM
                            for ( int __q2 = 0; __q2 < __q1; ++__q2 )
                            {
                                for ( int __m2 = 0; __m2 < M_model->mMaxM(__q2); ++__m2 )
                                {

                                    M_Cmm_pr[__q1][__m1][__q2][__m2].conservativeResize( __N, __N );
                                    M_Cmm_pr[__q2][__m2][__q1][__m1].conservativeResize( __N, __N );

                                    M_Mqm[__q2][__m2]->multVector(  __Y, __W );
                                    __W->scale( -1. );
                                    M_model->l2solve( __Z2, __W );
                                    double prod = M_model->scalarProduct( __Z1, __Z2 );
                                    M_Cmm_pr[ __q1][ __m1][ __q2][ __m2]( elem,__l ) = prod;
                                    M_Cmm_pr[ __q2][ __m2][ __q1][ __m1]( __l,elem ) = prod;
                                } // m2
                            } // q2
                            //diagonal elements
                            M_Mqm[__q1][__m1]->multVector(  __Y, __W );
                            __W->scale( -1. );
                            M_model->l2solve( __Z2, __W );
                            M_Cmm_pr[__q1][__m1][__q1][__m1].conservativeResize( __N, __N );
                            double prod = M_model->scalarProduct( __Z1, __Z2 );
                            M_Cmm_pr[ __q1][ __m1][ __q1][ __m1]( elem,__l ) = prod;

                        }// optimize
                        else
                        {

                            for ( int __q2 = 0; __q2 < __Qm; ++__q2 )
                            {
                                for ( int __m2 = 0; __m2 < M_model->mMaxM(__q2); ++__m2 )
                                {

                                    M_Cmm_pr[__q1][__m1][__q2][__m2].conservativeResize( __N, __N );

                                    M_Mqm[__q2][__m2]->multVector(  __Y, __W );
                                    __W->scale( -1. );
                                    M_model->l2solve( __Z2, __W );

                                    M_Cmm_pr[ __q1][ __m1][ __q2][ __m2]( elem,__l ) = M_model->scalarProduct( __Z1, __Z2 );
                                }// m2
                            }// q2
                        }//no optimize
                    }// __l
                }// m1
            } // q1
        }// elem

        for ( int __j = 0; __j < ( int )__N; ++__j )
        {
            *__X=M_model->rBFunctionSpace()->primalBasisElement(__j);
            for ( int __q1 = 0; __q1 < __Qm; ++__q1 )
            {
                for ( int __m1 = 0; __m1 < M_model->mMaxM(__q1); ++__m1 )
                {
                    M_Mqm[__q1][ __m1]->multVector(  __X, __W );
                    __W->scale( -1. );
                    M_model->l2solve( __Z1, __W );

                    for ( int elem=__N-number_of_added_elements; elem<__N; elem++ )
                    {
                        *__Y=M_model->rBFunctionSpace()->primalBasisElement(elem);

                        if( optimize )
                        {
                            //case we don't use EIM
                            for ( int __q2 = 0; __q2 < __q1; ++__q2 )
                            {
                                for ( int __m2 = 0; __m2 < M_model->mMaxM(__q2); ++__m2 )
                                {
                                    M_Mqm[__q2][__m2]->multVector(  __Y, __W );
                                    __W->scale( -1. );
                                    M_model->l2solve( __Z2, __W );
                                    double prod = M_model->scalarProduct( __Z1, __Z2 );
                                    M_Cmm_pr[ __q1][ __m1][ __q2][ __m2]( __j,elem ) = prod;
                                    M_Cmm_pr[ __q2][ __m2][ __q1][ __m1]( elem,__j ) = prod;
                                } // m2
                            } // q2
                            //diagonal elements
                            M_Mqm[__q1][__m1]->multVector(  __Y, __W );
                            __W->scale( -1. );
                            M_model->l2solve( __Z2, __W );
                            double prod = M_model->scalarProduct( __Z1, __Z2 );
                            M_Cmm_pr[ __q1][ __m1][ __q1][ __m1]( __j,elem ) = prod;

                        }// optimize
                        else
                        {
                            for ( int __q2 = 0; __q2 < __Qm; ++__q2 )
                            {
                                for ( int __m2 = 0; __m2 < M_model->mMaxM(__q2); ++__m2 )
                                {
                                    M_Mqm[__q2][ __m2]->multVector(  __Y, __W );
                                    __W->scale( -1. );
                                    M_model->l2solve( __Z2, __W );

                                    M_Cmm_pr[ __q1][ __m1][ __q2][ __m2]( __j,elem ) = M_model->scalarProduct( __Z1, __Z2 );

                                }// m2
                            }// q2
                        }// no optimize
                    }// elem
                } // m1
            }// q1
        } // __j

        if( this->worldComm().isMasterRank() )
            std::cout << "     o M_Cmm_pr updated in " << ti.elapsed() << "s\n";


        sparse_matrix_ptrtype Atq1 = M_model->newMatrix();
        sparse_matrix_ptrtype Atq2 = M_model->newMatrix();
        sparse_matrix_ptrtype Mtq1 = M_model->newMatrix();
        sparse_matrix_ptrtype Mtq2 = M_model->newMatrix();
        //
        // Dual
        //

        ti.start();

        LOG(INFO) << "[offlineResidual] Cmf_du Cma_du Cmm_du\n";
        for ( int elem=__N-number_of_added_elements; elem<__N; elem++ )
        {
            *__X=M_model->rBFunctionSpace()->dualBasisElement(elem);

            for ( int __q1 = 0; __q1 < __Qm; ++__q1 )
            {
                for ( int __m1 = 0; __m1 < M_model->mMaxM(__q1); ++__m1 )
                {

                    M_Mqm[__q1][__m1]->multVector(  __X, __W );
                    __W->close();
                    __W->scale( -1. );
                    M_model->l2solve( __Z1, __W );

                    for ( int __q2 = 0; __q2 < __QOutput; ++__q2 )
                    {
                        for ( int __m2 = 0; __m2 < M_model->mMaxF(M_output_index,__q2); ++__m2 )
                        {

                            M_Cmf_du[__q1][__m1][__q2][__m2].conservativeResize( __N );
                            M_Cmf_du_ini[__q1][__m1][__q2][__m2].conservativeResize( __N );

                            *__Fdu = *M_Fqm[M_output_index][__q2][__m2];
                            __Fdu->close();
                            __Fdu->scale( -1.0 );
                            M_model->l2solve( __Z2, __Fdu );

                            M_Cmf_du[ __q1][__m1][ __q2][__m2]( elem ) = 2.0*M_model->scalarProduct( __Z1, __Z2 );

                            *__Fdu = *M_Fqm[M_output_index][__q2][__m2];
                            __Fdu->close();
                            M_model->l2solve( __Z2, __Fdu );
                            M_Cmf_du_ini[ __q1][__m1][ __q2][__m2]( elem ) = 2.0*M_model->scalarProduct( __Z1, __Z2 );
                        }// m2
                    } // q2
                } // m1
            } // q1
        } // elem


        if( this->worldComm().isMasterRank() )
            std::cout << "     o M_Cmf_du updated in " << ti.elapsed() << "s\n";

        ti.start();

        for ( int elem=__N-number_of_added_elements; elem<__N; elem++ )
        {
            *__X=M_model->rBFunctionSpace()->dualBasisElement(elem);

            for ( int __q1 = 0; __q1 < __Qm; ++__q1 )
            {
                for ( int __m1 = 0; __m1 < M_model->mMaxM(__q1); ++__m1 )
                {

                    M_Mqm[__q1][__m1]->multVector(  __X, __W );
                    __W->scale( -1. );
                    M_model->l2solve( __Z1, __W );

                    for ( int __l = 0; __l < ( int )__N; ++__l )
                    {

                        *__Y=M_model->rBFunctionSpace()->dualBasisElement(__l);

                        for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
                        {
                            for ( int __m2 = 0; __m2 < M_model->mMaxA(__q2); ++__m2 )
                            {

                                if( boption(_prefix=M_prefix,_name="crb.use-symmetric-matrix") )
                                    Atq2 = M_Aqm[__q2][__m2];
                                else
                                    M_Aqm[__q2][__m2]->transpose( Atq2 );

                                M_Cma_du[__q1][__m1][__q2][__m2].conservativeResize( __N, __N );

                                Atq2->multVector(  __Y, __W );
                                __W->scale( -1. );
                                M_model->l2solve( __Z2, __W );

                                M_Cma_du[ __q1][__m1][__q2][ __m2]( elem,__l ) = 2.0*M_model->scalarProduct( __Z1, __Z2 );
                            }// m2
                        }// q2
                    }//__l
                }// m1
            }//q1
        }//elem

        for ( int __j = 0; __j < ( int )__N; ++__j )
        {
            *__X=M_model->rBFunctionSpace()->dualBasisElement(__j);
            for ( int __q1 = 0; __q1 < __Qm; ++__q1 )
            {
                for ( int __m1 = 0; __m1 < M_model->mMaxM(__q1); ++__m1 )
                {

                    M_Mqm[__q1][__m1]->multVector(  __X, __W );
                    __W->scale( -1. );
                    M_model->l2solve( __Z1, __W );

                    for ( int elem=__N-number_of_added_elements; elem<__N; elem++ )
                    {
                        *__Y=M_model->rBFunctionSpace()->dualBasisElement(elem);

                        for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
                        {
                            for ( int __m2 = 0; __m2 < M_model->mMaxA(__q2); ++__m2 )
                            {

                                if( boption(_prefix=M_prefix,_name="crb.use-symmetric-matrix") )
                                    Atq2 = M_Aqm[__q2][__m2];
                                else
                                    M_Aqm[__q2][__m2]->transpose( Atq2 );


                                Atq2->multVector(  __Y, __W );
                                __W->scale( -1. );
                                M_model->l2solve( __Z2, __W );

                                M_Cma_du[ __q1][__m1][__q2][ __m2]( __j,elem ) = 2.0*M_model->scalarProduct( __Z1, __Z2 );
                            }// m2
                        }// q2

                    }// elem
                } // m1
            }// q1
        } // __j


        if( this->worldComm().isMasterRank() )
            std::cout << "     o Cma_du updated in " << ti.elapsed() << "s\n";

        ti.start();


        for ( int elem=__N-number_of_added_elements; elem<__N; elem++ )
        {
            *__X=M_model->rBFunctionSpace()->dualBasisElement(elem);

            for ( int __q1 = 0; __q1 < __Qm; ++__q1 )
            {
                for ( int __m1 = 0; __m1 < M_model->mMaxM(__q1); ++__m1 )
                {

                    M_Mqm[__q1][__m1]->multVector(  __X, __W );
                    __W->scale( -1. );
                    M_model->l2solve( __Z1, __W );

                    for ( int __l = 0; __l < ( int )__N; ++__l )
                    {
                        *__Y=M_model->rBFunctionSpace()->dualBasisElement(__l);

                        if( optimize )
                        {
                            //case we don't use EIM
                            for ( int __q2 = 0; __q2 < __q1; ++__q2 )
                            {
                                for ( int __m2 = 0; __m2 < M_model->mMaxM(__q2); ++__m2 )
                                {

                                    M_Cmm_du[__q1][__m1][__q2][__m2].conservativeResize( __N, __N );
                                    M_Cmm_du[__q2][__m2][__q1][__m1].conservativeResize( __N, __N );

                                    M_Mqm[__q2][__m2]->multVector(  __Y, __W );
                                    __W->scale( -1. );
                                    M_model->l2solve( __Z2, __W );
                                    double prod = M_model->scalarProduct( __Z1, __Z2 );
                                    M_Cmm_du[ __q1][ __m1][ __q2][ __m2]( elem,__l ) = prod;
                                    M_Cmm_du[ __q2][ __m2][ __q1][ __m1]( __l,elem ) = prod;
                                } // m2
                            } // q2
                            //diagonal elements
                            M_Mqm[__q1][__m1]->multVector(  __Y, __W );
                            __W->scale( -1. );
                            M_model->l2solve( __Z2, __W );
                            M_Cmm_du[__q1][__m1][__q1][__m1].conservativeResize( __N, __N );
                            double prod = M_model->scalarProduct( __Z1, __Z2 );
                            M_Cmm_du[ __q1][ __m1][ __q1][ __m1]( elem,__l ) = prod;

                        }// optimize
                        else
                        {
                            for ( int __q2 = 0; __q2 < __Qm; ++__q2 )
                            {
                                for ( int __m2 = 0; __m2 < M_model->mMaxM(__q2); ++__m2 )
                                {
                                    M_Cmm_du[__q1][__m1][__q2][__m2].conservativeResize( __N, __N );

                                    M_Mqm[__q2][__m2]->multVector(  __Y, __W );
                                    __W->scale( -1. );
                                    M_model->l2solve( __Z2, __W );

                                    M_Cmm_du[ __q1][ __m1][ __q2][ __m2]( elem,__l ) = M_model->scalarProduct( __Z1, __Z2 );
                                }//m2
                            }//q2
                        }
                    }//__l
                }//m1
            }//q1
        }//elem

        for ( int __j = 0; __j < ( int )__N; ++__j )
        {
            *__X=M_model->rBFunctionSpace()->dualBasisElement(__j);

            for ( int __q1 = 0; __q1 < __Qm; ++__q1 )
            {
                for ( int __m1 = 0; __m1 < M_model->mMaxM(__q1); ++__m1 )
                {

                    M_Mqm[__q1][__m1]->multVector(  __X, __W );
                    __W->scale( -1. );
                    M_model->l2solve( __Z1, __W );

                    for ( int elem=__N-number_of_added_elements; elem<__N; elem++ )
                    {
                        *__Y=M_model->rBFunctionSpace()->dualBasisElement(elem);

                        if( optimize )
                        {
                            //case we don't use EIM
                            for ( int __q2 = 0; __q2 < __q1; ++__q2 )
                            {
                                for ( int __m2 = 0; __m2 < M_model->mMaxM(__q2); ++__m2 )
                                {
                                    M_Mqm[__q2][__m2]->multVector(  __Y, __W );
                                    __W->scale( -1. );
                                    M_model->l2solve( __Z2, __W );
                                    double prod = M_model->scalarProduct( __Z1, __Z2 );
                                    M_Cmm_du[ __q1][ __m1][ __q2][ __m2]( __j,elem ) = prod;
                                    M_Cmm_du[ __q2][ __m2][ __q1][ __m1]( elem,__j ) = prod;
                                } // m2
                            } // q2
                            //diagonal elements
                            M_Mqm[__q1][__m1]->multVector(  __Y, __W );
                            __W->scale( -1. );
                            M_model->l2solve( __Z2, __W );
                            double prod = M_model->scalarProduct( __Z1, __Z2 );
                            M_Cmm_du[ __q1][ __m1][ __q1][ __m1]( __j,elem ) = prod;

                        }// optimize
                        else
                        {
                            for ( int __q2 = 0; __q2 < __Qm; ++__q2 )
                            {
                                for ( int __m2 = 0; __m2 < M_model->mMaxM(__q2); ++__m2 )
                                {
                                    M_Mqm[__q2][__m2]->multVector(  __Y, __W );
                                    __W->scale( -1. );
                                    M_model->l2solve( __Z2, __W );

                                    M_Cmm_du[ __q1][ __m1][ __q2][ __m2]( __j,elem ) = M_model->scalarProduct( __Z1, __Z2 );
                                }// m2
                            }// q2
                        }//no optimize
                    }// elem
                } // m1
            }// q1
        } // __j


        if( this->worldComm().isMasterRank() )
            std::cout << "     o Cmm_du updated in " << ti.elapsed() << "s\n";
        ti.start();

    }//if( ! M_model_executed_in_steady_mode )

    LOG(INFO) << "[offlineResidual] Done.\n";

}



template<typename TruthModelType>
void
CRB<TruthModelType>::offlineResidualEim( int Ncur, mpl::bool_<true>, int number_of_added_elements )
{
    offlineResidualEim( Ncur, mpl::bool_<false>(), number_of_added_elements );
}

template<typename TruthModelType>
void
CRB<TruthModelType>::offlineResidual( int Ncur, mpl::bool_<false> , int number_of_added_elements )
{
    if ( ioption(_prefix=M_prefix,_name="crb.offline-residual-version" ) == 0 )
        offlineResidualV0( Ncur,mpl::bool_<false>(),number_of_added_elements );
    else
        offlineResidualV1( Ncur,mpl::bool_<false>(),number_of_added_elements );
}
template<typename TruthModelType>
void
CRB<TruthModelType>::offlineResidualV0( int Ncur, mpl::bool_<false> , int number_of_added_elements )
{
    Feel::Timer ti;
    int __QLhs = M_model->Qa();
    int __QRhs = M_model->Ql( 0 );
    int __QOutput = M_model->Ql( M_output_index );
    int __N = Ncur;

    bool use_ser = ioption(_prefix=M_prefix,_name="ser.rb-frequency");
    int added_elements = use_ser ? __N:number_of_added_elements;


    if( this->worldComm().isMasterRank() )
        std::cout << "     o N=" << Ncur << " QLhs=" << __QLhs
                  << " QRhs=" << __QRhs << " Qoutput=" << __QOutput << "\n";

    vector_ptrtype __X( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __Y( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __Fdu( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __Z1(  M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __Z2(  M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __W(  M_backend->newVector( M_model->functionSpace() ) );
    namespace ublas = boost::numeric::ublas;

#if 0
    std::vector< std::vector<sparse_matrix_ptrtype> > Aqm,Mqm;
    std::vector< std::vector<vector_ptrtype> > MFqm;
    std::vector< std::vector<std::vector<vector_ptrtype> > > Fqm,Lqm;
#endif
    //boost::tie( Mqm, Aqm, Fqm ) = M_model->computeAffineDecomposition();
    boost::tie( M_Mqm, M_Aqm, M_Fqm ) = M_model->computeAffineDecomposition();
    __X->zero();
    __X->add( 1.0 );

    if( this->worldComm().isMasterRank() )
        std::cout << "     o initialize offlineResidual in " << ti.elapsed() << "s\n";

    //std::cout << "measure of domain= " << M_model->scalarProduct( __X, __X ) << "\n";
#if 0
    ublas::vector<value_type> mu( P );

    for ( int q = 0; q < P; ++q )
    {
        mu[q] = mut( 0.0 );
    }

#endif


    // Primal
    // no need to recompute this term each time
    if ( Ncur == M_Nm || use_ser )
    {
        ti.start();
        LOG(INFO) << "[offlineResidual] Compute Primal residual data\n";
        LOG(INFO) << "[offlineResidual] C0_pr\n";

        // see above Z1 = C^-1 F and Z2 = F
        for ( int __q1 = 0; __q1 < __QRhs; ++__q1 )
        {
            for ( int __m1 = 0; __m1 < M_model->mMaxF(0,__q1); ++__m1 )
            {
                //LOG(INFO) << "__Fq->norm1=" << Fq[0][__q1][__m1]->l2Norm() << "\n";
                M_model->l2solve( __Z1, M_Fqm[0][__q1][__m1] );
                //for ( int __q2 = 0;__q2 < __q1;++__q2 )
                for ( int __q2 = 0; __q2 < __QRhs; ++__q2 )
                {
                    for ( int __m2 = 0; __m2 < M_model->mMaxF(0,__q2); ++__m2 )
                    {
                        //LOG(INFO) << "__Fq->norm 2=" << Fq[0][__q2][__m2]->l2Norm() << "\n";
                        M_model->l2solve( __Z2, M_Fqm[0][__q2][__m2] );
                        //M_C0_pr[__q1][__m1][__q2][__m2] = M_model->scalarProduct( __X, Fq[0][__q2][__m2] );
                        M_C0_pr[__q1][__m1][__q2][__m2] = M_model->scalarProduct( __Z1, __Z2 );
                        //M_C0_pr[__q2][__q1] = M_C0_pr[__q1][__q2];
                        //VLOG(1) <<"M_C0_pr[" << __q1 << "][" << __m1 << "]["<< __q2 << "][" << __m2 << "]=" << M_C0_pr[__q1][__m1][__q2][__m2] << "\n";
                        //LOG(INFO) << "M_C0_pr[" << __q1 << "][" << __m1 << "]["<< __q2 << "][" << __m2 << "]=" << M_C0_pr[__q1][__m1][__q2][__m2] << "\n";
                    }//end of loop __m2
                }//end of loop __q2
                //M_C0_pr[__q1][__q1] = M_model->scalarProduct( __X, __X );
            }//end of loop __m1
        }//end of loop __q1

        if( this->worldComm().isMasterRank() )
            std::cout << "     o M_C0_pr updated in " << ti.elapsed() << "s\n";

    }// Ncur==M_Nm


#if 0
    parameter_type const& mu = M_WNmu->at( 0 );
    //std::cout << "[offlineResidual] mu=" << mu << "\n";
    beta_vector_type betaAqm;
    std::vector<beta_vector_type> betaFqm;
    boost::tie( boost::tuples::ignore, betaAqm, betaFqm ) = M_model->computeBetaQm( mu );
    value_type __c0_pr = 0.0;

    for ( int __q1 = 0; __q1 < __QRhs; ++__q1 )
    {
        for ( int __m1 = 0; __m1 < M_model->mMaxF(0,__q1); ++__m1 )
        {
            value_type fq1 = betaFqm[0][__q1][__m1];

            for ( int __q2 = 0; __q2 < __QRhs; ++__q2 )
            {
                for ( int __m2 = 0; __m2 < M_model->mMaxF(M_output_index,__q2); ++__m2 )
                {
                    value_type fq2 = betaFqm[0][__q2][__m2];
                    __c0_pr += M_C0_pr[__q1][__m1][__q2][__m2]*fq1*fq2;
                }
            }
        }//end of loop __m1
    }//end of loop __q1


    //std::cout << "c0=" << __c0_pr << std::endl;

    std::vector< sparse_matrix_ptrtype > A;
    std::vector< std::vector<vector_ptrtype> > F;
    boost::tie( A, F ) = M_model->update( mu );
    M_model->l2solve( __X, F[0] );
    //std::cout << "c0 2 = " << M_model->scalarProduct( __X, __X ) << "\n";;
#endif

    ti.start();

    bool optimize = boption(_prefix=M_prefix,_name="crb.optimize-offline-residual") ;

    //
    //  Primal
    //
    LOG(INFO) << "[offlineResidual] Lambda_pr, Gamma_pr\n";

    for ( int elem=__N-added_elements; elem<__N; elem++ )
    {
        *__X=M_model->rBFunctionSpace()->primalBasisElement(elem);

        for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
        {
            for ( int __m1 = 0; __m1 < M_model->mMaxA(__q1) ; ++__m1 )
            {
                M_Aqm[__q1][__m1]->multVector(  __X, __W );
                __W->scale( -1. );
                //std::cout << "__W->norm=" << __W->l2Norm() << "\n";
                M_model->l2solve( __Z1, __W );
                for ( int __q2 = 0; __q2 < __QRhs; ++__q2 )
                {
                    for ( int __m2 = 0; __m2 < M_model->mMaxF(0, __q2); ++__m2 )
                    {
                        M_Lambda_pr[__q1][__m1][__q2][__m2].conservativeResize( __N );
                        //__Y = Fq[0][__q2];
                        //std::cout << "__Fq->norm=" << Fq[0][__q2]->l2Norm() << "\n";
                        M_model->l2solve( __Z2, M_Fqm[0][__q2][__m2] );
                        M_Lambda_pr[ __q1][ __m1][ __q2][ __m2]( elem ) = 2.0*M_model->scalarProduct( __Z1, __Z2 );
                        //VLOG(1) << "M_Lambda_pr[" << __q1 << "][" << __m1 << "][" << __q2 << "][" << __m2 << "][" << __j << "]=" <<  M_Lambda_pr[__q1][__m1][__q2][__m2][__j] << "\n";
                        //std::cout << "M_Lambda_pr[" << __q1 << "][" << __m1 << "][" << __q2 << "][" << __m2 << "][" << __j << "]=" << M_Lambda_pr[__q1][__m1][__q2][__m2][__j] << "\n";
                    }//m2
                }//end of __q2
            }//end of loop __m1
        }//end of __q1
    }//elem

    if( this->worldComm().isMasterRank() )
        std::cout << "     o Lambda_pr updated in " << ti.elapsed() << "s\n";

    ti.start();

    for ( int elem=__N-added_elements; elem<__N; elem++ )
    {
        *__X=M_model->rBFunctionSpace()->primalBasisElement(elem);
        for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
        {
            for ( int __m1 = 0; __m1 < M_model->mMaxA(__q1); ++__m1 )
            {
                M_Aqm[__q1][__m1]->multVector(  __X, __W );
                __W->scale( -1. );
                M_model->l2solve( __Z1, __W );

                for ( int __l = 0; __l < ( int )__N; ++__l )
                {
                    *__Y=M_model->rBFunctionSpace()->primalBasisElement(__l);

                    if( optimize )
                    {
                        //case we don't use EIM
                        for ( int __q2 = 0; __q2 < __q1; ++__q2 )
                        {
                            for ( int __m2 = 0; __m2 < M_model->mMaxA(__q2); ++__m2 )
                            {
                                M_Gamma_pr[__q1][__m1][__q2][__m2].conservativeResize( __N, __N );
                                M_Gamma_pr[__q2][__m2][__q1][__m1].conservativeResize( __N, __N );

                                M_Aqm[__q2][__m2]->multVector(  __Y, __W );
                                __W->scale( -1. );
                                M_model->l2solve( __Z2, __W );
                                double prod = M_model->scalarProduct( __Z1, __Z2 );
                                M_Gamma_pr[ __q1][ __m1][ __q2][ __m2]( elem,__l ) = prod;
                                M_Gamma_pr[ __q2][ __m2][ __q1][ __m1]( __l,elem ) = prod;
                            } // m2
                        } // q2
                        //diagonal elements
                        M_Aqm[__q1][__m1]->multVector(  __Y, __W );
                        __W->scale( -1. );
                        M_model->l2solve( __Z2, __W );
                        M_Gamma_pr[__q1][__m1][__q1][__m1].conservativeResize( __N, __N );
                        double prod = M_model->scalarProduct( __Z1, __Z2 );
                        M_Gamma_pr[ __q1][ __m1][ __q1][ __m1]( elem,__l ) = prod;
                    }//optimize
                    else
                    {
                        for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
                        {
                            for ( int __m2 = 0; __m2 < M_model->mMaxA(__q2); ++__m2 )
                            {
                                M_Gamma_pr[__q1][__m1][__q2][__m2].conservativeResize( __N, __N );

                                M_Aqm[__q2][__m2]->multVector(  __Y, __W );
                                __W->scale( -1. );
                                M_model->l2solve( __Z2, __W );
                                M_Gamma_pr[ __q1][ __m1][ __q2][ __m2]( elem,__l ) = M_model->scalarProduct( __Z1, __Z2 );
                            } // m2
                        } // q2
                    }//no optimize
                }//end of loop over l
            } // m1
        } // q1
    } // elem

    for ( int __j = 0; __j < ( int )__N; ++__j )
    {
        *__X=M_model->rBFunctionSpace()->primalBasisElement(__j);
        for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
        {
            for ( int __m1 = 0; __m1 < M_model->mMaxA(__q1); ++__m1 )
            {
                M_Aqm[__q1][__m1]->multVector(  __X, __W );
                __W->scale( -1. );
                M_model->l2solve( __Z1, __W );

                //column N-1
                //int __l = __N-1;
                for ( int elem=__N-added_elements; elem<__N; elem++ )
                {
                    *__Y=M_model->rBFunctionSpace()->primalBasisElement(elem);

                    if( optimize )
                    {
                        //case we don't use EIM
                        for ( int __q2 = 0; __q2 < __q1; ++__q2 )
                        {
                            for ( int __m2 = 0; __m2 < M_model->mMaxA(__q2); ++__m2 )
                            {
                                M_Aqm[__q2][__m2]->multVector(  __Y, __W );
                                __W->scale( -1. );
                                M_model->l2solve( __Z2, __W );
                                double prod = M_model->scalarProduct( __Z1, __Z2 );
                                M_Gamma_pr[ __q1][ __m1][ __q2][ __m2]( __j,elem ) = prod;
                                M_Gamma_pr[ __q2][ __m2][ __q1][ __m1]( elem,__j ) = prod;
                            } // m2
                        } // q2
                        //diagonal elements
                        M_Aqm[__q1][__m1]->multVector(  __Y, __W );
                        __W->scale( -1. );
                        M_model->l2solve( __Z2, __W );
                        double prod = M_model->scalarProduct( __Z1, __Z2 );
                        M_Gamma_pr[ __q1][ __m1][ __q1][ __m1]( __j,elem ) = prod;
                    }//optimize
                    else
                    {
                        for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
                        {
                            for ( int __m2 = 0; __m2 < M_model->mMaxA(__q2); ++__m2 )
                            {
                                M_Aqm[__q2][__m2]->multVector(  __Y, __W );
                                __W->scale( -1. );
                                M_model->l2solve( __Z2, __W );
                                M_Gamma_pr[ __q1][ __m1][ __q2][ __m2]( __j,elem ) = M_model->scalarProduct( __Z1, __Z2 );
                            } // m2
                        } // q2
                    }//no optimize
                }// end of loop elem
            }// m1
        }// q1
    }// end of loop __j

    if( this->worldComm().isMasterRank() )
        std::cout << "     o Gamma_pr updated in " << ti.elapsed() << "s\n";
    sparse_matrix_ptrtype Atq1 = M_model->newMatrix();
    sparse_matrix_ptrtype Atq2 = M_model->newMatrix();
    ti.start();

    //
    // Dual
    //
    // compute this only once
    if( M_solve_dual_problem )
    {
        if ( Ncur == M_Nm )
        {
            LOG(INFO) << "[offlineResidual] Compute Dual residual data\n";
            LOG(INFO) << "[offlineResidual] C0_du\n";

            for ( int __q1 = 0; __q1 < __QOutput; ++__q1 )
            {
                for ( int __m1 = 0; __m1 < M_model->mMaxF(M_output_index,__q1); ++__m1 )
                {
                    *__Fdu = *M_Fqm[M_output_index][__q1][__m1];
                    __Fdu->close();
                    __Fdu->scale( -1.0 );
                    M_model->l2solve( __Z1, __Fdu );

                    for ( int __q2 = 0; __q2 < __QOutput; ++__q2 )
                    {
                        for ( int __m2 = 0; __m2 < M_model->mMaxF(M_output_index,__q2); ++__m2 )
                        {
                            *__Fdu = *M_Fqm[M_output_index][__q2][__m2];
                            __Fdu->close();
                            __Fdu->scale( -1.0 );
                            M_model->l2solve( __Z2, __Fdu );
                            M_C0_du[__q1][__m1][__q2][__m2] = M_model->scalarProduct( __Z1, __Z2 );
                            //M_C0_du[__q2][__q1] = M_C0_du[__q1][__q2];
                            //M_C0_du[__q1][__q1] = M_model->scalarProduct( __Xdu, __Xdu );
                        }//end of loop __m2
                    }//end of loop __q2
                }//end of loop __m1
            }//end of loop __q1

            if( this->worldComm().isMasterRank() )
                std::cout << "     o C0_du updated in " << ti.elapsed() << "s\n";
            ti.start();
        }

        LOG(INFO) << "[offlineResidual] Lambda_du, Gamma_du\n";

        for ( int elem=__N-added_elements; elem<__N; elem++ )
        {
            *__X=M_model->rBFunctionSpace()->dualBasisElement(elem);

            for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
            {
                for ( int __m1 = 0; __m1 < M_model->mMaxA(__q1); ++__m1 )
                {

                    if( boption(_prefix=M_prefix,_name="crb.use-symmetric-matrix") )
                        Atq1 = M_Aqm[__q1][__m1];
                    else
                        M_Aqm[__q1][__m1]->transpose( Atq1 );

                    Atq1->multVector(  __X, __W );
                    __W->close();
                    __W->scale( -1. );
                    //std::cout << "__W->norm=" << __W->l2Norm() << "\n";
                    M_model->l2solve( __Z1, __W );

                    for ( int __q2 = 0; __q2 < __QOutput; ++__q2 )
                    {
                        for ( int __m2 = 0; __m2 < M_model->mMaxF(M_output_index,__q2) ; ++__m2 )
                        {
                            M_Lambda_du[__q1][__m1][__q2][__m2].conservativeResize( __N );

                            *__Fdu = *M_Fqm[M_output_index][__q2][__m2];
                            __Fdu->close();
                            __Fdu->scale( -1.0 );
                            M_model->l2solve( __Z2, __Fdu );
                            M_Lambda_du[ __q1][ __m1][ __q2][ __m2]( elem ) = 2.0*M_model->scalarProduct( __Z2, __Z1 );
                            //VLOG(1) << "M_Lambda_pr[" << __q1 << "][" << __m1 << "][" << __q2 << "][" << __m2 << "][" << __j << "]=" << M_Lambda_pr[__q1][__m1][__q2][__m2][__j] << "\n";
                            //std::cout << "M_Lambda_pr[" << __q1 << "][" << __m1 << "][" << __q2 << "][" << __m2 << "][" << __j << "]=" << M_Lambda_pr[__q1][__m1][__q2][__m2][__j] << "\n";
                        } // m2
                    } // q2
                } // m1
            } // q1
        }//elem

        if( this->worldComm().isMasterRank() )
            std::cout << "     o Lambda_du updated in " << ti.elapsed() << "s\n";
        ti.start();

        for ( int elem=__N-added_elements; elem<__N; elem++ )
        {
            //int __j = __N-1;
            *__X=M_model->rBFunctionSpace()->dualBasisElement(elem);

            for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
            {
                for ( int __m1 = 0; __m1 < M_model->mMaxA(__q1); ++__m1 )
                {
                    if( boption(_prefix=M_prefix,_name="crb.use-symmetric-matrix") )
                        Atq1=M_Aqm[__q1][__m1];
                    else
                        M_Aqm[__q1][__m1]->transpose( Atq1 );

                    Atq1->multVector(  __X, __W );
                    __W->scale( -1. );
                    M_model->l2solve( __Z1, __W );

                    for ( int __l = 0; __l < ( int )__N; ++__l )
                    {
                        *__Y=M_model->rBFunctionSpace()->dualBasisElement(__l);

                        if( optimize )
                        {
                            //case we don't use EIM
                            for ( int __q2 = 0; __q2 < __q1; ++__q2 )
                            {
                                for ( int __m2 = 0; __m2 < M_model->mMaxA(__q2); ++__m2 )
                                {
                                    if( boption(_prefix=M_prefix,_name="crb.use-symmetric-matrix") )
                                        Atq2 = M_Aqm[__q2][__m2];
                                    else
                                        M_Aqm[__q2][__m2]->transpose( Atq2 );

                                    M_Gamma_du[__q1][__m1][__q2][__m2].conservativeResize( __N, __N );
                                    M_Gamma_du[__q2][__m2][__q1][__m1].conservativeResize( __N, __N );

                                    Atq2->multVector(  __Y, __W );
                                    __W->scale( -1. );
                                    M_model->l2solve( __Z2, __W );
                                    double prod = M_model->scalarProduct( __Z1, __Z2 );
                                    M_Gamma_du[ __q1][ __m1][ __q2][ __m2]( elem,__l ) = prod;
                                    M_Gamma_du[ __q2][ __m2][ __q1][ __m1]( __l,elem ) = prod;
                                } // m2
                            } // q2
                            //diagonal elements
                            Atq1->multVector(  __Y, __W );
                            __W->scale( -1. );
                            M_model->l2solve( __Z2, __W );
                            M_Gamma_du[__q1][__m1][__q1][__m1].conservativeResize( __N, __N );
                            double prod = M_model->scalarProduct( __Z1, __Z2 );
                            M_Gamma_du[ __q1][ __m1][ __q1][ __m1]( elem,__l ) = prod;
                        }// optimize
                        else
                        {
                            for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
                            {
                                for ( int __m2 = 0; __m2 < M_model->mMaxA(__q2); ++__m2 )
                                {
                                    if( boption(_prefix=M_prefix,_name="crb.use-symmetric-matrix") )
                                        Atq2 = M_Aqm[__q2][__m2];
                                    else
                                        M_Aqm[__q2][__m2]->transpose( Atq2 );

                                    M_Gamma_du[__q1][__m1][__q2][__m2].conservativeResize( __N, __N );

                                    Atq2->multVector(  __Y, __W );
                                    __W->scale( -1. );
                                    M_model->l2solve( __Z2, __W );
                                    M_Gamma_du[ __q1][ __m1][ __q2][ __m2]( elem,__l ) = M_model->scalarProduct( __Z1, __Z2 );
                                } // m2
                            }// q2
                        } // no optimize
                    }//__l
                }//m1
            }//q1
        }//__elem

        // update column __N-1
        for ( int __j = 0; __j < ( int )__N; ++__j )
        {
            *__X=M_model->rBFunctionSpace()->dualBasisElement(__j);

            for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
            {
                for ( int __m1 = 0; __m1 < M_model->mMaxA(__q1); ++__m1 )
                {
                    if( boption(_prefix=M_prefix,_name="crb.use-symmetric-matrix") )
                        Atq1=M_Aqm[__q1][__m1];
                    else
                        M_Aqm[__q1][__m1]->transpose( Atq1 );

                    Atq1->multVector(  __X, __W );
                    __W->scale( -1. );
                    M_model->l2solve( __Z1, __W );

                    //int __l = __N-1;
                    for ( int elem=__N-added_elements; elem<__N; elem++ )
                    {
                        *__Y=M_model->rBFunctionSpace()->dualBasisElement(elem);

                        if( optimize )
                        {
                            //case we don't use EIM
                            for ( int __q2 = 0; __q2 < __q1; ++__q2 )
                            {
                                for ( int __m2 = 0; __m2 < M_model->mMaxA(__q2); ++__m2 )
                                {
                                    if( boption(_prefix=M_prefix,_name="crb.use-symmetric-matrix") )
                                        Atq2 = M_Aqm[__q2][__m2];
                                    else
                                        M_Aqm[__q2][__m2]->transpose( Atq2 );

                                    Atq2->multVector(  __Y, __W );
                                    __W->scale( -1. );
                                    M_model->l2solve( __Z2, __W );
                                    double prod = M_model->scalarProduct( __Z1, __Z2 );
                                    M_Gamma_du[ __q1][ __m1][ __q2][ __m2]( __j, elem ) = prod;
                                    M_Gamma_du[ __q2][ __m2][ __q1][ __m1]( elem, __j ) = prod;
                                } // m2
                            } // q2
                            //diagonal elements
                            Atq1->multVector(  __Y, __W );
                            __W->scale( -1. );
                            M_model->l2solve( __Z2, __W );
                            double prod = M_model->scalarProduct( __Z1, __Z2 );
                            M_Gamma_du[ __q1][ __m1][ __q1][ __m1]( __j,elem ) = prod;
                        }// optimize
                        else
                        {
                            for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
                            {
                                for ( int __m2 = 0; __m2 < M_model->mMaxA(__q2); ++__m2 )
                                {
                                    if( boption(_prefix=M_prefix,_name="crb.use-symmetric-matrix") )
                                        Atq2 = M_Aqm[__q2][__m2];
                                    else
                                        M_Aqm[__q2][__m2]->transpose( Atq2 );

                                    Atq2->multVector(  __Y, __W );
                                    __W->scale( -1. );
                                    M_model->l2solve( __Z2, __W );
                                    M_Gamma_du[ __q1][ __m1][ __q2][ __m2]( __j,elem ) = M_model->scalarProduct( __Z1, __Z2 );
                                }//m2
                            }// q2
                        }//no optimize
                    }// elem
                }// m1
            }// q1
        } // __j

        if( this->worldComm().isMasterRank() )
            std::cout << "     o Gamma_du updated in " << ti.elapsed() << "s\n";
        ti.start();
        LOG(INFO) << "[offlineResidual] Done.\n";

    }//end of if (solve_dual_problem)

}
template<typename TruthModelType>
void
CRB<TruthModelType>::offlineResidualV1( int Ncur, mpl::bool_<false> , int number_of_added_elements )
{
    Feel::Timer ti;
    int __QLhs = M_model->Qa();
    int __QRhs = M_model->Ql( 0 );
    int __QOutput = M_model->Ql( M_output_index );
    int __N = Ncur;

    if( Environment::isMasterRank() )
        std::cout << "     o N=" << Ncur << " QLhs=" << __QLhs
                  << " QRhs=" << __QRhs << " Qoutput=" << __QOutput << "\n";

    vector_ptrtype __X( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __Fdu( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __W( M_backend->newVector( M_model->functionSpace() ) );
    namespace ublas = boost::numeric::ublas;
#if 0
    std::vector< std::vector<sparse_matrix_ptrtype> > Aqm,Mqm;
    std::vector< std::vector<vector_ptrtype> > MFqm;
    std::vector< std::vector<std::vector<vector_ptrtype> > > Fqm,Lqm;
#endif
    //boost::tie( Mqm, Aqm, Fqm ) = M_model->computeAffineDecomposition();
    boost::tie( M_Mqm, M_Aqm, M_Fqm ) = M_model->computeAffineDecomposition();

    if( Environment::isMasterRank() )
        std::cout << "     o initialize offlineResidual in " << ti.elapsed() << "s\n";

    bool optimize = boption(_prefix=M_prefix,_name="crb.optimize-offline-residual") ;

    ti.start();

    //
    //  Primal
    //
    LOG(INFO) << "[offlineResidual] Lambda_pr, Gamma_pr\n";

    if ( M_precomputeResidualPrF.empty()/*Ncur == M_Nm*/ )
    {
        // update M_precomputeResidualPrF
        M_precomputeResidualPrF.resize(__QRhs);
        for ( int __q2 = 0; __q2 < __QRhs; ++__q2 )
        {
            M_precomputeResidualPrF[__q2].resize( M_model->mMaxF(0, __q2) );
            for ( int __m2 = 0; __m2 < M_model->mMaxF(0, __q2); ++__m2 )
            {
                auto myvec = M_backend->newVector( M_model->functionSpace() );
                M_model->l2solve( myvec, M_Fqm[0][__q2][__m2] );
                M_precomputeResidualPrF[__q2][__m2] = myvec;
            }
        }

        // update M_C0_pr
        for ( int __q1 = 0; __q1 < __QRhs; ++__q1 )
        {
            for ( int __m1 = 0; __m1 < M_model->mMaxF(0,__q1); ++__m1 )
            {
                auto myZ1 = M_precomputeResidualPrF[__q1][__m1];
                for ( int __q2 = 0; __q2 < __QRhs; ++__q2 )
                {
                    for ( int __m2 = 0; __m2 < M_model->mMaxF(0,__q2); ++__m2 )
                    {
                        auto myZ2 = M_precomputeResidualPrF[__q2][__m2];
                        M_C0_pr[__q1][__m1][__q2][__m2] = M_model->scalarProduct( myZ1, myZ2 );
                    }
                }
            }
        }
        if( Environment::isMasterRank() )
            std::cout << "     o M_C0_pr updated\n";
    }

    // update M_precomputeResidualPrA
    int nElementsToPrecompute = (M_precomputeResidualPrA.empty())? __N : number_of_added_elements;
    M_precomputeResidualPrA.resize(__QLhs);
    for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
    {
        M_precomputeResidualPrA[__q1].resize( M_model->mMaxA(__q1) );
        for ( int __m1 = 0; __m1 < M_model->mMaxA(__q1); ++__m1 )
        {
            M_precomputeResidualPrA[__q1][__m1].resize(__N);
        }
    }
    for ( int __j=__N-nElementsToPrecompute/*__N-number_of_added_elements*/; __j<__N; __j++ )
    {
        *__X=M_model->rBFunctionSpace()->primalBasisElement(__j);
        for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
        {
            for ( int __m1 = 0; __m1 < M_model->mMaxA(__q1); ++__m1 )
            {
                M_Aqm[__q1][__m1]->multVector(  __X, __W );
                __W->scale( -1. );
                auto myvec = M_backend->newVector( M_model->functionSpace() );
                M_model->l2solve( myvec, __W );
                M_precomputeResidualPrA[__q1][__m1][__j] = myvec;
            }
        }
    }
    // update M_Gamma_pr
    for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
    {
        for ( int __m1 = 0; __m1 < M_model->mMaxA(__q1); ++__m1 )
        {
            for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
            {
                for ( int __m2 = 0; __m2 < M_model->mMaxA(__q2); ++__m2 )
                {
                    M_Gamma_pr[__q1][__m1][__q2][__m2].conservativeResize( __N, __N );

                    for ( int elem=__N-number_of_added_elements; elem<__N; elem++ )
                    {
                        auto myZ1 = M_precomputeResidualPrA[__q1][__m1][elem];
                        auto myZ3 = M_precomputeResidualPrA[__q2][__m2][elem];
                        for ( int __l = 0; __l < __N; ++__l )
                        {
                            auto myZ2 = M_precomputeResidualPrA[__q2][__m2][__l];
                            auto myZ4 = M_precomputeResidualPrA[__q1][__m1][__l];
                            double scalProd = M_model->scalarProduct( myZ1, myZ2 );
                            M_Gamma_pr[ __q1][ __m1][ __q2][ __m2]( elem,__l ) = scalProd;
                            double scalProd2 = M_model->scalarProduct( myZ4, myZ3 );
                            M_Gamma_pr[ __q1][ __m1][ __q2][ __m2]( __l,elem ) = scalProd2;
                        }
                    }
                }
            }
        }
    }

    // update M_Lambda_pr
    for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
    {
        for ( int __m1 = 0; __m1 < M_model->mMaxA(__q1) ; ++__m1 )
        {
            for ( int __q2 = 0; __q2 < __QRhs; ++__q2 )
            {
                for ( int __m2 = 0; __m2 < M_model->mMaxF(0, __q2); ++__m2 )
                {
                    M_Lambda_pr[__q1][__m1][__q2][__m2].conservativeResize( __N );

                    auto myZ2 = M_precomputeResidualPrF[__q2][__m2];
                    for ( int elem=__N-number_of_added_elements; elem<__N; elem++ )
                    {
                        auto myZ1 = M_precomputeResidualPrA[__q1][__m1][elem];
                        M_Lambda_pr[ __q1][ __m1][ __q2][ __m2]( elem ) = 2.0*M_model->scalarProduct( myZ1, myZ2 );
                    }
                }
            }
        }
    }

    if( Environment::isMasterRank() )
        std::cout << "     o Gamma_pr and Lambda_pr updated in " << ti.elapsed() << "s\n";

    //
    // Dual
    //
    if( M_solve_dual_problem )
    {
        ti.start();

        sparse_matrix_ptrtype Atq1;
        bool opIsSym = boption(_prefix=M_prefix,_name="crb.use-symmetric-matrix");
        if ( !opIsSym )
            Atq1 = M_model->newMatrix();

        // compute this only once
        if ( M_precomputeResidualDuF.empty()/*Ncur == M_Nm*/ )
        {
            // update M_precomputeResidualDuF
            M_precomputeResidualDuF.resize(__QOutput);
            for ( int __q2 = 0; __q2 < __QOutput; ++__q2 )
            {
                M_precomputeResidualDuF[__q2].resize( M_model->mMaxF(M_output_index,__q2) );
                for ( int __m2 = 0; __m2 < M_model->mMaxF(M_output_index,__q2); ++__m2 )
                {
                    *__Fdu = *M_Fqm[M_output_index][__q2][__m2];
                    __Fdu->close();
                    __Fdu->scale( -1.0 );
                    auto myvec = M_backend->newVector( M_model->functionSpace() );
                    M_model->l2solve( myvec, __Fdu );
                    M_precomputeResidualDuF[__q2][__m2] = myvec;
                }
            }

            // update M_C0_du
            for ( int __q1 = 0; __q1 < __QOutput; ++__q1 )
            {
                for ( int __m1 = 0; __m1 < M_model->mMaxF(M_output_index,__q1); ++__m1 )
                {
                    auto myZ1 = M_precomputeResidualDuF[__q1][__m1];
                    for ( int __q2 = 0; __q2 < __QOutput; ++__q2 )
                    {
                        for ( int __m2 = 0; __m2 < M_model->mMaxF(M_output_index,__q2); ++__m2 )
                        {
                            auto myZ2 = M_precomputeResidualDuF[__q2][__m2];
                            M_C0_du[__q1][__m1][__q2][__m2] = M_model->scalarProduct( myZ1, myZ2 );
                        }
                    }
                }
            }
            if( Environment::isMasterRank() )
                std::cout << "     o M_C0_du updated\n";
        }

        // update M_precomputeResidualDuA
        int nElementsToPrecomputeDual = (M_precomputeResidualDuA.empty())? __N : number_of_added_elements;
        M_precomputeResidualDuA.resize(__QLhs);
        for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
        {
            M_precomputeResidualDuA[__q1].resize( M_model->mMaxA(__q1) );
            for ( int __m1 = 0; __m1 < M_model->mMaxA(__q1); ++__m1 )
            {
                M_precomputeResidualDuA[__q1][__m1].resize(__N);
            }
        }
        for ( int __j=__N-nElementsToPrecomputeDual/*__N-number_of_added_elements*/; __j<__N; __j++ )
        {
            *__X=M_model->rBFunctionSpace()->dualBasisElement(__j);
            for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
            {
                for ( int __m1 = 0; __m1 < M_model->mMaxA(__q1); ++__m1 )
                {
                    if( opIsSym )
                        Atq1=M_Aqm[__q1][__m1];
                    else
                        M_Aqm[__q1][__m1]->transpose( Atq1 );

                    Atq1->multVector(  __X, __W );
                    __W->scale( -1. );
                    auto myvec = M_backend->newVector( M_model->functionSpace() );
                    M_model->l2solve( myvec, __W );
                    M_precomputeResidualDuA[__q1][__m1][__j] = myvec;
                }
            }
        }

        // update M_Gamma_du
        for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
        {
            for ( int __m1 = 0; __m1 < M_model->mMaxA(__q1); ++__m1 )
            {
                for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
                {
                    for ( int __m2 = 0; __m2 < M_model->mMaxA(__q2); ++__m2 )
                    {
                        M_Gamma_du[__q1][__m1][__q2][__m2].conservativeResize( __N, __N );

                        for ( int elem=__N-number_of_added_elements; elem<__N; elem++ )
                        {
                            auto myZ1 = M_precomputeResidualDuA[__q1][__m1][elem];
                            auto myZ3 = M_precomputeResidualDuA[__q2][__m2][elem];
                            for ( int __l = 0; __l < __N; ++__l )
                            {
                                auto myZ2 = M_precomputeResidualDuA[__q2][__m2][__l];
                                auto myZ4 = M_precomputeResidualDuA[__q1][__m1][__l];
                                double scalProd = M_model->scalarProduct( myZ1, myZ2 );
                                M_Gamma_du[ __q1][ __m1][ __q2][ __m2]( elem,__l ) = scalProd;
                                double scalProd2 = M_model->scalarProduct( myZ4, myZ3 );
                                M_Gamma_du[ __q1][ __m1][ __q2][ __m2]( __l,elem ) = scalProd2;
                            }
                        }
                    }
                }
            }
        }

        // update M_Lambda_du
        for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
        {
            for ( int __m1 = 0; __m1 < M_model->mMaxA(__q1) ; ++__m1 )
            {
                for ( int __q2 = 0; __q2 < __QOutput; ++__q2 )
                {
                    for ( int __m2 = 0; __m2 < M_model->mMaxF(M_output_index,__q2) ; ++__m2 )
                    {
                        M_Lambda_du[__q1][__m1][__q2][__m2].conservativeResize( __N );

                        auto myZ2 = M_precomputeResidualDuF[__q2][__m2];
                        for ( int elem=__N-number_of_added_elements; elem<__N; elem++ )
                        {
                            auto myZ1 = M_precomputeResidualDuA[__q1][__m1][elem];
                            M_Lambda_du[ __q1][ __m1][ __q2][ __m2]( elem ) = 2.0*M_model->scalarProduct( myZ1, myZ2 );
                        }
                    }
                }
            }
        }

        if( Environment::isMasterRank() )
            std::cout << "     o Gamma_du and Lambda_du updated in " << ti.elapsed() << "s\n";


    }//end of if (M_solve_dual_problem)
}



template<typename TruthModelType>
void
CRB<TruthModelType>::offlineResidualEim( int Ncur, mpl::bool_<false> , int number_of_added_elements )
{
    Feel::Timer ti;
    int __QLhs = M_model->Qa();
    int __QRhs = M_model->Ql( 0 );
    int __QOutput = M_model->Ql( M_output_index );
    int __N = Ncur;

    if( this->worldComm().isMasterRank() )
        std::cout << "     o offline residual eim "<<std::endl;

    vector_ptrtype __X( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __Y( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __Fdu( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __Z1(  M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __Z2(  M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __W(  M_backend->newVector( M_model->functionSpace() ) );
    namespace ublas = boost::numeric::ublas;

#if 0
    std::vector< std::vector<sparse_matrix_ptrtype> > Aqm,Mqm;
    std::vector< std::vector<vector_ptrtype> > MFqm;
    std::vector< std::vector<std::vector<vector_ptrtype> > > Fqm,Lqm;
#endif
    //boost::tie( Mqm, Aqm, Fqm ) = M_model->computeAffineDecomposition();
    boost::tie( M_Mqm, M_Aqm, M_Fqm ) = M_model->computeAffineDecomposition();

    __X->zero();
    __X->add( 1.0 );

    auto all_eim_interpolation_errors = M_model->eimInterpolationErrorEstimation();
    auto eim_interpolation_errors_A = all_eim_interpolation_errors.template get<1>() ;
    auto eim_interpolation_errors_F = all_eim_interpolation_errors.template get<2>() ;

    auto endA = eim_interpolation_errors_A.end();
    auto endF = eim_interpolation_errors_F[0].end();

    // Primal
    // no need to recompute this term each time
    if ( Ncur == M_Nm )
    {
        ti.start();
        for ( int __q1 = 0; __q1 < __QRhs; ++__q1 )
        {
            auto itq1 = eim_interpolation_errors_F[0].find(__q1);
            if( itq1 != endF )
            {
                //remember that in C++ index begins at 0
                //so to have max+1, we call [max]
                int Mmaxq1 = M_model->mMaxF(0,__q1);
                M_model->l2solve( __Z1, M_Fqm[0][__q1][Mmaxq1] );
                for ( int __q2 = 0; __q2 < __QRhs; ++__q2 )
                {
                    int Mmaxq2 = M_model->mMaxF(0,__q2);
                    auto itq2 = eim_interpolation_errors_F[0].find(__q2);
                    if( itq2 != endF )
                    {
                        M_model->l2solve( __Z2, M_Fqm[0][__q2][Mmaxq2] );
                        M_C0_pr_eim[__q1][__q2] = M_model->scalarProduct( __Z1, __Z2 );
                    }
                    else
                    {
                        M_C0_pr_eim[__q1][__q2] = 0;
                    }
                }//end of loop __q2
            }
            else
            {
                //no eim error associated to M_Fqm[0][__q1]
                for ( int __q2 = 0; __q2 < __QRhs; ++__q2 )
                {
                    M_C0_pr_eim[__q1][__q2] = 0;
                }
            }
        }//end of loop __q1

        if( this->worldComm().isMasterRank() )
            std::cout << "     o M_C0_pr_eim updated in " << ti.elapsed() << "s\n";

    }// Ncur==M_Nm

    ti.start();

    //
    //  Primal
    //
    LOG(INFO) << "[offlineResidual] Lambda_pr, Gamma_pr\n";

    for ( int elem=__N-number_of_added_elements; elem<__N; elem++ )
    {
        *__X=M_model->rBFunctionSpace()->primalBasisElement(elem);

        for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
        {
            auto itq1 = eim_interpolation_errors_A.find(__q1);
            if( itq1 != endA )
            {
                int Mmaxq1 = M_model->mMaxA(__q1);
                M_Aqm[__q1][Mmaxq1]->multVector(  __X, __W );
                __W->scale( -1. );
                M_model->l2solve( __Z1, __W );
                for ( int __q2 = 0; __q2 < __QRhs; ++__q2 )
                {
                    int Mmaxq2 = M_model->mMaxF(0,__q2);
                    auto itq2 = eim_interpolation_errors_F[0].find(__q2);
                    M_Lambda_pr_eim[__q1][__q2].conservativeResize( __N );
                    if( itq2 != endF )
                    {
                        M_model->l2solve( __Z2, M_Fqm[0][__q2][Mmaxq2] );
                        M_Lambda_pr_eim[ __q1][ __q2]( elem ) = 2.0*M_model->scalarProduct( __Z1, __Z2 );
                    }
                    else
                    {
                        M_Lambda_pr_eim[ __q1][ __q2]( elem ) = 0;
                    }
                }
            }
            else
            {
                for ( int __q2 = 0; __q2 < __QRhs; ++__q2 )
                {
                    M_Lambda_pr_eim[__q1][__q2].conservativeResize( __N );
                    M_Lambda_pr_eim[ __q1][ __q2]( elem ) = 0;
                }
            }
        }//end of __q1
    }//elem

    if( this->worldComm().isMasterRank() )
        std::cout << "     o Lambda_pr_eim updated in " << ti.elapsed() << "s\n";

    ti.start();

    for ( int elem=__N-number_of_added_elements; elem<__N; elem++ )
    {
        *__X=M_model->rBFunctionSpace()->primalBasisElement(elem);
        for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
        {
            auto itq1 = eim_interpolation_errors_A.find(__q1);
            if( itq1 != endA )
            {
                int Mmaxq1 = M_model->mMaxA(__q1);
                M_Aqm[__q1][Mmaxq1]->multVector(  __X, __W );
                __W->scale( -1. );
                M_model->l2solve( __Z1, __W );
                for ( int __l = 0; __l < ( int )__N; ++__l )
                {
                    *__Y=M_model->rBFunctionSpace()->primalBasisElement(__l);
                    for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
                    {
                        auto itq2 = eim_interpolation_errors_A.find(__q2);
                        if( itq2 != endA )
                        {
                            int Mmaxq2 = M_model->mMaxA(__q2);
                            M_Aqm[__q2][Mmaxq2]->multVector(  __Y, __W );
                            M_Gamma_pr_eim[__q1][__q2].conservativeResize( __N, __N );
                            __W->scale( -1. );
                            M_model->l2solve( __Z2, __W );
                            M_Gamma_pr_eim[ __q1][ __q2]( elem,__l ) = M_model->scalarProduct( __Z1, __Z2 );
                        }
                        else
                        {
                            M_Gamma_pr_eim[__q1][__q2].conservativeResize( __N, __N );
                            M_Gamma_pr_eim[__q1][__q2]( elem,__l ) = 0;
                        }
                    }
                }// l
            }// if eim error associated to q1
            else
            {
                for ( int __l = 0; __l < ( int )__N; ++__l )
                {
                    for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
                    {
                        M_Gamma_pr_eim[__q1][__q2].conservativeResize( __N, __N );
                        M_Gamma_pr_eim[__q1][__q2]( elem,__l ) = 0;
                    }
                }
            }//if no eim associated to q1
        } // q1
    } // elem

    for ( int __j = 0; __j < ( int )__N; ++__j )
    {
        *__X=M_model->rBFunctionSpace()->primalBasisElement(__j);
        for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
        {
            auto itq1 = eim_interpolation_errors_A.find(__q1);
            if( itq1 != endA )
            {
                int Mmaxq1 = M_model->mMaxA(__q1);
                M_Aqm[__q1][Mmaxq1]->multVector(  __X, __W );
                __W->scale( -1. );
                M_model->l2solve( __Z1, __W );

                //column N-1
                //int __l = __N-1;
                for ( int elem=__N-number_of_added_elements; elem<__N; elem++ )
                {
                    *__Y=M_model->rBFunctionSpace()->primalBasisElement(elem);

                    for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
                    {
                        auto itq2 = eim_interpolation_errors_A.find(__q2);
                        if( itq2 != endA )
                        {
                            int Mmaxq2 = M_model->mMaxA(__q2);
                            M_Aqm[__q2][Mmaxq2]->multVector(  __Y, __W );
                            __W->scale( -1. );
                            M_model->l2solve( __Z2, __W );
                            M_Gamma_pr_eim[ __q1][ __q2]( __j,elem ) = M_model->scalarProduct( __Z1, __Z2 );
                        }
                        else
                        {
                            M_Gamma_pr_eim[ __q1][ __q2]( __j,elem ) = 0;
                        }
                    } // q2
                }// end of loop elem
            }//if eim error associated to q1
            else
            {
                for ( int elem=__N-number_of_added_elements; elem<__N; elem++ )
                {
                    for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
                    {
                        M_Gamma_pr_eim[ __q1][ __q2]( __j,elem ) = 0;
                    }
                }
            }
        }// q1
    }// end of loop __j

    if( this->worldComm().isMasterRank() )
        std::cout << "     o Gamma_pr_eim updated in " << ti.elapsed() << "s\n";
    sparse_matrix_ptrtype Atq1 = M_model->newMatrix();
    sparse_matrix_ptrtype Atq2 = M_model->newMatrix();
    ti.start();

    //
    // Dual
    //
    // compute this only once
    if( M_solve_dual_problem )
    {
        auto endFo = eim_interpolation_errors_F[M_output_index].end();

        if ( Ncur == M_Nm )
        {
            LOG(INFO) << "[offlineResidual] Compute Dual residual data\n";
            LOG(INFO) << "[offlineResidual] C0_du\n";

            for ( int __q1 = 0; __q1 < __QOutput; ++__q1 )
            {
                auto itq1 = eim_interpolation_errors_F[M_output_index].find(__q1);
                if( itq1 != endFo )
                {
                    int Mmaxq1 = M_model->mMaxF(M_output_index,__q1);
                    *__Fdu = *M_Fqm[M_output_index][__q1][Mmaxq1];
                    __Fdu->close();
                    __Fdu->scale( -1.0 );
                    M_model->l2solve( __Z1, __Fdu );
                    for ( int __q2 = 0; __q2 < __QOutput; ++__q2 )
                    {
                        auto itq2 = eim_interpolation_errors_F[M_output_index].find(__q2);
                        if( itq2 != endFo )
                        {
                            int Mmaxq2 = M_model->mMaxF(M_output_index,__q2);
                            *__Fdu = *M_Fqm[M_output_index][__q2][Mmaxq2];
                            __Fdu->close();
                            __Fdu->scale( -1.0 );
                            M_model->l2solve( __Z2, __Fdu );
                            M_C0_du_eim[__q1][__q2] = M_model->scalarProduct( __Z1, __Z2 );
                        }
                        else
                        {
                            M_C0_du_eim[__q1][__q2] = 0;
                        }
                    }//q2
                }//if eim error associated to q1
                else
                {
                    for ( int __q2 = 0; __q2 < __QOutput; ++__q2 )
                    {
                        M_C0_du_eim[__q1][__q2] = 0;
                    }
                }
            }//end of loop __q1

            if( this->worldComm().isMasterRank() )
                std::cout << "     o C0_du_eim updated in " << ti.elapsed() << "s\n";
            ti.start();
        }

        LOG(INFO) << "[offlineResidual] Lambda_du, Gamma_du\n";

        for ( int elem=__N-number_of_added_elements; elem<__N; elem++ )
        {
            *__X=M_model->rBFunctionSpace()->dualBasisElement(elem);

            for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
            {
                auto itq1 = eim_interpolation_errors_A.find(__q1);
                if( itq1 != endA )
                {
                    int Mmaxq1 = M_model->mMaxA(__q1);
                    if( boption(_prefix=M_prefix,_name="crb.use-symmetric-matrix") )
                        Atq1 = M_Aqm[__q1][Mmaxq1];
                    else
                        M_Aqm[__q1][Mmaxq1]->transpose( Atq1 );
                    Atq1->multVector(  __X, __W );
                    __W->scale( -1. );
                    M_model->l2solve( __Z1, __W );

                    for ( int __q2 = 0; __q2 < __QOutput; ++__q2 )
                    {
                        M_Lambda_du_eim[__q1][__q2].conservativeResize( __N );
                        auto itq2 = eim_interpolation_errors_F[M_output_index].find(__q2);
                        if( itq2 != endFo )
                        {
                            int Mmaxq2 = M_model->mMaxF(M_output_index,__q2);
                            *__Fdu = *M_Fqm[M_output_index][__q2][Mmaxq2];
                            __Fdu->scale( -1.0 );
                            M_model->l2solve( __Z2, __Fdu );
                            M_Lambda_du_eim[__q1][__q2]( elem ) = 2.0*M_model->scalarProduct( __Z2, __Z1 );
                        }
                        else
                        {
                            M_Lambda_du_eim[__q1][__q2]( elem ) = 0;
                        }
                    } // q2
                }
                else
                {
                    for ( int __q2 = 0; __q2 < __QOutput; ++__q2 )
                    {
                        M_Lambda_du_eim[__q1][__q2].conservativeResize( __N );
                        M_Lambda_du_eim[__q1][__q2]( elem ) = 0;
                    }
                }
            } // q1
        }//elem

        if( this->worldComm().isMasterRank() )
            std::cout << "     o Lambda_du_eim updated in " << ti.elapsed() << "s\n";
        ti.start();

        for ( int elem=__N-number_of_added_elements; elem<__N; elem++ )
        {
            //int __j = __N-1;
            *__X=M_model->rBFunctionSpace()->dualBasisElement(elem);

            for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
            {
                auto itq1 = eim_interpolation_errors_A.find(__q1);
                if( itq1 != endA )
                {
                    int Mmaxq1 = M_model->mMaxA(__q1);

                    if( boption(_prefix=M_prefix,_name="crb.use-symmetric-matrix") )
                        Atq1=M_Aqm[__q1][Mmaxq1];
                    else
                        M_Aqm[__q1][Mmaxq1]->transpose( Atq1 );

                    Atq1->multVector(  __X, __W );
                    __W->scale( -1. );
                    M_model->l2solve( __Z1, __W );

                    for ( int __l = 0; __l < ( int )__N; ++__l )
                    {
                        *__Y=M_model->rBFunctionSpace()->dualBasisElement(__l);
                        for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
                        {
                            M_Gamma_du_eim[__q1][__q2].conservativeResize( __N, __N );

                            auto itq2 = eim_interpolation_errors_A.find(__q2);
                            if( itq2 != endA )
                            {
                                int Mmaxq2 = M_model->mMaxA(__q2);
                                if( boption(_prefix=M_prefix,_name="crb.use-symmetric-matrix") )
                                    Atq2 = M_Aqm[__q2][Mmaxq2];
                                else
                                    M_Aqm[__q2][Mmaxq2]->transpose( Atq2 );

                                Atq2->multVector(  __Y, __W );
                                __W->scale( -1. );
                                M_model->l2solve( __Z2, __W );
                                M_Gamma_du_eim[__q1][__q2]( elem,__l ) = M_model->scalarProduct( __Z1, __Z2 );
                            }
                            else
                            {
                                M_Gamma_du_eim[__q1][__q2]( elem,__l ) = 0;
                            }
                        }// q2
                    }//__l
                }
                else
                {
                    for ( int __l = 0; __l < ( int )__N; ++__l )
                    {
                        for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
                        {
                            M_Gamma_du_eim[__q1][__q2].conservativeResize( __N, __N );
                            M_Gamma_du_eim[__q1][__q2]( elem,__l ) = 0;
                        }
                    }
                }
            }//q1
        }//__elem

        // update column __N-1
        for ( int __j = 0; __j < ( int )__N; ++__j )
        {
            *__X=M_model->rBFunctionSpace()->dualBasisElement(__j);

            for ( int __q1 = 0; __q1 < __QLhs; ++__q1 )
            {

                auto itq1 = eim_interpolation_errors_A.find(__q1);
                if( itq1 != endA )
                {
                    int Mmaxq1 = M_model->mMaxA(__q1);

                    if( boption(_prefix=M_prefix,_name="crb.use-symmetric-matrix") )
                        Atq1=M_Aqm[__q1][Mmaxq1];
                    else
                        M_Aqm[__q1][Mmaxq1]->transpose( Atq1 );

                    Atq1->multVector(  __X, __W );
                    __W->scale( -1. );
                    M_model->l2solve( __Z1, __W );

                    for ( int elem=__N-number_of_added_elements; elem<__N; elem++ )
                    {
                        *__Y=M_model->rBFunctionSpace()->dualBasisElement(elem);

                        for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
                        {
                            auto itq2 = eim_interpolation_errors_A.find(__q2);
                            if( itq2 != endA )
                            {
                                int Mmaxq2 = M_model->mMaxA(__q2);

                                if( boption(_prefix=M_prefix,_name="crb.use-symmetric-matrix") )
                                    Atq2 = M_Aqm[__q2][Mmaxq2];
                                else
                                    M_Aqm[__q2][Mmaxq2]->transpose( Atq2 );

                                Atq2->multVector(  __Y, __W );
                                __W->scale( -1. );
                                M_model->l2solve( __Z2, __W );
                                M_Gamma_du_eim[ __q1][ __q2]( __j,elem ) = M_model->scalarProduct( __Z1, __Z2 );
                            }
                            else
                            {
                                M_Gamma_du_eim[ __q1][ __q2]( __j,elem ) = 0;
                            }
                        }// q2
                    }// elem
                }
                else
                {
                    for ( int elem=__N-number_of_added_elements; elem<__N; elem++ )
                    {
                        for ( int __q2 = 0; __q2 < __QLhs; ++__q2 )
                        {
                            M_Gamma_du_eim[ __q1][ __q2]( __j,elem ) = 0;
                        }
                    }
                }
            }// q1
        } // __j

        if( this->worldComm().isMasterRank() )
            std::cout << "     o Gamma_du_eim updated in " << ti.elapsed() << "s\n";
        ti.start();
        LOG(INFO) << "[offlineResidual eim] Done.\n";

    }//end of if (M_solve_dual_problem)
}

template<typename TruthModelType>
void
CRB<TruthModelType>::printErrorsDuringRbConstruction( void )
{

    std::ofstream conv;
    std::string file_name = "crb-offline-error.dat";
    typedef convergence_type::left_map::const_iterator iterator;

    if( this->worldComm().isMasterRank() )
    {
        conv.open(file_name, std::ios::app);
        conv << "NbBasis" << "\t" << "output" << "\t" << "primal" << "\t" << "dual\n";

        for(iterator it = M_rbconv.left.begin(); it != M_rbconv.left.end(); ++it)
            conv<<it->first<<"\t"<<it->second.template get<0>()<<"\t"<<it->second.template get<1>()<<"\t"<<it->second.template get<2>()<<"\n";
        //for(iterator it = M_rbconv.left.begin(); it != M_rbconv.left.end(); ++it)
        //    LOG(INFO)<<"N : "<<it->first<<"  -  delta_du : "<<it->second.template get<2>()<<"\n";
    }
    conv.close();

}



template<typename TruthModelType>
void
CRB<TruthModelType>::printMuSelection( void )
{
    LOG(INFO)<<" List of parameter selectionned during the offline algorithm \n";
    for(int k=0;k<M_WNmu->size();k++)
    {
        std::cout<<" mu"<<k<<"= [ ";
        LOG(INFO)<<" mu"<<k<<"= [ ";
        parameter_type const& _mu = M_WNmu->at( k );
        for( int i=0; i<_mu.size()-1; i++ )
        {
            LOG(INFO)<<_mu(i)<<" , ";
            std::cout<<_mu(i)<<" , ";
        }
        LOG(INFO)<<_mu( _mu.size()-1 )<<" ] \n";
        std::cout<<_mu( _mu.size()-1 )<<" ] "<<std::endl;
    }
}

template<typename TruthModelType>
void
CRB<TruthModelType>::printRbPicardIterations()
{
    if ( this->worldComm().isMasterRank() )
    {
        std::ofstream file;
        std::string filename = "RB-offline_picard_summary.dat";
        file.open(filename, std::ios::out | std::ios::app);

        std::string name = fs::current_path().string() + "/" + filename;
        fs::path pathfile( name.c_str() );
        if( fs::is_empty(pathfile) )
            file << "NbBasis" << "\t " << "Nb_iter" << "\t" << "Increment \n";
        file << M_N << "\t" << offline_iterations_summary.first << "\t" << offline_iterations_summary.second << "\n";
        file.close();
    }
}

template<typename TruthModelType>
void
CRB<TruthModelType>::printOnlineRbPicardIterations(parameter_type const& mu) const
{
    if ( this->worldComm().isMasterRank() )
    {
        std::ofstream file;
        std::string filename = "RB-online_picard_summary.dat";
        file.open(filename, std::ios::out | std::ios::app);

        double rb_online_mean_iterations = 0;
        int rb_online_min_terations = 0;
        int rb_online_max_terations = 0;
        double rb_online_max_increments = 0;

        std::string name = fs::current_path().string() + "/" + filename;
        fs::path pathfile( name.c_str() );
        auto size = mu.size();
        if( fs::is_empty(pathfile) )
        {
            file << "mu" << "\t " << "Nb_iter" << "\t" << "Increment \n";
        }
        // file << M_N << "\t" << offline_iterations_summary.first << "\t" << offline_iterations_summary.second << "\n";

        file << "[ ";
        for ( int i=0; i<size-1; i++ )
            file << mu[i] <<" , ";
        file << mu[size-1]<<" ]" << "\t ";
        file << online_iterations_summary.first << "\t" << online_iterations_summary.second << "\n";
        file.close();
    }
}

template<typename TruthModelType>
typename CRB<TruthModelType>::element_type
CRB<TruthModelType>::runWithExpansion( parameter_type const& mu , int N , int time_index )
{
    int Nwn;

    if( N > 0 )
        Nwn = N;
    else
        Nwn = M_N;

    std::vector<vectorN_type> uN;
    std::vector<vectorN_type> uNdu;
    std::vector<vectorN_type> uNold;
    std::vector<vectorN_type> uNduold;

    auto o = lb( Nwn, mu, uN, uNdu , uNold, uNduold );
    int size = uN.size();

    element_type ucrb;
    if( time_index == -1 )
    {
        ucrb = expansion( uN[size-1], Nwn );
        //ucrb = Feel::expansion( M_model->rBFunctionSpace()->primalRB(), uN[size-1] , Nwn);
    }
    else
    {
        CHECK( time_index < size )<<" call crb::expansion with a wrong value of time index : "<<time_index<<" or size of uN vector is only "<<size;
        ucrb = expansion( uN[time_index], Nwn );
        //ucrb = Feel::expansion( M_model->rBFunctionSpace()->primalRB(), uN[time_index] , Nwn);
    }
    return ucrb;
}


template<typename TruthModelType>
typename CRB<TruthModelType>::element_type
CRB<TruthModelType>::expansion( vectorN_type const& u, int N, bool dual ) const
{
    auto WN = dual ? M_model->rBFunctionSpace()->dualRB() : M_model->rBFunctionSpace()->primalRB();
    int Nwn;

    if( N > 0 )
        Nwn = N;
    else
        Nwn = M_N;

    //FEELPP_ASSERT( M_WN.size() == u.size() )( M_WN.size() )( u.size() ).error( "invalid expansion size");
    FEELPP_ASSERT( Nwn <= WN.size() )( Nwn )( WN.size() ).error( "invalid expansion size ( Nwn and WN ) ");
    FEELPP_ASSERT( Nwn <= u.size() )( Nwn )( WN.size() ).error( "invalid expansion size ( Nwn and u ) ");
    //FEELPP_ASSERT( Nwn == u.size() )( Nwn )( u.size() ).error( "invalid expansion size");
    return Feel::expansion( WN, u, Nwn );
}


template<typename TruthModelType>
//typename boost::tuple<std::vector<double>,double, typename CRB<TruthModelType>::solutions_tuple, typename CRB<TruthModelType>::matrix_info_tuple,
//double, double, typename CRB<TruthModelType>::upper_bounds_tuple >
CRBResults
CRB<TruthModelType>::run( parameter_type const& mu, vectorN_type & time, double eps , int N, bool print_rb_matrix)
{
    //int Nwn = M_N;
    //int Nwn_max = ioption(_prefix=M_prefix,_name="crb.dimension-max");
#if 0
    if (  M_error_type!=CRB_EMPIRICAL )
    {
        auto lo = M_rbconv.right.range( boost::bimaps::unbounded, boost::bimaps::_key <= eps );

        for ( auto it = lo.first; it != lo.second; ++it )
        {
            std::cout << "rbconv[" << it->first <<"]=" << it->second << "\n";
        }

        auto it = M_rbconv.project_left( lo.first );
        Nwn = it->first;
        std::cout << "Nwn = "<< Nwn << " error = "<< it->second.template get<0>() << " eps=" << eps << "\n";
    }

    if ( boption(_prefix=M_prefix,_name="crb.check.residual") )
    {
        std::vector< std::vector<double> > primal_residual_coefficients = error_estimation.template get<1>();
        std::vector< std::vector<double> > dual_residual_coefficients = error_estimation.template get<2>();
        std::vector<element_ptrtype> Un,Unold,Undu,Unduold;
        buildFunctionFromRbCoefficients(Nwn, uN, M_WN, Un );
        buildFunctionFromRbCoefficients(Nwn, uNold, M_WN, Unold );
        buildFunctionFromRbCoefficients(Nwn, uNdu, M_WNdu, Undu );
        buildFunctionFromRbCoefficients(Nwn, uNduold, M_WNdu, Unduold );
        compareResidualsForTransientProblems(Nwn, mu , Un, Unold, Undu, Unduold, primal_residual_coefficients, dual_residual_coefficients );
    }
#endif
    M_model->countAffineDecompositionTerms();

    std::vector<vectorN_type> uN;
    std::vector<vectorN_type> uNdu;
    std::vector<vectorN_type> uNold;
    std::vector<vectorN_type> uNduold;

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
    Feel::Timer t1;
    auto o = lb( Nwn, mu, uN, uNdu , uNold, uNduold , print_rb_matrix);
    double time_prediction=t1.elapsed();
    auto output_vector=o.template get<0>();
    double output_vector_size=output_vector.size();
    double output = output_vector[output_vector_size-1];
    t1.start();
    auto error_estimation = delta( Nwn, mu, uN, uNdu , uNold, uNduold );
    double time_error_estimation=t1.elapsed();
    auto vector_output_upper_bound = error_estimation.template get<0>();
    double output_upper_bound = vector_output_upper_bound[0];
    auto matrix_info = o.template get<1>();
    auto primal_coefficients = error_estimation.template get<1>();
    auto dual_coefficients = error_estimation.template get<2>();

    //compute dual norm of primal/dual residual at final time
    int nb_dt = primal_coefficients.size();
    int final_time_index = nb_dt-1;
    double primal_residual_norm = 0;
    double dual_residual_norm = 0;

    if ( M_error_type == CRB_RESIDUAL || M_error_type == CRB_RESIDUAL_SCM )
    {
        int nb_coeff = primal_coefficients[final_time_index].size();
        for(int i=0 ; i<nb_coeff ; i++)
            primal_residual_norm += primal_coefficients[final_time_index][i] ;

        if( M_solve_dual_problem )
        {
            if ( M_model->isSteady() )
                dual_residual_norm =  math::abs( dual_coefficients[0][0]+dual_coefficients[0][1]+dual_coefficients[0][2] ) ;
            else
                dual_residual_norm =  math::abs( dual_coefficients[0][2]+dual_coefficients[0][4]+dual_coefficients[0][5] ) ;
        }

        primal_residual_norm = math::sqrt( math::abs(primal_residual_norm) );
        dual_residual_norm = math::sqrt( math::abs(dual_residual_norm) );
    }

    double delta_pr = error_estimation.template get<3>();
    double delta_du = error_estimation.template get<4>();

    time.resize(2);
    time(0)=time_prediction;
    time(1)=time_error_estimation;

    int size = uN.size();

    auto upper_bounds = boost::make_tuple(vector_output_upper_bound , delta_pr, delta_du , primal_coefficients , dual_coefficients );
    auto solutions = boost::make_tuple( uN , uNdu, uNold, uNduold);

    CRBResults r(boost::make_tuple( output_vector , Nwn , solutions, matrix_info , primal_residual_norm , dual_residual_norm, upper_bounds ));
    r.setParameter( mu );
    return r;
}


template<typename TruthModelType>
void
CRB<TruthModelType>::run( const double * X, unsigned long N, double * Y, unsigned long P )
{
    return run( X, N, Y, P, mpl::bool_<model_type::is_time_dependent>() );
}



template<typename TruthModelType>
void
CRB<TruthModelType>::run( const double * X, unsigned long N, double * Y, unsigned long P, mpl::bool_<true> )
{
    //std::cout<<"N = "<<N<<" et P = "<<P<<std::endl;

    for ( unsigned long p= 0; p < N-5; ++p ) std::cout<<"mu["<<p<<"] = "<<X[p]<<std::endl;


    parameter_type mu( M_Dmu );

    // the last parameter is the max error
    for ( unsigned long p= 0; p < N-5; ++p )
        mu( p ) = X[p];


    //double meshSize  = X[N-4];
    //M_model->setMeshSize(meshSize);
    // setOutputIndex( ( int )X[N-3] );
    // std::cout<<"output index = "<<X[N-3]<<std::endl;
    // int Nwn =  X[N-2];
    // std::cout<<"Nwn : "<<X[N-2]<<std::endl;
    // size_type maxerror = X[N-1];
    // std::cout<<"maxerror = "<<X[N-1]<<std::endl;
    //CRBErrorType errorType =(CRBErrorType)X[N-1];
    //setCRBErrorType(errorType);

    setOutputIndex( ( int )X[N-5] );
    //std::cout<<"output_index = "<<X[N-5]<<std::endl;
    int Nwn =  X[N-4];
    //std::cout<<" Nwn = "<<Nwn<<std::endl;
    int maxerror = X[N-3];
    //std::cout<<" maxerror = "<<maxerror<<std::endl;
    CRBErrorType errorType =( CRBErrorType )X[N-2];
    //std::cout<<"errorType = "<<X[N-2]<<std::endl;
    setCRBErrorType( errorType );
    M_compute_variance = X[N-1];
    //std::cout<<"M_compute_variance = "<<M_compute_variance<<std::endl;


#if 0

    if (  M_error_type!=CRB_EMPIRICAL )
    {
        auto lo = M_rbconv.right.range( boost::bimaps::unbounded,boost::bimaps::_key <= maxerror );

        for ( auto it = lo.first; it != lo.second; ++it )
        {
            std::cout << "rbconv[" << it->first <<"]=" << it->second << "\n";
        }

        auto it = M_rbconv.project_left( lo.first );
        Nwn = it->first;
    }

#endif

    std::vector<vectorN_type> uN;
    std::vector<vectorN_type> uNdu;
    std::vector<vectorN_type> uNold ;
    std::vector<vectorN_type> uNduold;

    FEELPP_ASSERT( P == 2 )( P ).warn( "invalid number of outputs" );
    auto o = lb( Nwn, mu, uN, uNdu , uNold, uNduold );
    auto e = delta( Nwn, mu, uN, uNdu , uNold, uNduold );
    auto output_vector=o.template get<0>();
    double output_vector_size=output_vector.size();
    double output = output_vector[output_vector_size-1];
    auto vector_err = e.template get<0>();
    int size=vector_err.size();
    double error = vector_err[size-1];
    Y[0]  = output;
    Y[1]  = error;
}


template<typename TruthModelType>
void
CRB<TruthModelType>::run( const double * X, unsigned long N, double * Y, unsigned long P, mpl::bool_<false> )
{

    parameter_type mu( M_Dmu );

    for ( unsigned long p= 0; p < N-5; ++p )
    {
        mu( p ) = X[p];
        std::cout << "mu( " << p << " ) = " << mu( p ) << std::endl;
    }

    // std::cout<<" list of parameters : [";
    // for ( unsigned long i=0; i<N-1; i++ ) std::cout<<X[i]<<" , ";
    // std::cout<<X[N-1]<<" ] "<<std::endl;


    setOutputIndex( ( int )X[N-5] );
    //std::cout<<"output_index = "<<X[N-5]<<std::endl;
    int Nwn =  X[N-4];
    //std::cout<<" Nwn = "<<Nwn<<std::endl;
    int maxerror = X[N-3];
    //std::cout<<" maxerror = "<<maxerror<<std::endl;
    CRBErrorType errorType =( CRBErrorType )X[N-2];
    //std::cout<<"errorType = "<<X[N-2]<<std::endl;
    setCRBErrorType( errorType );
    M_compute_variance = X[N-1];
    //std::cout<<"M_compute_variance = "<<M_compute_variance<<std::endl;

#if 0
    auto lo = M_rbconv.right.range( boost::bimaps::unbounded,boost::bimaps::_key <= maxerror );


    for ( auto it = lo.first; it != lo.second; ++it )
    {
        std::cout << "rbconv[" << it->first <<"]=" << it->second << "\n";
    }

    auto it = M_rbconv.project_left( lo.first );
    Nwn = it->first;
    //std::cout << "Nwn = "<< Nwn << " error = "<< it->second << " maxerror=" << maxerror << " ErrorType = "<<errorType <<"\n";
#endif
    //int Nwn = M_WN.size();
    std::vector<vectorN_type> uN;
    std::vector<vectorN_type> uNdu;
    std::vector<vectorN_type> uNold ;
    std::vector<vectorN_type> uNduold;

    FEELPP_ASSERT( P == 2 )( P ).warn( "invalid number of outputs" );
    auto o = lb( Nwn, mu, uN, uNdu , uNold, uNduold );
    auto e = delta( Nwn, mu, uN, uNdu , uNold, uNduold );
    auto output_vector=o.template get<0>();
    double output_vector_size=output_vector.size();
    double output = output_vector[output_vector_size-1];
    auto vector_err = e.template get<0>();
    int size=vector_err.size();
    double error = vector_err[size-1];
    Y[0]  = output;
    Y[1]  = error;
}






template<typename TruthModelType>
void
CRB<TruthModelType>::projectionOnPodSpace( const element_type & u , element_ptrtype& projection, const std::string& name_of_space )
{

    projection->zero();

    if ( name_of_space=="dual" )
    {
        if ( M_orthonormalize_dual )
            //in this case we can simplify because elements of reduced basis are orthonormalized
        {
            for( auto const& du : M_model->rBFunctionSpace()->dualRB() )
            {
                element_type e = unwrap_ptr( du ).functionSpace()->element();
                e = unwrap_ptr( du );
                double k =  M_model->scalarProduct( u, e );
                e.scale( k );
                projection->add( 1 , e );
            }
        }//end of if orthonormalize_dual

        else
        {
            matrixN_type MN ( ( int )M_N, ( int )M_N ) ;
            vectorN_type FN ( ( int )M_N );

            for ( size_type i=0; i<M_N; i++ )
            {
                for ( size_type j=0; j<i; j++ )
                {
                    MN( i,j ) = M_model->scalarProduct( M_model->rBFunctionSpace()->dualBasisElement(j) , M_model->rBFunctionSpace()->dualBasisElement(i) );
                    MN( j,i ) = MN( i,j );
                }

                MN( i,i ) = M_model->scalarProduct( M_model->rBFunctionSpace()->dualBasisElement(i) , M_model->rBFunctionSpace()->dualBasisElement(i) );
                FN( i ) = M_model->scalarProduct( u,M_model->rBFunctionSpace()->dualBasisElement(i) );
            }

            vectorN_type projectionN ( ( int ) M_N );
            projectionN = MN.lu().solve( FN );
            int index=0;
            for( auto const& du : M_model->rBFunctionSpace()->dualRB() )
            {
                element_type e = unwrap_ptr( du ).functionSpace()->element();
                e = unwrap_ptr( du );
                double k =  projectionN( index );
                e.scale( k );
                projection->add( 1 , e );
                index++;
            }
        }//end of if( ! orthonormalize_dual )
    }//end of projection on dual POD space

    else
    {
        if ( M_orthonormalize_primal )
        {
            for( auto const& pr : M_model->rBFunctionSpace()->primalRB() )
            {
                auto e = unwrap_ptr( pr ).functionSpace()->element();
                e = unwrap_ptr( pr );
                double k =  M_model->scalarProduct( u, e );
                e.scale( k );
                projection->add( 1 , e );
            }
        }//end of if orthonormalize_primal

        else
        {
            matrixN_type MN ( ( int )M_N, ( int )M_N ) ;
            vectorN_type FN ( ( int )M_N );

            for ( size_type i=0; i<M_N; i++ )
            {
                for ( size_type j=0; j<i; j++ )
                {
                    MN( i,j ) = M_model->scalarProduct( M_model->rBFunctionSpace()->primalBasisElement(j) , M_model->rBFunctionSpace()->primalBasisElement(i) );
                    MN( j,i ) = MN( i,j );
                }

                MN( i,i ) = M_model->scalarProduct( M_model->rBFunctionSpace()->primalBasisElement(i) , M_model->rBFunctionSpace()->primalBasisElement(i) );
                FN( i ) = M_model->scalarProduct( u, M_model->rBFunctionSpace()->primalBasisElement(i) );
            }

            vectorN_type projectionN ( ( int ) M_N );
            projectionN = MN.lu().solve( FN );
            int index=0;
            for( auto const& pr : M_model->rBFunctionSpace()->primalRB() )
            {
                element_type e = unwrap_ptr( pr ).functionSpace()->element();
                e = unwrap_ptr( pr );
                double k =  projectionN( index );
                e.scale( k );
                projection->add( 1 , e );
                index++;
            }
        }//end of if( ! orthonormalize_primal )
    }//end of projection on primal POD space

}

template<typename TruthModelType>
void
CRB<TruthModelType>::updateApeeBasisConstruction( int N , std::vector< parameter_type > const & new_primal_parameters, std::vector< parameter_type > const & new_dual_parameters ) const
{
    int index = M_primal_V.size(); // current size
    int primal_size=new_primal_parameters.size();
    int new_size = index + primal_size;
    for(int i=0; i<primal_size; i++)
    {
        auto mu = new_primal_parameters[i];
        computePrimalApeeBasis( mu );
    }

    int dual_size=new_dual_parameters.size();
    new_size = index + dual_size;
    for(int i=0; i<dual_size; i++)
    {
        auto mu = new_dual_parameters[i];
        computeDualApeeBasis( mu );
    }

}

template<typename TruthModelType>
void
CRB<TruthModelType>::updateApeeOfflineResidualComputation( int N , std::vector< parameter_type > const & new_primal_parameters, std::vector< parameter_type > const & new_dual_parameters )
{
    int index = M_primal_V.size(); // current size
    int primal_size=new_primal_parameters.size();
    int new_size = index + primal_size;
    M_primal_V.conservativeResize( new_size );
    for(int i=0; i<primal_size; i++)
    {
        auto mu = new_primal_parameters[i];
        auto u = offlineFixedPointPrimal( mu );
        M_primal_V(index) = computeSquareDualNormOfPrimalResidual( mu , u );
    }
    index = M_dual_V.size(); // current size

    element_ptrtype dual_initial_field( new element_type( M_model->functionSpace() ) );
    int dual_size=new_dual_parameters.size();
    new_size = index + dual_size;
    M_dual_V.conservativeResize( new_size );
    for(int i=0; i<dual_size; i++)
    {
        auto mu = new_dual_parameters[i];
        auto udu = offlineFixedPointDual( mu, dual_initial_field);
        M_dual_V(index) = computeSquareDualNormOfDualResidual( mu , udu );
    }

}


template<typename TruthModelType>
double
CRB<TruthModelType>::computeSquareDualNormOfPrimalResidual( parameter_type const& mu, element_type const & u )
{
    sparse_matrix_ptrtype A;
    std::vector<vector_ptrtype> F;

    bool is_linear=M_model->isLinear();
    if( is_linear )
        boost::tie( boost::tuples::ignore, A, F ) = M_model->update( mu );
    else
        boost::tie( boost::tuples::ignore, A, F ) = M_model->update( mu, u );

    vector_ptrtype Aun( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype Un( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype Frhs( M_backend->newVector( M_model->functionSpace() ) );
    *Un = u;
    A->multVector( Un, Aun );
    Aun->close();
    Aun->scale( -1 );
    *Frhs = *F[0];

    auto primal_residual = Frhs ;
    primal_residual->add( *Aun );

    vector_ptrtype __e_pr(  M_backend->newVector( M_model->functionSpace() ) );
    M_model->l2solve( __e_pr, primal_residual );

    double dual_norm_primal_residual = math::sqrt( M_model->scalarProduct( __e_pr,__e_pr ) );

    return dual_norm_primal_residual ;
}

template<typename TruthModelType>
double
CRB<TruthModelType>::computeSquareDualNormOfDualResidual( parameter_type const& mu, element_type const & udu )
{
    sparse_matrix_ptrtype A,At;
    std::vector<vector_ptrtype> F;

    boost::tie( boost::tuples::ignore, A, F ) = M_model->update( mu );

    At = M_model->newMatrix();
    if( boption(_prefix=M_prefix,_name="crb.use-symmetric-matrix") )
        At = A;
    else
        A->transpose( At );

    vector_ptrtype Atun( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype Undu( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype Lrhs( M_backend->newVector( M_model->functionSpace() ) );
    *Undu = udu;
    At->multVector( Undu, Atun );
    Atun->close();
    Atun->scale( -1 );
    *Lrhs = *F[M_output_index];

    auto dual_residual = Lrhs ;
    dual_residual->add( *Atun );

    vector_ptrtype __e_du(  M_backend->newVector( M_model->functionSpace() ) );
    M_model->l2solve( __e_du, dual_residual );

    double dual_norm_dual_residual = math::sqrt( M_model->scalarProduct( __e_du,__e_du ) );

    return dual_norm_dual_residual ;
}

template<typename TruthModelType>
void
CRB<TruthModelType>::assemblePrimalAndDualApeeBasisMatrices( int N ) const
{
    int Qa = M_model->Qa();
    int Qf = M_model->Ql( 0 );
    int Ql = M_model->Ql( M_output_index );
    int primal_size = Qf*Qf + N*Qa*Qf + N*N*Qa*Qa;
    int dual_size = Ql*Ql + N*Qa*Ql + N*N*Qa*Qa;

    LOG( INFO ) << "[assemblePrimalAndDualApeeBasisMatrices] primal matrix : "<<primal_size<<" * "<<primal_size<<" and dual matrix : "<<dual_size<<" * "<<dual_size;
    int primal_rows=M_primal_T.rows();
    int primal_cols=M_primal_T.cols();
    FEELPP_ASSERT( primal_rows == primal_cols )( primal_rows )( primal_cols ).error( "invalid size of matrix M_primal_T (should be a square matrix)");

    int number_of_elements_to_add = primal_size - primal_rows;
    FEELPP_ASSERT( number_of_elements_to_add  > 0 )( number_of_elements_to_add )(primal_size)(primal_rows).error( "assemblePrimalAndDualApeeBasisMatrices was called but there is no elements to add in M_primal_T ");

    M_primal_T.conservativeResize( primal_size , primal_size );
    for(int i=primal_rows; i<primal_size; i++)
    {
        M_primal_T.col( i ) = M_primal_apee_basis[ i ];
    }

    int dual_rows=M_dual_T.rows();
    int dual_cols=M_dual_T.cols();
    FEELPP_ASSERT( dual_rows == dual_cols )( dual_rows )( dual_cols ).error( "invalid size of matrix M_dual_T (should be a square matrix)");

    number_of_elements_to_add = dual_size - dual_rows;
    FEELPP_ASSERT( number_of_elements_to_add  > 0 )( number_of_elements_to_add )(dual_size)(dual_rows).error( "assemblePrimalAndDualApeeBasisMatrices was called but there is no elements to add in M_dual_T ");

    M_dual_T.conservativeResize( dual_size , dual_size );
    for(int i=dual_rows; i<dual_size; i++)
    {
        M_dual_T.col( i ) = M_dual_apee_basis[ i ];
    }

}

template<typename TruthModelType>
void
CRB<TruthModelType>::selectPrimalApeeParameters( int N, std::vector< parameter_type > & new_parameters )
{
    int Qa = M_model->Qa();
    int Qf = M_model->Ql( 0 );
    int primal_size = Qf*Qf + N*Qa*Qf + N*N*Qa*Qa;
    int current_size = M_primal_apee_mu->size();
    //number of parameters to add :
    int elem_to_add = primal_size - current_size;
    LOG( INFO ) << "[selectPrimalApeeParameters] elem_to_add : "<<elem_to_add;
    new_parameters.resize( elem_to_add );
    parameter_type mu( M_Dmu );
    size_type index=0;
    for(int i=0; i<elem_to_add; i++)
    {
        bool already_exist;
        do
        {
            //initialization
            already_exist=false;
            //pick randomly an element in parameter space
            mu = M_Dmu->element();
            //make sure that the new mu is not already is M_primal_apee_mu or in M_WNmu
            for( auto _mu : *M_primal_apee_mu )
            {
                if( mu == _mu )
                    already_exist=true;
            }
            for( auto _mu : *M_WNmu )
            {
                if( mu == _mu )
                    already_exist=true;
            }
        }
        while( already_exist );
        M_primal_apee_mu->push_back( mu , index );
        new_parameters[i]=mu;
        //M_primal_apee_mu_complement = M_primal_apee_mu->complement();
    }
}

template<typename TruthModelType>
void
CRB<TruthModelType>::selectDualApeeParameters( int N, std::vector< parameter_type > & new_parameters )
{
    int Qa = M_model->Qa();
    int Ql = M_model->Ql( M_output_index );
    int dual_size = Ql*Ql + N*Qa*Ql + N*N*Qa*Qa;
    int current_size = M_dual_apee_mu->size();
    //number of parameters to add :
    int elem_to_add = dual_size - current_size;
    LOG( INFO ) << "[selectDualApeeParameters] elem_to_add : "<<elem_to_add;
    new_parameters.resize( elem_to_add );
    parameter_type mu( M_Dmu );
    size_type index=0;
    for(int i=0; i<elem_to_add; i++)
    {
        bool already_exist;
        do
        {
            //initialization
            already_exist=false;
            //pick randomly an element in parameter space
            mu = M_Dmu->element();
            //make sure that the new mu is not already is M_dual_apee_mu or in M_WNmu
            for( auto _mu : *M_dual_apee_mu )
            {
                if( mu == _mu )
                    already_exist=true;
            }
            for( auto _mu : *M_WNmu )
            {
                if( mu == _mu )
                    already_exist=true;
            }
        }
        while( already_exist );
        M_dual_apee_mu->push_back( mu , index );
        new_parameters[i]=mu;
        //M_dual_apee_mu_complement = M_dual_apee_mu->complement();
    }
}

template<typename TruthModelType>
void
CRB<TruthModelType>::computePrimalApeeBasis( parameter_type const& mu ) const
{
    int N = M_WNmu->size();
    //compute RB solution
    std::vector< vectorN_type > uN_vector;
    std::vector< vectorN_type > uNdu_vector;
    std::vector< vectorN_type > uNold_vector;
    std::vector< vectorN_type > uNduold_vector;
    auto tuple = lb( N, mu, uN_vector, uNdu_vector, uNold_vector, uNduold_vector );

    vectorN_type uN;
    double time;
    if ( M_model->isSteady() )
    {
        uN = uNdu_vector[0];
        time=1e30;
    }
    else
        throw std::logic_error( "[CRB::computeApeeBasis] ERROR not yet implemented for transient case, if you want to run CRB on a transient problem, make sure that crb.use-accurate-apee=false" );

    auto Xpr = computeOnlinePrimalApeeVector( mu , uN );
    M_primal_apee_basis.push_back(Xpr);
}

template<typename TruthModelType>
typename CRB<TruthModelType>::vectorN_type
CRB<TruthModelType>::computeOnlinePrimalApeeVector( parameter_type const& mu , vectorN_type const & uN, vectorN_type const & uNold, double dt, double time ) const
{
    vectorN_type Xpr;
    int N = M_WNmu->size();
    LOG( INFO ) << "[computePrimalApeeBasis] N = "<<N;
    int Qa = M_model->Qa();
    int Qf = M_model->Ql( 0 );

    int primal_size = Qf*Qf + N*Qa*Qf + N*N*Qa*Qa;
    Xpr.resize(primal_size);

    beta_vector_type betaAqm;
    beta_vector_type betaMqm;
    std::vector<beta_vector_type> betaFqm;

    bool is_linear=M_model->isLinear();
    //get beta coefficients
    if( is_linear )
    {
        boost::tie( betaMqm, betaAqm, betaFqm ) = M_model->computeBetaQm( mu ,time );
    }
    else
    {
        //important note :
        //when lambda expressions will be totally operational
        //we will call computeBetaQm( uN, mu, tim )
        //and the test if( load_elements_db ) will disappear
        if( M_loadElementsDb )
            boost::tie( betaMqm, betaAqm, betaFqm ) = M_model->computeBetaQm( this->expansion( uN , N ), mu ,time );
        else
            boost::tie( betaMqm, betaAqm, betaFqm ) = M_model->computeBetaQm( mu ,time );
    }

    int idx=0;//start

    for ( int __q1 = 0; __q1 < Qf; ++__q1 )
    {
        for ( int __m1 = 0; __m1 < M_model->mMaxF(0,__q1); ++__m1 )
        {
            value_type fq1 = betaFqm[0][__q1][__m1];
            for ( int __q2 = 0; __q2 < Qf; ++__q2 )
            {
                for ( int __m2 = 0; __m2 < M_model->mMaxF(0,__q2); ++__m2 )
                {
                    value_type fq2 = betaFqm[0][__q2][__m2];
                    Xpr(idx) = fq1*fq2;
                    idx++;
                }//m2
            }//q2
        }//m1
    }//q1

    for ( int __q1 = 0; __q1 < Qa; ++__q1 )
    {
        for ( int __m1 = 0; __m1 < M_model->mMaxA(__q1); ++__m1 )
        {
            value_type aq1 = betaAqm[__q1][__m1];

            for ( int __q2 = 0; __q2 < Qf; ++__q2 )
            {
                for ( int __m2 = 0; __m2 < M_model->mMaxF(0,__q2); ++__m2 )
                {
                    value_type fq2 = betaFqm[0][__q2][__m2];
                    for(int n=0; n<N; n++)
                    {
                        Xpr(idx) = uN(n) *aq1 * fq2;
                        idx++;
                    }
                }//m2
            }//q2


            for ( int __q2 = 0; __q2 < Qa; ++__q2 )
            {
                for ( int __m2 = 0; __m2 < M_model->mMaxA(__q2); ++__m2 )
                {
                    value_type aq2 = betaAqm[__q2][__m2];
                    for(int n1=0; n1<N; n1++)
                    {
                        for(int n2=0; n2<N; n2++)
                        {
                            Xpr(idx) = uN(n1)*aq1*uN(n2)*aq2;
                            idx++;
                        }
                    }
                }//m2
            }//q2
        }//m1
    }//q1

    return Xpr;
}

template<typename TruthModelType>
void
CRB<TruthModelType>::computeDualApeeBasis( parameter_type const& mu ) const
{
    int N = M_WNmu->size();
    //compute RB solution
    std::vector< vectorN_type > uN_vector;
    std::vector< vectorN_type > uNdu_vector;
    std::vector< vectorN_type > uNold_vector;
    std::vector< vectorN_type > uNduold_vector;
    auto tuple = lb( N, mu, uN_vector, uNdu_vector, uNold_vector, uNduold_vector );

    vectorN_type uNdu;
    double time;

    if ( M_model->isSteady() )
    {
        uNdu = uNdu_vector[0];
        time=1e30;
    }
    else
        throw std::logic_error( "[CRB::computeApeeBasis] ERROR not yet implemented for transient case, if you want to run CRB on a transient problem, make sure that crb.use-accurate-apee=false" );

    auto Xdu = computeOnlineDualApeeVector( mu , uNdu );
    M_dual_apee_basis.push_back(Xdu);

}

template<typename TruthModelType>
typename CRB<TruthModelType>::vectorN_type
CRB<TruthModelType>::computeOnlineDualApeeVector( parameter_type const& mu , vectorN_type const & uNdu, vectorN_type const & uNold , double dt, double time ) const
{

    vectorN_type Xdu;
    int N = M_WNmu->size();
    int Qa = M_model->Qa();
    int Ql = M_model->Ql( M_output_index );

    int dual_size = Ql*Ql + N*Qa*Ql + N*N*Ql*Ql;
    Xdu.resize(dual_size);

    beta_vector_type betaAqm;
    beta_vector_type betaMqm;
    std::vector<beta_vector_type> betaFqm;

    bool is_linear=M_model->isLinear();
    //get beta coefficients
    if( is_linear )
    {
        boost::tie( betaMqm, betaAqm, betaFqm ) = M_model->computeBetaQm( mu ,time );
    }
    else
    {
        //important note :
        //when lambda expressions will be totally operational
        //we will call computeBetaQm( uN, mu, tim )
        //and the test if( load_elements_db ) will disappear
        if( M_loadElementsDb )
            boost::tie( betaMqm, betaAqm, betaFqm ) =  M_model->computeBetaQm( this->expansion( uNdu , N ), mu ,time );
        else
            boost::tie( betaMqm, betaAqm, betaFqm ) = M_model->computeBetaQm( mu ,time );
    }

    int idx=0;//start

    for ( int __q1 = 0; __q1 < Ql; ++__q1 )
    {
        for ( int __m1 = 0; __m1 < M_model->mMaxF(M_output_index,__q1); ++__m1 )
        {
            value_type fq1 = betaFqm[M_output_index][__q1][__m1];
            for ( int __q2 = 0; __q2 < Ql; ++__q2 )
            {
                for ( int __m2 = 0; __m2 < M_model->mMaxF(M_output_index,__q2); ++__m2 )
                {
                    value_type fq2 = betaFqm[M_output_index][__q2][__m2];
                    Xdu(idx) = fq1*fq2;
                    idx++;
                }//m2
            }//q2
        } //m1
    } //q1

    for ( int __q1 = 0; __q1 < Qa; ++__q1 )
    {
        for ( int __m1 = 0; __m1 < M_model->mMaxA(__q1); ++__m1 )
        {
            value_type aq1 = betaAqm[__q1][__m1];
            for ( int __q2 = 0; __q2 < Ql; ++__q2 )
            {
                for ( int __m2 = 0; __m2 < M_model->mMaxF(M_output_index,__q2); ++__m2 )
                {
                    value_type fq2 = betaFqm[M_output_index][__q2][__m2];
                    for(int n=0; n<N; n++)
                    {
                        Xdu(idx) = uNdu(n)*aq1*fq2;
                        idx++;
                    }
                }//m2
            }//q2

            for ( int __q2 = 0; __q2 < Qa; ++__q2 )
            {
                for ( int __m2 = 0; __m2 < M_model->mMaxA(__q2); ++__m2 )
                {
                    value_type aq2 = betaAqm[__q2][__m2];
                    for(int n1=0; n1<N; n1++)
                    {
                        for(int n2=0; n2<N; n2++)
                        {
                            Xdu(idx) = uNdu(n1)*aq1*uNdu(n2)*aq2;
                            idx++;
                        }
                    }
                }//m2
            }//q2
        }//m1
    }//q1

    return Xdu;
}

template<typename TruthModelType>
double
CRB<TruthModelType>::computeOnlinePrimalApee(  int N, parameter_type const& mu , vectorN_type const & uN, vectorN_type const & uNold , double dt, double time ) const
{
    int primal_rows = M_primal_T.rows();
    FEELPP_ASSERT( primal_rows >= N )( primal_rows )( N ).error( "invalid size of RB dimension or M_primal_T is not correctly filled");

    vectorN_type lambda;
    int Qa = M_model->Qa();
    int Qf = M_model->Ql( 0 );
    int primal_size = Qf*Qf + N*Qa*Qf + N*N*Qa*Qa; //d in the F.Casenave's paper
    lambda.resize(primal_size);

    if( M_computeMatrixInfo )
    {
        double cond = computeConditioning( M_primal_T );
        double det = M_primal_T.determinant();
        LOG( INFO ) << "[computeOnlinePrimalApee] matrix information :";
        LOG( INFO ) << "condition number = "<<cond;
        LOG( INFO ) << "determinant = "<<det;
    }

    auto X = computeOnlinePrimalApeeVector( mu , uN );
    lambda = M_primal_T.block( 0,0,N,N ).lu().solve( X );

    double result=0;
    for(int i=0; i<primal_size; i++)
    {
        result += lambda(i)*M_primal_V(i);
    }

    return result;
}

template<typename TruthModelType>
double
CRB<TruthModelType>::computeOnlineDualApee(  int N, parameter_type const& mu , vectorN_type const & uN, vectorN_type const & uNold , double dt, double time ) const
{
    int dual_rows = M_dual_T.rows();
    FEELPP_ASSERT( dual_rows >= N )( dual_rows )( N ).error( "invalid size of RB dimension or M_dual_T is not correctly filled");

    vectorN_type lambda;
    int Qa = M_model->Qa();
    int Ql = M_model->Ql( M_output_index );
    int dual_size = Ql*Ql + N*Qa*Ql + N*N*Ql*Ql; //d in the F.Casenave's paper
    lambda.resize(dual_size);

    if( M_computeMatrixInfo )
    {
        double cond = computeConditioning( M_dual_T );
        double det = M_dual_T.determinant();
        LOG( INFO ) << "[computeOnlineDualApee] matrix information :";
        LOG( INFO ) << "condition number = "<<cond;
        LOG( INFO ) << "determinant = "<<det;
    }

    auto X = computeOnlineDualApeeVector( mu , uN );
    lambda = M_dual_T.block( 0,0,N,N ).lu().solve( X );

    double result=0;
    for(int i=0; i<dual_size; i++)
    {
        result += lambda(i)*M_dual_V(i);
    }

    return result;
}

template<typename TruthModelType>
void
CRB<TruthModelType>::computationalTimeStatistics(std::string appname)
{

    double min=0,max=0,mean=0,standard_deviation=0;
    int n_eval = ioption(_prefix=M_prefix,_name="crb.computational-time-neval");

    vectorN_type time;
    Eigen::Matrix<double, Eigen::Dynamic, 1> time_crb_prediction;
    Eigen::Matrix<double, Eigen::Dynamic, 1> time_crb_error_estimation;
    time_crb_prediction.resize( n_eval );
    time_crb_error_estimation.resize( n_eval );

    sampling_ptrtype Sampling( new sampling_type( M_Dmu ) );
    Sampling->logEquidistribute( n_eval  );

    bool cvg = boption(_prefix=M_prefix,_name="crb.cvg-study");
    int dimension = ioption(_prefix=M_prefix,_name="crb.dimension");
    double tol = doption(_prefix=M_prefix,_name="crb.online-tolerance");

    int N=dimension;//by default we perform only one time statistics

    if( cvg ) //if we want to compute time statistics for every crb basis then we start at 1
        N=1;

    int proc_number =  this->worldComm().globalRank();
    int master =  this->worldComm().masterRank();

    //write on a file
    std::string file_name_prediction = "cvg-timing-crb-prediction.dat";
    std::string file_name_error_estimation = "cvg-timing-crb-error-estimation.dat";

    std::ofstream conv1;
    std::ofstream conv2;
    if( proc_number == master )
    {
        conv1.open(file_name_prediction, std::ios::app);
        conv1 << "NbBasis" << "\t" << "min" <<"\t"<< "max" <<"\t"<< "mean"<<"\t"<<"standard_deviation" << "\n";
        conv2.open(file_name_error_estimation, std::ios::app);
        conv2 << "NbBasis" << "\t" << "min" <<"\t"<< "max" <<"\t"<< "mean"<<"\t"<<"standard_deviation" << "\n";
    }

    //loop over basis functions (if cvg option)
    for(; N<=dimension; N++)
    {

        int mu_number = 0;
        for( auto mu : *Sampling )
        {
            //boost::mpi::timer tcrb;
            auto o = this->run( mu, time, tol , N );
            time_crb_prediction( mu_number ) = time(0) ;
            time_crb_error_estimation( mu_number ) = time(1) ;
            mu_number++;
        }

        auto stat_prediction = M_model->computeStatistics( time_crb_prediction , appname );
        auto stat_error_estimation = M_model->computeStatistics( time_crb_error_estimation , appname );
        min=stat_prediction(0);
        max=stat_prediction(1);
        mean=stat_prediction(2);
        standard_deviation=stat_prediction(3);
        if( proc_number == master )
            conv1 << N << "\t" << min << "\t" << max<< "\t"<< mean<< "\t"<< standard_deviation<<"\n";

        min=stat_error_estimation(0);
        max=stat_error_estimation(1);
        mean=stat_error_estimation(2);
        standard_deviation=stat_error_estimation(3);
        if( proc_number == master )
            conv2 << N << "\t" << min << "\t" << max<< "\t"<< mean<< "\t"<< standard_deviation<<"\n";

    }//loop over basis functions
    conv1.close();
    conv2.close();
}

#if 1
namespace detail
{
template<class Archive>
struct SaveRbSpaceCtxComposite
{
    SaveRbSpaceCtxComposite( Archive & ar, std::string const& rbctxNameInDb, boost::any const& ctxRbAny )
        :
        M_ar( ar ),
        M_rbctxNameInDb( rbctxNameInDb ),
        M_ctxRbAny( ctxRbAny ),
        M_currentSubSpaceIndex( 0 )
        {}
    template <typename T>
    void operator()( T & x ) const
        {
            typedef typename T::first_type key_type;
            typedef typename T::second_type::element_type rbsubspace_type;
            typedef typename rbsubspace_type::ctxrbset_ptrtype ctxrbset_ptrtype;
            if ( !boost::any_cast<ctxrbset_ptrtype>( &M_ctxRbAny ) )
            {
                ++M_currentSubSpaceIndex;
                return;
            }

            auto ctxRb = boost::any_cast<ctxrbset_ptrtype>( M_ctxRbAny );

            M_ar & boost::serialization::make_nvp( (M_rbctxNameInDb+"_spaceId").c_str(), M_currentSubSpaceIndex );
            M_ar & boost::serialization::make_nvp( M_rbctxNameInDb.c_str(), *ctxRb);
            ++M_currentSubSpaceIndex;
        }
    Archive & M_ar;
    std::string const& M_rbctxNameInDb;
    boost::any const& M_ctxRbAny;
    mutable int M_currentSubSpaceIndex;
};

template<typename TheRbSpaceType,class Archive>
void
saveRbSpaceCtx( TheRbSpaceType const& rbSpace, Archive & ar, std::string const& rbctxNameInDb, std::vector<boost::any> const& ctxRbAnyVec, mpl::false_ )
{
    typedef typename TheRbSpaceType::element_type::ctxrbset_ptrtype ctxrbset_ptrtype;
    CHECK( ctxRbAnyVec.size() <= 1 ) << "only one context in non composite";
    for ( boost::any const& ctxRbAny : ctxRbAnyVec )
    {
        CHECK( boost::any_cast<ctxrbset_ptrtype>( &ctxRbAny ) ) << "invalid rbctx cast in non composite";
        auto ctxRb = boost::any_cast<ctxrbset_ptrtype>( ctxRbAny );
        ar & boost::serialization::make_nvp( rbctxNameInDb.c_str(), *ctxRb );
    }
}
template<typename TheRbSpaceType,class Archive>
void
saveRbSpaceCtx( TheRbSpaceType const& rbSpace, Archive & ar, std::string const& rbctxNameInDb, std::vector<boost::any> const& ctxRbAnyVec, mpl::true_ )
{
    int numberOfCtx = ctxRbAnyVec.size();
    ar & boost::serialization::make_nvp( (rbctxNameInDb+"_number_of_context").c_str(), numberOfCtx );
    for ( int k=0;k<numberOfCtx;++k )
    {
        std::string rbctxNameInDb2 = fmt::format("{}_{}",rbctxNameInDb,k);
        boost::fusion::for_each( rbSpace->rbfunctionspaces(), SaveRbSpaceCtxComposite<Archive>( ar, rbctxNameInDb2, ctxRbAnyVec[k] ) );
    }
}
template<typename TheRbSpaceType,class Archive>
void
saveRbSpaceCtx( TheRbSpaceType const& rbSpace, Archive & ar, std::string const& rbctxNameInDb, std::vector<boost::any> const& ctxRbAnyVec )
{
    saveRbSpaceCtx( rbSpace,ar,rbctxNameInDb,ctxRbAnyVec, mpl::bool_<TheRbSpaceType::element_type::is_composite>() );
}


template<class Archive>
struct LoadRbSpaceCtxComposite
{
    LoadRbSpaceCtxComposite( Archive & ar, std::string const& rbctxNameInDb, int subspaceId, boost::any & ctxRbAny )
        :
        M_ar( ar ),
        M_rbctxNameInDb( rbctxNameInDb ),
        M_ctxRbAny( ctxRbAny ),
        M_subspaceId( subspaceId ),
        M_currentSubSpaceIndex( 0 )
        {}
    template <typename T>
    void operator()( T & x ) const
        {
            typedef typename T::first_type key_type;
            typedef typename T::second_type::element_type rbsubspace_type;
            typedef typename rbsubspace_type::ctxrbset_ptrtype ctxrbset_ptrtype;
            typedef typename rbsubspace_type::ctxrbset_type ctxrbset_type;
            if ( M_currentSubSpaceIndex != M_subspaceId )
            {
                ++M_currentSubSpaceIndex;
                return;
            }
            //std::cout << "load rbcomposite "<<M_rbctxNameInDb << " with index " << M_subspaceId << "\n";
            M_ctxRbAny = ctxrbset_ptrtype( new ctxrbset_type( x.second ) );
            auto ctxRb = boost::any_cast<ctxrbset_ptrtype>( M_ctxRbAny );
            M_ar & boost::serialization::make_nvp( M_rbctxNameInDb.c_str(), *ctxRb);
            ++M_currentSubSpaceIndex;
        }
    Archive & M_ar;
    std::string const& M_rbctxNameInDb;
    boost::any & M_ctxRbAny;
    int M_subspaceId;
    mutable int M_currentSubSpaceIndex;
};

template<typename TheRbSpaceType,class Archive>
void
loadRbSpaceCtx( TheRbSpaceType const& rbSpace, Archive & ar, std::string const& rbctxNameInDb, std::vector<boost::any> & ctxRbAny, mpl::false_ )
{
    typedef typename TheRbSpaceType::element_type::ctxrbset_ptrtype ctxrbset_ptrtype;
    typedef typename TheRbSpaceType::element_type::ctxrbset_type ctxrbset_type;
    ctxRbAny.resize( 1 );
    ctxRbAny[0] = ctxrbset_ptrtype( new ctxrbset_type( rbSpace ) );
    auto ctxRb = boost::any_cast<ctxrbset_ptrtype>( ctxRbAny[0] );
    ar & boost::serialization::make_nvp( rbctxNameInDb.c_str(), *ctxRb );
    ctxRb->setRbFunctionSpace( rbSpace );
}
template<typename TheRbSpaceType,class Archive>
void
loadRbSpaceCtx( TheRbSpaceType const& rbSpace, Archive & ar, std::string const& rbctxNameInDb, std::vector<boost::any> & ctxRbAny, mpl::true_ )
{
    int numberOfCtx = 0;
    ar & boost::serialization::make_nvp( (rbctxNameInDb+"_number_of_context").c_str(), numberOfCtx );
    ctxRbAny.resize( numberOfCtx );
    for ( int k=0;k<numberOfCtx;++k )
    {
        std::string rbctxNameInDb2 = fmt::format("{}_{}",rbctxNameInDb,k);
        int subspaceId = -1;
        ar & boost::serialization::make_nvp( (rbctxNameInDb2+"_spaceId").c_str(), subspaceId );
        boost::fusion::for_each( rbSpace->rbfunctionspaces(), LoadRbSpaceCtxComposite<Archive>( ar, rbctxNameInDb2, subspaceId, ctxRbAny[k] ) );
    }
}

template<typename TheRbSpaceType,class Archive>
void
loadRbSpaceCtx( TheRbSpaceType const& rbSpace, Archive & ar, std::string const& rbctxNameInDb, std::vector<boost::any> & ctxRbAny )
{
    loadRbSpaceCtx( rbSpace,ar,rbctxNameInDb,ctxRbAny, mpl::bool_<TheRbSpaceType::element_type::is_composite>() );
}

}
#endif


template<typename TruthModelType>
template<class Archive>
void
CRB<TruthModelType>::save( Archive & ar, const unsigned int version ) const
{
    int proc_number = this->worldComm().globalRank();

    LOG(INFO) <<"[CRB::save] version : "<<version<<std::endl;

    ar & boost::serialization::base_object<super>( *this );
    ar & BOOST_SERIALIZATION_NVP( M_output_index );
    ar & BOOST_SERIALIZATION_NVP( M_N );
    ar & BOOST_SERIALIZATION_NVP( M_rbconv );
    ar & BOOST_SERIALIZATION_NVP( M_error_type );
    ar & BOOST_SERIALIZATION_NVP( M_Xi );
    ar & BOOST_SERIALIZATION_NVP( M_WNmu );
    ar & BOOST_SERIALIZATION_NVP( M_Aqm_pr );
    //ar & BOOST_SERIALIZATION_NVP( M_hAqm );
    ar & BOOST_SERIALIZATION_NVP( M_Aqm_du );
    ar & BOOST_SERIALIZATION_NVP( M_Aqm_pr_du );
    ar & BOOST_SERIALIZATION_NVP( M_Fqm_pr );
    //ar & BOOST_SERIALIZATION_NVP( M_hFqm );
    ar & BOOST_SERIALIZATION_NVP( M_Fqm_du );
    ar & BOOST_SERIALIZATION_NVP( M_Lqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Lqm_du );

    ar & BOOST_SERIALIZATION_NVP( M_C0_pr );
    ar & BOOST_SERIALIZATION_NVP( M_C0_du );
    ar & BOOST_SERIALIZATION_NVP( M_Lambda_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Lambda_du );
    ar & BOOST_SERIALIZATION_NVP( M_Gamma_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Gamma_du );


    ar & BOOST_SERIALIZATION_NVP( M_Mqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Mqm_du );
    ar & BOOST_SERIALIZATION_NVP( M_Mqm_pr_du );

    if ( model_type::is_time_dependent )
    {

        ar & BOOST_SERIALIZATION_NVP( M_coeff_pr_ini_online );
        ar & BOOST_SERIALIZATION_NVP( M_coeff_du_ini_online );
        ar & BOOST_SERIALIZATION_NVP( M_Cmf_pr );
        ar & BOOST_SERIALIZATION_NVP( M_Cma_pr );
        ar & BOOST_SERIALIZATION_NVP( M_Cmm_pr );
        ar & BOOST_SERIALIZATION_NVP( M_Cmf_du );
        ar & BOOST_SERIALIZATION_NVP( M_Cmf_du_ini );
        ar & BOOST_SERIALIZATION_NVP( M_Cma_du );
        ar & BOOST_SERIALIZATION_NVP( M_Cmm_du );
    }

    ar & BOOST_SERIALIZATION_NVP( M_algebraicInnerProductPrimal );

    ar & BOOST_SERIALIZATION_NVP ( M_database_contains_variance_info );
    if( M_database_contains_variance_info )
        ar & BOOST_SERIALIZATION_NVP( M_variance_matrix_phi );

    ar & BOOST_SERIALIZATION_NVP( M_Fqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_InitialGuessV_pr );

    ar & BOOST_SERIALIZATION_NVP( M_current_mu );
    ar & BOOST_SERIALIZATION_NVP( M_no_residual_index );

    ar & BOOST_SERIALIZATION_NVP( M_maxerror );
    ar & BOOST_SERIALIZATION_NVP( M_use_newton );
    ar & BOOST_SERIALIZATION_NVP( M_Jqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Rqm_pr );
    //ar & BOOST_SERIALIZATION_NVP( M_hRqm );

    ar & BOOST_SERIALIZATION_NVP( M_primal_apee_basis );
    ar & BOOST_SERIALIZATION_NVP( M_dual_apee_basis );
    ar & BOOST_SERIALIZATION_NVP( M_primal_V );
    ar & BOOST_SERIALIZATION_NVP( M_dual_V );
    ar & BOOST_SERIALIZATION_NVP( M_primal_T );
    ar & BOOST_SERIALIZATION_NVP( M_dual_T );

    ar & BOOST_SERIALIZATION_NVP( M_model_executed_in_steady_mode );

    if( version > 0 )
    {
        ar & BOOST_SERIALIZATION_NVP( M_C0_pr_eim );
        ar & BOOST_SERIALIZATION_NVP( M_C0_du_eim );
        ar & BOOST_SERIALIZATION_NVP( M_Lambda_pr_eim );
        ar & BOOST_SERIALIZATION_NVP( M_Lambda_du_eim );
        ar & BOOST_SERIALIZATION_NVP( M_Gamma_pr_eim );
        ar & BOOST_SERIALIZATION_NVP( M_Gamma_du_eim );

        if ( model_type::is_time_dependent )
        {
            ar & BOOST_SERIALIZATION_NVP( M_Cmf_pr_eim );
            ar & BOOST_SERIALIZATION_NVP( M_Cma_pr_eim );
            ar & BOOST_SERIALIZATION_NVP( M_Cmm_pr_eim );
            ar & BOOST_SERIALIZATION_NVP( M_Cmf_du_eim );
            ar & BOOST_SERIALIZATION_NVP( M_Cmf_du_ini_eim );
            ar & BOOST_SERIALIZATION_NVP( M_Cma_du_eim );
            ar & BOOST_SERIALIZATION_NVP( M_Cmm_du_eim );
        }

    }

#if 1
    ar & BOOST_SERIALIZATION_NVP( M_hasRbSpaceContextEim );
    if ( M_hasRbSpaceContextEim )
    {
        auto rbSpaceContextEim = M_model->model()->rbSpaceContextEim();
        std::vector<std::string> rbSpaceContextEimNames;
        for ( auto const& rbSpaceContextPair : rbSpaceContextEim )
            rbSpaceContextEimNames.push_back( rbSpaceContextPair.first );
        ar & BOOST_SERIALIZATION_NVP( rbSpaceContextEimNames );

        for ( std::string const& eimName : rbSpaceContextEimNames )
        {
            std::string rbctxNameInDb = fmt::format("rbSpaceContextEim_{}",eimName);
            Feel::detail::saveRbSpaceCtx( M_model->rBFunctionSpace(), ar, rbctxNameInDb, rbSpaceContextEim.find( eimName )->second );
        }
    }
#endif

#if 0
    for(int i=0; i<M_N; i++)
        ar & BOOST_SERIALIZATION_NVP( M_WN[i] );
    for(int i=0; i<M_N; i++)
        ar & BOOST_SERIALIZATION_NVP( M_WNdu[i] );

    auto mesh = mesh_type::New();
    auto is_mesh_loaded = mesh->load( _name="mymesh",_path=this->dbLocalPath(),_type="binary" );

    if ( ! is_mesh_loaded )
    {
        auto first_element = M_WN[0];
        mesh = first_element.functionSpace()->mesh() ;
        mesh->save( _name="mymesh",_path=this->dbLocalPath(),_type="binary" );
    }

#endif

}


template<typename TruthModelType>
template<class Archive>
void
CRB<TruthModelType>::load( Archive & ar, const unsigned int version )
{

    //if( version <= 4 )
    //    throw std::logic_error( "[CRB::load] ERROR while loading the existing database, since version 5 there was many changes. Please use the option --crb.rebuild-database=true " );
    int proc_number = this->worldComm().globalRank();

    LOG(INFO) <<"[CRB::load] version"<< version <<"\n";

#if 0
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
#endif

    typedef boost::bimap< int, double > old_convergence_type;
    ar & boost::serialization::base_object<super>( *this );
    ar & BOOST_SERIALIZATION_NVP( M_output_index );
    ar & BOOST_SERIALIZATION_NVP( M_N );

	ar & BOOST_SERIALIZATION_NVP( M_rbconv );

    ar & BOOST_SERIALIZATION_NVP( M_error_type );
    ar & BOOST_SERIALIZATION_NVP( M_Xi );
    ar & BOOST_SERIALIZATION_NVP( M_WNmu );
    ar & BOOST_SERIALIZATION_NVP( M_Aqm_pr );
    //ar & BOOST_SERIALIZATION_NVP( M_hAqm );
    ar & BOOST_SERIALIZATION_NVP( M_Aqm_du );
    ar & BOOST_SERIALIZATION_NVP( M_Aqm_pr_du );
    ar & BOOST_SERIALIZATION_NVP( M_Fqm_pr );
    //ar & BOOST_SERIALIZATION_NVP( M_hFqm );
    ar & BOOST_SERIALIZATION_NVP( M_Fqm_du );
    ar & BOOST_SERIALIZATION_NVP( M_Lqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Lqm_du );
    ar & BOOST_SERIALIZATION_NVP( M_C0_pr );
    ar & BOOST_SERIALIZATION_NVP( M_C0_du );
    ar & BOOST_SERIALIZATION_NVP( M_Lambda_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Lambda_du );
    ar & BOOST_SERIALIZATION_NVP( M_Gamma_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Gamma_du );

    ar & BOOST_SERIALIZATION_NVP( M_Mqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Mqm_du );
    ar & BOOST_SERIALIZATION_NVP( M_Mqm_pr_du );

    if ( model_type::is_time_dependent )
    {
        ar & BOOST_SERIALIZATION_NVP( M_coeff_pr_ini_online );
        ar & BOOST_SERIALIZATION_NVP( M_coeff_du_ini_online );
        ar & BOOST_SERIALIZATION_NVP( M_Cmf_pr );
        ar & BOOST_SERIALIZATION_NVP( M_Cma_pr );
        ar & BOOST_SERIALIZATION_NVP( M_Cmm_pr );
        ar & BOOST_SERIALIZATION_NVP( M_Cmf_du );
        ar & BOOST_SERIALIZATION_NVP( M_Cmf_du_ini );
        ar & BOOST_SERIALIZATION_NVP( M_Cma_du );
        ar & BOOST_SERIALIZATION_NVP( M_Cmm_du );
    }

    ar & BOOST_SERIALIZATION_NVP( M_algebraicInnerProductPrimal );

    ar & BOOST_SERIALIZATION_NVP ( M_database_contains_variance_info );
    if( M_database_contains_variance_info )
        ar & BOOST_SERIALIZATION_NVP( M_variance_matrix_phi );

    ar & BOOST_SERIALIZATION_NVP( M_Fqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_InitialGuessV_pr );

    ar & BOOST_SERIALIZATION_NVP( M_current_mu );
    ar & BOOST_SERIALIZATION_NVP( M_no_residual_index );

    ar & BOOST_SERIALIZATION_NVP( M_maxerror );
    ar & BOOST_SERIALIZATION_NVP( M_use_newton );
    ar & BOOST_SERIALIZATION_NVP( M_Jqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Rqm_pr );
    //ar & BOOST_SERIALIZATION_NVP( M_hRqm );

    if( boption(_prefix=M_prefix,_name="crb.use-newton") != M_use_newton  )
    {
        if( M_use_newton )
        {
            if( this->worldComm().isMasterRank() )
                std::cout<<"[CRB::loadDB] WARNING in the database used the option use-newton=true and it's not the case in your options so make sure that crb.rebuild-database=true !" <<std::endl;
            LOG( INFO )<<"[CRB::loadDB] WARNING in the database used the option use-newton=true and it's not the case in your options so make sure that crb.rebuild-database=true !" ;
        }
        else
        {
            if( this->worldComm().isMasterRank() )
                std::cout<< "[CRB::loadDB] WARNING in the database used the option use-newton=false and it's not the case in your options so make sure that crb.rebuild-database=true !"<<std::endl;
            LOG( INFO )<< "[CRB::loadDB] WARNING in the database used the option use-newton=false and it's not the case in your options so make sure that crb.rebuild-database=true !";
        }
    }

    ar & BOOST_SERIALIZATION_NVP( M_primal_apee_basis );
    ar & BOOST_SERIALIZATION_NVP( M_dual_apee_basis );
    ar & BOOST_SERIALIZATION_NVP( M_primal_V );
    ar & BOOST_SERIALIZATION_NVP( M_dual_V );
    ar & BOOST_SERIALIZATION_NVP( M_primal_T );
    ar & BOOST_SERIALIZATION_NVP( M_dual_T );

    ar & BOOST_SERIALIZATION_NVP( M_model_executed_in_steady_mode );
    bool current_option=boption(_prefix=M_prefix,_name="crb.is-model-executed-in-steady-mode");
    if( M_model_executed_in_steady_mode != current_option )
    {
        if( M_model_executed_in_steady_mode && this->worldComm().isMasterRank() )
            std::cout<<"[CRB::loadDB] WARNING in the database used, the model was executed in steady mode but now you want to execute it in transient mode. make sure that --crb.rebuild-database=true"<<std::endl;
        LOG( INFO ) <<"[CRB::loadDB] WARNING in the database used, the model was executed in steady mode but now you want to execute it in transient mode. make sure that --crb.rebuild-database=true";
    }

    //For version == 0 there was no error estimation on EIM
    if( version > 0 )
    {
        ar & BOOST_SERIALIZATION_NVP( M_C0_pr_eim );
        ar & BOOST_SERIALIZATION_NVP( M_C0_du_eim );
        ar & BOOST_SERIALIZATION_NVP( M_Lambda_pr_eim );
        ar & BOOST_SERIALIZATION_NVP( M_Lambda_du_eim );
        ar & BOOST_SERIALIZATION_NVP( M_Gamma_pr_eim );
        ar & BOOST_SERIALIZATION_NVP( M_Gamma_du_eim );

        if ( model_type::is_time_dependent )
        {
            ar & BOOST_SERIALIZATION_NVP( M_Cmf_pr_eim );
            ar & BOOST_SERIALIZATION_NVP( M_Cma_pr_eim );
            ar & BOOST_SERIALIZATION_NVP( M_Cmm_pr_eim );
            ar & BOOST_SERIALIZATION_NVP( M_Cmf_du_eim );
            ar & BOOST_SERIALIZATION_NVP( M_Cmf_du_ini_eim );
            ar & BOOST_SERIALIZATION_NVP( M_Cma_du_eim );
            ar & BOOST_SERIALIZATION_NVP( M_Cmm_du_eim );
        }

    }// version > 0 => EIM error estimation


#if 1
    ar & BOOST_SERIALIZATION_NVP( M_hasRbSpaceContextEim );
    if ( M_hasRbSpaceContextEim )
    {
        // reload rbcontext only if rBFunctionSpace is init
        if ( M_model && M_model->rBFunctionSpace() )
        {
            //std::cout << "CRB reload rbContext with a rBFunctionSpace\n";
            std::vector<std::string> rbSpaceContextEimNames;
            ar & BOOST_SERIALIZATION_NVP( rbSpaceContextEimNames );

            std::map<std::string,std::vector<boost::any> > rbSpaceContextEim;
            for ( std::string const& eimName : rbSpaceContextEimNames )
            {
                std::string rbctxNameInDb = fmt::format("rbSpaceContextEim_{}",eimName);
                Feel::detail::loadRbSpaceCtx( M_model->rBFunctionSpace(), ar, rbctxNameInDb, rbSpaceContextEim[eimName] );
            }
            M_model->model()->setRbSpaceContextEim( rbSpaceContextEim );

        }
        else
            M_hasRbSpaceContextEim = false;
    }
#endif

#if 0
    std::cout << "[loadDB] output index : " << M_output_index << "\n"
              << "[loadDB] N : " << M_N << "\n"
              << "[loadDB] error type : " << M_error_type << "\n";

    for ( auto it = M_rbconv.begin(), en = M_rbconv.end();it != en; ++it )
        std::cout << "[loadDB] convergence: (" << it->left << ","  << it->right  << ")\n";

    element_type temp = Xh->element();

    M_WN.resize( M_N );
    M_WNdu.resize( M_N );

    for( int i = 0 ; i < M_N ; i++ )
    {
        temp.setName( fmt::format( "fem-primal-{}", i  );
        ar & BOOST_SERIALIZATION_NVP( temp );
        M_WN[i] = temp;
    }

    for( int i = 0 ; i < M_N ; i++ )
    {
        temp.setName( fmt::format( "fem-dual-{}", i ) );
        ar & BOOST_SERIALIZATION_NVP( temp );
        M_WNdu[i] = temp;
    }

#endif
    LOG(INFO) << "[CRB::load] end of load function" << std::endl;
}


template<typename TruthModelType>
bool
CRB<TruthModelType>::printErrorDuringOfflineStep()
{
    bool print = boption(_prefix=M_prefix,_name="crb.print-error-during-rb-construction");
    return print;
}
#if 0
template<typename TruthModelType>
bool
CRB<TruthModelType>::rebuildDB()
{
    bool rebuild_db = boption(_prefix=M_prefix,_name="crb.rebuild-database");
    int Nrestart = ioption(_prefix=M_prefix,_name="crb.restart-from-N");
    bool rebuild=false;
    //if ( rebuild_db && Nrestart < 1 )
    if ( rebuild_db || Nrestart == 0 )
    {
        std::cout << "rebuild_db = " << rebuild_db << std::endl;
        std::cout << "Nrestart = " << Nrestart << std::endl;
        rebuild=true;
    }
    if( M_N == 0 )
    {
        std::cout << "M_N = " << M_N << std::endl;
        rebuild=true;
    }
    return rebuild;
}
#endif

template<typename TruthModelType>
void
CRB<TruthModelType>::generateSuperSampling()
{
    LOG(INFO) << "[CRB::offline] compute super sampling\n";

    int proc_number = worldComm().globalRank();
    int total_proc = worldComm().globalSize();
    bool all_proc_same_sampling = boption(_prefix=M_prefix, _name="crb.all-procs-have-same-sampling" );
    int sampling_size = ioption(_prefix=M_prefix,_name="crb.sampling-size");
    std::string sampling_mode = soption(_prefix=M_prefix,_name="crb.sampling-mode");

    std::string file_name;
    if( all_proc_same_sampling )
        file_name = fmt::format("M_Xi_{}_{}", sampling_size, sampling_mode );
    else
        file_name = fmt::format("M_Xi_{}_{}-proc{}on{}", sampling_size, sampling_mode, proc_number, total_proc );

    std::ifstream file ( file_name );

    if ( !file )
    {
        std::string supersamplingname =fmt::format("Dmu-{}-generated-by-master-proc",sampling_size );
        if( sampling_mode == "random" )
            this->M_Xi->randomize( sampling_size , all_proc_same_sampling, supersamplingname, false );
        else if( sampling_mode == "log-random" )
            this->M_Xi->randomize( sampling_size , all_proc_same_sampling , supersamplingname, true );
        else if( sampling_mode == "log-equidistribute" )
            this->M_Xi->logEquidistribute( sampling_size , all_proc_same_sampling , supersamplingname );
        else if( sampling_mode == "equidistribute" )
            this->M_Xi->equidistribute( sampling_size , all_proc_same_sampling , supersamplingname );
        else
            throw std::logic_error( "[CRBSaddlePoint::offline] ERROR invalid option crb.sampling-mode, please select between log-random, log-equidistribute or equidistribute" );

        if ( all_proc_same_sampling )
            this->worldComm().barrier();
        if ( !all_proc_same_sampling || this->worldComm().isMasterRank() )
            M_Xi->writeOnFile(file_name);
    }
    else
    {
        this->M_Xi->clear();
        this->M_Xi->readFromFile(file_name);
    }

    this->M_WNmu->setSuperSampling( this->M_Xi );
} //generateSuperSampling()

template<typename TruthModelType>
bool
CRB<TruthModelType>::buildSampling()
{
    bool use_predefined_WNmu = boption(_prefix=M_prefix,_name="crb.use-predefined-WNmu");
    int N_log_equi = ioption(_prefix=M_prefix,_name="crb.use-logEquidistributed-WNmu");
    int N_equi = ioption(_prefix=M_prefix,_name="crb.use-equidistributed-WNmu");
    int N_random = ioption(_prefix=M_prefix,_name="crb.use-random-WNmu" );

    std::string file_name = fmt::format("SamplingWNmu");
    std::ifstream file ( file_name );

    if ( use_predefined_WNmu ) // In this case we want to read the sampling
    {
        this->M_WNmu->clear();
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
    else if ( this->M_error_type==CRB_NO_RESIDUAL )// We generate the sampling with choosen strategy
    {
        this->M_WNmu->clear();
        std::ifstream file ( file_name );
        if ( file && boption(_prefix=M_prefix,_name="crb.reload-last-sampling") )
        {
            Feel::cout << "[CRB::offline] Reload last sampling\n";
            this->M_WNmu->readFromFile(file_name);
        }
        else
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
                bool use_log = boption(_prefix=M_prefix,_name="crb.randomize.use-log");
                this->M_WNmu->randomize( N_random , true, "",use_log );
                if( Environment::isMasterRank() )
                    std::cout<<"[CRB::offline] Randomize WNmu ( sampling size : "
                             <<N_random<<" )"<<std::endl;
                LOG( INFO )<<"[CRB::offline] Randomize WNmu ( sampling size : "
                           <<N_random<<" )";
            }
            else // In this case we don't know what sampling to use
                throw std::logic_error( "[CRB::offline] ERROR : You have to choose an appropriate strategy for the offline sampling : random, equi, logequi or predefined" );

        }

        this->M_WNmu->writeOnFile(file_name);
        use_predefined_WNmu=true;
    } //build sampling

    return use_predefined_WNmu;
} //buildSampling()

template<typename TruthModelType>
void
CRB<TruthModelType>::addBasis( element_type& u, element_type& udu, parameter_type& mu )
{
    tic();
    M_model->rBFunctionSpace()->addPrimalBasisElement( u );
    toc("Add Primal Basis Function");
    tic();
    M_model->rBFunctionSpace()->addDualBasisElement( udu );
    toc("Add Dual Basis Function");
}

template<typename TruthModelType>
void
CRB<TruthModelType>::orthonormalizeBasis( int number_of_added_elements )
{
    bool POD_WN = boption(_prefix=M_prefix,_name="crb.apply-POD-to-WN") ;
    if(  POD_WN &&  ! M_model->isSteady() )
    {
        pod_ptrtype POD = pod_ptrtype( new pod_type() );
        POD->setModel( M_model );
        mode_set_type ModeSet;
        POD->setNm( M_N );
        bool use_solutions=false;
        bool is_primal=true;
        POD->pod( ModeSet, is_primal, M_model->rBFunctionSpace()->primalRB() , use_solutions );
        M_model->rBFunctionSpace()->setPrimalBasis( ModeSet );
        if( M_solve_dual_problem )
        {
            ModeSet.clear();
            POD->pod( ModeSet, false,  M_model->rBFunctionSpace()->dualRB() , use_solutions );
            M_model->rBFunctionSpace()->setDualBasis( ModeSet );
        }
    }
    else
    {
        double norm_max = doption(_prefix=M_prefix,_name="crb.orthonormality-tol");
        int max_iter = ioption(_prefix=M_prefix,_name="crb.orthonormality-max-iter");
        if ( M_orthonormalize_primal )
        {
            tic();
            double norm = norm_max+1;
            int iter=0;
            double old = 10;
            while( norm >= norm_max && iter < max_iter)
            {
                norm = orthonormalize( M_N, M_model->rBFunctionSpace()->primalRB(), number_of_added_elements );
                iter++;
                //if the norm doesn't change
                if( math::abs(old-norm) < norm_max )
                    norm=0;
                old=norm;
            }
            M_model->rBFunctionSpace()->updatePrimalBasisForUse();
            toc("Primal Orthonormalization");
        }
        if ( M_orthonormalize_dual && M_solve_dual_problem )
        {
            tic();
            double norm = norm_max+1;
            int iter=0;
            double old = 10;
            while( norm >= norm_max && iter < max_iter )
            {
                norm = orthonormalize( M_N, M_model->rBFunctionSpace()->dualRB() , number_of_added_elements );
                iter++;
                if( math::abs(old-norm) < norm_max )
                    norm=0;
                old=norm;
            }
            M_model->rBFunctionSpace()->updateDualBasisForUse();
            toc("Dual Orthonormalization");
        }
    }//orthonormalization
}


template<typename TruthModelType>
void
CRB<TruthModelType>::buildRbMatrix( int number_of_added_elements, parameter_type& mu, element_ptrtype dual_initial_field )
{
    tic();
    // we only compute the last line and last column of reduced matrices (last added elements)
    int number_of_elements_to_update = number_of_added_elements;
    // in the case of cobuild, we have to update all since affine decomposition has changed
    if( ioption(_prefix=M_prefix,_name="ser.rb-frequency") != 0 && !M_rebuild)
        number_of_elements_to_update = M_N;
    // In case of SER use + error estimation, we compute \hat{A}, \hat{F} (resp. \hat{R}) to compute norm of residual (Riesz)
    int ser_error_estimation = M_SER_errorEstimation;

    if( ! M_use_newton )
    {
        LOG(INFO) << "[CRB::offline] compute Aq_pr, Aq_du, Aq_pr_du" << "\n";

        M_hAqm.resize( M_model->Qa() );
        for  (size_type q = 0; q < M_model->Qa(); ++q )
        {
            M_hAqm[q].resize( M_model->mMaxA(q) );
            for( size_type m = 0; m < M_model->mMaxA(q); ++m )
            {
                M_Aqm_pr[q][m].conservativeResize( M_N, M_N );
                M_Aqm_du[q][m].conservativeResize( M_N, M_N );
                M_Aqm_pr_du[q][m].conservativeResize( M_N, M_N );

                for ( size_type i = M_N - number_of_elements_to_update; i < M_N; i++ )
                {
                    for ( size_type j = 0; j < M_N; ++j )
                    {
                        M_Aqm_pr[q][m]( i, j ) = M_model->Aqm(q , m , M_model->rBFunctionSpace()->primalBasisElement(j), M_model->rBFunctionSpace()->primalBasisElement(i) );//energy
                        M_Aqm_du[q][m]( i, j ) = M_model->Aqm( q , m , M_model->rBFunctionSpace()->dualBasisElement(i), M_model->rBFunctionSpace()->dualBasisElement(j), true );
                        M_Aqm_pr_du[q][m]( i, j ) = M_model->Aqm(q , m , M_model->rBFunctionSpace()->primalBasisElement(j), M_model->rBFunctionSpace()->dualBasisElement(i) );
                    }
                }

                for ( size_type j=M_N - number_of_elements_to_update; j < M_N; j++ )
                {
                    for ( size_type i = 0; i < M_N; ++i )
                    {
                        M_Aqm_pr[q][m]( i, j ) = M_model->Aqm(q , m , M_model->rBFunctionSpace()->primalBasisElement(j), M_model->rBFunctionSpace()->primalBasisElement(i) );
                        M_Aqm_du[q][m]( i, j ) = M_model->Aqm(q , m , M_model->rBFunctionSpace()->dualBasisElement(i), M_model->rBFunctionSpace()->dualBasisElement(j) , true );
                        M_Aqm_pr_du[q][m]( i, j ) = M_model->Aqm(q , m , M_model->rBFunctionSpace()->primalBasisElement(j), M_model->rBFunctionSpace()->dualBasisElement(i) );
                    }
                }

                if( ser_error_estimation )
                {
                    M_hAqm[q][m].resize( M_N );
                    auto Aqm_xi_n = backend()->newVector( M_model->functionSpace() );
                    auto xi_n = backend()->newVector( M_model->functionSpace() );
                    for( int n=0; n < M_N; n++ )
                    {
                        M_hAqm[q][m][n] = backend()->newVector( M_model->functionSpace() );
                        *xi_n = M_model->rBFunctionSpace()->primalBasisElement(n);
                        auto Aqm = M_model->Aqm(q,m);
                        Aqm->multVector( xi_n, Aqm_xi_n );
                        M_model->l2solve( M_hAqm[q][m][n], Aqm_xi_n );
                    }
                }
            }//loop over m
        }//loop over q

        LOG(INFO) << "[CRB::offline] compute Mq_pr, Mq_du, Mq_pr_du" << "\n";

        LOG(INFO) << "[CRB::offline] compute Fq_pr, Fq_du" << "\n";
        M_hFqm.resize( M_model->Ql( 0 ) );
        for ( size_type q = 0; q < M_model->Ql( 0 ); ++q )
        {
            M_hFqm[q].resize( M_model->mMaxF( 0, q ) );
            for( size_type m = 0; m < M_model->mMaxF( 0, q ); ++m )
            {
                M_Fqm_pr[q][m].conservativeResize( M_N );
                M_Fqm_du[q][m].conservativeResize( M_N );
                for ( size_type l = 1; l <= number_of_elements_to_update; ++l )
                {
                    int index = M_N-l;
                    M_Fqm_pr[q][m]( index ) = M_model->Fqm( 0, q, m, M_model->rBFunctionSpace()->primalBasisElement( index ) );
                    M_Fqm_du[q][m]( index ) = M_model->Fqm( 0, q, m, M_model->rBFunctionSpace()->dualBasisElement( index ) );
                }
                if( ser_error_estimation )
                {
                    M_hFqm[q][m] = backend()->newVector( M_model->functionSpace() );
                    auto fqm = M_model->Fqm(0, q, m);
                    M_model->l2solve( M_hFqm[q][m], fqm );
                }
            }//loop over m (>= mMaxF - cobuild_eim_freq)
        }//loop over q

    }//end of "if ! use_newton"

    if( M_use_newton )
    {
        LOG(INFO) << "[CRB::offline] compute Jq_pr " << "\n";

        for  (size_type q = 0; q < M_model->Qa(); ++q )
        {
            for( size_type m = 0; m < M_model->mMaxA(q); ++m )
            {
                M_Jqm_pr[q][m].conservativeResize( M_N, M_N );

                for ( size_type i = M_N - number_of_elements_to_update; i < M_N; i++ )
                {
                    for ( size_type j = 0; j < M_N; ++j )
                    {
                        M_Jqm_pr[q][m]( i, j ) = M_model->Jqm(q , m ,
                                                              M_model->rBFunctionSpace()->primalBasisElement(j),
                                                              M_model->rBFunctionSpace()->primalBasisElement(i) );
                        if (i!=j)
                            M_Jqm_pr[q][m]( j, i ) = M_model->Jqm(q , m ,
                                                                  M_model->rBFunctionSpace()->primalBasisElement(i),
                                                                  M_model->rBFunctionSpace()->primalBasisElement(j) );
                    }
                }
            }//loop over m
        }//loop over q

        if( ser_error_estimation )
        {
            auto RF_A = M_model->RF_Aqm();
            M_hAqm.resize( RF_A.size() );
            for(int q = 0; q < RF_A.size(); ++q )
            {
                M_hAqm[q].resize( RF_A[q].size() );
                for(int m = 0; m < RF_A[q].size(); ++m )
                {
                    M_hAqm[q][m].resize( M_N );
                    auto Aqm_xi_n = backend()->newVector( M_model->functionSpace() );
                    auto xi_n = backend()->newVector( M_model->functionSpace() );
                    for( int n=0; n < M_N; n++ )
                    {
                        M_hAqm[q][m][n] = backend()->newVector( M_model->functionSpace() );
                        *xi_n = M_model->rBFunctionSpace()->primalBasisElement(n);
                        auto Aqm = RF_A[q][m];
                        Aqm->multVector( xi_n, Aqm_xi_n );
                        M_model->l2solve( M_hAqm[q][m][n], Aqm_xi_n );
                    }
                }//m
            }//q
        }

        LOG(INFO) << "[CRB::offline] compute Rq_pr" << "\n";

        for ( size_type q = 0; q < M_model->Ql( 0 ); ++q )
        {
            for( size_type m = 0; m < M_model->mMaxF( 0, q ); ++m )
            {
                M_Rqm_pr[q][m].conservativeResize( M_N );

                for ( size_type l = 1; l <= number_of_elements_to_update; ++l )
                {
                    int index = M_N-l;
                    //M_Rqm_pr[q][m]( index ) = inner_product( *M_Rqm[0][q][m] , M_model->rBFunctionSpace()->primalBasisElement(index) );
                    M_Rqm_pr[q][m]( index ) = M_model->Fqm( 0, q, m, M_model->rBFunctionSpace()->primalBasisElement(index) );
                }
            }//loop over m
        }//loop over q

        if( ser_error_estimation )
        {
            auto RF_F = M_model->RF_Fqm();
            M_hFqm.resize( RF_F.size() );
            for(int q = 0; q < RF_F.size(); ++q )
            {
                M_hFqm[q].resize( RF_F[q].size() );
                for(int m = 0; m < RF_F[q].size(); ++m )
                {
                    M_hFqm[q][m] = backend()->newVector( M_model->functionSpace() );
                    auto fqm = RF_F[q][m];
                    M_model->l2solve( M_hFqm[q][m], fqm );
                }//m
            }//q
        }

    }//end if use_newton case


    if( !M_model->isLinear() )
    {
        int q_max = M_model->QInitialGuess();
        for ( size_type q = 0; q < q_max; ++q )
        {
            int m_max =M_model->mMaxInitialGuess(q);
            for( size_type m = 0; m < m_max; ++m )
            {
                M_InitialGuessV_pr[q][m].conservativeResize( M_N );
                for ( size_type l = 1; l <= number_of_elements_to_update; ++l )
                {
                    int index = M_N-l;
                    for( int idx = 0; idx<=index; idx++ )
                        M_InitialGuessV_pr[q][m]( idx ) = M_model->InitialGuessVqm( q, m, M_model->rBFunctionSpace()->primalBasisElement(idx) );
                }
            }
        }
    }


    for ( size_type q = 0; q < M_model->Qm(); ++q )
    {
        for( size_type m = 0; m < M_model->mMaxM(q); ++m )
        {
            M_Mqm_pr[q][m].conservativeResize( M_N, M_N );
            M_Mqm_du[q][m].conservativeResize( M_N, M_N );
            M_Mqm_pr_du[q][m].conservativeResize( M_N, M_N );

            for ( size_type i=M_N - number_of_elements_to_update; i < M_N; i++ )
            {
                for ( size_type j = 0; j < M_N; ++j )
                {
                    M_Mqm_pr[q][m]( i, j ) = M_model->Mqm(q, m,
                                                          M_model->rBFunctionSpace()->primalBasisElement(i),
                                                          M_model->rBFunctionSpace()->primalBasisElement(j) );
                    M_Mqm_du[q][m]( i, j ) = M_model->Mqm(q, m,
                                                          M_model->rBFunctionSpace()->dualBasisElement(i),
                                                          M_model->rBFunctionSpace()->dualBasisElement(j), true );
                    M_Mqm_pr_du[q][m]( i, j ) = M_model->Mqm( q, m,
                                                              M_model->rBFunctionSpace()->dualBasisElement(i),
                                                              M_model->rBFunctionSpace()->primalBasisElement(j) );
                }
            }
            for ( size_type j = M_N - number_of_elements_to_update; j < M_N ; j++ )
            {
                for ( size_type i = 0; i < M_N; ++i )
                {
                    M_Mqm_pr[q][m]( i, j ) = M_model->Mqm(q, m,
                                                          M_model->rBFunctionSpace()->primalBasisElement(i),
                                                          M_model->rBFunctionSpace()->primalBasisElement(j) );
                    M_Mqm_du[q][m]( i, j ) = M_model->Mqm(q, m,
                                                          M_model->rBFunctionSpace()->dualBasisElement(i),
                                                          M_model->rBFunctionSpace()->dualBasisElement(j), true );
                    M_Mqm_pr_du[q][m]( i, j ) = M_model->Mqm(q, m,
                                                             M_model->rBFunctionSpace()->dualBasisElement(i),
                                                             M_model->rBFunctionSpace()->primalBasisElement(j) );
                }
            }
        }//loop over m
    }//loop over q

    LOG(INFO) << "[CRB::offline] compute Lq_pr, Lq_du" << "\n";

    for ( size_type q = 0; q < M_model->Ql( M_output_index ); ++q )
    {
        for( size_type m = 0; m < M_model->mMaxF( M_output_index, q ); ++m )
        {
            M_Lqm_pr[q][m].conservativeResize( M_N );
            M_Lqm_du[q][m].conservativeResize( M_N );

            for ( size_type l = 1; l <= number_of_elements_to_update; ++l )
            {
                int index = M_N-l;
                M_Lqm_pr[q][m]( index ) = M_model->Fqm( M_output_index, q, m, M_model->rBFunctionSpace()->primalBasisElement(index) );
                M_Lqm_du[q][m]( index ) = M_model->Fqm( M_output_index, q, m, M_model->rBFunctionSpace()->dualBasisElement(index) );
            }
        }//loop over m
    }//loop over q

    LOG(INFO) << "compute coefficients needed for the initialization of unknown in the online step\n";
    if ( model_type::is_time_dependent || !M_model->isSteady() )
    {
        element_ptrtype primal_initial_field ( new element_type ( M_model->functionSpace() ) );
        M_model->initializationField( primal_initial_field, mu ); //fill initial_field
        M_coeff_pr_ini_online.conservativeResize( M_N );
        if ( M_orthonormalize_primal )
        {
            for ( size_type elem=M_N - number_of_elements_to_update; elem<M_N; elem++ )
            {
                //primal
                double k =  M_model->scalarProduct( *primal_initial_field, M_model->rBFunctionSpace()->primalBasisElement(elem) );
                M_coeff_pr_ini_online(elem)= k ;
            }
        }

        else if ( !M_orthonormalize_primal )
        {
            matrixN_type MN ( ( int )M_N, ( int )M_N ) ;
            vectorN_type FN ( ( int )M_N );

            //primal
            for ( size_type i=0; i<M_N; i++ )
            {
                for ( size_type j=0; j<i; j++ )
                {
                    MN( i,j ) = M_model->scalarProduct( M_model->rBFunctionSpace()->primalBasisElement(j),
                                                        M_model->rBFunctionSpace()->primalBasisElement(i) );
                    MN( j,i ) = MN( i,j );
                }

                MN( i,i ) = M_model->scalarProduct( M_model->rBFunctionSpace()->primalBasisElement(i),
                                                    M_model->rBFunctionSpace()->primalBasisElement(i) );
                FN( i ) = M_model->scalarProduct( *primal_initial_field, M_model->rBFunctionSpace()->primalBasisElement(i) );
            }

            vectorN_type projectionN ( ( int ) M_N );
            projectionN = MN.lu().solve( FN );

            for ( size_type i=M_N - number_of_elements_to_update; i<M_N; i++ )
            {
                M_coeff_pr_ini_online(i)= projectionN( i ) ;
            }
        }

        if ( M_solve_dual_problem )
        {
            M_coeff_du_ini_online.conservativeResize( M_N );

            if ( M_orthonormalize_dual )
            {
                for ( size_type elem=M_N - number_of_elements_to_update; elem<M_N; elem++ )
                {
                    double k =  M_model->scalarProduct( *dual_initial_field, M_model->rBFunctionSpace()->dualBasisElement(elem) );
                    M_coeff_du_ini_online(elem)= k ;
                }
            }
            else if ( !M_orthonormalize_dual )
            {
                matrixN_type MNdu ( ( int )M_N, ( int )M_N ) ;
                vectorN_type FNdu ( ( int )M_N );

                //dual
                for ( size_type i=0; i<M_N; i++ )
                {
                    for ( size_type j=0; j<i; j++ )
                    {
                        MNdu( i,j ) = M_model->scalarProduct( M_model->rBFunctionSpace()->dualBasisElement(j),
                                                              M_model->rBFunctionSpace()->dualBasisElement(i) );
                        MNdu( j,i ) = MNdu( i,j );
                    }

                    MNdu( i,i ) = M_model->scalarProduct( M_model->rBFunctionSpace()->dualBasisElement(i),
                                                          M_model->rBFunctionSpace()->dualBasisElement(i) );
                    FNdu( i ) = M_model->scalarProduct( *dual_initial_field, M_model->rBFunctionSpace()->dualBasisElement(i) );
                }

                vectorN_type projectionN ( ( int ) M_N );
                projectionN = MNdu.lu().solve( FNdu );

                for ( size_type i=M_N - number_of_elements_to_update; i<M_N; i++ )
                {
                    M_coeff_du_ini_online(i)= projectionN( i ) ;
                }
            }
        }
    }

    if ( true )
    {
        M_algebraicInnerProductPrimal.conservativeResize( M_N, M_N );
        for ( size_type i=M_N - number_of_elements_to_update; i<M_N; i++ )
        {
            for ( size_type j=0; j<M_N; j++ )
            {
                M_algebraicInnerProductPrimal(i,j) = inner_product( M_model->rBFunctionSpace()->primalBasisElement(i),
                                                                    M_model->rBFunctionSpace()->primalBasisElement(j) );
                if ( i!=j )
                    M_algebraicInnerProductPrimal(j,i) = M_algebraicInnerProductPrimal(i,j);
            }
        }
    }

    toc("Projection on reduced basis space");
}

template<typename TruthModelType>
void
CRB<TruthModelType>::saveRB()
{
    M_elements_database.setWn( boost::make_tuple( M_model->rBFunctionSpace()->primalRB() , M_model->rBFunctionSpace()->dualRB() ) );

    M_elements_database.saveNewElementToDB();
}


template<typename TruthModelType>
bool
CRB<TruthModelType>::showMuSelection()
{
    bool show = boption(_prefix=M_prefix,_name="crb.show-mu-selection");
    return show;
}

template<typename TruthModelType>
void
CRB<TruthModelType>::loadDB( std::string const& filename, crb::load l )
{
    auto fname = this->db( filename );

    if ( ( l == crb::load::all ) ||  (l == crb::load::fe ) )
        this->setLoadBasisFromDB( true );
    else
        this->setLoadBasisFromDB( false );
    this->loadJson( fname.string() );

    LOG(INFO) << "Loaded DB CRB " << fname;
}




template<typename TruthModelType>
void
CRB<TruthModelType>::loadJson( std::string const& filename, size_type loadingContext )
{
    // first load the model
    M_model->loadJson( filename, "crbmodel" );

    if ( !fs::exists( filename ) )
    {
        LOG(INFO) << "Could not find " << filename << std::endl;
        return;
    }

    std::string dbDir = fs::path( filename ).parent_path().string();

    auto json_str_wo_comments = removeComments(readFromFile(filename));
    //LOG(INFO) << "json file without comment:" << json_str_wo_comments;

    boost::property_tree::ptree ptree;
    std::istringstream istr( json_str_wo_comments );
    boost::property_tree::read_json( istr, ptree );
    this->setup( ptree, loadingContext, dbDir );
}
template<typename TruthModelType>
void
CRB<TruthModelType>::saveJson()
{
    if ( this->worldComm().isMasterRank() )
    {
        //std::string filenameJson = (this->dbLocalPath()/fs::path("crb.json")).string();
        std::string filenameJson = (this->dbLocalPath()/fs::path(this->jsonFilename())).string();
        std::cout << "saveDB: " << filenameJson << std::endl;

        boost::property_tree::ptree ptree;

        ptree.add( "uuid", this->idStr() );

        boost::property_tree::ptree ptreeCrbModel;
        M_model->updatePropertyTree( ptreeCrbModel );
        ptree.add_child( "crbmodel", ptreeCrbModel );

        boost::property_tree::ptree ptreeReducedBasisSpace;
        std::string meshFilename = fmt::format("{}_mesh_p{}.json",this->name(),this->worldComm().size());
        ptreeReducedBasisSpace.add( "mesh-filename",meshFilename );
        ptreeReducedBasisSpace.add( "database-filename", M_elements_database.dbFilename() );
        // add functionspace file (useful when nproc offline != nproc online)
        std::string functionspaceFilename = fmt::format("{}_functionSpace_p{}.json", this->name(), this->worldComm().size());
        ptreeReducedBasisSpace.add( "functionspace-filename", functionspaceFilename );
        
        ptreeReducedBasisSpace.add( "dimension", M_N );
        if ( M_model && M_model->rBFunctionSpace() && M_model->rBFunctionSpace()->functionSpace() )
        {
            ptreeReducedBasisSpace.add( "mesh-context",M_model->rBFunctionSpace()->functionSpace()->mesh()->components().context() );

            auto feSpace = M_model->rBFunctionSpace()->functionSpace();
            boost::property_tree::ptree ptreeFiniteElementSpace;
            ptreeFiniteElementSpace.add( "dimension", feSpace->nDof() );
            ptreeFiniteElementSpace.add( "basis-name", feSpace->basisName() );
            ptreeReducedBasisSpace.add_child( "finite-element-space", ptreeFiniteElementSpace );
        }
        ptree.add_child( "reduced-basis-space", ptreeReducedBasisSpace );

        boost::property_tree::ptree ptreeCrb;//Database;
        ptreeCrb.add( "dimension", M_N );
        ptreeCrb.add( "name", this->name() );
        // ptreeCrb.add( "database-filename",(this->dbLocalPath() / this->dbFilename()).string() );
        ptreeCrb.add( "database-filename", this->dbFilename() );
        ptreeCrb.add( "has-solve-dual-problem",M_solve_dual_problem );
        ptreeCrb.add( "output-index", M_output_index );
        ptreeCrb.add( "output-name", M_output_name );
        ptreeCrb.add( "error-type", M_error_type );
        ptree.add_child( "crb", ptreeCrb );

        if ( M_error_type == CRBErrorType::CRB_RESIDUAL_SCM )
        {
            boost::property_tree::ptree ptreeParameterScmA;
            M_scmA->updatePropertyTree( ptreeParameterScmA );
            ptree.add_child( "scmA", ptreeParameterScmA );
            if ( !M_model->isSteady() )
            {
                boost::property_tree::ptree ptreeParameterScmM;
                M_scmM->updatePropertyTree( ptreeParameterScmM );
                ptree.add_child( "scmM", ptreeParameterScmM );
            }
        }        
        write_json( filenameJson, ptree );
    }
}

template<typename TruthModelType>
void
CRB<TruthModelType>::setup( boost::property_tree::ptree const& ptree, size_type loadingContext, std::string const& dbDir )
{
    auto i = ptree.template get<std::string>( "uuid" );
    this->setId( boost::lexical_cast<uuids::uuid>( i ) );
    auto const& ptreeCrb = ptree.get_child( "crb" );
    M_N = ptreeCrb.template get<int>( "dimension" );
    M_output_index = ptreeCrb.template get<int>( "output-index" );
    int et = ptreeCrb.template get<int>( "error-type" );
    M_error_type = static_cast<CRBErrorType>( et );
    this->setName( ptreeCrb.template get<std::string>( "name" ) );
    M_solve_dual_problem = ptreeCrb.template get<bool>( "has-solve-dual-problem" );
    std::string dbname = ptreeCrb.template get<std::string>( "database-filename" );
    fs::path dbnamePath = fs::path( dbname );
    this->setDBFilename( dbnamePath.filename().string() );
    if ( dbnamePath.is_absolute() )
        this->setDBDirectory( fs::absolute( dbnamePath ).parent_path().string() );
    else if ( !dbDir.empty() )
        this->setDBDirectory( dbDir );
    this->setIsLoaded( false );
    CHECK( this->loadDB() ) << "crb load fails";
    if ( M_error_type == CRBErrorType::CRB_RESIDUAL_SCM )
    {
        auto const& ptreeScmA = ptree.get_child( "scmA" );
        M_scmA->setup( ptreeScmA, dbDir );
        auto ptreeOptionalScmM = ptree.get_child_optional( "scmM" );
        if ( ptreeOptionalScmM )
            M_scmM->setup( *ptreeOptionalScmM, dbDir );
    }

    auto const& ptreeReducedBasisSpace = ptree.get_child( "reduced-basis-space" );
    if ( (loadingContext == 1 || M_loadElementsDb ) && M_model )
    {
        // M_model->rBFunctionSpace()->setModel( M_model->model() );
        // //M_model->rBFunctionSpace()->setup( ptreeReducedBasisSpace );
        // M_elements_database.setModel( M_model );
        M_elements_database.setup( ptreeReducedBasisSpace, dbDir );
        M_loadElementsDb = true;
    }
    else
    {
        int rbdim = ptreeReducedBasisSpace.template get<size_type>( "dimension" );
        M_model->rBFunctionSpace()->setDimension( rbdim );
    }

}

template<typename TruthModelType>
void
CRB<TruthModelType>::saveDB()
{

    if ( this->worldComm().isMasterRank() )
    {
        fs::ofstream ofs( this->dbLocalPath() / this->dbFilename() );
        if ( ofs )
        {
            // boost::archive::text_oarchive oa( ofs );
            boost::archive::binary_oarchive oa( ofs );
            // write class instance to archive
            oa << *this;
            // archive and stream closed when destructors are called
        }
    }
    saveJson();
}

template<typename TruthModelType>
bool
CRB<TruthModelType>::loadDB()
{
    //if ( this->rebuildDB() )
    //    return false;

    if( this->isDBLoaded() )
    {
        return true;
    }

    fs::path db = this->lookForDB();

    if ( db.empty() )
    {
        return false;
    }
    if ( !fs::exists( db ) )
    {
        return false;
    }
    //std::cout << "Loading " << db << "...\n";
    fs::ifstream ifs( db );

    if ( ifs )
    {
        //boost::archive::text_iarchive ia( ifs );
        boost::archive::binary_iarchive ia( ifs );
        // write class instance to archive
        ia >> *this;
        //std::cout << "Loading " << db << " done...\n";
        this->setIsLoaded( true );
        // archive and stream closed when destructors are called
        return true;
    }

    return false;
}

template<typename TruthModelType>
void
CRB<TruthModelType>::setupOfflineFromDB()
{
    this->loadDB( (this->dbLocalPath()/fs::path(this->jsonFilename())).string(), crb::load::all );

    if( this->worldComm().isMasterRank() )
        std::cout << "Database CRB available and loaded with " << M_N <<" basis\n";
    LOG(INFO) << "Database CRB available and loaded with " << M_N <<" basis";

    M_elements_database.setMN( M_N );
    if( M_loadElementsDb )
    {
        if( M_elements_database.loadDB() )
        {
            if( this->worldComm().isMasterRank() )
            {
                if(M_elements_database.dbFileFormat()=="boost")
                {
                    std::cout<<"Database for basis functions " << M_elements_database.lookForDB() << " available and loaded\n";
                    LOG(INFO) << "Database for basis functions " << M_elements_database.lookForDB() << " available and loaded";
                }
                else
                {
                    std::cout<<"Database for basis functions in HDF5 available and loaded\n";
                    LOG(INFO) << "Database for basis functions in HDF5 available and loaded";
                }
            }
                
            auto basis_functions = M_elements_database.wn();
            M_model->rBFunctionSpace()->setBasis( basis_functions );
        }
        else
            M_N = 0;
    }
    if ( M_error_type == CRBErrorType::CRB_RESIDUAL_SCM )
    {
        if ( M_scmA->loadDB() )
        {
            if( this->worldComm().isMasterRank() )
                std::cout << "Database for SCM_A " << M_scmA->lookForDB() << " available and loaded\n";
            LOG( INFO ) << "Database for SCM_A " << M_scmA->lookForDB() << " available and loaded";
        }
        else
            M_N = 0;

        if ( !M_model->isSteady() )
        {
            if ( M_scmM->loadDB() )
            {
                if( this->worldComm().isMasterRank() )
                    std::cout << "Database for SCM_M " << M_scmM->lookForDB() << " available and loaded";
                LOG( INFO ) << "Database for SCM_M " << M_scmM->lookForDB() << " available and loaded";
            }
            else
                M_N = 0;
        }
    }

}

} // Feel

namespace boost
{
namespace serialization
{
template< typename T>
struct version< Feel::CRB<T> >
{
    // at the moment the version of the CRB DB is 1. if any changes is done
    // to the format it is mandatory to increase the version number below
    // and use the new version number of identify the new entries in the DB
    typedef mpl::int_<1> type;
    typedef mpl::integral_c_tag tag;
    static const unsigned int value = version::type::value;
};
template<typename T> const unsigned int version<Feel::CRB<T> >::value;
}
}
#endif /* __CRB_H */
