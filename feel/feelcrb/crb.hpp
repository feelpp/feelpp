/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-11-24

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
   \file crb.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-11-24
 */
#ifndef __CRB_H
#define __CRB_H 1

#include <boost/multi_array.hpp>
#include <boost/tuple/tuple.hpp>
#include "boost/tuple/tuple_io.hpp"
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/bimap.hpp>
#include <boost/bimap/support/lambda.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
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
#include <feel/feelcrb/pod.hpp>
#include <feel/feeldiscr/bdf2.hpp>


namespace Feel
{
  /**
   * CRBErrorType
   * Determine the type of error estimation used
   * - CRB_RESIDUAL : use the residual error estimation without algorithm SCM
   * - CRB_RESIDUAL_SCM : use the residual error estimation and also algorithm SCM
   * - CRB_NO_RESIDUAL : in this case we don't compute error estimation
   * - CRB_EMPIRICAL : compute |S_n - S_{n-1}| where S_n is the output obtained by using a reduced basis with n elements
   */
  enum CRBErrorType { CRB_RESIDUAL = 0, CRB_RESIDUAL_SCM=1, CRB_NO_RESIDUAL=2 , CRB_EMPIRICAL=3};

/**
 * \class CRB
 * \brief Certifed Reduced Basis class
 *
 * Implements the certified reduced basis method
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename TruthModelType>
class CRB : public CRBDB
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

    typedef double value_type;
    typedef boost::tuple<double,double> bounds_type;

    typedef ParameterSpace<TruthModelType::ParameterSpaceDimension> parameterspace_type;
    typedef boost::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef typename parameterspace_type::element_type parameter_type;
    typedef typename parameterspace_type::element_ptrtype parameter_ptrtype;
    typedef typename parameterspace_type::sampling_type sampling_type;
    typedef typename parameterspace_type::sampling_ptrtype sampling_ptrtype;

    typedef boost::tuple<double, parameter_type, size_type> relative_error_type;
    typedef relative_error_type max_error_type;

    typedef boost::bimap< int, double > convergence_type;
    typedef typename convergence_type::value_type convergence;


    //! scm
    typedef CRBSCM<truth_model_type> scm_type;
    typedef boost::shared_ptr<scm_type> scm_ptrtype;

    //! POD
    typedef POD<truth_model_type> pod_type;
    typedef boost::shared_ptr<pod_type> pod_ptrtype;



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
    typedef typename model_type::theta_vector_type theta_vector_type;


    typedef Eigen::VectorXd y_type;
    typedef std::vector<y_type> y_set_type;
    typedef std::vector<boost::tuple<double,double> > y_bounds_type;

    typedef std::vector<element_type> wn_type;

    typedef std::vector<double> vector_double_type;
    typedef boost::shared_ptr<vector_double_type> vector_double_ptrtype;

    typedef Eigen::VectorXd vectorN_type;
    typedef Eigen::MatrixXd matrixN_type;



    typedef std::vector<element_type> mode_set_type;
    typedef boost::shared_ptr<mode_set_type> mode_set_ptrtype;

    typedef boost::multi_array<value_type, 2> array_2_type;
    typedef boost::multi_array<vectorN_type, 2> array_3_type;
    typedef boost::multi_array<matrixN_type, 2> array_4_type;

    //! mesh type
    typedef typename model_type::mesh_type mesh_type;

    //! space type
    typedef typename model_type::space_type space_type;

    //! time discretization
	typedef Bdf<space_type>  bdf_type;
	typedef boost::shared_ptr<bdf_type> bdf_ptrtype;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    CRB()
        :
        super( "noname" ),
        M_model(),
        M_output_index( 0 ),
        M_tolerance( 1e-2),
        M_iter_max( 3 ),
        M_factor( -1 ),
        M_error_type( CRB_NO_RESIDUAL ),
        M_Dmu( new parameterspace_type ),
        M_Xi( new sampling_type( M_Dmu ) ),
        M_WNmu( new sampling_type( M_Dmu, 1, M_Xi ) ),
        M_WNmu_complement()
        {

        }

    //! constructor from command line options
    CRB( std::string  name,
         po::variables_map const& vm )
        :
        super( (boost::format( "%1%" ) % vm["crb.error-type"].template as<int>()).str(),
               name,
               (boost::format( "%1%-%2%-%3%" ) % name % vm["crb.output-index"].template as<int>() % vm["crb.error-type"].template as<int>() ).str(),
               vm ),

        M_model(),
        M_backend( backend_type::build( vm ) ),
        M_output_index( vm["crb.output-index"].template as<int>() ),
        M_tolerance( vm["crb.error-max"].template as<double>() ),
        M_iter_max( vm["crb.dimension-max"].template as<int>() ),
        M_factor(vm["crb.factor"].template as<int>() ),
        M_error_type( CRBErrorType(vm["crb.error-type"].template as<int>() ) ),
        M_Dmu( new parameterspace_type ),
        M_Xi( new sampling_type( M_Dmu ) ),
        M_WNmu( new sampling_type( M_Dmu, 1, M_Xi ) ),
        M_WNmu_complement(),
        M_scm( new scm_type( name, vm ) )
        {
            if ( this->loadDB() )
                std::cout << "Database " << this->lookForDB() << " available and loaded\n";


        }

    //! copy constructor
    CRB( CRB const & o )
        :
        super( o ),
        M_output_index( o.M_output_index ),
        M_tolerance( o.M_tolerance ),
        M_iter_max( o.M_iter_max ),
        M_factor( o.M_factor ),
        M_error_type( o.M_error_type ),
        M_Dmu( o.M_Dmu ),
        M_Xi( o.M_Xi ),
        M_WNmu( o.M_WNmu ),
        M_WNmu_complement( o.M_WNmu_complement ),
        M_C0_pr( o.M_C0_pr ),
        M_C0_du( o.M_C0_du ),
        M_Lambda_pr( o.M_Lambda_pr ),
        M_Lambda_du( o.M_Lambda_du ),
        M_Gamma_pr( o.M_Gamma_pr ),
        M_Gamma_du( o.M_Gamma_du )
        {}

    //! destructor
    ~CRB()
        {}


    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    CRB& operator=( CRB const & o)
        {
            if (this != &o )
            {
            }
            return *this;
        }
    //@}

    /** @name Accessors
     */
    //@{


    //! return factor
    int factor() const {return factor; }
    //! \return max iterations
    int maxIter() const { return M_iter_max; }

    //! \return the parameter space
    parameterspace_ptrtype Dmu() const { return M_Dmu; }

    //! \return the output index
    int outputIndex() const { return M_output_index; }

    //! \return the dimension of the reduced basis space
    int dimension() const { return M_N; }

    //! \return the train sampling used to generate the reduced basis space
    sampling_ptrtype trainSampling() const { return M_Xi; }

    //! \return the error type
    CRBErrorType errorType() const { return M_error_type; }

    //! \return the scm object
    scm_ptrtype scm() const { return M_scm; }


    //@}

    /** @name  Mutators
     */
    //@{

    //! set the output index
    void setOutputIndex( uint16_type oindex )
        {
            int M_prev_o = M_output_index;
            M_output_index = oindex;
            if ( M_output_index >= M_model->Nl()  )
                M_output_index = M_prev_o;
            //std::cout << " -- crb set output index to " << M_output_index << " (max output = " << M_model->Nl() << ")\n";
            this->setDBFilename( (boost::format( "%1%-%2%-%3%.crbdb" ) % this->name() % M_output_index % M_error_type ).str() );
            if ( M_output_index != M_prev_o )
                this->loadDB();
            //std::cout << "Database " << this->lookForDB() << " available and loaded\n";

        }

    //! set the crb error type
    void setCRBErrorType( CRBErrorType error ) { M_error_type = error; }

    //! set offline tolerance
    void setTolerance( double tolerance ) { M_tolerance = tolerance; }

    //! set the truth offline model
    void setTruthModel( truth_model_ptrtype const& model )
        {
            M_model = model;
            M_Dmu = M_model->parameterSpace();
            M_Xi = sampling_ptrtype( new sampling_type( M_Dmu ) );
            M_WNmu = sampling_ptrtype( new sampling_type( M_Dmu ) );

            M_scm->setTruthModel( M_model );
        }

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


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * generate offline the space \f$W_N\f$
     */
    void generate();

    /**
     * orthonormalize the basis
     */
    void orthonormalize( size_type N, wn_type& wn, int Nm = 1 );


    /*
     * check orthonormality
     */
    //void checkOrthonormality( int N, const wn_type& wn) const;
    void checkOrthonormality( int N, const wn_type& wn) const;

    /**
     * check the reduced basis space invariant properties
     * \param N dimension of \f$W_N\f$
     */
    void check( size_type N )  const;

    /**
     * compute effectivity indicator of the error estimation overall a given
     * parameter space
     *
     * \param max_ei : maximum efficiency indicator (output)
     * \param min_ei : minimum efficiency indicator (output)
     * \param Dmu (input) parameter space
     * \param N : sampling size (optional input with default value)
     */
    void computeErrorEstimationEfficiencyIndicator (parameterspace_ptrtype const& Dmu, double& max_ei, double& min_ei,int N);

    /**
     * Returns the lower bound of the output
     *
     * \param mu \f$ \mu\f$ the parameter at which to evaluate the output
     * \param N the size of the reduced basis space to use
     * \param uN primal solution
     * \param uNdu dual solution
     * \param K : index of time ( time = K*dt) at which we want to evaluate the output
     * Note : K as a default value for non time-dependent problems
     *
     *\return compute online the lower bound
     */
    value_type lb( size_type N, parameter_type const& mu, vectorN_type& uN, vectorN_type& uNdu , int K=0) const;
    value_type lb( size_type N, parameter_type const& mu, vectorN_type& uN, vectorN_type& uNdu, mpl::bool_<true> , int K=0 ) const;
    value_type lb( size_type N, parameter_type const& mu, vectorN_type& uN, vectorN_type& uNdu, mpl::bool_<false>, int K=0 ) const;

    /**
     * Returns the lower bound of the output
     *
     * \param mu \f$ \mu\f$ the parameter at which to evaluate the output
     * \param N the size of the reduced basis space to use
     * \param uN primal solution
     * \param uNdu dual solution
     *
     *\return compute online the lower bound
     */
    value_type lb( parameter_ptrtype const& mu, size_type N, vectorN_type& uN, vectorN_type& uNdu ) const
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
    value_type ub( size_type N, parameter_type const& mu, vectorN_type& uN, vectorN_type& uNdu ) const
        {
            return lb( N, mu, uN, uNdu ) + delta( N, mu, uN, uNdu );
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
    value_type delta( size_type N, parameter_ptrtype const& mu, vectorN_type const& uN, vectorN_type const& uNdu, int k=0 ) const
        {
            return delta( N, *mu, uN, uNdu );
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
    value_type delta( size_type N, parameter_type const& mu, vectorN_type const& uN, vectorN_type const& uNdu, int k=0) const;

    /**
     * Returns the upper bound of the output associed to \f$\mu\f$
     *
     * \param mu \f$ \mu\f$ the parameter at which to evaluate the output
     * \param N the dimension of \f$W_N\f$
     *
     *\return compute online the lower bound
     */
    value_type ub( size_type K, parameter_ptrtype const& mu, vectorN_type& uN, vectorN_type& uNdu ) const
        {
            return ub( K, *mu, uN, uNdu );
        }

    /**
     * Offline computation
     *
     * \return the convergence history (max error)
     */
    convergence_type offline();


    /**
     * offline computation with or without error estimation (called by offline function)
     *
     */
    void offlineNoErrorEstimation();
    void offlineNoErrorEstimation(mpl::bool_<true>);
    void offlineNoErrorEstimation(mpl::bool_<false>);
    void offlineWithErrorEstimation();
    void offlineWithErrorEstimation(mpl::bool_<true>);
    void offlineWithErrorEstimation(mpl::bool_<false>);



    /**
     * \brief Retuns maximum value of the relative error
     * \param N number of elements in the reduced basis <=> M_N
     */
    max_error_type maxErrorBounds( size_type N ) const;

    /**
     * evaluate online the residual
     */
    value_type N2Q2( int Ncur, parameter_type const& mu, vectorN_type const& Un, vectorN_type const& Undu ) const;
    value_type N2Q2( int Ncur, parameter_type const& mu, vectorN_type const& Un, vectorN_type const& Undu, mpl::bool_<true> ) const;
    value_type N2Q2( int Ncur, parameter_type const& mu, vectorN_type const& Un, vectorN_type const& Undu, mpl::bool_<false> ) const;

    /**
     * generate offline the residual
     */
    void generateN2Q2( int Ncur );
    void generateN2Q2( int Ncur, mpl::bool_<true> );
    void generateN2Q2( int Ncur, mpl::bool_<false> );

    /*
     * compute empirical error estimation, ie : |S_n - S{n-1}|
     * \param Ncur : number of elements in the reduced basis <=> M_N
     * \param mu : parameters value (one value per parameter)
     * \output : empirical error estimation
     */
    value_type empiricalErrorEstimation ( int Ncur, parameter_type const& mu, int k) const ;


    /**
     * run the certified reduced basis with P parameters and returns 1 output
     */
    boost::tuple<double,double,double> run( parameter_type const& mu, double eps = 1e-6 );

    /**
     * run the certified reduced basis with P parameters and returns 1 output
     */
    void run( const double * X, unsigned long N, double * Y, unsigned long P );

    /**
     * \return a random sampling
     */
    sampling_type randomSampling( int N ) { M_Xi->randomize(N); return *M_Xi; }

    /**
     * \return a equidistributed sampling
     */
    sampling_type equidistributedSampling( int N ) { M_Xi->equidistribute(N); return *M_Xi; }

    /**
     * save the CRB database
     */
    void saveDB();

    /**
     * load the CRB database
     */
    bool loadDB();

    /**
     *  do the projection on the POD space of u (for transient problems)
     *  \param u : the solution to project (input parameter)
     *  \param projection : the projection (output parameter)
     *  \param name_of_space : primal or dual
     */
    void projectionOnPodSpace(const element_ptrtype & u , element_ptrtype& projection ,const std::string& name_of_space="primal" );

    //@}

private:
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const;

    template<class Archive>
    void load(Archive & ar, const unsigned int version) ;

    BOOST_SERIALIZATION_SPLIT_MEMBER()
    private:

    truth_model_ptrtype M_model;

    int M_factor;

    backend_ptrtype M_backend;

    int M_output_index;
    double M_tolerance;
    int M_iter_max;


    CRBErrorType M_error_type;

    // parameter space
    parameterspace_ptrtype M_Dmu;

    // fine sampling of the parameter space
    sampling_ptrtype M_Xi;

    // sampling of parameter space to build WN
    sampling_ptrtype M_WNmu;
    sampling_ptrtype M_WNmu_complement;

    // reduced basis space
    wn_type M_WN;
    wn_type M_WNdu;


    size_type M_N;
    size_type M_Nm;

    bool orthonormalize_primal;
    bool orthonormalize_dual;

    convergence_type M_rbconv;

    //scm
    scm_ptrtype M_scm;


    //time
    bdf_ptrtype M_bdf_primal;
    bdf_ptrtype M_bdf_primal_save;
    bdf_ptrtype M_bdf_dual;

    // left hand side
    std::vector<matrixN_type> M_Aq_pr;
    std::vector<matrixN_type> M_Aq_du;
    std::vector<matrixN_type> M_Aq_pr_du;

    //mass matrix
    std::vector<matrixN_type> M_Mq_pr;
    std::vector<matrixN_type> M_Mq_du;
    std::vector<matrixN_type> M_Mq_pr_du;

    // right hand side
    std::vector<vectorN_type> M_Fq_pr;
    std::vector<vectorN_type> M_Fq_du;
    // output
    std::vector<vectorN_type> M_Lq_pr;
    std::vector<vectorN_type> M_Lq_du;

    array_2_type M_C0_pr;
    array_2_type M_C0_du;
    array_3_type M_Lambda_pr;
    array_3_type M_Lambda_du;
    array_4_type M_Gamma_pr;
    array_4_type M_Gamma_du;


};

po::options_description crbOptions( std::string const& prefix = "" );



template<typename TruthModelType>
typename CRB<TruthModelType>::convergence_type
CRB<TruthModelType>::offline()

{
    boost::timer ti;
    std::cout << "Offline CRB starts, this may take a while until Database is computed...\n";
    Log() << "[CRB::offline] Starting offline for output " << M_output_index << "\n";
    Log() << "[CRB::offline] initialize underlying finite element model\n";
    M_model->init();
    std::cout << " -- model init done in " << ti.elapsed() << "s\n"; ti.restart();

    orthonormalize_primal = this->vm()["crb.orthonormalize_primal"].template as<bool>() ;
    orthonormalize_dual = this->vm()["crb.orthonormalize_dual"].template as<bool>() ;
    M_Nm = this->vm()["crb.Nm"].template as<int>() ;


    //scm_ptrtype M_scm = scm_ptrtype( new scm_type( M_vm ) );
    //M_scm->setTruthModel( M_model );
    //    std::vector<boost::tuple<double,double,double> > M_rbconv2 = M_scm->offline();

    Log() << "[CRB::offline] compute random sampling\n";
    // random sampling
    M_Xi->randomize( this->vm()["crb.sampling-size"].template as<int>() );
    //M_Xi->equidistribute( this->vm()["crb.sampling-size"].template as<int>() );
    M_WNmu->setSuperSampling( M_Xi );

    std::cout << " -- sampling init done in " << ti.elapsed() << "s\n"; ti.restart();
    if ( M_error_type == CRB_RESIDUAL || M_error_type == CRB_RESIDUAL_SCM )
    {
        int __QLhs = M_model->Qa();
        int __QRhs = M_model->Ql(0);
        int __QOutput = M_model->Ql(M_output_index);

        typename array_2_type::extent_gen extents2;
        M_C0_pr.resize( extents2[__QRhs][__QRhs] );
        M_C0_du.resize( extents2[__QOutput][__QOutput] );

        typename array_3_type::extent_gen extents3;
        M_Lambda_pr.resize( extents3[__QLhs][__QRhs] );
        M_Lambda_du.resize( extents3[__QLhs][__QOutput] );

        typename array_4_type::extent_gen extents4;
        M_Gamma_pr.resize( extents4[__QLhs][__QLhs] );
        M_Gamma_du.resize( extents4[__QLhs][__QLhs] );
        std::cout << " -- residual data init done in " << ti.elapsed() << "\n"; ti.restart();
    }

    if(M_error_type == CRB_NO_RESIDUAL)
        offlineNoErrorEstimation();
    else
        offlineWithErrorEstimation();

    return M_rbconv;
}



template<typename TruthModelType>
void
CRB<TruthModelType>::offlineNoErrorEstimation()
{
    return offlineNoErrorEstimation(mpl::bool_<model_type::is_time_dependent>());
}


template<typename TruthModelType>
void
CRB<TruthModelType>::offlineNoErrorEstimation(mpl::bool_<true>)
{

    parameter_type mu( M_Dmu );

    double relative_error = 1e30;

    // empty sets
    M_WNmu->clear();

    size_type index;

    sparse_matrix_ptrtype M,A,Adu;
    std::vector<vector_ptrtype> F,L;


    const int Ndof = M_model->functionSpace()->nDof();//number of dofs used

    // dimension of reduced basis space
    M_N = 0;

    size_type Np = 1;

    Log() << "[CRB::offlineNoErrorEstimation] compute affine decomposition\n";
    std::vector<sparse_matrix_ptrtype> Aq;
    std::vector<sparse_matrix_ptrtype> Mq;
    std::vector<std::vector<vector_ptrtype> > Fq,Lq;
    boost::tie( Mq, Aq, Fq ) = M_model->computeAffineDecomposition();

    Log() << "[CRB::offlineNoErrorEstimation] allocate reduced basis data structures\n";
    M_Aq_pr.resize( M_model->Qa() );
    M_Aq_du.resize( M_model->Qa() );
    M_Aq_pr_du.resize( M_model->Qa() );

    M_Mq_pr.resize( M_model->Qm() );
    M_Mq_du.resize( M_model->Qm() );
    M_Mq_pr_du.resize( M_model->Qm() );

    M_Fq_pr.resize( M_model->Ql(0) );
    M_Fq_du.resize( M_model->Ql(0) );
    M_Lq_pr.resize( M_model->Ql(M_output_index) );
    M_Lq_du.resize( M_model->Ql(M_output_index) );

    element_ptrtype u( new element_type( M_model->functionSpace() ) );
    element_ptrtype uproj( new element_type( M_model->functionSpace() ) );
    element_ptrtype udu( new element_type( M_model->functionSpace() ) );
    vector_ptrtype U( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype Rhs( M_backend->newVector( M_model->functionSpace() ) );

    sampling_ptrtype Sampling;
    Sampling = sampling_ptrtype( new sampling_type( M_Dmu ) );

    int sampling_size = M_iter_max / M_Nm ;
    Sampling->logEquidistribute( sampling_size );
    Log() << "[CRB::offlineNoErrorEstimation] sampling parameter space with "<< sampling_size <<"elements done\n";

    Log() << "[CRB::offlineNoErrorEstimation] strategy "<< M_error_type <<"\n";
    std::cout << "[CRB::offlineNoErrorEstimation] strategy "<< M_error_type <<"\n";

    bool reuse_prec = this->vm()["crb.reuse_prec"].template as<bool>() ;

    for( size_type index = 0; index < Sampling->size(); ++index )
    {
        parameter_type const& mu = M_Xi->at( index );

        Log() <<"========================================"<<"\n";
        std::cout << "============================================================\n";
        std::cout << "N=" << M_N << "\n";
        Log() << "N=" << M_N << "\n";

        backend_ptrtype backend_primal_problem = backend_type::build( BACKEND_PETSC );
        backend_ptrtype backend_dual_problem = backend_type::build( BACKEND_PETSC );

        // for a given parameter \p mu assemble the left and right hand side

        u->setName( (boost::format( "fem-primal-%1%" ) % (M_N-1)).str() );

        M_bdf_primal = bdf( _space=M_model->functionSpace(), _vm=this->vm() , _name="bdf_primal");
        M_bdf_primal_save = bdf( _space=M_model->functionSpace(), _vm=this->vm() , _name="bdf_primal_save");


        //set parameters for time discretization
        if( M_model->isSteady() )
        {
            M_bdf_primal->setSteady();
            M_bdf_primal_save->setSteady();
        }
        else
        {
            M_bdf_primal->setTimeInitial( M_model->timeInitial() );
            M_bdf_primal->setTimeStep( M_model->timeStep() );
            M_bdf_primal->setTimeFinal( M_model->timeFinal() );
            M_bdf_primal->setOrder( M_model->timeOrder() );

            M_bdf_primal_save->setTimeInitial( M_model->timeInitial() );
            M_bdf_primal_save->setTimeStep( M_model->timeStep() );
            M_bdf_primal_save->setTimeFinal( M_model->timeFinal() );
            M_bdf_primal_save->setOrder( M_model->timeOrder() );

        }
        M_bdf_primal->start();
        M_bdf_primal_save->start();

        //initialization of unknown
        //double initialization_field=M_model->initializationField();
        M_model->initializationField( u );
        //u->setConstant(initialization_field);
        M_bdf_primal->initialize(*u);
        M_bdf_primal_save->initialize(*u);

        //direct problem
        double bdf_coeff = M_bdf_primal->polyDerivCoefficient(0);
        auto vec_bdf_poly = backend_primal_problem->newVector(M_model->functionSpace() );


        for ( M_bdf_primal->start() , M_bdf_primal_save->start() ;
              !M_bdf_primal->isFinished() && !M_bdf_primal_save->isFinished() ;
              M_bdf_primal->next() , M_bdf_primal_save->next() )
        {
            auto bdf_poly = M_bdf_primal->polyDeriv();
            boost::tie( M, A, F ) = M_model->update( mu , M_bdf_primal->time() );

            A->addMatrix( bdf_coeff, M);
            *Rhs = *F[0];
            *vec_bdf_poly = bdf_poly;
            Rhs->addVector( *vec_bdf_poly, *M);

            if( reuse_prec )
            {
                auto ret = backend_primal_problem->solve( _matrix=A, _solution=u, _rhs=Rhs, _reuse_prec=(M_bdf_primal->iteration() >=2) );
                if ( !ret.get<0>() )
                Log()<<"[CRB] WARNING : at time "<<M_bdf_primal->time()<<" we have not converged ( nb_it : "<<ret.get<1>()<<" and residual : "<<ret.get<2>() <<" ) \n";
            }
            else
            {
                auto ret = backend_primal_problem->solve( _matrix=A, _solution=u, _rhs=Rhs );
                if ( !ret.get<0>() )
                Log()<<"[CRB] WARNING : at time "<<M_bdf_primal->time()<<" we have not converged ( nb_it : "<<ret.get<1>()<<" and residual : "<<ret.get<2>() <<" ) \n";
            }

            M_bdf_primal->shiftRight( *u );

            if ( !M_model->isSteady() )
            {
                element_ptrtype projection( new element_type( M_model->functionSpace() ) );
                projectionOnPodSpace( u , projection, "primal" );
                *uproj = *u;
                uproj->add( -1 , *projection );
                M_bdf_primal_save->shiftRight( *uproj );
            }
        }

        std::cout<<"direct problem solved"<<std::endl;

        //dual problem

        M_bdf_dual = bdf( _space=M_model->functionSpace(), _vm=this->vm() , _name="bdf_dual" );

        //set parameters for time discretization
        if( M_model->isSteady() )
        {
            M_bdf_dual->setSteady();
        }
        else
        {
            M_bdf_dual->setTimeInitial( M_model->timeFinal() );
            M_bdf_dual->setTimeStep( -M_model->timeStep() );
            M_bdf_dual->setTimeFinal( M_model->timeInitial() );
            M_bdf_dual->setOrder( M_model->timeOrder() );
        }

        M_bdf_dual->start();

        Adu = M_model->newMatrix();

        //ini
        *udu = *u;
        M_bdf_dual->initialize(*udu);

#if 0

        bdf_coeff = M_bdf_dual->polyDerivCoefficient(0);
        for ( M_bdf_dual->start(); !M_bdf_dual->isFinished(); M_bdf_dual->next() )
        {

            auto bdf_poly = M_bdf_dual->polyDeriv();
            boost::tie(M, A, F ) = M_model->update( mu , M_bdf_dual->time() );
            A->addMatrix( bdf_coeff, M);
            A->transpose( Adu );
            *Rhs = *F[M_output_index];
            Rhs->scale(-1);
            *vec_bdf_poly = bdf_poly;
            Rhs->addVector( *vec_bdf_poly, *M);

            if( reuse_prec )
            {
                auto ret = backend_primal_problem->solve( _matrix=Adu, _solution=udu, _rhs=Rhs, _reuse_prec=(M_bdf_dual->iteration() >=2) );
                if ( !ret.get<0>() )
                Log()<<"[CRB] WARNING (adjoint model) : at time "<<M_bdf_dual->time()<<" we have not converged ( nb_it : "<<ret.get<1>()<<" and residual : "<<ret.get<2>() <<" ) \n";

            }
            else
            {
                auto ret = backend_primal_problem->solve( _matrix=Adu, _solution=udu, _rhs=Rhs );
                if ( !ret.get<0>() )
                Log()<<"[CRB] WARNING (adjoint model) : at time "<<M_bdf_dual->time()<<" we have not converged ( nb_it : "<<ret.get<1>()<<" and residual : "<<ret.get<2>() <<" ) \n";
            }

            if ( !M_model->isSteady() )
            {
                element_ptrtype projection( new element_type( M_model->functionSpace() ) );
                projectionOnPodSpace( udu , projection, "dual" );
                *uproj = *udu;
                uproj->add( -1 , *projection );
            }

            M_bdf_dual->shiftRight( *udu );
        }
        std::cout<<"dual problem solved"<<std::endl;

#endif

        for( int l = 0; l < M_model->Nl(); ++l )
            Log() << "u^T F[" << l << "]= " << inner_product( *u, *F[l] ) << "\n";
        Log() << "[CRB::offlineNoErrorEstimation] energy = " << A->energy( *u, *u ) << "\n";


        M_WNmu->push_back( mu, index );
        M_WNmu_complement = M_WNmu->complement();

        if( M_model->isSteady() )
        {
            M_WN.push_back( *u );
            M_WNdu.push_back( *udu );
        }//end of steady case
        else
        {
            //POD in time
            Log()<<"[CRB::offlineNoErrorEstimation] start of POD \n";

            pod_ptrtype POD = pod_ptrtype( new pod_type(  ) );

            POD->setNm(M_Nm);
            POD->setBdf( M_bdf_primal_save );
            POD->setModel(M_model);
            mode_set_type ModeSet;
            POD->pod(ModeSet);

            //now : loop over number modes per mu
            for(int i=0;i<M_Nm;i++)
            {
                M_WN.push_back( ModeSet[i] );
            }

            //and now the dual
            //POD->setBdf( M_bdf_dual );
            //mode_set_type ModeSetdu;
            //POD->pod(ModeSetdu);
            element_ptrtype mode_zero ( new element_type( M_model->functionSpace() ) );

            for(int i=0;i<M_Nm;i++)
            {
                //M_WNdu.push_back( ModeSetdu[i] ) ;
                M_WNdu.push_back( *mode_zero ) ;
            }
        }//end of transient case

        M_N+=M_Nm;

        if(orthonormalize_primal)
        {
            orthonormalize( M_N, M_WN, M_Nm );
            orthonormalize( M_N, M_WN, M_Nm );
        }
        if( orthonormalize_dual )
        {
            orthonormalize( M_N, M_WNdu, M_Nm );
            orthonormalize( M_N, M_WNdu, M_Nm );
        }


        Log() << "[CRB::offlineNoErrorEstimation] compute Aq_pr, Aq_du, Aq_pr_du" << "\n";
        for( int q = 0; q < M_model->Qa(); ++q )
        {
            M_Aq_pr[q].conservativeResize( M_N, M_N );
            M_Aq_du[q].conservativeResize( M_N, M_N );
            M_Aq_pr_du[q].conservativeResize( M_N, M_N );

            // only compute the last line and last column of reduced matrices
            for(int i = M_N-M_Nm; i < M_N; i++ )
            {
                for( int j = 0; j < M_N; ++j )
                {
                    M_Aq_pr[q]( i, j ) = Aq[q]->energy( M_WN[i], M_WN[j] );
                    M_Aq_du[q]( i, j ) = Aq[q]->energy( M_WNdu[i], M_WNdu[j], true );
                    M_Aq_pr_du[q]( i, j ) = Aq[q]->energy( M_WNdu[i], M_WN[j] );
                }
            }
            for(int j=M_N-M_Nm; j < M_N; j++ )
            {
                for( int i = 0; i < M_N; ++i )
                {
                    M_Aq_pr[q]( i, j ) = Aq[q]->energy( M_WN[i], M_WN[j] );
                    M_Aq_du[q]( i, j ) = Aq[q]->energy( M_WNdu[i], M_WNdu[j], true );
                    M_Aq_pr_du[q]( i, j ) = Aq[q]->energy( M_WNdu[i], M_WN[j] );
                }
            }
        }
        Log() << "[CRB::offlineNoErrorEstimation] compute Mq_pr, Mq_du, Mq_pr_du" << "\n";
        for( int q = 0; q < M_model->Qm(); ++q )
        {
            M_Mq_pr[q].conservativeResize( M_N, M_N );
            M_Mq_du[q].conservativeResize( M_N, M_N );
            M_Mq_pr_du[q].conservativeResize( M_N, M_N );

            // only compute the last line and last column of reduced matrices
            for( int i=M_N-M_Nm ; i < M_N; i++ )
            {
                for( int j = 0; j < M_N; ++j )
                {
                    M_Mq_pr[q]( i, j ) = Mq[q]->energy( M_WN[i], M_WN[j] );
                    M_Mq_du[q]( i, j ) = Mq[q]->energy( M_WNdu[i], M_WNdu[j], true );
                    M_Mq_pr_du[q]( i, j ) = Mq[q]->energy( M_WNdu[i], M_WN[j] );
                }
            }
            for(int j = M_N-M_Nm; j < M_N ; j++ )
            {
                for( int i = 0; i < M_N; ++i )
                {
                    M_Mq_pr[q]( i, j ) = Mq[q]->energy( M_WN[i], M_WN[j] );
                    M_Mq_du[q]( i, j ) = Mq[q]->energy( M_WNdu[i], M_WNdu[j], true );
                    M_Mq_pr_du[q]( i, j ) = Mq[q]->energy( M_WNdu[i], M_WN[j] );
                }
            }
        }

        Log() << "[CRB::offlineNoErrorEstimation] compute Fq_pr, Fq_du" << "\n";
        for( int q = 0; q < M_model->Ql(0); ++q )
        {
            M_Fq_pr[q].conservativeResize( M_N );
            M_Fq_du[q].conservativeResize( M_N );
            for( int l = 1; l <= M_Nm; ++l )
            {
                int index = M_N-l;
                M_Fq_pr[q]( index ) = M_model->Fq( 0, q, M_WN[index] );
                M_Fq_du[q]( index ) = M_model->Fq( 0, q, M_WNdu[index] );
            }
        }
        Log() << "[CRB::offlineNoErrorEstimation] compute Lq_pr, Lq_du" << "\n";
        for( int q = 0; q < M_model->Ql(M_output_index); ++q )
        {
            M_Lq_pr[q].conservativeResize( M_N );
            M_Lq_du[q].conservativeResize( M_N );
            for( int l = 1; l <= M_Nm; ++l )
            {
                int index = M_N-l;
                M_Lq_pr[q]( index ) = M_model->Fq( M_output_index, q, M_WN[index] );
                M_Lq_du[q]( index ) = M_model->Fq( M_output_index, q, M_WNdu[index] );
            }

        }

        check( M_WNmu->size() );
        double error=M_iter_max-M_N;
        M_rbconv.insert( convergence( M_N, error ) );

        std::cout << "============================================================\n";
        Log() <<"========================================"<<"\n";
        std::cout<<"number of elements in the reduced basis : "<<M_N<<std::endl;

    }
    this->saveDB();
    //return M_rbconv;
}


template<typename TruthModelType>
void
CRB<TruthModelType>::offlineNoErrorEstimation(mpl::bool_<false>)
{

    //in steady case M_Nm must be 1
    M_Nm = 1;

    parameter_type mu( M_Dmu );

    double relative_error = 1e30;

    // empty sets
    M_WNmu->clear();

    // start with M_C = { arg min mu, mu \in Xi }
    size_type index;
    boost::tie( mu, index ) = M_Xi->max();

    sparse_matrix_ptrtype A,At;
    std::vector<vector_ptrtype> F,L;


    //useful for POD in time
    int Nm = this->vm()["crb.Nm"].template as<int>() ;
    const int Ndof = M_model->functionSpace()->nDof();//number of dofs used
    matrixN_type M;
    matrixN_type Mdu;


    // dimension of reduced basis space
    M_N = 1;

    size_type Np = 1;

    Log() << "[CRB::offlineNoErrorEstimation] compute affine decomposition\n";
    std::vector<sparse_matrix_ptrtype> Aq;
    std::vector<std::vector<vector_ptrtype> > Fq,Lq;
    boost::tie( Aq, Fq ) = M_model->computeAffineDecomposition();

    Log() << "[CRB::offlineNoErrorEstimation] allocate reduced basis data structures\n";
    M_Aq_pr.resize( M_model->Qa() );
    M_Aq_du.resize( M_model->Qa() );
    M_Aq_pr_du.resize( M_model->Qa() );
    M_Fq_pr.resize( M_model->Ql(0) );
    M_Fq_du.resize( M_model->Ql(0) );
    M_Lq_pr.resize( M_model->Ql(M_output_index) );
    M_Lq_du.resize( M_model->Ql(M_output_index) );

    element_ptrtype u( new element_type( M_model->functionSpace() ) );
    element_ptrtype udu( new element_type( M_model->functionSpace() ) );
    vector_ptrtype U( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype Rhs( M_backend->newVector( M_model->functionSpace() ) );

    sampling_ptrtype Sampling;
    Sampling = sampling_ptrtype( new sampling_type( M_Dmu ) );

    Sampling->logEquidistribute( M_iter_max );
    Log() << "[CRB::offlineNoErrorEstimation] sampling parameter space with "<< M_iter_max <<"elements done\n";


    Log() << "[CRB::offlineNoErrorEstimation] strategy "<< M_error_type <<"\n";
    std::cout << "[CRB::offlineNoErrorEstimation] strategy "<< M_error_type <<"\n";

    for( size_type index = 0; index < Sampling->size(); ++index )
    {
        parameter_type const& mu = M_Xi->at( index );

        Log() <<"========================================"<<"\n";
        std::cout << "============================================================\n";
        std::cout << "N=" << M_N << "\n";
        Log() << "N=" << M_N << "\n";


        backend_ptrtype backendA = backend_type::build( BACKEND_PETSC );
        backend_ptrtype backendAt = backend_type::build( BACKEND_PETSC );

        // for a given parameter \p mu assemble the left and right hand side
        Log() << "[CRB::offlineNoErrorEstimation] update model for parameter" << "\n";
        boost::tie( A, F ) = M_model->update( mu );

        At = M_model->newMatrix();
        Log() << "[CRB::offlineNoErrorEstimation] transpose primal matrix" << "\n";
        A->transpose( At );
        u->setName( (boost::format( "fem-primal-%1%" ) % (M_N-1)).str() );
        udu->setName( (boost::format( "fem-dual-%1%" ) % (M_N-1)).str() );

        Log() << "[CRB::offlineNoErrorEstimation] solving primal" << "\n";
        backendA->solve( _matrix=A,  _solution=U, _rhs=F[0], _prec=A );
        //std::cout << "solving primal done" << std::endl;
        *u = *U;
        *Rhs = *F[M_output_index];
        Rhs->scale( -1 );
        Log() << "[CRB::offlineNoErrorEstimation] solving dual" << "\n";
        backendAt->solve( _matrix=At,  _solution=U, _rhs=Rhs, _prec=At );
        *udu = *U;
        //std::cout << "solving dual done" << std::endl;

        for( int l = 0; l < M_model->Nl(); ++l )
            Log() << "u^T F[" << l << "]= " << inner_product( *u, *F[l] ) << "\n";

        Log() << "[CRB::offlineNoErrorEstimation] energy = " << A->energy( *u, *u ) << "\n";

        M_WNmu->push_back( mu, index );
        M_WNmu_complement = M_WNmu->complement();

        Log() << "[CRB::offlineNoErrorEstimation] orthonormalize basis functions" << "\n";

        M_WN.push_back( *u );
        M_WNdu.push_back( *udu );

        if(orthonormalize_primal)
        {
            orthonormalize( M_N, M_WN );
            orthonormalize( M_N, M_WN );
        }
        if(orthonormalize_dual)
        {
            orthonormalize( M_N, M_WNdu );
            orthonormalize( M_N, M_WNdu );
        }

        Log() << "[CRB::offlineNoErrorEstimation] compute Aq_pr, Aq_du, Aq_pr_du" << "\n";
        for( int q = 0; q < M_model->Qa(); ++q )
            {
                M_Aq_pr[q].conservativeResize( M_N, M_N );
                M_Aq_du[q].conservativeResize( M_N, M_N );
                M_Aq_pr_du[q].conservativeResize( M_N, M_N );

                // only compute the last line and last column of reduced matrices
                int i = M_N-1;
                for( int j = 0; j < M_N; ++j )
                {
                    M_Aq_pr[q]( i, j ) = Aq[q]->energy( M_WN[i], M_WN[j] );
                    M_Aq_du[q]( i, j ) = Aq[q]->energy( M_WNdu[i], M_WNdu[j], true );
                    M_Aq_pr_du[q]( i, j ) = Aq[q]->energy( M_WNdu[i], M_WN[j] );
                }
                int j = M_N-1;
                for( int i = 0; i < M_N; ++i )
                {
                    M_Aq_pr[q]( i, j ) = Aq[q]->energy( M_WN[i], M_WN[j] );
                    M_Aq_du[q]( i, j ) = Aq[q]->energy( M_WNdu[i], M_WNdu[j], true );
                    M_Aq_pr_du[q]( i, j ) = Aq[q]->energy( M_WNdu[i], M_WN[j] );
                }
            }
            Log() << "[CRB::offlineNoErrorEstimation] compute Fq_pr, Fq_du" << "\n";
            for( int q = 0; q < M_model->Ql(0); ++q )
                {
                    M_Fq_pr[q].conservativeResize( M_N );
                    M_Fq_du[q].conservativeResize( M_N );
                    M_Fq_pr[q]( M_N-1 ) = M_model->Fq( 0, q, M_WN[M_N-1] );
                    M_Fq_du[q]( M_N-1 ) = M_model->Fq( 0, q, M_WNdu[M_N-1] );
                }
            Log() << "[CRB::offlineNoErrorEstimation] compute Lq_pr, Lq_du" << "\n";
            for( int q = 0; q < M_model->Ql(M_output_index); ++q )
                {
                    M_Lq_pr[q].conservativeResize( M_N );
                    M_Lq_du[q].conservativeResize( M_N );
                    M_Lq_pr[q]( M_N-1 ) = M_model->Fq( M_output_index, q, M_WN[M_N-1] );
                    M_Lq_du[q]( M_N-1 ) = M_model->Fq( M_output_index, q, M_WNdu[M_N-1] );
                }

        }

        if ( 0 )
        {
            u->zero();
            u->add( 1, M_WN[M_N-1] );
            udu->zero();
            udu->add( 1, M_WNdu[M_N-1] );
            vector_ptrtype Aun( M_backend->newVector( M_model->functionSpace() ) );
            vector_ptrtype Atun( M_backend->newVector( M_model->functionSpace() ) );
            vector_ptrtype Un( M_backend->newVector( M_model->functionSpace() ) );
            vector_ptrtype Undu( M_backend->newVector( M_model->functionSpace() ) );
            vector_ptrtype Frhs( M_backend->newVector( M_model->functionSpace() ) );
            vector_ptrtype Lrhs( M_backend->newVector( M_model->functionSpace() ) );
            *Un = *u;
            *Undu = *udu;
            A->multVector( Un, Aun );
            At->multVector( Undu, Atun );
            Aun->scale( -1 );
            Atun->scale( -1 );
            *Frhs = *F[0];
            *Lrhs = *F[M_output_index];
            Log() << "[CRB::offlineNoErrorEstimation] residual (f,f) " << M_N-1 << ":=" << M_model->scalarProduct( Frhs, Frhs ) << "\n";
            Log() << "[CRB::offlineNoErrorEstimation] residual (f,A) " << M_N-1 << ":=" << 2*M_model->scalarProduct( Frhs, Aun ) << "\n";
            Log() << "[CRB::offlineNoErrorEstimation] residual (A,A) " << M_N-1 << ":=" << M_model->scalarProduct( Aun, Aun ) << "\n";

            Log() << "[CRB::offlineNoErrorEstimation] residual (l,l) " << M_N-1 << ":=" << M_model->scalarProduct( Lrhs, Lrhs ) << "\n";
            Log() << "[CRB::offlineNoErrorEstimation] residual (l,At) " << M_N-1 << ":=" << 2*M_model->scalarProduct( Lrhs, Atun ) << "\n";
            Log() << "[CRB::offlineNoErrorEstimation] residual (At,At) " << M_N-1 << ":=" << M_model->scalarProduct( Atun, Atun ) << "\n";

            Aun->add( *Frhs );
            Lrhs->scale( -1 );
            Atun->add( *Lrhs );
            double err_primal = math::sqrt ( M_model->scalarProduct( Aun, Aun ) );
            double err_dual = math::sqrt ( M_model->scalarProduct( Atun, Atun ) );
            Log() << "[CRB::offlineNoErrorEstimation] primal residual for reduced basis function " << M_N-1 << ":=" << err_primal << "\n";
            Log() << "[CRB::offlineNoErrorEstimation] dual residual for reduced basis function " << M_N-1 << ":=" << err_dual << "\n";
        }

        check( M_WNmu->size() );
        double error=M_iter_max-M_N;
        M_rbconv.insert( convergence( M_N, error ) );

        ++M_N;

        std::cout << "============================================================\n";
        Log() <<"========================================"<<"\n";
        std::cout<<"number of elements in the reduced basis : "<<M_N<<std::endl;

    --M_N;

    this->saveDB();
    //return M_rbconv;
}





template<typename TruthModelType>
void
CRB<TruthModelType>::offlineWithErrorEstimation()
{
    return offlineWithErrorEstimation( mpl::bool_<model_type::is_time_dependent>() );
}

template<typename TruthModelType>
void
CRB<TruthModelType>::offlineWithErrorEstimation(mpl::bool_<true>)
{
    parameter_type mu( M_Dmu );

    double relative_error = 1e30;

    // empty sets
    M_WNmu->clear();

    // start with M_C = { arg min mu, mu \in Xi }
    size_type index;
    boost::tie( mu, index ) = M_Xi->max();

    std::cout << "  -- start with mu = " << mu << "\n";
    //std::cout << " -- WN size :  " << M_WNmu->size() << "\n";

    sparse_matrix_ptrtype M,A,Adu,At;
    std::vector<vector_ptrtype> F,L;


    //useful for POD in time
    const int Ndof = M_model->functionSpace()->nDof();//number of dofs used

    // dimension of reduced basis space
    M_N = 0;

    size_type Np = 1;

    Log() << "[CRB::offlineWithErrorEstimation] compute affine decomposition\n";
    std::vector<sparse_matrix_ptrtype> Aq;
    std::vector<sparse_matrix_ptrtype> Mq;
    std::vector<std::vector<vector_ptrtype> > Fq,Lq;
    boost::tie(Mq, Aq, Fq) = M_model->computeAffineDecomposition();

    // scm offline stage: build C_K
    if ( M_error_type == CRB_RESIDUAL_SCM )
        std::vector<boost::tuple<double,double,double> > M_rbconv2 = M_scm->offline();

    double maxerror = 1e10;
    //boost::tie( maxerror, mu, index ) = maxErrorBounds( N );

    Log() << "[CRB::offlineWithErrorEstimation] allocate reduced basis data structures\n";
    M_Aq_pr.resize( M_model->Qa() );
    M_Aq_du.resize( M_model->Qa() );
    M_Aq_pr_du.resize( M_model->Qa() );

    M_Mq_pr.resize( M_model->Qm() );
    M_Mq_du.resize( M_model->Qm() );
    M_Mq_pr_du.resize( M_model->Qm() );

    M_Fq_pr.resize( M_model->Ql(0) );
    M_Fq_du.resize( M_model->Ql(0) );
    M_Lq_pr.resize( M_model->Ql(M_output_index) );
    M_Lq_du.resize( M_model->Ql(M_output_index) );

    element_ptrtype u( new element_type( M_model->functionSpace() ) );
    element_ptrtype uproj( new element_type( M_model->functionSpace() ) );
    element_ptrtype udu( new element_type( M_model->functionSpace() ) );
    vector_ptrtype Rhs( M_backend->newVector( M_model->functionSpace() ) );


    Log() << "[CRB::offlineWithErrorEstimation] starting offline adaptive loop\n";
    Log() << "[CRB::offlineNoErrorEstimation] strategy "<< M_error_type <<"\n";
    std::cout << "[CRB::offlineNoErrorEstimation] strategy "<< M_error_type <<"\n";

    bool reuse_prec = this->vm()["crb.reuse_prec"].template as<bool>() ;

    while ( maxerror > M_tolerance && M_N < M_iter_max )
    {
        boost::timer timer, timer2;
        Log() <<"========================================"<<"\n";
        std::cout << "============================================================\n";
        std::cout << "N=" << M_N << "/"  << M_iter_max << " maxerror=" << maxerror << " / "  << M_tolerance << "\n";
        Log() << "N=" << M_N << "/"  << M_iter_max << " maxerror=" << maxerror << " / "  << M_tolerance << "\n";


        backend_ptrtype backend_primal_problem = backend_type::build( BACKEND_PETSC );
        backend_ptrtype backend_dual_problem = backend_type::build( BACKEND_PETSC );

        // for a given parameter \p mu assemble the left and right hand side

        u->setName( (boost::format( "fem-primal-%1%" ) % (M_N-1)).str() );

        if(M_model->isSteady() )
        {
            boost::tie( M, A, F ) = M_model->update( mu , 1e30 );
            std::cout << "  -- updated model for parameter in " << timer2.elapsed() << "s\n"; timer2.restart();
            Log() << "[CRB::offlineWithErrorEstimation] transpose primal matrix" << "\n";
            At = M_model->newMatrix();
            A->transpose( At );
            u->setName( (boost::format( "fem-primal-%1%" ) % (M_N-1)).str() );
            udu->setName( (boost::format( "fem-dual-%1%" ) % (M_N-1)).str() );

            Log() << "[CRB::offlineWithErrorEstimation] solving primal" << "\n";
            backend_primal_problem->solve( _matrix=A,  _solution=u, _rhs=F[0], _prec=A );
            //std::cout << "solving primal done" << std::endl;
            std::cout << "  -- primal problem solved in " << timer2.elapsed() << "s\n"; timer2.restart();
            *Rhs = *F[M_output_index];
            Rhs->scale( -1 );
            backend_dual_problem->solve( _matrix=At,  _solution=udu, _rhs=Rhs, _prec=At );
        }
        else
        {

            M_bdf_primal = bdf( _space=M_model->functionSpace(), _vm=this->vm() , _name="bdf_primal");
            M_bdf_primal_save = bdf( _space=M_model->functionSpace(), _vm=this->vm() , _name="bdf_primal_save");

            //set parameters for time discretization
            if( M_model->isSteady() )
            {
                M_bdf_primal->setSteady();
                M_bdf_primal_save->setSteady();
            }
            else
            {
               M_bdf_primal->setTimeInitial( M_model->timeInitial() );
               M_bdf_primal->setTimeStep( M_model->timeStep() );
               M_bdf_primal->setTimeFinal( M_model->timeFinal() );
               M_bdf_primal->setOrder( M_model->timeOrder() );

               M_bdf_primal_save->setTimeInitial( M_model->timeInitial() );
               M_bdf_primal_save->setTimeStep( M_model->timeStep() );
               M_bdf_primal_save->setTimeFinal( M_model->timeFinal() );
               M_bdf_primal_save->setOrder( M_model->timeOrder() );
            }

            M_bdf_primal->start();
            M_bdf_primal_save->start();

            //initialization of unknown
            M_model->initializationField( u );
            M_bdf_primal->initialize(*u);
            M_bdf_primal_save->initialize(*u);

            //direct problem
            double bdf_coeff = M_bdf_primal->polyDerivCoefficient(0);

            auto vec_bdf_poly = backend_primal_problem->newVector(M_model->functionSpace() );
            for ( M_bdf_primal->start(),M_bdf_primal_save->start();
                  !M_bdf_primal->isFinished() , !M_bdf_primal_save->isFinished();
                  M_bdf_primal->next() , M_bdf_primal_save->next() )
            {
                auto bdf_poly = M_bdf_primal->polyDeriv();

                boost::tie(M, A, F ) = M_model->update( mu , M_bdf_primal->time() );

                A->addMatrix( bdf_coeff, M);
                *Rhs = *F[0];
                *vec_bdf_poly = bdf_poly;
                Rhs->addVector( *vec_bdf_poly, *M);

                A->close();

                if( reuse_prec )
                {
                    auto ret = backend_primal_problem->solve( _matrix=A, _solution=u, _rhs=Rhs, _reuse_prec=(M_bdf_primal->iteration() >=2) );
                    if ( !ret.get<0>() )
                        Log()<<"[CRB] WARNING : at time "<<M_bdf_primal->time()<<" we have not converged ( nb_it : "<<ret.get<1>()<<" and residual : "<<ret.get<2>() <<" ) \n";
                }
                else
                {
                    auto ret = backend_primal_problem->solve( _matrix=A, _solution=u, _rhs=Rhs );
                    if ( !ret.get<0>() )
                        Log()<<"[CRB] WARNING : at time "<<M_bdf_primal->time()<<" we have not converged ( nb_it : "<<ret.get<1>()<<" and residual : "<<ret.get<2>() <<" ) \n";
                }

                M_bdf_primal->shiftRight( *u );

                if( ! M_model->isSteady() )
                {
                    element_ptrtype projection ( new element_type (M_model->functionSpace() ) );
                    projectionOnPodSpace ( u , projection, "primal" );
                    *uproj=*u;
                    M_bdf_primal_save->shiftRight( *uproj );
                }


                Log() << "u^T Rhs "<< inner_product( *u, *Rhs ) <<  " at time : "<<M_bdf_primal->time()<<"\n";
                for( int l = 0; l < M_model->Nl(); ++l )
                    Log() << "u^T F[" << l << "]= " << inner_product( *u, *F[l] ) << " at time : "<<M_bdf_primal->time()<<"\n";

                Log() << "[CRB::offlineWithErrorEstimation] energy = " << A->energy( *u, *u ) << "\n";
            }

            std::cout<<"direct problem solved"<<std::endl;

            //dual problem

            M_bdf_dual = bdf( _space=M_model->functionSpace(), _vm=this->vm() , _name="bdf_dual" );

#if 0
            //set parameters for time discretization
            if( M_model->isSteady() )
            {
                M_bdf_primal->setSteady();
            }
            else
            {
                M_bdf_dual->setTimeInitial( M_model->timeFinal() );
                M_bdf_dual->setTimeStep( -M_model->timeStep() );
                M_bdf_dual->setTimeFinal( M_model->timeInitial() );
                M_bdf_dual->setOrder( M_model->timeOrder() );
            }

            Adu = M_model->newMatrix();
            M_bdf_dual->start();

            //initialization
            *udu = *u;
            M_bdf_dual->initialize(*udu);

            bdf_coeff = M_bdf_dual->polyDerivCoefficient(0);
            for ( M_bdf_dual->start(); !M_bdf_dual->isFinished(); M_bdf_dual->next() )
            {
                auto bdf_poly = M_bdf_dual->polyDeriv();
                boost::tie(M, A, F ) = M_model->update( mu , M_bdf_dual->time() );
                A->addMatrix( bdf_coeff, M);
                A->transpose( Adu );
                *Rhs = *F[M_output_index];
                Rhs->scale(-1);
                *vec_bdf_poly = bdf_poly;
                Rhs->addVector( *vec_bdf_poly, *M);

                if( reuse_prec )
                {
                    auto ret = backend_primal_problem->solve( _matrix=Adu, _solution=udu, _rhs=Rhs, _reuse_prec=(M_bdf_dual->iteration() >=2) );
                    if ( !ret.get<0>() )
                        Log()<<"[CRB] WARNING (adjoint model) : at time "<<M_bdf_dual->time()<<" we have not converged ( nb_it : "<<ret.get<1>()<<" and residual : "<<ret.get<2>() <<" ) \n";
                }
                else
                {
                    auto ret = backend_primal_problem->solve( _matrix=Adu, _solution=udu, _rhs=Rhs );
                    if ( !ret.get<0>() )
                        Log()<<"[CRB] WARNING (adjoint model) : at time "<<M_bdf_dual->time()<<" we have not converged ( nb_it : "<<ret.get<1>()<<" and residual : "<<ret.get<2>() <<" ) \n";
                }

                if( ! M_model->isSteady() )
                {
                    element_ptrtype projection ( new element_type (M_model->functionSpace() ) );
                    projectionOnPodSpace ( udu , projection, "dual" );
                    *uproj=*udu;
                    M_bdf_primal_save->shiftRight( *uproj );
                }


                M_bdf_dual->shiftRight( *udu );
            }
            std::cout<<"dual problem solved"<<std::endl;

#endif

        }//end of transient case

        for( int l = 0; l < M_model->Nl(); ++l )
            Log() << "u^T F[" << l << "]= " << inner_product( *u, *F[l] ) << "\n";


        Log() << "[CRB::offlineNWithErrorEstimation] energy = " << A->energy( *u, *u ) << "\n";

        M_WNmu->push_back( mu, index );
        M_WNmu_complement = M_WNmu->complement();

        if( M_model->isSteady() )
        {
            M_WN.push_back( *u );
            M_WNdu.push_back( *udu );
        }//end of steady case
        else
        {
            //POD in time
            Log()<<"[CRB::offlineNoErrorEstimation] start of POD \n";

            pod_ptrtype POD = pod_ptrtype( new pod_type(  ) );

            POD->setNm(M_Nm);
            POD->setBdf( M_bdf_primal_save );
            POD->setModel(M_model);
            mode_set_type ModeSet;
            POD->pod(ModeSet);
            //now : loop over number modes per mu
            for(int i=0;i<M_Nm;i++)
            {
                M_WN.push_back( ModeSet[i] );
            }


            //and now the dual
            //POD->setBdf( M_bdf_dual );
            //mode_set_type ModeSetdu;
            //POD->pod(ModeSetdu);
            element_ptrtype mode_zero ( new element_type( M_model->functionSpace() ) );

            for(int i=0;i<M_Nm;i++)
            {
                //M_WNdu.push_back( ModeSetdu[i] ) ;
                M_WNdu.push_back( *mode_zero ) ;
            }

        }//end of transient case

        M_N+=M_Nm;


        if(orthonormalize_primal)
        {
            orthonormalize( M_N, M_WN, M_Nm );
            orthonormalize( M_N, M_WN, M_Nm );
        }
        if(orthonormalize_dual)
        {
            orthonormalize( M_N, M_WNdu, M_Nm );
            orthonormalize( M_N, M_WNdu, M_Nm );
        }

        Log() << "[CRB::offlineWithErrorEstimation] compute Aq_pr, Aq_du, Aq_pr_du" << "\n";
        for( int q = 0; q < M_model->Qa(); ++q )
        {
            M_Aq_pr[q].conservativeResize( M_N, M_N );
            M_Aq_du[q].conservativeResize( M_N, M_N );
            M_Aq_pr_du[q].conservativeResize( M_N, M_N );

            // only compute the last line and last column of reduced matrices
            for(int i = M_N-M_Nm; i < M_N; i++ )
            {
                for( int j = 0; j < M_N; ++j )
                {
                    M_Aq_pr[q]( i, j ) = Aq[q]->energy( M_WN[i], M_WN[j] );
                    M_Aq_du[q]( i, j ) = Aq[q]->energy( M_WNdu[i], M_WNdu[j], true );
                    M_Aq_pr_du[q]( i, j ) = Aq[q]->energy( M_WNdu[i], M_WN[j] );
                }
            }
            for(int j=M_N-M_Nm; j < M_N; j++ )
            {
                for( int i = 0; i < M_N; ++i )
                {
                    M_Aq_pr[q]( i, j ) = Aq[q]->energy( M_WN[i], M_WN[j] );
                    M_Aq_du[q]( i, j ) = Aq[q]->energy( M_WNdu[i], M_WNdu[j], true );
                    M_Aq_pr_du[q]( i, j ) = Aq[q]->energy( M_WNdu[i], M_WN[j] );
                }
            }
        }
        Log() << "[CRB::offlineWithErrorEstimation] compute Mq_pr, Mq_du, Mq_pr_du" << "\n";
        for( int q = 0; q < M_model->Qm(); ++q )
        {
            M_Mq_pr[q].conservativeResize( M_N, M_N );
            M_Mq_du[q].conservativeResize( M_N, M_N );
            M_Mq_pr_du[q].conservativeResize( M_N, M_N );

            // only compute the last line and last column of reduced matrices
            for( int i=M_N-M_Nm ; i < M_N; i++ )
            {
                for( int j = 0; j < M_N; ++j )
                {
                    M_Mq_pr[q]( i, j ) = Mq[q]->energy( M_WN[i], M_WN[j] );
                    M_Mq_du[q]( i, j ) = Mq[q]->energy( M_WNdu[i], M_WNdu[j], true );
                    M_Mq_pr_du[q]( i, j ) = Mq[q]->energy( M_WNdu[i], M_WN[j] );
                }
            }
            for(int j = M_N-M_Nm; j < M_N ; j++ )
            {
                for( int i = 0; i < M_N; ++i )
                {
                    M_Mq_pr[q]( i, j ) = Mq[q]->energy( M_WN[i], M_WN[j] );
                    M_Mq_du[q]( i, j ) = Mq[q]->energy( M_WNdu[i], M_WNdu[j], true );
                    M_Mq_pr_du[q]( i, j ) = Mq[q]->energy( M_WNdu[i], M_WN[j] );
                }
            }
        }

        Log() << "[CRB::offlineWithErrorEstimation] compute Fq_pr, Fq_du" << "\n";
        for( int q = 0; q < M_model->Ql(0); ++q )
        {
            M_Fq_pr[q].conservativeResize( M_N );
            M_Fq_du[q].conservativeResize( M_N );
            for( int l = 1; l <= M_Nm; ++l )
            {
                int index = M_N-l;
                M_Fq_pr[q]( index ) = M_model->Fq( 0, q, M_WN[index] );
                M_Fq_du[q]( index ) = M_model->Fq( 0, q, M_WNdu[index] );
            }
        }
        Log() << "[CRB::offlineWithErrorEstimation] compute Lq_pr, Lq_du" << "\n";
        for( int q = 0; q < M_model->Ql(M_output_index); ++q )
        {
            M_Lq_pr[q].conservativeResize( M_N );
            M_Lq_du[q].conservativeResize( M_N );
            for( int l = 1; l <= M_Nm; ++l )
            {
                int index = M_N-l;
                M_Lq_pr[q]( index ) = M_model->Fq( M_output_index, q, M_WN[index] );
                M_Lq_du[q]( index ) = M_model->Fq( M_output_index, q, M_WNdu[index] );
            }

        }


        timer2.restart();
        if(M_error_type==CRB_RESIDUAL || M_error_type == CRB_RESIDUAL_SCM)
            {
                std::cout << "  -- N2Q2 update starts\n";
                generateN2Q2( M_N );
                Log()<<"[CRB::offlineWithErrorEstimation] end of call generateN2Q2 and M_N = "<< M_N <<"\n";
                std::cout << "  -- N2Q2 updated in " << timer2.elapsed() << "s\n"; timer2.restart();
            }


        boost::tie( maxerror, mu, index ) = maxErrorBounds( M_N );
        std::cout << "  -- max error bounds computed in " << timer2.elapsed() << "s\n"; timer2.restart();

        M_rbconv.insert( convergence( M_N, maxerror ) );

        //mu = M_Xi->at( M_N );//M_WNmu_complement->min().get<0>();

        check( M_WNmu->size() );


        if ( this->vm()["crb.check.rb"].template as<int>() == 1 )std::cout << "  -- check reduced basis done in " << timer2.elapsed() << "s\n"; timer2.restart();

        std::cout << "time: " << timer.elapsed() << "\n";
        Log() << "time: " << timer.elapsed() << "\n";
        std::cout << "============================================================\n";
        Log() <<"========================================"<<"\n";
    }
    std::cout<<"number of elements in the reduced basis : "<<M_N<<std::endl;

    this->saveDB();
    std::cout << "Offline CRB is done\n";


}

template<typename TruthModelType>
void
CRB<TruthModelType>::offlineWithErrorEstimation(mpl::bool_<false>)
{

    parameter_type mu( M_Dmu );

    double relative_error = 1e30;

    // empty sets
    M_WNmu->clear();


    // start with M_C = { arg min mu, mu \in Xi }
    size_type index;
    boost::tie( mu, index ) = M_Xi->max();

    std::cout << "  -- start with mu = " << mu << "\n";
    //std::cout << " -- WN size :  " << M_WNmu->size() << "\n";

    sparse_matrix_ptrtype A,At;
    std::vector<vector_ptrtype> F,L;

    // dimension of reduced basis space
    M_N = 1;

    size_type Np = 1;

    Log() << "[CRB::offlineWithErrorEstimation] compute affine decomposition\n";
    std::vector<sparse_matrix_ptrtype> Aq;
    std::vector<std::vector<vector_ptrtype> > Fq,Lq;
    boost::tie( Aq, Fq) = M_model->computeAffineDecomposition();


    // scm offline stage: build C_K
    if ( M_error_type == CRB_RESIDUAL_SCM )
        std::vector<boost::tuple<double,double,double> > M_rbconv2 = M_scm->offline();

    double maxerror = 1e10;
    //boost::tie( maxerror, mu, index ) = maxErrorBounds( N );

    Log() << "[CRB::offlineWithErrorEstimation] allocate reduced basis data structures\n";
    M_Aq_pr.resize( M_model->Qa() );
    M_Aq_du.resize( M_model->Qa() );
    M_Aq_pr_du.resize( M_model->Qa() );
    M_Fq_pr.resize( M_model->Ql(0) );
    M_Fq_du.resize( M_model->Ql(0) );
    M_Lq_pr.resize( M_model->Ql(M_output_index) );
    M_Lq_du.resize( M_model->Ql(M_output_index) );

    element_ptrtype u( new element_type( M_model->functionSpace() ) );
    element_ptrtype udu( new element_type( M_model->functionSpace() ) );
    vector_ptrtype U( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype Rhs( M_backend->newVector( M_model->functionSpace() ) );

    Log() << "[CRB::offlineWithErrorEstimation] starting offline adaptive loop\n";
    Log() << "[CRB::offlineNoErrorEstimation] strategy "<< M_error_type <<"\n";
    std::cout << "[CRB::offlineNoErrorEstimation] strategy "<< M_error_type <<"\n";

    while ( maxerror > M_tolerance && M_N <= M_iter_max )
    {
        boost::timer timer, timer2;
        Log() <<"========================================"<<"\n";
        std::cout << "============================================================\n";
        std::cout << "N=" << M_N << "/"  << M_iter_max << " maxerror=" << maxerror << " / "  << M_tolerance << "\n";
        Log() << "N=" << M_N << "/"  << M_iter_max << " maxerror=" << maxerror << " / "  << M_tolerance << "\n";


        backend_ptrtype backendA = backend_type::build( BACKEND_PETSC );
        backend_ptrtype backendAt = backend_type::build( BACKEND_PETSC );

        // for a given parameter \p mu assemble the left and right hand side
        Log() << "[CRB::offlineWithErrorEstimation] update model for parameter" << "\n";

        boost::tie( A, F ) = M_model->update( mu );
        std::cout << "  -- updated model for parameter in " << timer2.elapsed() << "s\n"; timer2.restart();

        Log() << "[CRB::offlineWithErrorEstimation] transpose primal matrix" << "\n";
        At = M_model->newMatrix();
        A->transpose( At );
        u->setName( (boost::format( "fem-primal-%1%" ) % (M_N-1)).str() );
        udu->setName( (boost::format( "fem-dual-%1%" ) % (M_N-1)).str() );

        Log() << "[CRB::offlineWithErrorEstimation] solving primal" << "\n";
        backendA->solve( _matrix=A,  _solution=U, _rhs=F[0], _prec=A );
        //std::cout << "solving primal done" << std::endl;
        *u = *U;
        std::cout << "  -- primal problem solved in " << timer2.elapsed() << "s\n"; timer2.restart();
        *Rhs = *F[M_output_index];
        Rhs->scale( -1 );
        backendAt->solve( _matrix=At,  _solution=U, _rhs=Rhs, _prec=At );
        *udu = *U;
        //std::cout << "solving dual done" << std::endl;
        std::cout << "  -- dual problem solved in " << timer2.elapsed() << "s\n"; timer2.restart();

        for( int l = 0; l < M_model->Nl(); ++l )
            Log() << "u^T F[" << l << "]= " << inner_product( *u, *F[l] ) << "\n";

        Log() << "[CRB::offlineWithErrorEstimation] energy = " << A->energy( *u, *u ) << "\n";

        M_WNmu->push_back( mu, index );
        M_WNmu_complement = M_WNmu->complement();

        M_WN.push_back( *u );
        M_WNdu.push_back( *udu );

        if(orthonormalize_primal)
        {
            orthonormalize( M_N, M_WN );
            orthonormalize( M_N, M_WN );
        }
        if(orthonormalize_dual)
        {
            orthonormalize( M_N, M_WNdu );
            orthonormalize( M_N, M_WNdu );
        }

        Log() << "[CRB::offlineWithErrorEstimation] compute Aq_pr, Aq_du, Aq_pr_du" << "\n"; timer2.restart();
        for( int q = 0; q < M_model->Qa(); ++q )
            {
                M_Aq_pr[q].conservativeResize( M_N, M_N );
                M_Aq_du[q].conservativeResize( M_N, M_N );
                M_Aq_pr_du[q].conservativeResize( M_N, M_N );

                // only compute the last line and last column of reduced matrices
                int i = M_N-1;
                for( int j = 0; j < M_N; ++j )
                    {
                        M_Aq_pr[q]( i, j ) = Aq[q]->energy( M_WN[i], M_WN[j] );
                        M_Aq_du[q]( i, j ) = Aq[q]->energy( M_WNdu[i], M_WNdu[j], true );
                        M_Aq_pr_du[q]( i, j ) = Aq[q]->energy( M_WNdu[i], M_WN[j] );
                    }
                int j = M_N-1;
                for( int i = 0; i < M_N; ++i )
                    {
                        M_Aq_pr[q]( i, j ) = Aq[q]->energy( M_WN[i], M_WN[j] );
                        M_Aq_du[q]( i, j ) = Aq[q]->energy( M_WNdu[i], M_WNdu[j], true );
                        M_Aq_pr_du[q]( i, j ) = Aq[q]->energy( M_WNdu[i], M_WN[j] );
                    }
            }
        std::cout << "  -- Aq computed in " << timer2.elapsed() << "s\n"; timer2.restart();
        Log() << "[CRB::offlineWithErrorEstimation] compute Fq_pr, Fq_du" << "\n";
        for( int q = 0; q < M_model->Ql(0); ++q )
            {
                M_Fq_pr[q].conservativeResize( M_N );
                M_Fq_du[q].conservativeResize( M_N );
                M_Fq_pr[q]( M_N-1 ) = M_model->Fq( 0, q, M_WN[M_N-1] );
                M_Fq_du[q]( M_N-1 ) = M_model->Fq( 0, q, M_WNdu[M_N-1] );
            }
        Log() << "[CRB::offlineWithErrorEstimation] compute Lq_pr, Lq_du" << "\n";
        for( int q = 0; q < M_model->Ql(M_output_index); ++q )
            {
                M_Lq_pr[q].conservativeResize( M_N );
                M_Lq_du[q].conservativeResize( M_N );
                M_Lq_pr[q]( M_N-1 ) = M_model->Fq( M_output_index, q, M_WN[M_N-1] );
                M_Lq_du[q]( M_N-1 ) = M_model->Fq( M_output_index, q, M_WNdu[M_N-1] );
            }
        std::cout << "  -- Fq/Lq computed in " << timer2.elapsed() << "s\n"; timer2.restart();



#if 0
        Log() << "[CRB::offlineWithErrorEstimation] orthonormalize basis functions" << "\n";
        // orthonormalize twice
        if(orthonormalize_primal)
        {
            orthonormalize( M_N, M_WN );
            orthonormalize( M_N, M_WN );
        }
        if(orthonormalize_dual)
        {
            orthonormalize( M_N, M_WNdu );
            orthonormalize( M_N, M_WNdu );
        }
        std::cout << "  -- orthonormalization done in " << timer2.elapsed() << "s\n"; timer2.restart();
#endif// 0

        if (  0 )
        {
            u->zero();
            u->add( 1, M_WN[M_N-1] );
            udu->zero();
            udu->add( 1, M_WNdu[M_N-1] );
            vector_ptrtype Aun( M_backend->newVector( M_model->functionSpace() ) );
            vector_ptrtype Atun( M_backend->newVector( M_model->functionSpace() ) );
            vector_ptrtype Un( M_backend->newVector( M_model->functionSpace() ) );
            vector_ptrtype Undu( M_backend->newVector( M_model->functionSpace() ) );
            vector_ptrtype Frhs( M_backend->newVector( M_model->functionSpace() ) );
            vector_ptrtype Lrhs( M_backend->newVector( M_model->functionSpace() ) );
            *Un = *u;
            *Undu = *udu;
            A->multVector( Un, Aun );
            At->multVector( Undu, Atun );
            Aun->scale( -1 );
            Atun->scale( -1 );
            *Frhs = *F[0];
            *Lrhs = *F[M_output_index];
            Log() << "[CRB::offlineWithErrorEstimation] residual (f,f) " << M_N-1 << ":=" << M_model->scalarProduct( Frhs, Frhs ) << "\n";
            Log() << "[CRB::offlineWithErrorEstimation] residual (f,A) " << M_N-1 << ":=" << 2*M_model->scalarProduct( Frhs, Aun ) << "\n";
            Log() << "[CRB::offlineWithErrorEstimation] residual (A,A) " << M_N-1 << ":=" << M_model->scalarProduct( Aun, Aun ) << "\n";

            Log() << "[CRB::offlineWithErrorEstimation] residual (l,l) " << M_N-1 << ":=" << M_model->scalarProduct( Lrhs, Lrhs ) << "\n";
            Log() << "[CRB::offlineWithErrorEstimation] residual (l,At) " << M_N-1 << ":=" << 2*M_model->scalarProduct( Lrhs, Atun ) << "\n";
            Log() << "[CRB::offlineWithErrorEstimation] residual (At,At) " << M_N-1 << ":=" << M_model->scalarProduct( Atun, Atun ) << "\n";

            Aun->add( *Frhs );
            Lrhs->scale( -1 );
            Atun->add( *Lrhs );
            double err_primal = math::sqrt ( M_model->scalarProduct( Aun, Aun ) );
            double err_dual = math::sqrt ( M_model->scalarProduct( Atun, Atun ) );
            Log() << "[CRB::offlineWithErrorEstimation] primal residual for reduced basis function " << M_N-1 << ":=" << err_primal << "\n";
            Log() << "[CRB::offlineWithErrorEstimation] dual residual for reduced basis function " << M_N-1 << ":=" << err_dual << "\n";
        }

        timer2.restart();
        if(M_error_type==CRB_RESIDUAL || M_error_type == CRB_RESIDUAL_SCM)
            {
                std::cout << "  -- N2Q2 update starts\n";
                generateN2Q2( M_N );
                Log()<<"[CRB::offlineWithErrorEstimation] end of call generateN2Q2 and M_N = "<< M_N <<"\n";
                std::cout << "  -- N2Q2 updated in " << timer2.elapsed() << "s\n"; timer2.restart();
            }

        boost::tie( maxerror, mu, index ) = maxErrorBounds( M_N );
        std::cout << "  -- max error bounds computed in " << timer2.elapsed() << "s\n"; timer2.restart();

        M_rbconv.insert( convergence( M_N, maxerror ) );

        //mu = M_Xi->at( M_N );//M_WNmu_complement->min().get<0>();

        check( M_WNmu->size() );
        if ( this->vm()["crb.check.rb"].template as<int>() == 1 )std::cout << "  -- check reduced basis done in " << timer2.elapsed() << "s\n"; timer2.restart();
        ++M_N;


        std::cout << "time: " << timer.elapsed() << "\n";
        Log() << "time: " << timer.elapsed() << "\n";
        std::cout << "============================================================\n";
        Log() <<"========================================"<<"\n";
    }
    --M_N;
    std::cout<<"number of elements in the reduced basis : "<<M_N<<std::endl;

    this->saveDB();
    std::cout << "Offline CRB is done\n";
}






template<typename TruthModelType>
void
CRB<TruthModelType>::check( size_type N ) const
{
    if ( this->vm()["crb.check.rb"].template as<int>() == 0 )
        return;

    std::cout << "  -- check reduced basis\n";


   Log() << "----------------------------------------------------------------------\n";
    // check that for each mu associated to a basis function of \f$W_N\f$
   for( int k = std::max(0,(int)N-2); k < N; ++k )
    {
        Log() << "**********************************************************************\n";
        parameter_type const& mu = M_WNmu->at( k );
        vectorN_type uN( N );
        vectorN_type uNdu( N );
        double s = lb( N, mu, uN, uNdu );
        double err = delta( N, mu, uN, uNdu );

#if 0
        //if (  err > 1e-5 )
        // {
        std::cout << "[check] error bounds are not < 1e-10\n";
        std::cout << "[check] k = " << k << "\n";
        std::cout << "[check] mu = " << mu << "\n";
        std::cout << "[check] delta = " <<  err << "\n";
        std::cout << "[check] uN( " << k << " ) = " << uN(k) << "\n";
#endif
        // }
        double sfem = M_model->output(M_output_index, mu);

        int size = mu.size();
        std::cout<<"    o mu = [ ";
        for(int i=0;i<size-1;i++) std::cout<< mu[i] <<" , ";
        std::cout<< mu[size-1]<<" ]"<<std::endl;

        Log() << "[check] s= " << s << " +- " << err  << " | sfem= " << sfem << " | abs(sfem-srb) =" << math::abs( sfem - s ) << "\n";
        std::cout << "      + s= " << s << " +- " << err  << " | sfem= " << sfem << " | abs(sfem-srb) =" << math::abs( sfem - s ) << "\n";
    }
    Log() << "----------------------------------------------------------------------\n";

}

template<typename TruthModelType>
void
CRB<TruthModelType>::computeErrorEstimationEfficiencyIndicator (parameterspace_ptrtype const& Dmu, double& max_ei, double& min_ei,int N=4)
{
    vectorN_type uN( N );
    vectorN_type uNdu( N );

    //sampling of parameter space Dmu
    sampling_ptrtype Sampling;
    Sampling = sampling_ptrtype( new sampling_type( M_Dmu ) );
    //Sampling->equidistribute( N );
    Sampling->randomize( N );
    //Sampling->logEquidistribute( N );

    int RBsize = M_WNmu->size();

    y_type ei( Sampling->size() );
    for( size_type k = 0; k < Sampling->size(); ++k ){
        parameter_type const& mu = M_Xi->at( k );
        double s = lb( RBsize, mu, uN, uNdu );//output
        double sfem = M_model->output(M_output_index, mu); //true ouput
        double error_estimation = delta( RBsize,mu,uN,uNdu);
        ei(k) = error_estimation/math::abs(sfem-s);
        std::cout<<" efficiency indicator = "<<ei(k)<<" for parameters {";
        for(int i=0;i<mu.size();i++) std::cout<<mu[i]<<" ";
        std::cout<<"}  --  |sfem - s| = "<<math::abs(sfem-s)<<" and error estimation = "<<error_estimation<<std::endl;
    }

    Eigen::MatrixXf::Index index_max_ei;
    max_ei = ei.array().abs().maxCoeff( &index_max_ei );
    Eigen::MatrixXf::Index index_min_ei;
    min_ei = ei.array().abs().minCoeff( &index_min_ei );

    std::cout<<"[EI] max_ei = "<<max_ei<<" and min_ei = "<<min_ei<<min_ei<<" sampling size = "<<N<<std::endl;

}//end of computeErrorEstimationEfficiencyIndicator





template<typename TruthModelType>
typename CRB<TruthModelType>::value_type
CRB<TruthModelType>::lb( size_type N, parameter_type const& mu, vectorN_type & uN, vectorN_type & uNdu , int K) const
{
    return lb(N,mu,uN,uNdu,mpl::bool_<model_type::is_time_dependent>() , K);
}


template<typename TruthModelType>
typename CRB<TruthModelType>::value_type
CRB<TruthModelType>::lb( size_type N, parameter_type const& mu, vectorN_type & uN, vectorN_type & uNdu, mpl::bool_<true> , int K  ) const
{


    //if K>0 then the time at which we want to evaluate output is defined by
    //time_for_output = K * time_step
    //else it's the default value and in this case we take final time
    double time_for_output;
    double time_step=M_model->timeStep();


    if( K > 0)
    {
        time_for_output = K * time_step;
    }
    else
    {
        int model_K = M_model->computeNumberOfSnapshots();
        time_for_output = model_K * time_step;
    }

    if ( N > M_N ) N = M_N;

    uN.resize( N );

    theta_vector_type theta_aq;
    theta_vector_type theta_mq;
    std::vector<theta_vector_type> theta_fq, theta_lq;

    matrixN_type A ( (int)N, (int)N ) ;
    vectorN_type F ( (int)N );
    vectorN_type L ( (int)N );

    vectorN_type u_old ( (int)N );


    //-- initialization step
    element_ptrtype initial_field ( new element_type ( M_model->functionSpace() ) );
    element_ptrtype projection    ( new element_type ( M_model->functionSpace() ) );
    M_model->initializationField( initial_field ); //fill initial_field

    std::vector<double> coeffs;
    if(M_WN.size()==0)
    {
        throw std::logic_error( "[CRB::lb] ERROR : size of M_WN is zero so the projection of initial condition on the reduced basis is impossible" );
    }

    if( orthonormalize_primal)
    {
        BOOST_FOREACH( auto pr, M_WN )
        {
            double k =  M_model->scalarProduct(*initial_field, pr);
            coeffs.push_back(k);
        }
    }
    else
    {
        matrixN_type MN ( (int)M_N, (int)M_N ) ;
        vectorN_type FN ( (int)M_N );
        for(int i=0; i<M_N; i++)
        {
            for(int j=0; j<i; j++)
            {
                MN(i,j) = M_model->scalarProduct( M_WN[j] , M_WN[i] );
                MN(j,i) = MN(i,j);
            }
            MN(i,i) = M_model->scalarProduct( M_WN[i] , M_WN[i] );
            FN(i) = M_model->scalarProduct(*initial_field,M_WN[i]);
        }
        vectorN_type projectionN ((int) M_N);
        projectionN = MN.lu().solve( FN );

        for(int i=0;i <projectionN.size(); i++)
        {
            coeffs.push_back(projectionN(i));
        }
    }

    int index=0;
    for (std::vector<double>::iterator it = coeffs.begin(); it != coeffs.end(); ++it)
    {
        u_old(index) = *it;
        index++;
    }
    //-- end of initialization step

    int Qm;

    if( M_model->isSteady() )
    {
        time_step = 1e30;
        time_for_output = 1e30;
        Qm = 0;
    }
    else
    {
        Qm = M_model->Qm();
    }

    //vector containing oututs from time=time_step until time=time_for_output
    std::vector<double>output_vector;
    double output;
    for(double time=time_step; time<=time_for_output; time+=time_step)
    {
        boost::tie( theta_mq, theta_aq, theta_fq ) = M_model->computeThetaq( mu ,time);

        A.setZero(N,N);
        for(int q = 0;q < M_model->Qa(); ++q)
        {
            A += theta_aq[q]*M_Aq_pr[q].block(0,0,N,N);
        }

        F.setZero(N);
        for(int q = 0;q < M_model->Ql(0); ++q)
        {
            F += theta_fq[0][q]*M_Fq_pr[q].head(N);
        }

        for(int q = 0;q < Qm; ++q)
        {
            A += theta_mq[q]*M_Mq_pr[q].block(0,0,N,N)/time_step;
            F += theta_mq[q]*M_Mq_pr[q].block(0,0,N,N)*u_old/time_step;
        }

        uN = A.lu().solve( F );

        u_old = uN;

        L.setZero(N);
        for(int q = 0;q < M_model->Ql(M_output_index); ++q)
        {
            L += theta_fq[M_output_index][q]*M_Lq_pr[q].head(N);
        }

        output = L.dot( uN );
        output_vector.push_back(output);

    }

    double s_wo_correction = L.dot( uN );
    return s_wo_correction;
}


template<typename TruthModelType>
typename CRB<TruthModelType>::value_type
CRB<TruthModelType>::lb( size_type N, parameter_type const& mu, vectorN_type & uN, vectorN_type & uNdu, mpl::bool_<false> , int K) const
{
    if ( N > M_N ) N = M_N;
    vectorN_type vN( N );
    //std::cout << "[lb] mu= " << mu << "\n";
    //Log() << "[lb] <compliant> N        = " << N << "\n";
    //Log() << "[lb] <compliant> N        = " << M_N << "\n";
    //Log() << "[lb] <compliant> QLhs     = " << M_model->Qa() << "\n";
    //Log() << "[lb] <compliant> QRhs     = " << M_model->Ql(0) << "\n";

    theta_vector_type theta_aq;
    std::vector<theta_vector_type> theta_fq, theta_lq;
    boost::tie( theta_aq, theta_fq ) = M_model->computeThetaq( mu );

    //std::cout << "theta_aq = " << theta_aq << "\n";
    //std::cout << "theta_fq[0] = " << theta_fq[0] << "\n";
    //std::cout << "theta_fq[" << M_output_index << "] = " << theta_fq[M_output_index] << "\n";
    matrixN_type A( (int)N, (int)N);
    matrixN_type Adu( (int)N, (int)N);
    matrixN_type Aprdu( (int)N, (int)N);
    A.setZero( N, N );
    Adu.setZero( N, N );
    Aprdu.setZero( N, N );
    //std::cout << "[online] A[" << -1 << "]=" << A << "\n";
    for(int q = 0;q < M_model->Qa(); ++q)
    {
        //std::cout << "[online] Aq[" << q << "]=" << M_Aq_pr[q] << "\n";
        A += theta_aq[q]*M_Aq_pr[q].block(0,0,N,N);
        Adu += theta_aq[q]*M_Aq_du[q].block(0,0,N,N);
        Aprdu += theta_aq[q]*M_Aq_pr_du[q].block(0,0,N,N);
        //std::cout << "[lb] A[" << q << "]=" << A << "\n";
        //std::cout << "[lb] Adu[" << q << "]=" << Adu << "\n";
        //std::cout << "[lb] Aprdu[" << q << "]=" << Aprdu << "\n";
    }
    //std::cout << "[lb] A=" << A << "\n";
    //std::cout << "[lb] Adu=" << Adu << "\n";
    //std::cout << "[lb] Aprdu=" << Aprdu << "\n";
    vectorN_type F( (int)N);
    vectorN_type Fdu( (int)N);
    vectorN_type L( (int)N);
    vectorN_type Ldu( (int)N);
    F.setZero(N);
    Fdu.setZero(N);
    L.setZero(N);
    Ldu.setZero(N);
    for(int q = 0;q < M_model->Ql(0); ++q)
    {
        F += theta_fq[0][q]*M_Fq_pr[q].head(N);
        Fdu += theta_fq[0][q]*M_Fq_du[q].head(N);
    }
    for(int q = 0;q < M_model->Ql(M_output_index); ++q)
    {
        L += theta_fq[M_output_index][q]*M_Lq_pr[q].head(N);
        Ldu += theta_fq[M_output_index][q]*M_Lq_du[q].head(N);
    }
    //std::cout << "[lb] F=" << F << "\n";
    //std::cout << "[lb] Fdu=" << Fdu << "\n";
    //std::cout << "[lb] L=" << L << "\n";
    //std::cout << "[lb] Ldu=" << Ldu << "\n";
    // solve A uN = F
    uN.resize( N );
    uNdu.resize( N );
    uN = A.lu().solve( F );
    uNdu = Adu.lu().solve( -Ldu );
    //std::cout << "[lb] un = " << uN << "\n";
    //std::cout << "[lb] undu = " << uNdu << "\n";
    //std::cout << "[lb] energy = " << uN.dot(A*uN) << "\n";
    //std::cout << "[lb] duAprun = " << uNdu.dot(A*uN) << "\n";
    double s = L.dot( uN ) - ( Fdu.dot( uNdu ) - uNdu.dot(Aprdu*uN) );
#if 0
    double s_wo_correction = L.dot( uN );

    double e = N2Q2( N, mu, uN, uNdu );
    ////std::cout << "s = " << s << "\n";
    double sfem = M_model->output(M_output_index, mu);
    Log() << "s = " << s << " swocorr=" << s_wo_correction << " fdotu=" << F.dot( uN )
          << " sfem=" << sfem
          << " error=" << math::abs(sfem-s)
          << " error rb=" << e
          << "\n";
#endif



        return s;


}


template<typename TruthModelType>
typename CRB<TruthModelType>::value_type
CRB<TruthModelType>::delta( size_type N,
                            parameter_type const& mu,
                            vectorN_type const& uN,
                            vectorN_type const& uNdu, int k) const
{
    if ( M_error_type == CRB_NO_RESIDUAL )
      return -1;
    else if( M_error_type == CRB_EMPIRICAL )
        return empiricalErrorEstimation ( N, mu , k);
    else
      return N2Q2( N, mu, uN, uNdu );
}
template<typename TruthModelType>
typename CRB<TruthModelType>::max_error_type
CRB<TruthModelType>::maxErrorBounds( size_type N ) const
{

    vectorN_type uN( N );
    vectorN_type uNdu( N );

    y_type err( M_Xi->size() );

    if( M_error_type == CRB_EMPIRICAL )
    {
      if(N==M_Nm)
      {
          parameter_type mu (M_Dmu);
          size_type id;
          boost::tie (mu, id) = M_Xi->max();
          return boost::make_tuple( 1e5, M_Xi->at( id ), id );
      }
      else{
          y_type err1( M_WNmu_complement->size() );
          for( size_type k = 0; k < M_WNmu_complement->size(); ++k )
              {
                  parameter_type const& mu = M_WNmu_complement->at( k );
                  double _err = delta( N, mu, uN, uNdu, k);
                  err1(k) = _err;
              }
          Eigen::MatrixXf::Index index;
          double maxerr = err1.array().abs().maxCoeff( &index );
          Log() << "[maxErrorBounds] WNmu_complement N=" << N << " max Error = " << maxerr << " at index = " << index << "\n";
          parameter_type  mu = M_WNmu_complement->at( index );
          return boost::make_tuple( maxerr, mu, M_WNmu_complement->indexInSuperSampling( index ) );
      }

    }//end of if ( M_error_type == CRB_EMPIRICAL)
    else
    {
      for( size_type k = 0; k < M_Xi->size(); ++k )
          {
              //std::cout << "--------------------------------------------------\n";
              parameter_type const& mu = M_Xi->at( k );
              //std::cout << "[maxErrorBounds] mu=" << mu << "\n";
              double o = lb( N, mu, uN, uNdu );
              //std::cout << "[maxErrorBounds] output=" << o << "\n";
              double _err = delta( N, mu, uN, uNdu );
              //std::cout << "[maxErrorBounds] error=" << _err << "\n";
              err(k) = _err;
          }
    }//else

    Eigen::MatrixXf::Index index;
    double maxerr = err.array().abs().maxCoeff( &index );
    Log() << "[maxErrorBounds] N=" << N << " max Error = " << maxerr << " at index = " << index << "\n";
    return boost::make_tuple( maxerr, M_Xi->at( index ), index );
}
template<typename TruthModelType>
void
CRB<TruthModelType>::orthonormalize( size_type N, wn_type& wn, int Nm )
{
    std::cout << "  -- orthonormalization (Gram-Schmidt)\n";
    Debug (12000) << "[CRB::orthonormalize] orthonormalize basis for N=" << N << "\n";
    Debug (12000) << "[CRB::orthonormalize] orthonormalize basis for WN="
                  << wn.size() << "\n";
    Debug (12000) << "[CRB::orthonormalize] starting ...\n";

    for( size_type i = 0;i < N;++i )
    {
        for( size_type j = std::max(i+1,N-1); j < N; ++j )
        {
            value_type __rij_pr = M_model->scalarProduct(  wn[i], wn[ j ] );
            wn[j].add( -__rij_pr, wn[i] );
        }

    }
    // normalize
    for( int i =N-Nm;i < N;++i )
    {
        value_type __rii_pr = math::sqrt( M_model->scalarProduct(  wn[i], wn[i] ) );
        wn[i].scale( 1./__rii_pr );
    }

    Debug (12000) << "[CRB::orthonormalize] finished ...\n";
    Debug (12000) << "[CRB::orthonormalize] copying back results in basis\n";

    if ( this->vm()["crb.check.gs"].template as<int>() )
        checkOrthonormality( N , wn);

}

template <typename TruthModelType>
void
CRB<TruthModelType>::checkOrthonormality ( int N, const wn_type& wn) const
{

    if(wn.size()==0)
    {
        throw std::logic_error( "[CRB::checkOrthonormality] ERROR : size of wn is zero" );
    }

    if(orthonormalize_primal*orthonormalize_dual==0)
    {
        std::cout<<"Warning : calling checkOrthonormality is called but ";
        std::cout<<" orthonormalize_dual = "<<orthonormalize_dual;
        std::cout<<" and orthonormalize_primal = "<<orthonormalize_primal<<std::endl;
    }

    matrixN_type A, I;
    A.setZero( N, N );
    I.setIdentity( N, N );
    for( int i = 0;i < N;++i )
    {
        for( int j = 0;j < N;++j )
        {
            A( i, j ) = M_model->scalarProduct(  wn[i], wn[j] );
        }
    }
    A -= I;
    Debug(12000) << "orthonormalization: " << A.norm() << "\n";
    std::cout << "    o check : " << A.norm() << " (should be 0)\n";
    //FEEL_ASSERT( A.norm() < 1e-14 )( A.norm() ).error( "orthonormalization failed.");
}



//error estimation only to build reduced basis
template <typename TruthModelType>
typename CRB<TruthModelType>::value_type
CRB<TruthModelType>::empiricalErrorEstimation ( int Nwn, parameter_type const& mu , int second_index) const
{
  vectorN_type Un( Nwn );
  vectorN_type Undu( Nwn );
  double sn = lb( Nwn, mu, Un, Undu  );

  int nb_element = Nwn/M_factor*(M_factor>0) + (Nwn+M_factor)*(M_factor<0 && Nwn>(-M_factor)) +1*(M_factor<0 && Nwn<=(-M_factor));
  vectorN_type Un2( nb_element );
  vectorN_type Undu2( nb_element );
  double output_smaller_basis = lb(nb_element, mu, Un2, Undu2);
  double error_estimation = math::abs(sn-output_smaller_basis);

  return error_estimation;

}


template<typename TruthModelType>
typename CRB<TruthModelType>::value_type
CRB<TruthModelType>::N2Q2( int Ncur, parameter_type const& mu,vectorN_type const& Un,vectorN_type const& Undu ) const
{
    return N2Q2( Ncur, mu, Un, Undu, mpl::bool_<model_type::is_time_dependent>() );
}

template<typename TruthModelType>
typename CRB<TruthModelType>::value_type
CRB<TruthModelType>::N2Q2( int Ncur,parameter_type const& mu, vectorN_type const& Un,vectorN_type const& Undu, mpl::bool_<true> ) const
{
    //std::cout<<"[CRB::N2Q2] ERROR : N2Q2 need to be implemented for time dependent models"<<std::endl;
    //WARNING works only for steady case
    double time=1e30;

    theta_vector_type theta_aq;
    theta_vector_type theta_mq;
    std::vector<theta_vector_type> theta_fq, theta_lq;
    boost::tie( theta_mq, theta_aq, theta_fq ) = M_model->computeThetaq( mu, time );

    int __QLhs = M_model->Qa();
    int __QRhs = M_model->Ql(0);
    int __QOutput = M_model->Ql(M_output_index);
    int __size = Un.size();
    int __N = Ncur;

    // primal terms
    value_type __c0_pr = 0.0;
    for ( int __q1 = 0;__q1 < __QRhs;++__q1 )
    {
        value_type fq1 = theta_fq[0][__q1];
        for ( int __q2 = 0;__q2 < __QRhs;++__q2 )
        {
            value_type fq2 = theta_fq[0][__q2];
            __c0_pr += M_C0_pr[__q1][__q2]*fq1*fq2;
        }
    }

    value_type __lambda_pr = 0.0;
    value_type __gamma_pr = 0.0;
    for ( int __q1 = 0;__q1 < __QLhs;++__q1 )
    {
        value_type a_q1 = theta_aq[__q1];
        for ( int __q2 = 0;__q2 < __QRhs;++__q2 )
        {
            value_type f_q2 = theta_fq[0][__q2];
            __lambda_pr += a_q1*f_q2*M_Lambda_pr[__q1][__q2].dot( Un );
        }
        for ( int __q2 = 0;__q2 < __QLhs;++__q2 )
        {
            value_type a_q2 = theta_aq[__q2];
            __gamma_pr += a_q1 * a_q2 * Un.transpose()*(M_Gamma_pr[__q1][__q2]*Un);
        }
    }

    // dual terms
    value_type __c0_du = 0.0;
    for ( int __q1 = 0;__q1 < __QOutput;++__q1 )
    {
        value_type fq1 = theta_fq[M_output_index][__q1];
        for ( int __q2 = 0;__q2 < __QOutput;++__q2 )
        {
            value_type fq2 = theta_fq[M_output_index][__q2];
            __c0_du += M_C0_du[__q1][__q2]*fq1*fq2;
        }
    }

    value_type __lambda_du = 0.0;
    value_type __gamma_du = 0.0;
    for ( int __q1 = 0;__q1 < __QLhs;++__q1 )
    {
        value_type a_q1 = theta_aq[__q1];
        for ( int __q2 = 0;__q2 < __QOutput;++__q2 )
        {
            value_type a_q2 = theta_fq[M_output_index][__q2]*a_q1;
            __lambda_du += a_q2 * M_Lambda_du[__q1][__q2].dot(Undu);
        }
        for ( int __q2 = 0;__q2 < __QLhs;++__q2 )
        {
            value_type a_q2 = theta_aq[__q2]*a_q1;
            __gamma_du += a_q2*Undu.transpose()*M_Gamma_du[ __q1][ __q2]*Undu;
        }
    }

    //std::cout << "c0= " << __c0_du << std::endl;
    //std::cout << "lambda= " << __lambda_du << std::endl;
    //std::cout << "gamma= " << __gamma_du << std::endl;
    value_type delta_pr = math::sqrt( math::abs(__c0_pr+__lambda_pr+__gamma_pr) );
    value_type delta_du = math::sqrt( math::abs(__c0_du+__lambda_du+__gamma_du) );
    //std::cout << "delta_pr=" << delta_pr << std::endl;
    //std::cout << "delta_du=" << delta_du << std::endl;

    double alpha = 1;
    if ( M_error_type == CRB_RESIDUAL_SCM )
    {
        double alpha_up, lbti;
        boost::tie( alpha, lbti ) = M_scm->lb( mu );
        boost::tie( alpha_up, lbti ) = M_scm->ub( mu );
        std::cout << "alpha_lo = " << alpha << " alpha_hi = " << alpha_up << "\n";
    }
    return (delta_pr*delta_du)/alpha;
}

template<typename TruthModelType>
typename CRB<TruthModelType>::value_type
CRB<TruthModelType>::N2Q2( int Ncur,parameter_type const& mu, vectorN_type const& Un,vectorN_type const& Undu, mpl::bool_<false> ) const
{
    theta_vector_type theta_aq;
    std::vector<theta_vector_type> theta_fq, theta_lq;
    boost::tie( theta_aq, theta_fq ) = M_model->computeThetaq( mu );

    int __QLhs = M_model->Qa();
    int __QRhs = M_model->Ql(0);
    int __QOutput = M_model->Ql(M_output_index);
    int __size = Un.size();
    int __N = Ncur;

    // primal terms
    value_type __c0_pr = 0.0;
    for ( int __q1 = 0;__q1 < __QRhs;++__q1 )
    {
        value_type fq1 = theta_fq[0][__q1];
        for ( int __q2 = 0;__q2 < __QRhs;++__q2 )
        {
            value_type fq2 = theta_fq[0][__q2];
            __c0_pr += M_C0_pr[__q1][__q2]*fq1*fq2;
        }
    }

    value_type __lambda_pr = 0.0;
    value_type __gamma_pr = 0.0;
    for ( int __q1 = 0;__q1 < __QLhs;++__q1 )
    {
        value_type a_q1 = theta_aq[__q1];
        for ( int __q2 = 0;__q2 < __QRhs;++__q2 )
        {
            value_type f_q2 = theta_fq[0][__q2];
            __lambda_pr += a_q1*f_q2*M_Lambda_pr[__q1][__q2].dot( Un );
        }
        for ( int __q2 = 0;__q2 < __QLhs;++__q2 )
        {
            value_type a_q2 = theta_aq[__q2];
            __gamma_pr += a_q1 * a_q2 * Un.transpose()*(M_Gamma_pr[__q1][__q2]*Un);
        }
    }

    // dual terms
    value_type __c0_du = 0.0;
    for ( int __q1 = 0;__q1 < __QOutput;++__q1 )
    {
        value_type fq1 = theta_fq[M_output_index][__q1];
        for ( int __q2 = 0;__q2 < __QOutput;++__q2 )
        {
            value_type fq2 = theta_fq[M_output_index][__q2];
            __c0_du += M_C0_du[__q1][__q2]*fq1*fq2;
        }
    }

    value_type __lambda_du = 0.0;
    value_type __gamma_du = 0.0;
    for ( int __q1 = 0;__q1 < __QLhs;++__q1 )
    {
        value_type a_q1 = theta_aq[__q1];
        for ( int __q2 = 0;__q2 < __QOutput;++__q2 )
        {
            value_type a_q2 = theta_fq[M_output_index][__q2]*a_q1;
            __lambda_du += a_q2 * M_Lambda_du[__q1][__q2].dot(Undu);
        }
        for ( int __q2 = 0;__q2 < __QLhs;++__q2 )
        {
            value_type a_q2 = theta_aq[__q2]*a_q1;
            __gamma_du += a_q2*Undu.transpose()*M_Gamma_du[ __q1][ __q2]*Undu;
        }
    }

    //std::cout << "c0= " << __c0_du << std::endl;
    //std::cout << "lambda= " << __lambda_du << std::endl;
    //std::cout << "gamma= " << __gamma_du << std::endl;
    value_type delta_pr = math::sqrt( math::abs(__c0_pr+__lambda_pr+__gamma_pr) );
    value_type delta_du = math::sqrt( math::abs(__c0_du+__lambda_du+__gamma_du) );
    //std::cout << "delta_pr=" << delta_pr << std::endl;
    //std::cout << "delta_du=" << delta_du << std::endl;

    double alpha = 1;
    if ( M_error_type == CRB_RESIDUAL_SCM )
    {
        double alpha_up, lbti;
        boost::tie( alpha, lbti ) = M_scm->lb( mu );
        boost::tie( alpha_up, lbti ) = M_scm->ub( mu );
        std::cout << "alpha_lo = " << alpha << " alpha_hi = " << alpha_up << "\n";
    }
    return (delta_pr*delta_du)/alpha;
}


template<typename TruthModelType>
void
CRB<TruthModelType>::generateN2Q2( int Ncur )
{
    return generateN2Q2( Ncur, mpl::bool_<model_type::is_time_dependent>() );
}

template<typename TruthModelType>
void
CRB<TruthModelType>::generateN2Q2( int Ncur, mpl::bool_<true> )
{
    //std::cout<<"[CRB::generateN2Q2]  ERROR : need to be implemented for time dependent problems "<<std::endl;
    //WARNING works only for steady case
    double time=1e30;

    boost::timer ti;
    int __QLhs = M_model->Qa();
    int __QRhs = M_model->Ql(0);
    int __QOutput = M_model->Ql(M_output_index);
    int __size = Ncur;
    int __N = Ncur;
    std::cout << "     o N=" << Ncur << " QLhs=" << __QLhs
              << " QRhs=" << __QRhs << " Qoutput=" << __QOutput << "\n";
    vector_ptrtype __X( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __Fdu( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __Xdu( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __Y( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __Ydu( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __Z_pr(  M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __Z2_pr(  M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __W_pr(  M_backend->newVector( M_model->functionSpace() ) );
    namespace ublas = boost::numeric::ublas;

    std::vector<sparse_matrix_ptrtype> Aq;
    std::vector<sparse_matrix_ptrtype> Mq;
    std::vector<std::vector<vector_ptrtype> > Fq,Lq;
    boost::tie( Mq, Aq, Fq ) = M_model->computeAffineDecomposition();
    __X->zero();
    __X->add( 1.0 );
    //std::cout << "measure of domain= " << M_model->scalarProduct( __X, __X ) << "\n";
#if 0
    ublas::vector<value_type> mu(P);
    for ( int q = 0; q < P; ++q )
    {
        mu[q] = mut(0.0);
    }
#endif

    // Primal
    // no need to recompute this term each time
    if ( Ncur == 1 )
    {
        Log() << "[generateN2Q2] Compute Primal residual data\n";
        Log() << "[generateN2Q2] C0_pr\n";
        // see above X = C^-1 F and Y = F
        for ( int __q1 = 0;__q1 < __QRhs;++__q1 )
        {
            //Log() << "__Fq->norm1=" << Fq[0][__q1]->l2Norm() << "\n";
            M_model->l2solve( __X, Fq[0][__q1] );
            //for ( int __q2 = 0;__q2 < __q1;++__q2 )
            for ( int __q2 = 0;__q2 < __QRhs;++__q2 )
            {
                //Log() << "__Fq->norm 2=" << Fq[0][__q2]->l2Norm() << "\n";
                M_model->l2solve( __Y, Fq[0][__q2] );
                //M_C0_pr[__q1][__q2] = M_model->scalarProduct( __X, Fq[0][__q2] );
                M_C0_pr[__q1][__q2] = M_model->scalarProduct( __X, __Y );
                //M_C0_pr[__q2][__q1] = M_C0_pr[__q1][__q2];
                //Debug() << "M_C0_pr[" << __q1 << "][" << __q2 << "]=" << M_C0_pr[__q1][__q2] << "\n";
                //Log() << "M_C0_pr[" << __q1 << "][" << __q2 << "]=" << M_C0_pr[__q1][__q2] << "\n";
            }
            //M_C0_pr[__q1][__q1] = M_model->scalarProduct( __X, __X );
        }
    }
    std::cout << "     o initialize generateN2Q2 in " << ti.elapsed() << "s\n"; ti.restart();

#if 0
    parameter_type const& mu = M_WNmu->at( 0 );
    //std::cout << "[generateN2Q2] mu=" << mu << "\n";
    theta_vector_type theta_aq;
    std::vector<theta_vector_type> theta_fq;
    boost::tie( theta_aq, theta_fq ) = M_model->computeThetaq( mu );
    value_type __c0_pr = 0.0;
    for ( int __q1 = 0;__q1 < __QRhs;++__q1 )
    {
        value_type fq1 = theta_fq[0][__q1];
        for ( int __q2 = 0;__q2 < __QRhs;++__q2 )
        {
            value_type fq2 = theta_fq[0][__q2];
            __c0_pr += M_C0_pr[__q1][__q2]*fq1*fq2;
        }
    }

    //std::cout << "c0=" << __c0_pr << std::endl;

    sparse_matrix_ptrtype A;
    sparse_matrix_ptrtype M;
    std::vector<vector_ptrtype> F;
    boost::tie( M, A, F ) = M_model->update( mu, time );
    M_model->l2solve( __X, F[0] );
    //std::cout << "c0 2 = " << M_model->scalarProduct( __X, __X ) << "\n";;
#endif

    //
    //  Primal
    //
    Log() << "[generateN2Q2] Lambda_pr, Gamma_pr\n";
    for ( int __q1 = 0;__q1 < __QLhs;++__q1 )
    {
        for ( int __q2 = 0; __q2 < __QRhs;++__q2 )
        {
            M_Lambda_pr[__q1][__q2].conservativeResize( __N );

            *__X=M_WN[__N-1];
            Aq[__q1]->multVector(  __X, __W_pr );
            __W_pr->scale( -1. );
            //std::cout << "__W_pr->norm=" << __W_pr->l2Norm() << "\n";
            M_model->l2solve( __Z_pr, __W_pr );

            //__Y = Fq[0][__q2];
            //std::cout << "__Fq->norm=" << Fq[0][__q2]->l2Norm() << "\n";
            M_model->l2solve( __Y, Fq[0][__q2] );
            M_Lambda_pr[ __q1][ __q2](__N-1) = 2.0*M_model->scalarProduct( __Y, __Z_pr );
            //Debug() << "M_Lambda_pr[" << __q1 << "][" << __q2 << "][" << __j << "]=" << M_Lambda_pr[__q1][__q2][__j] << "\n";
            //std::cout << "M_Lambda_pr[" << __q1 << "][" << __q2 << "][" << __j << "]=" << M_Lambda_pr[__q1][__q2][__j] << "\n";
        }
    }
    std::cout << "     o Lambda_pr updated in " << ti.elapsed() << "s\n"; ti.restart();
    for ( int __q1 = 0;__q1 < __QLhs;++__q1 )
    {
        for ( int __q2 = 0; __q2 < __QLhs;++__q2 )
        {
            M_Gamma_pr[__q1][__q2].conservativeResize( __N, __N );

            // line N-1
            int __j=__N-1;
            *__X=M_WN[__j];
            Aq[__q1]->multVector(  __X, __W_pr );
            __W_pr->scale( -1. );
            //std::cout << "__W_pr->norm=" << __W_pr->l2Norm() << "\n";
            M_model->l2solve( __Z_pr, __W_pr );

            //Aq[__q2]->multVector(  __Z_pr, __W_pr );
            for ( int __l = 0; __l < ( int )__N;++__l )
            {
                *__X=M_WN[__l];
                Aq[__q2]->multVector(  __X, __W_pr );
                __W_pr->scale( -1. );
                //std::cout << "__W2_pr->norm=" << __W_pr->l2Norm() << "\n";
                M_model->l2solve( __Z2_pr, __W_pr );
                M_Gamma_pr[ __q1][ __q2](__j,__l) = M_model->scalarProduct( __Z_pr, __Z2_pr );
                //M_Gamma_pr[ __q2][ __q1][ __j ][__l] = M_Gamma_pr[ __q1][ __q2][ __j ][__l];
                //Debug() << "M_Gamma_pr[" << __q1 << "][" << __q2 << "][" << __j << "][" << __l << "]=" << M_Gamma_pr[__q1][__q2][__j][__l] << "\n";
                //std::cout << "M_Gamma_pr[" << __q1 << "][" << __q2 << "][" << __j << "][" << __l << "]=" << M_Gamma_pr[__q1][__q2][__j][__l] << "\n";
            }

            for ( int __j = 0;__j < ( int )__N;++__j )
            {
                *__X=M_WN[__j];
                Aq[__q1]->multVector(  __X, __W_pr );
                __W_pr->scale( -1. );
                //std::cout << "__W_pr->norm=" << __W_pr->l2Norm() << "\n";
                M_model->l2solve( __Z_pr, __W_pr );

                //Aq[__q2]->multVector(  __Z_pr, __W_pr );
                //column N-1
                int __l = __N-1;
                {
                    *__X=M_WN[__l];
                    Aq[__q2]->multVector(  __X, __W_pr );
                    __W_pr->scale( -1. );
                    //std::cout << "__W2_pr->norm=" << __W_pr->l2Norm() << "\n";
                    M_model->l2solve( __Z2_pr, __W_pr );
                    M_Gamma_pr[ __q1][ __q2](__j,__l) = M_model->scalarProduct( __Z_pr, __Z2_pr );
                    //M_Gamma_pr[ __q2][ __q1][ __j ][__l] = M_Gamma_pr[ __q1][ __q2][ __j ][__l];
                    //Debug() << "M_Gamma_pr[" << __q1 << "][" << __q2 << "][" << __j << "][" << __l << "]=" << M_Gamma_pr[__q1][__q2][__j][__l] << "\n";
                    //std::cout << "M_Gamma_pr[" << __q1 << "][" << __q2 << "][" << __j << "][" << __l << "]=" << M_Gamma_pr[__q1][__q2][__j][__l] << "\n";
                }
            }

        }// on N1
    } // on q1
    std::cout << "     o Gamma_pr updated in " << ti.elapsed() << "s\n"; ti.restart();
    sparse_matrix_ptrtype Atq1 = M_model->newMatrix();
    sparse_matrix_ptrtype Atq2 = M_model->newMatrix();
    //
    // Dual
    //
    // compute this only once
    if ( Ncur == 1 )
    {
        Log() << "[generateN2Q2] Compute Dual residual data\n";
        Log() << "[generateN2Q2] C0_du\n";
        for ( int __q1 = 0;__q1 < __QOutput;++__q1 )
        {
            *__Fdu = *Fq[M_output_index][__q1];
            __Fdu->scale(-1.0);
            M_model->l2solve( __Xdu, __Fdu );
            for ( int __q2 = 0;__q2 < __QOutput;++__q2 )
            {
                *__Fdu = *Fq[M_output_index][__q2];
                __Fdu->scale(-1.0);
                M_model->l2solve( __Ydu, __Fdu );
                M_C0_du[__q1][__q2] = M_model->scalarProduct( __Xdu, __Ydu );
                //M_C0_du[__q2][__q1] = M_C0_du[__q1][__q2];
            }
            //M_C0_du[__q1][__q1] = M_model->scalarProduct( __Xdu, __Xdu );
        }
        std::cout << "     o C0_du updated in " << ti.elapsed() << "s\n"; ti.restart();
    }

    Log() << "[generateN2Q2] Lambda_du, Gamma_du\n";
    for ( int __q1 = 0;__q1 < __QLhs;++__q1 )
    {
        Aq[__q1]->transpose( Atq1 );
        for ( int __q2 = 0; __q2 < __QOutput;++__q2 )
        {
            M_Lambda_du[__q1][__q2].conservativeResize( __N );

            *__Xdu=M_WNdu[__N-1];
            Atq1->multVector(  __Xdu, __W_pr );
            __W_pr->scale( -1. );
            //std::cout << "__W_pr->norm=" << __W_pr->l2Norm() << "\n";
            M_model->l2solve( __Z_pr, __W_pr );

            *__Fdu = *Fq[M_output_index][__q2];
            __Fdu->scale(-1.0);
            M_model->l2solve( __Y, __Fdu );
            M_Lambda_du[ __q1][ __q2](__N-1) = 2.0*M_model->scalarProduct( __Y, __Z_pr );
            //Debug() << "M_Lambda_pr[" << __q1 << "][" << __q2 << "][" << __j << "]=" << M_Lambda_pr[__q1][__q2][__j] << "\n";
            //std::cout << "M_Lambda_pr[" << __q1 << "][" << __q2 << "][" << __j << "]=" << M_Lambda_pr[__q1][__q2][__j] << "\n";
        } // q2
    } // q2
    std::cout << "     o Lambda_du updated in " << ti.elapsed() << "s\n"; ti.restart();
    for ( int __q1 = 0;__q1 < __QLhs;++__q1 )
    {
        Aq[__q1]->transpose( Atq1 );

        for ( int __q2 = 0; __q2 < __QLhs;++__q2 )
        {
            Aq[__q2]->transpose( Atq2 );
            M_Gamma_du[__q1][__q2].conservativeResize( __N, __N );
            // update line N-1
            {
                int __j = __N-1;
                *__Xdu=M_WNdu[__j];
                Atq1->multVector(  __Xdu, __W_pr );
                __W_pr->scale( -1. );
                //std::cout << "__W_pr->norm=" << __W_pr->l2Norm() << "\n";
                M_model->l2solve( __Z_pr, __W_pr );

                //Aq[__q2]->multVector(  __Z_pr, __W_pr );

                for ( int __l = 0; __l < ( int )__N;++__l )
                {
                    *__X=M_WNdu[__l];
                    Atq2->multVector(  __X, __W_pr );
                    __W_pr->scale( -1. );
                    //std::cout << "__W2_pr->norm=" << __W_pr->l2Norm() << "\n";
                    M_model->l2solve( __Z2_pr, __W_pr );
                    M_Gamma_du[ __q1][ __q2](__j,__l) = M_model->scalarProduct( __Z_pr, __Z2_pr );
                    //M_Gamma_pr[ __q2][ __q1][ __j ][__l] = M_Gamma_pr[ __q1][ __q2][ __j ][__l];
                    //Debug() << "M_Gamma_pr[" << __q1 << "][" << __q2 << "][" << __j << "][" << __l << "]=" << M_Gamma_pr[__q1][__q2][__j][__l] << "\n";
                    //std::cout << "M_Gamma_pr[" << __q1 << "][" << __q2 << "][" << __j << "][" << __l << "]=" << M_Gamma_pr[__q1][__q2][__j][__l] << "\n";
                }
            }

            // update column __N-1
            for ( int __j = 0;__j < ( int )__N;++__j )
            {
                *__Xdu=M_WNdu[__j];
                Atq1->multVector(  __Xdu, __W_pr );
                __W_pr->scale( -1. );
                //std::cout << "__W_pr->norm=" << __W_pr->l2Norm() << "\n";
                M_model->l2solve( __Z_pr, __W_pr );

                //Aq[__q2]->multVector(  __Z_pr, __W_pr );
                int __l = __N-1;
                {
                    *__X=M_WNdu[__l];
                    Atq2->multVector(  __X, __W_pr );
                    __W_pr->scale( -1. );
                    //std::cout << "__W2_pr->norm=" << __W_pr->l2Norm() << "\n";
                    M_model->l2solve( __Z2_pr, __W_pr );
                    M_Gamma_du[ __q1][ __q2](__j,__l) = M_model->scalarProduct( __Z_pr, __Z2_pr );
                    //M_Gamma_pr[ __q2][ __q1][ __j ][__l] = M_Gamma_pr[ __q1][ __q2][ __j ][__l];
                    //Debug() << "M_Gamma_pr[" << __q1 << "][" << __q2 << "][" << __j << "][" << __l << "]=" << M_Gamma_pr[__q1][__q2][__j][__l] << "\n";
                    //std::cout << "M_Gamma_pr[" << __q1 << "][" << __q2 << "][" << __j << "][" << __l << "]=" << M_Gamma_pr[__q1][__q2][__j][__l] << "\n";
                }
            }

        }// on N1
    } // on q1
    std::cout << "     o Gamma_du updated in " << ti.elapsed() << "s\n"; ti.restart();
    Log() << "[generateN2Q2] Done.\n";


}

template<typename TruthModelType>
void
CRB<TruthModelType>::generateN2Q2( int Ncur, mpl::bool_<false> )
{
    boost::timer ti;
    int __QLhs = M_model->Qa();
    int __QRhs = M_model->Ql(0);
    int __QOutput = M_model->Ql(M_output_index);
    int __size = Ncur;
    int __N = Ncur;
    std::cout << "     o N=" << Ncur << " QLhs=" << __QLhs
              << " QRhs=" << __QRhs << " Qoutput=" << __QOutput << "\n";
    vector_ptrtype __X( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __Fdu( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __Xdu( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __Y( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __Ydu( M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __Z_pr(  M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __Z2_pr(  M_backend->newVector( M_model->functionSpace() ) );
    vector_ptrtype __W_pr(  M_backend->newVector( M_model->functionSpace() ) );
    namespace ublas = boost::numeric::ublas;

    std::vector<sparse_matrix_ptrtype> Aq;
    std::vector<std::vector<vector_ptrtype> > Fq,Lq;
    boost::tie( Aq, Fq ) = M_model->computeAffineDecomposition();
    __X->zero();
    __X->add( 1.0 );
    //std::cout << "measure of domain= " << M_model->scalarProduct( __X, __X ) << "\n";
#if 0
    ublas::vector<value_type> mu(P);
    for ( int q = 0; q < P; ++q )
    {
        mu[q] = mut(0.0);
    }
#endif

    // Primal
    // no need to recompute this term each time
    if ( Ncur == 1 )
    {
        Log() << "[generateN2Q2] Compute Primal residual data\n";
        Log() << "[generateN2Q2] C0_pr\n";
        // see above X = C^-1 F and Y = F
        for ( int __q1 = 0;__q1 < __QRhs;++__q1 )
        {
            //Log() << "__Fq->norm1=" << Fq[0][__q1]->l2Norm() << "\n";
            M_model->l2solve( __X, Fq[0][__q1] );
            //for ( int __q2 = 0;__q2 < __q1;++__q2 )
            for ( int __q2 = 0;__q2 < __QRhs;++__q2 )
            {
                //Log() << "__Fq->norm 2=" << Fq[0][__q2]->l2Norm() << "\n";
                M_model->l2solve( __Y, Fq[0][__q2] );
                //M_C0_pr[__q1][__q2] = M_model->scalarProduct( __X, Fq[0][__q2] );
                M_C0_pr[__q1][__q2] = M_model->scalarProduct( __X, __Y );
                //M_C0_pr[__q2][__q1] = M_C0_pr[__q1][__q2];
                //Debug() << "M_C0_pr[" << __q1 << "][" << __q2 << "]=" << M_C0_pr[__q1][__q2] << "\n";
                //Log() << "M_C0_pr[" << __q1 << "][" << __q2 << "]=" << M_C0_pr[__q1][__q2] << "\n";
            }
            //M_C0_pr[__q1][__q1] = M_model->scalarProduct( __X, __X );
        }
    }
    std::cout << "     o initialize generateN2Q2 in " << ti.elapsed() << "s\n"; ti.restart();

#if 0
    parameter_type const& mu = M_WNmu->at( 0 );
    //std::cout << "[generateN2Q2] mu=" << mu << "\n";
    theta_vector_type theta_aq;
    std::vector<theta_vector_type> theta_fq;
    boost::tie( theta_aq, theta_fq ) = M_model->computeThetaq( mu );
    value_type __c0_pr = 0.0;
    for ( int __q1 = 0;__q1 < __QRhs;++__q1 )
    {
        value_type fq1 = theta_fq[0][__q1];
        for ( int __q2 = 0;__q2 < __QRhs;++__q2 )
        {
            value_type fq2 = theta_fq[0][__q2];
            __c0_pr += M_C0_pr[__q1][__q2]*fq1*fq2;
        }
    }

    //std::cout << "c0=" << __c0_pr << std::endl;

    sparse_matrix_ptrtype A;
    std::vector<vector_ptrtype> F;
    boost::tie( A, F ) = M_model->update( mu );
    M_model->l2solve( __X, F[0] );
    //std::cout << "c0 2 = " << M_model->scalarProduct( __X, __X ) << "\n";;
#endif

    //
    //  Primal
    //
    Log() << "[generateN2Q2] Lambda_pr, Gamma_pr\n";
    for ( int __q1 = 0;__q1 < __QLhs;++__q1 )
    {
        for ( int __q2 = 0; __q2 < __QRhs;++__q2 )
        {
            M_Lambda_pr[__q1][__q2].conservativeResize( __N );

            *__X=M_WN[__N-1];
            Aq[__q1]->multVector(  __X, __W_pr );
            __W_pr->scale( -1. );
            //std::cout << "__W_pr->norm=" << __W_pr->l2Norm() << "\n";
            M_model->l2solve( __Z_pr, __W_pr );

            //__Y = Fq[0][__q2];
            //std::cout << "__Fq->norm=" << Fq[0][__q2]->l2Norm() << "\n";
            M_model->l2solve( __Y, Fq[0][__q2] );
            M_Lambda_pr[ __q1][ __q2](__N-1) = 2.0*M_model->scalarProduct( __Y, __Z_pr );
            //Debug() << "M_Lambda_pr[" << __q1 << "][" << __q2 << "][" << __j << "]=" << M_Lambda_pr[__q1][__q2][__j] << "\n";
            //std::cout << "M_Lambda_pr[" << __q1 << "][" << __q2 << "][" << __j << "]=" << M_Lambda_pr[__q1][__q2][__j] << "\n";
        }
    }
    std::cout << "     o Lambda_pr updated in " << ti.elapsed() << "s\n"; ti.restart();
    for ( int __q1 = 0;__q1 < __QLhs;++__q1 )
    {
        for ( int __q2 = 0; __q2 < __QLhs;++__q2 )
        {
            M_Gamma_pr[__q1][__q2].conservativeResize( __N, __N );

            // line N-1
            int __j=__N-1;
            *__X=M_WN[__j];
            Aq[__q1]->multVector(  __X, __W_pr );
            __W_pr->scale( -1. );
            //std::cout << "__W_pr->norm=" << __W_pr->l2Norm() << "\n";
            M_model->l2solve( __Z_pr, __W_pr );

            //Aq[__q2]->multVector(  __Z_pr, __W_pr );
            for ( int __l = 0; __l < ( int )__N;++__l )
            {
                *__X=M_WN[__l];
                Aq[__q2]->multVector(  __X, __W_pr );
                __W_pr->scale( -1. );
                //std::cout << "__W2_pr->norm=" << __W_pr->l2Norm() << "\n";
                M_model->l2solve( __Z2_pr, __W_pr );
                M_Gamma_pr[ __q1][ __q2](__j,__l) = M_model->scalarProduct( __Z_pr, __Z2_pr );
                //M_Gamma_pr[ __q2][ __q1][ __j ][__l] = M_Gamma_pr[ __q1][ __q2][ __j ][__l];
                //Debug() << "M_Gamma_pr[" << __q1 << "][" << __q2 << "][" << __j << "][" << __l << "]=" << M_Gamma_pr[__q1][__q2][__j][__l] << "\n";
                //std::cout << "M_Gamma_pr[" << __q1 << "][" << __q2 << "][" << __j << "][" << __l << "]=" << M_Gamma_pr[__q1][__q2][__j][__l] << "\n";
            }

            for ( int __j = 0;__j < ( int )__N;++__j )
            {
                *__X=M_WN[__j];
                Aq[__q1]->multVector(  __X, __W_pr );
                __W_pr->scale( -1. );
                //std::cout << "__W_pr->norm=" << __W_pr->l2Norm() << "\n";
                M_model->l2solve( __Z_pr, __W_pr );

                //Aq[__q2]->multVector(  __Z_pr, __W_pr );
                //column N-1
                int __l = __N-1;
                {
                    *__X=M_WN[__l];
                    Aq[__q2]->multVector(  __X, __W_pr );
                    __W_pr->scale( -1. );
                    //std::cout << "__W2_pr->norm=" << __W_pr->l2Norm() << "\n";
                    M_model->l2solve( __Z2_pr, __W_pr );
                    M_Gamma_pr[ __q1][ __q2](__j,__l) = M_model->scalarProduct( __Z_pr, __Z2_pr );
                    //M_Gamma_pr[ __q2][ __q1][ __j ][__l] = M_Gamma_pr[ __q1][ __q2][ __j ][__l];
                    //Debug() << "M_Gamma_pr[" << __q1 << "][" << __q2 << "][" << __j << "][" << __l << "]=" << M_Gamma_pr[__q1][__q2][__j][__l] << "\n";
                    //std::cout << "M_Gamma_pr[" << __q1 << "][" << __q2 << "][" << __j << "][" << __l << "]=" << M_Gamma_pr[__q1][__q2][__j][__l] << "\n";
                }
            }

        }// on N1
    } // on q1
    std::cout << "     o Gamma_pr updated in " << ti.elapsed() << "s\n"; ti.restart();
    sparse_matrix_ptrtype Atq1 = M_model->newMatrix();
    sparse_matrix_ptrtype Atq2 = M_model->newMatrix();
    //
    // Dual
    //
    // compute this only once
    if ( Ncur == 1 )
    {
        Log() << "[generateN2Q2] Compute Dual residual data\n";
        Log() << "[generateN2Q2] C0_du\n";
        for ( int __q1 = 0;__q1 < __QOutput;++__q1 )
        {
            *__Fdu = *Fq[M_output_index][__q1];
            __Fdu->scale(-1.0);
            M_model->l2solve( __Xdu, __Fdu );
            for ( int __q2 = 0;__q2 < __QOutput;++__q2 )
            {
                *__Fdu = *Fq[M_output_index][__q2];
                __Fdu->scale(-1.0);
                M_model->l2solve( __Ydu, __Fdu );
                M_C0_du[__q1][__q2] = M_model->scalarProduct( __Xdu, __Ydu );
                //M_C0_du[__q2][__q1] = M_C0_du[__q1][__q2];
            }
            //M_C0_du[__q1][__q1] = M_model->scalarProduct( __Xdu, __Xdu );
        }
        std::cout << "     o C0_du updated in " << ti.elapsed() << "s\n"; ti.restart();
    }

    Log() << "[generateN2Q2] Lambda_du, Gamma_du\n";
    for ( int __q1 = 0;__q1 < __QLhs;++__q1 )
    {
        Aq[__q1]->transpose( Atq1 );
        for ( int __q2 = 0; __q2 < __QOutput;++__q2 )
        {
            M_Lambda_du[__q1][__q2].conservativeResize( __N );

            *__Xdu=M_WNdu[__N-1];
            Atq1->multVector(  __Xdu, __W_pr );
            __W_pr->scale( -1. );
            //std::cout << "__W_pr->norm=" << __W_pr->l2Norm() << "\n";
            M_model->l2solve( __Z_pr, __W_pr );

            *__Fdu = *Fq[M_output_index][__q2];
            __Fdu->scale(-1.0);
            M_model->l2solve( __Y, __Fdu );
            M_Lambda_du[ __q1][ __q2](__N-1) = 2.0*M_model->scalarProduct( __Y, __Z_pr );
            //Debug() << "M_Lambda_pr[" << __q1 << "][" << __q2 << "][" << __j << "]=" << M_Lambda_pr[__q1][__q2][__j] << "\n";
            //std::cout << "M_Lambda_pr[" << __q1 << "][" << __q2 << "][" << __j << "]=" << M_Lambda_pr[__q1][__q2][__j] << "\n";
        } // q2
    } // q2
    std::cout << "     o Lambda_du updated in " << ti.elapsed() << "s\n"; ti.restart();
    for ( int __q1 = 0;__q1 < __QLhs;++__q1 )
    {
        Aq[__q1]->transpose( Atq1 );

        for ( int __q2 = 0; __q2 < __QLhs;++__q2 )
        {
            Aq[__q2]->transpose( Atq2 );
            M_Gamma_du[__q1][__q2].conservativeResize( __N, __N );
            // update line N-1
            {
                int __j = __N-1;
                *__Xdu=M_WNdu[__j];
                Atq1->multVector(  __Xdu, __W_pr );
                __W_pr->scale( -1. );
                //std::cout << "__W_pr->norm=" << __W_pr->l2Norm() << "\n";
                M_model->l2solve( __Z_pr, __W_pr );

                //Aq[__q2]->multVector(  __Z_pr, __W_pr );

                for ( int __l = 0; __l < ( int )__N;++__l )
                {
                    *__X=M_WNdu[__l];
                    Atq2->multVector(  __X, __W_pr );
                    __W_pr->scale( -1. );
                    //std::cout << "__W2_pr->norm=" << __W_pr->l2Norm() << "\n";
                    M_model->l2solve( __Z2_pr, __W_pr );
                    M_Gamma_du[ __q1][ __q2](__j,__l) = M_model->scalarProduct( __Z_pr, __Z2_pr );
                    //M_Gamma_pr[ __q2][ __q1][ __j ][__l] = M_Gamma_pr[ __q1][ __q2][ __j ][__l];
                    //Debug() << "M_Gamma_pr[" << __q1 << "][" << __q2 << "][" << __j << "][" << __l << "]=" << M_Gamma_pr[__q1][__q2][__j][__l] << "\n";
                    //std::cout << "M_Gamma_pr[" << __q1 << "][" << __q2 << "][" << __j << "][" << __l << "]=" << M_Gamma_pr[__q1][__q2][__j][__l] << "\n";
                }
            }

            // update column __N-1
            for ( int __j = 0;__j < ( int )__N;++__j )
            {
                *__Xdu=M_WNdu[__j];
                Atq1->multVector(  __Xdu, __W_pr );
                __W_pr->scale( -1. );
                //std::cout << "__W_pr->norm=" << __W_pr->l2Norm() << "\n";
                M_model->l2solve( __Z_pr, __W_pr );

                //Aq[__q2]->multVector(  __Z_pr, __W_pr );
                int __l = __N-1;
                {
                    *__X=M_WNdu[__l];
                    Atq2->multVector(  __X, __W_pr );
                    __W_pr->scale( -1. );
                    //std::cout << "__W2_pr->norm=" << __W_pr->l2Norm() << "\n";
                    M_model->l2solve( __Z2_pr, __W_pr );
                    M_Gamma_du[ __q1][ __q2](__j,__l) = M_model->scalarProduct( __Z_pr, __Z2_pr );
                    //M_Gamma_pr[ __q2][ __q1][ __j ][__l] = M_Gamma_pr[ __q1][ __q2][ __j ][__l];
                    //Debug() << "M_Gamma_pr[" << __q1 << "][" << __q2 << "][" << __j << "][" << __l << "]=" << M_Gamma_pr[__q1][__q2][__j][__l] << "\n";
                    //std::cout << "M_Gamma_pr[" << __q1 << "][" << __q2 << "][" << __j << "][" << __l << "]=" << M_Gamma_pr[__q1][__q2][__j][__l] << "\n";
                }
            }

        }// on N1
    } // on q1
    std::cout << "     o Gamma_du updated in " << ti.elapsed() << "s\n"; ti.restart();
    Log() << "[generateN2Q2] Done.\n";

}
template<typename TruthModelType>
boost::tuple<double,double,double>
CRB<TruthModelType>::run( parameter_type const& mu, double eps )
{

    int Nwn = M_N;
    if( M_error_type!=CRB_EMPIRICAL )
    {
        auto lo = M_rbconv.right.range( boost::bimaps::unbounded, boost::bimaps::_key <= eps );
        for( auto it = lo.first; it != lo.second; ++it )
           {
               std::cout << "rbconv[" << it->first <<"]=" << it->second << "\n";
           }
        auto it = M_rbconv.project_left( lo.first );
        Nwn = it->first;
        std::cout << "Nwn = "<< Nwn << " error = "<< it->second << " eps=" << eps << "\n";
    }
    vectorN_type uN( Nwn );
    vectorN_type uNdu( Nwn );
    uN.setZero( Nwn );
    uNdu.setZero( Nwn );
    double o = lb( Nwn, mu, uN, uNdu  );
    double e = delta( Nwn, mu, uN, uNdu );
    return boost::make_tuple( o , e, Nwn );
}

template<typename TruthModelType>
void
CRB<TruthModelType>::run( const double * X, unsigned long N, double * Y, unsigned long P )
{

    parameter_type mu( M_Dmu );
    // the last parameter is the max error
    for( int p= 0; p < N-4; ++p )
        mu(p) = X[p];

#if 0
    std::cout<<" list of parameters : [";
    for(int i=0;i<N-1;i++) std::cout<<X[i]<<" , ";
    std::cout<<X[N-1]<<" ] "<<std::endl;
#endif
    //double meshSize  = X[N-4];
    //M_model->setMeshSize(meshSize);
    setOutputIndex( (int)X[N-4] );
    int Nwn =  X[N-3];
    int maxerror = X[N-2];
    CRBErrorType errorType =(CRBErrorType)X[N-1];
    setCRBErrorType(errorType);

    auto lo = M_rbconv.right.range( boost::bimaps::unbounded,boost::bimaps::_key <= maxerror );
#if 0
   for( auto it = lo.first; it != lo.second; ++it )
        {
            std::cout << "rbconv[" << it->first <<"]=" << it->second << "\n";
        }
#endif
    auto it = M_rbconv.project_left( lo.first );
    Nwn = it->first;
    //std::cout << "Nwn = "<< Nwn << " error = "<< it->second << " maxerror=" << maxerror << " ErrorType = "<<errorType <<"\n";

    //int Nwn = M_WN.size();
    vectorN_type uN( Nwn );
    vectorN_type uNdu( Nwn );
    uN.setZero( Nwn );
    uNdu.setZero( Nwn );


    FEEL_ASSERT( P == 2 )( P ).warn( "invalid number of outputs" );
    Y[0] = lb( Nwn, mu, uN, uNdu  );
    Y[1] = delta( Nwn, mu, uN, uNdu );

}




template<typename TruthModelType>
void
CRB<TruthModelType>::projectionOnPodSpace( const element_ptrtype & u , element_ptrtype& projection, const std::string& name_of_space)
{

    projection->zero();

    if( name_of_space=="dual")
    {
        if( orthonormalize_dual )
        //in this case we can simplify because elements of reduced basis are orthonormalized
        {
            BOOST_FOREACH( auto du, M_WNdu )
            {
                element_type e;
                e = du;
                double k =  M_model->scalarProduct(*u, e);
                e.scale(k);
                projection->add( 1 , e );
            }
        }//end of if orthonormalize_dual
        else
        {
            matrixN_type MN ( (int)M_N, (int)M_N ) ;
            vectorN_type FN ( (int)M_N );
            for(int i=0; i<M_N; i++)
            {
                for(int j=0; j<i; j++)
                {
                    MN(i,j) = M_model->scalarProduct( M_WNdu[j] , M_WNdu[i] );
                    MN(j,i) = MN(i,j);
                }
                MN(i,i) = M_model->scalarProduct( M_WNdu[i] , M_WNdu[i] );
                FN(i) = M_model->scalarProduct(*u,M_WNdu[i]);
            }
            vectorN_type projectionN ((int) M_N);
            projectionN = MN.lu().solve( FN );
            int index=0;
            BOOST_FOREACH( auto du, M_WNdu )
            {
                element_type e;
                e = du;
                double k =  projectionN(index);
                e.scale(k);
                projection->add( 1 , e );
                index++;
            }
       }//end of if( ! orthonormalize_dual )
    }//end of projection on dual POD space
    else
    {
        if( orthonormalize_primal)
        {
            BOOST_FOREACH( auto pr, M_WN )
            {
                element_type e;
                e = pr;
                double k =  M_model->scalarProduct(*u, e);
                e.scale(k);
                projection->add( 1 , e );
            }
        }//end of if orthonormalize_primal
        else
       {
           matrixN_type MN ( (int)M_N, (int)M_N ) ;
           vectorN_type FN ( (int)M_N );
           for(int i=0; i<M_N; i++)
           {
                for(int j=0; j<i; j++)
                {
                    MN(i,j) = M_model->scalarProduct( M_WN[j] , M_WN[i] );
                    MN(j,i) = MN(i,j);
                }
                MN(i,i) = M_model->scalarProduct( M_WN[i] , M_WN[i] );
                FN(i) = M_model->scalarProduct(*u,M_WN[i]);
            }
            vectorN_type projectionN ((int) M_N);
            projectionN = MN.lu().solve( FN );
            int index=0;
            BOOST_FOREACH( auto pr, M_WN )
            {
                element_type e;
                e = pr;
                double k =  projectionN(index);
                e.scale(k);
                projection->add( 1 , e );
                index++;
            }
       }//end of if( ! orthonormalize_primal )
    }//end of projection on primal POD space

}





template<typename TruthModelType>
template<class Archive>
void
CRB<TruthModelType>::save(Archive & ar, const unsigned int version) const
{
    ar & boost::serialization::base_object<super>(*this);
    ar & BOOST_SERIALIZATION_NVP( M_output_index );
    ar & BOOST_SERIALIZATION_NVP( M_N );
    ar & BOOST_SERIALIZATION_NVP( M_rbconv );
    ar & BOOST_SERIALIZATION_NVP( M_error_type );
    ar & BOOST_SERIALIZATION_NVP( M_Xi );
    ar & BOOST_SERIALIZATION_NVP( M_WNmu );
    ar & BOOST_SERIALIZATION_NVP( M_Aq_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Aq_du );
    ar & BOOST_SERIALIZATION_NVP( M_Aq_pr_du );
    ar & BOOST_SERIALIZATION_NVP( M_Fq_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Fq_du );
    ar & BOOST_SERIALIZATION_NVP( M_Lq_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Lq_du );

    ar & BOOST_SERIALIZATION_NVP( M_C0_pr );
    ar & BOOST_SERIALIZATION_NVP( M_C0_du );
    ar & BOOST_SERIALIZATION_NVP( M_Lambda_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Lambda_du );
    ar & BOOST_SERIALIZATION_NVP( M_Gamma_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Gamma_du );

    if( model_type::is_time_dependent )
    {
        ar & BOOST_SERIALIZATION_NVP( M_Mq_pr );
        ar & BOOST_SERIALIZATION_NVP( M_Mq_du );
        ar & BOOST_SERIALIZATION_NVP( M_Mq_pr_du );
    }

}

template<typename TruthModelType>
template<class Archive>
void
CRB<TruthModelType>::load(Archive & ar, const unsigned int version)
{
    ar & boost::serialization::base_object<super>(*this);
    ar & BOOST_SERIALIZATION_NVP( M_output_index );
    ar & BOOST_SERIALIZATION_NVP( M_N );
    ar & BOOST_SERIALIZATION_NVP( M_rbconv );
    ar & BOOST_SERIALIZATION_NVP( M_error_type );
    ar & BOOST_SERIALIZATION_NVP( M_Xi );
    ar & BOOST_SERIALIZATION_NVP( M_WNmu );
    ar & BOOST_SERIALIZATION_NVP( M_Aq_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Aq_du );
    ar & BOOST_SERIALIZATION_NVP( M_Aq_pr_du );
    ar & BOOST_SERIALIZATION_NVP( M_Fq_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Fq_du );
    ar & BOOST_SERIALIZATION_NVP( M_Lq_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Lq_du );

    ar & BOOST_SERIALIZATION_NVP( M_C0_pr );
    ar & BOOST_SERIALIZATION_NVP( M_C0_du );
    ar & BOOST_SERIALIZATION_NVP( M_Lambda_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Lambda_du );
    ar & BOOST_SERIALIZATION_NVP( M_Gamma_pr );
    ar & BOOST_SERIALIZATION_NVP( M_Gamma_du );

    if( model_type::is_time_dependent )
    {
        ar & BOOST_SERIALIZATION_NVP( M_Mq_pr );
        ar & BOOST_SERIALIZATION_NVP( M_Mq_du );
        ar & BOOST_SERIALIZATION_NVP( M_Mq_pr_du );
    }


#if 0
    std::cout << "[loadDB] output index : " << M_output_index << "\n"
              << "[loadDB] N : " << M_N << "\n"
              << "[loadDB] error type : " << M_error_type << "\n";
    for( auto it = M_rbconv.begin(), en = M_rbconv.end();
         it != en; ++it )
    {
        std::cout << "[loadDB] convergence: (" << it->left << ","  << it->right  << ")\n";
    }
#endif
}

template<typename TruthModelType>
void
CRB<TruthModelType>::saveDB()
{
    fs::ofstream ofs( this->dbLocalPath() / this->dbFilename() );
    if ( ofs )
    {
        boost::archive::text_oarchive oa(ofs);
        // write class instance to archive
        oa << *this;
        // archive and stream closed when destructors are called
    }
}
template<typename TruthModelType>
bool
CRB<TruthModelType>::loadDB()
{
    fs::path db = this->lookForDB();
    if ( db.empty() )
        return false;
    if ( !fs::exists( db ) )
        return false;
    //std::cout << "Loading " << db << "...\n";
    fs::ifstream ifs( db );

    if ( ifs )
    {
        boost::archive::text_iarchive ia(ifs);
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
#endif /* __CRB_H */
