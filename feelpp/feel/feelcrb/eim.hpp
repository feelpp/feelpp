/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-05-02

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

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
   \file eim.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Cecile Daversin <daversin@math.unistra.fr>
   \author Stephane Veys
   \date 2012-05-02
 */
#ifndef _FEELPP_EIM_HPP
#define _FEELPP_EIM_HPP 1

#include <limits>
#include <numeric>
#include <string>
#include <iostream>
#include <fstream>

#include <boost/ref.hpp>
#include <boost/next_prior.hpp>
#include <boost/type_traits.hpp>
#include <boost/tuple/tuple.hpp>
#if BOOST_VERSION >= 104700
//#include <boost/math/special_functions/nonfinite_num_facets.hpp>
#endif
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/base_object.hpp>

#include <feel/feelcrb/crbdb.hpp>
#include <feel/feelcrb/parameterspace.hpp>

//#include <feel/feelvf/vf.hpp>

#include <feel/feelcrb/crb.hpp>
#include <feel/feelcrb/crbmodel.hpp>
#include <feel/feeldiscr/reducedbasisspace.hpp>
#include <feel/feeldiscr/geometricspace.hpp>

//#include <Eigen/Core>

namespace Feel
{
class ModelCrbBaseBase;
class EimFunctionNoSolveBase {};

/**
  \class EIM
  \brief Empirical interpolation of a function to obtain an affine decomposition

  We are given a function \f$g (\: \cdot\: ; \mu) \in L^{\infty}
  (\Omega)\f$ of sufficient regularity.


  \section offline Offline Stage

  To begin, we choose \f$\mu^g_1\f$,
  and define \f$S^g_1 = \{ \mu^g_1 \}\f$, \f$\xi_1 \equiv g (x ; \mu^g_1)\f$, and
  \f$W^g_1 = {\rm span} \: \{\xi_1 \}\f$; we assume that \f$\xi_1  \neq 0\f$.
  Then, for \f$M \geq 2\f$, we set \f$\mu^g_M = \arg \max_{\mu \in
  \Xi^{^g}}\inf _{z \in W^g_{M-1}} \|g (\: \cdot \: ; \mu) - z
  \|_{L^{\infty} (\Omega)}\f$, where \f$\Xi^g\f$ is a suitably fine parameter
  sample over \f${\mathcal{D}}\f$. We then set \f$S^g_M = S^g_{M-1} \cup \mu^g_M\f$,
  \f$\xi_M = g (x;\mu^g_M)\f$, and
  \f$W^g_M = {\rm span} \: \{ \xi_m, \: 1 \leq m \leq M \}\f$. Note that,
  thanks to our truth
  approximation, \f$\mu^g_M\f$ is the solution of a <em>standard linear
  program</em>.

  We suppose that \f$M_{\max}\f$ is chosen such that the dimension of \f$\{g
  (\: \cdot \: ;
  \mu) \: | \: \mu \in
  \mathcal{D}\} \f$ exceeds \f$M_{\max}\f$.

  We now construct nested sets of interpolation points \f$T_M = \{ t_1,
  \ldots, t_M \}\f$, \f$1 \leq M
  \leq M_{\max}\f$.  We first set \f$t_1 = \arg \: {\rm ess} \: \sup_{x \in
  \Omega} | \xi_1 (x)|\f$,
  \f$q_1 = \xi_1 (x) / \xi_1 (t_1) \f$, \f$B^1_{11} = 1\f$.  Then for \f$M = 2,
  \ldots, M_{\max}\f$, we solve
  the linear system
  \f$ \sum^{M-1}_{j = 1} \: \sigma^{M-1}_j \: q_j(t_i) = \xi_M (t_i)\f$,
  \f$ 1 \leq i \leq M-1\f$, and set \f$r_M (x) = \xi_M (x) - \sum^{M-1}_{j =
  1}\: \sigma^{M-1}_j \:
  q_j (x)\f$, \f$t_M = \arg \: {\rm ess} \:
  \sup_{x \in \Omega} |r_M (x)|\f$, \f$q_M (x) = r_M (x) /r_M (t_M) \f$, and
  \f$B^M_{i \: j} = q_j (t_i)\f$,
  \f$1 \leq i,j \leq M\f$.

  \section online Online Stage

  Our coefficient function approximation is the interpolant of
  \f$g\f$ over \f$T_M\f$ :
  \f$g_M (x ;
  \mu) = \sum^M_{m = 1} \beta_m
  (\mu) \: q_m (x)\f$, where
  \f$\sum^M_{j = 1} \: B^M_{i \: j} \: \beta_j (\mu) = g (t_i ; \mu)\f$, \f$
  1 \leq i \leq M\f$.
  We define \f$\varepsilon_M (\mu) \equiv \| g (\: \cdot \: ; \mu) - g_M
  (\: \cdot \: ; \mu)\|_{L^{\infty} (\Omega)}\f$.

  \todo make it truly mesh independent.
  \todo find a generic solution for coordinates type \c node_type

  @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  @see
*/
template<typename ModelType>
class EIM
{

public:

    static const uint16_type nDim = ModelType::nDim;
    /** @name Typedefs
     */
    //@{
    typedef ModelType model_type;
    typedef typename ModelType::value_type value_type;

    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_type;
    typedef typename matrix_type::ColXpr column_type;
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector_type;
    typedef typename ModelType::node_type node_type;

    typedef typename ModelType::functionspace_type functionspace_type;
    typedef std::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef typename functionspace_type::element_type element_type;
    typedef typename functionspace_type::element_ptrtype element_ptrtype;

    typedef typename ModelType::model_functionspace_type model_functionspace_type;
    typedef std::shared_ptr<model_functionspace_type> model_functionspace_ptrtype;
    typedef typename model_functionspace_type::element_type model_element_type;
    typedef typename model_functionspace_type::element_type model_solution_type;

    typedef typename ModelType::parameterspace_type parameterspace_type;
    typedef typename ModelType::parameter_type parameter_type;
    typedef typename parameterspace_type::sampling_ptrtype sampling_ptrtype;
    typedef typename parameterspace_type::sampling_type sampling_type;

    typedef boost::tuple<double,Eigen::Matrix<double,nDim,1> > space_residual_type;
    typedef boost::tuple<double,parameter_type> parameter_residual_type;

    typedef Eigen::VectorXd vectorN_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    EIM( std::string const& name= "default", WorldComm const& worldComm = Environment::worldComm(), std::string const& prefix = "", std::string const& prefixModel = "" )
        :
        M_prefix( prefix ),
        M_prefixModel( prefixModel ),
        M_worldComm( worldComm ),
        M_is_read( false ),
        M_is_written( false ),
        M_name( name ),
        M_M( 1 ),
        M_offline_done( false ),
        M_offline_step( false ),
        M_adapt_SER( false ),
        M_tol( 1e-8 ),
        M_index_max(),
        M_model( 0 ),
        M_correct_RB_SER( false )
        {}
    EIM( model_type* model, sampling_ptrtype sampling, double __tol = 1e-8, bool offline_done=false, std::string const& prefix = "", std::string const& prefixModel = "" )
        :
        M_prefix( prefix ),
        M_prefixModel( prefixModel ),
        M_worldComm( model->worldComm() ),
        M_is_read( false ),
        M_is_written( false ),
        M_name( model->name() ),
        M_trainset( sampling ),
        M_M( 1 ),
        M_offline_done( offline_done ),
        M_offline_step( false ),
        M_adapt_SER( false ),
        M_tol( __tol ),
        M_index_max(),
        M_model( model ),
        M_correct_RB_SER( false )
        {

            int user_max = ioption(_prefix=this->M_prefix,_name="eim.dimension-max");
            int max_built = M_model->maxQ();
            bool enrich_database = boption(_prefix=this->M_prefix,_name="eim.enrich-database");
            bool cobuild = ( ioption(_prefix=this->M_prefixModel,_name="ser.eim-frequency") != 0 );

            bool do_offline=false;
            M_restart=false;

            if( (user_max+1) > max_built && M_offline_done )
            {
                if( enrich_database || cobuild)
                    do_offline=true;
                else
                    do_offline=false;
            }

            if( !M_offline_done || boption(_prefix=this->M_prefix,_name="eim.rebuild-database") )
            {
                do_offline=true;
                M_restart = true;
            }

            if ( do_offline )
            {
                LOG(INFO) << "construct EIM approximation...\n";

                if( enrich_database )
                {
                    if( this->worldComm().isMasterRank() )
                        std::cout<<model->name()<<" enrich the existing database..."<<std::endl;
                }
                else if( cobuild )
                {
                    if( this->worldComm().isMasterRank() )
                        std::cout<<model->name()<<" continue co-building process..."<<std::endl;
                }
                if( M_restart )
                {
                    M_model->initializeDataStructures();
                }

                this->setOfflineStep( true );
                offline();
            }
        }

    EIM( EIM const & __bbf )
        :
        M_is_read( __bbf.M_is_read ),
        M_is_written( __bbf.M_is_written ),
        M_name( __bbf.M_name ),
        M_M( __bbf.M_M ),
        M_offline_done( __bbf.M_offline_done ),
        M_tol( __bbf.M_tol ),
        M_index_max( __bbf.M_index_max ),
        M_model( __bbf.M_model ),
        M_correct_RB_SER( __bbf.M_correct_RB_SER )
        {}
    ~EIM()
        {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    //! \return the mpi communicators
    WorldComm const& worldComm() const { return M_worldComm; }

    /**
       \return the number of DOF in space discretization
    */
    size_type nDOF() const {  FEELPP_ASSERT( M_model != 0 ).error( "Invalid EIM model" ); return M_model->functionSpace()->nLocalDof(); }

    /**
     * return the set of reduced basis functions associated with the eim
     */
    std::vector<element_type> const& q() const {  return M_model->q(); }


    /**
       Check whether the offline stage has been executed and the database created.
       \return \c true if the offline has been executed and the DB saved, \c false otherwise
    */
    bool isOfflineDone() const
        {
            return M_offline_done;
        }

    //@}

    void setOfflineStep( bool b )
        {
            M_offline_step = b;
        }
    bool offlineStep() const
        {
            return M_offline_step;
        }
    void setAdaptationSER( bool b )
        {
            M_adapt_SER = b;
        }
    bool adaptationSER() const
        {
            return M_adapt_SER;
        }
    void setRbCorrection( bool b )
        {
            M_correct_RB_SER = b;
        }
    bool rbCorrection() const
        {
            return M_correct_RB_SER;
        }
    void setRestart( bool b )
        {
            M_restart = b;
        }
    bool restart() const
        {
            return M_restart;
        }

    /** @name  Mutators
     */
    //@{

    void setTrainSet( sampling_ptrtype pset ) { M_trainset = pset; }

    //@}

    /** @name  Methods
     */
    //@{

    /**
       \brief Online stage of the coefficient-function interpolation.

       Our coefficient function approximation is the interpolant of
       \f$g\f$ over \f$T_M\f$ :
       \f$g_M (x ;
       \mu) = \sum^M_{m = 1} \beta_m
       (\mu) \: q_m (x)\f$, where
       \f$\sum^M_{j = 1} \: B^M_{i \: j} \: \beta_j (\mu) = g (t_i ; \mu)\f$, \f$
       1 \leq i \leq M\f$.
       We define \f$\varepsilon_M (\mu) \equiv \| g (\: \cdot \: ; \mu) - g_M
       (\: \cdot \: ; \mu)\|_{L^{\infty} (\Omega)}\f$.

       Note that \f$ \Omega \f$ is given and $D^\mu$ is handled by the \c parameterset_type
       data structure.
    */
    vector_type beta( parameter_type const& mu  ) const { return M_model->beta( mu, M_model->mMax() ); }
    vector_type beta( parameter_type const& mu, model_solution_type const& T  ) const { return M_model->beta( mu, T, M_model->mMax() ); }
    vector_type beta( parameter_type const& mu, size_type M  ) const { return M_model->beta( mu , M ); }
    vector_type beta( parameter_type const& mu, model_solution_type const& T, size_type M  ) const {return M_model->beta( mu , T , M ); }

    void studyConvergence( parameter_type const & mu, model_solution_type & solution , std::vector< std::string > all_file_name ) const;
    boost::tuple<double,element_type> interpolationErrorEstimation( parameter_type const & mu, model_solution_type const& solution , int M ) const ;
    double errorEstimationLinf( parameter_type const & mu, model_solution_type const& solution , int M ) const ;

    void computationalTimeStatistics( std::string appname )  { return M_model->computationalTimeStatistics(); }
    element_type residual ( size_type M ) const;

    parameter_residual_type computeBestFit( sampling_ptrtype trainset, int __M);

    element_type operator()( parameter_type const& mu , int N) const { return expansion( M_model->q(), M_model->beta( mu, N ) , N); }
    element_type operator()( parameter_type const& mu, model_solution_type const& T , int N ) const { return expansion( M_model->q(), M_model->beta( mu, T, N ) , N ); }

    /**
       orthonormalize
    */
    void orthonormalize( std::vector<element_type>& );

    void offline();
    //@}
protected:

    std::string M_prefix;
    std::string M_prefixModel;

    //! mpi communicators
    WorldComm const& M_worldComm;

    mutable bool M_is_read;
    mutable bool M_is_written;

    std::string M_name;
    sampling_ptrtype M_trainset;
    size_type M_M;

    size_type M_max_q;//size of vector M_q ( to save/load )

    mutable bool M_offline_done;
    mutable bool M_offline_step;
    mutable bool M_adapt_SER;

    double M_tol;

    std::vector<size_type> M_index_max;

    model_type* M_model;

    mutable bool M_correct_RB_SER;
    bool M_restart;
    double M_greedy_maxerr;
    double M_greedy_rbmaxerr; // stores maximum value of rb approximation error indicator (only for SER with error estimation)
    bool M_criterion;
protected:

private:

    /**
       \brief Offline stage of the coefficient-function interpolation.

       We are given a function \f$g (\: \cdot\: ; \mu) \in L^{\infty}
       (\Omega)\f$ of sufficient regularity.  To begin, we choose \f$\mu^g_1\f$,
       and define \f$S^g_1 = \{ \mu^g_1 \}\f$, \f$\xi_1 \equiv g (x ; \mu^g_1)\f$, and
       \f$W^g_1 = {\rm span} \: \{\xi_1 \}\f$; we assume that \f$\xi_1  \neq 0\f$.
       Then, for \f$M \geq 2\f$, we set \f$\mu^g_M = \arg \max_{\mu \in
       \Xi^{^g}}\inf _{z \in W^g_{M-1}} \|g (\: \cdot \: ; \mu) - z
       \|_{L^{\infty} (\Omega)}\f$, where \f$\Xi^g\f$ is a suitably fine parameter
       sample over \f${\mathcal{D}}\f$. We then set \f$S^g_M = S^g_{M-1} \cup \mu^g_M\f$,
       \f$\xi_M = g (x;\mu^g_M)\f$, and
       \f$W^g_M = {\rm span} \: \{ \xi_m, \: 1 \leq m \leq M \}\f$. Note that,
       thanks to our truth
       approximation, \f$\mu^g_M\f$ is the solution of a <em>standard linear
       program</em>.

       We suppose that \f$M_{\max}\f$ is chosen such that the dimension of \f$\{g
       (\: \cdot \: ;
       \mu) \: | \: \mu \in
       \mathcal{D}\} \f$ exceeds \f$M_{\max}\f$.

       We now construct nested sets of interpolation points \f$T_M = \{ t_1,
       \ldots, t_M \}\f$, \f$1 \leq M
       \leq M_{\max}\f$.  We first set \f$t_1 = \arg \: {\rm ess} \: \sup_{x \in
       \Omega} | \xi_1 (x)|\f$,
       \f$q_1 = \xi_1 (x) / \xi_1 (t_1) \f$, \f$B^1_{11} = 1\f$.  Then for \f$M = 2,
       \ldots, M_{\max}\f$, we solve
       the linear system
       \f$ \sum^{M-1}_{j = 1} \: \sigma^{M-1}_j \: q_j(t_i) = \xi_M (t_i)\f$,
       \f$ 1 \leq i \leq M-1\f$, and set \f$r_M (x) = \xi_M (x) - \sum^{M-1}_{j =
       1}\: \sigma^{M-1}_j \:
       q_j (x)\f$, \f$t_M = \arg \: {\rm ess} \:
       \sup_{x \in \Omega} |r_M (x)|\f$, \f$q_M (x) = r_M (x) /r_M (t_M) \f$, and
       \f$B^M_{i \: j} = q_j (t_i)\f$,
       \f$1 \leq i,j \leq M\f$.
    */
    //void offline();
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template <typename ModelType>
typename boost::tuple< double, typename EIM<ModelType>::element_type >
EIM<ModelType>::interpolationErrorEstimation( parameter_type const & mu, model_solution_type const& solution , int M ) const
{
    double max = M_model->mMax();
    CHECK( M <= max ) << "Invalid number M for errorEstimation: " << M << " Mmax : " << max << "\n";
    //interpolation point associated to (M+1) basis function is index by M
    //remember that interpolationPoint(0) is the first interpolation point
    //auto t = M_model->interpolationPoint( M );
    //auto projected_expression = M_model->operator()( solution , t,  mu );

    auto projected_expression_at_points = M_model->evaluateExpressionAtInterpolationPoints( solution ,  mu , M+1);
    double projected_expression = projected_expression_at_points(M);
    auto eim_approximation = this->operator()(mu, solution, M);
    //double eim = eim_approximation(t)(0,0,0);
    auto eim_approximation_at_points = M_model->evaluateElementAtInterpolationPoints( eim_approximation , M+1 );
    double eim = eim_approximation_at_points(M);

    double coeff = math::abs( projected_expression - eim );

    auto result = M_model->q( M );
    result.scale( coeff );
    return boost::make_tuple(coeff,result);
}


template <typename ModelType>
double
EIM<ModelType>::errorEstimationLinf( parameter_type const & mu, model_solution_type const& solution , int M ) const
{
    double max = M_model->mMax();
    CHECK( M <= max ) << "Invalid number M for errorEstimation: " << M << " Mmax : " << max << "\n";
    auto projected_expression = M_model->operator()( solution  , mu );
    auto eim_approximation = this->operator()(mu, solution, M);
    auto diff = idv( projected_expression ) - idv( eim_approximation );
    auto norm = normLinf( _range=M_model->functionSpace()->template rangeElements<0>(), _pset=_Q<0>(), _expr= diff );
    double error = norm.template get<0>();
    return error;
}

template<typename ModelType>
void
EIM<ModelType>::orthonormalize( std::vector<element_type> & __Z )
{
    size_type __M = __Z.size();
    for ( size_type __i = 0;__i < __M-1; ++__i )
    {
        value_type __s = inner_product(__Z[__i],__Z[__M-1]);
        __Z[__M-1].add(- __s,__Z[__i]);
    }
    __Z[__M-1].scale(inner_product(__Z[__M-1],__Z[__M-1]));
}


template<typename ModelType>
typename EIM<ModelType>::parameter_residual_type
EIM<ModelType>::computeBestFit( sampling_ptrtype trainset, int __M)
{
    boost::mpi::timer timer;
    DVLOG(2) << "compute best fit  for m=" << __M
             << " and trainset of size " << trainset->size() << "...\n";
    using namespace vf;
    int proc_number =  this->worldComm().globalRank();

    vector_type maxerr( trainset->size() );
    maxerr.setZero();
    int index = 0;
    int index_criterion = 0;
    DVLOG(2) << "Compute best fit M=" << __M << "\n";
    vector_type rhs( __M );

    model_solution_type solution = M_model->modelFunctionSpace()->element();
    double time_solve = 0;
    double time_exp = 0;

    // Create the associated trainset (according to ser.error-estimation)
    sampling_ptrtype subtrainset( new sampling_type( M_model->parameterSpace() ) );
    bool ser_error_estimation = boption(_prefix=this->M_prefixModel,_name="ser.error-estimation");
    int subtrainset_method = ioption(_prefix=this->M_prefixModel,_name="ser.eim-subtrainset-method");

    if( ser_error_estimation && subtrainset_method != 0 )
        subtrainset = M_model->createSubTrainset( trainset, subtrainset_method );
    else
        subtrainset = trainset;

    std::vector<bool> error_criterion(subtrainset->size(), true);
    M_criterion=true;
    double max_rb_error = 0.0;
    bool ser = (ioption(_prefix=this->M_prefixModel,_name="ser.eim-frequency") != 0);
    double rtol = doption(_prefix=this->M_prefixModel,_name="ser.eim-greedy-rtol");
    std::vector<vectorN_type> uN; //eventually contains reduced basis approx.
    tic();
    for( auto const& mu : *subtrainset )
    {
#if !defined(NDEBUG)
        DVLOG(2) << "compute best fit check mu...\n";
        mu.check();
#endif

        //LOG_EVERY_N(INFO, 1 ) << " (every 10 mu) compute fit at mu="<< mu <<"\n" ;
        timer.restart();

        if ( M_model->modelUseSolve() )
        {
            if ( ser ) //Use SER
            {
                if( ser_error_estimation && doption(_prefix=this->M_prefixModel,_name="ser.eim-greedy-rtol")!=0 )
                {
                    // Compute the error indicator with mu
                    auto riesz = M_model->RieszResidualNorm( mu );
                    double error_indicator = riesz.template get<0>();
                    uN = riesz.template get<1>();
                    if( error_indicator > max_rb_error ) // Update the max error indicator to be used at the next step
                        max_rb_error = error_indicator;
                    if( M_greedy_rbmaxerr != 0 ) // Takes whole sampling with the roughest RB approx (greedy_maxerr uninitialized)
                        error_criterion[index] = error_indicator/M_greedy_rbmaxerr < rtol;
                }

                if( boption(_prefix=this->M_prefixModel,_name="ser.use-rb-in-eim-mu-selection") && error_criterion[index] )
                    solution = M_model->computeRbExpansion( mu, uN ); //RB
                else
                    solution = M_model->computePfem( mu ); //PFEM
            }
            else
                solution = M_model->solve( mu ); //FEM
            time_solve += timer.elapsed();
            timer.restart();
        }

        rhs = M_model->computeExpansionCoefficients( mu , solution, __M );
        auto z = expansion( M_model->q(), rhs, __M );
        auto resmax = M_model->computeMaximumOfResidual( mu , solution, z );
        time_exp += timer.elapsed();
        //DCHECK( rhs.size() == __M ) << "Invalid size rhs: " << rhs.size() << " M=" << __M  << " rhs = " << rhs << "\n";

        DCHECK( index < subtrainset->size() ) << "Invalid index " << index << " should be less than trainset size = " << subtrainset->size() << "\n";
        maxerr( index ) = resmax.template get<0>();
        if( error_criterion[index] )
            index_criterion++;
        index++;
    }
    toc("ComputeBestFit Greedy Algorithm - "+ M_model->name(), FLAGS_v>0);

    if( doption(_prefix=this->M_prefixModel,_name="ser.eim-greedy-rtol")!=0 )
    {
        M_greedy_rbmaxerr = max_rb_error; // Update the maximum error indicator value for the next basis
        if( this->worldComm().isMasterRank() )
            std::cout << "[computeBestFit] updated M_greedy_rbmaxerr = " << M_greedy_rbmaxerr << std::endl;
    }

    if( this->worldComm().isMasterRank() )
    {
        std::cout << "-- Mean time to solve(mu) : " << time_solve/trainset->size() << std::endl;
        std::cout << "-- Mean time to compute resmax : " << time_exp/trainset->size() << std::endl;
        if( doption(_prefix=this->M_prefixModel,_name="ser.eim-greedy-rtol")!=0 )
        {
            std::cout << "-- Number of parameters selected to compute resmax : " << index << "/" << subtrainset->size() << "\n";
            std::cout << "-- Number of parameters which satistify the criterion : " << index_criterion << "/" << subtrainset->size() << "\n";
        }
    }

    //LOG_ASSERT( index == subtrainset->size() ) << "Invalid index " << index << " should be equal to trainset size = " << subtrainset->size() << "\n";
    LOG_ASSERT( index <= subtrainset->size() ) << "Invalid index " << index << " should be inferior to trainset size = " << subtrainset->size() << "\n";
    auto err = maxerr.array().abs().maxCoeff( &index );

    if( this->worldComm().isMasterRank() )
        std::cout << "err=" << err << " reached at index " << index << " and mu=" << subtrainset->at(index) << "\n";
    if( ! error_criterion[index] )
        M_criterion = false;

    return boost::make_tuple( err, subtrainset->at(index) );
}

template<typename ModelType>
void
EIM<ModelType>::offline()
{
    using namespace vf;

    if ( this->worldComm().isMasterRank() )
    {
        if ( !fs::exists( M_model->dbLocalPath() ) )
            fs::create_directories( M_model->dbLocalPath() );
    }
    this->worldComm().barrier();

    bool expression_expansion = M_model->computeExpansionOfExpression();

    int max_z=0;
    int max_solution=0;
    int max_g=0;
    auto mu = M_model->parameterSpace()->element();
    node_type t;
    auto solution = M_model->modelFunctionSpace()->element();
    double time=0;
    double time_=0;
    boost::mpi::timer timer,timer2,timer3;
    double rb_online_mean_iterations, rb_online_max_increments;

    if( M_restart )
    {
        tic();
        if( this->worldComm().isMasterRank() )
        {
            std::cout<<" ******************offline EIM for expression "<<M_model->name()<<" start "<<std::endl;
        }

        DVLOG(2) << "[offline] starting offline stage ...\n";
        M_M = 1;
        M_max_q=0;
        if ( !M_trainset )
            M_trainset = M_model->parameterSpace()->sampling();
        if ( M_trainset->empty() )
        {
            int sampling_size = ioption(_prefix=this->M_prefix,_name="eim.sampling-size");
            std::string file_name = ( boost::format("eim_trainset_%1%") % sampling_size ).str();
            bool all_procs_have_same_sampling=true;
            std::string sampling_mode = "log-equidistribute";// random, log-random, log-equidistribute, equidistribute
            std::ifstream file ( file_name );
            if( ! file )
            {
                M_trainset->sample( sampling_size, sampling_mode, all_procs_have_same_sampling, file_name );
            }
            else
            {
                M_trainset->clear();
                M_trainset->readFromFile(file_name);
            }
        }

        DVLOG(2) << "[offline] create mu_1...\n";

        // min element in Dmu to start with (in // each proc have the same element)
        mu = M_model->parameterSpace()->max();
        mu.check();

        DVLOG( 2 ) << "mu ( of size "<<mu.size()<<"): \n"<<mu;

        //store this value
        M_model->clearParameterSampling();
        M_model->addParameter( mu );

        DVLOG( 2 ) <<" parameter added";
        timer2.restart();

        solution = M_model->solve( mu );
        time=timer2.elapsed();
        if( this->worldComm().isMasterRank() )
        {
            std::cout<<" -- model solution computed in "<<time<<"s"<<std::endl;
        }
        DVLOG( 2 ) << "solution computed";

        if( !expression_expansion )
        {
            M_model->addExpressionEvaluation( M_model->operator()( solution , mu ) );
            max_g++;
        }
        else
        {
            M_model->addSolution( solution );
            max_solution++;
        }

        if( this->worldComm().isMasterRank() )
            std::cout << "compute finite element solution at mu_1 done";
        VLOG(2) << "compute finite element solution at mu_1 done";

        DVLOG(2) << "compute T^" << 0 << "...\n";
        // Build T^0

        timer2.restart();
        auto zmax = M_model->computeMaximumOfExpression( mu , solution );
        // store space coordinate where max absolute value occurs
        t = zmax.template get<1>();
        time=timer2.elapsed();
        if( this->worldComm().isMasterRank() )
        {
            std::cout<<" -- maximum of expression computed in "<<time<<"s"<<std::endl;
        }

        M_model->addInterpolationPoint( t );
        DVLOG( 2 )<<"add the interpolation point : \n"<<t;
        DVLOG( 2 ) << "norm Linf = " << zmax.template get<0>() << " at " << zmax.template get<1>() << "\n";

        auto zero = vf::project( _space=M_model->functionSpace(), _range=M_model->functionSpace()->template rangeElements<0>(), _expr=cst(0.) );
        if( expression_expansion )
        {
            M_model->addZ( zero );
            max_z++;
        }

        timer2.restart();
        auto q = M_model->Residual(0,zero);
        time=timer2.elapsed();
        if( this->worldComm().isMasterRank() )
        {
            std::cout<<" -- expression evaluated in mu in "<<time<<"s"<<std::endl;
        }
        M_max_q++;
        DVLOG( 2 ) << "max-q : "<<M_max_q;
        M_model->addBasis( q );
        DVLOG( 2 ) << "basis q added";

        M_model->setMax(M_M, M_max_q, max_g, max_z, max_solution);
        ++M_M;
        M_model->clearOfflineError();
        timer2.restart();
        M_model->fillInterpolationMatrixFirstTime( );
        time=timer2.elapsed();
        time_=timer3.elapsed();
        if( this->worldComm().isMasterRank() )
        {
            std::cout<<" -- interpolation matrix filled in "<<time<<"s"<<std::endl;
            std::cout<<" -- time for this basis : "<<time_<<"s"<<std::endl;
        }

        M_greedy_rbmaxerr = 0;
        toc("EIM offline initialization (M=1) - " + M_model->name(), FLAGS_v>0);
    }//if M_restart
    else
    {
        M_M = M_model->mMax()+1;
        M_max_q = M_model->maxQ();
        max_g = M_model->maxG();
        max_z = M_model->maxZ();
        max_solution = M_model->maxSolution();

        int cobuild_freq = ioption(_prefix=this->M_prefixModel,_name="ser.eim-frequency");
        int use_rb = boption(_prefix=this->M_prefixModel,_name="ser.use-rb-in-eim-mu-selection") || boption(_prefix=this->M_prefixModel,_name="ser.use-rb-in-eim-basis-build");

        if( cobuild_freq != 0 && use_rb )
        {
            // Cobuild : If the first group (FEM solve) of EIM basis has already been built, go to loadDB (crb)
            if( M_model->mMax()-1 >= cobuild_freq && !M_model->rbBuilt() )
            {
                if( this->worldComm().isMasterRank() )
                    std::cout << "First group of EIM has already been built, start to load rb..." << std::endl;
                return;
            }
        }

    }//if ! M_restart
    /**
       \par build \f$W^g_M\f$
    */
    double err = 1;
    int cobuild_freq = ioption(_prefix=this->M_prefixModel,_name="ser.eim-frequency");
    LOG(INFO) << "[eim] cobuild frequency = " << cobuild_freq << "\n";
    int user_max = ioption(_prefix=this->M_prefix,_name="eim.dimension-max");
    int Mmax;
    int restart = M_restart ? 1 : 0;

    if( ioption(_prefix=this->M_prefixModel,_name="ser.eim-frequency") != 0 ) // SER
    {
        if( M_restart ) // First EIM basis with SER : do only initialization step (require only one FEM approx)
            Mmax = 0;
        else if ( M_model->adaptationSER() && M_model->mMax() < user_max ) // R-adaptation
            Mmax = M_model->mMax() + 1;
        else if ( M_model->mMax() + cobuild_freq  <= user_max )
            Mmax = M_model->mMax() + cobuild_freq;
        else
            Mmax = user_max;
    }
    else
        Mmax = user_max;

    //to deal with error estimation we need to build an "extra" basis function
    if( Mmax == user_max  )
        Mmax++;

    if( this->worldComm().isMasterRank() )
    {
        std::cout << "M_M = " << M_M << ", Mmax = " << Mmax << std::endl;
        std::cout << "RB correction = " << this->rbCorrection() << std::endl;
    }

    // Print maxerror (greedy) to file
    std::string eim_greedy_file_name = "cvg-eim-"+M_model->name()+"-Greedy-max-error.dat";
    std::string eim_greedy_inc_file_name = "cvg-eim-"+M_model->name()+"-Greedy-max-error-inc.dat";
    std::ofstream greedy_maxerr, greedy_maxerr_inc;
    if( this->worldComm().isMasterRank() )
    {
        greedy_maxerr.open(eim_greedy_file_name, std::ios::app);
        if( doption(_prefix=this->M_prefixModel,_name="ser.radapt-eim-rtol") != 0 )
            greedy_maxerr_inc.open(eim_greedy_inc_file_name, std::ios::app);
        if( M_M == 2 )
        {
            greedy_maxerr << "M" << "\t" << "maxErr" <<"\n";
            if( doption(_prefix=this->M_prefixModel,_name="ser.radapt-eim-rtol") != 0 )
                greedy_maxerr_inc << "M" << "\t" << "Increment" <<"\n";
        }
    }

    tic();
    LOG(INFO) << "start greedy algorithm...\n";
    for( ; M_M <=Mmax; ++M_M ) //err >= this->M_tol ) //Mmax == 1 : the basis has already been built at init step
    {
        if( !this->offlineStep() )
            break;

        timer3.restart();
        //LOG(INFO) << "M=" << M_M << "...\n";
        if( this->worldComm().isMasterRank() )
        {
                std::cout<<" ================================ "<<std::endl;
        }

        DVLOG(2) << "compute best fit error...\n";
        timer2.restart();
        // compute mu = arg max inf ||G(.;mu)-z||_infty
        tic();
        auto bestfit = computeBestFit( M_trainset, this->M_M-1 );
        toc("ComputeBestFit - M_M= " + std::to_string(M_M) + " - " + M_model->name(), FLAGS_v>0);

        // Print summary of EIM iterations
        if( boption(_prefix=this->M_prefixModel,_name="ser.use-rb-in-eim-mu-selection") && boption(_prefix=this->M_prefixModel,_name="ser.print-rb-iterations_info") )
            M_model->printRbIterationsSER( M_M-1 );

        time=timer2.elapsed();
        double error=bestfit.template get<0>();
        if( this->worldComm().isMasterRank() )
        {
            std::cout<<" -- best fit computed in "<<time<<"s -- absolute associated error : "<<error<<std::endl;
            greedy_maxerr << M_M << "\t" << error <<"\n";
        }

        M_model->addOfflineError(error);
        mu = bestfit.template get<1>();
        timer2.restart();

        if( ioption(_prefix=this->M_prefixModel,_name="ser.eim-frequency") != 0  ) // SER
        {
            // SER : choose between RB and PFEM approx
            tic();
            if( boption(_prefix=this->M_prefixModel,_name="ser.use-rb-in-eim-basis-build") && M_criterion )
                solution = M_model->computeRbExpansion( mu ); // Use current RB approx
            else
                solution = M_model->computePfem( mu ); // Use parametric FE with current affine decomposition
            toc("Compute solution (SER) - "+ M_model->name(), FLAGS_v>0);

            if( M_M > 2 ) // Ensure M_greedy_maxerr (error for previous EIM approx.) is initialized
            {
                // Need adapt group size (r-adaptation) ?
                double increment = math::abs( error - M_greedy_maxerr );
                double inc_relative = increment/math::abs( M_greedy_maxerr );
                if( this->worldComm().isMasterRank() && doption(_prefix=this->M_prefixModel,_name="ser.radapt-eim-rtol")!=0 )
                {
                    std::cout << " -- Absolute error (Greedy) relative increment = " << inc_relative
                              << ", rtol = " << doption(_prefix=this->M_prefixModel,_name="ser.radapt-eim-rtol") << std::endl;
                    greedy_maxerr_inc << M_M << "\t" << inc_relative <<"\n";
                }

                // Adapt only if user has given tol in option (default : 0)
                // increment > 1e-10 and inc_relative > 0 avoid to take into account EIM used for constants
                if( increment > 1e-10 && inc_relative > 0 && inc_relative < doption(_prefix=this->M_prefixModel,_name="ser.radapt-eim-rtol") )
                    this->setAdaptationSER( true );

                // Use corrected RB if increment is not sufficient
                this->setRbCorrection( false ); //Re-init to false
                if( increment > 1e-10 && inc_relative > 0 && inc_relative < doption(_prefix=this->M_prefixModel,_name="ser.corrected-rb-rtol") )
                {
                    if( this->worldComm().isMasterRank() )
                        std::cout << " -- Relative increment < tol : RB approx used for next EIM will be corrected " << std::endl;
                    this->setRbCorrection( true );
                }
            }

            // If criterion is satisfied,
            if( !adaptationSER() )
                M_greedy_maxerr = error;
        }
        else
        {
            tic();
            solution = M_model->solve( mu ); //No use of SER : use FE model since we don't have affine decomposition yet
            toc("Compute solution (No SER) - "+ M_model->name(), FLAGS_v>0);
        }

        time=timer2.elapsed();
        if( this->worldComm().isMasterRank() )
        {
            std::cout<<" -- model solution computed in "<<time<<"s"<<std::endl;
        }

        timer2.restart();
        time=timer2.elapsed();
        if( this->worldComm().isMasterRank() )
        {
            std::cout<<" -- maximum of expression computed in "<<time<<"s"<<std::endl;
        }

        //if we want to impose the use of dimension-max functions, we don't want to stop here
        //if ( (bestfit.template get<0>()/gmax.template get<0>()) < doption(_name="eim.error-max") &&  ! boption(_name="eim.use-dimension-max-functions") )
        if ( bestfit.template get<0>() < doption(_prefix=this->M_prefix,_name="eim.error-max") &&  ! boption(_prefix=this->M_prefix,_name="eim.use-dimension-max-functions") )
        {
            M_M--;
            this->setOfflineStep(false);
            break;
        }

        /**
         * we have a new \f$\mu^g_M\f$, insert it in \f$S^g_M\f$ and the
         * corresponding \f$z \in W^g_M\f$
         */
        DVLOG(2) << "[offline] S(" << this->M_M-1 << ") = " << bestfit.template get<1>() << "\n";

        if( ! expression_expansion )
        {
            timer2.restart();
            //M_model->addExpressionEvaluation( M_model->operator()( mu ) ); //projection de l'expression
            M_model->addExpressionEvaluation( M_model->operator()( solution, mu ) ); //projection de l'expression
            time=timer2.elapsed();
            if( this->worldComm().isMasterRank() )
            {
                std::cout<<" -- expression evaluated in mu in "<<time<<"s"<<std::endl;
            }
            max_g++;
        }
        else
        {
            M_model->addSolution( solution );
            max_solution++;
        }

        // build T^m such that T^m-1 \subset T^m
        DVLOG(2) << "[offline] compute residual M="<< M_M << "..." <<"\n";
        //res = this->residual(M_M-1);

        //LOG(INFO) << "residual = " << res << "\n";
        LOG(INFO) << "[offline] compute arg sup |residual|..." <<"\n";
        //auto resmax = normLinf( _range=elements(M_model->mesh()), _pset=_Q<5>(), _expr=idv(res) );
        auto coeff = M_model->computeExpansionCoefficients( mu ,  solution , M_M-1 );
        auto z = expansion( M_model->q(), coeff , M_M-1 );

        if( expression_expansion )
        {
            M_model->addZ( z );
            max_z++;
        }

        timer2.restart();
        tic();
        auto resmax = M_model->computeMaximumOfResidual( mu, solution , z );
        toc("computeMaximumOfResidual - "+ M_model->name(), FLAGS_v>0);
        time=timer2.elapsed();
        if( this->worldComm().isMasterRank() )
        {
            std::cout<<" -- Maximum of residual computed in "<<time<<"s"<<std::endl;
        }

        t = resmax.template get<1>();

        // store space coordinate where max absolute value occurs
        DVLOG(2) << "[offline] store coordinates where max absolute value is attained : \n" << resmax.template get<1>();

        M_model->addParameter( mu );
        M_model->addInterpolationPoint( t ); // Store interpolation point
        element_type res = M_model->Residual( M_M-1, z );

        //DVLOG(2) << "[offline] scale new basis function by " << 1./resmax.template get<0>() << "..." <<"\n";
        //std::cout << "[offline] eval (scaling) = " << res( t )(0,0,0) << std::endl;
        //res.scale( 1./res( t )(0,0,0) );
        DVLOG(2) << "store new basis function..." <<"\n";

        //if we want to impose the use of dimension-max functions, we don't want to stop here
        if ( resmax.template get<0>() < doption(_prefix=this->M_prefix,_name="eim.error-max") &&  ! boption(_prefix=this->M_prefix,_name="eim.use-dimension-max-functions") )
        {
            M_M--;
            this->setOfflineStep(false);
            break;
        }

        //M_model->addInterpolationPoint( t ); // Store interpolation point
        M_model->addBasis( res );
        M_max_q++;
        M_model->setMax(M_M, M_max_q, max_g, max_z, max_solution);
        timer2.restart();
        M_model->fillInterpolationMatrix( );
        time=timer2.elapsed();
        time_=timer3.elapsed();
        if( this->worldComm().isMasterRank() )
        {
            std::cout<<" -- interpolation matrix filled in "<<time<<"s"<<std::endl;
            std::cout<<" -- time for this basis : "<<time_<<"s"<<std::endl;
            std::cout<<" M_M : "<<M_M<<std::endl;
        }

        VLOG(2) << "================================================================================\n";

    }
    toc("EIM offline (M > 2) - " + M_model->name(), FLAGS_v>0);

    if( this->worldComm().isMasterRank() )
    {
        greedy_maxerr.close();
        if( doption(_prefix=this->M_prefixModel,_name="ser.radapt-eim-rtol") != 0 )
            greedy_maxerr_inc.close();
    }

    time=timer.elapsed();
    if( this->worldComm().isMasterRank() )
        std::cout<<"Total time for offline step of EIM "<<M_model->name()<<" : "<<time<<"s\n"<<std::endl;
    DVLOG(2) << "[offline] M_max = " << M_M << "...\n";

    this->M_offline_done = true;
    if( Mmax >= user_max ) //dimension-max has been reached
        this->setOfflineStep(false);
}

template<typename ModelType>
void
EIM<ModelType>::studyConvergence( parameter_type const & mu , model_solution_type & solution , std::vector<std::string> all_file_name ) const
{
    LOG(INFO) << " Convergence study \n";
    int proc_number =  this->worldComm().globalRank();

    std::vector<double> l2ErrorVec(M_model->mMax(), 0.0);

    //here are files containing l2error and others quantities for all runs performed
    int number_of_files=all_file_name.size();
    std::ofstream fileL2 ( all_file_name[0] ,std::ios::out | std::ios::app );
    std::ofstream fileL2estimated (  all_file_name[1] , std::ofstream::out | std::ofstream::app );
    std::ofstream fileL2ratio (  all_file_name[2] , std::ofstream::out | std::ofstream::app );
    std::ofstream fileLINF (  all_file_name[3] , std::ofstream::out | std::ofstream::app );
    std::ofstream fileLINFestimated (  all_file_name[4] , std::ofstream::out | std::ofstream::app );
    std::ofstream fileLINFratio (  all_file_name[5] , std::ofstream::out | std::ofstream::app );
    std::string head = "#Nbasis \t  L2_error \t L2_estimated \t ratio_l2 \t linf_error \t linf_estimated \t ratio_linf \n";

    int max = M_model->mMax();
    int Nmax=0;

    double relative_l2_error;
    double relative_linf_error;
    double absolute_l2_error;
    double interpolation_error;
    double absolute_linf_error_estimated ;
    double absolute_l2_error_estimated;
    double relative_l2_error_estimated;
    double relative_linf_error_estimated;
    double relative_ratio_l2;
    double relative_ratio_linf;
    //double absolute_ratio_linf;

    //As we print error estimation, we stop at max-1
    //because we need to access to the max^th basis function
    if( this->worldComm().isMasterRank() )
    {
        Nmax = max;
        fileL2 << Nmax<< "\t";
        fileL2estimated << Nmax-1 <<"\t";
        fileL2ratio << Nmax-1  <<"\t" ;
        fileLINF << Nmax  <<"\t" ;
        fileLINFestimated << Nmax-1  <<"\t" ;
        fileLINFratio << Nmax-1  <<"\t" ;
    }
    for(int N=1; N<=max; N++)
    {
        double exprl2norm = 0 , diffl2norm = 0 ;
        exprl2norm =M_model->projExpressionL2Norm( solution , mu );
        auto eim_approximation = this->operator()(mu , solution, N);
        diffl2norm = M_model->projDiffL2Norm( solution , mu , eim_approximation );
        double absolute_linf_error = M_model->projDiffLinfNorm( solution , mu , eim_approximation );
        double exprlinfnorm = M_model->projExprLinfNorm( solution , mu );

        if( N < max )
        {
            relative_l2_error = diffl2norm / exprl2norm ;
            relative_linf_error = absolute_linf_error / exprlinfnorm ;
            absolute_l2_error = diffl2norm ;
            //interpolation error : || projection_g - g ||_L2
            interpolation_error = M_model->interpolationError( solution , mu );
            absolute_linf_error_estimated = this->errorEstimationLinf( mu , solution, N );
            relative_linf_error_estimated = absolute_linf_error_estimated/exprlinfnorm;
            auto tuple = this->interpolationErrorEstimation( mu , solution, N );
            auto error_estimation_element = tuple.template get<1>();
            absolute_l2_error_estimated = error_estimation_element.l2Norm();
            relative_l2_error_estimated = absolute_l2_error_estimated/exprl2norm;
            relative_ratio_l2 = math::abs( relative_l2_error_estimated / relative_l2_error );
            //absolute_ratio_linf = math::abs( absolute_linf_error_estimated / absolute_linf_error );
            relative_ratio_linf = math::abs( relative_linf_error_estimated / relative_linf_error );
        }

        std::string str = "\t";
        if( N == Nmax ) str = "\n";

        if( this->worldComm().isMasterRank() )
        {
            fileL2            << relative_l2_error            <<str;
            fileLINF          << relative_linf_error          <<str;

            if( N == Nmax-1 )
                str = "\n";
            else
                str= "\t";
            if( N < max )
            {
                fileL2estimated   << relative_l2_error_estimated  <<str;
                fileL2ratio       << relative_ratio_l2            <<str;
                fileLINFestimated << relative_linf_error_estimated<<str;
                if( relative_linf_error < 1e-14 )
                    fileLINFratio     << 0          <<str;
                else
                    fileLINFratio     << relative_ratio_linf          <<str;
            }
        }

    }//loop over basis functions

    fileL2.close();
    fileL2estimated.close();
    fileL2ratio.close();
    fileLINF.close();
    fileLINFestimated.close();
    fileLINFratio.close();

#if 0
    bool use_expression = boption(_prefix=this->M_prefix,_name="eim.compute-error-with-truth-expression");

        if( use_expression )
        {
            exprl2norm =M_model->expressionL2Norm( solution , mu );
            auto eim_approximation = this->operator()(mu , solution, N);
            diffl2norm = M_model->diffL2Norm( solution , mu , eim_approximation );
        }

        int size = mu.size();
        LOG(INFO)<<" mu = [ ";
        for ( int i=0; i<size-1; i++ ) LOG(INFO)<< mu[i] <<" , ";
        LOG(INFO)<< mu[size-1]<<" ]\n";

        //old version
        // Compute expression
        auto expression = M_model->operator()(mu);
        // Compute eim expansion
        auto eim_approximation = this->operator()(mu , N);
        //Compute l2error
        LOG(INFO) << "EIM name : " << M_model->name() << "\n";
        double norm_l2_expression = expression.l2Norm();
        LOG(INFO) << "norm_l2 expression = " << norm_l2_expression << "\n";
        double norm_l2_approximation = eim_approximation.l2Norm();
        LOG(INFO) << "norm_l2_approximation = " << norm_l2_approximation << "\n";
        auto l2_error = math::abs( norm_l2_expression - norm_l2_approximation ) / norm_l2_expression;
        LOG(INFO) << "norm l2 error = " << l2_error << "\n";

    std::string mu_str;
    for ( int i=0; i<mu.size(); i++ )
        mu_str= mu_str + ( boost::format( "_%1%" ) %mu(i) ).str() ;

    std::string file_name = "cvg-eim-"+M_model->name()+"-"+mu_str+".dat";
    if( std::ifstream( file_name ) )
        std::remove( file_name.c_str() );

    std::ofstream conv;
    if( proc_number == this->worldComm().masterRank() )
        {
            conv.open(file_name, std::ios::app);
            conv << "#Nbasis" << "\t" << "L2_error \t L2_estimated \t ratio_l2 \t linf_error \t linf_estimated \t ratio_linf \t interpolation_error " <<"\n";
        }

#endif

}

template<typename SpaceType, typename ModelSpaceType, typename ParameterSpaceType>
class EIMFunctionBase : public CRBDB

{
    typedef CRBDB super_type;
public:

    typedef SpaceType functionspace_type;
    typedef std::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef typename functionspace_type::element_ptrtype element_ptrtype;
    typedef typename functionspace_type::element_type solution_type;
    typedef typename functionspace_type::element_type element_type;

    typedef typename functionspace_type::mesh_type mesh_type;
    typedef typename functionspace_type::mesh_ptrtype mesh_ptrtype;
    typedef typename functionspace_type::value_type value_type;
    typedef typename functionspace_type::Context context_type;

    typedef ModelSpaceType model_functionspace_type;
    typedef std::shared_ptr<model_functionspace_type> model_functionspace_ptrtype;
    typedef typename model_functionspace_type::element_type model_element_type;
    typedef typename model_functionspace_type::element_type model_solution_type;
    typedef typename model_functionspace_type::element_ptrtype model_element_ptrtype;
    typedef typename model_functionspace_type::Context model_solution_context_type;

    static const uint16_type nDim = mesh_type::nDim;

    typedef ParameterSpaceType parameterspace_type;
    typedef std::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef typename parameterspace_type::element_type parameter_type;
    typedef Eigen::Matrix<double, nDim, 1> node_type;
    //typedef Eigen::Matrix<double, SpaceType::basis_type::nLocalDof, Eigen::Dynamic> matrix_basis_pc_type;

    typedef Eigen::VectorXd vectorN_type;
    typedef Eigen::MatrixXd matrixN_type;

    typedef typename parameterspace_type::sampling_ptrtype sampling_ptrtype;
    typedef typename parameterspace_type::sampling_type sampling_type;

    typedef EIM<EIMFunctionBase<SpaceType,model_functionspace_type, ParameterSpaceType> > eim_type;
    typedef typename eim_type::vector_type vector_type;

    typedef std::shared_ptr<eim_type> eim_ptrtype;
    //typedef typename eim_type::betam_type betam_type;
    //typedef typename eim_type::qm_type qm_type;

    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_type;
    typedef std::shared_ptr<matrix_type> matrix_ptrtype;
    typedef boost::tuple<double,parameter_type> parameter_residual_type;

    enum ResidualNormType { Linfty = 0, L2 = 1, LinftyVec = 2 };

    typedef GeometricSpace<mesh_type> geometricspace_type;
    typedef std::shared_ptr<geometricspace_type> geometricspace_ptrtype;
    typedef typename geometricspace_type::Context geometricspace_context_type;
    typedef std::shared_ptr<geometricspace_context_type> geometricspace_context_ptrtype;

    typedef CRBBase<model_functionspace_type,parameterspace_type> crb_type;
    typedef std::shared_ptr<crb_type> crb_ptrtype;

    EIMFunctionBase( parameterspace_ptrtype const& pspace,
                     sampling_ptrtype const& sampling,
                     std::string const& modelname,
                     std::string const& name,
                     uuids::uuid const& uid,
                     std::string const& prefix)
        :
        super_type( name, "eim", uid ),
        M_prefix( prefix ),
        M_fspace(),
        M_pspace( pspace ),
        M_trainset( sampling ),
        M_modelname( modelname ),
        M_name( name ),
        M_computeExpansionOfExpression( boption(_prefix=this->M_prefix,_name="eim.compute-expansion-of-expression") ),
        M_normUsedForResidual( ResidualNormType::Linfty)
        {
            this->setDBDirectory( modelname,uid );
            std::string norm_used = soption(_prefix=this->M_prefix,_name="eim.norm-used-for-residual");
            if( norm_used == "Linfty" )
                M_normUsedForResidual = ResidualNormType::Linfty;
            else if( norm_used == "L2" )
                M_normUsedForResidual = ResidualNormType::L2;
            else if( norm_used == "LinftyVec" )
                M_normUsedForResidual = ResidualNormType::LinftyVec;
            else
                CHECK( false ) <<"[EIM options] The name of the norm "<<norm_used<<" is not known\n";

            LOG(INFO)<< "EimFunctionBase constructor\n";
        }
    EIMFunctionBase( functionspace_ptrtype const& fspace,
                     parameterspace_ptrtype const& pspace,
                     sampling_ptrtype const& sampling,
                     std::string const& modelname,
                     std::string const& name,
                     uuids::uuid const& uid,
                     std::string const& prefix )
        :
        EIMFunctionBase( pspace, sampling, modelname, name, uid, prefix )
        {
            M_fspace = fspace;
        }
    ~EIMFunctionBase() override
        {}
    std::string const& modelName() const { return M_modelname; }
    void setModelName( std::string const& name ) { M_modelname = name; }

    functionspace_ptrtype const& functionSpace() const { return M_fspace; }
    functionspace_ptrtype functionSpace() { return M_fspace; }
    parameterspace_ptrtype const& parameterSpace() const { return M_pspace; }
    parameterspace_ptrtype parameterSpace()  { return M_pspace; }
    sampling_ptrtype trainSet() { return M_trainset; }
    sampling_ptrtype const& trainSet() const { return M_trainset; }
    virtual void setTrainSet( sampling_ptrtype tset ) { M_trainset = tset; }

    bool computeExpansionOfExpression() const { return M_computeExpansionOfExpression; }

    void addOnlineTime( const double time )
    {
        int size = M_online_time.size();
        M_online_time.conservativeResize( size+1 );
        M_online_time( size ) = time;
    }
    Eigen::VectorXd onlineTime() const { return M_online_time; }

    mesh_ptrtype const& mesh() const { return M_fspace->mesh(); }
    mesh_ptrtype mesh()  { return M_fspace->mesh(); }

    virtual bool modelUseSolve() const = 0;

    virtual element_type operator()( parameter_type const& ) = 0;
    virtual element_type operator()( model_solution_type const& , parameter_type const& ) = 0;
    virtual element_type interpolant( parameter_type const& ) = 0;
    virtual element_type interpolant( parameter_type const& , model_solution_type const & , int ) = 0;

    value_type operator()( node_type const& x, parameter_type const& mu )
        {
            VLOG(2) << "calling EIMFunctionBase::operator()( x=" << x << ", mu=" << mu << ")\n";
            element_type v = this->operator()( mu );
            value_type res = v(x)(0,0,0);
            VLOG(2) << "EIMFunctionBase::operator() v(x)=" << res << "\n";
            return res;
        }

    // evaluate eim expansion at interpolation points in space and mu in parameter where T provides the coefficient
    //value_type operator()( vector_type const& T, parameter_type const& mu ) = 0;

    value_type operator()( model_solution_type const& T, node_type const& x, parameter_type const& mu )
        {
            VLOG(2) << "calling EIMFunctionBase::operator()( x=" << x << ", mu=" << mu << ")\n";
            element_type v = this->operator()( T, mu );
            value_type res = v(x)(0,0,0);
            VLOG(2) << "EIMFunctionBase::operator() v(x)=" << res << "\n";
            return res;
        }

#if 0
    virtual vector_type operator()( model_solution_context_type/*context_type*/ const& ctx, parameter_type const& mu , int M=-1) = 0;
    virtual vector_type operator()( model_solution_type const& T, model_solution_context_type/*context_type*/ const& ctx, parameter_type const& mu , int M=-1) = 0;
#endif
    virtual element_type const& q( int m )  const = 0;
    virtual std::vector<element_type> const& q() const = 0;
    virtual vector_type  beta( parameter_type const& mu ) /*const*/ = 0;
    virtual vector_type  beta( parameter_type const& mu, model_solution_type const& T ) /*const*/ = 0;
    virtual vector_type  beta( parameter_type const& mu, vectorN_type const& T ) /*const*/ = 0;
    virtual vector_type  beta( parameter_type const& mu , size_type M )  = 0;
    virtual vector_type  beta( parameter_type const& mu, model_solution_type const& T , size_type M)  = 0;
    virtual vector_type  beta( parameter_type const& mu, vectorN_type const& T , size_type M)  = 0;

    virtual size_type  mMax(bool & error) const = 0;
    virtual size_type  mMax() const = 0;

    virtual boost::tuple<double,element_type> interpolationErrorEstimation ( parameter_type const& mu, model_solution_type const& solution , int M) const = 0;
    virtual double errorEstimationLinf( parameter_type const & mu, model_solution_type const& solution , int M ) const=0 ;
    virtual node_type interpolationPoint( int position ) const = 0;

    virtual void studyConvergence( parameter_type const & mu , model_solution_type & solution, std::vector< std::string > all_file_name ) const = 0;
    virtual void computationalTimeStatistics( std::string appname )  = 0;
    virtual double expressionL2Norm( model_solution_type const& T , parameter_type const& mu ) const = 0;
    virtual double diffL2Norm( model_solution_type const& T , parameter_type const& mu , element_type const& eim_expansion ) const = 0;
    virtual double projExpressionL2Norm( model_solution_type const& T , parameter_type const& mu ) const = 0;
    virtual double projDiffL2Norm( model_solution_type const& T , parameter_type const& mu , element_type const& eim_expansion ) const = 0;
    virtual double projDiffLinfNorm( model_solution_type const& T , parameter_type const& mu , element_type const& eim_expansion ) const = 0;
    virtual double projExprLinfNorm( model_solution_type const& T , parameter_type const& mu ) const = 0;
    virtual double interpolationError( model_solution_type const& T , parameter_type const& mu ) const = 0;
    virtual model_solution_type solve( parameter_type const&  mu ) = 0;
    virtual element_type residual( element_type const& z, element_type const& g ) = 0;
    virtual vector_type computeExpansionCoefficients( parameter_type const& mu, model_solution_type const& solution,  int M)=0;
    virtual boost::tuple<double,node_type> computeMaximumOfExpression( parameter_type const& mu, model_solution_type const& solution )=0;
    virtual boost::tuple<double,node_type> computeMaximumOfResidual( parameter_type const& mu, model_solution_type const& solution  , element_type const& z)=0;
    virtual void fillInterpolationMatrix( )=0;
    virtual void fillInterpolationMatrixFirstTime( )=0;

    virtual model_functionspace_ptrtype modelFunctionSpace() const = 0;
    virtual model_functionspace_ptrtype modelFunctionSpace() = 0;
    virtual void addInterpolationPoint( node_type const& t ) = 0;

    virtual void updateRbSpaceContext( boost::any const& rbspacebase ) = 0;
    virtual bool hasRbSpaceContext() const = 0;
    virtual bool hasRbSpaceContext2() const = 0;
    virtual boost::any rbSpaceContext() const = 0;
    virtual boost::any rbSpaceContext2() const = 0;
    virtual void setRbSpaceContext( boost::any const& rbCtxBase ) = 0;
    virtual void setRbSpaceContext2( boost::any const& rbCtxBase ) = 0;

    virtual void setMax(int m, int max_q, int max_g, int max_z, int max_solution) = 0;
    virtual int maxQ() const = 0;
    virtual int maxG() const = 0;
    virtual int maxZ() const = 0;
    virtual int maxSolution() const = 0;

    virtual bool offlineStep() const = 0;
    virtual bool adaptationSER() const = 0;
    virtual void setAdaptationSER(bool b)=0;
    virtual bool rbCorrection() const = 0;
    virtual void setRbCorrection(bool b)=0;
    virtual void setRestart(bool b)=0;
    virtual void setRB( crb_ptrtype rb )=0;
    virtual bool rbBuilt() const = 0;
    virtual void offline()=0;

    virtual void addBasis( element_type const &q ) = 0;
    virtual void addExpressionEvaluation( element_type const &g ) = 0;
    virtual void addZ( element_type const &z ) = 0;
    virtual void addSolution( model_solution_type const &solutions ) = 0;
    virtual void addParameter( parameter_type const& mu ) = 0;
    virtual void clearParameterSampling() = 0;
    virtual boost::any buildBasisFunction( int m ) = 0;
    virtual element_type Residual( int m, element_type const& z) = 0;

    virtual void printInterpolationPointsSelection() const=0;
    virtual void printMuSelection() const=0;
    virtual void printOfflineError() const=0;
    virtual void printRbIterationsSER( int M ) const=0;
    virtual void addOfflineError(double error) =0;
    virtual void clearOfflineError() =0;

    virtual void initializeDataStructures() = 0;

    virtual vector_type evaluateExpressionAtInterpolationPoints(model_solution_type const &solution, parameter_type const& mu, int M)=0;
    virtual vector_type evaluateElementAtInterpolationPoints(element_type const & element, int M)=0;
    virtual model_solution_type computeRbExpansion( parameter_type const& mu)=0;
    virtual model_solution_type computeRbExpansion( parameter_type const& mu, std::vector<vectorN_type> uN)=0;
    virtual model_solution_type computePfem( parameter_type const& mu )=0;
    virtual boost::tuple<double,std::vector<vectorN_type> > RieszResidualNorm( parameter_type const& mu )=0;
    virtual sampling_ptrtype createSubTrainset( sampling_ptrtype const& trainset, int method )=0;

protected :
    std::string M_prefix;
    functionspace_ptrtype M_fspace;
    parameterspace_ptrtype M_pspace;
    sampling_ptrtype M_trainset;
    std::string M_modelname;
    std::string M_name;
    Eigen::VectorXd M_online_time;//contains online computational time

    bool M_computeExpansionOfExpression;
    ResidualNormType M_normUsedForResidual;
};

template<typename EimSpaceType,typename ModelType>
class EIMFunctionFeSpaceDb
{
public :
    typedef EimSpaceType eim_functionspace_type;
    typedef std::shared_ptr<eim_functionspace_type> eim_functionspace_ptrtype;
    typedef typename eim_functionspace_type::element_type eim_element_type;

    typedef ModelType model_type;
    typedef std::shared_ptr<model_type> model_ptrtype;
    typedef std::weak_ptr<model_type> model_weakptrtype;

    typedef typename model_type::functionspace_type model_functionspace_type;
    typedef std::shared_ptr<model_functionspace_type> model_functionspace_ptrtype;
    typedef typename model_functionspace_type::element_type model_element_type;

    EIMFunctionFeSpaceDb()
        :
        M_max_q( 0 ),
        M_max_g( 0 ),
        M_max_z( 0 ),
        M_max_solution( 0 )
        {}

    eim_functionspace_ptrtype const& eimFunctionSpace() { return M_eimFunctionSpace; }
    model_functionspace_ptrtype const modelFunctionSpace() const { return M_model.lock()->functionSpace(); }

    int maxQ() const { return M_max_q; }
    int maxG() const { return M_max_g; }
    int maxZ() const { return M_max_z; }
    int maxSolution() const { return M_max_solution; }

    std::vector<eim_element_type> const& q() const { return M_q_vector; }
    eim_element_type const& q( int k ) const { return M_q_vector[k]; }
    std::vector<eim_element_type> const& g() const { return M_g_vector; }
    eim_element_type const& g( int k ) const { return M_g_vector[k]; }
    std::vector<eim_element_type> const& z() const { return M_z_vector; }
    eim_element_type const& z( int k ) const { return M_z_vector[k]; }
    std::vector<model_element_type> const& solutions() const { return M_solution_vector; }
    model_element_type const& solution( int k ) const { return M_solution_vector[k]; }

    void setEimFunctionSpace( eim_functionspace_ptrtype const& space ) { M_eimFunctionSpace = space; }
    void setModel( model_ptrtype const& model ) { M_model = model; }

    void setMax( int max_q, int max_g, int max_z, int max_solution )
    {
        M_max_q = max_q;
        M_max_g = max_g;
        M_max_z = max_z;
        M_max_solution = max_solution;
    }

    void addBasis( eim_element_type const &q )
    {
        M_q_vector.push_back( q );
    }
    void addExpressionEvaluation( eim_element_type const &g )
    {
        M_g_vector.push_back( g );
    }
    void addZ( eim_element_type const &z )
    {
        M_z_vector.push_back( z );
    }
    void addSolution( model_element_type const &solution )
    {
        M_solution_vector.push_back( solution );
    }

    void clear()
    {
        M_q_vector.clear();
        M_g_vector.clear();
        M_z_vector.clear();
        M_solution_vector.clear();
        this->setMax( 0,0,0,0 );
    }

    bool load( std::string const& dbName, std::string const& dbDir, WorldComm const& worldComm )
    {
        if  ( !this->eimFunctionSpace() )
        {
            if ( !M_model.lock() )
                return false;
            else if ( !this->modelFunctionSpace() )
                return false;
        }

        std::string feDbFilename = ( boost::format( "%1%_p%2%.crbdb" ) %dbName %worldComm.globalRank() ).str();
        fs::ifstream ifs( fs::path(dbDir) / feDbFilename );
        if ( ifs )
        {
            //boost::archive::text_iarchive ia( ifs );
            boost::archive::binary_iarchive ia( ifs );
            //write class instance to archive
            ia >> *this;
            return true;
        }
        return false;
    }

    void save( std::string const& dbName, std::string const& dbDir, WorldComm const& worldComm )
    {
        std::string feDbFilename = ( boost::format( "%1%_p%2%.crbdb" ) %dbName %worldComm.globalRank() ).str();
        fs::ofstream ofs( fs::path(dbDir) / feDbFilename );
        if ( ofs )
        {
            boost::archive::binary_oarchive oa( ofs );
            oa << *this;
        }
    }


    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & __ar, const unsigned int __version )
    {
        //max of eim basis function projected on function space
        __ar & BOOST_SERIALIZATION_NVP( M_max_q );
        DVLOG(2) << "max-q saved/loaded";
        __ar & BOOST_SERIALIZATION_NVP( M_max_g );
        DVLOG(2) << "max-g saved/loaded";
        __ar & BOOST_SERIALIZATION_NVP( M_max_z );
        DVLOG(2) << "max-z saved/loaded\n";
        __ar & BOOST_SERIALIZATION_NVP( M_max_solution );
        DVLOG(2) << "max-solution saved/loaded\n";

        if ( Archive::is_loading::value )
        {
            if( M_q_vector.size() == 0 && this->eimFunctionSpace() )
            {
                for(int i = 0; i < M_max_q; i++ )
                    M_q_vector.push_back( this->eimFunctionSpace()->element() );
                for( int i = 0; i < M_max_q; ++ i )
                    __ar & BOOST_SERIALIZATION_NVP( M_q_vector[i] );
            }
            DVLOG(2) << "q saved/loaded ( "<<M_max_g<<" elements ) ";

            if( M_g_vector.size() == 0 && this->eimFunctionSpace() )
            {
                for(int i = 0; i < M_max_g; i++ )
                    M_g_vector.push_back( this->eimFunctionSpace()->element() );
                for( int i = 0; i < M_max_g; ++ i )
                    __ar & BOOST_SERIALIZATION_NVP( M_g_vector[i] );
            }
            DVLOG(2) << "g saved/loaded ( "<<M_max_g<<" elements ) ";

            if( M_z_vector.size() == 0 && this->eimFunctionSpace() )
            {
                for(int i = 0; i < M_max_z; i++ )
                    M_z_vector.push_back( this->eimFunctionSpace()->element() );
                for( int i = 0; i < M_max_z; ++ i )
                    __ar & BOOST_SERIALIZATION_NVP( M_z_vector[i] );
            }
            DVLOG(2) << "z saved/loaded ( "<<M_max_z<<" elements ) ";

            if( M_solution_vector.size() == 0 && this->modelFunctionSpace() )
            {
                for(int i = 0; i < M_max_solution; i++ )
                    M_solution_vector.push_back( this->modelFunctionSpace()->element() );
                for( int i = 0; i < M_max_solution; ++ i )
                    __ar & BOOST_SERIALIZATION_NVP( M_solution_vector[i] );
            }
            DVLOG(2) << "vector of solutions saved/loaded ( "<<M_max_solution<<" elements ) ";
        }
        else
        {
            for( int i = 0; i < M_max_q; ++ i )
                __ar & BOOST_SERIALIZATION_NVP( M_q_vector[i] );
            for( int i = 0; i < M_max_g; ++ i )
                __ar & BOOST_SERIALIZATION_NVP( M_g_vector[i] );
            for( int i = 0; i < M_max_z; ++ i )
                __ar & BOOST_SERIALIZATION_NVP( M_z_vector[i] );
            for( int i = 0; i < M_max_solution; ++ i )
                __ar & BOOST_SERIALIZATION_NVP( M_solution_vector[i] );
        }
    }

private :
    eim_functionspace_ptrtype M_eimFunctionSpace;
    model_functionspace_ptrtype M_modelFunctionSpace;
    model_weakptrtype M_model;

    int M_max_q;
    int M_max_g;
    int M_max_z;
    int M_max_solution;

    std::vector< eim_element_type > M_q_vector ;
    std::vector< eim_element_type > M_g_vector ;
    std::vector< eim_element_type > M_z_vector ; // we need to store this to be able to compute expression of eim basis functions and the user can use them in integrals
    std::vector< model_element_type > M_solution_vector ;
};


template<typename ModelType, typename SpaceType, typename ExprType, int SubSpaceId = -1, int SubSpaceId2 = -1>
class EIMFunction
    : public EIMFunctionBase<SpaceType, typename ModelType::functionspace_type, typename ModelType::parameterspace_type>
{
    typedef EIMFunctionBase<SpaceType, typename ModelType::functionspace_type, typename ModelType::parameterspace_type> super;
public:
    typedef ModelType model_type;
    //typedef ModelType* model_ptrtype;
    typedef std::shared_ptr<model_type> model_ptrtype;
    typedef std::weak_ptr<model_type> model_weakptrtype;

    static const bool model_use_nosolve = boost::is_base_of<EimFunctionNoSolveBase,model_type>::type::value;
    static const bool model_use_solve = !model_use_nosolve;
    //static const bool model_use_solve = boost::is_base_of<ModelCrbBaseBase,model_type>::type::value;
    static const bool model_is_modelcrbbase = boost::is_base_of<ModelCrbBaseBase,model_type>::type::value;

    typedef SpaceType functionspace_type;
    typedef std::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef typename super::element_type element_type;
    typedef typename super::element_ptrtype element_ptrtype;
    typedef typename super::model_element_type model_element_type;
    typedef typename super::model_element_ptrtype model_element_ptrtype;
    typedef typename super::value_type value_type;

    typedef typename super::model_solution_type model_solution_type;
    typedef typename super::solution_type solution_type;

    typedef typename super::model_functionspace_type model_functionspace_type;
    typedef std::shared_ptr<model_functionspace_type> model_functionspace_ptrtype;
    typedef typename super::model_solution_context_type model_solution_context_type;
    static const bool use_subspace_element = (SubSpaceId != -1 && model_functionspace_type::is_composite);

    struct fake_element_composite_type
    {
        template<int i>
        struct sub_element
        {
            typedef typename functionspace_type::element_type type;
        };
    };
    typedef typename mpl::if_< mpl::bool_<model_functionspace_type::is_composite>,
                               typename model_functionspace_type::element_type,
                               fake_element_composite_type >::type model_element_composite_type;
    typedef typename mpl::if_< mpl::bool_<use_subspace_element>,
                               mpl::identity<typename model_element_composite_type::template sub_element<SubSpaceId>::type>,
                               mpl::identity<typename model_functionspace_type::element_type> >::type::type model_element_expr_type;
    typedef typename mpl::if_< mpl::bool_<use_subspace_element>,
                               mpl::identity<typename model_element_composite_type::template sub_element<SubSpaceId2>::type>,
                               mpl::identity<typename model_functionspace_type::element_type> >::type::type model_element2_expr_type;

    typedef typename model_element_expr_type::functionspace_type  model_element_expr_functionspace_type;
    typedef typename model_element_expr_functionspace_type::Context model_element_expr_context_type;
    typedef std::shared_ptr<model_element_expr_context_type> model_element_expr_context_ptrtype;

    typedef typename SpaceType::mesh_type mesh_type;
    static const uint16_type nDim = mesh_type::nDim;

    typedef typename super::parameterspace_type parameterspace_type;
    typedef typename super::parameter_type parameter_type;

    typedef typename parameterspace_type::sampling_ptrtype sampling_ptrtype;
    typedef typename parameterspace_type::sampling_type sampling_type;

    typedef typename super::geometricspace_type geometricspace_type;
    typedef typename super::geometricspace_ptrtype geometricspace_ptrtype;
    typedef typename super::geometricspace_context_type geometricspace_context_type;
    typedef typename super::geometricspace_context_ptrtype geometricspace_context_ptrtype;

    // reduced basis space
    //typedef ReducedBasisSpace<model_type> rbfunctionspace_type;
    // typedef typename model_type::rbfunctionspace_type rbfunctionspace_type;
    template <typename TheModelType,typename TheModelHasRbSpaceType>
    struct RbSpaceFromModel
    {
        typedef ReducedBasisSpace<typename TheModelType::functionspace_type> type;
    };
    template <typename TheModelType>
    struct RbSpaceFromModel<TheModelType,mpl::bool_<true> >
    {
        typedef typename TheModelType::rbfunctionspace_type type;
    };
    typedef typename RbSpaceFromModel<model_type,mpl::bool_<model_is_modelcrbbase> >::type rbfunctionspace_type;
    typedef std::shared_ptr<rbfunctionspace_type> rbfunctionspace_ptrtype;

    template <typename TheRbSpaceType,int TheSubSpaceId>
    struct RbSpaceCtxFromRbSpace
    {
        typedef typename TheRbSpaceType::template sub_rbfunctionspace<TheSubSpaceId>::type::ctxrbset_type type;
    };
    template <typename TheRbSpaceType>
    struct RbSpaceCtxFromRbSpace<TheRbSpaceType,-1>
    {
        typedef typename TheRbSpaceType::ctxrbset_type type;
    };
    typedef typename RbSpaceCtxFromRbSpace<rbfunctionspace_type,SubSpaceId>::type rbfunctionspace_context_type;
    typedef std::shared_ptr<rbfunctionspace_context_type> rbfunctionspace_context_ptrtype;
    typedef typename RbSpaceCtxFromRbSpace<rbfunctionspace_type,SubSpaceId2>::type rbfunctionspace_context2_type;
    typedef std::shared_ptr<rbfunctionspace_context2_type> rbfunctionspace_context2_ptrtype;

    typedef ExprType expr_type;
    typedef std::shared_ptr<expr_type> expr_ptrtype;

    typedef typename super::eim_type eim_type;
    typedef typename super::eim_ptrtype eim_ptrtype;
    typedef typename super::context_type context_type;
    typedef typename super::vector_type vector_type;

    typedef CRBModel<ModelType> crbmodel_type;
    typedef std::shared_ptr<crbmodel_type> crbmodel_ptrtype;
    typedef CRBBase<typename model_type::functionspace_type,parameterspace_type> crb_type;
    typedef std::shared_ptr<crb_type> crb_ptrtype;

    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_type;
    typedef std::shared_ptr<matrix_type> matrix_ptrtype;

    typedef boost::tuple<double,parameter_type> parameter_residual_type;
    typedef typename super::node_type node_type;

    typedef Eigen::VectorXd vectorN_type;

    EIMFunction( model_ptrtype model,
                 functionspace_ptrtype space,
                 model_element_expr_type &u,
                 model_element2_expr_type &u2,
                 parameter_type& mu,
                 expr_type& expr,
                 sampling_ptrtype sampling,
                 std::string const& name,
                 std::string const& dbfilename,
                 std::string const& dbdirectory,
                 std::string const& prefix)
        :
        super( space, model->parameterSpace(), sampling, model->modelName(), name, model->uuid(), prefix ),
        M_model( model ),
        M_prefixModel( model->prefix() ),
        M_expr( expr ),
        M_u( &u ),
        M_u2( &u2 ),
        M_mu( mu ),
        M_crb_built( false ),
        M_mu_sampling( new sampling_type ( model->parameterSpace() , 0 , sampling ) ),
        M_t(),
        M_B(),
        M_offline_error(),
        M_eim(),
        M_write_nl_solutions( boption(_prefix=this->M_prefix, _name="eim.elements.write") ),
        M_write_nl_directory( soption(_prefix=this->M_prefix, _name="eim.elements.directory") )
        {
            if ( model )
                M_eimFeSpaceDb.setModel( model );
            if ( this->functionSpace() )
            {
                M_eimFeSpaceDb.setEimFunctionSpace( this->functionSpace() );

                if ( model_use_nosolve && this->functionSpace()->mesh() )
                {
                    geometricspace_ptrtype geospace( new geometricspace_type( this->functionSpace()->mesh() ) );
                    M_ctxGeoEim.reset( new geometricspace_context_type( geospace ) );
                }
                M_ctxFeBasisEim.reset( new context_type( this->functionSpace() ) );
                M_internalModelFeFunc = this->functionSpace()->elementPtr();
            }

            this->initExprFeContex( mpl::bool_<use_subspace_element>() );

            // update ouput path of database
            if ( dbfilename.empty() )
            {
                //this->setDBFilename( ( boost::format( "%1%.crbdb" ) %this->name() ).str() );
                //this->addDBSubDirectory( "EIMFunction_"+model->modelName() );
                this->addDBSubDirectory( "eim" );
            }
            else
            {
                fs::path dbfilenamePath = dbfilename;
                this->setDBFilename( dbfilenamePath.filename().string() );
                if ( dbfilenamePath.is_relative() )
                {
                    std::string mysubDir;
                    if ( !dbfilenamePath.parent_path().empty() )
                        mysubDir = dbfilenamePath.parent_path().string();

                    fs::path dbdirectoryPath = dbdirectory;
                    if ( !dbdirectoryPath.is_absolute() )
                        dbdirectoryPath = fs::absolute( dbdirectoryPath );
                    this->setDBDirectory( dbdirectoryPath.string() );

                    dbfilenamePath = (dbdirectoryPath/dbfilenamePath).string();

                    if ( !mysubDir.empty() )
                        this->addDBSubDirectory( mysubDir );
                }
                else
                {
                    this->setDBDirectory( dbfilenamePath.parent_path().string() );
                }
            }

            // reload maybe a database
            bool hasLoadedDb = this->loadDB();
            if ( !hasLoadedDb )
                LOG(INFO) << "No EIMFunction database ";
            else
                LOG( INFO ) << "EIMFunction loaded";

            // build eim basis
            if ( this->functionSpace() )
                M_eim.reset( new eim_type( this, sampling , 1e-8, hasLoadedDb, this->M_prefix ) );

            if ( M_write_nl_solutions )
            {
                if ( this->worldComm().isMasterRank() )
                {
                    boost::filesystem::path dir( M_write_nl_directory );
                    if ( boost::filesystem::exists(dir) && boption(_prefix=this->M_prefix, _name="eim.elements.clean-directory" ) )
                    {
                        boost::filesystem::remove_all(dir);
                        boost::filesystem::create_directory(dir);
                    }
                    else if ( !boost::filesystem::exists(dir) )
                    {
                        boost::filesystem::create_directory(dir);
                    }
                }
            }
            this->worldComm().barrier();

        }

    /**
     * storing of shared_ptr model can be necessary (as EimNoSolve)
     */
    void attachModel( model_ptrtype const& model ) { M_modelAttached = model; }

    void initExprFeContex( mpl::false_ /**/ )
        {
            if ( this->model()->functionSpace() )
                M_ctxFeModelSolution.reset( new model_element_expr_context_type( this->model()->functionSpace() ) );
        }
    void initExprFeContex( mpl::true_ /**/ )
        {
            if ( this->model()->functionSpace() )
                M_ctxFeModelSolution.reset( new model_element_expr_context_type( this->model()->functionSpace()->template functionSpace<SubSpaceId>() ) );
        }

    model_ptrtype model() const { return M_model.lock(); }

    bool modelUseSolve() const override { return model_use_solve;}

    auto expr( parameter_type const& mu ) const
    {
        M_mu = mu;
#if !defined(NDEBUG)
        M_mu.check();
#endif
        return M_expr( idv(*M_u), gradv(*M_u2) );
    }

    auto expr( parameter_type const& mu, model_solution_type const& u ) const
    {
        M_mu = mu;
#if !defined(NDEBUG)
        M_mu.check();
#endif
        return this->expr( u,mpl::bool_<use_subspace_element>() );
    }
    auto expr( model_solution_type const& u, mpl::false_ /**/ ) const
    {
        *M_u = u;
        return M_expr( idv(u) );
    }
    auto expr( model_solution_type const& u, mpl::true_ /**/ ) const
    {
        auto const& uexpr = u.template element<SubSpaceId>();
        *M_u = uexpr;
        auto const& uexpr2 = u.template element<SubSpaceId2>();
        if ( SubSpaceId != SubSpaceId2 )
            *M_u2 = uexpr2;
        //return M_expr( idv(uexpr), idv(uexpr2) );
        return M_expr( idv(uexpr), gradv(uexpr2) );
    }

    auto expr( parameter_type const& mu, typename rbfunctionspace_context_type::rbspace_type::element_type const& urb ) const
    {
        M_mu = mu;
#if !defined(NDEBUG)
        M_mu.check();
#endif
        // return M_expr( idv(urb) );
        return M_expr( idv(urb), gradv(urb) );
    }
    auto expr( parameter_type const& mu, typename rbfunctionspace_context_type::rbspace_type::element_type const& urb,
               typename rbfunctionspace_context2_type::rbspace_type::element_type const& urb2 ) const
    {
        M_mu = mu;
#if !defined(NDEBUG)
        M_mu.check();
#endif
        //return M_expr( idv(urb) );
        return M_expr( idv(urb), gradv(urb2) );
    }

    //!
    //!
    //!
    void loadDB( std::string const& filename, crb::load l ) override {}

    void saveDB() override
    {
        if ( this->worldComm().isMasterRank() )
        {
            fs::ofstream ofs( this->dbLocalPath() / this->dbFilename() );
            if ( ofs )
            {
                //boost::archive::text_oarchive oa( ofs );
                boost::archive::binary_oarchive oa( ofs );
                // write class instance to archive
                oa << *this;
                // archive and stream closed when destructors are called
            }
        }

        M_eimFeSpaceDb.save( this->name() + "_fespace", this->dbDirectory(), this->worldComm() );
    }

    /**
     * load the EIM database
     */
    bool loadDB() override
    {
        if( M_eimFeSpaceDb.q().size() > 0 && M_t.size() > 0 )
        {
            return true;
        }

        fs::path db = this->lookForDB();
        if ( db.empty() || !fs::exists( db ) )
        {
            return false;
        }

        // if ( !fs::exists( db ) )
        //     return false;

        fs::ifstream ifs( db );
        if ( ifs )
        {
            //boost::archive::text_iarchive ia( ifs );
            boost::archive::binary_iarchive ia( ifs );
            //write class instance to archive
            ia >> *this;
            //std::cout << "Loading " << db << " done...\n";

            // add interpolation points to the context
            if( M_ctxFeBasisEim && M_ctxFeBasisEim->nPoints()==0 && M_ctxFeModelSolution )
            {
                typename Feel::node<value_type>::type no(nDim);
                for(int i=0 ; i<M_t.size(); i++)
                {
                    node_type t = M_t[i];
                    for(int dim =0; dim < nDim; ++dim )
                        no(dim) = t(dim);
                    M_ctxFeBasisEim->add( no );
                    M_ctxFeModelSolution->add( no );
                }
            }

            // load eim fe space
            bool hasLoadedFeSpaceDb = M_eimFeSpaceDb.load( this->name() + "_fespace", this->dbDirectory(), this->worldComm() );
            if ( this->functionSpace() )
                CHECK( hasLoadedFeSpaceDb ) << "loading of fespace database fails";

            this->setIsLoaded( true );
            // archive and stream closed when destructors are called
            return true;
        }

        return false;
    }

    void initializeDataStructures() override
    {
        M_mu_sampling->clear();
        if ( M_ctxGeoEim )
            M_ctxGeoEim->removeCtx();
        M_ctxFeBasisEim->removeCtx();
        M_ctxFeModelSolution->removeCtx();
        M_t.clear();
        M_M_max=0;
        M_B.resize(0,0);
        M_eimFeSpaceDb.clear();
    }

    vector_type beta( parameter_type const& mu ) /*const*/ override { return this->beta( mu, this->mMax() ); }
    vector_type beta( parameter_type const& mu, model_solution_type const& T ) /*const*/ override { return this->beta( mu, T, this->mMax() ); }
    vector_type beta( parameter_type const& mu, vectorN_type const& urb ) /*const*/ override { return this->beta( mu, urb, this->mMax() ); }

    vector_type
    beta( parameter_type const& mu, size_type __M ) override
    {

        vector_type __beta( __M );
        if ( model_use_nosolve && M_ctxGeoEim )
            __beta = this->operator()( *M_ctxGeoEim, mu , __M );
        else
        {
            DCHECK( M_ctxFeModelSolution ) << "no fe context";
            __beta = this->operator()( *M_ctxFeModelSolution, mu , __M );
        }
        DCHECK( __beta.size() == __M ) << "Invalid size beta: " << __beta.size() << " M=" << __M  << " beta = " << __beta << "\n";
        this->M_B.block(0,0,__M,__M).template triangularView<Eigen::UnitLower>().solveInPlace(__beta);

        return __beta;
    }

    vector_type
    beta( parameter_type const& mu, model_solution_type const& T, size_type __M ) override
    {
        DCHECK( M_ctxFeModelSolution ) << "no fe context";
        // beta=B_M\g(Od(indx),mut(i))'
        vector_type __beta( __M );
        __beta = this->operator()( T, *M_ctxFeModelSolution, mu , __M );
        DCHECK( __beta.size() == __M ) << "Invalid size beta: " << __beta.size() << " M=" << __M  << " beta = " << __beta << "\n";

        this->M_B.block(0,0,__M,__M).template triangularView<Eigen::UnitLower>().solveInPlace(__beta);
        return __beta;
    }

    vector_type
    beta( parameter_type const& mu, vectorN_type const& urb, size_type __M ) override
    {
        if ( !M_ctxRbModelSolution )
            return this->beta(mu,__M );

        DCHECK( M_ctxRbModelSolution ) << "no rbspace context";
        vector_type __beta( __M );
        if ( SubSpaceId != SubSpaceId2 && M_ctxRbModelSolution2 )
            __beta = this->operator()( urb, *M_ctxRbModelSolution, *M_ctxRbModelSolution2, mu , __M );
        else
            __beta = this->operator()( urb, *M_ctxRbModelSolution, mu , __M );
        this->M_B.block(0,0,__M,__M).template triangularView<Eigen::UnitLower>().solveInPlace(__beta);
        return __beta;
    }


    model_functionspace_ptrtype modelFunctionSpace() const override { return this->model()->functionSpace();}
    model_functionspace_ptrtype modelFunctionSpace() override { return this->model()->functionSpace();}

    void addInterpolationPoint( node_type const& t ) override
    {
        M_t.push_back( t );
        typename Feel::node<value_type>::type no(nDim);
        for(int i =0;i < nDim; ++i ) no(i) = M_t.back()(i);
        // add in precompute object the last magic point
        if ( M_ctxGeoEim )
            M_ctxGeoEim->add( no );
        M_ctxFeBasisEim->add( no );
        M_ctxFeModelSolution->add( no );
        std::for_each( M_t.begin(), M_t.end(), []( node_type const& t ) { DVLOG(2) << "t=" << t << "\n"; } );
    }

    bool hasRbSpaceContext() const override
    {
        return M_ctxRbModelSolution.use_count() > 0;
    }
    bool hasRbSpaceContext2() const override
    {
        return M_ctxRbModelSolution2.use_count() > 0;
    }
    boost::any rbSpaceContext() const override
    {
        return M_ctxRbModelSolution;
    }
    boost::any rbSpaceContext2() const override
    {
        return M_ctxRbModelSolution2;
    }

    void updateRbSpaceContext( boost::any const& rbspacebase ) override
    {
        // no update if model is not a modelcrbbase
        if ( !model_is_modelcrbbase )
            return;

        this->updateRbSpaceContext( rbspacebase, mpl::bool_<use_subspace_element>() );
    }
    void updateRbSpaceContext( boost::any const& rbspacebase, mpl::false_ )
    {
        if ( !boost::any_cast<rbfunctionspace_ptrtype>( &rbspacebase ) )
        {
            std::cout << "[EIMFunction::updateRbSpaceContext] cast fails with rbfunctionspace_ptrtype\n";
            return;
        }
        //std::cout << "[EIMFunction::updateRbSpaceContext] cast ok\n";
        rbfunctionspace_ptrtype rbspace = boost::any_cast<rbfunctionspace_ptrtype>( rbspacebase );

        M_ctxRbModelSolution.reset( new rbfunctionspace_context_type( rbspace ) );
        typename Feel::node<value_type>::type no(nDim);
        for (auto const& pt : M_t )
        {
            for(int i =0;i < nDim; ++i ) no(i) = pt(i);
            M_ctxRbModelSolution->add( no );
        }
        M_ctxRbModelSolution->update();
    }
    void updateRbSpaceContext( boost::any const& rbspacebase, mpl::true_ )
    {
        if ( !boost::any_cast<rbfunctionspace_ptrtype>( &rbspacebase ) )
        {
            std::cout << "[EIMFunction::updateRbSpaceContext (composite)] cast fails with rbfunctionspace_ptrtype\n";
            return;
        }

        rbfunctionspace_ptrtype rbspace = boost::any_cast<rbfunctionspace_ptrtype>( rbspacebase );

        M_ctxRbModelSolution.reset( new rbfunctionspace_context_type( rbspace->template rbFunctionSpace<SubSpaceId>() ) );
        typename Feel::node<value_type>::type no(nDim);
        for (auto const& pt : M_t )
        {
            for(int i =0;i < nDim; ++i ) no(i) = pt(i);
            M_ctxRbModelSolution->add( no );
        }
        M_ctxRbModelSolution->update();

        if ( SubSpaceId != SubSpaceId2 )
        {
            M_ctxRbModelSolution2.reset( new rbfunctionspace_context2_type( rbspace->template rbFunctionSpace<SubSpaceId2>() ) );
            typename Feel::node<value_type>::type no(nDim);
            for (auto const& pt : M_t )
            {
                for(int i =0;i < nDim; ++i ) no(i) = pt(i);
                M_ctxRbModelSolution2->add( no );
            }
            M_ctxRbModelSolution2->update();
        }

    }
    void setRbSpaceContext( boost::any const& rbCtxBase ) override
    {
        if ( !boost::any_cast<rbfunctionspace_context_ptrtype>( &rbCtxBase ) )
        {
            std::cout << "[EIMFunction::setRbSpaceContext] cast fails with rbfunctionspace_context_ptrtype\n";
            return;
        }
        //std::cout << "[EIMFunction::setRbSpaceContext] cast ok\n";
        M_ctxRbModelSolution = boost::any_cast<rbfunctionspace_context_ptrtype>( rbCtxBase );
    }
    void setRbSpaceContext2( boost::any const& rbCtxBase ) override
    {
        if ( !boost::any_cast<rbfunctionspace_context2_ptrtype>( &rbCtxBase ) )
        {
            std::cout << "[EIMFunction::setRbSpaceContext2] cast fails with rbfunctionspace_context2_ptrtype\n";
            return;
        }
        //std::cout << "[EIMFunction::setRbSpaceContext2] cast ok\n";
        M_ctxRbModelSolution2 = boost::any_cast<rbfunctionspace_context2_ptrtype>( rbCtxBase );
    }

    node_type interpolationPoint( int position ) const override
    {
        int size = M_t.size();
        DCHECK( position < size ) << "Invalid point position: " << position << " M_t.size() =" << M_t.size() << "\n";
        return M_t[position];
    }

    void addBasis( element_type const &q ) override
    {
        M_eimFeSpaceDb.addBasis( q );
    }
    void addExpressionEvaluation( element_type const &g ) override
    {
        M_eimFeSpaceDb.addExpressionEvaluation( g );
    }
    void addZ( element_type const &z ) override
    {
        M_eimFeSpaceDb.addZ( z );
    }
    void addSolution( model_solution_type const &solution ) override
    {
        M_eimFeSpaceDb.addSolution( solution );
    }
    void addParameter( parameter_type const &mu ) override
    {
        M_mu_sampling->addElement( mu );
    }
    void clearParameterSampling() override
    {
        M_mu_sampling->clear();
    }
    void addOfflineError( double error ) override
    {
        M_offline_error.push_back( error );
    }
    void clearOfflineError() override
    {
        M_offline_error.resize(0);
    }

    model_solution_type solve( parameter_type const&  mu ) override
    {
        M_mu = mu;
#if !defined(NDEBUG)
        M_mu.check();
#endif
        model_solution_type u = this->model()->functionSpace()->element();
        bool need_solve = true;

        if ( M_write_nl_solutions )
        {
            need_solve = !u.load( _path=M_write_nl_directory,
                                  _suffix=std::to_string(mu.key()), _type="hdf5" );
            if ( need_solve )
                LOG(INFO) << "EIM : Unable to load nl solution in direcotry "
                          << M_write_nl_directory << ", for parameter : " << mu.toString()
                          <<" / " << mu.key()<< ". Solve function will be called.";
            else
                LOG(INFO) << "EIM : NL solution loaded in direcotry "
                          << M_write_nl_directory << ", for parameter : " << mu.toString()
                          <<" / " << mu.key();
        }

        if ( need_solve )
        {
            LOG(INFO) << "EIM : calling solve function for parameter " << mu.toString()
                      <<" / " << mu.key();
            u = this->model()->solve(mu);

            if ( M_write_nl_solutions )
            {
                LOG(INFO) << "EIM : Wrting solution on disk in directory "
                          << M_write_nl_directory << ", for parameter : " << mu.toString()
                          <<" / " << mu.key();
                u.save( _path=M_write_nl_directory,
                        _suffix=std::to_string(mu.key()), _type="hdf5" );
            }
        }

        return u;
    }

    element_type operator()( parameter_type const&  mu ) override
        {
            model_solution_type sol;
            if( ioption(_prefix=this->M_prefixModel,_name="ser.eim-frequency") != 0 ) //Use SER
            {
                if( boption(_prefix=this->M_prefixModel,_name="ser.use-rb-in-eim-basis-build") )
                    sol = this->computeRbExpansion( mu ); //RB
                else
                    sol = this->computePfem( mu ); //PFEM
            }
            else
                sol = this->solve( mu ); //FEM
            auto eimexpr = this->expr( mu, sol );
            //LOG(INFO) << "operator() mu=" << mu << "\n" << "sol=" << M_u << "\n";
            return vf::project( _space=this->functionSpace(), _range=this->functionSpace()->template rangeElements<0>(), _expr=eimexpr );
        }
    element_type operator()( model_solution_type const& T, parameter_type const&  mu ) override
        {
            auto eimexpr = this->expr( mu, T );
            return vf::project( _space=this->functionSpace(), _range=this->functionSpace()->template rangeElements<0>(), _expr=eimexpr );
        }

    vector_type operator()( geometricspace_context_type const& ctx, parameter_type const& mu , int M )
        {
            auto eimexpr = this->expr( mu );
            bool applyProjection = false;//true;
            bool doMpiComm = false;//true;//false;
            return evaluateFromContext( _context=ctx, _expr=eimexpr, _max_points_used=M,
                                        _mpi_communications=doMpiComm, _projection=applyProjection );
        }

    vector_type operator()( model_element_expr_context_type const& ctx, parameter_type const& mu , int M)
        {
            model_solution_type sol;
            if( ioption(_prefix=this->M_prefixModel,_name="ser.eim-frequency") != 0 ) //Use SER
            {
                if( boption(_prefix=this->M_prefixModel,_name="ser.use-rb-in-eim-basis-build") )
                    sol = this->computeRbExpansion( mu ); //RB
                else
                    sol = this->computePfem( mu ); //PFEM
            }
            else
                sol = this->solve( mu ); //FEM

            auto eimexpr = this->expr( mu, sol );

            //auto projected_expr = vf::project( _space=this->functionSpace(), _expr=M_expr );
            //return evaluateFromContext( _context=ctx, _expr=idv(projected_expr) , _max_points_used=M );
            bool applyProjection = false;//true;
            bool doMpiComm = true;//false;
            return evaluateFromContext( _context=ctx, _expr=eimexpr, _max_points_used=M,
                                        _mpi_communications=doMpiComm, _projection=applyProjection );
        }
    vector_type operator()( model_solution_type const& T, model_element_expr_context_type const& ctx, parameter_type const& mu , int M)
        {
            auto eimexpr = this->expr( mu,T );
            //auto projected_expr = vf::project( _space=this->functionSpace(), _expr=M_expr );
            //return evaluateFromContext( _context=ctx, _expr=idv(projected_expr) , _max_points_used=M );
            bool applyProjection= false;//true;
            bool doMpiComm = true;//false;
            return evaluateFromContext( _context=ctx, _expr=eimexpr , _max_points_used=M,
                                        _mpi_communications=doMpiComm, _projection=applyProjection );
        }

    vector_type operator()( vectorN_type const& urb, rbfunctionspace_context_type const& ctx, parameter_type const& mu , int M )
        {
            auto urbelt = ctx.rbFunctionSpace()->element();
            int dimRb = std::min((int)urbelt.size(),(int)urb.size());
            for ( int k=0; k<dimRb; ++k )
                urbelt(k) = urb(k);
            auto eimexpr = this->expr( mu,urbelt );
            return evaluateFromContext( _context=ctx, _expr=eimexpr , _max_points_used=M,
                                        _mpi_communications=false, _projection=false );
        }

    vector_type operator()( vectorN_type const& urb, rbfunctionspace_context_type const& ctx, rbfunctionspace_context2_type const& ctx2, parameter_type const& mu , int M )
        {
            auto urbelt = ctx.rbFunctionSpace()->element();
            int dimRb = std::min((int)urbelt.size(),(int)urb.size());
            for ( int k=0; k<dimRb; ++k )
                urbelt(k) = urb(k);
            auto urbelt2 = ctx2.rbFunctionSpace()->element();
            int dimRb2 = std::min((int)urbelt2.size(),(int)urb.size());
            for ( int k=0; k<dimRb2; ++k )
                urbelt2(k) = urb(k);
            auto eimexpr2 = this->expr( mu,urbelt,urbelt2 );
            return evaluateFromContext( _context=ctx, _context2=ctx2, _expr=eimexpr2 , _max_points_used=M,
                                        _mpi_communications=false, _projection=false );
        }

    vector_type computeExpansionCoefficients( parameter_type const& mu, model_solution_type const& solution,  int M) override
    {
        vector_type rhs( M );

        auto eimexpr = this->expr( mu, solution );

        if( 0 )// boption(_name="eim.compute-expansion-of-expression") )
        {
            rhs = evaluateFromContext( _context=*M_ctxFeModelSolution, _expr=eimexpr );
        }
        else
        {
            //auto proj_g = vf::project( _space=this->functionSpace(), _expr=M_expr );
            //rhs = proj_g.evaluate( M_ctx );
            rhs = evaluateFromContext( _context=*M_ctxFeModelSolution, _expr=eimexpr , _max_points_used=M, _projection=true );
        }

        M_B.block(0,0,M,M).template triangularView<Eigen::UnitLower>().solveInPlace(rhs);
        DVLOG(2) << "solve B sol = rhs with rhs = \n" << rhs <<"\n";

        return rhs;
    }


    vector_type evaluateExpressionAtInterpolationPoints(model_solution_type const &solution, parameter_type const& mu, int M) override
    {
        auto eimexpr = this->expr( mu, solution );
        return evaluateFromContext( _context=*M_ctxFeModelSolution, _expr=eimexpr , _max_points_used=M, _projection=true );
    }
    vector_type evaluateElementAtInterpolationPoints(element_type const & element, int M) override
    {
        return evaluateFromContext( _context=*M_ctxFeBasisEim, _expr=idv(element) , _max_points_used=M, _projection=true );
    }

    // compute the maximum of the residual using either real expression
    // or its projection on the functionspace
    // mu : parameter
    // z  : expansion in the EIM basis
    // g  : projection of the expression
    boost::tuple<double,node_type> computeMaximumOfResidual( parameter_type const& mu, model_solution_type const& solution, element_type const& z) override
    {
        double max=0;
        node_type node(mesh_type::nDim);
        for(int d=0; d<nDim; d++) node(d)=0;

        auto eimexpr = this->expr( mu, solution );
        auto residual_expr = eimexpr - idv(z);

        // auto proj_g = vf::project( _space=this->functionSpace(), _expr=M_expr );
        // auto residual_projected_expr = idv(proj_g)-idv(z);
        if ( !this->M_computeExpansionOfExpression || ( this->M_normUsedForResidual == super::ResidualNormType::LinftyVec ) )
        {
            this->M_internalModelFeFunc->on(_range=this->functionSpace()->template rangeElements<0>(),_expr=residual_expr );
            // this->M_internalModelFeFunc->on(_range=elements(this->mesh()),_expr=M_expr );
            // this->M_internalModelFeFunc->add(-1., z );
        }
        auto residual_projected_expr = idv( this->M_internalModelFeFunc );

        switch ( this->M_normUsedForResidual )
        {
        case super::ResidualNormType::Linfty:
        {
            if( this->M_computeExpansionOfExpression )
                boost::tie( max, node ) = normLinf( _range=this->functionSpace()->template rangeElements<0>(), _pset=_Q<0>(), _expr= residual_expr );
            else
                boost::tie( max, node ) = normLinf( _range=this->functionSpace()->template rangeElements<0>(), _pset=_Q<0>(), _expr= residual_projected_expr );
        }
        break;
        case super::ResidualNormType::L2:
        {
            if( this->M_computeExpansionOfExpression )
                max = normL2( _range=this->functionSpace()->template rangeElements<0>(), _expr=residual_expr );
            else
                max = normL2( _range=this->functionSpace()->template rangeElements<0>(), _expr=residual_projected_expr );
        }
        break;
        case super::ResidualNormType::LinftyVec:
        {
            max = this->M_internalModelFeFunc->linftyNorm();
        }
        break;
        }

        return boost::make_tuple( max , node );
    }

    boost::tuple<double,node_type> computeMaximumOfExpression( parameter_type const& mu, model_solution_type const& solution ) override
    {
        double max=0;
        node_type node(mesh_type::nDim);
        for(int d=0; d<nDim; d++) node(d)=0;

        auto eimexpr = this->expr( mu, solution );
#if 1
        if ( !this->M_computeExpansionOfExpression || ( this->M_normUsedForResidual == super::ResidualNormType::LinftyVec ) )
            this->M_internalModelFeFunc->on(_range=this->functionSpace()->template rangeElements<0>(),_expr=eimexpr );
        auto projected_expr = idv( this->M_internalModelFeFunc );

        switch ( this->M_normUsedForResidual )
        {
        case super::ResidualNormType::Linfty:
        {
            if( this->M_computeExpansionOfExpression )
                boost::tie( max, node ) = normLinf( _range=this->functionSpace()->template rangeElements<0>(), _pset=_Q<0>(), _expr= eimexpr );
            else
                boost::tie( max, node ) = normLinf( _range=this->functionSpace()->template rangeElements<0>(), _pset=_Q<0>(), _expr= projected_expr );
        }
        break;
        case super::ResidualNormType::L2:
        {
            if( this->M_computeExpansionOfExpression )
                max = normL2( _range=this->functionSpace()->template rangeElements<0>(), _expr=eimexpr );
            else
                max = normL2( _range=this->functionSpace()->template rangeElements<0>(), _expr=projected_expr );
        }
        break;
        case super::ResidualNormType::LinftyVec:
        {
            max = this->M_internalModelFeFunc->linftyNorm();
        }
        break;
        }

#else
        auto proj_g = vf::project( _space=this->functionSpace(), _range=this->functionSpace()->template rangeElements<0>(),_expr=eimexpr/*M_expr*/ );
        auto projected_expr = idv( proj_g );

        std::string norm_used = soption(_prefix=this->M_prefix,_name="eim.norm-used-for-residual");
        bool check_name_norm = false;
        DVLOG( 2 ) << "[computeMaximumOfExpression] norm used : "<<norm_used;
        if( norm_used == "Linfty" )
        {
            check_name_norm=true;
            if( boption(_prefix=this->M_prefix,_name="eim.compute-expansion-of-expression") )
            {
                auto exprmax = normLinf( _range=this->functionSpace()->template rangeElements<0>(), _pset=_Q<0>(), _expr= eimexpr );
                max = exprmax.template get<0>();
                node = exprmax.template get<1>();
            }
            else
            {
                auto exprmax = normLinf( _range=this->functionSpace()->template rangeElements<0>(), _pset=_Q<0>(), _expr=projected_expr);
                max = exprmax.template get<0>();
                node = exprmax.template get<1>();
            }
        }
        if( norm_used == "L2" )
        {
            check_name_norm=true;
            if( boption(_prefix=this->M_prefix,_name="eim.compute-expansion-of-expression") )
            {
                double norm = math::sqrt( integrate( _range=this->functionSpace()->template rangeElements<0>() ,_expr=eimexpr*eimexpr).evaluate()( 0,0 ) );
                max = norm;
            }
            else
            {
                double norm = math::sqrt( integrate( _range=this->functionSpace()->template rangeElements<0>() ,_expr=projected_expr*projected_expr).evaluate()( 0,0 ) );
                max = norm;
            }
        }
        if( norm_used == "LinftyVec" )
        {
            check_name_norm=true;
            if( boption(_prefix=this->M_prefix,_name="eim.compute-expansion-of-expression") )
            {
                auto projection = vf::project( _space=this->functionSpace(), _range=this->functionSpace()->template rangeElements<0>(),_expr=eimexpr );
                double norm = projection.linftyNorm();
                max = norm ;
            }
            else
            {
                auto projection = vf::project( _space=this->functionSpace(), _range=this->functionSpace()->template rangeElements<0>(),_expr=projected_expr );
                double norm = projection.linftyNorm();
                max = norm ;
            }
        }
        CHECK( check_name_norm ) <<"[EIM options] The name of the norm "<<norm_used<<" is not known\n";
#endif

        return boost::make_tuple( max , node );
    }

    element_type residual( element_type const& z, element_type const& g ) override
    {
        auto residual_projected_expr = idv(g) - idv(z);

        //element_type projection;
        //projection = vf::project( _space=this->functionSpace(),_expr=residual_projected_expr );
        element_type projection = this->functionSpace()->element( residual_projected_expr );

        return projection;
    }


    void fillInterpolationMatrixFirstTime(  ) override
    {
        // update interpolation matrix
        // TODO: update only the new line and eventually the new column rather than recomputing everything
        M_B.resize(1 , 1 );
        M_B(0,0) = 1;

        CHECK( M_ctxFeBasisEim->nPoints() == 1 );

        if( this->M_computeExpansionOfExpression )
        {
            M_mu = M_mu_sampling->at(0);
            auto eimexpr = this->expr( M_mu, M_eimFeSpaceDb.solution(0) );

            auto expression_evaluated = evaluateFromContext( _context=*M_ctxFeModelSolution, _expr=eimexpr );
            double eval = expression_evaluated(0);
            auto expression_q = eimexpr / eval ; // normalization
            auto q_evaluated = evaluateFromContext( _context=*M_ctxFeModelSolution, _expr=expression_q );
            eval = q_evaluated(0);
            CHECK( math::abs(eval-1) < 1e-10 )
                << "q[0](t[0])) != 1 " << "q[0] = " << eval << " t[0] = \n "<<M_ctxFeBasisEim->node(0)<<" \n" ;
        }
        else
        {
            auto eimexpr = this->expr( M_mu );//, *M_u );
            auto projected_g = vf::project( _space=this->functionSpace(), _range=this->functionSpace()->template rangeElements<0>(),_expr=eimexpr );

            auto projected_g_evaluated = evaluateFromContext( _context=*M_ctxFeBasisEim, _expr=idv( projected_g ) );
            double eval = projected_g_evaluated(0);
            projected_g.scale( 1./eval );
            auto q_evaluated = evaluateFromContext( _context=*M_ctxFeBasisEim, _expr=idv( projected_g ) );
            eval = q_evaluated(0);
            CHECK( math::abs(eval-1) < 1e-10 )
                << "q[0](t[0])) != 1 " << "q[0] = " << eval << " t[0] = \n "<<M_ctxFeBasisEim->node(0)<<" \n" ;
        }

        DVLOG( 2 )<<" M_B : \n "<<M_B;
        saveDB();

    }//fillInterpolationMatrixFirstTime

    //build eim m^th eim basis function
    //m : index of the eim basis
    //q : vector of projected on functionspace eim basis functions
    boost::any buildBasisFunction( int m ) override
    {
        //rank of the current processor
        int proc_number = this->worldComm().globalRank();

        int npoints = M_ctxFeBasisEim->nPoints();

        CHECK( m < npoints ) << "m = "<<m<<" and there is "<<npoints<<" interpolation points" ;
        if( this->M_computeExpansionOfExpression )
        {
            auto interpolation_point = this->functionSpace()->context();
            int proc_having_the_point ;

            //we want to normalize the expression with its evaluation at
            //the interpolation point
            interpolation_point.clear();
            proc_having_the_point = M_ctxFeBasisEim->processorHavingPoint( m );
            if( proc_number == proc_having_the_point )
            {
                auto basis = M_ctxFeBasisEim->at( m );
                interpolation_point.addCtx( basis , proc_number );
            }

            vector_type expression_evaluated;

            M_mu = M_mu_sampling->at(m);
            auto eimexpr = this->expr( M_mu, M_eimFeSpaceDb.solution(m) );
            auto residual = eimexpr - idv( M_eimFeSpaceDb.z(m) );
            expression_evaluated = evaluateFromContext( _context=interpolation_point , _expr=residual );
            double eval = expression_evaluated( 0 );
            auto expression_q  = residual / eval; // __j^th normalized basis function

            return expression_q;
        }
        else
        {
            //in that case we have already the eim basis functions
            return M_eimFeSpaceDb.q(m);
        }

        return 0;
    }

    element_type Residual( int m, element_type const& z) override
    {
        if( this->M_computeExpansionOfExpression )
        {
            M_mu = M_mu_sampling->at(m);
            auto eimexpr = this->expr( M_mu, M_eimFeSpaceDb.solution(m) );
            auto residual =  eimexpr - idv(z);

            int proc_number = this->worldComm().globalRank();
            auto interpolation_point = this->functionSpace()->context();
            int proc_having_the_point ;
            //we want to normalize the expression with its evaluation at
            //the interpolation point
            interpolation_point.clear();
            proc_having_the_point = M_ctxFeBasisEim->processorHavingPoint( m );
            if( proc_number == proc_having_the_point )
            {
                auto basis = M_ctxFeBasisEim->at( m );
                interpolation_point.addCtx( basis , proc_number );
            }
            vector_type expression_evaluated;
            element_type res;
            if( m == 0 )
            {
                expression_evaluated = evaluateFromContext( _context=interpolation_point , _expr=eimexpr );
                double eval = expression_evaluated( 0 );
                res = vf::project( _space=this->functionSpace(), _range=this->functionSpace()->template rangeElements<0>(), _expr=eimexpr/eval );
            }
            else
            {
                expression_evaluated = evaluateFromContext( _context=interpolation_point , _expr=residual );
                double eval = expression_evaluated( 0 );
                res = vf::project( _space=this->functionSpace(), _range=this->functionSpace()->template rangeElements<0>(), _expr=residual/eval );
            }
            return res;
        }
        else
        {
            auto residual = vf::project(_space=this->functionSpace(), _range=this->functionSpace()->template rangeElements<0>(), _expr=idv(M_eimFeSpaceDb.g(m)) - idv( z ) );
            auto t = M_t[m];
            residual.scale( 1./residual(t)(0,0,0) );
            return residual;
        }
    }


    void fillInterpolationMatrix() override
    {
        //rank of the current processor
        int proc_number = this->worldComm().globalRank();

        int size = M_ctxFeBasisEim->nPoints();
        M_B.conservativeResize( size , size );

        DVLOG( 2 ) << "solution.size() : "<<M_eimFeSpaceDb.solutions().size();
        DVLOG( 2 ) << "q.size() : "<< M_eimFeSpaceDb.q().size();

        auto eimexpr = this->expr( M_mu );//, *M_u );

        auto element_zero = project( _space=this->functionSpace(), _range=this->functionSpace()->template rangeElements<0>(), _expr=cst(0) );
        auto example_of_q = ( eimexpr - idv(element_zero) )/1.0;

        typedef decltype( example_of_q ) expression_type ;

        for( int __j = 0; __j < size; ++__j )
        {
            auto any_q = buildBasisFunction( __j );
            if( this->M_computeExpansionOfExpression )
            {

                if (!boost::any_cast<expression_type>(&any_q))
                    throw std::logic_error( "[EIM::fillInterpolationMatrix] not cast possible for eim basis function expression" );
                auto q = boost::any_cast<expression_type>( any_q );
                M_B.col( __j ) = evaluateFromContext( _context=*M_ctxFeBasisEim, _expr=q );
            }
            else
            {
                if (!boost::any_cast<element_type>(&any_q))
                    throw std::logic_error( "[EIM::fillInterpolationMatrix] not cast possible for eim basis function projected in function space" );
                auto q = boost::any_cast<element_type>( any_q );
                M_B.col( __j ) = q.evaluate( *M_ctxFeBasisEim );
            }
        }

        DVLOG( 2 )<<" M_B : \n "<<M_B;
        //google::FlushLogFiles(google::GLOG_INFO);
        saveDB();

    }//fillInterpolationMatrix

    /*
     computeRbExpansion returns RB solution expansion to be used with SER (return FE solution is RB space is not available)
     If option ser.corrected-rb=true, return the RB solution corrected with Riesz representation of RB residual
     */
    model_solution_type computeRbExpansion( parameter_type const& mu ) override
    {
        if( this->rbBuilt() )
            return computeRbExpansion( mu, typename boost::is_base_of<ModelCrbBaseBase,model_type>::type() );
        else
            return this->solve( mu );
    }
    model_solution_type computeRbExpansion( parameter_type const& mu, boost::mpl::bool_<false>)
    {
        return this->model()->functionSpace()->element();
    }
    model_solution_type computeRbExpansion( parameter_type const& mu, boost::mpl::bool_<true>)
    {
        //Compute RB approx (online)
        std::vector<vectorN_type> uN, uNdu, uNold, uNduold;
        auto o = M_crb->lb( M_crb->dimension(), mu, uN, uNdu , uNold, uNduold );
        return computeRbExpansion( mu, uN );
    }

    model_solution_type computeRbExpansion( parameter_type const& mu, std::vector<vectorN_type> uN) override
    {
        if( this->rbBuilt() )
            return computeRbExpansion( mu, uN, typename boost::is_base_of<ModelCrbBaseBase,model_type>::type() );
        else
            return this->solve( mu );
    }
    model_solution_type computeRbExpansion( parameter_type const& mu, std::vector<vectorN_type> uN, boost::mpl::bool_<false>)
    {
        return this->model()->functionSpace()->element();
    }
    model_solution_type computeRbExpansion( parameter_type const& mu, std::vector<vectorN_type> uN, boost::mpl::bool_<true>)
    {
        int size = uN.size();
        if( size!= 0 )
        {
            // Compute RB expansion from uN
            model_solution_type sol = M_crb->expansion( uN[size-1], M_crb->dimension(), false );

            this->M_rb_online_iterations.push_back( M_crb->online_iterations().first );
            this->M_rb_online_increments.push_back( M_crb->online_iterations().second );
            return sol;
        }
        else
            return computeRbExpansion( mu );
    }

    /* computePfem returns PFEM solution (FE with affine decomposition) */
    model_solution_type computePfem( parameter_type const& mu ) override
    {
        if( M_crb )
            return computePfem( mu, typename boost::is_base_of<ModelCrbBaseBase,model_type>::type() );
        else
            return this->solve( mu );
    }
    model_solution_type computePfem( parameter_type const& mu, boost::mpl::bool_<false> )
    {
        return this->model()->functionSpace()->element();
    }
    model_solution_type computePfem( parameter_type const& mu, boost::mpl::bool_<true> )
    {
        return M_crb->solveFemModelUsingAffineDecomposition( mu );
    }

    boost::tuple<double,std::vector<vectorN_type> > RieszResidualNorm( parameter_type const& mu ) override
    {
        if( this->rbBuilt() )
            return RieszResidualNorm( mu, typename boost::is_base_of<ModelCrbBaseBase,model_type>::type() );
        else
            return boost::make_tuple(0,std::vector<vectorN_type>());
    }
    boost::tuple<double,std::vector<vectorN_type> > RieszResidualNorm( parameter_type const& mu, boost::mpl::bool_<false>)
    {
        return boost::make_tuple(0,std::vector<vectorN_type>());
    }
    boost::tuple<double,std::vector<vectorN_type> > RieszResidualNorm( parameter_type const& mu, boost::mpl::bool_<true>)
    {
        //Compute RB approx (online)
        std::vector<vectorN_type> uN, uNdu, uNold, uNduold;
        auto o = M_crb->lb( M_crb->dimension(), mu, uN, uNdu , uNold, uNduold );
        auto error = M_crb->computeRieszResidualNorm( mu, uN );
        return boost::make_tuple( error, uN );
    }

    sampling_ptrtype createSubTrainset( sampling_ptrtype const& trainset, int method ) override
    {
        sampling_ptrtype subtrainset( new sampling_type( this->model()->parameterSpace() ) );
        std::vector<parameter_type> subvector;

        // Compute min/max of residual (Riesz) norm
        std::vector<double> norm( trainset->size() );
        int i=0;
        for( auto const& mu : *trainset )
        {
            norm[i] = RieszResidualNorm( mu ).template get<0>();
            i++;
        }
        auto min = min_element(norm.begin(), norm.end());
        auto max = max_element(norm.begin(), norm.end());

        auto min_idx = std::distance(norm.begin(), min);
        auto max_idx = std::distance(norm.begin(), max);

        // Deduce the tolerance used to select inputs of the sub-trainset
        // TODO : use other expr for tol(min,max) ?
        double tol = ( *min + *max )/2.;

        // Proc zero build the vector of the selected inputs
        if( this->worldComm().isMasterRank() )
        {
            std::cout << "[eim::createSubTrainset] error indicator min = " << *min << " reached at mu = " << trainset->at(min_idx) << std::endl;
            std::cout << "[eim::createSubTrainset] error indicator max = " << *max << " reached at mu = " << trainset->at(max_idx) << std::endl;
            std::cout << "[eim::createSubTrainset] error indicator tol(mean) = " << tol << std::endl;

            for(int i=0; i<trainset->size(); i++)
            {
                auto mu = trainset->at(i);
                if ( method == 1 && norm[i] < tol )
                    subvector.push_back( mu );
                else if ( method == 2 && norm[i] > tol )
                    subvector.push_back( mu );
            }
        }
        // All procs have the same sub-trainset
        boost::mpi::broadcast( this->worldComm() , subvector, this->worldComm().masterRank() );
        subtrainset->setElements( subvector );
        if( this->worldComm().isMasterRank() )
            std::cout << "[createSubTrainset] Original trainset size = " << trainset->size() << ", sub-trainset size = " << subtrainset->size() << "\n";
        if( subtrainset->size() != 0 )
            return subtrainset;
        else
        {
            if( this->worldComm().isMasterRank() )
                std::cout << "[createSubTrainset] No elements added in subtrainset, consider full trainset \n";
            return trainset;
        }
    }

    //Let g the expression that we want to have an eim expasion
    //here is computed its l2 norm : || g ||_L2
    double expressionL2Norm( model_solution_type const& T , parameter_type const& mu ) const override
    {
        auto eimexpr = this->expr( mu, T );
        auto mesh = this->functionSpace()->mesh();
        return normL2( _range=this->functionSpace()->template rangeElements<0>(), _expr=eimexpr );
    }

    //Let geim the eim expansion of g
    //here is computed || g - geim ||_L2
    double diffL2Norm(  model_solution_type const& T , parameter_type const& mu , element_type const & eim_expansion ) const override
    {
        auto eimexpr = this->expr( mu, T );
        auto mesh = this->modelFunctionSpace()->mesh();
        auto difference = eimexpr - idv(eim_expansion);
        return normL2( _range=this->functionSpace()->template rangeElements<0>(), _expr=difference );
    }


    //Let \pi_g the projection of g on the function space
    //here is computed || \pi_g ||_L2
    double projExpressionL2Norm( model_solution_type const& T , parameter_type const& mu ) const override
    {
        auto eimexpr = this->expr( mu, T );
        auto mesh = this->functionSpace()->mesh();
        auto pi_g = vf::project( _space=this->functionSpace(), _range=this->functionSpace()->template rangeElements<0>(), _expr=eimexpr );
        return normL2( _range=this->functionSpace()->template rangeElements<0>(), _expr=idv(pi_g) );
    }

    //here is computed || \pi_g - geim ||_L2
    double projDiffL2Norm( model_solution_type const& T , parameter_type const& mu , element_type const& eim_expansion ) const override
    {
        auto eimexpr = this->expr( mu, T );
        auto mesh = this->functionSpace()->mesh();
        auto pi_g = vf::project( _space=this->functionSpace(), _range=this->functionSpace()->template rangeElements<0>(), _expr=eimexpr );
        auto diff = pi_g - eim_expansion;
        return normL2( _range=this->functionSpace()->template rangeElements<0>(), _expr=idv(diff) );
    }

    //here is computed || \pi_g - geim ||_Linf
    double projDiffLinfNorm( model_solution_type const& T , parameter_type const& mu , element_type const& eim_expansion ) const override
    {
        auto eimexpr = this->expr( mu, T );
        auto mesh = this->modelFunctionSpace()->mesh();
        auto pi_g = vf::project( _space=this->functionSpace(), _range=this->functionSpace()->template rangeElements<0>(), _expr=eimexpr );
        auto diff = pi_g - eim_expansion;
        auto linf = normLinf( _range=this->functionSpace()->template rangeElements<0>(), _pset=_Q<5>(), _expr=idv(diff) );
        return linf.template get<0>();
    }

    //here is computed || \pi_g  ||_Linf
    double projExprLinfNorm( model_solution_type const& T , parameter_type const& mu ) const override
    {
        auto eimexpr = this->expr( mu, T );
        auto mesh = this->functionSpace()->mesh();
        auto pi_g = vf::project( _space=this->functionSpace(), _range=this->functionSpace()->template rangeElements<0>(), _expr=eimexpr );
        auto linf = normLinf( _range=this->functionSpace()->template rangeElements<0>(), _pset=_Q<5>(), _expr=idv(pi_g) );
        return linf.template get<0>();
    }

    double interpolationError(model_solution_type const& T , parameter_type const& mu ) const override
    {
        auto eimexpr = this->expr( mu, T );
        auto mesh = this->modelFunctionSpace()->mesh();
        auto pi_g = vf::project( _space=this->functionSpace(), _range=this->functionSpace()->template rangeElements<0>(), _expr=eimexpr );
        auto difference = eimexpr - idv(pi_g);
        return normL2( _range=this->functionSpace()->template rangeElements<0>(), _expr=difference );
    }

    void computationalTimeStatistics( std::string appname ) override
        {
            computationalTimeStatistics( appname, typename boost::is_base_of<ModelCrbBaseBase,model_type>::type() );
        }
    void computationalTimeStatistics( std::string appname, boost::mpl::bool_<false> )
        {}
    void computationalTimeStatistics( std::string appname, boost::mpl::bool_<true> )
    {
#if 0 // this function is removed for now since it need M_crb

        //auto crbmodel = crbmodel_ptrtype( new crbmodel_type( M_model , CRBModelMode::CRB ) );
        if( !this->modelBuilt() )
            M_crbmodel = crbmodel_ptrtype( new crbmodel_type( this->model(), crb::stage::offline/*M_model*/ ) );
        //make sure that the CRB DB is already build
        if( !this->rbBuilt() )
            M_crb = crb_ptrtype( new crb_type( appname,
                                               M_crbmodel,
                                               crb::stage::offline ) );

        if ( !M_crb->isDBLoaded() || M_crb->rebuild() )
        {
            if( this->worldComm().isMasterRank() )
                LOG( INFO ) << "No CRB DB available, do crb offline computations...";
            M_crb->setOfflineStep( true );
            M_crb->offline();
        }

        int n_eval = ioption(_prefix=this->M_prefix,_name="eim.computational-time-neval");

        Eigen::Matrix<double, Eigen::Dynamic, 1> time_crb;
        Eigen::Matrix<double, Eigen::Dynamic, 1> time;
        time_crb.resize( n_eval );

        typename crb_type::sampling_ptrtype Sampling( new typename crb_type::sampling_type( this->model()->parameterSpace() ) );
        Sampling->logEquidistribute( n_eval  );

        //dimension
        int N =  ioption(_prefix=this->M_prefixModel,_name="crb.dimension");
        //reduced basis approximation space
        auto WN = M_crb->wn();
        int mu_number = 0;
        for( auto mu : *Sampling )
        {
            //LOG( INFO ) << "[computational] mu = \n"<<mu;

            boost::mpi::timer tcrb;
            auto o = M_crb->run( mu, time, doption(_prefix=this->M_prefixModel,_name="crb.online-tolerance") , N);
            auto solutions=o.template get<2>();
            auto uN = solutions.template get<0>();//vector of solutions ( one solution at each time step )

            int size=uN.size();
            auto u_crb = M_crb->expansion( uN[size-1] , N , false );

            boost::mpi::timer teim;
            this->beta( mu , u_crb );
            this->addOnlineTime( teim.elapsed() );
            time_crb( mu_number ) = tcrb.elapsed() ;
            mu_number++;
        }

        this->model()->computeStatistics( super::onlineTime() , super::name() );
        this->model()->computeStatistics( time_crb , super::name()+" - global crb timing" );
#endif
    }

    bool offlineStep() const override { return M_eim->offlineStep(); }
    bool adaptationSER() const override { return M_eim->adaptationSER(); }
    void setAdaptationSER(bool b) override {M_eim->setAdaptationSER(b);}
    bool rbCorrection() const override { return M_eim->rbCorrection(); }
    void setRbCorrection(bool b) override {M_eim->setRbCorrection(b);}

    void offline() override { M_eim->offline(); }
    void setRestart(bool b) override { M_eim->setRestart(b);}

    void setRB( crb_ptrtype rb ) override
    {
        M_crb = rb;
        if ( M_crb )
            M_crb_built=true;
    }
    bool rbBuilt() const override { return M_crb_built; }

    void setTrainSet( sampling_ptrtype tset ) override { M_eim->setTrainSet( tset ); }
    element_type interpolant( parameter_type const& mu ) override
    {
        auto beta = this->beta( mu );
        return expansion( M_eimFeSpaceDb.q(), this->beta( mu ) , M_M_max);
    }
    element_type interpolant( parameter_type const& mu , model_solution_type const & solution , int M) override
    {
        return expansion( M_eimFeSpaceDb.q(), this->beta( mu , solution , M) , M );
    }

    //return M_eim->operator()( mu , M_eim->mMax() ); }

    //element_type const& q( int m ) const { return M_eim->q( m ); }
    element_type const& q( int m ) const override
    {
        return M_eimFeSpaceDb.q( m );
    }

    void printInterpolationPointsSelection() const override
    {
        if ( this->worldComm().isMasterRank() )
        {
            std::ofstream file;
            std::string filename = "InterpolationPointsSelection-"+super::name()+".dat";
            file.open(filename);
            if( M_t[0].size() == 1 )
                file << "NbBasis" << "\t " << "x\n";
            if( M_t[0].size() == 2 )
                file << "NbBasis" << "\t " << "x\t y\n";
            if( M_t[0].size() == 3 )
                file << "NbBasis" << "\t " << "x\t y\t z\n";
            int size=M_t.size();
            for(int i=0 ; i<size; i++)
            {
                file << i <<"\t ";
                node_type t = M_t[i];
                int sizet = t.size();
                for( int j=0; j<sizet-1; j++)
                    file<< t(j)<<"\t ";
                file<< t(sizet-1) << "\n";
            }
        }
    }

    void printMuSelection() const override
    {
        if ( this->worldComm().isMasterRank() )
        {
            std::ofstream file;
            std::string filename = "MuSelection-"+super::name()+".dat";
            file.open(filename);
            int number=0;
            for( auto mu : *M_mu_sampling )
            {
                int sizemu=mu.size();
                if( number == 0 )
                {
                    file << "NbBasis\t ";

                    for(int i=0; i<sizemu-1; i++)
                        file<<(boost::format("Comp%1%") %i ).str()<<"\t ";
                    file<<(boost::format("Comp%1%") %(sizemu-1) ).str()<<"\n ";
                }
                file<< number <<"\t ";
                for(int i=0; i<sizemu-1; i++)
                    file << mu(i)<<"\t ";
                file<< mu(sizemu-1) << "\n" ;
                number++;
            }
            file.close();
        }
    }

    void printOfflineError() const override
    {
        if ( this->worldComm().isMasterRank() )
        {
            std::ofstream file;
            std::string filename = "OfflineError-"+super::name()+".dat";
            file.open(filename);
            file << "NbBasis" << "\t " << "Error\n";
            int size=M_offline_error.size();
            for(int i=0 ; i<size; i++)
            {
                file << i+1 <<"\t ";
                file<< M_offline_error[i]<<"\n ";
            }
            file.close();
        }
    }

    void printRbIterationsSER( int M ) const override
    {
        if ( this->worldComm().isMasterRank() )
        {
            std::ofstream file;
            std::string filename = "EIM-" + super::name() + "-rb_online_greedy_summary.dat";
            file.open(filename, std::ios::out | std::ios::app);

            double rb_online_mean_iterations = 0;
            int rb_online_min_terations = 0;
            int rb_online_max_terations = 0;
            double rb_online_max_increments = 0;
            if( M_rb_online_iterations.size() != 0 )
            {
                rb_online_min_terations = *std::min_element( M_rb_online_iterations.begin(), M_rb_online_iterations.end() );
                rb_online_max_terations = *std::max_element( M_rb_online_iterations.begin(), M_rb_online_iterations.end() );
                rb_online_mean_iterations = std::accumulate( M_rb_online_iterations.begin(), M_rb_online_iterations.end(), 0.0 );
                rb_online_mean_iterations /= M_rb_online_iterations.size();
                rb_online_max_increments = *std::max_element( M_rb_online_increments.begin(), M_rb_online_increments.end() );
            }

            std::string name = fs::current_path().string() + "/" + filename;
            fs::path pathfile( name.c_str() );
            if( fs::is_empty(pathfile) )
                file << "NbBasis" << "\t " << "Min_iter" << "\t" << "Mean_iter_sup" << "\t" << "Max_iter" << "\t" << "Max_increment \n";
            file << M << "\t"
                 << rb_online_min_terations << "\t"
                 << std::ceil(rb_online_mean_iterations) << "\t"
                 << rb_online_max_terations << "\t"
                 << rb_online_max_increments << "\n";
            file.close();
        }
    }

    std::vector<element_type> const& q() const override { return M_eimFeSpaceDb.q(); }


    void studyConvergence( parameter_type const & mu , model_solution_type & solution, std::vector< std::string > all_file_name ) const override { return M_eim->studyConvergence( mu , solution , all_file_name ) ; }
    boost::tuple<double,element_type> interpolationErrorEstimation( parameter_type const & mu , model_solution_type const& solution, int M ) const override { return M_eim->interpolationErrorEstimation(mu , solution, M) ; }
    double errorEstimationLinf( parameter_type const & mu, model_solution_type const& solution , int M ) const override { return M_eim->errorEstimationLinf(mu, solution, M) ; }
    //size_type mMax() const { return M_eim->mMax(); }
    size_type mMax( bool & error) const override
    {
        int max=0;
        int user_max = ioption(_prefix=this->M_prefix,_name="eim.dimension-max");
        // int built = M_max_q;
        //if the user wants to enrich the database we return M_M_max or
        //if the eim expansion contains less terms that expected by the user then there is no error.
        //But if there is already enough basis functions then we return M_M_max-1
        //to deal with error estimation
        if( (user_max+1) > M_eimFeSpaceDb.maxQ() )
        {
            max = M_M_max;
            //in that case if the eim expansion is finished there is no error associated
            //or if the user wants to enrich the DB it is to soon to known if it will be error
            error=false;
        }
        else
        {
            max = M_M_max-1;
            error=true;
        }
        return max;
    }
    size_type mMax() const override
    {
        int max=0;
        int user_max = ioption(_prefix=this->M_prefix,_name="eim.dimension-max");
        // int built = M_max_q;
        //if the user wants to enrich the database we return M_M_max or
        //if the eim expansion contains less terms that expected by the user then there is no error.
        //But if there is already enough basis functions then we return M_M_max-1
        //to deal with error estimation
        if( (user_max+1) > M_eimFeSpaceDb.maxQ() )
        {
            max = M_M_max;
        }
        else
        {
            max = M_M_max-1;
        }
        return max;
    }

    void setMax(int m, int max_q, int max_g, int max_z, int max_solution ) override
    {
        M_M_max = m;
        M_eimFeSpaceDb.setMax( max_q, max_g, max_z, max_solution );
    }

    int maxQ() const override
    {
        return M_eimFeSpaceDb.maxQ();
    }
    int maxG() const override
    {
        return M_eimFeSpaceDb.maxG();
    }
    int maxZ() const override
    {
        return M_eimFeSpaceDb.maxZ();
    }
    int maxSolution() const override
    {
        return M_eimFeSpaceDb.maxSolution();
    }

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & __ar, const unsigned int __version )
    {
        __ar & BOOST_SERIALIZATION_NVP( M_M_max );
        DVLOG(2) << "M saved/loaded\n";

        // save t
        __ar & BOOST_SERIALIZATION_NVP( M_t );
        DVLOG(2) << "t saved/loaded\n";

        __ar & BOOST_SERIALIZATION_NVP( M_offline_error );

        __ar & BOOST_SERIALIZATION_NVP( M_mu_sampling );

        // save B
        __ar & BOOST_SERIALIZATION_NVP( M_B );
        DVLOG(2) << "B saved/loaded\n";

        // load/save geometricspace context
        if ( Archive::is_loading::value )
        {
            bool hasCtxGeoEim = false;
            __ar & BOOST_SERIALIZATION_NVP( hasCtxGeoEim );
            // std::cout << "EIM load geoctx " << hasCtxGeoEim <<"\n";
            if ( hasCtxGeoEim )
            {
                auto geospace = std::make_shared<geometricspace_type>( this->worldCommPtr() );
                if ( this->functionSpace() && this->functionSpace()->mesh() )
                    geospace->setMesh( this->functionSpace()->mesh() );
                M_ctxGeoEim = std::make_shared<geometricspace_context_type>( geospace );
                __ar & BOOST_SERIALIZATION_NVP( *M_ctxGeoEim );
            }
        }
        else
        {
            bool hasCtxGeoEim = ( M_ctxGeoEim )? true : false;
            // std::cout << "EIM save geoctx " << hasCtxGeoEim <<"\n";
            __ar & BOOST_SERIALIZATION_NVP( hasCtxGeoEim );
            if ( hasCtxGeoEim )
                __ar & BOOST_SERIALIZATION_NVP( *M_ctxGeoEim );
        }
    }

private:
    EIMFunctionFeSpaceDb<functionspace_type,model_type> M_eimFeSpaceDb;
    // model_ptrtype M_model;
    model_weakptrtype M_model;
    std::string M_prefixModel;
    model_ptrtype M_modelAttached;
    expr_type M_expr;
    // model_solution_type& M_u;
    model_element_expr_type * M_u;
    model_element2_expr_type * M_u2;
    parameter_type& M_mu;
    crb_ptrtype M_crb;
    bool M_crb_built;
    sampling_ptrtype M_mu_sampling ;
    std::vector< double > M_evaluation_vector; // we need to store this to be able to compute expression of eim basis functions and the user can use them in integrals
    int M_M_max;
    std::vector<node_type> M_t;
    geometricspace_context_ptrtype M_ctxGeoEim;
    std::shared_ptr<context_type> M_ctxFeBasisEim;
    std::shared_ptr<model_element_expr_context_type> M_ctxFeModelSolution;
    rbfunctionspace_context_ptrtype M_ctxRbModelSolution;
    rbfunctionspace_context2_ptrtype M_ctxRbModelSolution2;
    matrix_type M_B;
    std::vector<double> M_offline_error;
    eim_ptrtype M_eim;

    bool  M_write_nl_solutions;
    std::string M_write_nl_directory;

    std::vector<int> M_rb_online_iterations;
    std::vector<double> M_rb_online_increments;

    std::shared_ptr<element_type> M_internalModelFeFunc;
};

namespace detail
{

template<typename ArgModelType,typename ArgExprType,typename ArgSpaceType,typename ArgElementType,typename ArgElement2Type>
struct compute_eim_return
{
    // typedef typename boost::remove_reference<typename boost::remove_pointer<typename parameter::binding<Args, tag::model>::type>::type>::type::element_type model1_type;
    // typedef typename boost::remove_const<typename boost::remove_pointer<model1_type>::type>::type model_type;
    // typedef typename boost::remove_reference<typename parameter::binding<Args, tag::expr>::type>::type expr_type;
    // typedef typename boost::remove_reference<typename parameter::binding<Args, tag::space>::type>::type::element_type space_type;
    // typedef typename boost::remove_reference<typename parameter::binding<Args, tag::element>::type>::type element_type;
    // typedef typename boost::remove_reference<typename parameter::binding<Args, tag::element2,element_type>::type>::type element2_type;
    using model_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<ArgModelType>>>;
    using expr_type = std::decay_t<ArgExprType>;
    using space_type = std::decay_t<ArgSpaceType>;
    using element_type = std::decay_t<ArgElementType>;
    using element2_type = std::decay_t<ArgElement2Type>;
    static const int subspaceid = mpl::if_< mpl::bool_< model_type::functionspace_type::is_composite>,
                                            mpl::int_<element_type::functionspace_type::basis_type::TAG>,// must be improve! loop on subspaces and detect the same
                                            mpl::int_<-1> >::type::value;
    static const int subspaceid2 = mpl::if_< mpl::bool_< model_type::functionspace_type::is_composite>,
                                            mpl::int_<element2_type::functionspace_type::basis_type::TAG>,// must be improve! loop on subspaces and detect the same
                                            mpl::int_<-1> >::type::value;

    typedef EIMFunction<model_type, space_type, expr_type, subspaceid, subspaceid2 > type;
    typedef std::shared_ptr<type> ptrtype;
};
}

#if 0
BOOST_PARAMETER_FUNCTION(
    ( typename Feel::detail::compute_eim_return<Args>::ptrtype ), // 1. return type
    eim,                        // 2. name of the function template
    tag,                                        // 3. namespace of tag types
    ( required
      ( in_out(model),          * )
      ( in_out(element),        * )
      ( in_out(parameter),        * )
      ( in_out(expr),          * )
      ( name, * )
      ( space, *)
        ) // required
    ( optional
      ( in_out(element2),        *, element )
      //( space, *( boost::is_convertible<mpl::_,std::shared_ptr<FunctionSpaceBase> > ), model->functionSpace() )
      //( space, *, model->functionSpace() )
      ( sampling, *, model->parameterSpace()->sampling() )
      ( verbose, (int), 0 )
      ( filename, *( boost::is_convertible<mpl::_,std::string> ), "" )
      ( directory, *( boost::is_convertible<mpl::_,std::string> ), "" )
      ( prefix, (std::string), "")
        ) // optionnal
)
#endif
template <typename ... Ts>
auto eim( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && model = args.get(_model);
    auto && element = args.get(_element);
    auto && parameter = args.get(_parameter);
    auto && expr = args.get(_expr);
    std::string const& name = args.get(_name);
    auto && space = args.get(_space);
    auto && element2 = args.get_else(_element2,element);
    auto && sampling = args.get_else(_sampling,model->parameterSpace()->sampling() );
    std::string const& filename = args.get_else(_filename,"");
    std::string const& directory = args.get_else(_directory,"");
    std::string const& prefix = args.get_else(_prefix,"");

    using eim_helper_args_type =  Feel::detail::compute_eim_return<decltype(model),decltype(expr),decltype(space),decltype(element),decltype(element2)>;
    using eim_type = typename eim_helper_args_type::type;
    using eim_ptrtype = std::shared_ptr<eim_type>;

    auto eimFunc = std::make_shared<eim_type>( model, space, element, element2, parameter, expr, sampling, name, filename, directory, prefix );
    if ( eim_type::model_use_nosolve )
        eimFunc->attachModel( model );
    return eimFunc;
}


template<typename ModelType>
struct EimFunctionNoSolve : public EimFunctionNoSolveBase
{
    typedef typename ModelType::functionspace_type functionspace_type;
    typedef typename ModelType::functionspace_ptrtype functionspace_ptrtype;
    typedef typename functionspace_type::element_type element_type;
    typedef typename functionspace_type::element_ptrtype element_ptrtype;
    typedef typename ModelType::parameterspace_type parameterspace_type;
    typedef typename ModelType::parameterspace_ptrtype parameterspace_ptrtype;
    typedef typename parameterspace_type::element_type parameter_type;
    typedef typename parameterspace_type::sampling_type sampling_type;
    typedef typename parameterspace_type::sampling_ptrtype sampling_ptrtype;
    typedef typename functionspace_type::value_type value_type;

    // reduced basis space
    //typedef typename ModelType::rbfunctionspace_type rbfunctionspace_type;
    //typedef std::shared_ptr<rbfunctionspace_type> rbfunctionspace_ptrtype;

    typedef std::shared_ptr<ModelType> model_ptrtype;
    typedef std::weak_ptr<ModelType> model_weakptrtype;

    static const int nb_spaces = functionspace_type::nSpaces;
    typedef typename mpl::if_< boost::is_same< mpl::int_<nb_spaces> , mpl::int_<2> > , fusion::vector< mpl::int_<0>, mpl::int_<1> >  ,
                       typename mpl::if_ < boost::is_same< mpl::int_<nb_spaces> , mpl::int_<3> > , fusion::vector < mpl::int_<0> , mpl::int_<1> , mpl::int_<2> > ,
                                  typename mpl::if_< boost::is_same< mpl::int_<nb_spaces> , mpl::int_<4> >, fusion::vector< mpl::int_<0>, mpl::int_<1>, mpl::int_<2>, mpl::int_<3> >,
                                                     fusion::vector< mpl::int_<0>, mpl::int_<1>, mpl::int_<2>, mpl::int_<3>, mpl::int_<4> >
                                                     >::type >::type >::type index_vector_type;

    EimFunctionNoSolve( model_ptrtype const& model, std::string const& prefix = "" )
        :
        M_model( model ),
        M_prefix( prefix )
        {
            if ( model->functionSpace() )
            {
                M_elt = model->functionSpace()->elementPtr();
                M_elt->setConstant( boost::lexical_cast<value_type>("inf") );
            }
        }
    element_type const& solve( parameter_type const& mu ) const
    {
        DVLOG(2) << "no solve required\n";
        // static const bool is_composite = functionspace_type::is_composite;
        // return solve( mu , mpl::bool_< is_composite >() );
        CHECK( M_elt ) << " element not init";
        return *M_elt;
    }
#if 0
    element_type solve( parameter_type const& mu , mpl::bool_<false> )
    {
        //value_type x = boost::lexical_cast<value_type>("inf");
        //M_elt = vf::project( _space=M_model->functionSpace(), _expr=cst(x) );
        //M_elt.setConstant( x );
        return M_elt;
    }
    element_type solve( parameter_type const& mu , mpl::bool_<true> )
    {
#if 0
        ProjectInfCompositeCase project_inf_composite_case( M_elt );
        index_vector_type index_vector;
        fusion::for_each( index_vector, project_inf_composite_case );
        return project_inf_composite_case.element();
#else
        return M_elt;
#endif
    }
#endif

    std::string /*const&*/ modelName() const { return M_model.lock()->modelName(); }
    std::string prefix() const { return M_prefix; }
    uuids::uuid uuid() const { return M_model.lock()->uuid(); }
    functionspace_ptrtype const& functionSpace() const { return M_model.lock()->functionSpace(); }
    parameterspace_ptrtype const& parameterSpace() const { return M_model.lock()->parameterSpace(); }
#if 0
    struct ProjectInfCompositeCase
    {
        ProjectInfCompositeCase( element_type & composite_element)
            :
            M_element( composite_element )
        {}

        template< typename T >
        void
        operator()( const T& t ) const
        {
            auto view = M_element.template element< T::value >();
            view.setConstant( boost::lexical_cast<value_type>("inf") );
            //auto space = view.functionSpace();
            //view = vf::project( _space=space, _expr=cst( boost::lexical_cast<value_type>("inf") ) );
        }

        element_type element()
        {
            return M_element;
        }

        element_type M_element;

    }; //struct ProjectInfOnSubspace
#endif
    // model_ptrtype M_model;
    model_weakptrtype M_model;
    std::string M_prefix;
    element_ptrtype M_elt;
};

template<typename ModelType>
std::shared_ptr<EimFunctionNoSolve<ModelType>>
eim_no_solve( std::shared_ptr<ModelType> model )
{
    return std::shared_ptr<EimFunctionNoSolve<ModelType>>( new EimFunctionNoSolve<ModelType>( model ) );
}

template< typename ExprType , typename EimType >
Expr< vf_div< Expr< vf_sub< ExprType, Expr< OpId< typename EimType::element_type , __VALUE> > > > ,  Cst<double>  > >
eimBasisExpression(int m, ExprType const& expr, EimType const& eim)
{
    auto any_type = eim->buildBasisFunction( m );
    typedef typename EimType::element_type element_type;

    typedef
        Expr< vf_div <
            Expr< vf_sub< ExprType , Expr< OpId< element_type , __VALUE> >  > > ,
            Cst< double >
        > >basis_type;

    if (!boost::any_cast<basis_type>(&any_type))
    {
        throw std::logic_error( "[qExpression] not cast possible for eim basis function expression" );
    }
    return boost::any_cast<basis_type>(any_type);
}

}
#endif /* _FEELPP_EIM_HPP */
