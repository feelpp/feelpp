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
#include <boost/math/special_functions/nonfinite_num_facets.hpp>
#endif
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/base_object.hpp>

#include <feel/feelcrb/crbdb.hpp>
#include <feel/feelcrb/parameterspace.hpp>

#include <feel/feelvf/vf.hpp>

#include <feel/feelcrb/crb.hpp>
#include <feel/feelcrb/crbmodel.hpp>

#include <Eigen/Core>

namespace Feel
{
class ModelCrbBaseBase {};

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
class EIM : public CRBDB
{
    typedef  CRBDB super;

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
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef typename functionspace_type::element_type element_type;
    typedef typename functionspace_type::element_ptrtype element_ptrtype;

    typedef typename ModelType::solution_type solution_type;

    typedef typename ModelType::parameterspace_type parameterspace_type;
    typedef typename ModelType::parameter_type parameter_type;
    typedef typename parameterspace_type::sampling_ptrtype sampling_ptrtype;

    typedef boost::tuple<double,Eigen::Matrix<double,nDim,1> > space_residual_type;
    typedef boost::tuple<double,parameter_type> parameter_residual_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    EIM()
        :
        super(),
        M_vm(),
        M_is_read( false ),
        M_is_written( false ),
        M_name( "default" ),
        M_M( 1 ),
        M_M_max( 1 ),
        M_WN(),
        M_offline_done( false ),
        M_tol( 1e-8 ),
        M_q(),
        M_B(),
        M_t(),
        M_index_max(),
        M_model( 0 )
        {}
    EIM( po::variables_map const& vm, model_type* model, sampling_ptrtype sampling, double __tol = 1e-8 )
        :
        super(model->modelName(), model->name(), model->name(), vm ),
        M_vm( vm ),
        M_is_read( false ),
        M_is_written( false ),
        M_name( model->name() ),
        M_trainset( sampling ),
        M_M( 1 ),
        M_M_max( 1 ),
        M_WN( M_vm["eim.dimension-max"].template as<int>() ),
        M_offline_done( false ),
        M_tol( __tol ),
        M_q(),
        M_B(),
        M_t(),
        M_index_max(),
        M_model( model )
        {
            if ( !loadDB() || M_vm["eim.rebuild-database"].template as<bool>() )
            {
                LOG(INFO) << "construct EIM approximation...\n";
                M_offline_done = false;
            }

            if ( !this->isOfflineDone() )
                offline();
#if 0
            //old version of convergence
            if( M_vm["eim.study-cvg"].template as<bool>() )
            {
                M_WN=1;
                do{
                    offline();
                    studyConvergence();

                    M_offline_done = false;
                    M_WN++;
                }while( M_WN < M_vm["eim.dimension-max"].template as<int>() );
            }
#endif

        }

    EIM( EIM const & __bbf )
        :
    super(__bbf),
        M_is_read( __bbf.M_is_read ),
        M_is_written( __bbf.M_is_written ),
        M_name( __bbf.M_name ),
        M_M( __bbf.M_M ),
        M_M_max (__bbf.M_M_max ),
        M_WN(__bbf.M_WN ),
        M_offline_done( __bbf.M_offline_done ),
        M_tol( __bbf.M_tol ),
        M_q( __bbf.M_q ),
        M_B( __bbf.M_B ),
        M_t( __bbf.M_t ),
        M_index_max( __bbf.M_index_max ),
        M_model( __bbf.M_model )
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

    /**
       \return the number of DOF in space discretization
    */
    size_type nDOF() const {  FEELPP_ASSERT( M_model != 0 ).error( "Invalid EIM model" ); return M_model->functionSpace()->nLocalDof(); }

    /**
     * return the set of reduced basis functions associated with the eim
     */
    std::vector<element_type> const& q() const {  return M_q; }

    /**
     * return the m-th reduced basis function associated with the eim
     */
    element_type const&
    q( size_type __m ) const
        {
            FEELPP_ASSERT( __m >= 0 && __m < M_M )( __m )( M_M ).error( "out of bounds access" );

            return M_q[ __m ];
        }

    size_type mMax() const { return M_M_max; }


    /**
       Check whether the offline stage has been executed and the database created.
       \return \c true if the offline has been executed and the DB saved, \c false otherwise
    */
    bool isOfflineDone() const
        {
            return M_offline_done;
        }

    //@}




    /** @name  Mutators
     */
    //@{

    void setTrainSet( sampling_ptrtype pset ) { M_trainset = pset; }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * save the CRB database
     */
    void saveDB()
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

    /**
     * load the CRB database
     */
    bool loadDB()
        {
            //if ( this->rebuildDB() )
            //return false;

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
    vector_type beta( parameter_type const& mu  ) const { return beta( mu, this->mMax() ); }
    vector_type beta( parameter_type const& mu, solution_type const& T  ) const { return beta( mu, T, this->mMax() ); }
    vector_type beta( parameter_type const& mu, size_type M  ) const;
    vector_type beta( parameter_type const& mu, solution_type const& T, size_type M  ) const;

    std::vector<double> studyConvergence( parameter_type const & mu, solution_type & solution ) const;

    void computationalTimeStatistics( std::string appname )  { return M_model->computationalTimeStatistics(); }
    element_type residual ( size_type M ) const;

    parameter_residual_type computeBestFit( sampling_ptrtype trainset, int __M );

    element_type operator()( parameter_type const& mu , int N) const { return expansion( M_q, beta( mu ) , N); }
    element_type operator()( parameter_type const& mu, solution_type const& T , int N ) const { return expansion( M_q, beta( mu, T ) , N ); }
    /**
       orthonormalize
    */
    void orthonormalize( std::vector<element_type>& );

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & __ar, const unsigned int __version )
        {
            LOG(INFO) << "serializing...\n";

            __ar & BOOST_SERIALIZATION_NVP( M_name );

            LOG(INFO) << "name saved/loaded...\n";

            __ar & BOOST_SERIALIZATION_NVP( M_offline_done );
            LOG(INFO) << "offline status...\n";

            __ar & BOOST_SERIALIZATION_NVP( M_M_max );
            LOG(INFO) << "M saved/loaded\n";
            M_M = M_M_max;

            // save index
            __ar & BOOST_SERIALIZATION_NVP( M_index_max );
            LOG(INFO) << "index saved/loaded\n";

            // save t
            __ar & BOOST_SERIALIZATION_NVP( M_t );
            LOG(INFO) << "t saved/loaded\n";


            if ( Archive::is_loading::value )
            {
                for( int i = 0; i < M_M_max; ++ i )
                {
                    M_q.push_back( M_model->functionSpace()->element() );
                }
                for( int i = 0; i < M_M_max; ++ i )
                    __ar & BOOST_SERIALIZATION_NVP( M_q[i] );
                // save q
                LOG(INFO) << "q saved/loaded\n";
            }
            else
            {
                for( int i = 0; i < M_M_max; ++ i )
                    __ar & BOOST_SERIALIZATION_NVP( M_q[i] );
            }


            // save B
            __ar & BOOST_SERIALIZATION_NVP( M_B );
            LOG(INFO) << "B saved/loaded\n";

        }
    //@}


protected:

    po::variables_map M_vm;
    mutable bool M_is_read;
    mutable bool M_is_written;

    std::string M_name;
    sampling_ptrtype M_trainset;
    size_type M_M;
    size_type M_M_max;

    size_type M_WN;

    mutable bool M_offline_done;

    double M_tol;

    std::vector<element_type> M_g;

    std::vector<element_type> M_q;

    matrix_type M_B;

    std::vector<node_type> M_t;

    std::vector<size_type> M_index_max;

    model_type* M_model;
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
    void offline();
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<typename ModelType>
typename EIM<ModelType>::vector_type
EIM<ModelType>::beta( parameter_type const& mu, size_type __M ) const
{
    // beta=B_M\g(Od(indx),mut(i))'
    vector_type __beta( __M );
    for ( size_type __m = 0;__m < __M;++__m )
    {
        __beta[__m] = M_model->operator()( this->M_t[__m], mu );
    }
    this->M_B.block(0,0,__M,__M).template triangularView<Eigen::UnitLower>().solveInPlace(__beta);

    return __beta;
}
template<typename ModelType>
typename EIM<ModelType>::vector_type
EIM<ModelType>::beta( parameter_type const& mu, solution_type const& T, size_type __M ) const
{
    // beta=B_M\g(Od(indx),mut(i))'
    vector_type __beta( __M );
    for ( size_type __m = 0;__m < __M;++__m )
    {
        __beta[__m] = M_model->operator()( T, this->M_t[__m], mu );
    }
    this->M_B.block(0,0,__M,__M).template triangularView<Eigen::UnitLower>().solveInPlace(__beta);

    return __beta;
}

template<typename ModelType>
typename EIM<ModelType>::element_type
EIM<ModelType>::residual( size_type __M ) const
{
    LOG(INFO) << "compute residual for m=" << __M << "...\n";
    vector_type rhs( __M );

    //LOG(INFO) << "g[" << __M << "]=" << M_g[__M] << "\n";
    for ( size_type __m = 0;__m < __M;++__m )
    {
        LOG(INFO) << "t[" << __m << "]=" << M_t[__m] << "\n";
        rhs[__m]= M_g[__M]( M_t[__m] )(0,0,0);
    }
    this->M_B.block(0,0,__M,__M).template triangularView<Eigen::UnitLower>().solveInPlace(rhs);
    LOG(INFO) << "solve B sol = rhs with rhs = " << rhs <<"\n";

    // res(:,i)=M_g(:,i)-q(:,0:i)*sigma
    LOG(INFO) << "compute residual..." <<"\n";
    using namespace vf;
    auto z = expansion( M_q, rhs, __M );
    LOG(INFO) << "return residual..." <<"\n";
    return vf::project( _space=M_model->functionSpace(),
                        _expr=idv(M_g[__M])-idv( z ) );
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
EIM<ModelType>::computeBestFit( sampling_ptrtype trainset, int __M )
{
    LOG(INFO) << "compute best fit  for m=" << __M
              << " and trainset of size " << trainset->size() << "...\n";
    using namespace vf;
    parameter_type mu = M_model->parameterSpace()->element();

    vector_type maxerr( trainset->size() );
    maxerr.setZero();
    int index = 0;
    LOG(INFO) << "Compute best fit M=" << __M << "\n";
    BOOST_FOREACH( mu, *trainset )
    {
        LOG(INFO) << "compute best fit check mu...\n";
        mu.check();
        LOG_EVERY_N(INFO, 1 ) << " (every 10 mu) compute fit at mu="<< mu <<"\n" ;
        // evaluate model at mu

        auto Z = M_model->operator()( mu );

        vector_type rhs( __M );
        for ( size_type __m = 0;__m < __M;++__m )
        {
            rhs[__m]= Z( M_t[__m] )(0,0,0);
        }
        this->M_B.block(0,0,__M,__M).template triangularView<Eigen::UnitLower>().solveInPlace(rhs);
        auto res = vf::project( _space=M_model->functionSpace(),
                                _expr=idv(Z)-idv( expansion( M_q, rhs, __M ) ) );

        std::string norm_used = option(_name="eim.norm-used-for-residual").template as<std::string>();
        bool check_name_norm = false;
        LOG( INFO ) << "[computeBestFit] norm used : "<<norm_used;
        if( norm_used == "Linfty" )
        {
            check_name_norm=true;
            auto resmax = normLinf( _range=elements(M_model->mesh()), _pset=_Q<5>(), _expr=idv(res) );
            LOG_ASSERT( index < trainset->size() ) << "Invalid index " << index << " should be less than trainset size = " << trainset->size() << "\n";
            maxerr( index++ ) = resmax.template get<0>();
        }
        if( norm_used == "L2" )
        {
            check_name_norm=true;
            double norm = math::sqrt( integrate( _range=elements( M_model->mesh() ) ,_expr=idv(res)*idv(res) ).evaluate()( 0,0 ) );
            LOG_ASSERT( index < trainset->size() ) << "Invalid index " << index << " should be less than trainset size = " << trainset->size() << "\n";
            maxerr( index++ ) = norm;
        }
        if( norm_used == "LinftyVec" )
        {
            check_name_norm=true;
            double norm = res.linftyNorm();
            LOG_ASSERT( index < trainset->size() ) << "Invalid index " << index << " should be less than trainset size = " << trainset->size() << "\n";
            maxerr( index++ ) = norm ;
        }
        CHECK( check_name_norm ) <<"[EIM] The name of the norm "<<norm_used<<" is not known\n";

        int index2;
        auto err = maxerr.array().abs().maxCoeff( &index2 );
        LOG_EVERY_N(INFO, 1 ) << " (every 10 mu) maxerr=" <<  err << " at index = " << index2 << " at mu = " << trainset->at(index2) << "\n";
    }
    LOG_ASSERT( index == trainset->size() ) << "Invalid index " << index << " should be equal to trainset size = " << trainset->size() << "\n";
    auto err = maxerr.array().abs().maxCoeff( &index );
    LOG(INFO)<< "err=" << err << " reached at index " << index << " and mu=" << trainset->at(index) << "\n";
    return boost::make_tuple( err, trainset->at(index) );
}
template<typename ModelType>
void
EIM<ModelType>::offline(  )
{
    using namespace vf;

    if ( this->isOfflineDone() )
        return;
    LOG(INFO) << "[offline] starting offline stage ...\n";
    M_M = 1;
    LOG(INFO) << "[offline] create mu_1...\n";
    // min element in Dmu to start with (in // each proc have the same element)
    auto mu = M_model->parameterSpace()->min();

    if ( !M_trainset )
        M_trainset = M_model->parameterSpace()->sampling();
    if ( M_trainset->empty() )
    {
        int sampling_size = M_vm["eim.sampling-size"].template as<int>();
        std::string file_name = ( boost::format("eim_trainset_%1%") % sampling_size ).str();
        std::ifstream file ( file_name );
        if( ! file )
        {
            M_trainset->randomize( sampling_size  );
            M_trainset->writeOnFile(file_name);
        }
        else
        {
            M_trainset->clear();
            M_trainset->readFromFile(file_name);
        }
    }

    // store residual
    auto res = M_model->functionSpace()->element();

    //if( !M_vm["eim.study-cvg"].template as<bool>() || M_WN == 1 )
    // {
    LOG(INFO) << "compute finite element solution at mu_1...\n";
    M_g.push_back( M_model->operator()( mu ) );

    LOG(INFO) << "compute T^" << 0 << "...\n";
    // Build T^0
    auto zmax = normLinf( _range=elements(M_model->mesh()), _pset=_Q<5>(), _expr=idv(M_g[0]) );
    // store space coordinate where max absolute value occurs
    M_t.push_back( zmax.template get<1>() );
    LOG(INFO) << "norm Linf = " << zmax.template get<0>() << " at " << zmax.template get<1>() << "\n";
    //LOG(INFO) << "g = " << M_g[0] << "\n";

    LOG(INFO) << "compute and insert q_0...\n";
    // insert first element
    auto q = M_g[0];
    //q.scale( 1./zmax.template get<0>() );
    q.scale( 1./ M_g[0]( M_t[0] )( 0, 0, 0 ) );
    M_q.push_back( q );

    LOG(INFO) << "compute entry (0,0) of interpolation matrix...\n";
    this->M_B.resize( 1, 1 );
    this->M_B( 0, 0 ) = 1;
    CHECK( math::abs( M_q[0]( M_t[0] )( 0, 0, 0 ) - 1 ) < 1e-10 )
        << "q[0](t[0] != 1 " << "q[0] = " << M_q[0]( M_t[0] )( 0, 0, 0 )
        << "  t[0] = "<< M_t[0] << "\n";

    ++M_M;
    //}
    //else
    //M_M = M_WN - 1;

    /**
       \par build \f$W^g_M\f$
    */
    double err = 1;

   //for residual storage
   //auto res = M_model->functionSpace()->element();

    LOG(INFO) << "start greedy algorithm...\n";
    //for(  ; M_M < M_WN; ++M_M ) //err >= this->M_tol )
    for(  ; M_M < M_vm["eim.dimension-max"].template as<int>(); ++M_M ) //err >= this->M_tol )
    {

        LOG(INFO) << "M=" << M_M << "...\n";

        LOG(INFO) << "compute best fit error...\n";
        // compute mu = arg max inf ||G(.;mu)-z||_infty
        auto bestfit = computeBestFit( M_trainset, this->M_M-1 );

        auto g_bestfit = M_model->operator()( bestfit.template get<1>() );
        auto gmax = normLinf( _range=elements(M_model->mesh()), _pset=_Q<5>(), _expr=idv(g_bestfit) );

        DVLOG(2) << "best fit max error = " << bestfit.template get<0>() << " relative error = " << bestfit.template get<0>()/gmax.template get<0>() << " at mu = "
                 << bestfit.template get<1>() << "  tolerance=" << M_vm["eim.error-max"].template as<double>() << "\n";

        //if we want to impose the use of dimension-max functions, we don't want to stop here
        if ( (bestfit.template get<0>()/gmax.template get<0>()) < M_vm["eim.error-max"].template as<double>() &&  ! M_vm["eim.use-dimension-max-functions"].template as<bool>() )
            break;

        /**
         * we have a new \f$\mu^g_M\f$, insert it in \f$S^g_M\f$ and the
         * corresponding \f$z \in W^g_M\f$
         */
        LOG(INFO) << "[offline] S(" << this->M_M-1 << ") = " << bestfit.template get<1>() << "\n";

        // update M_g(:,M-1)
        M_g.push_back( g_bestfit );

        //orthonormalize( M_g );

        // build T^m such that T^m-1 \subset T^m
        LOG(INFO) << "[offline] compute residual M="<< M_M << "..." <<"\n";
        res = this->residual(M_M-1);
        //LOG(INFO) << "residual = " << res << "\n";
        LOG(INFO) << "[offline] compute arg sup |residual|..." <<"\n";
        auto resmax = normLinf( _range=elements(M_model->mesh()), _pset=_Q<5>(), _expr=idv(res) );
        LOG(INFO) << "[offline] store coordinates where max absolute value is attained " << resmax.template get<1>() << "..." <<"\n";
        // store space coordinate where max absolute value occurs
        M_t.push_back( resmax.template get<1>() );

        LOG(INFO) << "[offline] scale new basis function by " << 1./resmax.template get<0>() << "..." <<"\n";
        res.scale( 1./res(resmax.template get<1>())(0,0,0) );
        LOG(INFO) << "store new basis function..." <<"\n";
        M_q.push_back( res );

        std::for_each( M_t.begin(), M_t.end(), []( node_type const& t ) { DVLOG(2) << "t=" << t << "\n"; } );
        // update interpolation matrix
        M_B.conservativeResize( M_M, M_M );
        for( int __i = 0; __i < M_M; ++__i )
        {
            for( int __j = 0; __j < M_M; ++__j )
            {
                this->M_B( __i, __j ) = M_q[__j]( M_t[__i] )(0,0,0);

            }
        }
        DVLOG(2) << "[offline] Interpolation matrix: M_B = " << this->M_B <<"\n";
#if 0
        for( int __i = 0; __i < M_M; ++__i )
        {
            for( int __j = __i+1; __j < M_M; ++__j )
            {
                LOG_ASSERT( math::abs( M_q[__j]( M_t[__i] )(0,0,0)) < 1e-10 )
                    << "q[j](t_i) when j > i should be 0, = "
                    << math::abs( M_q[__j]( M_t[__i] )(0,0,0)) << "\n";
            }
            //this->M_B(__i,__i)=1;
            //this->M_B(__i,__i)=;
            LOG_ASSERT( math::abs( M_q[__i]( M_t[__i] )( 0, 0, 0 ) - 1 ) < 1e-10 )
                << "q[ " << __i << "](t[" << __i << "] != 1 "
                << "t[" << __i << "] = "<< M_t[__i] << "\n";
        }
#endif


        VLOG(2) << "================================================================================\n";
        //if we want to impose the use of dimension-max functions, we don't want to stop here
        if ( resmax.template get<0>() < M_vm["eim.error-max"].template as<double>() &&  ! M_vm["eim.use-dimension-max-functions"].template as<bool>() )
        {
            ++M_M;
            break;
        }
    }

    LOG(INFO) << "[offline] M_max = " << M_M-1 << "...\n";
    this->M_M_max = this->M_M-1;

    this->M_offline_done = true;
    LOG(INFO) << "[offline] saving DB...\n";
    saveDB();
    LOG(INFO) << "[offline] done with offline stage ...\n";
}

template<typename ModelType>
std::vector<double>
EIM<ModelType>::studyConvergence( parameter_type const & mu , solution_type & solution) const
{
    LOG(INFO) << " Convergence study \n";
    int proc_number =  Environment::worldComm().globalRank();

    std::vector<double> l2ErrorVec(this->mMax(), 0.0);

    std::string mu_str;
    for ( int i=0; i<mu.size(); i++ )
        mu_str= mu_str + ( boost::format( "_%1%" ) %mu(i) ).str() ;

    std::string file_name = "cvg-eim-"+M_model->name()+"-"+mu_str+".dat";
    if( std::ifstream( file_name ) )
        std::remove( file_name.c_str() );

    std::ofstream conv;
    if( proc_number == Environment::worldComm().masterRank() )
        {
            conv.open(file_name, std::ios::app);
            conv << "#Nb_basis" << "\t" << "L2_error" << "\n";
        }

    bool use_expression = option(_name="eim.compute-error-with-truth-expression").template as<bool>();
    int max = this->mMax();
    for(int N=1; N<=max; N++)
    {

        double exprl2norm = 0 , diffl2norm = 0 ;

        if( use_expression )
        {
            exprl2norm =M_model->expressionL2Norm( solution , mu );
            auto eim_approximation = this->operator()(mu , solution, N);
            diffl2norm = M_model->diffL2Norm( solution , mu , eim_approximation );
        }
        else
        {
            exprl2norm =M_model->projExpressionL2Norm( solution , mu );
            auto eim_approximation = this->operator()(mu , solution, N);
            diffl2norm = M_model->projDiffL2Norm( solution , mu , eim_approximation );
        }

        double l2_error = diffl2norm / exprl2norm ;
        //interpolation error : || projection_g - g ||_L2
        double interpolation_error = M_model->interpolationError( solution , mu );
#if 0
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
#endif
        l2ErrorVec[N-1] = l2_error; // /!\ l2ErrorVec[i] represents error with i+1 bases

        if( proc_number == Environment::worldComm().masterRank() )
            conv << N << "\t" << l2_error << "\t" << interpolation_error << "\n";

    }//loop over basis functions

    conv.close();
    return l2ErrorVec;
}

template<typename SpaceType, typename ModelSpaceType, typename ParameterSpaceType>
class EIMFunctionBase
{
public:

    typedef SpaceType functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef typename functionspace_type::element_type element_type;
    typedef typename functionspace_type::element_ptrtype element_ptrtype;
    typedef typename functionspace_type::mesh_type mesh_type;
    typedef typename functionspace_type::mesh_ptrtype mesh_ptrtype;
    typedef typename functionspace_type::value_type value_type;

    typedef ModelSpaceType model_functionspace_type;
    typedef typename model_functionspace_type::element_type solution_type;

    static const uint16_type nDim = mesh_type::nDim;

    typedef ParameterSpaceType parameterspace_type;
    typedef boost::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef typename parameterspace_type::element_type parameter_type;
    typedef typename parameterspace_type::sampling_ptrtype sampling_ptrtype;
    typedef Eigen::Matrix<double, nDim, 1> node_type;

    typedef EIM<EIMFunctionBase<SpaceType,model_functionspace_type, ParameterSpaceType> > eim_type;
    typedef typename eim_type::vector_type vector_type;

    typedef boost::shared_ptr<eim_type> eim_ptrtype;
    //typedef typename eim_type::betam_type betam_type;
    //typedef typename eim_type::qm_type qm_type;


    EIMFunctionBase( po::variables_map const& vm,
                     functionspace_ptrtype fspace,
                     parameterspace_ptrtype pspace,
                     sampling_ptrtype sampling,
                     std::string const& modelname,
                     std::string const& name )
        :
        M_vm( vm ),
        M_fspace( fspace ),
        M_pspace( pspace ),
        M_trainset( sampling ),
        M_modelname( modelname ),
        M_name( name )
        {
            LOG(INFO)<< "EimFunctionBase constructor\n";
        }
    virtual ~EIMFunctionBase()
        {}
    std::string const& name() const { return M_name; }
    void setName( std::string const& name ) { M_name = name; }
    std::string modelName() const { return M_modelname; }
    void setModelName( std::string const& name ) { M_modelname = name; }

    functionspace_ptrtype functionSpace() const { return M_fspace; }
    functionspace_ptrtype functionSpace()  { return M_fspace; }
    parameterspace_ptrtype parameterSpace() const { return M_pspace; }
    parameterspace_ptrtype parameterSpace()  { return M_pspace; }
    sampling_ptrtype trainSet() { return M_trainset; }
    sampling_ptrtype trainSet() const { return M_trainset; }
    virtual void setTrainSet( sampling_ptrtype tset ) { M_trainset = tset; }

    void addOnlineTime( const double time )
    {
        int size = M_online_time.size();
        M_online_time.conservativeResize( size+1 );
        M_online_time( size ) = time;
    }
    Eigen::VectorXd onlineTime() const { return M_online_time; }

    mesh_ptrtype mesh() const { return M_fspace->mesh(); }
    mesh_ptrtype mesh()  { return M_fspace->mesh(); }

    virtual element_type operator()( parameter_type const& ) = 0;
    virtual element_type operator()( solution_type const& T, parameter_type const& ) = 0;
    virtual element_type interpolant( parameter_type const& ) = 0;
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

    value_type operator()( solution_type const& T, node_type const& x, parameter_type const& mu )
        {
            VLOG(2) << "calling EIMFunctionBase::operator()( x=" << x << ", mu=" << mu << ")\n";
            element_type v = this->operator()( T, mu );
            value_type res = v(x)(0,0,0);
            VLOG(2) << "EIMFunctionBase::operator() v(x)=" << res << "\n";
            return res;
        }

    virtual element_type const& q( int m )  const = 0;
    virtual vector_type  beta( parameter_type const& mu ) const = 0;
    virtual vector_type  beta( parameter_type const& mu, solution_type const& T ) const = 0;
    virtual size_type  mMax() const = 0;

    virtual std::vector<double> studyConvergence( parameter_type const & mu , solution_type & solution ) const = 0;
    virtual void computationalTimeStatistics( std::string appname )  = 0;
    virtual double expressionL2Norm( solution_type const& T , parameter_type const& mu ) const = 0;
    virtual double diffL2Norm( solution_type const& T , parameter_type const& mu , element_type const& eim_expansion ) const = 0;
    virtual double projExpressionL2Norm( solution_type const& T , parameter_type const& mu ) const = 0;
    virtual double projDiffL2Norm( solution_type const& T , parameter_type const& mu , element_type const& eim_expansion ) const = 0;
    virtual double interpolationError( solution_type const& T , parameter_type const& mu ) const = 0;

    po::variables_map M_vm;
    functionspace_ptrtype M_fspace;
    parameterspace_ptrtype M_pspace;
    sampling_ptrtype M_trainset;
    std::string M_modelname;
    std::string M_name;
    Eigen::VectorXd M_online_time;//contains online computational time
};



template<typename ModelType, typename SpaceType, typename ExprType>
class EIMFunction
    : public EIMFunctionBase<SpaceType, typename ModelType::functionspace_type, typename ModelType::parameterspace_type>
{
    typedef EIMFunctionBase<SpaceType, typename ModelType::functionspace_type, typename ModelType::parameterspace_type> super;
public:
    typedef ModelType model_type;
    //typedef ModelType* model_ptrtype;
    typedef boost::shared_ptr<model_type> model_ptrtype;

    typedef SpaceType functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef typename super::element_type element_type;
    typedef typename super::element_ptrtype element_ptrtype;

    typedef typename ModelType::element_type solution_type;

    typedef typename super::parameterspace_type parameterspace_type;
    typedef typename super::parameter_type parameter_type;
    typedef typename super::sampling_ptrtype sampling_ptrtype;

    typedef ExprType expr_type;
    typedef boost::shared_ptr<expr_type> expr_ptrtype;

    typedef typename super::eim_type eim_type;
    typedef typename super::eim_ptrtype eim_ptrtype;

    typedef typename super::vector_type vector_type;

    typedef CRBModel<ModelType> crbmodel_type;
    typedef boost::shared_ptr<crbmodel_type> crbmodel_ptrtype;
    typedef CRB<crbmodel_type> crb_type;
    typedef boost::shared_ptr<crb_type> crb_ptrtype;

    EIMFunction( po::variables_map const& vm,
                 model_ptrtype model,
                 functionspace_ptrtype space,
                 solution_type& u,
                 parameter_type& mu,
                 expr_type& expr,
                 sampling_ptrtype sampling,
                 std::string const& name )
        :
        super( vm, space, model->parameterSpace(), sampling, model->modelName(), name ),
        M_vm( vm ),
        M_model( model ),
        M_expr( expr ),
        M_u( u ),
        M_mu( mu ),
        M_eim( new eim_type( vm, this, sampling ) )
        {
        }
    element_type operator()( parameter_type const&  mu )
        {
            M_mu = mu;
#if !defined(NDEBUG)
            M_mu.check();
#endif
            M_u = M_model->solve( mu );
            //LOG(INFO) << "operator() mu=" << mu << "\n" << "sol=" << M_u << "\n";
            return vf::project( _space=this->functionSpace(), _expr=M_expr );
        }
    element_type operator()( solution_type const& T, parameter_type const&  mu )
        {
            M_mu = mu;
#if !defined(NDEBUG)
            M_mu.check();
#endif
            // no need to solve we have already an approximation (typically from
            // an nonlinear iteration procedure)
            M_u = T;
            return vf::project( _space=this->functionSpace(), _expr=M_expr );
        }

    //Let g the expression that we want to have an eim expasion
    //here is computed its l2 norm : || g ||_L2
    double expressionL2Norm( solution_type const& T , parameter_type const& mu ) const
    {
        M_mu = mu;
#if !defined(NDEBUG)
        M_mu.check();
#endif
        M_u = T;
        auto mesh = this->functionSpace()->mesh();
        return math::sqrt( integrate( _range=elements( mesh ), _expr=M_expr*M_expr ).evaluate()( 0,0 ) );
    }

    //Let geim the eim expansion of g
    //here is computed || g - geim ||_L2
    double diffL2Norm(  solution_type const& T , parameter_type const& mu , element_type const & eim_expansion ) const
    {
        M_mu = mu;
#if !defined(NDEBUG)
        M_mu.check();
#endif
        M_u = T;
        auto mesh = this->functionSpace()->mesh();
        auto difference = M_expr - idv(eim_expansion);
        return math::sqrt( integrate( _range=elements( mesh ), _expr=difference*difference ).evaluate()( 0,0 ) );
    }


    //Let \pi_g the projection of g on the function space
    //here is computed || \pi_g ||_L2
    double projExpressionL2Norm( solution_type const& T , parameter_type const& mu ) const
    {
        M_mu = mu;
#if !defined(NDEBUG)
        M_mu.check();
#endif
        M_u = T;
        auto mesh = this->functionSpace()->mesh();
        auto pi_g = vf::project( _space=this->functionSpace(), _expr=M_expr );
        return math::sqrt( integrate( _range=elements( mesh ), _expr=idv(pi_g)*idv(pi_g) ).evaluate()( 0,0 ) );
    }

    //here is computed || \pi_g - geim ||_L2
    double projDiffL2Norm( solution_type const& T , parameter_type const& mu , element_type const& eim_expansion ) const
    {
        M_mu = mu;
#if !defined(NDEBUG)
        M_mu.check();
#endif
        M_u = T;
        auto mesh = this->functionSpace()->mesh();
        auto pi_g = vf::project( _space=this->functionSpace(), _expr=M_expr );
        auto diff = pi_g - eim_expansion;
        return math::sqrt( integrate( _range=elements( mesh ), _expr=idv(diff)*idv(diff) ).evaluate()( 0,0 ) );
    }

    double interpolationError( solution_type const& T , parameter_type const& mu ) const
    {
        M_mu = mu;
#if !defined(NDEBUG)
        M_mu.check();
#endif
        M_u = T;
        auto mesh = this->functionSpace()->mesh();
        auto pi_g = vf::project( _space=this->functionSpace(), _expr=M_expr );
        auto difference = M_expr - idv(pi_g);
        return math::sqrt( integrate( _range=elements( mesh ), _expr=difference*difference ).evaluate()( 0,0 ) );
    }

    void computationalTimeStatistics( std::string appname )
        {
            computationalTimeStatistics( appname, typename boost::is_base_of<ModelCrbBaseBase,model_type>::type() );
        }
    void computationalTimeStatistics( std::string appname, boost::mpl::bool_<false> )
        {}
    void computationalTimeStatistics( std::string appname, boost::mpl::bool_<true> )
    {
        //auto crbmodel = crbmodel_ptrtype( new crbmodel_type( M_vm , CRBModelMode::CRB ) );
        auto crbmodel = crbmodel_ptrtype( new crbmodel_type( M_model , CRBModelMode::CRB ) );
        //make sure that the CRB DB is already build
        M_crb = crb_ptrtype( new crb_type( appname,
                                           M_vm ,
                                           crbmodel ) );

        if ( !M_crb->isDBLoaded() || M_crb->rebuildDB() )
        {
            if( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
                LOG( INFO ) << "No CRB DB available, do crb offline computations...";
            M_crb->setOfflineStep( true );
            M_crb->offline();
        }

        int n_eval = option(_name="eim.computational-time-neval").template as<int>();

        Eigen::Matrix<double, Eigen::Dynamic, 1> time_crb;
        time_crb.resize( n_eval );

        typename crb_type::sampling_ptrtype Sampling( new typename crb_type::sampling_type( M_model->parameterSpace() ) );
        Sampling->logEquidistribute( n_eval  );

        //dimension
        int N =  option(_name="crb.dimension").template as<int>();
        //reduced basis approximation space
        auto WN = M_crb->wn();
        int mu_number = 0;
        BOOST_FOREACH( auto mu, *Sampling )
        {
            //LOG( INFO ) << "[computational] mu = \n"<<mu;
            boost::mpi::timer tcrb;
            auto o = M_crb->run( mu,  option(_name="crb.online-tolerance").template as<double>() , N);
            auto solutions=o.template get<2>();
            auto uN = solutions.template get<0>();

            auto u_crb = M_crb->expansion( uN , N , WN );

            boost::mpi::timer teim;
            this->beta( mu , u_crb );
            this->addOnlineTime( teim.elapsed() );
            time_crb( mu_number ) = tcrb.elapsed() ;
            mu_number++;
        }

        M_model->computeStatistics( super::onlineTime() , super::name() );
        M_model->computeStatistics( time_crb , super::name()+" - global crb timing" );

    }

    void setTrainSet( sampling_ptrtype tset ) { M_eim->setTrainSet( tset ); }
    element_type interpolant( parameter_type const& mu ) { return M_eim->operator()( mu , M_eim->mMax() ); }

    element_type const& q( int m ) const { return M_eim->q( m ); }

    vector_type  beta( parameter_type const& mu ) const { return M_eim->beta( mu ); }
    vector_type  beta( parameter_type const& mu, solution_type const& T ) const { return M_eim->beta( mu, T ); }

    std::vector<double> studyConvergence( parameter_type const & mu , solution_type & solution ) const { return M_eim->studyConvergence( mu , solution ) ; };

    size_type mMax() const { return M_eim->mMax(); }

private:
    po::variables_map M_vm;
    model_ptrtype M_model;
    expr_type M_expr;
    solution_type& M_u;
    parameter_type& M_mu;
    eim_ptrtype M_eim;
    crb_ptrtype M_crb;
};

namespace detail
{

template<typename Args>
struct compute_eim_return
{
    typedef typename boost::remove_reference<typename boost::remove_pointer<typename parameter::binding<Args, tag::model>::type>::type>::type::element_type model1_type;
    typedef typename boost::remove_const<typename boost::remove_pointer<model1_type>::type>::type model_type;
    typedef typename boost::remove_reference<typename parameter::binding<Args, tag::expr>::type>::type expr_type;
    typedef typename boost::remove_reference<typename parameter::binding<Args, tag::space>::type>::type::element_type space_type;
    typedef EIMFunction<model_type, space_type, expr_type> type;
    typedef boost::shared_ptr<EIMFunction<model_type, space_type, expr_type> > ptrtype;
};
}

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
      ( options, *, Environment::vm())
      //( space, *( boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> > ), model->functionSpace() )
      //( space, *, model->functionSpace() )
      ( sampling, *, model->parameterSpace()->sampling() )
      ( verbose, (int), 0 )
        ) // optionnal
)
{
    Feel::detail::ignore_unused_variable_warning( args );
    typedef typename Feel::detail::compute_eim_return<Args>::type eim_type;
    typedef typename Feel::detail::compute_eim_return<Args>::ptrtype eim_ptrtype;
    return  eim_ptrtype(new eim_type( options, model, space, element, parameter, expr, sampling, name ) );
} // eim

template<typename ModelType>
struct EimFunctionNoSolve
{
    typedef typename ModelType::functionspace_type functionspace_type;
    typedef typename ModelType::functionspace_ptrtype functionspace_ptrtype;
    typedef typename functionspace_type::element_type element_type;
    typedef typename ModelType::parameterspace_type parameterspace_type;
    typedef typename ModelType::parameterspace_ptrtype parameterspace_ptrtype;
    typedef typename parameterspace_type::element_type parameter_type;
    typedef typename parameterspace_type::sampling_type sampling_type;
    typedef typename parameterspace_type::sampling_ptrtype sampling_ptrtype;
    typedef typename functionspace_type::value_type value_type;

    typedef boost::shared_ptr<ModelType> model_ptrtype;

    static const int nb_spaces = functionspace_type::nSpaces;
    typedef typename mpl::if_< boost::is_same< mpl::int_<nb_spaces> , mpl::int_<2> > , fusion::vector< mpl::int_<0>, mpl::int_<1> >  ,
                       typename mpl::if_ < boost::is_same< mpl::int_<nb_spaces> , mpl::int_<3> > , fusion::vector < mpl::int_<0> , mpl::int_<1> , mpl::int_<2> > ,
                                  typename mpl::if_< boost::is_same< mpl::int_<nb_spaces> , mpl::int_<4> >, fusion::vector< mpl::int_<0>, mpl::int_<1>, mpl::int_<2>, mpl::int_<3> >,
                                                     fusion::vector< mpl::int_<0>, mpl::int_<1>, mpl::int_<2>, mpl::int_<3>, mpl::int_<4> >
                                                     >::type >::type >::type index_vector_type;

    EimFunctionNoSolve( model_ptrtype model )
        :
        M_model( model ),
        M_elt( M_model->functionSpace()->element() )
        {}

    element_type solve( parameter_type const& mu )
    {
        DVLOG(2) << "no solve required\n";
        static const bool is_composite = functionspace_type::is_composite;
        return solve( mu , mpl::bool_< is_composite >() );
    }
    element_type solve( parameter_type const& mu , mpl::bool_<false> )
    {
        value_type x = boost::lexical_cast<value_type>("inf");
        M_elt = vf::project( _space=M_model->functionSpace(), _expr=cst(x) );
        return M_elt;
    }
    element_type solve( parameter_type const& mu , mpl::bool_<true> )
    {
        ProjectInfCompositeCase project_inf_composite_case( M_elt );
        index_vector_type index_vector;
        fusion::for_each( index_vector, project_inf_composite_case );
        return project_inf_composite_case.element();
    }

    std::string modelName() const { return M_model->modelName(); }
    functionspace_ptrtype functionSpace() { return M_model->functionSpace(); }
    parameterspace_ptrtype parameterSpace() { return M_model->parameterSpace(); }

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
            auto space = view.functionSpace();
            view = vf::project( _space=space, _expr=cst( boost::lexical_cast<value_type>("inf") ) );
        }

        element_type element()
        {
            return M_element;
        }

        element_type M_element;

    }; //struct ProjectInfOnSubspace

    model_ptrtype M_model;
    element_type M_elt;
};

template<typename ModelType>
boost::shared_ptr<EimFunctionNoSolve<ModelType>>
eim_no_solve( boost::shared_ptr<ModelType> model )
{
    return boost::shared_ptr<EimFunctionNoSolve<ModelType>>( new EimFunctionNoSolve<ModelType>( model ) );
}


po::options_description eimOptions( std::string const& prefix ="");
}
#endif /* _FEELPP_EIM_HPP */
