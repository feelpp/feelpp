/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2012-05-02
 */
#ifndef _FEELPP_EIM_HPP
#define _FEELPP_EIM_HPP 1

#include <limits>
#include <string>

#include <boost/ref.hpp>
#include <boost/next_prior.hpp>
#include <boost/type_traits.hpp>
#include <boost/tuple/tuple.hpp>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/base_object.hpp>

#include <feel/feelcrb/parameterspace.hpp>

#include <feel/feelvf/vf.hpp>

#include <Eigen/Core>

namespace Feel
{
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

  @author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
    typedef typename ModelType::functionspace_ptrtype functionspace_ptrtype;
    typedef typename functionspace_type::element_type element_type;
    typedef typename functionspace_type::element_ptrtype element_ptrtype;

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
        M_is_read( false ),
        M_is_written( false ),
        M_name( "default" ),
        M_M( 1 ),
        M_M_max( 1 ),
        M_offline_done( false ),
        M_tol( 1e-8 ),
        M_q(),
        M_B(),
        M_t(),
        M_index_max(),
        M_model( 0 )
        {}
    EIM( model_type* model, double __tol = 1e-8 )
        :
        M_is_read( false ),
        M_is_written( false ),
        M_name( model->name() ),
        M_M( 1 ),
        M_M_max( 1 ),
        M_offline_done( false ),
        M_tol( __tol ),
        M_q(),
        M_B(),
        M_t(),
        M_index_max(),
        M_model( model )
        {
            LOG(INFO) << "construct EIM approximation...\n";
            if (  !this->isOfflineDone() )
                offline();
        }

    EIM( EIM const & __bbf )
        :
        M_is_read( __bbf.M_is_read ),
        M_is_written( __bbf.M_is_written ),
        M_name( __bbf.M_name ),
        M_M( __bbf.M_M ),
        M_M_max (__bbf.M_M_max ),
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
    size_t nDOF() const {  FEELPP_ASSERT( M_model != 0 ).error( "Invalid EIM model" ); return M_model->functionSpace()->nLocalDof(); }

    matrix_type const& q() const {  return M_q; }

    column_type
    qm( size_t __m ) const
        {
            FEELPP_ASSERT( __m >= 0 && __m < M_M )( __m )( M_M ).error( "out of bounds access" );

            return M_q.col( __m );
        }
    value_type qxm( size_t __i, size_t __m ) const
        {
            FEELPP_ASSERT( __m >= 0 && __m < M_M )( __m )( M_M ).error( "out of bounds access" );
            FEELPP_ASSERT( __i >= 0 && __i < M_q.rows() )( __i )( M_q.rows() ).error( "out of bounds access" );

            return M_q( __i, __m );
        }


    size_t mMax() const { return M_M_max; }


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
    vector_type beta( parameter_type const& mu, size_t M  ) const;

    element_type residual ( size_t M ) const;

    parameter_residual_type computeBestFit( sampling_ptrtype trainset, int __M );

    /**
       orthonormalize
    */
    void orthonormalize( std::vector<element_type>& );

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & __ar, const unsigned int __version )
        {
            Debug( 7450 ) << "serializing...\n";

            __ar & BOOST_SERIALIZATION_NVP( M_name );

            Debug( 7450 ) << "name saved/loaded...\n";

            __ar & BOOST_SERIALIZATION_NVP( M_offline_done );
            Debug( 7450 ) << "offline status...\n";

            __ar & BOOST_SERIALIZATION_NVP( M_M_max );
            Debug( 7450 ) << "M saved/loaded\n";
            M_M = M_M_max;

            // save index
            __ar & BOOST_SERIALIZATION_NVP( M_index_max );
            Debug( 7450 ) << "index saved/loaded\n";

            // save t
            __ar & BOOST_SERIALIZATION_NVP( M_t );
            Debug( 7450 ) << "t saved/loaded\n";

            // save q
            __ar & BOOST_SERIALIZATION_NVP( M_q );
            Debug( 7450 ) << "q saved/loaded\n";

            // save B
            __ar & BOOST_SERIALIZATION_NVP( M_B );
            Debug( 7450 ) << "B saved/loaded\n";

        }
    //@}

    /**
       Given the \f$\mu\f$ in parameter space,
       compute \f$ q_m( x ) \forall x \in \Omega\f$
    */
    void blackboxQ( vector_type& __g, int __m )
        {
            FEELPP_ASSERT( M_t.size() ).error( "t size is 0" );
            FEELPP_ASSERT( M_q.rows() && M_q.cols() ).error( "q size is 0" );
            FEELPP_ASSERT( M_B.rows() && M_B.cols() ).error( "B size is 0" );
            if ( __g.size() != M_q.rows() )
                __g.resize( M_q.rows() );
            __g = qm( __m );
        }

    /**
       Given the \f$\mu\f$ in parameter space,
       compute \f$ g_M( x, \mu ) \forall x \in \Omega\f$
    */
    void blackbox( vector_type& __g )
        {
            FEELPP_ASSERT( M_t.size() ).error( "t size is 0" );
            FEELPP_ASSERT( M_q.rows() && M_q.cols() ).error( "q size is 0" );
            FEELPP_ASSERT( M_B.rows() && M_B.cols() ).error( "B size is 0" );

            vector_type __b;
            beta( __b, M_M_max );

            __g.conservativeResize( M_q.rows() );

            __g = M_q.block(0,0,M_q.rows(),M_M_max)*__b;
        }
    //@}

protected:

    mutable bool M_is_read;
    mutable bool M_is_written;

    std::string M_name;

    size_t M_M;
    size_t M_M_max;

    mutable bool M_offline_done;

    double M_tol;

    std::vector<element_type> M_g;

    std::vector<element_type> M_q;

    matrix_type M_B;

    std::vector<node_type> M_t;

    std::vector<size_t> M_index_max;

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


};

template<typename ModelType>
typename EIM<ModelType>::vector_type
EIM<ModelType>::beta( parameter_type const& mu, size_t __M ) const
{
    // beta=B_M\g(Od(indx),mut(i))'
    vector_type __beta( __M );
    for ( size_t __m = 0;__m < __M;++__m )
    {
        __beta[__m] = M_model->operator()( this->M_t[__m], mu );
    }
    this->M_B.block(0,0,__M,__M).template triangularView<Eigen::UnitLower>().solveInPlace(__beta);
    return __beta;
}
template<typename ModelType>
typename EIM<ModelType>::element_type
EIM<ModelType>::residual( size_t __M ) const
{
    LOG(INFO) << "compute residual for m=" << __M << "...\n";
    vector_type rhs( __M );
    for ( size_t __m = 0;__m < __M;++__m )
    {
        rhs[__m]= M_g[__M]( M_t[__m] )(0,0,0);
    }
    this->M_B.block(0,0,__M,__M).template triangularView<Eigen::UnitLower>().solveInPlace(rhs);
    LOG(INFO) << "solve B sol = rhs with rhs = " << rhs <<"\n";

    // res(:,i)=M_g(:,i)-q(:,0:i)*sigma
    LOG(INFO) << "compute residual..." <<"\n";
    using namespace vf;
    auto z = expansion( M_q, rhs );
    LOG(INFO) << "return residual..." <<"\n";
    return vf::project( _space=M_model->functionSpace(),
                        _expr=idv(M_g[__M])-idv( z ) );
}

template<typename ModelType>
void
EIM<ModelType>::orthonormalize( std::vector<element_type> & __Z )
{
    FEELPP_ASSERT( __Z.rows() > 0 && __Z.cols() > 0 )( __Z.rows() )( __Z.cols() ).error( "invalid number of rows or columns" );
    // here we use the norm2 but it might be a good one
    size_t __M = __Z.cols();
    for ( size_t __i = 0;__i < __M-1; ++__i )
    {
        value_type __s = __Z.col(__i)*__Z.col(__M-1);
        __Z.col(__M-1) -= __s*__Z.col(__i);
    }
    value_type __normZ = __Z.col(__M-1).norm();;
    __Z.col(__M-1).normalize();
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
        LOG_EVERY_N(INFO, 10 ) << " (every 10 mu) compute fit at mu="<< mu <<"\n" ;
        // evaluate model at mu
        auto Z = M_model->operator()( mu );

        vector_type rhs( __M );
        for ( size_t __m = 0;__m < __M;++__m )
        {
            rhs[__m]= Z( M_t[__m] )(0,0,0);
        }
        this->M_B.block(0,0,__M,__M).template triangularView<Eigen::UnitLower>().solveInPlace(rhs);
        auto res = vf::project( _space=M_model->functionSpace(),
                                _expr=idv(Z)-idv( expansion( M_q, rhs ) ) );
        auto resmax = normLinf( _range=elements(M_model->mesh()), _pset=_Q<5>(), _expr=idv(res) );
        maxerr( index ) = resmax.template get<0>();
        int index2;
        auto err = maxerr.array().abs().maxCoeff( &index2 );
        LOG_EVERY_N(INFO, 10 ) << " (every 10 mu) maxerr=" <<  err << " at index = " << index2 << " at mu = " << trainset->at(index2) << "\n";
    }
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

    LOG(INFO) << "create training set...\n";
    auto trainset = M_model->parameterSpace()->sampling();
    trainset->randomize( 30 );

    LOG(INFO) << "create mu_1...\n";
    // random element in Dmu to start with
    auto mu = M_model->parameterSpace()->element();

    LOG(INFO) << "compute finite element solution at mu_1...\n";
    M_g.push_back( M_model->operator()( mu ) );

    LOG(INFO) << "compute T^" << 0 << "...\n";
    // Build T^0
    auto zmax = normLinf( _range=elements(M_model->mesh()), _pset=_Q<5>(), _expr=idv(M_g[0]) );
    // store space coordinate where max absolute value occurs
    M_t.push_back( zmax.template get<1>() );
    LOG(INFO) << "norm Linf = " << zmax.template get<0>() << " at " << zmax.template get<1>() << "\n";

    LOG(INFO) << "compute and insert q_0...\n";
    // insert first element
    auto q = M_g[0];
    q.scale( zmax.template get<0>() );
    M_q.push_back( q );

    LOG(INFO) << "compute entry (0,0) of interpolation matrix...\n";
    this->M_B.resize( 1, 1 );
    this->M_B( 0, 0 ) = 1;

    /**
       \par build \f$W^g_M\f$
    */
    double err = 1;

    LOG(INFO) << "start greedy algorithm...\n";
    while (  err >= this->M_tol )
    {
        LOG(INFO) << "M=" << M_M << "...\n";

        LOG(INFO) << "compute best fit error...\n";
        // compute mu = arg max inf ||G(.;mu)-z||_infty
        auto bestfit = computeBestFit( trainset, this->M_M );

        if ( this->M_M == 1 && bestfit.template get<0>() < 1e-12 )
            break;

        /**
         * we have a new \f$\mu^g_M\f$, insert it in \f$S^g_M\f$ and the
         * corresponding \f$z \in W^g_M\f$
         */
        ++this->M_M;

        LOG(INFO) << "S(" << this->M_M-1 << ") = " << bestfit.template get<1>() << "\n";

        // update M_g(:,M-1)
        M_g.push_back( M_model->operator()( bestfit.template get<1>() ) );

        //orthonormalize( M_g );

        // build T^m such that T^m-1 \subset T^m
        M_B.conservativeResize( M_M, M_M );

        for( int __i = 0; __i < M_M; ++__i )
        {
            for( int __j = 0; __j < __i; ++__j )
            {
                this->M_B( __i, __j ) = M_q[__j]( M_t[__i] )(0,0,0);

            }
            this->M_B(__i,__i)=1;
        }
        LOG(INFO) << "Interpolation matrix: M_B = " << this->M_B <<"\n";

        LOG(INFO) << "compute residual M="<< M_M-1 << "..." <<"\n";
        auto res = this->residual(M_M-1);

        LOG(INFO) << "compute arg sup |residual|..." <<"\n";
        auto resmax = normLinf( _range=elements(M_model->mesh()), _pset=_Q<5>(), _expr=idv(res) );
        LOG(INFO) << "store coordinates where max absolute value is attained " << resmax.template get<1>() << "..." <<"\n";
        // store space coordinate where max absolute value occurs
        M_t.push_back( resmax.template get<1>() );

        LOG(INFO) << "store new basis function..." <<"\n";
        res.scale( 1./resmax.template get<0>() );
        // insert new q
        M_q.push_back( res );

        if ( this->M_M > 20 )
            break;
    }

    LOG(INFO) << "M_max = " << M_M << "...\n";
    this->M_M_max = this->M_M;

    this->M_offline_done = true;

}

template<typename SpaceType, typename ParameterSpaceType>
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
    static const uint16_type nDim = mesh_type::nDim;

    typedef ParameterSpaceType parameterspace_type;
    typedef boost::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef typename parameterspace_type::element_type parameter_type;

    typedef EIM<EIMFunctionBase<SpaceType,ParameterSpaceType> > eim_type;
    typedef boost::shared_ptr<eim_type> eim_ptrtype;
    //typedef typename eim_type::betam_type betam_type;
    //typedef typename eim_type::qm_type qm_type;

    typedef Eigen::Matrix<double, nDim, 1> node_type;

    EIMFunctionBase( functionspace_ptrtype fspace,
                     parameterspace_ptrtype pspace,
                     std::string const& name )
        :
        M_fspace( fspace ),
        M_pspace( pspace ),
        M_name( name )
        {
            LOG(INFO)<< "EimFunctionBase constructor\n";
        }
    virtual ~EIMFunctionBase()
        {}
    std::string const& name() const { return M_name; }
    void setName( std::string const& name ) { M_name = name; }

    functionspace_ptrtype functionSpace() const { return M_fspace; }
    functionspace_ptrtype functionSpace()  { return M_fspace; }
    parameterspace_ptrtype parameterSpace() const { return M_pspace; }
    parameterspace_ptrtype parameterSpace()  { return M_pspace; }
    mesh_ptrtype mesh() const { return M_fspace->mesh(); }
    mesh_ptrtype mesh()  { return M_fspace->mesh(); }

    virtual element_type operator()( parameter_type const& ) = 0;
    value_type operator()( node_type const& x, parameter_type const& mu )
        {
            LOG(INFO) << "calling EIMFunctionBase::operator()( x=" << x << ", mu=" << mu << ")\n";
            element_type v = this->operator()( mu );
            value_type res = v(x);
            LOG(INFO) << "EIMFunctionBase::operator() v(x)=" << res << "\n";
            return res;
        }
#if 0
    qm_type q( parameter_type const& ) { return M_eim->blackboxQ( mu ); }
    _type q( parameter_type const& ) { return M_eim->blackboxQ( mu ); }
#endif

    functionspace_ptrtype M_fspace;
    parameterspace_ptrtype M_pspace;
    std::string M_name;
    eim_ptrtype M_eim;
};



template<typename ModelType, typename ExprType>
class EIMFunction
    : public EIMFunctionBase<typename ModelType::functionspace_type, typename ModelType::parameterspace_type>
{
    typedef EIMFunctionBase<typename ModelType::functionspace_type, typename ModelType::parameterspace_type> super;
public:
    typedef ModelType model_type;
    typedef ModelType* model_ptrtype;

    typedef typename super::functionspace_type functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef typename super::element_type element_type;
    typedef typename super::element_ptrtype element_ptrtype;

    typedef typename super::parameterspace_type parameterspace_type;
    typedef typename super::parameter_type parameter_type;

    typedef ExprType expr_type;
    typedef boost::shared_ptr<expr_type> expr_ptrtype;

    typedef typename super::eim_type eim_type;
    typedef typename super::eim_ptrtype eim_ptrtype;

    EIMFunction( model_ptrtype model,
                 element_type& u,
                 parameter_type& mu,
                 expr_type& expr,
                 std::string const& name )
        :
        super( model->functionSpace(), model->parameterSpace(), name ),
        M_model( model ),
        M_expr( expr ),
        M_u( u ),
        M_mu( mu ),
        M_eim( new eim_type( this ) )
        {

        }
    element_type operator()( parameter_type const&  mu )
        {
            M_mu = mu;
            M_u = M_model->solve( mu );
            return vf::project( _space=M_model->functionSpace(), _expr=M_expr );
        }

private:
    model_ptrtype M_model;
    expr_type M_expr;
    element_type& M_u;
    parameter_type& M_mu;
    eim_ptrtype M_eim;
};

namespace detail
{

template<typename Args>
struct compute_eim_return
{
    typedef typename boost::remove_reference<typename boost::remove_pointer<typename parameter::binding<Args, tag::model>::type>::type>::type model1_type;
    typedef typename boost::remove_const<typename boost::remove_pointer<model1_type>::type>::type model_type;
    typedef typename boost::remove_reference<typename parameter::binding<Args, tag::expr>::type>::type expr_type;
    typedef EIMFunction<model_type, expr_type> type;
    typedef boost::shared_ptr<EIMFunction<model_type, expr_type> > ptrtype;
};
}

BOOST_PARAMETER_FUNCTION(
    ( typename detail::compute_eim_return<Args>::ptrtype ), // 1. return type
    eim,                        // 2. name of the function template
    tag,                                        // 3. namespace of tag types
    ( required
      ( in_out(model),          * )
      ( in_out(element),        * )
      ( in_out(parameter),        * )
      ( in_out(expr),          * )
      ( name, * )
        ) // required
    ( optional
      ( verbose, (int), 0 )
        ) // optionnal
)
{
    Feel::detail::ignore_unused_variable_warning( args );
    typedef typename detail::compute_eim_return<Args>::type eim_type;
    typedef typename detail::compute_eim_return<Args>::ptrtype eim_ptrtype;
    return  eim_ptrtype(new eim_type( model, element, parameter, expr, name ) );
} // eim



}
#endif /* _FEELPP_EIM_HPP */

