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

#include <feel/feelvf/vf.hpp>

#include <Eigen/Core>
#include <glpk.h>

namespace Feel
{
template<typename T>
class EIMBase
{
public:

    /** @name Typedefs
     */
    //@{

    typedef T value_type;

    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_type;
    typedef typename matrix_type::ColXpr column_type;
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector_type;
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> node_type;

    //@}

    EIMBase()
        :
        _M_is_read( false ),
        _M_is_written( false ),
        _M_name( "default" ),
        _M_M( 1 ),
        _M_M_max( 1 ),
        _M_offline_done( false ),
        _M_tol( 1e-8 ),
        _M_q(),
        _M_B(),
        _M_t(),
        _M_index_max()
        {
        }
    EIMBase( std::string const& __name, double __tol = 1e-8 )
        :
        _M_is_read( false ),
        _M_is_written( false ),
        _M_name( __name ),
        _M_M( 1 ),
        _M_M_max( 1 ),
        _M_offline_done( false ),
        _M_tol( __tol ),
        _M_q(),
        _M_B(),
        _M_t(),
        _M_index_max()
        {
        }
    EIMBase( EIMBase const& __bbf )
        :
        _M_is_read( __bbf._M_is_read ),
        _M_is_written( __bbf._M_is_written ),
        _M_name( __bbf._M_name ),
        _M_M( __bbf._M_M ),
        _M_M_max (__bbf._M_M_max ),
        _M_offline_done( __bbf._M_offline_done ),
        _M_tol( __bbf._M_tol ),
        _M_q( __bbf._M_q ),
        _M_B( __bbf._M_B ),
        _M_t( __bbf._M_t ),
        _M_index_max( __bbf._M_index_max )
        {}
    virtual ~EIMBase() {}

    /** @name Accessors
     */
    //@{

    /**
       \return the number of DOF in space discretization
    */
    size_t nbDOF() const {  return _M_q.rows(); }

    matrix_type const& q() const {  return _M_q; }

    column_type
    qm( size_t __m ) const
        {
            FEELPP_ASSERT( __m >= 0 && __m < _M_M )( __m )( _M_M ).error( "out of bounds access" );

            return _M_q.col( __m );
        }
    value_type qxm( size_t __i, size_t __m ) const
        {
            FEELPP_ASSERT( __m >= 0 && __m < _M_M )( __m )( _M_M ).error( "out of bounds access" );
            FEELPP_ASSERT( __i >= 0 && __i < _M_q.rows() )( __i )( _M_q.rows() ).error( "out of bounds access" );

            return _M_q( __i, __m );
        }


    size_t Mmax() const { return _M_M_max; }


    //@}

    /** @name  Methods
     */
    //@{

    /**
       orthonormalize
    */
    void orthonormalize( matrix_type& );

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & __ar, const unsigned int __version )
        {
            Debug( 7450 ) << "serializing...\n";

            __ar & BOOST_SERIALIZATION_NVP( _M_name );

            Debug( 7450 ) << "name saved/loaded...\n";

            __ar & BOOST_SERIALIZATION_NVP( _M_offline_done );
            Debug( 7450 ) << "offline status...\n";

            __ar & BOOST_SERIALIZATION_NVP( _M_M_max );
            Debug( 7450 ) << "M saved/loaded\n";
            _M_M = _M_M_max;

            // save index
            __ar & BOOST_SERIALIZATION_NVP( _M_index_max );
            Debug( 7450 ) << "index saved/loaded\n";

            // save t
            __ar & BOOST_SERIALIZATION_NVP( _M_t );
            Debug( 7450 ) << "t saved/loaded\n";

            // save q
            __ar & BOOST_SERIALIZATION_NVP( _M_q );
            Debug( 7450 ) << "q saved/loaded\n";

            // save B
            __ar & BOOST_SERIALIZATION_NVP( _M_B );
            Debug( 7450 ) << "B saved/loaded\n";

        }

    /**
       offline stage
    */
    virtual void offline(  )
        {}

    /**
       Check whether the offline stage has been executed and the database created.
       \return \c true if the offline has been executed and the DB saved, \c false otherwise
    */
    bool isOfflineDone() const
        {
            return _M_offline_done;
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
    virtual void beta( vector_type& __beta, size_t M  ) const
        {}

    virtual void beta( vector_type& __beta, value_type& __err, size_t M ) const
        {}

    /**
       Given the \f$\mu\f$ in parameter space,
       compute \f$ q_m( x ) \forall x \in \Omega\f$
    */
    void blackboxQ( vector_type& __g, int __m )
        {
            FEELPP_ASSERT( _M_t.size() ).error( "t size is 0" );
            FEELPP_ASSERT( _M_q.rows() && _M_q.cols() ).error( "q size is 0" );
            FEELPP_ASSERT( _M_B.rows() && _M_B.cols() ).error( "B size is 0" );
            if ( __g.size() != _M_q.rows() )
                __g.resize( _M_q.rows() );
            __g = qm( __m );
        }

    /**
       Given the \f$\mu\f$ in parameter space,
       compute \f$ g_M( x, \mu ) \forall x \in \Omega\f$
    */
    void blackbox( vector_type& __g )
        {
            FEELPP_ASSERT( _M_t.size() ).error( "t size is 0" );
            FEELPP_ASSERT( _M_q.rows() && _M_q.cols() ).error( "q size is 0" );
            FEELPP_ASSERT( _M_B.rows() && _M_B.cols() ).error( "B size is 0" );

            vector_type __b;
            beta( __b, _M_M_max );

            __g.conservativeResize( _M_q.rows() );

            __g = _M_q.block(0,0,_M_q.rows(),_M_M_max)*__b;
        }
    //@}

protected:

    mutable bool _M_is_read;
    mutable bool _M_is_written;

    std::string _M_name;

    size_t _M_M;
    size_t _M_M_max;

    mutable bool _M_offline_done;

    double _M_tol;


    matrix_type _M_q;

    matrix_type _M_B;

    std::vector<node_type> _M_t;

    std::vector<size_t> _M_index_max;

};

template<typename T>
void
EIMBase<T>::orthonormalize( matrix_type& __Z )
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
    :
        public EIMBase<typename ModelType::value_type>
{
public:


    /** @name Typedefs
     */
    //@{
    typedef typename ModelType::value_type value_type;
    typedef EIMBase<value_type> super;

    typedef typename super::matrix_type matrix_type;
    typedef typename super::vector vector_type;
    typedef typename ModelType::space_type space_type;
    typedef typename ModelType::space_ptrtype space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef typename space_type::element_ptrtype element_ptrtype;
    typedef ModelType f_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    EIM( space_ptrtype const& Xh, double __tol = 1e-8 )
        :
        super( f_type().name(), __tol ),
        M_Xh( Xh )
        {
            if (  !this->isOfflineDone() )
                offline( M_Xh );
        }

    EIM( EIM const & __bbf )
        :
        super( __bbf )
        {}
    ~EIM()
        {}

    //static EIMBase<value_type>* create() { return new EIM<value_type, Dmu, F>; }

    //@}

    /** @name Operator overloads
     */
    //@{


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
    void beta( vector_type& __beta, size_t M  ) const;

    void beta( vector_type& __beta, value_type& __err, size_t M ) const;

    template<class Archive>
    void
    serialize( Archive& __ar, const unsigned int __version )
        {
            __ar & BOOST_SERIALIZATION_NVP( boost::serialization::base_object<super>(*this) );
        }
    //@}


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


private:

    mutable f_type _M_f;
    space_ptrtype M_Xh;
    std::vector<element_type> _M_q;
};

template<typename ModelType>
void
EIM<ModelType>::beta( vector_type& __beta, size_t __M ) const
{
    if ( __M <= 0 )
        __M = 1;
    if ( __M >= this->Mmax() )
        __M = this->Mmax();
    // beta=B_M\g(Od(indx),mut(i))'
    __beta.resize( __M );
    for ( size_t __m = 0;__m < __M;++__m )
    {
        //TODO: __beta[__m] = _M_f( this->_M_t[__m], Dmu::instance().current() );
    }


    //TODO: this->_M_B.block(0,0,__M,__M).triangularView<Eigen::UnitLower>().solveInPlace(__beta);
}
template<typename ModelType>
void
EIM<ModelType>::beta( vector_type& __beta, value_type& __err, size_t __M ) const
{
    if ( __M <= 0 )
        __M = 1;
    if ( __M >= this->Mmax() )
        __M = this->Mmax()-1;
    if (  this->Mmax() == 1 )// affine case
    {
        beta( __beta, __M );
        __err = 0;
        return;
    }

    beta( __beta, __M );
    //TODO: value_type gmu = _M_f( this->_M_t[__M], Dmu::instance().current() ); // g(t_{M+1})
    value_type gMmu = 0; // gM(t_{M+1})
    for ( size_t __m = 0; __m < __M;++__m )
        gMmu += __beta[ __m ] * qxm( this->_M_index_max[__M], __m );
    //TODO: __err = std::abs( gmu - gMmu );
}

template<typename ModelType>
void
EIM<ModelType>::offline(  )
{
#if 0
    if ( this->isOfflineDone() )
        return;


    matrix_type __Z( M_Xh->nLocalDof(), this->_M_M );

    //TODO: __Z.col( this->_M_M-1 ) = expr
    //TODO: __f.assign_special( gmm::mat_col( __Z, this->_M_M-1 ) );

    orthonormalize( __Z );

    //std::cout << "Z(" << this->_M_M-1 << ") = " << __Z.col( this->_M_M-1 ) << "\n";
    //std::cout << "\n";

    /**
       \par build \f$W^g_M\f$


    */
    double err = 1;

    while (  err >= this->_M_tol )
    {
        /**
           build the \f$ M\times M\f$ matrix \f$ A_z = Z^T * Z\f$
        */
        matrix_type __A_z(this->_M_M, this->_M_M);
        __A_z = __Z.transpose()*__Z;
        //std::cout << "Az = " << __A_z << "\n";

        /**
           linear program to find the \f$ sup_\mu h(\mu) \f$ where \f$ h(mu) = inf_z ||f-z||_\infty\f$
        */
        //TODO: Dmu::instance().setRandom();

#if 0
        boost::tie( mu, err )  = lp( _data = data );

        // push mu into parameter set
#endif

        if ( this->_M_M == 1 && err < 1e-12 )
            break;

        /** we have a new \f$\mu^g_M\f$, insert it in \f$S^g_M\f$ and the corresponding \f$z \in W^g_M\f$
         */
        ++this->_M_M;

        //TODO: std::cout << "S(" << this->_M_M-1 << ") = " << __S << "\n";

        // update Z(:,M-1)
        //TODO: __Z.col(this->_M_M-1) = XXXXX;

        orthonormalize( __Z );

        //std::cout << "Z(" << this->_M_M-1 << ") = " << gmm::mat_col( __Z, this->_M_M-1 ) << "\n";
        //std::cout << "\n";

        if ( this->_M_M > 20 )
            break;
    }


    // Build T^M
    this->_M_t.conservativeResize( this->_M_M );

    auto __z_col = __Z.col(0);
    int index_max;
    value_type __z_max = __z_col.maxCoeff(&index_max);

    this->_M_index_max.push_back( index_max);
    std::cout << "index = [ ";
    std::copy( this->_M_index_max.begin(), this->_M_index_max.end(), std::ostream_iterator<size_t>( std::cout, " " ) );
    std::cout << " ]\n";

    // store t_0
    this->_M_t[0] = M_Xh->dofPoint( this->_M_index_max[0] ).template get<0>();

    this->_M_q.resize( this->_M_M );
    // copy Z(:,1) in q(:,1)
    this->_M_q[0] = __z_col;
    // divide by max(Z(:,1))
    this->_M_q[0]/=__z_max;

    // residual : res(:,1)=Z(:,1)
    matrix_type __res( M_Xh->nLocalDof(), this->_M_M );
    __res.col(0) = __z_col;

    this->_M_B.resize( 1, 1 );
    this->_M_B( 0, 0 ) = 1;

    for ( size_t __i = 1;__i < this->_M_M;++__i )
    {
        this->_M_B.conservativeResize( __i, __i );
        matrix_type __b( __i, 1 );

        this->_M_B = this->_M_q.block( this->_M_index_max, 0, 1, __i );
        __b = __Z.block( this->_M_index_max, __i, 1, 1 );
        //std::cout << "M_B " << __i-1 << " = " << this->_M_B[__i-1] <<"\n";

        vector_type __sigma( __i );
        __sigma = __b.col(0);


        //TODO: this->_M_B.triangularView<Eigen::UnitLower>().solveInPlace(__sigma);

        // res(:,i)=Z(:,i)-q(:,0:i)*sigma
        __res.col(__i) = __Z.col(__i) -__sigma*this->_M_q.block( 0, 0, this->_M_q.rows(), __i );

        auto __res_col = __res.col(__i);
        int index_max;
        auto __res_max = __res_col.maxCoeff(&index_max);
        this->_M_index_max.push_back( index_max );

        std::cout << "index = [ ";
        std::copy( this->_M_index_max.begin(), this->_M_index_max.end(), std::ostream_iterator<size_t>( std::cout, " " ) );
        std::cout << " ]\n";

        // store t_i: TODO
        //this->_M_t[__i] = __mesh.point_of_dof( this->_M_index_max[__i] );


        _M_q.col(__i) = __res.col(__i)/boost::get<0>( __res_max );
    }
    if (  this->_M_M > 1 )
    {
        this->_M_B.conservativeResize( this->_M_M, this->_M_M );
        this->_M_B = _M_q.block(this->_M_index_max,0,1,this->_M_M);
    }
    std::cout << "M_B = " << this->_M_B <<"\n";
    //std::cout << "q(" << this->_M_M-1 << ") = " << gmm::mat_col( _M_q, this->_M_M-1 ) << "\n";
    //std::cout << "\n";

    this->_M_M_max = this->_M_M;
#endif
    this->_M_offline_done = true;

}

template<typename SpaceType, typename ParameterSpaceType>
class EIMFunctionBase
{
public:
    typedef SpaceType space_type;
    typedef typename space_type::element_type element_type;
    typedef ParameterSpaceType parameterspace_type;
    typedef typename parameterspace_type::element_type parameter_type;

    virtual element_type operator()( parameter_type const& ) = 0;
};


template<typename ModelType, typename ExprType>
class EIMFunction
    : public EIMFunctionBase<typename ModelType::functionspace_type, typename ModelType::parameterspace_type>
{
    typedef EIMFunctionBase<typename ModelType::functionspace_type, typename ModelType::parameterspace_type> super;
public:
    typedef ModelType model_type;
    typedef boost::shared_ptr<ModelType> model_ptrtype;
    typedef typename ModelType::functionspace_type space_type;
    typedef typename ModelType::functionspace_ptrtype space_ptrtype;


    typedef typename space_type::element_type element_type;
    typedef typename space_type::element_ptrtype element_ptrtype;
    typedef typename ModelType::parameterspace_type parameterspace_type;
    typedef typename parameterspace_type::element_type parameter_type;
    typedef ExprType expr_type;
    typedef boost::shared_ptr<expr_type> expr_ptrtype;

    EIMFunction( model_ptrtype model, element_ptrtype u, expr_ptrtype expr )
        :
        super( model->functionSpace(), model->parameterSpace() ),
        M_model( model ),
        M_u( u ),
        M_expr( expr )
        {}
    element_type operator()( parameter_type const&  mu )
        {
            *M_u = M_model->solve( mu );
            return vf::project( _space=M_model->functionSpace(), _expr=*M_expr );
        }
private:
    model_ptrtype M_model;
    expr_ptrtype M_expr;
    element_ptrtype M_u;
};



}
#endif /* _FEELPP_EIM_HPP */

