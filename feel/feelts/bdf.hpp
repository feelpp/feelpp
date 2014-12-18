/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2006-12-30

   Copyright (C) 2006-2008 Université Joseph Fourier (Grenoble)

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 3.0 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file bdf.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-12-30
*/
#ifndef _BDF_H
#define _BDF_H

#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <boost/timer.hpp>
#include <boost/shared_array.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/utility.hpp>

#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <boost/parameter.hpp>
#include <feel/feelcore/parameter.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelts/tsbase.hpp>
#include <feel/feeldiscr/functionspace.hpp>

namespace Feel
{
namespace ublas = boost::numeric::ublas;
namespace fs = boost::filesystem;

enum BDFTimeScheme { BDF_ORDER_ONE=1, BDF_ORDER_TWO, BDF_ORDER_THREE, BDF_ORDER_FOUR, BDF_MAX_ORDER = 4 };

/**
 * \class Bdf
 * \ingroup SpaceTime
 * \brief Backward differencing formula time discretization
 *
 * A differential equation of the form
 *
 * \f$ M u' = A u + f \f$
 *
 * is discretized in time as
 *
 * \f$ M p'(t_{k+1}) = A u_{k+1} + f_{k+1} \f$
 *
 * where p denotes the polynomial of order n in t that interpolates
 * \f$ (t_i,u_i) \f$ for \f$ i = k-n+1,...,k+1\f$.
 *
 * The approximative time derivative \f$ p'(t_{k+1}) \f$ is a linear
 * combination of state vectors \f$u_i\f$:
 *
 * \f$ p'(t_{k+1}) = \frac{1}{\Delta t} (\alpha_0 u_{k+1} - \sum_{i=0}^n \alpha_i u_{k+1-i} )\f$
 *
 * Thus we have
 *
 * \f$ \frac{\alpha_0}{\Delta t} M u_{k+1} = A u_{k+1} + f + M \bar{p} \f$
 *
 * with
 *
 * \f$ \bar{p} = \frac{1}{\Delta t} \sum_{i=1}^n \alpha_i u_{k+1-i} \f$
 *
 * This class stores the n last state vectors in order to be able to
 * calculate \f$ \bar{p} \f$. It also provides \f$ \alpha_i \f$
 * and can extrapolate the new state from the n last states with a
 * polynomial of order n-1:
 *
 * \f$ u_{k+1} \approx \sum_{i=0}^{n-1} \beta_i u_{k-i} \f$
 */
template<typename SpaceType>
class Bdf : public TSBase
{
    friend class boost::serialization::access;
    typedef TSBase super;
public:
    typedef Bdf<SpaceType> bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;
    typedef SpaceType space_type;
    typedef boost::shared_ptr<space_type>  space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef typename space_type::return_type return_type;
    typedef typename element_type::value_type value_type;
    typedef boost::shared_ptr<element_type> unknown_type;
    typedef std::vector< unknown_type > unknowns_type;
    typedef typename node<value_type>::type node_type;

    typedef typename super::time_iterator time_iterator;

    /**
     * Constructor
     *
     * @param space approximation space
     * @param n order of the BDF
     */
    Bdf( po::variables_map const& vm, space_ptrtype const& space, std::string const& name, std::string const& prefix="" );

    /**
     * Constructor
     * @param space approximation space
     * @param name name of the BDF
     */
    Bdf( space_ptrtype const& space, std::string const& name );

    //! copy operator
    Bdf( Bdf const& b )
        :
        super( b ),
        M_order( b.M_order ),
        M_strategyHighOrderStart( b.M_strategyHighOrderStart ),
        M_order_cur( b.M_order_cur ),
        M_last_iteration_since_order_change( b.M_last_iteration_since_order_change ),
        M_iterations_between_order_change( b.M_iterations_between_order_change ),
        M_space( b.M_space ),
        M_unknowns( b.M_unknowns ),
        M_alpha( b.M_alpha ),
        M_beta( b.M_beta ),
        M_prefix( b.M_prefix )
    {}

    ~Bdf();


    //! return a deep copy of the bdf object
    bdf_ptrtype deepCopy() const
    {
        auto b = bdf_ptrtype( new bdf_type( *this ) );

        for ( auto it = b->M_unknowns.begin(), en = b->M_unknowns.end(); it != en; ++ it )
        {
            *it = unknown_type( new element_type( M_space ) );
        }

        return b;
    }


    //! return the curent order used at current time time
    int timeOrder() const
    {
        return M_order_cur;
    }

    //! return the order in time
    int bdfOrder() const
    {
        return M_order;
    }

    void setOrder( int order )
    {
        M_order = order;
    }

    //!return the prefix
    std::string bdfPrefix() const
    {
        return M_prefix;
    }

    //! return the number of iterations between order change
    int numberOfIterationsBetweenOrderChange() const
    {
        return M_iterations_between_order_change;
    }

    //! return the number of iterations since last order change
    int numberOfIterationsSinceOrderChange() const
    {
        return this->iteration()-M_last_iteration_since_order_change;
    }


    //! return a vector of the times prior to timeInitial() (included)
    std::map<int,double> priorTimes() const
        {
            std::map<int,double> prior;
            for( int i = 0; i < this->M_order; ++i )
                prior[i]=timeInitial()-i*timeStep();
            return prior;
        }

    /**
       Initialize all the entries of the unknown vector to be derived with the
       vector u0 (duplicated)
    */
    void initialize( element_type const& u0 );

    /**
       Initialize all the entries of the unknown vector to be derived with a
       set of vectors uv0
    */
    void initialize( unknowns_type const& uv0 );

    /**
       start the bdf
    */
    double start();
    double start( element_type const& u0 );
    double start( unknowns_type const& uv0 );

    /**
       restart the bdf
    */
    double restart();


    /**
       Update the vectors of the previous time steps by shifting on the right
       the old values.
       @param u_curr current (new) value of the state vector
    */
    template<typename container_type>
    void shiftRight( typename space_type::template Element<value_type, container_type> const& u_curr );

    double next() const
    {
        double tcur = super::next();

        if ( M_strategyHighOrderStart == 1 )
        {
            if ( ( ( M_iteration - M_last_iteration_since_order_change ) == M_iterations_between_order_change ) &&
                 M_order_cur < M_order )
            {
                M_last_iteration_since_order_change = M_iteration;
                ++M_order_cur;
            }
        }

        return tcur;
    }

    template<typename container_type>
    double
    next( typename space_type::template Element<value_type, container_type> const& u_curr )
    {
        this->shiftRight( u_curr );

        return this->next();
    }

		/**
		 * Return \f$ \alpha_i \f$
		 */
    double polyCoefficient( int i ) const
    {
        CHECK( i >=0 && i < BDF_MAX_ORDER-1 ) <<  "[BDF] invalid index " << i;
        return M_beta[this->timeOrder()-1][i];
    }  
		/**
		 * Return \f$ \frac{\alpha_i}{\Delta t} \f$
		 */
    double polyDerivCoefficient( int i ) const
    {
        CHECK( i >=0 && i <= BDF_MAX_ORDER ) << "[BDF] invalid index " << i;
        return M_alpha[this->timeOrder()-1][i]/math::abs( this->timeStep() );
        //return M_alpha[this->timeOrder()-1][i]/this->timeStep();
    }

    //! Returns the right hand side \f$ \bar{p} \f$ of the time derivative
    //! formula
    element_type polyDeriv() const;

    //! Compute the polynomial extrapolation approximation of order n-1 of
    //! u^{n+1} defined by the n stored state vectors
    element_type poly() const;

    //! Return a vector with the last n state vectors
    unknowns_type const& unknowns() const;

    //! Return a vector with the last n state vectors
    element_type& unknown( int i );

    template<typename container_type>
    void setUnknown( int i,  typename space_type::template Element<value_type, container_type> const& e )
    {
        *M_unknowns[i] = e;
    }

    void showMe( std::ostream& __out = std::cout ) const;

    void loadCurrent();

    void print() const
    {
        LOG(INFO) << "============================================================\n";
        LOG(INFO) << "BDF Information\n";
        LOG(INFO) << "   time step : " << this->timeStep() << "\n";
        LOG(INFO) << "time initial : " << this->timeInitial() << "\n";
        LOG(INFO) << "  time final : " << this->timeFinal() << "\n";
        LOG(INFO) << "  time order : " << this->timeOrder() << "\n";
    }


private:
    void init();


    void saveCurrent();

    //! save/load Bdf metadata
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
        DVLOG(2) << "[BDF::serialize] saving/loading archive\n";
        ar & boost::serialization::base_object<TSBase>( *this );
    }

    //! compute BDF coefficients
    void computeCoefficients();

private:

    //! bdf order
    int M_order;

    /**
     * strategy to start with high order scheme :
     * 0 start with order given
     * 1 start with order 1 and increase step by step up to order given
     */
    int M_strategyHighOrderStart;

    //! bdf order used at the current time
    mutable int M_order_cur;

    mutable int M_last_iteration_since_order_change;
    int M_iterations_between_order_change;

    //! space
    space_ptrtype M_space;

    //! Last n state vectors
    unknowns_type M_unknowns;

    //! Coefficients \f$ \alpha_i \f$ of the time bdf discretization
    std::vector<ublas::vector<double> > M_alpha;

    //! Coefficients \f$ \beta_i \f$ of the extrapolation
    std::vector<ublas::vector<double> > M_beta;

    std::string M_prefix;
};

template <typename SpaceType>
Bdf<SpaceType>::Bdf( po::variables_map const& vm,
                     space_ptrtype const& __space,
                     std::string const& name,
                     std::string const& prefix )
    :
    super( vm, name, prefix, __space->worldComm() ),
    M_order( vm[prefixvm( prefix, "bdf.order" )].as<int>() ),
    M_strategyHighOrderStart( vm[prefixvm( prefix, "bdf.strategy-high-order-start" )].as<int>() ),
    M_order_cur( M_order ),
    M_iterations_between_order_change( vm[prefixvm( prefix, "bdf.iterations-between-order-change" )].as<int>() ),
    M_space( __space ),
    M_alpha( BDF_MAX_ORDER ),
    M_beta( BDF_MAX_ORDER ),
    M_prefix( prefix )
{
    M_unknowns.resize( BDF_MAX_ORDER );

    for ( uint8_type __i = 0; __i < ( uint8_type )BDF_MAX_ORDER; ++__i )
    {
        M_unknowns[__i] = unknown_type( new element_type( M_space ) );
        M_unknowns[__i]->zero();
    }
    computeCoefficients();

}

template <typename SpaceType>
Bdf<SpaceType>::Bdf( space_ptrtype const& __space,
                     std::string const& name  )
    :
    super( name, __space->worldComm() ),
    M_order( 1 ),
    M_strategyHighOrderStart( 0 ),
    M_order_cur( 1 ),
    M_iterations_between_order_change( 1 ),
    M_space( __space ),
    M_alpha( BDF_MAX_ORDER ),
    M_beta( BDF_MAX_ORDER ),
    M_prefix( "" )
{
    M_unknowns.resize( BDF_MAX_ORDER );

    for ( uint8_type __i = 0; __i < ( uint8_type )BDF_MAX_ORDER; ++__i )
    {
        M_unknowns[__i] = unknown_type( new element_type( M_space ) );
        M_unknowns[__i]->zero();
    }
    computeCoefficients();
}

template <typename SpaceType>
void
Bdf<SpaceType>::computeCoefficients()
{
    for ( int i = 0; i < BDF_MAX_ORDER; ++i )
    {
        M_alpha[ i ].resize( i+2 );
        M_beta[ i ].resize( i+1 );
    }

    for ( int i = 0; i < BDF_MAX_ORDER; ++i )
    {
        if (  i == 0 ) // BDF_ORDER_ONE:
        {
            M_alpha[i][ 0 ] = 1.; // Backward Euler
            M_alpha[i][ 1 ] = 1.;
            M_beta[i][ 0 ] = 1.; // u^{n+1} \approx u^n
        }

        else if ( i == 1 ) // BDF_ORDER_TWO:
        {
            M_alpha[i][ 0 ] = 3. / 2.;
            M_alpha[i][ 1 ] = 2.;
            M_alpha[i][ 2 ] = -1. / 2.;
            M_beta[i][ 0 ] = 2.;
            M_beta[i][ 1 ] = -1.;
        }

        else if ( i == 2 ) // BDF_ORDER_THREE:
        {
            M_alpha[i][ 0 ] = 11. / 6.;
            M_alpha[i][ 1 ] = 3.;
            M_alpha[i][ 2 ] = -3. / 2.;
            M_alpha[i][ 3 ] = 1. / 3.;
            M_beta[i][ 0 ] = 3.;
            M_beta[i][ 1 ] = -3.;
            M_beta[i][ 2 ] = 1.;
        }

        else if ( i == 3 ) /// BDF_ORDER_FOUR:
        {
            M_alpha[i][ 0 ] = 25. / 12.;
            M_alpha[i][ 1 ] = 4.;
            M_alpha[i][ 2 ] = -3.;
            M_alpha[i][ 3 ] = 4. / 3.;
            M_alpha[i][ 4 ] = -1. / 4.;
            M_beta[i][ 0 ] = 4.;
            M_beta[i][ 1 ] = -6.;
            M_beta[i][ 2 ] = 4.;
            M_beta[i][ 3 ] = -1.;
        }
    }
}
template <typename SpaceType>
void
Bdf<SpaceType>::init()
{
    this->setPathSave( (boost::format("%3%bdf_o_%1%_dt_%2%")
                        %this->bdfOrder()
                        %this->timeStep()
                        %this->bdfPrefix()  ).str() );

    super::init();

    if ( this->isRestart() )
    {
        for ( int p = 0; p < std::min( M_order, M_iteration+1 ); ++p )
        {
            // create and open a character archive for output
            std::ostringstream ostr;

            if( M_rankProcInNameOfFiles )
                ostr << M_name << "-" << M_iteration-p<<"-proc"<<this->worldComm().globalRank()<<"on"<<this->worldComm().globalSize();
            else
                ostr << M_name << "-" << M_iteration-p;

            DVLOG(2) << "[Bdf::init()] load file: " << ostr.str() << "\n";

            fs::ifstream ifs;

            if ( this->restartPath().empty() ) ifs.open( this->path()/ostr.str() );

            else ifs.open( this->restartPath()/this->path()/ostr.str() );

            //fs::ifstream ifs (this->restartPath() / this->path() / ostr.str(), std::ios::binary);

            // load data from archive
            boost::archive::binary_iarchive ia( ifs );
            ia >> *M_unknowns[p];
        }
    }
}


template <typename SpaceType>
Bdf<SpaceType>::~Bdf()
{}


template <typename SpaceType>
void
Bdf<SpaceType>::initialize( element_type const& u0 )
{
    M_time_values_map.clear();
    std::ostringstream ostr;

    if( M_rankProcInNameOfFiles )
        ostr << M_name << "-" << 0<<"-proc"<<this->worldComm().globalRank()<<"on"<<this->worldComm().globalSize();
    else
        ostr << M_name << "-" << 0;
    //M_time_values_map.insert( std::make_pair( 0, boost::make_tuple( 0, ostr.str() ) ) );
    //M_time_values_map.push_back( 0 );
    M_time_values_map.push_back( M_Ti );
    std::for_each( M_unknowns.begin(), M_unknowns.end(), *boost::lambda::_1 = u0 );
    this->saveCurrent();
}

template <typename SpaceType>
void
Bdf<SpaceType>::initialize( unknowns_type const& uv0 )
{
    M_time_values_map.clear();
    std::ostringstream ostr;

    if( M_rankProcInNameOfFiles )
        ostr << M_name << "-" << 0<<"-proc"<<this->worldComm().globalRank()<<"on"<<this->worldComm().globalSize();
    else
        ostr << M_name << "-" << 0;
    //M_time_values_map.insert( std::make_pair( 0, boost::make_tuple( 0, ostr.str() ) ) );
    //M_time_values_map.push_back( 0);
    M_time_values_map.push_back( M_Ti );

    std::copy( uv0.begin(), uv0.end(), M_unknowns.begin() );
    this->saveCurrent();
}

template <typename SpaceType>
double
Bdf<SpaceType>::start()
{
    this->init();
    double ti = super::start();
    M_last_iteration_since_order_change = 1;
    switch ( M_strategyHighOrderStart )
    {
    default :
    case 0 : M_order_cur = M_order; break;
    case 1 : M_order_cur = 1; break;
    }
    return ti;
}

template <typename SpaceType>
double
Bdf<SpaceType>::start( element_type const& u0 )
{
    this->init();
    this->initialize( u0 );
    double ti = super::start();
    M_last_iteration_since_order_change = 1;
    switch ( M_strategyHighOrderStart )
    {
    default :
    case 0 : M_order_cur = M_order; break;
    case 1 : M_order_cur = 1; break;
    }
    return ti;
}

template <typename SpaceType>
double
Bdf<SpaceType>::start( unknowns_type const& uv0 )
{
    this->init();
    this->initialize( uv0 );
    double ti = super::start();
    M_last_iteration_since_order_change = 1;
    switch ( M_strategyHighOrderStart )
    {
    default :
    case 0 : M_order_cur = M_order; break;
    case 1 : M_order_cur = 1; break;
    }
    return ti;
}

template <typename SpaceType>
double
Bdf<SpaceType>::restart()
{
    this->init();

    double ti = super::restart();
    M_last_iteration_since_order_change = 1;

    switch ( M_strategyHighOrderStart )
    {
    default :
    case 0 :
    {
        M_order_cur = M_order;
    }
    break;
    case 1 :
    {
        M_order_cur = 1;
        for ( int i = 2; i<=M_iteration; ++i )
        {
            if ( ( ( i - M_last_iteration_since_order_change ) == M_iterations_between_order_change ) &&
                 M_order_cur < M_order )
            {
                M_last_iteration_since_order_change = i;
                ++M_order_cur;
            }
        }
    }
    break;
    }

    return ti;
}

template <typename SpaceType>
const
typename Bdf<SpaceType>::unknowns_type&
Bdf<SpaceType>::unknowns() const
{
    return M_unknowns;
}

template <typename SpaceType>
typename Bdf<SpaceType>::element_type&
Bdf<SpaceType>::unknown( int i )
{
    DVLOG(2) << "[Bdf::unknown] id: " << i << " l2norm = " << M_unknowns[i]->l2Norm() << "\n";
    return *M_unknowns[i];
}


template <typename SpaceType>
void
Bdf<SpaceType>::saveCurrent()
{
    if (!this->saveInFile()) return;

    bool doSave=false;
    for ( uint8_type i = 0; i < this->timeOrder() && !doSave; ++i )
        {
            int iterTranslate = M_iteration + this->timeOrder()-(i+1);
            if (iterTranslate % this->saveFreq()==0) doSave=true;
        }

    if (!doSave) return;

    TSBaseMetadata bdfsaver( *this );
    bdfsaver.save();

    {
        std::ostringstream ostr;

        if( M_rankProcInNameOfFiles )
            ostr << M_name << "-" << M_iteration<<"-proc"<<this->worldComm().globalRank()<<"on"<<this->worldComm().globalSize();
        else
            ostr << M_name << "-" << M_iteration;

        fs::ofstream ofs( M_path_save / ostr.str() );


        // load data from archive
        boost::archive::binary_oarchive oa( ofs );
        oa << *M_unknowns[0];
    }
}

template <typename SpaceType>
void
Bdf<SpaceType>::loadCurrent()
{
    //TSBaseMetadata bdfsaver( *this );
    //bdfsaver.save();

    {
        std::ostringstream ostr;

        if( M_rankProcInNameOfFiles )
            ostr << M_name << "-" << M_iteration<<"-proc"<<this->worldComm().globalRank()<<"on"<<this->worldComm().globalSize();
        else
            ostr << M_name << "-" << M_iteration;

        fs::ifstream ifs( M_path_save / ostr.str() );

        // load data from archive
        boost::archive::binary_iarchive ia( ifs );
        ia >> *M_unknowns[0];
    }
}

template <typename SpaceType>
template<typename container_type>
void
Bdf<SpaceType>::shiftRight( typename space_type::template Element<value_type, container_type> const& __new_unk )
{
    DVLOG(2) << "shiftRight: inserting time " << this->time() << "s\n";
    super::shiftRight();

    // shift all previously stored bdf data
    using namespace boost::lambda;
    typename unknowns_type::reverse_iterator __it = boost::next( M_unknowns.rbegin() );
    std::for_each( M_unknowns.rbegin(), boost::prior( M_unknowns.rend() ),
                   ( *lambda::_1 = *( *lambda::var( __it ) ), ++lambda::var( __it ) ) );
    // u(t^{n}) coefficient is in M_unknowns[0]
    *M_unknowns[0] = __new_unk;
    int i = 0;
    BOOST_FOREACH( boost::shared_ptr<element_type>& t, M_unknowns  )
    {
        DVLOG(2) << "[Bdf::shiftright] id: " << i << " l2norm = " << t->l2Norm() << "\n";
        ++i;
    }

    // save newly stored bdf data
    this->saveCurrent();
}


template <typename SpaceType>
typename Bdf<SpaceType>::element_type
Bdf<SpaceType>::polyDeriv() const
{
    element_type __t( M_space );
    __t.zero();

    FEELPP_ASSERT( __t.size() == M_space->nDof() )( __t.size() )( M_space->nDof() ).error( "invalid space element size" );
    FEELPP_ASSERT( __t.size() == M_unknowns[0]->size() )( __t.size() )( M_unknowns[0]->size() ).error( "invalid space element size" );

    for ( uint8_type i = 0; i < this->timeOrder(); ++i )
        __t.add( this->polyDerivCoefficient( i+1 ), *M_unknowns[i] );

    return __t;
}

template <typename SpaceType>
typename Bdf<SpaceType>::element_type
Bdf<SpaceType>::poly() const
{
    element_type __t( M_space );
    __t.zero();

    FEELPP_ASSERT( __t.size() == M_space->nDof() )( __t.size() )( M_space->nDof() ).error( "invalid space element size" );
    FEELPP_ASSERT( __t.size() == M_unknowns[0]->size() )( __t.size() )( M_unknowns[0]->size() ).error( "invalid space element size" );

    for ( uint8_type i = 0; i < this->timeOrder(); ++i )
        __t.add(  this->polyCoefficient( i ),  *M_unknowns[ i ] );

    return __t;
}



BOOST_PARAMETER_FUNCTION(
    ( boost::shared_ptr<Bdf<typename meta::remove_all<typename parameter::binding<Args, tag::space>::type>::type::element_type> > ),
    bdf, tag,
    ( required
      ( space,*( boost::is_convertible<mpl::_,boost::shared_ptr<Feel::FunctionSpaceBase> > ) ) )
    ( optional
      ( vm,*, Environment::vm() )
      ( prefix,*,"" )
      ( name,*,"bdf" )
      ( order,*( boost::is_integral<mpl::_> ),vm[prefixvm( prefix,"bdf.order" )].template as<int>() )
      ( initial_time,*( boost::is_floating_point<mpl::_> ),vm[prefixvm( prefix,"bdf.time-initial" )].template as<double>() )
      ( final_time,*( boost::is_floating_point<mpl::_> ),vm[prefixvm( prefix,"bdf.time-final" )].template as<double>() )
      ( time_step,*( boost::is_floating_point<mpl::_> ),vm[prefixvm( prefix,"bdf.time-step" )].template as<double>() )
      ( strategy,*( boost::is_integral<mpl::_> ),vm[prefixvm( prefix,"bdf.strategy" )].template as<int>() )
      ( steady,*( bool ),vm[prefixvm( prefix,"bdf.steady" )].template as<bool>() )
      ( restart,*( boost::is_integral<mpl::_> ),vm[prefixvm( prefix,"bdf.restart" )].template as<bool>() )
      ( restart_path,*,vm[prefixvm( prefix,"bdf.restart.path" )].template as<std::string>() )
      ( restart_at_last_save,*( boost::is_integral<mpl::_> ),vm[prefixvm( prefix,"bdf.restart.at-last-save" )].template as<bool>() )
      ( save,*( boost::is_integral<mpl::_> ),vm[prefixvm( prefix,"bdf.save" )].template as<bool>() )
      ( freq,*(boost::is_integral<mpl::_> ),vm[prefixvm( prefix,"bdf.save.freq" )].template as<int>() )
      ( rank_proc_in_files_name,*( boost::is_integral<mpl::_> ),vm[prefixvm( prefix,"bdf.rank-proc-in-files-name" )].template as<bool>() )
    ) )
{
    typedef typename meta::remove_all<space_type>::type::element_type _space_type;
    auto thebdf = boost::shared_ptr<Bdf<_space_type> >( new Bdf<_space_type>( vm,space,name,prefix ) );
    thebdf->setTimeInitial( initial_time );
    thebdf->setTimeFinal( final_time );
    thebdf->setTimeStep( time_step );
    thebdf->setOrder( order );
    thebdf->setSteady( steady );
    thebdf->setRestart( restart );
    thebdf->setRestartPath( restart_path );
    thebdf->setRestartAtLastSave( restart_at_last_save );
    thebdf->setSaveInFile( save );
    thebdf->setSaveFreq( freq );
    thebdf->setRankProcInNameOfFiles( rank_proc_in_files_name );
    return thebdf;
}


}
#endif
