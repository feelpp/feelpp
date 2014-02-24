/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
   Date: 2014-02-24

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
   \file newmark.hpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2014-02-24
*/
#ifndef _NEWMARK_H
#define _NEWMARK_H

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


/**
 * \class Newmark
 * \ingroup SpaceTime
 * \brief Newmark discretization
 *
 */
template<typename SpaceType>
class Newmark : public TSBase
{
    friend class boost::serialization::access;
    typedef TSBase super;
public:
    typedef Newmark<SpaceType> newmark_type;
    typedef boost::shared_ptr<newmark_type> newmark_ptrtype;
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
    Newmark( po::variables_map const& vm, space_ptrtype const& space, std::string const& name, std::string const& prefix="" );

    /**
     * Constructor
     * @param space approximation space
     * @param name name of the BDF
     */
    Newmark( space_ptrtype const& space, std::string const& name );

    //! copy operator
    Newmark( Newmark const& b )
        :
        super( b ),
        M_space( b.M_space ),
        M_unknowns( b.M_unknowns ),
        M_currentVel( b.M_currentVel ),
        M_currentAcc( b.M_currentAcc ),
        M_previousVel( b.M_previousVel ),
        M_previousAcc( b.M_previousAcc ),
        M_polyDeriv( b.M_polyDeriv ),
        M_gamma( b.M_gamma ),
        M_beta( b.M_beta )
    {}

    ~Newmark();

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

    double next() const { return super::next(); }

    template<typename container_type>
    double
    next( typename space_type::template Element<value_type, container_type> const& u_curr )
        {
            this->shiftRight( u_curr );
            return super::next();
        }

    template<typename container_type>
    void updateVelAcc( typename space_type::template Element<value_type, container_type> const& u_curr );


    double polyDerivCoefficient() const
    {
      return 1.0/(M_beta*std::pow(this->timeStep(),2));
    }

    //! Returns the right hand side \f$ \bar{p} \f$ of the time derivative
    //! formula
    element_type const& polyDeriv() const;

    //! Return a vector with the last n state vectors
    unknowns_type const& unknowns() const;

    //! Return a vector with the last n state vectors
    element_type& unknown( int i );

    element_type const& previousDisp( int i=0 ) const;
    element_type const& previousVel() const;
    element_type const& previousAcc() const;


    template<typename container_type>
    void setUnknown( int i,  typename space_type::template Element<value_type, container_type> const& e )
    {
        *M_unknowns[i] = e;
    }

    void showMe( std::ostream& __out = std::cout ) const;

    void loadCurrent();


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


    //! compute second derivative poly of rhs
    void computePolyDeriv();

private:

    //! space
    space_ptrtype M_space;

    //! Last n state vectors
    unknowns_type M_unknowns;

    unknown_type M_currentVel, M_currentAcc;
    unknown_type M_previousVel, M_previousAcc;
    unknown_type M_polyDeriv;

    //! Coefficients \f$ \gamma \f$ and \f$ \beta \f$  of the time newmark discretization
    double M_gamma, M_beta;

};

template <typename SpaceType>
Newmark<SpaceType>::Newmark( po::variables_map const& vm,
                             space_ptrtype const& __space,
                             std::string const& name,
                             std::string const& prefix )
    :
    super( vm, name, prefix, __space->worldComm() ),
    M_space( __space ),
    M_currentVel( unknown_type( new element_type( M_space ) ) ),
    M_currentAcc( unknown_type( new element_type( M_space ) ) ),
    M_previousVel( unknown_type( new element_type( M_space ) ) ),
    M_previousAcc( unknown_type( new element_type( M_space ) ) ),
    M_polyDeriv( unknown_type( new element_type( M_space ) ) ),
    M_gamma( 0.5 ),
    M_beta( 0.25 )
{
    M_unknowns.resize( 2 );

    for ( uint8_type __i = 0; __i < 2; ++__i )
    {
        M_unknowns[__i] = unknown_type( new element_type( M_space ) );
        //M_unknowns[__i]->zero();
    }

    //computeCoefficients();
}

template <typename SpaceType>
Newmark<SpaceType>::Newmark( space_ptrtype const& __space,
                             std::string const& name  )
    :
    super( name, __space->worldComm() ),
    M_space( __space ),
    M_currentVel( unknown_type( new element_type( M_space ) ) ),
    M_currentAcc( unknown_type( new element_type( M_space ) ) ),
    M_previousVel( unknown_type( new element_type( M_space ) ) ),
    M_previousAcc( unknown_type( new element_type( M_space ) ) ),
    M_polyDeriv( unknown_type( new element_type( M_space ) ) ),
    M_gamma( 0.5 ),
    M_beta( 0.25 )
{
    M_unknowns.resize( 2 );

    for ( uint8_type __i = 0; __i < 2; ++__i )
    {
        M_unknowns[__i] = unknown_type( new element_type( M_space ) );
        //M_unknowns[__i]->zero();
    }

    //computeCoefficients();
}

template <typename SpaceType>
void
Newmark<SpaceType>::computeCoefficients()
{
}

template <typename SpaceType>
void
Newmark<SpaceType>::init()
{
    super::init();

    if ( this->isRestart() )
    {
#warning TODO :  RESTART IN NEWMARK 
#if 0
        for ( int p = 0; p < std::min( M_order, M_iteration+1 ); ++p )
        {
            // create and open a character archive for output
            std::ostringstream ostr;

            if( M_rankProcInNameOfFiles )
                ostr << M_name << "-" << M_iteration-p<<"-proc"<<this->worldComm().globalRank()<<"on"<<this->worldComm().globalSize();
            else
                ostr << M_name << "-" << M_iteration-p;

            DVLOG(2) << "[Newmark::init()] load file: " << ostr.str() << "\n";

            fs::ifstream ifs;

            if ( this->restartPath().empty() ) ifs.open( this->path()/ostr.str() );

            else ifs.open( this->restartPath()/this->path()/ostr.str() );

            //fs::ifstream ifs (this->restartPath() / this->path() / ostr.str(), std::ios::binary);

            // load data from archive
            boost::archive::binary_iarchive ia( ifs );
            ia >> *M_unknowns[p];
        }
#endif
    }
}


template <typename SpaceType>
Newmark<SpaceType>::~Newmark()
{}


template <typename SpaceType>
void
Newmark<SpaceType>::initialize( element_type const& u0 )
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
Newmark<SpaceType>::initialize( unknowns_type const& uv0 )
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
Newmark<SpaceType>::start()
{
    this->init();

    return super::start();
}

template <typename SpaceType>
double
Newmark<SpaceType>::start( element_type const& u0 )
{
    this->init();
    this->initialize( u0 );
    auto res = super::start();
    return res;
}

template <typename SpaceType>
double
Newmark<SpaceType>::start( unknowns_type const& uv0 )
{
    this->init();
    this->initialize( uv0 );
    auto res = super::start();
    return res;
}

template <typename SpaceType>
double
Newmark<SpaceType>::restart()
{
    this->init();

    return super::restart();
}

template <typename SpaceType>
const
typename Newmark<SpaceType>::unknowns_type&
Newmark<SpaceType>::unknowns() const
{
    return M_unknowns;
}

template <typename SpaceType>
typename Newmark<SpaceType>::element_type&
Newmark<SpaceType>::unknown( int i )
{
    DVLOG(2) << "[Newmark::unknown] id: " << i << " l2norm = " << M_unknowns[i]->l2Norm() << "\n";
    return *M_unknowns[i];
}

template <typename SpaceType>
typename Newmark<SpaceType>::element_type const&
Newmark<SpaceType>::previousDisp( int i ) const
{
    return *M_unknowns[i];
}

template <typename SpaceType>
typename Newmark<SpaceType>::element_type const&
Newmark<SpaceType>::previousVel() const
{
  return *M_previousVel;
}

template <typename SpaceType>
typename Newmark<SpaceType>::element_type const&
Newmark<SpaceType>::previousAcc() const
{
  return *M_previousAcc;
}

template <typename SpaceType>
void
Newmark<SpaceType>::saveCurrent()
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
Newmark<SpaceType>::loadCurrent()
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
Newmark<SpaceType>::shiftRight( typename space_type::template Element<value_type, container_type> const& __new_unk )
{
    DVLOG(2) << "shiftRight: inserting time " << this->time() << "s\n";
    super::shiftRight();

    // update M_currentVel and M_currentAcc with new disp
    this->updateVelAcc(__new_unk);

    *M_previousVel=*M_currentVel;
    *M_previousAcc=*M_currentAcc;

    // shift all previously stored bdf data
#if 0
    using namespace boost::lambda;
    typename unknowns_type::reverse_iterator __it = boost::next( M_unknowns.rbegin() );
    std::for_each( M_unknowns.rbegin(), boost::prior( M_unknowns.rend() ),
                   ( *lambda::_1 = *( *lambda::var( __it ) ), ++lambda::var( __it ) ) );
    // u(t^{n}) coefficient is in M_unknowns[0]
#endif
    *M_unknowns[1] = *M_unknowns[0];
    *M_unknowns[0] = __new_unk;

    int i = 0;
    BOOST_FOREACH( boost::shared_ptr<element_type>& t, M_unknowns  )
    {
        DVLOG(2) << "[Newmark::shiftright] id: " << i << " l2norm = " << t->l2Norm() << "\n";
        ++i;
    }

    // save newly stored bdf data
    this->saveCurrent();

    // compute second derivative poly of rhs
    this->computePolyDeriv();
}


template <typename SpaceType>
typename Newmark<SpaceType>::element_type const&
Newmark<SpaceType>::polyDeriv() const
{
    return *M_polyDeriv;
}

template <typename SpaceType>
void
Newmark<SpaceType>::computePolyDeriv()
{
    double coef_disp = 1.0/(M_beta*std::pow(this->timeStep(),2));
    double coef_vel  = 1.0/(M_beta*this->timeStep());
    double coef_acc  = 1.0/(2*M_beta) - 1.0;
    M_polyDeriv->zero();
    M_polyDeriv->add(coef_disp, this->previousDisp());
    M_polyDeriv->add(coef_vel, this->previousVel());
    M_polyDeriv->add(coef_acc, this->previousAcc());
}

template <typename SpaceType>
template<typename container_type>
void
Newmark<SpaceType>::updateVelAcc( typename space_type::template Element<value_type, container_type> const& __new_unk )
{
  double deltaT = this->timeStep();

  double cst_vel_displ = M_gamma/(M_beta*deltaT);
  double cst_vel_vel = (M_gamma/M_beta)-1;
  double cst_vel_acc = deltaT*((M_gamma/(2*M_beta))-1 );
  double cst_acc_displ = 1./(M_beta*std::pow(deltaT,2));
  double cst_acc_vel = 1./(M_beta*deltaT);
  double cst_acc_acc = ((1./(2*M_beta))-1);

  M_currentVel->zero();
  M_currentVel->add(cst_vel_displ, __new_unk );
  M_currentVel->add(-cst_vel_displ, this->previousDisp());
  M_currentVel->add(-cst_vel_vel, this->previousVel());
  M_currentVel->add(-cst_vel_acc, this->previousAcc());

  M_currentAcc->zero();
  M_currentAcc->add(cst_acc_displ, __new_unk );
  M_currentAcc->add(-cst_acc_displ, this->previousDisp());
  M_currentAcc->add(-cst_acc_vel, this->previousVel());
  M_currentAcc->add(-cst_acc_acc, this->previousAcc());
}



BOOST_PARAMETER_FUNCTION(
    ( boost::shared_ptr<Newmark<typename meta::remove_all<typename parameter::binding<Args, tag::space>::type>::type::element_type> > ),
    newmark, tag,
    ( required
      ( space,*( boost::is_convertible<mpl::_,boost::shared_ptr<Feel::FunctionSpaceBase> > ) ) )
    ( optional
      ( vm,*, Environment::vm() )
      ( prefix,*,"" )
      ( name,*,"newmark" )
      //( order,*( boost::is_integral<mpl::_> ),vm[prefixvm( prefix,"bdf.order" )].template as<int>() )
      ( initial_time,*( boost::is_floating_point<mpl::_> ),vm[prefixvm( prefix,"ts.time-initial" )].template as<double>() )
      ( final_time,*( boost::is_floating_point<mpl::_> ),vm[prefixvm( prefix,"ts.time-final" )].template as<double>() )
      ( time_step,*( boost::is_floating_point<mpl::_> ),vm[prefixvm( prefix,"ts.time-step" )].template as<double>() )
      //( strategy,*( boost::is_integral<mpl::_> ),vm[prefixvm( prefix,"ts.strategy" )].template as<int>() )
      ( steady,*( bool ),vm[prefixvm( prefix,"ts.steady" )].template as<bool>() )
      ( restart,*( boost::is_integral<mpl::_> ),vm[prefixvm( prefix,"ts.restart" )].template as<bool>() )
      ( restart_path,*,vm[prefixvm( prefix,"ts.restart.path" )].template as<std::string>() )
      ( restart_at_last_save,*( boost::is_integral<mpl::_> ),vm[prefixvm( prefix,"ts.restart.at-last-save" )].template as<bool>() )
      ( save,*( boost::is_integral<mpl::_> ),vm[prefixvm( prefix,"ts.save" )].template as<bool>() )
      ( freq,*(boost::is_integral<mpl::_> ),vm[prefixvm( prefix,"ts.save.freq" )].template as<int>() )
      ( rank_proc_in_files_name,*( boost::is_integral<mpl::_> ),vm[prefixvm( prefix,"ts.rank-proc-in-files-name" )].template as<bool>() )
    ) )
{
    typedef typename meta::remove_all<space_type>::type::element_type _space_type;
    auto thenewmark = boost::shared_ptr<Newmark<_space_type> >( new Newmark<_space_type>( vm,space,name,prefix ) );
    thenewmark->setTimeInitial( initial_time );
    thenewmark->setTimeFinal( final_time );
    thenewmark->setTimeStep( time_step );
    thenewmark->setOrder( 2/*order*/ );
    thenewmark->setSteady( steady );
    thenewmark->setRestart( restart );
    thenewmark->setRestartPath( restart_path );
    thenewmark->setRestartAtLastSave( restart_at_last_save );
    thenewmark->setSaveInFile( save );
    thenewmark->setSaveFreq( freq );
    thenewmark->setRankProcInNameOfFiles( rank_proc_in_files_name );
    return thenewmark;
}


}
#endif
