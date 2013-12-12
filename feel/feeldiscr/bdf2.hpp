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
#include <feel/feeldiscr/functionspace.hpp>

namespace Feel
{
namespace ublas = boost::numeric::ublas;
namespace fs = boost::filesystem;

enum BDFState { BDF_UNITIALIZED = 0, BDF_RUNNING, BDF_STOPPED };

enum BDFTimeScheme { BDF_ORDER_ONE=1, BDF_ORDER_TWO, BDF_ORDER_THREE, BDF_ORDER_FOUR, BDF_MAX_ORDER = 4 };
enum BDFStragegy { BDF_STRATEGY_DT_CONSTANT,
                   BDF_STRATEGY_DT_ADAPTATIVE
                 };

class BdfBase
{
    friend class boost::serialization::access;
public:
    typedef std::vector<double>::iterator time_iterator;
    typedef std::vector<double>::const_iterator time_const_iterator;
    typedef std::vector<double> time_values_map_type;

    BdfBase()
        :
        M_order( 1 ),
        M_order_cur( 1 ),
        M_name( "bdf" ),
        M_time( 0.0 ),
        M_iteration( 0 ),
        M_Ti( 0.0 ),
        M_Tf( 1.0 ),
        M_dt( 1.0 ),
        M_strategy( BDF_STRATEGY_DT_CONSTANT ),
        M_state( BDF_UNITIALIZED ),
        M_n_restart( 0 ),
        M_restart( false ),
        M_restartPath( "" ),
        M_restartAtLastSave( false ),
        M_alpha( BDF_MAX_ORDER ),
        M_beta( BDF_MAX_ORDER ),
        M_saveInFile( true ),
        M_saveFreq( 1 ),
        M_rankProcInNameOfFiles( false ),
        M_worldComm( Environment::worldComm() )
    {}

#if 0
    template <class ArgumentPack>
    BdfBase( ArgumentPack const& args )
        :
        M_order( args[_order | 1] ),
        M_name( args[_name | "bdf"] ),
        M_time( args[_initial_time | 0] ),
        M_iteration( 0 ),
        M_Ti( args[_initial_time | 0] ),
        M_Tf( args[_final_time | 1] ),
        M_dt( args[_time_step | 0.1] ),
        M_iterations_between_order_change( 1 ),
        M_strategy( ( BDFStragegy )args[_strategy | BDF_STRATEGY_DT_CONSTANT] ),
        M_state( BDF_UNITIALIZED ),
        M_n_restart( 0 ),
        M_restart( args[_restart | false] ),
        M_restartPath( args[_restart_path | ""] ),
        M_restartAtLastSave( args[_restart_at_last_save | false] )
        M_alpha( BDF_MAX_ORDER ),
        M_beta( BDF_MAX_ORDER ),
        M_saveInFile( true ),
        M_saveFreq( 1 ),
        M_rankProcInNameOfFiles( false ),
        M_worldComm( Environment::worldComm() )
    {

    }
#endif
    BdfBase( po::variables_map const& vm, std::string name, std::string const& prefix, WorldComm const& worldComm )
        :
        M_order( vm[prefixvm( prefix, "bdf.order" )].as<int>() ),
        M_order_cur( M_order ),
        M_name( name ),
        M_time( vm[prefixvm( prefix, "bdf.time-initial" )].as<double>() ),
        M_iteration( 0 ),
        M_iterations_between_order_change( vm[prefixvm( prefix, "bdf.iterations-between-order-change" )].as<int>() ),
        M_Ti( vm[prefixvm( prefix, "bdf.time-initial" )].as<double>() ),
        M_Tf( vm[prefixvm( prefix, "bdf.time-final" )].as<double>() ),
        M_dt( vm[prefixvm( prefix, "bdf.time-step" )].as<double>() ),
        M_strategy( ( BDFStragegy )vm[prefixvm( prefix, "bdf.strategy" )].as<int>() ),
        M_state( BDF_UNITIALIZED ),
        M_n_restart( 0 ),
        M_restart( vm[prefixvm( prefix, "bdf.restart" )].as<bool>() ),
        M_restartPath( vm[prefixvm( prefix, "bdf.restart.path" )].as<std::string>() ),
        M_restartAtLastSave( vm[prefixvm( prefix, "bdf.restart.at-last-save" )].as<bool>() ),
        M_alpha( BDF_MAX_ORDER ),
        M_beta( BDF_MAX_ORDER ),
        M_saveInFile( vm[prefixvm( prefix, "bdf.save" )].as<bool>() ),
        M_saveFreq( vm[prefixvm( prefix, "bdf.save.freq" )].as<int>() ),
        M_rankProcInNameOfFiles( vm[prefixvm( prefix, "bdf.rank-proc-in-files-name" )].as<bool>() ),
        M_worldComm( worldComm )
    {
    }
    BdfBase( std::string name, WorldComm const& worldComm )
        :
        M_order( 1 ),
        M_order_cur( 1 ),
        M_name( name ),
        M_time( 0. ),
        M_iteration( 0 ),
        M_iterations_between_order_change( 1 ),
        M_Ti( 0. ),
        M_Tf( 1.0 ),
        M_dt( 1.0 ),
        M_strategy( ( BDFStragegy )0 ),
        M_state( BDF_UNITIALIZED ),
        M_n_restart( 0 ),
        M_restart( false ),
        M_restartPath( "" ),
        M_restartAtLastSave( false ),
        M_alpha( BDF_MAX_ORDER ),
        M_beta( BDF_MAX_ORDER ),
        M_saveInFile( true ),
        M_saveFreq( 1 ),
        M_worldComm( worldComm )
    {
    }

    BdfBase( BdfBase const& b )
        :
        M_order( b.M_order ),
        M_order_cur( b.M_order_cur ),
        M_name( b.M_name ),
        M_time( b.M_time ),
        M_iteration( b.M_iteration ),
        M_iterations_between_order_change( b.M_iterations_between_order_change ),
        M_Ti( b.M_Ti ),
        M_Tf( b.M_Tf ),
        M_dt( b.M_dt ),
        M_strategy( b.M_strategy ),
        M_state( b.M_state ),
        M_n_restart( b.M_n_restart ),
        M_restart( b.M_restart ),
        M_restartPath( b.M_restartPath ),
        M_restartAtLastSave( b.M_restartAtLastSave ),
        M_time_values_map( b.M_time_values_map ),
        M_alpha( b.M_alpha ),
        M_beta( b.M_beta ),
        M_saveInFile( b.M_saveInFile ),
        M_saveFreq( b.M_saveFreq ),
        M_rankProcInNameOfFiles( b.M_rankProcInNameOfFiles ),
        M_worldComm( b.M_worldComm )
    {}

    virtual ~BdfBase() {}

    double polyCoefficient( int i ) const
    {
        FEELPP_ASSERT( i >=0 && i < BDF_MAX_ORDER-1 ).error( "[BDF] invalid index" );
        return M_beta[this->timeOrder()-1][i];
    }
    double polyDerivCoefficient( int i ) const
    {
        FEELPP_ASSERT( i >=0 && i <= BDF_MAX_ORDER ).error( "[BDF] invalid index" );
        return M_alpha[this->timeOrder()-1][i]/math::abs( this->timeStep() );
        //return M_alpha[this->timeOrder()-1][i]/this->timeStep();
    }

    BdfBase& operator=( BdfBase const& b )
    {
        if ( this != &b )
        {
            M_order = b.M_order;
            M_order_cur = b.M_order_cur;
            M_name = b.M_name;
            M_time = b.M_time;
            M_iteration = b.M_iteration;
            M_Ti = b.M_Ti;
            M_Tf = b.M_Tf;
            M_dt = b.M_dt;
            M_iterations_between_order_change = b.M_iterations_between_order_change;
            M_n_restart = b.M_n_restart;
            M_restart = b.M_restart;
            M_restartPath = b.M_restartPath;
            M_restartAtLastSave = b.M_restartAtLastSave;
            M_strategy = b.M_strategy;
            M_state = b.M_state;

            M_alpha = b.M_alpha;
            M_beta = b.M_beta;
            M_saveInFile = b.M_saveInFile;
            M_rankProcInNameOfFiles = b.M_rankProcInNameOfFiles;

            M_time_values_map = b.M_time_values_map;
            M_worldComm = b.M_worldComm;
        }

        return *this;
    }

    //! save/load Bdf metadata
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
        DVLOG(2) << "[BDF::serialize] serialize BDFBase\n";
#if 0
        ar & M_order;
        ar & M_name;
        ar & M_time;

        ar & M_n_restart;

        ar & M_Tf;
#endif
        //ar & M_time_orders;
        ar & boost::serialization::make_nvp( "time_values", M_time_values_map );

        DVLOG(2) << "[BDF::serialize] time orders size: " << M_time_orders.size() << "\n";
        DVLOG(2) << "[BDF::serialize] time values size: " << M_time_values_map.size() << "\n";

        for ( auto it = M_time_values_map.begin(), en = M_time_values_map.end(); it!=en; ++it )
        {
            //LOG(INFO) << "[Bdf] order " << i << "=" << M_time_orders[i] << "\n";
            DVLOG(2) << "[Bdf::serialize] value " << *it << "\n";

        }

        DVLOG(2) << "[BDF::serialize] serialize BDFBase done\n";
    }

    //! return the order in time
    int timeOrder() const
    {
        return M_order_cur;
    }

    //! return the initial time
    double timeInitial() const
    {
        return M_Ti;
    }

    //! return the final time
    double timeFinal() const
    {
        return M_Tf;
    }

    //! return the final step
    double timeStep() const
    {
        return M_dt;
    }

    //! return the BDF strategy
    BDFStragegy strategy() const
    {
        return M_strategy;
    }

    //! return the number of iterations between order change
    int numberOfIterationsBetweenOrderChange() const
    {
        return M_iterations_between_order_change;
    }

    //! return the number of iterations since last order change
    int numberOfIterationsSinceOrderChange() const
    {
        return M_iteration-M_last_iteration_since_order_change;
    }

    //! return the number of restarts
    int nRestart() const
    {
        return M_n_restart;
    }

    //! return the value of the bool restart
    bool isRestart() const
    {
        return M_restart;
    }

    //! return value of do restart at last save
    bool doRestartAtLastSave() const
    {
        return M_restartAtLastSave;
    }

    //! return the current time
    double time() const
    {
        return M_time;
    }

    //! return the iteration number
    int iteration() const
    {
        return M_iteration;
    }

    //! return the real time in seconds spent in the iteration
    double realTimePerIteration() const
    {
        FEELPP_ASSERT( state() == BDF_RUNNING ).error( "invalid BDF state" );
        M_real_time_per_iteration = M_timer.elapsed();
        return M_real_time_per_iteration;
    }

    //! return vector of time
    time_values_map_type timeValues() const
    {
        return M_time_values_map;
    }

    //! start the bdf
    double start()
    {
        M_state = BDF_RUNNING;
        M_timer.restart();
        M_iteration = 1;
        M_time = M_Ti+this->timeStep();
        // warning: this needs to be fixed wrt restart
        M_last_iteration_since_order_change = 1;
        M_order_cur = 1;
        return M_Ti;
    }

    //! restart the bdf
    double restart()
    {
        M_state = BDF_RUNNING;
        M_timer.restart();
        M_time = M_Ti+this->timeStep();
        M_last_iteration_since_order_change = 1;
        M_order_cur = 1;
        ++M_iteration;

        for ( int i = 2; i<=M_iteration; ++i )
        {
            if ( ( ( i - M_last_iteration_since_order_change ) == M_iterations_between_order_change )&&
                    M_order_cur < M_order )
            {
                M_last_iteration_since_order_change = i;
                ++M_order_cur;
            }
        }

        return M_Ti;
    }


    //! return true if Bdf is finished, false otherwise
    bool isFinished() const
    {
        bool finished=false;

        if ( M_Tf>M_Ti )
        {
            if ( M_time > M_Tf )
            {
                M_state = BDF_STOPPED;
                finished=true;
            }
        }

        else
        {
            if ( M_time < M_Tf )
            {
                M_state = BDF_STOPPED;
                finished=true;
            }
        }

        return finished;
    }

    /**
     * advance in time
     *
     * \todo implements the strategies here, at the moment constant
     * time step
     */
    virtual double next() const
    {
        FEELPP_ASSERT( state() == BDF_RUNNING ).error( "invalid BDF state" );
        M_real_time_per_iteration = M_timer.elapsed();
        M_timer.restart();
        M_time += M_dt;
        ++M_iteration;

        if ( ( ( M_iteration - M_last_iteration_since_order_change ) == M_iterations_between_order_change )&&
                M_order_cur < M_order )
        {
            M_last_iteration_since_order_change = M_iteration;
            ++M_order_cur;
        }

        return M_time;
    }

    virtual void shiftRight()
    {
        // create and open a character archive for output
        std::ostringstream ostr;

        if( M_rankProcInNameOfFiles )
            ostr << M_name << "-" << M_iteration<<"-proc"<< this->worldComm().globalRank()<<"on"<<this->worldComm().globalSize();
        else
            ostr << M_name << "-" << M_iteration;

        DVLOG(2) << "[BdfBase::shiftRight] solution name " << ostr.str() << "\n";

        //M_time_values_map.insert( std::make_pair( M_iteration, this->time() ) );
        M_time_values_map.push_back( this->time() );
        //update();
    }
    virtual void update()
    {
        double tn = M_time_values_map[M_iteration];
        double tn1 = M_time_values_map[M_iteration-1];
        double tn2= 0;

        if ( M_iteration >= 2 )
            tn2 = M_time_values_map[M_iteration-2];

        // order 3 and 4 not yet done
#if 0
        double tn3= 0;

        if ( M_iteration >= 3 )
            tn3 = M_time_values_map[M_iteration-3];

#endif

        for ( int i = 0; i < BDF_MAX_ORDER; ++i )
        {
            if (  i == 0 ) // BDF_ORDER_ONE:
            {
                M_alpha[i][ 0 ] = 1./( tn-tn1 ); // Backward Euler
                M_alpha[i][ 1 ] = 1./( tn-tn1 );
                M_beta[i][ 0 ] = 1.; // u^{n+1} \approx u^n
            }

            if (  i == 1 ) // BDF_ORDER_TWO:
            {
                double denom = ( tn-tn1 )*( tn1-tn2 );
                M_alpha[i][ 0 ] = 0.5*( tn1-2*tn2+tn )/denom;
                M_alpha[i][ 1 ] = ( tn2-tn )/denom;
                M_alpha[i][ 2 ] = 0.5/( tn1-tn2 );
                M_beta[i][ 0 ] = 1+denom;
                M_beta[i][ 1 ] = denom;
            }
        }
    }
    //! return the state of the BDF
    BDFState state() const
    {
        return M_state;
    }

    //! return the relative path where the bdf data is stored
    fs::path path()
    {
        return M_path_save;
    }
    fs::path restartPath()
    {
        return M_restartPath;
    }

    bool saveInFile() const
    {
        return M_saveInFile;
    }

    int saveFreq() const
    {
        return M_saveFreq;
    }

    bool rankProcInNameOfFiles() const
    {
        return  M_rankProcInNameOfFiles ;
    }

    WorldComm const& worldComm() const
    {
        return M_worldComm;
    }

    void setOrder( int order )
    {
        M_order = order;
    }
    void setTimeInitial( double ti )
    {
        M_Ti = ti;
    }
    void setTimeStep( double dt )
    {
        M_dt = dt;
    }
    void setTimeFinal( double T )
    {
        M_Tf = T;
    }
    void setStrategy( int strategy )
    {
        M_strategy = ( BDFStragegy )strategy;
    }
    void setSteady( bool steady = true )
    {
        if ( steady )
        {
            M_dt=1e30;
            M_Tf=1e30;
        }
    }
    void setRestart( bool doRestart )
    {
        M_restart=doRestart;
    }
    void setRestartPath( std::string s )
    {
        M_restartPath=s;
    }
    void setRestartAtLastSave( bool b )
    {
        M_restartAtLastSave=b;
    }
    void setSaveInFile( bool b )
    {
        M_saveInFile = b;
    }
    void setSaveFreq(int v)
    {
        M_saveFreq=v;
    }
    void setRankProcInNameOfFiles( bool b )
    {
        M_rankProcInNameOfFiles = b;
    }

    void print() const
    {
        LOG(INFO) << "============================================================\n";
        LOG(INFO) << "BDF Information\n";
        LOG(INFO) << "   time step : " << this->timeStep() << "\n";
        LOG(INFO) << "time initial : " << this->timeInitial() << "\n";
        LOG(INFO) << "  time final : " << this->timeFinal() << "\n";
        LOG(INFO) << "  time order : " << this->timeOrder() << "\n";
    }
protected:
    //! time order
    int M_order;
    mutable int M_order_cur;

    //! name of the file holding the bdf data
    std::string M_name;

    //! time
    mutable double M_time;
    mutable int M_iteration;
    mutable int M_last_iteration_since_order_change;
    int M_iterations_between_order_change;

    //! initial time to start
    double M_Ti;

    //! final time
    double M_Tf;

    //! timestep
    double M_dt;

    //! BDF strategy  (constant timestep or adaptative timestep)
    BDFStragegy M_strategy;

    //! Bdf state
    mutable BDFState M_state;

    int M_n_restart;

    //! do a restart
    bool M_restart;

    //! restart path
    fs::path M_restartPath;

    //! do restart with ti the last save
    bool M_restartAtLastSave;

    //! timer for real time per iteration
    mutable boost::timer M_timer;

    //! real time spent per iteration (in seconds)
    mutable double M_real_time_per_iteration;

    //! vector that holds the time order at each bdf step
    std::vector<int> M_time_orders;

    //! vector that holds the time value at each bdf step
    time_values_map_type M_time_values_map;

    //! path to the directory to store the functions
    fs::path M_path_save;

    //! Coefficients \f$ \alpha_i \f$ of the time bdf discretization
    std::vector<ublas::vector<double> > M_alpha;

    //! Coefficients \f$ \beta_i \f$ of the extrapolation
    std::vector<ublas::vector<double> > M_beta;

    //! Save solutions in file after each step
    bool M_saveInFile;

    //! frequence of save solutions
    int M_saveFreq;

    //! put the rank of the processor in generated files
    bool M_rankProcInNameOfFiles;

    //!  mpi communicator tool
    WorldComm M_worldComm;

protected:
    void init()
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

        std::ostringstream ostr;
        ostr << "bdf_o_" << M_order << "_dt_" << M_dt;
        //ostr << "bdf/" << M_name << "/o_" << M_order << "/dt_" << M_dt;
        M_path_save = ostr.str();

        // if directory does not exist, create it
        if ( !fs::exists( M_path_save ) && this->saveInFile() && this->worldComm().isMasterRank() )
            fs::create_directories( M_path_save );
        // be sure that all process can find the path after
        this->worldComm().barrier();

        if ( M_restart )
        {
            //fs::ifstream ifs;
            fs::path thepath;

            if ( this->restartPath().empty() ) thepath=this->path()/"metadata";//   ifs.open(this->path()/"metadata");

            else thepath=this->restartPath()/this->path()/"metadata"; //ifs.open(this->restartPath()/this->path()/"metadata");

            // read the saved bdf data
            if ( fs::exists( thepath/* this->restartPath() / this->path() / "metadata" )*/ ) )
            {
                DVLOG(2) << "[Bdf] loading metadata from " << M_path_save.string() << "\n";

                //fs::ifstream ifs( this->restartPath() / this->path() / "metadata")
                fs::ifstream ifs( thepath );


                boost::archive::text_iarchive ia( ifs );
                ia >> BOOST_SERIALIZATION_NVP( *this );
                DVLOG(2) << "[Bdf::init()] metadata loaded\n";
                //BdfBaseMetadata bdfloader( *this );
                //bdfloader.load();

                // modify Ti with last saved time
                if (this->doRestartAtLastSave()) M_Ti = M_time_values_map.back();

                M_iteration = 0;
                // look for M_ti in the time values
                bool found = false;
                BOOST_FOREACH( auto time, M_time_values_map )
                {
                    if ( math::abs( time-M_Ti ) < 1e-10 )
                    {
                        //M_iteration = time.first;
                        //std::cout << "time found " << time << std::endl;
                        found = true;
                        break;
                    }

                    ++M_iteration;
                }
                //std::cout << "M_iteration " << M_iteration << std::endl;

                if ( !found )
                {
                    DVLOG(2) << "[Bdf] intial time " << M_Ti << " not found\n";
                    M_Ti = 0.0;
                    M_iteration = 0;
                    M_time_values_map.clear();
                    return;
                }

                else
                {
                    if (this->saveFreq()==1)
                    {
                        M_time_values_map.resize( M_iteration+1 );
                    }
                    else
                    {
                        int nItBack = M_iteration % this->saveFreq();
                        M_iteration-=nItBack;
                        M_time_values_map.resize( M_iteration+1 );
                        M_Ti = M_time_values_map.back();
                    }
                }

                DVLOG(2) << "[Bdf] initial time is Ti=" << M_Ti << "\n";

                DVLOG(2) << "[Bdf::init()] file index: " << M_iteration << "\n";
            }

            else
            {
                M_Ti = 0.0;
                M_time_values_map.clear();
            }
        }

    } // init
};
class BdfBaseMetadata
{
public:
    BdfBaseMetadata( BdfBase& bdf )
        :
        M_bdf( bdf )
    {
    }

    void load()
    {
        fs::ifstream ifs;

        if ( M_bdf.restartPath().empty() ) ifs.open( M_bdf.path()/"metadata" );

        else ifs.open( M_bdf.restartPath()/M_bdf.path()/"metadata" );

        //fs::ifstream ifs( M_bdf.path() / "metadata");

        boost::archive::text_iarchive ia( ifs );
        ia >> BOOST_SERIALIZATION_NVP( M_bdf );
        DVLOG(2) << "[Bdf::init()] metadata loaded\n";
    }

    void save()
    {
        if ( !M_bdf.saveInFile() ) return;

        // only master process write
        if ( M_bdf.worldComm().isMasterRank() )
        {
            fs::ofstream ofs( M_bdf.path() / "metadata" );

            boost::archive::text_oarchive oa( ofs );
            oa << BOOST_SERIALIZATION_NVP( ( BdfBase const& )M_bdf );
            DVLOG(2) << "[Bdf::init()] metadata saved\n";
        }
        // to be sure that all process can read the metadata file
        M_bdf.worldComm().barrier();

    }

private:
    BdfBase& M_bdf;
};
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
 * (t_i,u_i) for i = k-n+1,...,k+1.
 *
 * The approximative time derivative \f$ p'(t_{k+1}) \f$ is a linear
 * combination of state vectors u_i:
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
 * calculate \f$ \bar{p} \f$. It also provides alpha_i
 * and can extrapolate the the new state from the n last states with a
 * polynomial of order n-1:
 *
 * \f$ u_{k+1} \approx \sum_{i=0}^{n-1} \beta_i u_{k-i} \f$
 */
template<typename SpaceType>
class Bdf : public BdfBase
{
    friend class boost::serialization::access;
    typedef BdfBase super;
public:
    typedef Bdf<SpaceType> bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;
    typedef SpaceType space_type;
    typedef boost::shared_ptr<space_type>  space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef typename space_type::return_type return_type;
    typedef typename element_type::value_type value_type;
    //typedef boost::numeric::ublas::vector< element_type > unknowns_type;
    typedef boost::shared_ptr<element_type> unknown_type;
    typedef std::vector< unknown_type > unknowns_type;
    typedef typename node<value_type>::type node_type;

    typedef typename super::time_iterator time_iterator;
#if 0
    BOOST_PARAMETER_CONSTRUCTOR(
        Bdf, ( BdfBase ), tag,
        ( required ( final_time,* ) )
        ( optional
          ( name,* )
          ( order,* )
          ( initial_time,* )
          ( time_step,* )
          ( strategy,* ) ) )
#endif
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
        M_space( b.M_space ),
        M_unknowns( b.M_unknowns )
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
        M_unknowns[i]->assign( e );
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
        ar & boost::serialization::base_object<BdfBase>( *this );
    }

private:

    //! space
    space_ptrtype M_space;

    //! Last n state vectors
    unknowns_type M_unknowns;
};

template <typename SpaceType>
Bdf<SpaceType>::Bdf( po::variables_map const& vm,
                     space_ptrtype const& __space,
                     std::string const& name,
                     std::string const& prefix )
    :
    super( vm, name, prefix, __space->worldComm() ),
    M_space( __space )
{
    M_unknowns.resize( BDF_MAX_ORDER );

    for ( uint8_type __i = 0; __i < ( uint8_type )BDF_MAX_ORDER; ++__i )
    {
        M_unknowns[__i] = unknown_type( new element_type( M_space ) );
        M_unknowns[__i]->zero();
    }

}

template <typename SpaceType>
Bdf<SpaceType>::Bdf( space_ptrtype const& __space,
                     std::string const& name  )
    :
    super( name, __space->worldComm() ),
    M_space( __space )
{
    M_unknowns.resize( BDF_MAX_ORDER );

    for ( uint8_type __i = 0; __i < ( uint8_type )BDF_MAX_ORDER; ++__i )
    {
        M_unknowns[__i] = unknown_type( new element_type( M_space ) );
        M_unknowns[__i]->zero();
    }

}

template <typename SpaceType>
void
Bdf<SpaceType>::init()
{
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

    return super::start();
}

template <typename SpaceType>
double
Bdf<SpaceType>::start( element_type const& u0 )
{
    this->init();
    this->initialize( u0 );
    auto res = super::start();
    return res;
}

template <typename SpaceType>
double
Bdf<SpaceType>::start( unknowns_type const& uv0 )
{
    this->init();
    this->initialize( uv0 );
    auto res = super::start();
    return res;
}

template <typename SpaceType>
double
Bdf<SpaceType>::restart()
{
    this->init();

    return super::restart();
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

    BdfBaseMetadata bdfsaver( *this );
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
    //BdfBaseMetadata bdfsaver( *this );
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
    thebdf->setStrategy( strategy );
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
