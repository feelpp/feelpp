/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-12-20

  Copyright (C) 2013 Feel++ Consortium

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
   \file tsbase.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-20
 */
#ifndef FEELPP_TSBASE_HPP
#define FEELPP_TSBASE_HPP 1

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

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/typetraits.hpp>


namespace Feel
{
namespace ublas = boost::numeric::ublas;
namespace fs = boost::filesystem;

enum TSState { TS_UNITIALIZED = 0, TS_RUNNING, TS_STOPPED };
enum TSStragegy { TS_STRATEGY_DT_CONSTANT,TS_STRATEGY_DT_ADAPTATIVE};

class TSBase
{
    friend class boost::serialization::access;
public:
    typedef std::vector<double>::iterator time_iterator;
    typedef std::vector<double>::const_iterator time_const_iterator;
    typedef std::vector<double> time_values_map_type;

    TSBase()
        :
        M_order( 1 ),
        M_order_cur( 1 ),
        M_name( "bdf" ),
        M_time( 0.0 ),
        M_iteration( 0 ),
        M_Ti( 0.0 ),
        M_Tf( 1.0 ),
        M_dt( 1.0 ),
        M_state( TS_UNITIALIZED ),
        M_n_restart( 0 ),
        M_restart( false ),
        M_restartPath( "" ),
        M_restartStepBeforeLastSave( 0 ),
        M_restartAtLastSave( false ),
        M_saveInFile( true ),
        M_saveFreq( 1 ),
        M_rankProcInNameOfFiles( false ),
        M_worldComm( Environment::worldComm() )
    {}

    TSBase( po::variables_map const& vm, std::string name, std::string const& prefix, WorldComm const& worldComm )
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
        M_state( TS_UNITIALIZED ),
        M_n_restart( 0 ),
        M_restart( vm[prefixvm( prefix, "bdf.restart" )].as<bool>() ),
        M_restartPath( vm[prefixvm( prefix, "bdf.restart.path" )].as<std::string>() ),
        M_restartStepBeforeLastSave( vm[prefixvm( prefix, "bdf.restart.step-before-last-save" )].as<int>() ),
        M_restartAtLastSave( vm[prefixvm( prefix, "bdf.restart.at-last-save" )].as<bool>() ),
        M_saveInFile( vm[prefixvm( prefix, "bdf.save" )].as<bool>() ),
        M_saveFreq( vm[prefixvm( prefix, "bdf.save.freq" )].as<int>() ),
        M_rankProcInNameOfFiles( vm[prefixvm( prefix, "bdf.rank-proc-in-files-name" )].as<bool>() ),
        M_worldComm( worldComm )
    {
    }
    TSBase( std::string name, WorldComm const& worldComm )
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
        M_state( TS_UNITIALIZED ),
        M_n_restart( 0 ),
        M_restart( false ),
        M_restartPath( "" ),
        M_restartStepBeforeLastSave( 0 ),
        M_restartAtLastSave( false ),
        M_saveInFile( true ),
        M_saveFreq( 1 ),
        M_worldComm( worldComm )
    {
    }

    TSBase( TSBase const& b )
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
        M_state( b.M_state ),
        M_n_restart( b.M_n_restart ),
        M_restart( b.M_restart ),
        M_restartPath( b.M_restartPath ),
        M_restartStepBeforeLastSave( b.M_restartStepBeforeLastSave ),
        M_restartAtLastSave( b.M_restartAtLastSave ),
        M_time_values_map( b.M_time_values_map ),
        M_saveInFile( b.M_saveInFile ),
        M_saveFreq( b.M_saveFreq ),
        M_rankProcInNameOfFiles( b.M_rankProcInNameOfFiles ),
        M_worldComm( b.M_worldComm )
    {}

    virtual ~TSBase() {}

    TSBase& operator=( TSBase const& b )
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
            M_state = b.M_state;
            M_restart = b.M_restart;
            M_restartPath = b.M_restartPath;
            M_restartStepBeforeLastSave = b.M_restartStepBeforeLastSave;
            M_restartAtLastSave = b.M_restartAtLastSave;

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

    //! reload metadata to get all time step already done
    void reloadMetaData()
    {
        CHECK(  M_restart ) <<" bdf is not in restart mode ";
        this->init();
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

    //! return the step index before the last used to restart
    int restartStepBeforeLastSave() const
    {
        return M_restartStepBeforeLastSave;
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
        CHECK( state() == TS_RUNNING ) << "invalid Time Stepping state, it should be " << TS_RUNNING << " (TS_RUNNING) and it is " << state();
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
        M_state = TS_RUNNING;
        M_timer.restart();
        M_iteration = 1;
        M_time = M_Ti+this->timeStep();
        // warning: this needs to be fixed wrt restart
        M_last_iteration_since_order_change = 1;
        //M_order_cur = 1;
        M_order_cur = M_order;
        return M_Ti;
    }

    //! restart the bdf
    double restart()
    {
        M_state = TS_RUNNING;
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
                M_state = TS_STOPPED;
                finished=true;
            }
        }

        else
        {
            if ( M_time < M_Tf )
            {
                M_state = TS_STOPPED;
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
        CHECK( state() == TS_RUNNING ) << "invalid Time Stepping state, it should be " << TS_RUNNING << " (TS_RUNNING) and it is " << state();
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

        DVLOG(2) << "[TSBase::shiftRight] solution name " << ostr.str() << "\n";

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

    }
    //! return the state of the BDF
    TSState state() const
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

    TSStragegy strategy() const
    {
        //return M_strategy;
        return TSStragegy::TS_STRATEGY_DT_CONSTANT;
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
    void setRestartPath( std::string const& s )
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

    //! state of the time stepping algorithm
    mutable TSState M_state;

    //! restart
    int M_n_restart;

    //! do a restart
    bool M_restart;

    //! restart path
    fs::path M_restartPath;

    //! do restart with the ieme step before the last save
    int M_restartStepBeforeLastSave;

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
                //TSBaseMetadata bdfloader( *this );
                //bdfloader.load();

                // modify Ti with last saved time
                if ( this->doRestartAtLastSave() )
                    M_Ti = M_time_values_map.back();
                // modify Ti with ieme step before the last
                if ( this->restartStepBeforeLastSave() > 0 )
                {
                    const int nbTimeStep = M_time_values_map.size();
                    CHECK( nbTimeStep-1-this->restartStepBeforeLastSave() >=0 ) << "error with restartStepBeforeLastSave " << this->restartStepBeforeLastSave()
                                                                               << "must be less to " << nbTimeStep << std::endl;
                    M_Ti =  M_time_values_map[ nbTimeStep-1-this->restartStepBeforeLastSave() ];
                }

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
class TSBaseMetadata
{
public:
    TSBaseMetadata( TSBase& bdf )
        :
        M_ts( bdf )
    {
    }

    void load()
    {
        fs::ifstream ifs;

        if ( M_ts.restartPath().empty() ) ifs.open( M_ts.path()/"metadata" );

        else ifs.open( M_ts.restartPath()/M_ts.path()/"metadata" );

        //fs::ifstream ifs( M_ts.path() / "metadata");

        boost::archive::text_iarchive ia( ifs );
        ia >> BOOST_SERIALIZATION_NVP( M_ts );
        DVLOG(2) << "[Bdf::init()] metadata loaded\n";
    }

    void save()
    {
        if ( !M_ts.saveInFile() ) return;

        // only master process write
        if ( M_ts.worldComm().isMasterRank() )
        {
            fs::ofstream ofs( M_ts.path() / "metadata" );

            boost::archive::text_oarchive oa( ofs );
            oa << BOOST_SERIALIZATION_NVP( ( TSBase const& )M_ts );
            DVLOG(2) << "[Bdf::init()] metadata saved\n";
        }
        // to be sure that all process can read the metadata file
        M_ts.worldComm().barrier();

    }

private:
    TSBase& M_ts;
};

} // namespace Feel
#endif
