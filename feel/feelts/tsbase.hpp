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
#include <feel/feelcore/worldcomm.hpp>


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

    TSBase();
    TSBase( po::variables_map const& vm, std::string name, std::string const& prefix, WorldComm const& worldComm );
    TSBase( std::string name, WorldComm const& worldComm );
    TSBase( TSBase const& b );

    virtual ~TSBase() {}

    TSBase& operator=( TSBase const& b );

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

        //DVLOG(2) << "[BDF::serialize] time orders size: " << M_time_orders.size() << "\n";
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
        // if initiliaze has been called M_iteration start to 1 else 0
        M_iteration = M_time_values_map.size();//1;
        M_time = M_Ti+this->timeStep();
        return M_Ti;
    }

    //! restart the bdf
    double restart()
    {
        M_state = TS_RUNNING;
        M_timer.restart();
        M_time = M_Ti+this->timeStep();
        ++M_iteration;
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
        return M_time;
    }

    virtual void shiftRight()
    {
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
    bool isSteady() const { return M_steady; }

    void setPathSave( std::string s )
    {
        M_path_save = s;
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
        M_steady = steady;
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

    virtual void print() const
    {
        LOG(INFO) << "============================================================\n";
        LOG(INFO) << "TS Information\n";
        LOG(INFO) << "   time step : " << this->timeStep() << "\n";
        LOG(INFO) << "time initial : " << this->timeInitial() << "\n";
        LOG(INFO) << "  time final : " << this->timeFinal() << "\n";
    }
protected:

    //! name of the file holding the bdf data
    std::string M_name;

    //! time
    mutable double M_time;

    //! iteration
    mutable int M_iteration;

    //! is steady
    bool M_steady;
    
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
    //std::vector<int> M_time_orders;

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
    void init();
};
class TSBaseMetadata
{
public:
    TSBaseMetadata( TSBase& bdf )
        :
        M_ts( bdf )
    {
    }

    void load();

    void save();

private:
    TSBase& M_ts;
};

} // namespace Feel
#endif
