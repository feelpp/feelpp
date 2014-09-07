/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2014-03-17

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
   \file tsbase.cpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2014-03-17
 */

#include <feel/feelcore/environment.hpp>
#include <feel/feelts/tsbase.hpp>

namespace Feel
{

TSBase::TSBase()
    :
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

TSBase::TSBase( po::variables_map const& vm, std::string name, std::string const& prefix, WorldComm const& worldComm )
    :
    M_name( name ),
    M_time( vm[prefixvm( prefix, "bdf.time-initial" )].as<double>() ),
    M_iteration( 0 ),
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
{}
TSBase::TSBase( std::string name, WorldComm const& worldComm )
    :
    M_name( name ),
    M_time( 0. ),
    M_iteration( 0 ),
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
{}

TSBase::TSBase( TSBase const& b )
    :
    M_name( b.M_name ),
    M_time( b.M_time ),
    M_iteration( b.M_iteration ),
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

TSBase&
TSBase::operator=( TSBase const& b )
{
    if ( this != &b )
    {
        M_name = b.M_name;
        M_time = b.M_time;
        M_iteration = b.M_iteration;
        M_Ti = b.M_Ti;
        M_Tf = b.M_Tf;
        M_dt = b.M_dt;
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


void
TSBase::init()
{
#if 0
    std::ostringstream ostr;
    ostr << "bdf_o_" << M_order << "_dt_" << M_dt;
    //ostr << "bdf/" << M_name << "/o_" << M_order << "/dt_" << M_dt;
    M_path_save = ostr.str();
#endif
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
                    //std::cout << "[Bdf::init()] file index: " << M_iteration << "\n";
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







void
TSBaseMetadata::load()
{
    fs::ifstream ifs;

    if ( M_ts.restartPath().empty() ) ifs.open( M_ts.path()/"metadata" );

    else ifs.open( M_ts.restartPath()/M_ts.path()/"metadata" );

    //fs::ifstream ifs( M_ts.path() / "metadata");

    boost::archive::text_iarchive ia( ifs );
    ia >> BOOST_SERIALIZATION_NVP( M_ts );
    DVLOG(2) << "[Bdf::init()] metadata loaded\n";
}

void
TSBaseMetadata::save()
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




} // namespace Feel
