/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2014-03-17

  Copyright (C) 2013-2016 Feel++ Consortium

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
    M_reverse( false ),
    M_reverseLoad( false ),
    M_n_restart( 0 ),
    M_restart( false ),
    M_restartPath( "" ),
    M_restartStepBeforeLastSave( 0 ),
    M_restartAtLastSave( false ),
    M_saveInFile( true ),
    M_saveFreq( 1 ),
    M_rankProcInNameOfFiles( false ),
    M_fileFormat( "binary" ),
    M_worldComm( Environment::worldComm() ),
    M_prefix(),
    M_displayStats( false )
{}

TSBase::TSBase( std::string name, std::string const& prefix, WorldComm const& worldComm, po::variables_map const& vm )
    :
    M_name( name ),
    M_time( doption(_prefix=prefix,_name="ts.time-initial",_vm=vm) ),
    M_iteration( 0 ),
    M_Ti( doption(_prefix=prefix,_name="ts.time-initial",_vm=vm) ),
    M_Tf( doption(_prefix=prefix,_name="ts.time-final",_vm=vm) ),
    M_dt( doption(_prefix=prefix,_name="ts.time-step",_vm=vm) ),
    M_state( TS_UNITIALIZED ),
    M_reverse( false ),
    M_reverseLoad( false ),
    M_n_restart( 0 ),
    M_restart( boption(_prefix=prefix,_name="ts.restart",_vm=vm) ),
    M_restartPath( soption(_prefix=prefix,_name="ts.restart.path",_vm=vm) ),
    M_restartStepBeforeLastSave( ioption(_prefix=prefix,_name="ts.restart.step-before-last-save",_vm=vm) ),
    M_restartAtLastSave( boption(_prefix=prefix,_name="ts.restart.at-last-save",_vm=vm) ),
    M_saveInFile( boption(_prefix=prefix,_name="ts.save",_vm=vm) ),
    M_saveFreq( ioption(_prefix=prefix,_name="ts.save.freq",_vm=vm) ),
    M_rankProcInNameOfFiles( boption(_prefix=prefix,_name="ts.rank-proc-in-files-name",_vm=vm) ),
    M_fileFormat( soption(_prefix=prefix,_name="ts.file-format",_vm=vm) ),
    M_worldComm( worldComm ),
    M_prefix( prefix ),
    M_displayStats( boption(_prefix=prefix,_name="ts.display-stats",_vm=vm) )
{}

TSBase::TSBase( std::string name, std::string const& prefix, WorldComm const& worldComm, po::variables_map const& vm,
                double ti, double tf, double dt, bool steady, bool reverse, bool restart, std::string const& restart_path, bool restart_at_last_save,
                bool save, int freq, bool rank_proc_in_files_name, std::string const& format )
    :
    M_name( name ),
    M_time( ti ),
    M_iteration( 0 ),
    M_steady( false ),
    M_Ti( ti ),
    M_Tf( tf ),
    M_dt( dt ),
    M_state( TS_UNITIALIZED ),
    M_reverse( false ),
    M_reverseLoad( false ),
    M_n_restart( 0 ),
    M_restart( restart ),
    M_restartPath( restart_path ),
    M_restartStepBeforeLastSave( ioption(_prefix=prefix,_name="ts.restart.step-before-last-save",_vm=vm) ),
    M_restartAtLastSave( restart_at_last_save ),
    M_saveInFile( save ),
    M_saveFreq( freq ),
    M_rankProcInNameOfFiles( rank_proc_in_files_name ),
    M_fileFormat( format ),
    M_worldComm( worldComm ),
    M_prefix( prefix ),
    M_displayStats( boption(_prefix=prefix,_name="ts.display-stats",_vm=vm) )
{
    this->setSteady( steady );
    this->setReverse( reverse );
}

TSBase::TSBase( std::string name, WorldComm const& worldComm )
    :
    M_name( name ),
    M_time( 0. ),
    M_iteration( 0 ),
    M_Ti( 0. ),
    M_Tf( 1.0 ),
    M_dt( 1.0 ),
    M_state( TS_UNITIALIZED ),
    M_reverse( false ),
    M_reverseLoad( false ),
    M_n_restart( 0 ),
    M_restart( false ),
    M_restartPath( "" ),
    M_restartStepBeforeLastSave( 0 ),
    M_restartAtLastSave( false ),
    M_saveInFile( true ),
    M_saveFreq( 1 ),
    M_fileFormat( "binary" ),
    M_worldComm( worldComm ),
    M_prefix(),
    M_displayStats( false)
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
    M_reverse( b.M_reverse ),
    M_reverseLoad( b.M_reverseLoad ),
    M_n_restart( b.M_n_restart ),
    M_restart( b.M_restart ),
    M_restartPath( b.M_restartPath ),
    M_restartStepBeforeLastSave( b.M_restartStepBeforeLastSave ),
    M_restartAtLastSave( b.M_restartAtLastSave ),
    M_time_values_map( b.M_time_values_map ),
    M_saveInFile( b.M_saveInFile ),
    M_saveFreq( b.M_saveFreq ),
    M_rankProcInNameOfFiles( b.M_rankProcInNameOfFiles ),
    M_fileFormat( b.M_fileFormat ),
    M_worldComm( b.M_worldComm ),
    M_prefix( b.M_prefix ),
    M_displayStats( b.M_displayStats )
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
        M_reverse = b.M_reverse;
        M_reverseLoad = b.M_reverseLoad;
        M_restart = b.M_restart;
        M_restartPath = b.M_restartPath;
        M_restartStepBeforeLastSave = b.M_restartStepBeforeLastSave;
        M_restartAtLastSave = b.M_restartAtLastSave;

        M_saveInFile = b.M_saveInFile;
        M_rankProcInNameOfFiles = b.M_rankProcInNameOfFiles;
        M_fileFormat = b.M_fileFormat;

        M_time_values_map = b.M_time_values_map;
        M_worldComm = b.M_worldComm;

        M_displayStats = b.M_displayStats;
    }

    return *this;
}

void
TSBase::init()
{
    CHECK( M_fileFormat == "binary" || M_fileFormat == "hdf5" ) << "invalid file format " << M_fileFormat;
#if 0
    std::ostringstream ostr;
    ostr << "bdf_o_" << M_order << "_dt_" << M_dt;
    //ostr << "bdf/" << M_name << "/o_" << M_order << "/dt_" << M_dt;
    M_path_save = ostr.str();
#endif

    if ( this->saveInFile() )
    {
        // If directory does not exist, create it.
        if ( this->worldComm().isMasterRank() && !fs::exists( M_path_save ) )
            fs::create_directories( M_path_save );
        // Be sure that all process can find the path after.
        this->worldComm().barrier();
    }

    if ( M_restart )
    {
        fs::path thepath;
        if ( not this->restartPath().empty() )
        {   // Default metadata path.
            thepath=this->restartPath()/this->path()/"metadata";
        }
        else
        {   // Restart metadata path.
            thepath=this->path()/"metadata";
        }

        // Read the saved bdf metadata.
        if ( fs::exists( thepath ) )
        {
            DVLOG(2) << "[TSBase::init()] loading metadata from " << M_path_save.string() << "\n";

            std::ifstream ifs( thepath );
            boost::archive::text_iarchive ia( ifs );
            ia >> BOOST_SERIALIZATION_NVP( *this );

            DVLOG(2) << "[TSBase::init()] metadata loaded\n";

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

            int iteration = 0;
            // Determine if Ti is found in map.
            bool found = false;

            // Determine current iteration.
            for( auto time : M_time_values_map )
            {
                if ( math::abs( time-M_Ti ) < 1e-10 )
                {
                    found = true;
                    break;
                }
                ++iteration;
            }

            if ( found )
            {
                M_iteration = iteration;
                if( this->isReverse() )
                {
                    M_iteration = this->iterationNumber() - iteration;
                }
                if (this->saveFreq()==1)
                {
                    M_time_values_map.resize( M_iteration+1 );
                }
                else // saveFreq != 1.
                {
                    int nItBack = M_iteration % this->saveFreq();
                    M_iteration-=nItBack;

                    M_time_values_map.resize( M_iteration+1 );
                    M_Ti = M_time_values_map.back();
                }
            }
            else // M_Ti not found.
            {
                DVLOG(2) << "[TSBase::init()] initial time " << M_Ti << " not found\n";
                M_time_values_map.clear();
            }
            DVLOG(2) << "[TSBase::init()] initial time is Ti=" << M_Ti << "\n";
            DVLOG(2) << "[TSBase::init()] file index: " << M_iteration << "\n";
        }
        else // Metadata path does not exist.
        {
            M_time_values_map.clear();
        }
    } // restart
} // init

void
TSBaseMetadata::load()
{
    std::ifstream ifs;

    if ( M_ts.restartPath().empty() )
    {   // Default metadata path.
        ifs.open( M_ts.path()/"metadata" );
    }
    else
    {   // Restart metadata path.
        ifs.open( M_ts.restartPath()/M_ts.path()/"metadata" );
    }

    boost::archive::text_iarchive ia( ifs );
    ia >> BOOST_SERIALIZATION_NVP( M_ts );
    DVLOG(2) << "[TSBaseMetadata::init()] metadata loaded\n";
}

void
TSBaseMetadata::save()
{
    if ( !M_ts.saveInFile() ) return;

    // only master process write
    if ( M_ts.worldComm().isMasterRank() )
    {
        std::ofstream ofs( M_ts.path() / "metadata" );

        boost::archive::text_oarchive oa( ofs );
        oa << BOOST_SERIALIZATION_NVP( ( TSBase const& )M_ts );
        DVLOG(2) << "[TSBaseMetadata::init()] metadata saved\n";
    }
    // to be sure that all process can read the metadata file
    M_ts.worldComm().barrier();
}


} // namespace Feel
