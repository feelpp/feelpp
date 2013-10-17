#include <cstdlib>

#include <unistd.h>
#include <sstream>

#include <petscsys.h>

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/pslogger.hpp>

namespace Feel
{

PsLogger::PsLogger( std::string fileName, WorldComm const& worldComm, std::string format,  bool all_procs)
    :
    M_worldComm( worldComm ),
    M_fileName( fileName + (boost::format("-%1%_%2%") %this->worldComm().globalSize() %this->worldComm().globalRank()).str() )
{
    std::stringstream command;

    bool only_one_proc=(! all_procs) && ( Environment::worldComm().globalRank()==Environment::worldComm().masterRank() );
    bool write_a_file = all_procs || only_one_proc ;
    if( write_a_file )
    {
        command << "echo logging output of ps, format: > " << this->fileName() << std::ends;
        system( command.str().c_str() );
        command.str( "" );
        command << "echo " << format << " >> " << this->fileName() << std::ends;
        system( command.str().c_str() );
        command.str( "" );
#if defined( __APPLE__ )
        command << "ps -p " << getpid()
                << " -o \"" << format << "\" >> " << this->fileName() << std::ends;
#else
        command << "ps --no-header -p " << getpid()
                << " -o \"" << format << "\" >> " << this->fileName() << std::ends;
#endif
        M_command = command.str();
    }//write a file
}

void PsLogger::log( std::string logMessage , bool all_procs )
{
    bool only_one_proc=(! all_procs) && ( Environment::worldComm().globalRank()==Environment::worldComm().masterRank() );
    bool write_a_file = all_procs || only_one_proc ;

    if ( logMessage.length() > 0 && write_a_file )
    {
        std::stringstream command;
        command << "echo " << logMessage << " >> " << this->fileName() << std::ends;
        system( command.str().c_str() );
    }

    Environment::logMemoryUsage( logMessage );
    system( M_command.c_str() );

}

}
