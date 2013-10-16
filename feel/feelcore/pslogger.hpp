#ifndef _PSLOGGER_HPP_
#define _PSLOGGER_HPP_

#include <string>
#include <feel/feelcore/environment.hpp>

namespace Feel
{

/** Writes output of the system command ps to a logfile.
    @author Christoph Winkelmann, 2002
*/
class PsLogger
{

public:

    /** Constructor.
        @param fileName name of the logfile
        @param format the format of the output of ps. By default, memory
        consumption (rss), the percentage of real memory used by this process (pmem)
        and cpu load (pcpu) are logged. See man ps for more
        information about format specifiers.
        @param all_procs : write on all processors if true
    */
    PsLogger( std::string fileName, WorldComm const& worldComm=Environment::worldComm(), std::string format="rss pmem pcpu" , bool all_procs=true);

    /** writes the log message and the output of ps to the logfile.
        @param logMessage the log message to write into the logfile
        @param all_procs : write on all processors if true
    */
    void log( std::string logMessage , bool all_procs=true);

    std::string fileName() const { return M_fileName; }

    WorldComm const& worldComm() const { return M_worldComm; }

private:

    WorldComm M_worldComm;
    std::string M_fileName;
    std::string M_command;

};

} // namespace Feel

#endif // _PSLOGGER_HPP_
