#ifndef _PSLOGGER_HPP_
#define _PSLOGGER_HPP_

#include <string>

/** Writes output of the system command ps to a logfile.
    @author Christoph Winkelmann, 2002
*/
class PsLogger
{

public:

    /** Constructor.
        @param fileName name of the logfile
        @param format the format of the output of ps. By default, memory
        consumption (rss) and cpu load (pcpu) are logged. See man ps for more
        information about format specifiers.
    */
    PsLogger( std::string fileName, std::string format="rss pcpu" );

    /** writes the log message and the output of ps to the logfile.
        @param logMessage the log message to write into the logfile
    */
    void log( std::string logMessage );

private:

    std::string M_fileName;
    std::string M_command;

};

#endif // _PSLOGGER_HPP_
