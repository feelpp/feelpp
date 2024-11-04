#ifndef FEELPP_MPI_LOG_SINK_HPP
#define FEELPP_MPI_LOG_SINK_HPP

#include <glog/logging.h>
#include <fmt/core.h>
#include <fstream>
#include <string>
#include <thread>
#include <chrono>

namespace Feel 
{

/**
 * @class MpiLogSink
 * @brief Custom logging sink for MPI applications to direct logs to separate files or console.
 */
class MpiLogSink : public google::LogSink 
{
public:
    /**
     * @brief Constructor for MpiLogSink
     * @param rank The MPI rank of the process
     * @param log_option Option to control logging behavior ("master" for rank 0, "all" for all ranks)
     * @param output_option Option to control console output ("stdout", "stderr", or "none")
     * @param log_dir Directory to store log files (default is current working directory)
     * @param log_memory Flag to enable logging memory usage (default is false)
     */
    MpiLogSink(int rank, const std::string& log_option, const std::string& output_option, const std::string& log_dir = "", bool log_memory = false);

    /**
     * @brief Destructor for MpiLogSink
     */
    ~MpiLogSink() override;

    /**
     * @brief Sends the log message to the appropriate output (file or console)
     * @param severity Log severity level
     * @param full_filename Full path of the source file generating the log message
     * @param base_filename Base name of the source file generating the log message
     * @param line Line number where the log message is generated
     * @param tm_time Pointer to time structure for the log message
     * @param message Log message content
     * @param message_len Length of the log message content
     */
    void send(google::LogSeverity severity, const char* full_filename,
              const char* base_filename, int line,
              const struct ::tm* tm_time, const char* message, size_t message_len) override;

private:
    enum class LogOption { Master, All, None };
    enum class OutputOption { None, Stdout, Stderr };

    int rank_;                      ///< MPI rank
    LogOption log_option_;        ///< Logging option: "master" or "all"
    OutputOption output_option_;     ///< Console output option: "stdout", "stderr", or "none"
    std::ofstream log_file_;        ///< File stream for log output
    bool log_memory_;               ///< Flag to enable logging memory usage

    // Helper functions for converting strings to enum values
    static LogOption logOptionFromString(const std::string& option);
    static OutputOption outputOptionFromString(const std::string& option);
};

} // namespace Feel

#endif // FEELPP_MPI_LOG_SINK_HPP