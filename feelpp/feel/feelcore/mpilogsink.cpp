#include <feel/feelcore/mpilogsink.hpp>
#include <glog/logging.h>
#include <mpi.h>
#include <fmt/core.h>
#include <fmt/chrono.h>
#include <fmt/color.h>
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <sstream>
#include <thread>

namespace Feel 
{


MpiLogSink::MpiLogSink(int rank, const std::string& log_option, const std::string& output_option, const std::string& log_dir, bool log_memory)
    : rank_(rank),
      log_option_(logOptionFromString(log_option)),
      output_option_(outputOptionFromString(output_option)),
      log_memory_(log_memory)
{
    if ((log_option_ == LogOption::Master && rank_ == 0) || log_option_ == LogOption::All) 
    {
        std::string log_filename = fmt::format("{}_{}.log", log_dir, rank_);
        log_file_.open(log_filename, std::ios::out);

        if (!log_file_.is_open()) 
        {
            throw std::runtime_error("Failed to open log file: " + log_filename);
        }
    }
}

MpiLogSink::~MpiLogSink() 
{
    this->close();
}

void
MpiLogSink::close()
{
    if (log_file_.is_open()) 
    {
        log_file_.close();
    }
}
double getMemoryUsageInGB() 
{
    std::ifstream file("/proc/self/status");
    std::string line;
    double memory_kb = 0.0;

    while (std::getline(file, line)) {
        if (line.find("VmRSS:") == 0) {  // Resident Set Size (physical memory)
            std::istringstream iss(line);
            std::string key, unit;
            iss >> key >> memory_kb >> unit;
            break;
        }
    }
    return memory_kb / (1024 * 1024);  // Convert from KB to GB
}

std::string formatForConsole(google::LogSeverity severity, const std::string& message) 
{
    const char* color_code;
    switch (severity) {
        case google::INFO: color_code = "\033[32m"; break;    // Green for INFO
        case google::WARNING: color_code = "\033[33m"; break; // Yellow for WARNING
        case google::ERROR: color_code = "\033[31m"; break;   // Red for ERROR
        case google::FATAL: color_code = "\033[35m"; break;   // Magenta for FATAL
        default: color_code = "\033[0m"; break;
    }
    return fmt::format("{}{}{}", color_code, message, "\033[0m");
}

void MpiLogSink::send(google::LogSeverity severity, const char* full_filename,
                      const char* base_filename, int line,
                      const struct ::tm* tm_time, const char* message, size_t message_len) 
{


    std::ostringstream oss;
    oss << std::this_thread::get_id();
    std::string thread_id_str = oss.str();

    auto now = std::chrono::system_clock::now();
    auto now_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()).count() % 1000;

    //std::string severity_name = fmt::format("{:<7}", google::GetLogSeverityName(severity));

    std::string memory_usage = log_memory_ ? fmt::format(" [Mem: {:.2f} GB] ", getMemoryUsageInGB()) : " ";

    const bool do_log =( log_option_ == LogOption::All ) || ( log_option_ == LogOption::Master && rank_ == 0 );
    if (do_log && log_file_.is_open() )
    {
        std::string log_message = fmt::format("[{}]: [{}] [{:%Y-%m-%d %H:%M:%S}.{:03}]{}[{}:{}]: {}\n",
                                              rank_, //thread_id_str, 
                                              google::GetLogSeverityName(severity)[0], 
                                              *tm_time, now_ms, 
                                              memory_usage,
                                              base_filename, line,
                                              std::string(message, message_len));
        log_file_ << log_message;
        log_file_.flush();
    }

    if (output_option_ != OutputOption::None || severity == google::ERROR || severity == google::FATAL)  
    {
        // Determine color based on severity
        fmt::text_style severity_style;
        switch (severity) {
            case google::INFO:
                severity_style = fmt::fg(fmt::color::green) | fmt::emphasis::bold;
                break;
            case google::WARNING:
                severity_style = fmt::fg(fmt::color::yellow) | fmt::emphasis::bold;
                break;
            case google::ERROR:
                severity_style = fmt::fg(fmt::color::red) | fmt::emphasis::bold;
                break;
            case google::FATAL:
                severity_style = fmt::fg(fmt::color::magenta) | fmt::emphasis::bold;
                break;
            default:
                severity_style = fmt::fg(fmt::color::white);
                break;
        }

        // Create log message with severity-specific styling
        std::string log_message = fmt::format("{}: {} [{:%Y-%m-%d %H:%M:%S}.{:03}]{}[{}]: {}\n",
                                      fmt::styled(fmt::format("[{}]", rank_), fmt::fg(fmt::color::blue) | fmt::emphasis::bold),
                                      fmt::styled(fmt::format("[{}]", google::GetLogSeverityName(severity)[0]), severity_style),
                                      *tm_time, now_ms,
                                      memory_usage,
                                      fmt::styled(fmt::format("{}:{}", base_filename, line), fmt::emphasis::underline),
                                      std::string(message, message_len));

        std::string console_message = formatForConsole(severity, log_message);
        if (output_option_ == OutputOption::Stdout && do_log )
        {
            std::cout << console_message;
        } 
        else if ( (output_option_ == OutputOption::Stderr) && do_log && ( severity != google::ERROR ) && ( severity != google::FATAL ) )
        {
            std::cerr << console_message;
        }
        if ( severity == google::ERROR || severity == google::FATAL) 
        {
            std::cerr << console_message;
        }
        if (severity == google::FATAL) 
        {
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
}



MpiLogSink::LogOption MpiLogSink::logOptionFromString(const std::string& option) 
{
    if (option == "master") return LogOption::Master;
    if (option == "all") return LogOption::All;
    return LogOption::None;
}

MpiLogSink::OutputOption MpiLogSink::outputOptionFromString(const std::string& option) 
{
    if (option == "stdout") return OutputOption::Stdout;
    if (option == "stderr") return OutputOption::Stderr;
    return OutputOption::None;
}

} // namespace Feel