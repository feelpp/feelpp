#ifndef FEELPP_CMDLINEOPTIONS_HPP
#define FEELPP_CMDLINEOPTIONS_HPP 1

#include <optional>
#include <functional>
#include <boost/program_options.hpp>

namespace Feel {

namespace po = boost::program_options;

struct CommandLineOptions
{
    using init_function_type = std::function<void(po::options_description const&,po::variables_map &)>;
    CommandLineOptions() = default;
    explicit CommandLineOptions( po::options_description const& _options, init_function_type func={} );
    explicit CommandLineOptions( po::variables_map const& vm );
    CommandLineOptions( CommandLineOptions const& ) = default;
    CommandLineOptions( CommandLineOptions && ) = default;

    po::variables_map const& vm() const;

private :
    std::optional<po::variables_map> M_vm;
};

} // ns Feel

#endif // FEELPP_CMDLINEOPTIONS_HPP
