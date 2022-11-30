#include <feel/feelcore/cmdlineoptions.hpp>
#include <feel/feelcore/environment.hpp>

namespace Feel {

CommandLineOptions::CommandLineOptions( po::variables_map const& vm ):
    M_vm( vm ) // TODO check if equal to Environment::vm()
{}

CommandLineOptions::CommandLineOptions( po::options_description const& _options, init_function_type func )
{
    M_vm.emplace();
    auto mycmdparser = Environment::commandLineParser();
    po::parsed_options parsed = mycmdparser.options( _options ).
        style(po::command_line_style::allow_long | po::command_line_style::long_allow_adjacent | po::command_line_style::long_allow_next).
        allow_unregistered().run();
    po::store(parsed,*M_vm);
    for ( auto & configFile : Environment::configFiles() )
    {
        std::istringstream & iss = std::get<1>( configFile );
        po::store(po::parse_config_file(iss, _options,true), *M_vm);
    }
    if ( func )
        std::invoke(func, _options, *M_vm );
    po::notify(*M_vm);
}

po::variables_map const& CommandLineOptions::vm() const
{
    if ( M_vm.has_value() )
        return *M_vm;
    else
        return Environment::vm();
}

} // ns Feel
