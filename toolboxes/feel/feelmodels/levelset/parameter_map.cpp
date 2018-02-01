#include <feel/feelmodels/levelset/parameter_map.hpp>

namespace Feel {

int parameter_map::iget( std::string const& param, int defaultValue) const
{
    return this->get<int>( param, int(defaultValue) );
}

double parameter_map::dget( std::string const& param, double defaultValue ) const
{
    return this->get<double>( param, double(defaultValue) );
}

std::string parameter_map::sget( std::string const& param, std::string const& defaultValue ) const
{
    return this->get<std::string>( param, defaultValue );
}

std::pair<bool, int> 
parameter_map::iparam( std::string const& param, int defaultValue ) const
{
    return this->param<int>( param, int(defaultValue) );
}

std::pair<bool, double> 
parameter_map::dparam( std::string const& param, double defaultValue ) const
{
    return this->param<double>( param, double(defaultValue) );
}

std::pair<bool, std::string> 
parameter_map::sparam( std::string const& param, std::string const& defaultValue ) const
{
    return this->param<std::string>( param, defaultValue );
}

} // namespace Feel
