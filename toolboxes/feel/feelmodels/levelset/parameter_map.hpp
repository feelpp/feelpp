#ifndef _PARAMETERS_MAP_HPP
#define _PARAMETERS_MAP_HPP 1

#include <boost/any.hpp>
#include <map>

namespace Feel {

class parameter_map
    : public std::map<std::string, boost::any>
{
public:
    template<typename T>
    T get( std::string const& param, T const& defaultValue ) const;

    int iget( std::string const& param, int defaultValue = 0) const;
    double dget( std::string const& param, double defaultValue = 0.) const;
    std::string sget( std::string const& param, std::string const& defaultValue = "") const;

    template<typename T>
    std::pair<bool, T> param( std::string const& param, T const& defaultValue ) const;

    std::pair<bool, int> iparam( std::string const& param, int defaultValue = 0) const;
    std::pair<bool, double> dparam( std::string const& param, double defaultValue = 0.) const;
    std::pair<bool, std::string> sparam( std::string const& param, std::string const& defaultValue = "") const;
};

template<typename T>
T parameter_map::get( std::string const& param, T const& defaultValue ) const
{
    auto p = this->find( param );
    if( p != this->end() )
    {
        return boost::any_cast<T>(p->second);
    }

    // default
    return defaultValue;
}

template<typename T>
std::pair<bool, T>
parameter_map::param( std::string const& param, T const& defaultValue ) const
{
    auto p = this->find( param );
    if( p != this->end() )
    {
        return std::make_pair( true, boost::any_cast<T>(p->second) );
    }

    // default
    return std::make_pair( false, defaultValue );
}

} // namespace Feel

#endif
