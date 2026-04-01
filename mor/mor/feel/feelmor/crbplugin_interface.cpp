//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 11 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!
#include <algorithm>
#include <boost/dll.hpp>
#include <feel/feelmor/crbplugin_interface.hpp>



namespace Feel {
namespace dll = boost::dll;

namespace
{
std::string
normalizePluginStem( std::string libname, std::string const& pluginname )
{
    if ( libname.empty() )
        libname = fmt::format( "feelpp_mor_plugin_{}", pluginname );

    auto const last_separator = libname.find_last_of( "/\\" );
    if ( last_separator != std::string::npos )
        libname = libname.substr( last_separator + 1 );

    auto const so_pos = libname.find( ".so" );
    if ( so_pos != std::string::npos )
        libname = libname.substr( 0, so_pos );

    if ( libname.rfind( "lib", 0 ) == 0 )
        libname = libname.substr( 3 );

    return libname;
}

std::vector<fs::path>
resolveVersionedPluginCandidates( fs::path const& dirname, std::string const& stem )
{
    std::vector<fs::path> candidates;
    if ( !fs::exists( dirname ) || !fs::is_directory( dirname ) )
        return candidates;

    auto const versioned_prefix = fmt::format( "lib{}.so.", stem );
    auto const bare_name = fmt::format( "lib{}", stem );
    auto const soname = fmt::format( "lib{}.so", stem );
    for ( auto const& entry : fs::directory_iterator( dirname ) )
    {
        auto const filename = entry.path().filename().string();
        if ( filename == bare_name || filename == soname || filename.rfind( versioned_prefix, 0 ) == 0 )
            candidates.push_back( entry.path() );
    }

    std::sort( candidates.begin(), candidates.end(),
               []( fs::path const& lhs, fs::path const& rhs )
               {
                   auto const& lhs_name = lhs.filename().string();
                   auto const& rhs_name = rhs.filename().string();
                   if ( lhs_name.size() != rhs_name.size() )
                       return lhs_name.size() < rhs_name.size();
                   return lhs_name < rhs_name;
               } );
    candidates.erase( std::unique( candidates.begin(), candidates.end() ), candidates.end() );
    return candidates;
}
}

namespace detail
{
class CRBPluginManagerImpl:
        public std::map<std::string, crbpluginapi_create_ft >,
        public boost::noncopyable
{
public:
    using value_type = crbpluginapi_create_ft;
    using key_type = std::string;

    using crbplugin_manager_type = std::map<key_type, value_type>;

};
struct CRBPluginManager : public  Feel::Singleton<CRBPluginManagerImpl> {
    ~CRBPluginManager() {}
};
}

std::shared_ptr<CRBPluginAPI>
factoryCRBPlugin( std::string const& pluginname, std::string const& pluginlibname, std::string const& dirname )
{
    auto p = Feel::detail::CRBPluginManager::instance().find( pluginname );
    if ( p != Feel::detail::CRBPluginManager::instance().end() )
    {
        return p->second();
    }
    else
    {
        auto const plugin_symbol = "create_crbplugin_" + pluginname;
        auto const stem = normalizePluginStem( pluginlibname, pluginname );
        auto const dirpath = fs::path( dirname );
        std::string last_error;

        LOG( INFO ) << fmt::format( "[feelpp.mor.factoryCRBPlugin] plugin name: {} plugin libname: {} dirname: {}", pluginname, stem, dirname );

        auto try_load = [&]( fs::path const& libpath, dll::load_mode::type mode ) -> std::shared_ptr<CRBPluginAPI>
        {
            auto const boost_libpath = boost::dll::fs::path( libpath.string() );
            Feel::detail::CRBPluginManager::instance().operator[]( pluginname ) =
                boost::dll::import_alias<crbpluginapi_create_t>( boost_libpath, plugin_symbol, mode );
            auto loaded = Feel::detail::CRBPluginManager::instance().find( pluginname );
            auto plugin = loaded->second();
            LOG( INFO ) << fmt::format( "[feelpp.mor.factoryCRBPlugin] loaded plugin: {}", libpath.string() );
            Logger::flush();
            return plugin;
        };

        try
        {
            auto decorated_path = ( dirpath / stem ).make_preferred();
            LOG( INFO ) << fmt::format( "[feelpp.mor.factoryCRBPlugin] loading plugin: {}...", decorated_path.string() );
            Logger::flush();
            return try_load( decorated_path, dll::load_mode::append_decorations );
        }
        catch ( std::exception const& err )
        {
            last_error = err.what();
            LOG( WARNING ) << fmt::format( "[feelpp.mor.factoryCRBPlugin] decorated load failed for {}: {}", stem, last_error );
        }

        for ( auto const& candidate : resolveVersionedPluginCandidates( dirpath, stem ) )
        {
            try
            {
                LOG( INFO ) << fmt::format( "[feelpp.mor.factoryCRBPlugin] trying versioned plugin candidate: {}", candidate.string() );
                Logger::flush();
                return try_load( candidate, dll::load_mode::default_mode );
            }
            catch ( std::exception const& err )
            {
                last_error = err.what();
            }
        }

        throw std::runtime_error( fmt::format( "[feelpp.mor.factoryCRBPlugin] unable to load plugin {} from {}: {}", stem, dirname, last_error ) );
    }
}
}
