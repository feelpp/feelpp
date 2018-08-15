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
#include <boost/dll.hpp>
#include <feel/feelcrb/crbplugin_interface.hpp>



namespace Feel {
namespace dll = boost::dll;

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
struct CRBPluginManager : public  Feel::Singleton<CRBPluginManagerImpl> {};
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
#if defined( __APPLE__ )
        std::string libext = ".dylib";
#else
        std::string libext = ".so";
#endif
        std::string libname = pluginlibname;
        if ( libname.empty() )
            libname = ("libfeelpp_crb_" + pluginname + libext);
        fs::path pname = fs::path(dirname) / libname;
        //std::cout << "loading " << pname.string() << std::endl;

        Feel::detail::CRBPluginManager::instance().operator[]( pluginname ) = 
            boost::dll::import_alias<crbpluginapi_create_t>(pname,
                                                            "create_crbplugin_"+pluginname,
                                                            dll::load_mode::append_decorations );
        auto p = Feel::detail::CRBPluginManager::instance().find( pluginname );
        auto plugin = p->second();
        //std::cout << "Loaded the plugin " << plugin->name().c_str() << std::endl;
        return plugin;
    }
}
}
