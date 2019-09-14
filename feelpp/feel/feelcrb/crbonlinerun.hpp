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
//! @date 02 Sep 2019
//! @copyright 2019 Feel++ Consortium
//!
#ifndef FEELPP_CRBONLINERUN_HPP
#define FEELPP_CRBONLINERUN_HPP 1

#include <feel/feelcore/environment.hpp>
#include <feel/feelcrb/crbplugin_interface.hpp>


namespace Feel {

std::shared_ptr<Feel::CRBPluginAPI>
loadPlugin( std::string const& name, std::string const& id,
            std::string const& dirname = Environment::expand( soption(_name="plugin.dir") ),
            std::string const& libname = "",
            bool loadFiniteElementDatabase = false,
            int loadLast = 2 );


std::shared_ptr<Feel::CRBPluginAPI>
loadPlugin( std::string const& pluginname = Environment::expand( soption(_name="plugin.name") ),
            std::string const& dirname = Environment::expand( soption(_name="plugin.dir") ),
            std::string const& libname = "",
            bool loadFiniteElementDatabase = false,
            int loadLast = 2 );

std::string
loadModelName( std::string const& filename );

bool
runCrbOnline( std::vector<std::shared_ptr<Feel::CRBPluginAPI>> plugin );
    
void runCrbOnlineQuery( std::string const& dirname, std::string const& libname, bool loadFiniteElementDatabase, int loadLast );
void runCrbOnlineList();
void runCrbOnlineCompare( std::string const& dirname, std::string const& libname, bool loadFiniteElementDatabase, int loadLast );
}

#endif
