/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 22 Nov 2019

 Copyright (C) 2019 Feel++ Consortium

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#include <feel/feelcore/environment.hpp>
#include <map>
#include <fmt/core.h>
#include <fmt/ostream.h>

int main( int argc, char** argv )
{
    using namespace Feel;
    using Feel::cout;
    Environment env( _argc = argc, _argv = argv,
                     _about = about( _name = "info",
                                     _author = "Feel++ Consortium",
                                     _email = "feelpp-devel@feelpp.org" ) );
    std::cout << fmt::format("       Feel++ version: {}", Info::version() ) << std::endl;
    std::cout << fmt::format(" Feel++ version major: {}", Info::versionMajor() ) << std::endl;
    std::cout << fmt::format(" Feel++ version minor: {}", Info::versionMinor() ) << std::endl;
    std::cout << fmt::format(" Feel++ version micro: {}", Info::versionMicro() ) << std::endl;
    std::cout << fmt::format("Feel++ version string: {}", Info::versionString() ) << std::endl;
    std::cout << fmt::format("       Feel++ buildid: {}", Info::buildId() ) << std::endl;
    std::cout << fmt::format("        Feel++ prefix: {}", Info::prefix() ) << std::endl;
    std::cout << fmt::format("        Feel++ libdir: {}", Info::libdir() ) << std::endl;
    std::cout << fmt::format("     Feel++ plugindir: {}", Info::plugindir() ) << std::endl;
    std::cout << fmt::format("       Feel++ datadir: {}", Info::datadir() ) << std::endl;
    
    if ( Info::gitMetadata().populated())
    {
        std::cout << fmt::format( "******* Git ******* " ) << std::endl;
        std::cout << fmt::format("      Feel++ git commit author: {} ",Info::gitMetadata().authorName()) << std::endl;
        std::cout << fmt::format("Feel++ git commit author email: {} ",Info::gitMetadata().authorEmail()) << std::endl;
    }
    else
    {
        std::cout << "Feel++ git metadata not populated" << std::endl;
    }
    return 0;
}
