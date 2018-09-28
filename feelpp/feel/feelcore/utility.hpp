/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 27 sept. 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
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
#ifndef FEELPP_UTILITY_HPP
#define FEELPP_UTILITY_HPP 1

#include <feel/feelcore/feelmacros.hpp>
#include <string>

// the following snippet of code detects the current OS and
// defines the appropriate macro that is used to wrap some
// platform specific things
#if defined(_WIN32) || defined(_WIN64)
#define FEELPP_OS_WINDOWS
#elif defined(__APPLE__)
#define FEELPP_OS_MACOS
#elif defined(__unix__) || defined(__unix)
#define FEELPP_OS_LINUX
#else
#error unsupported platform
#endif

namespace Feel {

/**
 * helper function that read a file \p filename and store its contents into a \c
 * std::string
 */
FEELPP_EXPORT std::string readFromFile( std::string const& filename );


//! Get a unique char (no return cariage requrired).
//! This function is a custom getch.
FEELPP_EXPORT int getOneChar();

//! Ask the user to fill a password.
//! The password is hidden with asterisks.
//! \param message a message to display to the user
//! \return The (clear) password in a string.
//
//! \note: This function is unix specific. 
FEELPP_EXPORT std::string askPassword( std::string const& message );

}
#endif
