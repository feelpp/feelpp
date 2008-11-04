/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-02-19

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006 Université Joseph Fourier

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
/**
   \file info.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-02-19
 */
#include <sstream>
#include <lifeconfig.h>
#include <life/lifecore/info.hpp>

#define stringize2(x) #x
#define stringize(x) stringize2(x)

namespace Life
{
unsigned long long
Info::buildId()
{
    return LIFE_BUILDID;
}

unsigned long long
Info::revision()
{
    return LIFE_REVISION;
}

unsigned long long
Info::version()
{
    return LIFE_VERSION;
}

unsigned int
Info::versionMajor()
{
    return LIFE_VERSION_MAJOR;
}

unsigned int
Info::versionMinor()
{
    return LIFE_VERSION_MINOR;
}

unsigned int
Info::versionMicro()
{
    return LIFE_VERSION_MICRO;
}

char const*
Info::versionString()
{
    static bool _created = false;
    static std::ostringstream ostr;
    if ( !_created )
        {
            ostr << LIFE_VERSION_MAJOR << "."
                 << LIFE_VERSION_MINOR << "."
                 << LIFE_VERSION_MICRO << "-r"
                 << LIFE_REVISION << "-"
                 << LIFE_BUILDID;
            _created = true;
        }
    return ostr.str().c_str();
}

char const*
Info::prefix()
{
    return stringize(LIFE_PREFIX);
}

char const*
Info::datadir()
{
    return stringize(LIFE_DATADIR);
}

}

