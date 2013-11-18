/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-02-19

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2008,2009,2010 Université de Grenoble 1


  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-02-19
 */
#include <sstream>

#include <feel/feelconfig.h>
#include <feel/feelinfo.h>
#include <feel/feelcore/info.hpp>

#define stringize2(x) #x
#define stringize(x) stringize2(x)

namespace Feel
{
char const*
Info::buildId()
{
    return stringize(FEELPP_BUILDID);
}

char const*
Info::revision()
{
    return stringize(FEELPP_REVISION);
}

unsigned long long
Info::version()
{
    return FEELPP_VERSION;
}

unsigned int
Info::versionMajor()
{
    return FEELPP_VERSION_MAJOR;
}

unsigned int
Info::versionMinor()
{
    return FEELPP_VERSION_MINOR;
}

unsigned int
Info::versionMicro()
{
    return FEELPP_VERSION_MICRO;
}

char const*
Info::versionString()
{
    return stringize( FEELPP_VERSION_MAJOR ) "." stringize( FEELPP_VERSION_MINOR ) "." stringize( FEELPP_VERSION_MICRO )  stringize(FEELPP_VERSION_PRERELEASE) stringize(FEELPP_VERSION_METADATA) stringize(FEELPP_BUILDID);
}

char const*
Info::prefix()
{
    return stringize( FEELPP_PREFIX );
}

char const*
Info::datadir()
{
    return stringize( FEELPP_DATADIR );
}

}
