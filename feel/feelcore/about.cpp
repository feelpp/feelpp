/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
    kate: space-indent on; indent-width 4; mixedindent off; indent-mode cstyle; encoding utf-8;
   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2005-03-17

   Copyright (C) 2009 Université de Grenoble 1
   Copyright (C) 2005,2006 EPFL

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
   \file about.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-03-17
 */
#include <cstring>
#include <iostream>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/about.hpp>

namespace Feel
{
std::string
AboutPerson::name() const
{
    return _M_Name;
}

std::string
AboutPerson::task() const
{
    return _M_Task;
}

std::string
AboutPerson::emailAddress() const
{
    return _M_EmailAddress;
}


std::string
AboutPerson::webAddress() const
{
    return _M_WebAddress;
}


class AboutDataPrivate
{
public:
    AboutDataPrivate()
    {}
    AboutDataPrivate( AboutDataPrivate const& /* d */ )
    {}
    ~AboutDataPrivate()
    {
    }
};

AboutData::AboutData( const char*  appName,
                      const char*  programName,
                      const char*  version,
                      const char*  shortDescription,
                      int licenseType,
                      const char*  copyrightStatement,
                      const char*  text,
                      const char*  homePageAddress,
                      const char*  bugsEmailAddress
                    ) :
    _M_ProgramName( programName ),
    _M_Version( version ),
    _M_ShortDescription( shortDescription ),
    _M_LicenseKey( licenseType ),
    _M_CopyrightStatement( copyrightStatement ),
    _M_OtherText( text ),
    _M_HomepageAddress( homePageAddress ),
    _M_BugEmailAddress( bugsEmailAddress ),
    _M_LicenseText (),
    d( new AboutDataPrivate )
{
    if ( appName )
    {
        const char *p = ::strrchr( appName, '/' );

        if ( p )
            _M_AppName = p+1;

        else
            _M_AppName = appName;
    }

    else
        _M_AppName = "";
}
AboutData::AboutData( AboutData const& ad )
    :
    _M_AppName( ad._M_AppName ),
    _M_ProgramName( ad._M_ProgramName ),
    _M_Version( ad._M_Version ),
    _M_ShortDescription( ad._M_ShortDescription ),
    _M_LicenseKey( ad._M_LicenseKey ),
    _M_CopyrightStatement( ad._M_CopyrightStatement ),
    _M_OtherText( ad._M_OtherText ),
    _M_HomepageAddress( ad._M_HomepageAddress ),
    _M_BugEmailAddress( ad._M_BugEmailAddress ),
    _M_AuthorList( ad._M_AuthorList ),
    _M_CreditList( ad._M_CreditList ),
    _M_LicenseText ( ad._M_LicenseText ),
    d( new AboutDataPrivate( *ad.d ) )
{
}
AboutData::~AboutData()
{
    delete d;
}

void
AboutData::addAuthor( std::string const & name, std::string const & task,
                      std::string const & emailAddress, std::string const & webAddress )
{
    _M_AuthorList.push_back( AboutPerson( name,task,emailAddress,webAddress ) );

}

void
AboutData::addCredit( std::string const & name, std::string const & task,
                      std::string const & emailAddress, std::string const & webAddress )
{
    _M_CreditList.push_back( AboutPerson( name,task,emailAddress,webAddress ) );
}

void
AboutData::setLicenseText( std::string const & licenseText )
{
    _M_LicenseText = licenseText;
    _M_LicenseKey = License_Custom;
}

void
AboutData::setAppName( std::string const & appName )
{
    _M_AppName = appName;
}

void
AboutData::setProgramName( const char* programName )
{
    _M_ProgramName = programName;
}

void
AboutData::setVersion( const char* version )
{
    _M_Version = version;
}

void
AboutData::setShortDescription( std::string const & shortDescription )
{
    _M_ShortDescription = shortDescription;
}

void
AboutData::setLicense( LicenseKey licenseKey )
{
    _M_LicenseKey = licenseKey;
}

void
AboutData::setCopyrightStatement( std::string const & copyrightStatement )
{
    _M_CopyrightStatement = copyrightStatement;
}

void
AboutData::setOtherText( std::string const & otherText )
{
    _M_OtherText = otherText;
}

void
AboutData::setHomepage( std::string const & homepage )
{
    _M_HomepageAddress = homepage;
}

void
AboutData::setBugAddress( std::string const & bugAddress )
{
    _M_BugEmailAddress = bugAddress;
}

void
AboutData::setProductName( std::string const & productName )
{
    _M_ProductName = productName;
}

std::string
AboutData::appName() const
{
    return _M_AppName;
}

std::string
AboutData::productName() const
{
    if ( !_M_ProductName.empty() )
        return _M_ProductName;

    else
        return appName();
}

std::string
AboutData::programName() const
{
    return _M_ProgramName;
}

std::string
AboutData::version() const
{
    return _M_Version;
}

std::string
AboutData::shortDescription() const
{
    return _M_ShortDescription;
}

std::string
AboutData::homepage() const
{
    return _M_HomepageAddress;
}

std::string
AboutData::bugAddress() const
{
    return _M_BugEmailAddress;
}

const std::vector<AboutPerson>&
AboutData::authors() const
{
    return _M_AuthorList;
}

const std::vector<AboutPerson>&
AboutData::credits() const
{
    return _M_CreditList;
}
std::string
AboutData::otherText() const
{
    return _M_OtherText;
}


std::string
AboutData::license() const
{
    std::string result;

    if ( !copyrightStatement().empty() )
        result = copyrightStatement() + "\n\n";

    std::string l;
    std::string f;

    switch ( _M_LicenseKey )
    {
    case License_GPL_V2:
        l = "GPL v2";
        break;

    case License_LGPL_V2:
        l = "LGPL v2";
        break;

    case License_BSD:
        l = "BSD License";
        break;

    case License_Artistic:
        l = "Artistic License";
        break;

    case License_QPL_V1_0:
        l = "QPL v1.0";
        break;

    case License_Custom:
        if ( !_M_LicenseText.empty() )
            return( _M_LicenseText );

        // fall through
    default:
        result += "No licensing terms for this program have been specified.\n"
                  "Please check the documentation or the source for any\n"
                  "licensing terms.\n";
        return result;
    }

    if ( !l.empty() )
        result += "This program is distributed under the terms of the " + l;

    return result;
}

std::string
AboutData::copyrightStatement() const
{
    return _M_CopyrightStatement;
}

std::ostream&
operator<<( std::ostream& os, AboutData const& /* about */ )
{
    return os;
}
}
