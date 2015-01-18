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

#include <boost/parameter.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/parameter.hpp>
#include <feel/feelcore/about.hpp>

namespace Feel
{
std::string
AboutPerson::name() const
{
    return M_Name;
}

std::string
AboutPerson::task() const
{
    return M_Task;
}

std::string
AboutPerson::emailAddress() const
{
    return M_EmailAddress;
}


std::string
AboutPerson::webAddress() const
{
    return M_WebAddress;
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

AboutData::AboutData( std::string const & appName,
                      std::string const & programName,
                      std::string const & version,
                      std::string const & shortDescription,
                      int licenseType,
                      std::string const & copyrightStatement,
                      std::string const & text,
                      std::string const & homePageAddress,
                      std::string const & bugsEmailAddress
                    ) :
    M_ProgramName( programName ),
    M_Version( version ),
    M_ShortDescription( shortDescription ),
    M_LicenseKey( licenseType ),
    M_CopyrightStatement( copyrightStatement ),
    M_OtherText( text ),
    M_HomepageAddress( homePageAddress ),
    M_BugEmailAddress( bugsEmailAddress ),
    M_LicenseText ()//,
    //d( new AboutDataPrivate )
{
    if ( appName.size() > 0 )
    {
        size_t found = appName.find_last_of("/\\");

        if ( found != std::string::npos )
            M_AppName = appName.substr(found + 1);

        else
            M_AppName = appName;
    }

    else
        M_AppName = "";
}
AboutData::AboutData( AboutData const& ad )
    :
    M_AppName( ad.M_AppName ),
    M_ProgramName( ad.M_ProgramName ),
    M_Version( ad.M_Version ),
    M_ShortDescription( ad.M_ShortDescription ),
    M_LicenseKey( ad.M_LicenseKey ),
    M_CopyrightStatement( ad.M_CopyrightStatement ),
    M_OtherText( ad.M_OtherText ),
    M_HomepageAddress( ad.M_HomepageAddress ),
    M_BugEmailAddress( ad.M_BugEmailAddress ),
    M_AuthorList( ad.M_AuthorList ),
    M_CreditList( ad.M_CreditList ),
    M_LicenseText ( ad.M_LicenseText )//,
    //d( new AboutDataPrivate( *ad.d ) )
{
}
AboutData::~AboutData()
{
#if 0
    if ( d )
    {
        delete d;
        d = 0;
    }
#endif
}

void
AboutData::addAuthor( std::string const & name, std::string const & task,
                      std::string const & emailAddress, std::string const & webAddress )
{
    M_AuthorList.push_back( AboutPerson( name,task,emailAddress,webAddress ) );

}

void
AboutData::addCredit( std::string const & name, std::string const & task,
                      std::string const & emailAddress, std::string const & webAddress )
{
    M_CreditList.push_back( AboutPerson( name,task,emailAddress,webAddress ) );
}

void
AboutData::setLicenseText( std::string const & licenseText )
{
    M_LicenseText = licenseText;
    M_LicenseKey = License_Custom;
}

void
AboutData::setAppName( std::string const & appName )
{
    M_AppName = appName;
}

void
AboutData::setProgramName( std::string const & programName )
{
    M_ProgramName = programName;
}

void
AboutData::setVersion( std::string const & version )
{
    M_Version = version;
}

void
AboutData::setShortDescription( std::string const & shortDescription )
{
    M_ShortDescription = shortDescription;
}

void
AboutData::setLicense( LicenseKey licenseKey )
{
    M_LicenseKey = licenseKey;
}

void
AboutData::setCopyrightStatement( std::string const & copyrightStatement )
{
    M_CopyrightStatement = copyrightStatement;
}

void
AboutData::setOtherText( std::string const & otherText )
{
    M_OtherText = otherText;
}

void
AboutData::setHomepage( std::string const & homepage )
{
    M_HomepageAddress = homepage;
}

void
AboutData::setBugAddress( std::string const & bugAddress )
{
    M_BugEmailAddress = bugAddress;
}

void
AboutData::setProductName( std::string const & productName )
{
    M_ProductName = productName;
}

std::string
AboutData::appName() const
{
    return M_AppName;
}

std::string
AboutData::productName() const
{
    if ( !M_ProductName.empty() )
        return M_ProductName;

    else
        return appName();
}

std::string
AboutData::programName() const
{
    return M_ProgramName;
}

std::string
AboutData::version() const
{
    return M_Version;
}

std::string
AboutData::shortDescription() const
{
    return M_ShortDescription;
}

std::string
AboutData::homepage() const
{
    return M_HomepageAddress;
}

std::string
AboutData::bugAddress() const
{
    return M_BugEmailAddress;
}

const std::vector<AboutPerson>&
AboutData::authors() const
{
    return M_AuthorList;
}

const std::vector<AboutPerson>&
AboutData::credits() const
{
    return M_CreditList;
}
std::string
AboutData::otherText() const
{
    return M_OtherText;
}


std::string
AboutData::license() const
{
    std::string result;

    if ( !copyrightStatement().empty() )
        result = copyrightStatement() + "\n\n";

    std::string l;
    std::string f;

    switch ( M_LicenseKey )
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
        if ( !M_LicenseText.empty() )
            return( M_LicenseText );

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
    return M_CopyrightStatement;
}

std::ostream&
operator<<( std::ostream& os, AboutData const& /* about */ )
{
    return os;
}
}
