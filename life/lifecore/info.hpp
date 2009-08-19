/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-02-19

  Copyright (C) 2008,2009 Université de Grenoble 1
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
   \file info.hpp

   The file was created from KDE/kdelibs/kdecore/kdeversion.hpp and
   accomodated to Life needs.

   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-02-19
 */
#ifndef _LIFE_VERSION_H_
#define _LIFE_VERSION_H_

#define LIFE_MAKE_VERSION( a,b,c ) (((a) << 16) | ((b) << 8) | (c))
#define LIFE_IS_VERSION(a,b,c) ( LIFE_VERSION >= LIFE_MAKE_VERSION(a,b,c) )

/**
 * Namespace for general LIFE functions.
 */
namespace Life
{
/**
 * \class Info
 * \brief information provider for versioning and installation directories
 * \author Christophe Prud'homme
 *
 */
class Info
{
public:
    /**
     * Build id is the Life compilation date/time (e.g 2005/12/07 - 18:18:09)
     *
     * This is a unique id for each release.
     * It permits Life to check if a more recent version of itself exist.
     * Format: year/month/day hour:minutes:seconds
     */
    static unsigned long long buildId();

    /**
     * Revision number (Subversion revision).
     */
    static unsigned long long revision();

    /**
     * Returns the encoded number of LIFE's version, see the LIFE_VERSION macro.
     * In contrary to that macro this function returns the number of the actually
     * installed LIFE version, not the number of the LIFE version that was
     * installed when the program was compiled.
     * @return the version number, encoded in a single unsigned long long
     * @since 0.7
     */
    static unsigned long long version();

    /**
     * Returns the major number of LIFE's version, e.g.
     * 0 for LIFE 0.7
     * @return the major version number
     * @since 0.7
     */
    static unsigned int versionMajor();

    /**
     * Returns the minor number of LIFE's version, e.g.
     * 7 for LIFE 0.7.0
     * @return the minor version number
     * @since 0.7
     */
    static unsigned int versionMinor();

    /**
     * Returns the micro number of LIFE's version, e.g.
     * 0 for LIFE 0.7.0
     * @return the extra information
     * @since 0.7
     */
    static unsigned int versionMicro();

    /**
     * Returns the LIFE version as string, e.g. "0.7.0".
     * @return the LIFE version. You can keep the string forever
     * @since 0.7
     */
    static char const* versionString();

    /**
     * \brief prefix directory
     *
     * A prefix used in constructing the default values of the
     * variables listed below. The default value of prefix should be
     * /usr/local. When building the complete GNU system, the prefix
     * will be empty and /usr will be a symbolic link to /.
     */
    static char const* prefix();

    /**
     * \brief datadir directory
     *
     * The directory for installing idiosyncratic read-only
     * architecture-independent data files for this program
     */
    static char const* datadir();

private:
    Info();
    Info( Info const& );
    Info& operator=( Info const& );
}; // Info

}

#endif // _LIFE_VERSION_H_
