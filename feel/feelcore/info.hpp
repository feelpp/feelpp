/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \file info.hpp

   The file was created from KDE/kdelibs/kdecore/kdeversion.hpp and
   accomodated to Feel needs.

   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-02-19
 */
#ifndef _FEELPP_VERSION_H_
#define _FEELPP_VERSION_H_

#define FEELPP_MAKE_VERSION( a,b,c ) (((a) << 16) | ((b) << 8) | (c))
#define FEELPP_IS_VERSION(a,b,c) ( FEELPP_VERSION >= FEELPP_MAKE_VERSION(a,b,c) )

/**
 * Namespace for general FEEL functions.
 */
namespace Feel
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
     * Build id is the Feel compilation date/time (e.g 2005/12/07 - 18:18:09)
     *
     * This is a unique id for each release.
     * It permits Feel to check if a more recent version of itself exist.
     * Format: year/month/day hour:minutes:seconds
     */
    static unsigned long long buildId();

    /**
     * Revision number (Subversion revision).
     */
    static unsigned long long revision();

    /**
     * Returns the encoded number of FEEL's version, see the FEELPP_VERSION macro.
     * In contrary to that macro this function returns the number of the actually
     * installed FEEL version, not the number of the FEEL version that was
     * installed when the program was compiled.
     * @return the version number, encoded in a single unsigned long long
     * @since 0.7
     */
    static unsigned long long version();

    /**
     * Returns the major number of FEEL's version, e.g.
     * 0 for FEEL 0.7
     * @return the major version number
     * @since 0.7
     */
    static unsigned int versionMajor();

    /**
     * Returns the minor number of FEEL's version, e.g.
     * 7 for FEEL 0.7.0
     * @return the minor version number
     * @since 0.7
     */
    static unsigned int versionMinor();

    /**
     * Returns the micro number of FEEL's version, e.g.
     * 0 for FEEL 0.7.0
     * @return the extra information
     * @since 0.7
     */
    static unsigned int versionMicro();

    /**
     * Returns the FEEL version as string, e.g. "0.7.0".
     * @return the FEEL version. You can keep the string forever
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

#endif // _FEELPP_VERSION_H_
