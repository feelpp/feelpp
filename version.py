# -*- Mode: python coding: latin-1 -*-
###  version.py ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2006-05-25
#
#  Copyright (C) 2006 Université Joseph Fourier
#
# Distributed under the GPL(GNU Public License):
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#
import os, time, sys, re, types, shutil, stat, httplib, ftplib, urllib, urllib2, zipfile, errno, tarfile, popen2

import string, fnmatch, glob

def VCReadLines(infolines):
        """
        Get the info from the revision system be it svk or svn
        """

        for line in infolines:
                line = line.split(':')
                if len(line) > 1:
                        revLine = line[0].strip().lower()
                        if revLine[2:] == 'vision':
                                return eval(line[1].strip())
        return 0

def LifeGetRevision(env):
        """
        Gets the revision number using the shell command 'svn info'

        You need the official Subversion command-line client installed and referenced inside your PATH.
        Under Windows TortoiseSVN is not enough.

        Example of a 'svn info' output:
        Path: .
        URL: http://tanguy_k@dev.openwengo.com/svn/openwengo/trunk
        Repository UUID: 30a43799-04e7-0310-8b2b-ea0d24f86d0e
        Revision: 3448
        Node Kind: directory
        Schedule: normal
        Last Changed Author: tanguy_k
        Last Changed Rev: 3447
        Last Changed Date: 2005-12-09 20:24:12 +0100 (Fri, 09 Dec 2005)

        @rtype integer
        @return Subversion revision number or 0 if an error occured
        """

#        infocmd = ' info ' + '\"' + env['root_build_dir'] + '\"'
        infocmd = ' info ' + env.GetLaunchDir()
        info = os.popen( 'svn'+infocmd )
        infolines=info.readlines()
        if len(infolines):
                return VCReadLines( infolines )
        else:
                info = os.popen( 'svk'+infocmd )
                infolines=info.readlines()
                if len(infolines):
                        return VCReadLines( infolines )
                else:
                        print "Invalid repository(must be svn or svk)\n"
        return 0

def LifeGetCurrentDateTime(env):
        """
        Gets the current date-time using the format YYYYMMDDHHMMSS.

        Example: 20051210125200 (meaning 2005-12-10 12:52:00)

        @rtype integer
        @return current date-time
        """

        return eval(time.strftime('%Y%m%d%H%M%S', time.gmtime()))

def LifeGetSourcePath(env):
        return os.path.abspath(os.path.join(env['root_build_dir'], env.Dir('.').srcnode().path))


def LifeGetBuildId(env):
    return {
        'build_id': LifeGetCurrentDateTime(env),
        'version_major': 1,
        'version_minor': 0,
        'version_micro': 0,
        'version_string': '1.0.0pre2',
        'version' :  (((1) << 16) | ((0) << 8) | (0)),
        'revision': LifeGetRevision(env),
        'package': 'life',
        }

def make_str(strlist, prefix):
    thestr=''
    for i in strlist:
            if i[0]!='#':
                    thestr = thestr + ' ' + prefix+i
    return thestr

def make_version(env):
    LifeBuildId=LifeGetBuildId(env)
    return LifeBuildId['version_string']+'-r'+str(LifeBuildId['revision'])+'-'+str(LifeBuildId['build_id'])
