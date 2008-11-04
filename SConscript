# -*- Mode: python coding: latin-1 -*-
###  SConscript ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
#       Date: 2006-05-25
#
#  Copyright (C) 2006 EPFL
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
import os,fnmatch,string,glob
import doxygen

Import('env')



projects=[
    'boost/libs/log/src',
    'life',
    'life/lifecore',
    'life/lifealg',
    'life/lifemesh',
    'life/lifefilters',
    'life/lifediscr',
    ]

# compile eventually the testsuite
testsuite=[]
if env['testsuite'] == 'execute' or env['testsuite'] == 'compile':
    testsuite=[
        'testsuite/lifecore',
        'testsuite/lifealg',
        'testsuite/lifemesh',
        'testsuite/lifepoly',
        'testsuite/lifediscr',
#        'testsuite/lifefilters',
        'testsuite/',
        ]
    if env['testsuite'] == 'execute':
        print "LifeV Testsuite..."
        print "LifeV Testsuite...done."

# compile eventually the examples
examples=[]
if env['examples']:
    examples=[
        'examples',
        ]
SConscript(dirs=projects+testsuite+examples)

#
# Install
#
def SGlob(pattern):
    path = string.replace(GetBuildPath('SConscript'), 'SConscript', '')

    result = []
    for i in glob.glob(path + pattern):
        result.append(string.replace(i, path, ''))
    return result

def Glob( pattern = '*.*', thedir = '.' ):
    import os, fnmatch
    files = []

    for file in os.listdir( Dir(thedir).abspath ):
        if fnmatch.fnmatch(file, pattern) :
            files.append( os.path.join( thedir, file ) )

    return files

if 'install' in COMMAND_LINE_TARGETS:
    # install gmm
    env.Install( os.path.join( env['prefix'], 'include/gmm'), SGlob('contrib/gmm/include/*.h') )
    env.Alias('install', os.path.join( env['prefix'], 'include/gmm') )

    # install boost

    cmdline='find ' + Dir('.').srcnode().abspath +'/boost/boost -type d -print'
    info = os.popen( cmdline  )
    infolines=info.readlines()
    for boostdir in infolines:
        boostdir=string.strip(boostdir)
        boostinstdir=string.replace(boostdir,Dir('.').srcnode().abspath +'/boost/', '')
        env.Install( os.path.join( env['prefix'], 'include/'+boostinstdir), SGlob('boost/'+boostinstdir+'/*.h*') )
        env.Alias('install', os.path.join( env['prefix'], 'include/'+boostinstdir) )

    env.Install( os.path.join( env['prefix'], 'lib'), Glob('*.a',env['build_dir']+'/boost/libs/log/src') )
    env.Alias('install', os.path.join( env['prefix'], 'lib') )

    # install life
    for i in ['life',
              'life/lifecore',
              'life/lifealg',
              'life/lifemesh',
              'life/lifepoly',
              'life/lifetime',
              'life/lifediscr',
              'life/lifefilters',
              'life/lifevf',]:
        env.Install( os.path.join( env['prefix'], 'include/'+i), SGlob(i+'/*.hpp') )
        env.Alias('install', os.path.join( env['prefix'], 'include/'+i))
        if os.path.exists(env['build_dir']+'/'+i):
            env.Install( os.path.join( env['prefix'], 'lib'), Glob('*.a',env['build_dir']+'/'+i) )
            env.Alias('install', os.path.join( env['prefix'], 'lib'))

    env.Install( os.path.join( env['prefix'], 'share/pkgconfig'), Glob('*.pc',env['build_dir']+'/life') )
    env.Alias('install', os.path.join( env['prefix'], 'share/pkgconfig'))
    env.Install( os.path.join( env['prefix'], 'include/'), SGlob('*.h') )
    env.Alias('install', os.path.join( env['prefix'], 'include/'))


#
# Test
if 'test' in COMMAND_LINE_TARGETS:
    env.Alias( 'test','testsuite/lifecore/test.lifecore.passed' )
    env.Alias( 'test','testsuite/lifemesh/test.lifemesh.passed' )
    env.Alias( 'test','testsuite/lifepoly/test.passed' )
    env.Alias( 'test','testsuite/lifediscr/test.passed' )

if 'testcore' in COMMAND_LINE_TARGETS:
    env.Alias( 'test','testsuite/lifecore/test.lifecore.passed' )

if 'testalg' in COMMAND_LINE_TARGETS:
    env.Alias( 'test','testsuite/lifealg/test.passed' )

if 'testmesh' in COMMAND_LINE_TARGETS:
    env.Alias( 'test','testsuite/lifemesh/test.lifemesh.passed' )

if 'testpoly' in COMMAND_LINE_TARGETS:
    env.Alias( 'test','testsuite/lifepoly/test.passed' )

if 'testfem' in COMMAND_LINE_TARGETS:
    env.Alias( 'test','testsuite/lifediscr/test.passed' )

#
# Doxygen
#env.Doxygen("Doxyfile")
