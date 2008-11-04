# -*- Mode: python -*-
###  SConstruct --- top level

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
from boost import *
from version import *
import os, copy

def default_dir( try_this, default_should_exists = '/usr/include'):
     if os.path.exists( try_this ):
          return try_this
     else :
          return default_should_exists
#
# get the environment
#
opts = Options(['soptions.cache', 'SOptions'])
opts.AddOptions(

    # scons flags
    BoolOption('verbose', 'scons in verbose mode', 0),

    # compilation flags
    BoolOption('warnings', 'compilation with -Wall and similar', 1),
    BoolOption('shared', 'enable shared libraries compilation', 0),
    BoolOption('static', 'enable static libraries compilation', 1),
    BoolOption('nobuildid', 'compilation without build id', 0),

    EnumOption('debug', 'debug output and symbols', 'no',
               allowed_values=('yes', 'no', 'full'),
               map={}, ignorecase=0),  # case sensitive

    EnumOption('opt', 'Optimized output and symbols', '0',
               allowed_values=('-1', '0', '1', '2', '3'),
               map={}, ignorecase=0),  # case sensitive

    # directories
    PathOption('prefix', 'prefix directory where lifev should be installed', '/usr'),

    # boost
    PathOption('boost_inc_dir', 'where the root of Boost headers is installed', '/usr/include'),
    PathOption('boost_lib_dir', 'where the Boost libraries are installed', '/usr/lib'),

    # qd
    PathOption('qd_inc_dir', 'where the QD headers are installed', default_dir('/usr/include/qd')),
    PathOption('qd_lib_dir', 'where the QD libraries are installed', '/usr/lib'),

    # mpi
    BoolOption('disable_mpi', 'Disable MPI', 0 ),
    PathOption('mpi_inc_dir', 'where the MPI headers are installed', default_dir('/usr/include/mpi')),
    PathOption('mpi_lib_dir', 'where the MPI libraries are installed', '/usr/lib'),

    # metis/parmetis
    PathOption('metis_inc_dir', 'where the Metis headers are installed', default_dir('/usr/include/metis')),
    PathOption('metis_lib_dir', 'where the Metis libraries are installed', '/usr/lib'),
    PathOption('parmetis_inc_dir', 'where the ParMetis headers are installed', default_dir('/usr/include/parmetis')),
    PathOption('parmetis_lib_dir', 'where the ParMetis libraries are installed', '/usr/lib'),

    # Superlu
    PathOption('superlu_inc_dir', 'where the Superlu headers are installed', default_dir('/usr/include/superlu')),
    PathOption('superlu_lib_dir', 'where the Superlu libraries are installed', '/usr/lib'),

    # ufsparse
    PathOption('ufsparse_inc_dir', 'where the UFSparse headers are installed', default_dir('/usr/include/ufsparse')),
    PathOption('ufsparse_lib_dir', 'where the UFSparse libraries are installed', '/usr/lib'),

    # vtk
    PathOption('vtk_inc_dir', 'where the VTK headers are installed', default_dir('/usr/include/vtk')),
    PathOption('vtk_lib_dir', 'where the VTK libraries are installed', '/usr/lib'),

    # petsc
    PathOption('petsc_inc_dir', 'where the Petsc headers are installed', default_dir('/usr/include/petsc')),
    PathOption('petsc_lib_dir', 'where the Petsc libraries are installed', '/usr/lib'),

    # trilinos
    BoolOption('disable_trilinos', 'Disable Trilinos', 0 ),
    PathOption('trilinos_inc_dir', 'where the Trilinos headers are installed', default_dir('/usr/include/trilinos/6.0')),
    PathOption('trilinos_lib_dir', 'where the Trilinos libraries are installed', '/usr/lib'),

    # buildir
    PathOption('root_build_dir', 'the root build directory', os.getcwd()),

    # enable/disable life parts
    EnumOption('testsuite', 'enable testsuite compilation and execution', 'compile',
               allowed_values=('compile', 'execute', 'no'),
               map={}, ignorecase=0),  # case sensitive

    BoolOption('examples', 'enable examples compilation', 0 ),
)
# the doxygen.py file needs to be in toolpath

env = Environment(tools = ["default", "doxygen"], toolpath = '.',
                  CPPPATH=['#', '#/boost', '#/contrib/gmm/include'],
                  LIBPATH = '',
                  options=opts)
theenv = {}
if os.environ.has_key('CXX'):
     env['CXX'] = os.environ['CXX']
if os.environ.has_key('CC'):
     env['CC'] = os.environ['CC']
if os.environ.has_key('HOME'):
     env['ENV']['HOME'] = os.environ['HOME']
if os.environ.has_key('PATH'):
     env['ENV']['PATH'] = os.environ['PATH']

Help(opts.GenerateHelpText(env))
opts.Save('soptions.cache', env)

def CheckPKGConfig(context, version):
     context.Message( 'Checking for pkg-config... ' )
     ret = context.TryAction('pkg-config --atleast-pkgconfig-version=%s' % version)[0]
     context.Result( ret )
     return ret

def CheckPKG(context, name):
     context.Message( 'Checking for %s... ' % name )
     ret = context.TryAction('pkg-config --exists \'%s\'' % name)[0]
     context.Result( ret )
     return ret

#
# Start configure checks
#
conf = Configure( env,
                  config_h='lifeconfig.h',
                  custom_tests = { 'CheckBoost' : CheckBoost,
                                   'CheckPKGConfig' : CheckPKGConfig,
                                   'CheckPKG' : CheckPKG
                                   }
                  )
if not conf.CheckPKGConfig('0.15.0'):
     print 'pkg-config >= 0.15.0 not found.'
#     Exit(1)

#
# Required libs
#

# Blas/Lapack
conf.CheckLib('blas','dgemm_',language='c++')
conf.CheckLib('lapack','dgesvd_',language='c++')

# Boost
boost_libs=['signals','program_options', 'date_time','filesystem','serialization','test']
if not conf.CheckBoost(boost_libs):
     print "Boost Libraries ",boost_libs," must be installed"
     print "They are mandatory for LifeV"
     print "See http://www.boost.org for more details about boost"
     Exit(1)

#
# Optional libs
#

# MPI
if not env['disable_mpi']:
     lastCPPPATH= copy.copy(conf.env['CPPPATH'])
     conf.env.Append( CPPPATH = conf.env['mpi_inc_dir'] )
     if conf.CheckCXXHeader('mpi.h'):
          conf.env.AppendUnique( CPPDEFINES = ['HAVE_MPI'] );
     else:
          conf.env.Replace( CPPPATH = lastCPPPATH )
          lastLIBPATH= copy.copy(conf.env['LIBPATH'])
          conf.env.Append( LIBPATH = conf.env['mpi_lib_dir'] )
          if conf.CheckLib('mpi','MPI_Init',language='c++'):
               conf.env.AppendUnique(LIBS = ['mpi++','mpi'])
          else:
               conf.env.Replace( LIBPATH = lastLIBPATH )


# Metis
lastCPPPATH= copy.copy(conf.env['CPPPATH'])
conf.env.Append( CPPPATH = conf.env['metis_inc_dir'] )
if not conf.CheckCXXHeader('metis.h'):
    conf.env.Replace( CPPPATH = lastCPPPATH )
lastLIBPATH= copy.copy(conf.env['LIBPATH'])
conf.env.Append( LIBPATH = conf.env['metis_lib_dir'] )
if conf.CheckLib('metis','METIS_PartGraphRecursive',language='c++'):
    conf.env.AppendUnique(LIBS = ['metis'])
else:
     conf.env.Replace( LIBPATH = lastLIBPATH )

# ParMetis
lastCPPPATH= copy.copy(conf.env['CPPPATH'])
conf.env.Append( CPPPATH = conf.env['parmetis_inc_dir'] )
if not conf.CheckCXXHeader('parmetis.h'):
    conf.env.Replace( CPPPATH = lastCPPPATH )
lastLIBPATH= copy.copy(conf.env['LIBPATH'])
conf.env.Append( LIBPATH = conf.env['parmetis_lib_dir'] )
if conf.CheckLib('parmetis','ParMETIS_V3_PartGeomKway',language='c++'):
    conf.env.AppendUnique(LIBS = ['parmetis'])
else:
     conf.env.Replace( LIBPATH = lastLIBPATH )


# QD
lastLIBPATH= copy.copy(conf.env['LIBPATH'])
lastCPPPATH= copy.copy(conf.env['CPPPATH'])
conf.env.Append( LIBPATH = conf.env['qd_lib_dir'] )
if not conf.CheckLib('qd','fpu_fix_start',language='c++'):
     conf.env.Replace( LIBPATH = lastLIBPATH )

conf.env.Append( CPPPATH = conf.env['qd_inc_dir']  )
if not conf.CheckCXXHeader('qd.h'):
    conf.env.Replace( CPPPATH = lastCPPPATH )

# Superlu check
lastLIBPATH= copy.copy(conf.env['LIBPATH'])
lastCPPPATH= copy.copy(conf.env['CPPPATH'])
conf.env.Append( LIBPATH = conf.env['superlu_lib_dir'] )
conf.env.Append( CPPPATH = conf.env['superlu_inc_dir']  )
if conf.CheckCXXHeader('slu_Cnames.h') and conf.CheckLib('superlu','superlu_free',language='c++'):
#     conf.env['umfpack_libs']=['amd', 'umfpack', 'lapack', 'blas']
     conf.env.AppendUnique( CPPDEFINES = ['GMM_USES_SUPERLU'] );
else :
     conf.env.Replace( CPPPATH = lastCPPPATH )
     conf.env.Replace( LIBPATH = lastLIBPATH )

# UMFPACK check
lastLIBPATH= copy.copy(conf.env['LIBPATH'])
lastCPPPATH= copy.copy(conf.env['CPPPATH'])
conf.env.Append( LIBPATH = conf.env['ufsparse_lib_dir'] )
conf.env.Append( CPPPATH = conf.env['ufsparse_inc_dir']  )
conf.CheckLib('amd','amd_postorder',language='c++')
if conf.CheckCXXHeader('umfpack.h') and conf.CheckLib('umfpack','umfpack_di_free_symbolic',language='c++'):
     conf.env['umfpack_libs']=['amd', 'umfpack', 'lapack', 'blas']
     conf.env.AppendUnique( CPPDEFINES = ['HAVE_UMFPACK'] );
else :
     conf.env.Replace( CPPPATH = lastCPPPATH )
     conf.env.Replace( LIBPATH = lastLIBPATH )

# VTK check
lastLIBPATH= copy.copy(conf.env['LIBPATH'])
lastCPPPATH= copy.copy(conf.env['CPPPATH'])
conf.env.Append( LIBPATH = conf.env['vtk_lib_dir'] )
conf.env.Append( CPPPATH = conf.env['vtk_inc_dir']  )
if conf.CheckLibWithHeader('vtkCommon','vtkConfigure.h','c++'):
#     conf.env.AppendUnique(LIBS = ['vtkCommon'])
     conf.env['vtk_libs']=['vtkRendering', 'vtkGraphics', 'vtkImaging', 'vtkCommon']
     conf.env.AppendUnique( CPPDEFINES = ['HAVE_VTK'] );
else:
     conf.env.Replace( CPPPATH = lastCPPPATH )
     conf.env.Replace( LIBPATH = lastLIBPATH )

# Petsc check
lastLIBPATH= copy.copy(conf.env['LIBPATH'])
lastCPPPATH= copy.copy(conf.env['CPPPATH'])
conf.env.Append( LIBPATH = conf.env['petsc_lib_dir'] )
conf.env.Append( CPPPATH = conf.env['petsc_inc_dir']  )

if not conf.CheckCXXHeader('petsc.h'):
    conf.env.Replace( CPPPATH = lastCPPPATH )

if (conf.CheckLib('petsc','PetscInitialize') and
    conf.CheckLib('petscvec','VecCreate') and
    conf.CheckLib('petscmat','MatCreate') and
    conf.CheckLib('petscksp','KSPCreate') and
    conf.CheckLib('petscdm','DMInitializePackage') ):
#     conf.env.AppendUnique(LIBS = ['vtkCommon'])
     conf.env['petsc_libs']=['petsc']
#     conf.env.AppendUnique( CPPDEFINES = ['HAVE_PETSC'] );
else:
     conf.env.Replace( CPPPATH = lastCPPPATH )
     conf.env.Replace( LIBPATH = lastLIBPATH )

# Trilinos check
if not env['disable_trilinos']:
	lastLIBPATH= copy.copy(conf.env['LIBPATH'])
	lastCPPPATH= copy.copy(conf.env['CPPPATH'])
	conf.env.Append( LIBPATH = conf.env['trilinos_lib_dir'] )
	conf.env.Append( CPPPATH = conf.env['trilinos_inc_dir']  )

	if conf.CheckLibWithHeader('epetra','Epetra_config.h','c++', autoadd=0):
		conf.env.AppendUnique( CPPDEFINES = ['HAVE_TRILINOS_EPETRA'] );
		conf.env.PrependUnique( LIBS = ['epetra'] )
	else:
     		conf.env.Replace( CPPPATH = lastCPPPATH )
     		conf.env.Replace( LIBPATH = lastLIBPATH )

	if conf.CheckLibWithHeader('epetraext','EpetraExt_ConfigDefs.h','c++', autoadd=0):
     		conf.env.PrependUnique( LIBS = ['epetraext'] )
     		conf.env.AppendUnique( CPPDEFINES = ['HAVE_TRILINOS_EPETRAEXT'] );
	else:
     		conf.env.Replace( CPPPATH = lastCPPPATH )
     		conf.env.Replace( LIBPATH = lastLIBPATH )
	if conf.CheckLibWithHeader('triutils','Triutils_config.h','c++', autoadd=0):
     		conf.env.PrependUnique( LIBS = ['triutils'] )
     		conf.env.AppendUnique( CPPDEFINES = ['HAVE_TRILINOS_TRIUTILS'] );
	else:
     		conf.env.Replace( CPPPATH = lastCPPPATH )
     		conf.env.Replace( LIBPATH = lastLIBPATH )

	if conf.CheckLibWithHeader('teuchos','Teuchos_config.h','c++', autoadd=0):
     		conf.env.PrependUnique( LIBS = ['teuchos'] )
     		conf.env.AppendUnique( CPPDEFINES = ['HAVE_TRILINOS_TEUCHOS'] );
	else:
     		conf.env.Replace( CPPPATH = lastCPPPATH )
     		conf.env.Replace( LIBPATH = lastLIBPATH )

	if conf.CheckLibWithHeader('aztecoo','AztecOO_config.h','c++', autoadd=0):
     		conf.env.PrependUnique( LIBS = ['aztecoo'] )
     		conf.env.AppendUnique( CPPDEFINES = ['HAVE_TRILINOS_AZTECOO','HAVE_AZTECOO_TEUCHOS'] );
	else:
     		conf.env.Replace( CPPPATH = lastCPPPATH )
     		conf.env.Replace( LIBPATH = lastLIBPATH )
	if conf.CheckLibWithHeader('amesos','Amesos_config.h','c++', autoadd=0):
     		conf.env.PrependUnique( LIBS = ['amesos'] )
     		conf.env.AppendUnique( CPPDEFINES = ['HAVE_TRILINOS_AMESOS'] );
	else:
     		conf.env.Replace( CPPPATH = lastCPPPATH )
     		conf.env.Replace( LIBPATH = lastLIBPATH )

	if conf.CheckLibWithHeader('ifpack','Ifpack_config.h','c++', autoadd=0):
                conf.env.PrependUnique( LIBS = ['ifpack'] )
        	conf.env.AppendUnique( CPPDEFINES = ['HAVE_TRILINOS_IFPACK'] );
	else:
     		conf.env.Replace( CPPPATH = lastCPPPATH )
     		conf.env.Replace( LIBPATH = lastLIBPATH )


	if conf.CheckLibWithHeader('anasazi','Anasazi_config.h','c++', autoadd=0):
     		conf.env.PrependUnique( LIBS = ['anasazi'] )
     		conf.env.AppendUnique( CPPDEFINES = ['HAVE_TRILINOS_ANASAZI'] );
	else:
     		conf.env.Replace( CPPPATH = lastCPPPATH )
     		conf.env.Replace( LIBPATH = lastLIBPATH )

	if conf.CheckLibWithHeader('ml','ml_config.h','c++', autoadd=0):
		conf.env.PrependUnique( LIBS = ['ml'] )
		conf.env.AppendUnique( CPPDEFINES = ['HAVE_TRILINOS_ML'] );
	else:
		conf.env.Replace( CPPPATH = lastCPPPATH )
		conf.env.Replace( LIBPATH = lastLIBPATH )

# we are finished with configuration
env = conf.Finish()
#print "LIBS=", env['LIBS']

# unit test builder
def builder_unit_test(target, source, env):
     if 'test' in COMMAND_LINE_TARGETS:
          for i in range(0,len(source)):
               print "================================================================================"
               app = str(source[i].abspath)
               print "Running Test ", app
               if os.spawnl(os.P_WAIT, app, app)==0:
                    open(str(target[0]),'w').write("PASSED\n")
               else:
                    return 1

# Create a builder for tests
bld = Builder(action = builder_unit_test)
env.Append(BUILDERS = {'Test' :  bld})


#
# export environment
#
Export("env")

#examine options
if env["warnings"]:
        env.Append( CXXFLAGS   = " -Wall" );

profile='sconsstd'
if env["debug"] == 'yes':
        env.AppendUnique( CXXFLAGS   = ['-g3', '-O0'] );
        env.AppendUnique( CPPDEFINES = ['DEBUG=1']);
        profile='sconsdebug'
elif env["opt"] != '-1':
        env.AppendUnique( CXXFLAGS   = ['-O'+env["opt"]] );
        if env["opt"] == '3':
            env.AppendUnique( CXXFLAGS   = ['-ftree-vectorize'] );
        env.AppendUnique( CPPDEFINES = ['NDEBUG']);
        profile='sconsopt'+env["opt"]
else :
        env.AppendUnique( CXXFLAGS   = ['-g', '-O2'] );
        env.AppendUnique( CPPDEFINES = ['NDEBUG']);

#print "build_dir=",env['root_build_dir'] + "/" + profile
env['os_sysname']=os.uname()[0]
env['os_release']=os.uname()[2]
env['os_machine']=os.uname()[4]
env['mach_build_dir']=env['os_sysname']+os.sep+env['os_machine']
env['build_dir']=env['root_build_dir'] + os.sep + profile + os.sep + env['mach_build_dir']

# make sure build_dir actually exists
if not os.path.exists(env['build_dir']):
     os.makedirs(env['build_dir'])

if env['verbose']:
     print "Summary of build variables"
     print "=========================="
     print " o-         CC:", env['CC']
     print " o-    CCFLAGS:", env['CCFLAGS']
     print " o-        CXX:", env['CXX']
     print " o-   CXXFLAGS:", env['CXXFLAGS']
     print " o-    CPPPATH:", env['CPPPATH']
     print " o- CPPDEFINES:", env['CPPDEFINES']
     print " o-    LIBPATH:", env['LIBPATH']
     print " o-       LIBS:", env['LIBS']
     print "=========================="

     print "Summary of environment"
     print "======================"
     print " o-        profile: ", profile
     print " o- root_build_dir: ", env['root_build_dir']
     print " o-      build_dir: ", env['build_dir']
     print "======================"

# Generate Doxyfile with proper version
#doxyfiletmpl=open('Doxyfile.tmpl','r')
#doxy=doxyfiletmpl.read()+'PROJECT_NUMBER='+make_version(env)+'\n'
#doxyfile=open('Doxyfile','w')
#doxyfile.write(doxy)
#doxyfile.close();

SConscript('SConscript', build_dir=env['build_dir'], duplicate=0)



