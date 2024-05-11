import os
import re
import sys
import platform
import subprocess

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = [
            '-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON',
            '-DCMAKE_CXX_COMPILER=@CMAKE_CXX_COMPILER@',
            '-DCMAKE_C_COMPILER=@CMAKE_C_COMPILER@',
            '-DCMAKE_INSTALL_PREFIX=@CMAKE_INSTALL_PREFIX@',
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
            '-DPython3_EXECUTABLE=' + sys.executable,
            '-DFEELPP_DIR=@FEELPP_DIR@',
            '-DFEELPP_TOOLBOXES_DIR=@FEELPP_TOOLBOXES_DIR@']

        if @PYFEELPP_SETUP_HAS_PARAVIEW_CMAKE_ARGS@ == 1 :
            cmake_args += @PYFEELPP_SETUP_PARAVIEW_CMAKE_ARGS@

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        env['PYTHONPATH'] = '@PYFEELPP_TOOLBOXES_PYTHONPATH@'

        self.build_temp=self.build_temp+"-"+ext.name
        print(self.build_temp)
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)

setup(
    name='pyfeelpptoolboxes',
    version='@FEELPP_VERSION_STRING@',
    author='Christophe Prud\'homme',
    author_email='christophe.prudhomme@feelpp.org',
    description='PyFeel++-Toolboxes: Python bindings for Feel++ Toolboxes',
    long_description='',
    package_dir={ 'pyfeelpptoolboxes': '@CMAKE_CURRENT_SOURCE_DIR@/pyfeelpptoolboxes' },
    packages=['pyfeelpptoolboxes',
              'pyfeelpptoolboxes.modelcore','pyfeelpptoolboxes.modelmesh',
              #'pyfeelpptoolboxes.advection',
              'pyfeelpptoolboxes.fluid','pyfeelpptoolboxes.solid',
              'pyfeelpptoolboxes.electric',
              'pyfeelpptoolboxes.heat',
              #'pyfeelpptoolboxes.maxwell',
              'pyfeelpptoolboxes.hdg',
              #'pyfeelpptoolboxes.levelset','pyfeelpptoolboxes.multifluid',
              'pyfeelpptoolboxes.thermoelectric',
              #'pyfeelpptoolboxes.fsi',
              
    ],
    install_requires=['pyfeelpp'],
    ext_modules=[CMakeExtension('_modelcore','@CMAKE_CURRENT_SOURCE_DIR@/pyfeelpptoolboxes/modelcore'),
                 CMakeExtension('_modelmesh','@CMAKE_CURRENT_SOURCE_DIR@/pyfeelpptoolboxes/modelmesh'),
                 #CMakeExtension('_advection','@CMAKE_CURRENT_SOURCE_DIR@/pyfeelpptoolboxes/advection'),
                 CMakeExtension('_fluid','@CMAKE_CURRENT_SOURCE_DIR@/pyfeelpptoolboxes/fluid'),
                 CMakeExtension('_solid','@CMAKE_CURRENT_SOURCE_DIR@/pyfeelpptoolboxes/solid'),
                 CMakeExtension('_electric','@CMAKE_CURRENT_SOURCE_DIR@/pyfeelpptoolboxes/electric'),
                 CMakeExtension('_heat','@CMAKE_CURRENT_SOURCE_DIR@/pyfeelpptoolboxes/heat'),
                 #CMakeExtension('_maxwell','@CMAKE_CURRENT_SOURCE_DIR@/pyfeelpptoolboxes/maxwell'),
                 CMakeExtension('_hdg','@CMAKE_CURRENT_SOURCE_DIR@/pyfeelpptoolboxes/hdg'),
                 #CMakeExtension('_levelset','@CMAKE_CURRENT_SOURCE_DIR@/pyfeelpptoolboxes/levelset'),
                 #CMakeExtension('_multifluid','@CMAKE_CURRENT_SOURCE_DIR@/pyfeelpptoolboxes/multifluid'),
                 CMakeExtension('_thermoelectric','@CMAKE_CURRENT_SOURCE_DIR@/pyfeelpptoolboxes/thermoelectric'),
                 #CMakeExtension('_fsi','@CMAKE_CURRENT_SOURCE_DIR@/pyfeelpptoolboxes/fsi'),
                 
    ],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
)
