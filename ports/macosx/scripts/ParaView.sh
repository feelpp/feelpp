#!/bin/bash
# ParaView building script for OS X
# author: Alexandre Ancel <alexandre.ancel@cemosis.fr>
# date: 01/09/2015

fail_on_error()
{
    echo "Error: $1"
    exit 1
}

echo "This script will build a ParaView version for OS X"
echo "This version will be compatible with Catalyst and ParaView plugins"
echo "You need several dependencies for ParaView, check http://www.paraview.org/Wiki/ParaView:Build_And_Install for more info"
echo ""
echo "For example, to install the dependencies with homebrew:"
echo "brew install openmpi"
echo "brew install --only-dependencies paraview"
echo ""

if [[ ! -d "./ParaView-v4.3.1-source" ]]; then
    tar zxvf ./ParaView-v4.3.1-source.tar.gz || fail_on_error "Please download ParaView from http://www.paraview.org/download/"
else
    echo "ParaView sources already extracted"
fi

cd ./ParaView-v4.3.1-source
mkdir -p build
mkdir -p install

SRCDIR=`pwd`
BUILDDIR=`pwd`/build
INSTALLDIR=`pwd`/install

echo "Building in $BUILDDIR, Installing in $INSTALLDIR"

cd $BUILDDIR

# if we didn't already built ParaView we do so
if [[ ! -d "$BUILDDIR/CMakeFiles/__macos_install" ]]; then
    # Configure the app
    # PARAVIEW_INSTALL_DEVELOPMENT_FILES is disabled in OS X, but when installing the headers are however installed and accessible "somewhere"
    cmake $SRCDIR -DPARAVIEW_ENABLE_CATALYST=ON -DPARAVIEW_ENABLE_PYTHON=ON -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON -DPARAVIEW_USE_MPI=ON -DMACOSX_APP_INSTALL_PREFIX=$INSTALLDIR || fail_on_error "Error during CMake step"
    # Do the initial build. It will take some time
    make -j 8 || fail_on_error "Error during make step"
    # Do the installation, even if the headers won't be installed
    # It will force the creation of the $BUILDDIR/CMakeFiles/__macos_install directory which will contain everything
    make -j 8 install || fail_on_error "Error during install step"
else
    echo "ParaView has already been built. If you want to rebuilt it, remove $BUILDDIR/CMakeFiles/__macos_install"
fi

echo "ParaView is built. To use Feel++ with the development headers"
echo "Export the following variables:"
#echo "export CMAKE_PREFIX_PATH=$BUILDDIR/CMakeFiles/__macos_install:\${CMAKE_PREFIX_PATH}"
#echo "export PATH=$BUILDDIR/CMakeFiles/__macos_install:\${PATH}"
#echo "export LD_LIBRARY_PATH=$BUILDDIR/CMakeFiles/__macos_install:\${LD_LIBRARY_PATH}"
echo "export PARAVIEW_DIR=$BUILDDIR/CMakeFiles/__macos_install"
