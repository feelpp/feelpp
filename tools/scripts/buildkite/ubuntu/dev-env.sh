#!/bin/bash

set -e

NPROCS=10
echo '--- apt-get update'
apt-get -qq update

echo '--- install apt-add-repository'
apt-get -y --force-yes install software-properties-common python-software-properties

echo '--- added gcc toolchain and clang repo'
apt-add-repository 'deb http://ppa.launchpad.net/ubuntu-toolchain-r/test/ubuntu trusty main'
apt-add-repository 'deb http://llvm.org/apt/trusty/ llvm-toolchain-trusty-3.7 main'
apt-cache search g++-4.9
apt-cache search clang-3.7

echo '--- added openturns repo'
curl http://ubuntu.openturns.org/openturns.org-repo.key | apt-key add -

echo '--- apt-get install'
apt-get -y --force-yes install \
        xauth cmake flex g++-4.9 clang-3.7 gfortran git ipython openmpi-bin pkg-config \
            wget bison sudo \
            libbz2-dev \
            libboost-all-dev libboost-mpi-dev \
            automake autoconf libtool \
            libopenblas-dev libcln-dev libcppunit-dev \
            libeigen3-dev liblapack-dev libmpfr-dev \
            libslepc3.4.2-dev \
            libopenmpi-dev libann-dev libglpk-dev \
            python-dev libhwloc-dev libvtk6-dev libpcre3-dev \
            libhdf5-openmpi-dev libeigen3-dev libcgal-dev \
            python-numpy python-vtk6 python-six python-ply \
            python-h5py python-urllib3 xterm tmux screen python-openturns

echo '--- apt-get clean'
apt-get clean
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*



export CXX=clang++-3.7
export CC=clang-3.7
export FEELPP_DEP_INSTALL_PREFIX=/usr/local

# Boost
export BOOST_VERSION=1.59.0
export BOOST_DIR=boost_1_59_0
echo '--- compiling/installing boost $BOOST_VERSION'

if ! [ -f ${FEELPP_DEP_INSTALL_PREFIX}/boost-${BOOST_VERSION} ]; then
wget http://sourceforge.net/projects/boost/files/boost/${BOOST_VERSION}/${BOOST_DIR}.tar.bz2/download -O ${BOOST_DIR}.tar.bz2 \
    && tar xjf ${BOOST_DIR}.tar.bz2 \
    && cd ${BOOST_DIR} \
    && echo "using mpi ;" >> user-config.jam \
    && echo "" >> user-config.jam \
    && ./bootstrap.sh \
    && ./bjam -j$NPROCS install \
      --layout=tagged \
      --prefix=${FEELPP_DEP_INSTALL_PREFIX} \
      --user-config=user-config.jam \
      variant=release \
      threading=multi \
      link=shared
touch ${FEELPP_DEP_INSTALL_PREFIX}/boost-${BOOST_VERSION}
fi

# Install PETSc from source
export PETSC_VERSION=3.6.3
echo '--- compiling/installing PETSc $PETSC_VERSION'
if ! [ -f ${FEELPP_DEP_INSTALL_PREFIX}/petsc-${PETSC_VERSION} ]; then
   cd /tmp && \
    wget -nc http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-${PETSC_VERSION}.tar.gz && \
    tar -xf petsc-${PETSC_VERSION}.tar.gz && \
    cd petsc-${PETSC_VERSION} && \
    ./configure --COPTFLAGS="-O2" \
                --CXXOPTFLAGS="-O2" \
                --FOPTFLAGS="-O2" \
                --with-blas-lib=/usr/lib/libopenblas.a --with-lapack-lib=/usr/lib/liblapack.a \
                --with-c-support \
                --with-debugging=0 \
                --with-shared-libraries \
                --download-suitesparse \
                --download-scalapack \
                --download-metis \
                --download-parmetis \
                --download-ptscotch \
                --download-hypre \
                --download-mumps \
                --download-blacs \
                --download-spai \
                --download-ml \
                --prefix=${FEELPP_DEP_INSTALL_PREFIX} && \
     make && \
     make install && \
     touch ${FEELPP_DEP_INSTALL_PREFIX}/petsc-${PETSC_VERSION}
   rm -rf /tmp/*
fi

# Install SLEPc from source
export SLEPC_VERSION=3.6.3
echo '--- compiling/installing SLEPc $SLEPC_VERSION'
if ! [ -f ${FEELPP_DEP_INSTALL_PREFIX}/slepc-${SLEPC_VERSION} ]; then
cd /tmp && \
    export PETSC_DIR=${FEELPP_DEP_INSTALL_PREFIX} && \
    wget -nc http://www.grycap.upv.es/slepc/download/download.php?filename=slepc-${SLEPC_VERSION}.tar.gz -O slepc-${SLEPC_VERSION}.tar.gz && \
    tar -xf slepc-${SLEPC_VERSION}.tar.gz && \
    cd slepc-${SLEPC_VERSION} && \
    ./configure --prefix=${FEELPP_DEP_INSTALL_PREFIX} && \
    make && \
    make install && \
    rm -rf /tmp/*
touch ${FEELPP_DEP_INSTALL_PREFIX}/slepc-${SLEPC_VERSION}
fi
export SLEPC_DIR=${FEELPP_DEP_INSTALL_PREFIX}
export PETSC_DIR=${FEELPP_DEP_INSTALL_PREFIX}

# Gmsh
# check to see if protobuf folder is empty
export GMSH_VERSION=2.12.0
echo '--- compiling/installing GMSH ${GMSH_VERSION}'
if ! [ -f ${FEELPP_DEP_INSTALL_PREFIX}/gmsh-${GMSH_VERSION} ]; then
cd /tmp \
    && wget http://www.geuz.org/gmsh/src/gmsh-${GMSH_VERSION}-source.tgz \
    && tar xvzf gmsh-${GMSH_VERSION}-source.tgz \
    && cd gmsh-${GMSH_VERSION}-source \
    && mkdir build \
    && cd build \
    && cmake \
        -DCMAKE_CXX_COMPILER=${CXX} \
        -DCMAKE_C_COMPILER=${CC} \
        -DCMAKE_INSTALL_PREFIX=${FEELPP_DEP_INSTALL_PREFIX} \
        -DCMAKE_BUILD_TYPE=release \
        -DENABLE_BUILD_LIB=ON \
        -DENABLE_BUILD_SHARED=ON \
        -DENABLE_BUILD_DYNAMIC=ON \
        -DENABLE_MPI=OFF \
        -DENABLE_MUMPS=OFF \
        -DENABLE_OPENMP=ON  \
        .. \
    && make -j$NPROCS \
    && make install \
    && rm -rf /tmp/*
touch ${FEELPP_DEP_INSTALL_PREFIX}/gmsh-${GMSH_VERSION}
fi

# ParaView
export PARAVIEW_VERSION=4.4.0
echo '--- compiling/installing PARAVIEW $PARAVIEW_VERSION'
if ! [ -f ${FEELPP_DEP_INSTALL_PREFIX}/paraview-${PARAVIEW_VERSION} ]; then
cd /tmp \
    && wget "http://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v4.4&type=source&os=all&downloadFile=ParaView-v${PARAVIEW_VERSION}-source.tar.gz" -O ParaView-v${PARAVIEW_VERSION}-source.tar.gz \
    && tar xvzf ParaView-v${PARAVIEW_VERSION}-source.tar.gz \
    && cd ParaView-v${PARAVIEW_VERSION}-source \
    && mkdir build \
    && cd build \
    && cmake \
        -DCMAKE_CXX_COMPILER=${CXX} \
        -DCMAKE_C_COMPILER=${CC} \
        -DCMAKE_INSTALL_PREFIX=${FEELPP_DEP_INSTALL_PREFIX} \
        -DCMAKE_BUILD_TYPE=Release \
        -DPARAVIEW_BUILD_QT_GUI=OFF \
        -DPARAVIEW_ENABLE_CATALYST=ON \
        -DPARAVIEW_ENABLE_PYTHON=ON \
        -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON \
        -DPARAVIEW_USE_MPI=ON \
        .. \
    && make -j$NPROCS \
    && make install \
    && rm -rf /tmp/*
touch ${FEELPP_DEP_INSTALL_PREFIX}/paraview-${PARAVIEW_VERSION}
fi
