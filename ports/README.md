# Description of Feel++ dependencies and how to build them

The compilation order indicated on this page must be ensured for the dependencies to work. The dependencies are indicated to help you decided which modules you have to load on clusters or how to set up your environment for the next dependency.

## Compilers

If you need to recompile the compiler, here are adequate configure commands.

### GCC 

Depends on: gcc (System)

```

# Add the following flags if required:
# --with-gmp=...
# --with-mpfr=...
# --with-mpc=...

./configure --prefix=${INSTALL_DIR}/gcc/${VERSION} \
  --enable-bootstrap \
  --enable-shared \
  --enable-threads=posix \
  --enable-checking=release \
  --with-system-zlib \
  --enable-__cxa_atexit \
  --disable-libunwind-exceptions \
  --enable-gnu-unique-object \
  --enable-languages=c,c++,fortran \
  --disable-dssi \
  --enable-libgcj-multifile \
  --with-ppl \
  --with-cloog \
  --with-tune=generic \
  --disable-multilib \
  --with-arch_32=i686 \
```

### Clang 

Depends on: gcc (built for Feel++ or system)

```
cmake ../llvm -DCMAKE_BUILD_TYPE=Release -DGCC_INSTALL_PREFIX=${FEELPP_GCC_PATH} \
  -DCMAKE_C_COMPILER=${FEELPP_GCC_PATH}/bin/gcc -DCMAKE_CXX_COMPILER=${FEELPP_GCC_PATH}/bin/g++ \ 
  -DBUILD_SHARED_LIBS=ON -DLLVM_BUILD_TOOLS=OFF -DCLANG_INCLUDE_DOCS=OFF -DCLANG_INCLUDE_TESTS=OFF \ -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}/clang/${CLANG_VERSION}/gcc-${FEELPP_GCC_VERSION}
```

## MPI

### OpenMPI 

Depends on: gcc (built for Feel++ or system)

```
openmpiDir=${INSTALL_DIR}/openmpi/${OPENMPI_VERSION}/gcc-${FEELPP_GCC_VERSION}

# Enable fortran binding: --enable-mpi-fortran=all
./configure CFLAGS=-m64 CXXFLAGS=-m64 FFLAGS=-m64 FCFLAGS=-m64 --prefix=$openmpiDir

make -j all install
export PATH=$openmpiDir/bin:$PATH
export LD_LIBRARY_PATH=$openmpiDir/lib:$LD_LIBRARY_PATH
```

## Libraries

### Boost

Depends on: gcc (built for Feel++ or system), OpenMPI

```
boostDir=${INSTALL_DIR}/boost/${BOOST_VERSION}/gcc-${FEELPP_GCC_VERSION}/openmpi-${FEELPP_OPENMPI_VERSION}

rm user-config.jam
touch user-config.jam
echo "using mpi ;" >> user-config.jam
echo "" >> user-config.jam

./bootstrap.sh
./bjam -j64 install --layout=tagged --prefix=$boostDir --user-config=user-config.jam variant=release threading=single,multi link=static,shared

export BOOSTROOT=$boostDir
export PATH=$BOOSTROOT/include:$PATH
export LD_LIBRARY_PATH=$BOOSTROOT/lib:$LD_LIBRARY_PATH
```

### Gmsh 

Depends on: gcc (built for Feel++ or system)
(We do not use OpenMPI nor PETSc as dependencies here)

```
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}/Gmsh/${GMSH_VERSION}/gcc-${FEELPP_GCC_VERSION} \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_MPI=OFF -DENABLE_BUILD_LIB=ON -DENABLE_BUILD_SHARED=ON \ 
  -DCMAKE_BUILD_TYPE=Release -DENABLE_PETSC=OFF ..
  
make -j$64 lib
make -j64 shared
make -j64 install

export GMSH_DIR=${INSTALL_DIR}/Gmsh/${GMSH_VERSION}/gcc-${FEELPP_GCC_VERSION}
```

## PETSc

Depends on: gcc (built for Feel++ or system), OpenMPI

```
# Useful option to get the list of packages to download if you are on a cluster with no internet access (firewall ...), 
# but seems to be buggy for installation. So remove it when you have all the required packages
# See https://bitbucket.org/petsc/petsc/issues/133/install-issue-no-cached-configure-in-rdict
#  --with-packages-dir=../packages
./configure \ 
   --with-shared-libraries=1  --with-debugging=0  --COPTFLAGS='-O3' --CXXOPTFLAGS='-O3' --FOPTFLAGS='-O3' \ --prefix=${INSTALL_DIR}/PETSc/${PETSC_VERSION}/gcc-${FEELPP_GCC_VERSION}/openmpi-${FEELPP_OPENMPI_VERSION} \
   --with-cc=`which mpicc`  --with-cxx=`which mpic++` --with-fc=`which mpif90` --with-mpiexec=`which mpiexec` \  --download-suitesparse=1  --download-ml  --download-metis \
   --download-parmetis  --download-blacs  --download-scalapack \ 
   --download-fblaslapack  --download-mumps  --download-hypre \ 
   --download-ptscotch --download-elemental --download-elemental-shared=1 \
   --download-fftw=1
 
 make PETSC_DIR=`pwd` PETSC_ARCH=arch-linux2-c-opt all
 make PETSC_DIR=`pwd` PETSC_ARCH=arch-linux2-c-opt install
 
 export PETSC_DIR=${INSTALL_DIR}/PETSc/${PETSC_VERSION}/gcc-${FEELPP_GCC_VERSION}/openmpi-${FEELPP_OPENMPI_VERSION}
```
## Slepc

Depends on: gcc (built for Feel++ or system), OpenMPI, PETSc

```
 ./configure --prefix=${INSTALL_DIR}/SLEPc/${SLEPC_VERSION}/gcc-${FEELPP_GCC_VERSION}/openmpi-${FEELPP_OPENMPI_VERSION}/PETSc-${FEELPP_PETSC_VERSION}
 
 make SLEPC_DIR=`pwd` PETSC_DIR=${FEELPP_PETSC_DIR}
 make SLEPC_DIR=`pwd` PETSC_DIR=${FEELPP_PETSC_DIR} install
 // Optionnal
 make SLEPC_DIR=`pwd` PETSC_DIR=${FEELPP_PETSC_DIR} PETSC_ARCH="" test
 
 export SLEPC_DIR=/data/software/install/slepc-3.6.1/openmpi-1.10.0
 ```
 
## hdf5

Depends on: gcc (built for Feel++ or system), OpenMPI

```
# For codes using fortran, add FC=`which mpif90` and the "--enable-fortran --enable-fortran2003" options
CC=`which mpicc` CXX=`which mpic++` ./configure --enable-parallel --prefix=${INSTALL_DIR}/hdf5/${HDF5_VERSION}/gcc-${FEELPP_GCC_VERSION}/openmpi-${FEELPP_OPENMPI_VERSION} --enable-build-all --enable-production

make install
```

## ParaView 

Depends on: gcc (built for Feel++ or system), OpenMPI, Python, Qt (if you use `-DPARAVIEW_BUILD_QT_GUI=ON`)

* Build for using DISPLAY export over SSH (see https://github.com/feelpp/feelpp/blob/develop/ports/linux/ParaView.sh for more information)
```
cmake ${PARAVIEW_SOURCE_DIR} \
-DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF \
-DVTK_RENDERING_BACKEND=OpenGL -DPARAVIEW_ENABLE_CATALYST=ON -DPARAVIEW_ENABLE_PYTHON=ON \
-DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON -DPARAVIEW_USE_MPI=ON \
-DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}/ParaView/${PARAVIEW_VERSION}/gcc-${FEELPP_GCC_VERSION}/openmpi-${FEELPP_OPENMPI_VERSION}

make -j $3
 ```
 
 * To use offscreen rendering (only available for NVIDIA with recent drivers > 355), please refer to this [wiki page](https://github.com/aancel/admin/wiki/Compile-ParaView-with-EGL-support-on-Ubuntu-14.04)
 
## VTK (not needed if ParaView is installed)

Depends on: gcc (built for Feel++ or system), OpenMPI, Python

```
cmake ${VTK_SOURCE_DIR} \
-DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}/VTK/${VTK_VERSION}/gcc-${FEELPP_GCC_VERSION}/openmpi-${FEELPP_OPENMPI_VERSION} 
-DCMAKE_BUILD_TYPE=Release -DVTK_USE_PARALLEL=ON -DBUILD_SHARED_LIBS=ON -DVTK_WRAP_PYTHON=ON -DVTK_USE_MPI=ON

make install
```

<!--
## To Export
Be carreful, there is a conflict :
```
/usr/bin/ld: warning: libmpi.so.1, needed by /usr/lib/libvtkParallel.so.5.8.0, may conflict with libmpi.so.12
```
Do not compile Feel++ with VTK for that configuration, we need to build the module.

```
export PATH=/data/software/install/openmpi-1.10.0/bin:$PATH
export LD_LIBRARY_PATH=/data/software/install/openmpi-1.10.0/lib:$LD_LIBRARY_PATH
export BOOSTROOT=/data/software/install/boost-1.59.0
export PATH=$BOOSTROOT/include:$PATH
export LD_LIBRARY_PATH=$BOOSTROOT/lib:$LD_LIBRARY_PATH
export GMSH_DIR=/data/software/install/gmsh-2.10.1
export PETSC_DIR=/data/software/install/petsc-3.6.1/openmpi-1.10.0/
export SLEPC_DIR=/data/software/install/slepc-3.6.1/openmpi-1.10.0
```
-->
