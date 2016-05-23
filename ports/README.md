L'ordre de compilation doit être respecté

Gmsh ne dépend ni de Petsc ni d'openMPI !

## OpenMPI 
```
openmpiDir=/data/software/install/openmpi-1.10.0
# Enable fortran binding: --enable-mpi-fortran=all
./configure CFLAGS=-m64 CXXFLAGS=-m64 FFLAGS=-m64 FCFLAGS=-m64 --prefix=$openmpiDir
make -j all install
export PATH=$openmpiDir/bin:$PATH
export LD_LIBRARY_PATH=$openmpiDir/lib:$LD_LIBRARY_PATH
```

## Boost
```
boostDir=/data/software/install/boost-1.59.0
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

## Gmsh

```
cmake -DCMAKE_INSTALL_PREFIX=/data/software/install/gmsh-2.10.1 -DENABLE_MPI=OFF -DENABLE_BUILD_LIB=ON -DENABLE_BUILD_SHARED=ON -DCMAKE_BUILD_TYPE=release -DENABLE_PETSC=OFF ..
make -j$64 lib
make -j64 shared
make -j64 install
export GMSH_DIR=/data/software/install/gmsh-2.10.1 
```

## PETSc
```
 ./configure --with-shared-libraries=1  --with-debugging=0  --COPTFLAGS='-O3' --CXXOPTFLAGS='-O3' --FOPTFLAGS='-O3'  --prefix=/data/software/install/petsc-3.6.1/openmpi-1.10.0/  --with-cc=`which mpicc`  --with-cxx=`which mpic++` --with-fc=`which mpif90` --with-mpiexec=`which mpiexec`  --download-suitesparse=1  --download-ml  --download-metis  --download-parmetis  --download-blacs  --download-scalapack  --download-fblaslapack  --download-mumps  --download-hypre  --download-ptscotch --download-elemental --download-elemental-shared=1 --with-cxx-dialect=C++11
 make PETSC_DIR=/data/software/src/petsc-3.6.1 PETSC_ARCH=arch-linux2-c-opt all
 make PETSC_DIR=/data/software/src/petsc-3.6.1 PETSC_ARCH=arch-linux2-c-opt install
 export PETSC_DIR=/data/software/install/petsc-3.6.1/openmpi-1.10.0/
```
## Slepc
```
 ./configure --prefix=/data/software/install/slepc-3.6.1/openmpi-1.10.0/
 make SLEPC_DIR=$PWD PETSC_DIR=/data/software/install/petsc-3.6.1/openmpi-1.10.0/
 make SLEPC_DIR=/data/software/src/slepc-3.6.1 PETSC_DIR=/data/software/install/petsc-3.6.1/openmpi-1.10.0/ install
 // Optionnal
 make SLEPC_DIR=/data/software/install/slepc-3.6.1/openmpi-1.10.0 PETSC_DIR=/data/software/install/petsc-3.6.1/openmpi-1.10.0/ PETSC_ARCH="" test
 export SLEPC_DIR=/data/software/install/slepc-3.6.1/openmpi-1.10.0
 ```
 
## hdf5
```
# For codes using fortran, add FC=`which mpif90` and the "--enable-fortran --enable-fortran2003" options
CC=`which mpicc` CXX=`which mpic++` ./configure --enable-parallel --prefix=/data/software/install/hdf5/1.8.15-patch1/gcc-4.9.0/openmpi-1.10 --enable-build-all --enable-production
make install
```

## ParaView 
* Build for using DISPLAY export over SSH (see https://github.com/feelpp/feelpp/blob/develop/ports/linux/ParaView.sh for more information)
```
cmake /data/software/src/ParaView/ParaView-v5.0.1-source -DBUILD_TESTING=OFF -DVTK_RENDERING_BACKEND=OpenGL -DPARAVIEW_ENABLE_CATALYST=ON -DPARAVIEW_ENABLE_PYTHON=ON -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON -DPARAVIEW_USE_MPI=ON -DCMAKE_INSTALL_PREFIX=/data/software/install/ParaView/5.0.1/gcc-4.9.0/openmpi-1.10
make -j $3
 ```
 
 * To use offscreen rendering (only available for NVIDIA with recent drivers > 355), please refer to this [wiki page](https://github.com/aancel/admin/wiki/Compile-ParaView-with-EGL-support-on-Ubuntu-14.04)
 
## VTK (not needed if ParaView is installed)
```
cmake /data/software/src/VTK/VTK5.10.1 -DCMAKE_INSTALL_PREFIX=/data/software/install/VTK/5.10.1/gcc-4.9.0/openmpi-1.10 -DCMAKE_BUILT_TYPE=Release -DVTK_USE_PARALLEL=ON -DBUILD_SHARED_LIBS=ON -DVTK_WRAP_PYTHON=ON -DVTK_USE_MPI=ON
make install
```

## FFTW
```
wget http://fftw.org/fftw-3.3.4.tar.gz && tar zxvf fftw-3.3.4.tar.gz && cd fftw-3.3.4
./configure --enable-mpi --enable-threads --enable-openmp --enable-shared --prefix=/data/software/install/fftw/3.3.4/gcc-4.9.0/openmpi-1.10
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
