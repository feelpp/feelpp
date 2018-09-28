#! /bin/bash

PETSC_VERSION=${1:-3.8.4}
PREFIX=${2:-/opt/petsc/${PETSC_VERSION}}
wget -nc http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-${PETSC_VERSION}.tar.gz 
tar -xf petsc-${PETSC_VERSION}.tar.gz 
cd petsc-${PETSC_VERSION} 
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
            --prefix=${PREFIX} 
make 
make install 
                    
# prefix=/opt/petsc
# export PETSC_DIR=`pwd`
# ./configure --prefix=$1 \
#     --with-shared-libraries=1 \
#     --with-cc=/usr/bin/mpicc \
#     --with-cxx=/usr/bin/mpic++ \
#     --with-fc=/usr/bin/mpif90 \
#     --enable-debugging=0\
#     --COPTFLAGS="-O3"\
#     --CXXOPTFLAGS="-O3"\
#     --FOPTFLAGS="-O3"\
#     --download-metis \
#     --download-parmetis \
#     --download-blacs \
#     --download-scalapack \
#     --download-fblaslapack \
#     --download-mumps \
#     --download-ptscotch \
#     --download-ml \
#     --download-hypre \
#     --download-suitesparse\
#     --download-superlu

# #    --download-fftw \
