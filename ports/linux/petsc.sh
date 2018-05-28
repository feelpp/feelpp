#! /bin/bash

prefix=/opt/petsc
export PETSC_DIR=`pwd`
./configure --prefix=$1 \
    --with-shared-libraries=1 \
    --with-cc=/usr/bin/mpicc \
    --with-cxx=/usr/bin/mpic++ \
    --with-fc=/usr/bin/mpif90 \
    --enable-debugging=0\
    --COPTFLAGS="-O3"\
    --CXXOPTFLAGS="-O3"\
    --FOPTFLAGS="-O3"\
    --download-metis \
    --download-parmetis \
    --download-blacs \
    --download-scalapack \
    --download-fblaslapack \
    --download-mumps \
    --download-ptscotch \
    --download-ml \
    --download-hypre \
    --download-suitesparse\
    --download-superlu

#    --download-fftw \
