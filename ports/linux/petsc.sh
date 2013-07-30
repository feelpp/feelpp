#! /bin/bash

prefix=/usr
export PETSC_DIR=`pwd`
./configure --prefix=/opt/petsc/$1 \
    --with-shared-libraries=1 \
    --with-cc=/usr/bin/mpicc \
    --with-cxx=/usr/bin/mpic++ \
    --with-fc=/usr/bin/mpif90 \
    --download-metis \
    --download-parmetis \
    --download-blacs \
    --download-scalapack \
    --download-f-blas-lapack \
    --download-mumps \
    --download-ptscotch
