#! /bin/bash

wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.4.3.tar.gz
tar xzf petsc-lite-3.4.3.tar.gz
cd petsc-3.4.3
export PETSC_DIR=`pwd`
./configure --download-f-blas-lapack --download-mpich\
    --with-shared-libraries=1 \
    --download-blacs \
    --download-scalapack \
    --download-mumps \
    --download-ptscotch
