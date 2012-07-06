#!/bin/sh

module load cmake/2.8.7   
module load gcc/4.6.3  
module load openmpi/gcc/1.4.5  

export prefix=/softs/cemracs/petsc/3.3
export MPI=/softs/openmpi/gcc/1.4.5/

./configure --with-python --with-debugging=0 \
    --with-c-support=1 --with-c++-support=1  \
    --with-shared-libraries=1 --with-mpi=1 --PETSC_ARCH=linux \
    --prefix=${destroot}${prefix}/lib/petsc \
    --with-cc=${MPI}/bin/mpicc \
    --with-cxx=${MPI}/bin/mpic++ \
    --with-mpiexec=${MPI}/bin/mpiexec \
        --download-umfpack=1 \
        --download-ml \
        --download-metis \
        --download-parmetis \
        --download-blacs \
        --download-scalapack \
        --download-f-blas-lapack \
        --download-mumps \
        --download-pastix \
        --download-ptscotch \
    --with-fc=${MPI}/bin/mpif90 \
    --LIBS=-lstdc++

