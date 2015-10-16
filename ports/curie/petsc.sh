#! /bin/bash

unset PETSC_DIR
unset PETSC_ARCH

./configure \
    PETSC_ARCH=linux-gnu-intel --prefix=$WORKDIR/packages-install/petsc-3.6.2/ \
    --LDFLAGS="-lrt -L/usr/local/intel-14.0.3.174/14.0.3.174/mkl/lib/intel64 -lmkl_rt -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -liomp5 -lm" \
    --with-clean=1 --with-shared-libraries=1 --with-debugging=0 \
    --with-cc=$MPI_ROOT/bin/mpicc --with-fc=$MPI_ROOT/bin/mpif90 --with-cxx=$MPI_ROOT/bin/mpic++ \
    --with-blas-lapack-dir=/usr/local/intel-14.0.3.174/14.0.3.174/mkl/lib/intel64 \
    --with-scalapack=1 --with-scalapack-include=/usr/local/intel-14.0.3.174/14.0.3.174/mkl/include --with-scalapack-lib="[-L/usr/local/intel-14.0.3.174/14.0.3.174/mkl/lib/intel64 -lmkl_rt -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -liomp5 -lm]" \
    --with-blacs=1 --with-blacs-include=/usr/local/intel-14.0.3.174/14.0.3.174/mkl/include --with-blacs-lib="[-L/usr/local/intel-14.0.3.174/14.0.3.174/mkl/lib/intel64 -lmkl_rt -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -liomp5 -lm]" \
    --with-hwloc=1 --with-hwloc-dir=/usr/local/hwloc-1.7.1 \
    --with-ptscotch=1 --with-ptscotch-dir=/usr/local/petsc-3.5.2/ --with-ml=1 --with-ml-include=/usr/local/petsc-3.5.2//include --with-ml-lib="[/usr/local/petsc-3.5.2//lib/libml.a,-lmpi_cxx,-lstdc++]" --with-hypre=1 --with-hypre-include=/usr/local/petsc-3.5.2//include --with-hypre-lib="[-L/usr/local/petsc-3.5.2//lib,-lHYPRE]" \
    --with-mumps=1 --with-mumps-include=/usr/local/petsc-3.5.2//include --with-mumps-lib="[-L/usr/local/petsc-3.5.2//lib,-lcmumps,-ldmumps,-lsmumps,-lzmumps,-lmumps_common,-lparmetis,-lmetis,-lpord,-lptesmumps,-lptscotch,-lptscotcherr]" \
    --with-metis=1 --with-metis-include=/usr/local/petsc-3.5.2//include --with-metis-lib="[-L/usr/local/petsc-3.5.2//lib,-lmetis]" --with-parmetis=1 --with-parmetis-include=/usr/local/petsc-3.5.2//include --with-parmetis-lib="[-L/usr/local/petsc-3.5.2//lib,-lparmetis,-lmetis]" \
    --with-suitesparse=1 --with-suitesparse-include=/usr/local/petsc-3.5.2//include --with-suitesparse-lib="[-L/usr/local/petsc-3.5.2//lib,-lumfpack,-lklu,-lcholmod,-lbtf,-lccolamd,-lcolamd,-lcamd,-lamd,-lmetis,-lsuitesparseconfig,rt]" \
    --with-hdf5=1  --with-hdf5-dir=$HDF5_ROOT
