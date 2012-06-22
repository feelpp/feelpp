#! /bin/sh

module load cmake/2.8.3
module load atlas/3.9.72 
cmake \
 -DBLAS_atlas_LIBRARY=/usr/local/atlas-3.9.72/lib/libatlas.a\
 -DBLAS_cblas_LIBRARY=/usr/local/atlas-3.9.72/lib/libcblas.a\
 -DBLAS_f77blas_LIBRARY=/usr/local/atlas-3.9.72/lib/libf77blas.a\
 -DCMAKE_INSTALL_PREFIX=$WORKDIR/local-gcc45 \
 $SCRATCHDIR/Gmsh/gmsh-2.6.0-source/
