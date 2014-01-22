#PETSC
cd $petscDir/src
wget -c http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.4.3.tar.gz
tar xzf petsc-3.4.3.tar.gz
cd petsc-3.4.3
./configure --with-shared-libraries=1 \
  --with-debugging=0 \
  COPTFLAGS='-O3 -march=p4 -mtune=p4' FOPTFLAGS='-O3 -qarch=p4 -qtune=p4' \
  --prefix=$petscDir \
  --download-umfpack=1 \
  --download-ml \
  --download-metis \
  --download-parmetis \
  --download-blacs \
  --download-scalapack \
  --download-f-blas-lapack \
  --download-mumps \
  --download-pastix \
  --download-ptscotch
make all
make test
make install
