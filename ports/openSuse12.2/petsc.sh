#PETSC
mkdir $workdir/_petsc
mkdir $petscDir
mkdir $PETSC_DIR
cd $workdir/_petsc
wget -c http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.5.2.tar.gz
if [ ! -d "petsc-3.5.2" ]; then
  tar xzf petsc-3.5.2.tar.gz
fi
cd petsc-3.5.2
./configure --with-shared-libraries=1 \
  --with-debugging=0 \
  COPTFLAGS='-O3 -march=p4 -mtune=p4' FOPTFLAGS='-O3 -qarch=p4 -qtune=p4' \
  --prefix=$petscDir \
  --with-cc=${openmpiDir}/bin/mpicc \
  --with-cxx=${openmpiDir}/bin/mpic++ \
  --with-mpiexec=${openmpiDir}/bin/mpiexec \
  --with-fc=${openmpiDir}/bin/mpifort \
  --download-suitesparse=1 \
  --download-ml \
  --download-metis \
  --download-parmetis \
  --download-blacs \
  --download-scalapack \
  --download-fblaslapack \
  --download-mumps \
  --download-pastix \
  --download-ptscotch
make all
make test
make install
