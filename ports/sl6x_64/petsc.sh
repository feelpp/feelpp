#PETSC
export PETSC_ARCH="arch-linux2-c-opt"
export PETSC_DIR="$workdir/petsc-3.3-p5"
cd $workdir
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.3-p5.tar.gz
tar xzf petsc-3.3-p5.tar.gz
cd $PETSC_DIR
./configure --with-shared-libraries=1 --with-debugging=0 COPTFLAGS='-O3 -march=p4 -mtune=p4' FOPTFLAGS='-O3 -qarch=p4 -qtune=p4'
make all
make test
