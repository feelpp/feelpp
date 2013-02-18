#PETSC
export PETSC_ARCH="arch-linux2-c-opt"
export PETSC_DIR="$petscDir/petsc-3.3-p5"
mkdir $workdir/_petsc
cd $workdir/_petsc
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.3-p5.tar.gz
tar xzf petsc-3.3-p5.tar.gz
cd petsc-3.3-p5
./configure --with-shared-libraries=1 \
	--with-debugging=0 \
	COPTFLAGS='-O3 -march=p4 -mtune=p4' FOPTFLAGS='-O3 -qarch=p4 -qtune=p4' \
	--prefix=$PETSC_DIR
make all
make test
make install
