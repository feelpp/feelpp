# On "allume" les dépots
sudo rm /etc/yum.repos.d/rpm*
sudo rm /etc/yum.repos.d/slc6-cernonly*
sudo sed -i "s/enabled=0/enabled=1/g" -e /etc/yum.repos.d
sudo yum upgrade

# Install cd cmake28
sudo yum install cmake28
export cmake=cmake28

#Installation de différentes lib requises
sudo yum install texinfo glibc-devel.i686 icu.x86_64 libicu-devel.x86_64 libmpc.x86_64 mpfr.x86_64 gmp-devel.x86_64 mpfr-devel.x86_64 libmpc-devel.x86_64

#Config général: tout est installé dans $home/work
export workdir=$HOME/work
mkdir $workdir

# Installation de gcc 4.7
mkdir $workdir/_gcc
cd $workdir/_gcc
wget ftp://ftp.mpi-sb.mpg.de/pub/gnu/mirror/gcc.gnu.org/pub/gcc/releases/gcc-4.7.2/gcc-4.7.2.tar.bz2
tar xjf gcc-4.7.2.tar.bz2
mkdir BUILD_DIR
cd BUILD_DIR
../gcc-4.7.2/configure --prefix=$workdir/gcc
make -j4
make -j4 install
#if gcc is set after $PATH, gcc --version -> 4.4
#if gcc is set before $PATH, gcc --version -> 4.7
export LD_LIBRARY_PATH=$workdir/gcc
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/vhuber/work/gcc/lib/../lib64
export CC=$workdir/gcc/bin/gcc
export CXX=$workdir/gcc/bin/g++
export FC=$workdir/gcc/bin/gfortran
export F90=$workdir/gcc/bin/gfortran

#installation d'openmPI
mkdir $workdir/_openmpi
cd $workdir/_openmpi
wget http://www.open-mpi.org/software/ompi/v1.6/downloads/openmpi-1.6.3.tar.bz2
tar xjf openmpi-1.6.3.tar.bz2
cd openmpi-1.6.3
./configure CFLAGS=-m64 CXXFLAGS=-m64 FFLAGS=-m64 FCFLAGS=-m64 --prefix=$workdir/openmpi
make -j all install
export PATH=$workdir/openmpi/bin:$PATH
export LD_LIBRARY_PATH=$workdir/openmpi/libmake:$LD_LIBRARY_PATH

#Installation de BOOST
mkdir $workdir/_boost
cd $workdir/_boost
wget http://ignum.dl.sourceforge.net/project/boost/boost/1.49.0/boost_1_49_0.tar.bz2
tar xjf boost_1_49_0.tar.bz2
cd boost_1_49_0
rm user-config.jam
echo "using mpi ;" >> user-config.jam
echo "" >> user-config.jam
./bootstrap.sh
./bjam -j4 install \
      --layout=tagged \
      --prefix=$workdir/boost\
      --user-config=user-config.jam \
      variant=release \
      threading=single,multi \
      link=static,shared
export Boost_DIR=$workdir/boost


#PETSC
export  PETSC_ARCH="arch-linux2-c-opt"
export  PETSC_DIR="$workdir/petsc-3.3-p5"
cd $workdir
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.3-p5.tar.gz
tar xzf petsc-3.3-p5.tar.gz
cd $PETSC_DIR
./configure --with-shared-libraries=1 --with-debugging=0 COPTFLAGS='-O3 -march=p4 -mtune=p4' FOPTFLAGS='-O3 -qarch=p4 -qtune=p4'
make all
make test

#GMSH
export GMSH_DIR=$workdir/gmsh
mkdir $GMSH_DIR
cd $workdir
wget http://geuz.org/gmsh/src/gmsh-2.6.1-source.tgz
tar xzf gmsh-2.6.1-source.tgz
cd gmsh-2.6.1-source
mkdir BuildDir
cd BuildDir
cmake -DCMAKE_INSTALL_PREFIX=$GMSH_DIR -DCMAKE_BUILD_TYPE=release -DCMAKE_PREFIX_PATH=$workdir/openmpi/include ..
make lib
make shared
make -j8 install
export PATH=$PATH:$GMSH_DIR/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GMSH_DIR/lib

# On clone le dépot GIT
cd $workdir
git clone https://github.com/feelpp/feelpp.git
mkdir feelppLib
mkdir feelpp_build_dir 
cd feelpp_build_dir
cmake $workdir/feelpp \
      -DCMAKE_BUILD_TYPE=release \
      -DCMAKE_INSTALL_PREFIX=$workdir/feelppLib
make -j4 install
