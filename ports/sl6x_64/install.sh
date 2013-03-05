# On "allume" les dépots (les deux premiers ne fonctionnent pas)
sudo rm /etc/yum.repos.d/rpm*
sudo rm /etc/yum.repos.d/slc6-cernonly*
sudo sed -i "s/enabled=0/enabled=1/g" -e /etc/yum.repos.d
sudo yum upgrade

# Install de quelques paquets nécessaires.
sudo yum -y update
sudo yum -y update yum
sudo yum -y install cmake28 texinfo glibc-devel.i686 icu.x86_64 libicu-devel.x86_64 libmpc.x86_64 mpfr.x86_64 gmp-devel.x86_64 mpfr-devel.x86_64 libmpc-devel.x86_64 zip unzip

#Pour installer ailleurs: changer la variable $workdir (s'assurer d'avoir les droits)
#espace nécessaire: 10Go (possibilité de minimiser beaucoup en supprimant les repertoires de compilation)
alias cmake=`cmake28`
export workdir=/usr/local/feelpp

export gccVersion=4.6.3
export gccDir=$workdir/gcc-$gccVersion
export openmpiDir=$workdir/openmpi/gcc-$gccVersion/
export boostDir=$workdir/boost/gcc-$gccVersion/
export gmshDir=$workdir/gmsh/gcc-$gccVersion/
export petscDir=$workdir/petsc/gcc-$gccVersion/
export feelppDir=$workdir/feelppDir/gcc-$gccVersion/

module load libs/mkl13

mkdir $workdir
mkdir $workdir/openmpi
mkdir $workdir/boost
mkdir $workdir/gmsh
mkdir $workdir/petsc
mkdir $workdir/feelppDir
mkdir $openmpiDir
mkdir $boostDir
mkdir $gmshDir
mkdir $petscDir
mkdir $feelppDir

export nbProc=`cat /proc/cpuinfo |grep processor | tail -n 1 | awk '{print $3}'|awk '{print "("$1"+1)/2"}' | bc`

#télécharge, compile et install gcc dans $workdir/gcc
./gcc47.sh
export PATH=$gccDir/bin:$PATH
export LD_LIBRARY_PATH=$gccDir/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$gccDir/lib64:$LD_LIBRARY_PATH
export CC=$gccDir/bin/gcc
export CXX=$gccDir/bin/g++
export FC=$gccDir/bin/gfortran
export F90=$gccDir/bin/gfortran
#télécharge, compile et install openmpi dans $workdir/openmpi
./openmpi.sh
export PATH=$openmpiDir/bin:$PATH
export LD_LIBRARY_PATH=$openmpiDir/lib:$LD_LIBRARY_PATH
#télécharge, compile et install boost dans $workdir/boost
./boost.sh
export BOOSTROOT=$boostDir
export PATH=$BOOSTROOT/include:$PATH
export LD_LIBRARY_PATH=$BOOSTROOT/lib:$LD_LIBRARY_PATH
#télécharge, compile et install gmsh dans $workdir/gmsh
export GMSH_DIR=$gmshDir
./gmsh.sh
export PATH=$GMSH_DIR/bin:$PATH
export LD_LIBRARY_PATH=$GMSH_DIR/lib:$LD_LIBRARY_PATH
#télécharge, compile et install petsc dans $workdir/petsc-version
./petsc.sh
export PETSC_DIR=$workdir/petsc/gcc-$gccVersion
#télécharge, compile et install feelpp dans $workdir/feelppLib
./feelpp.sh
