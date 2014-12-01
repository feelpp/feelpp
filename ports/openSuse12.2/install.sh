# On "allume" les dépots (les deux premiers ne fonctionnent pas)
export workdir=/opt/feelpp

export openmpiDir=$workdir/openmpi
export boostDir=$workdir/boost
export gmshDir=$workdir/gmsh
export petscDir=$workdir/petsc

mkdir $workdir
mkdir $openmpiDir
mkdir $boostDir
mkdir $gmshDir
mkdir $petscDir

export nbProc=`cat /proc/cpuinfo |grep processor | tail -n 1 | awk '{print $3}'|awk '{print "("$1"+1)/2"}' | bc`

#télécharge, compile et install openmpi dans $workdir/openmpi
./openmpi.sh
export PATH=$openmpiDir/bin:$PATH
export LD_LIBRARY_PATH=$openmpiDir/lib64:$LD_LIBRARY_PATH
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
export PETSC_DIR=$workdir/petsc
#télécharge, compile et install feelpp dans $workdir/feelppLib
#./feelpp.sh
