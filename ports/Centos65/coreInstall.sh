#!/bin/bash
#Pour installer ailleurs: changer la variable $workdir (s'assurer d'avoir les droits)
#espace nécessaire: 10Go (possibilité de minimiser beaucoup en supprimant les repertoires de compilation)
export basedir=$HOME
export workdir=libs

export boostDir=$basedir/$workdir/boost
export gmshDir=$basedir/$workdir/gmsh
export petscDir=$basedir/$workdir/petsc

mkdir -p $boostDir/src
mkdir -p $gmshDir/src
mkdir -p $petscDir/src

#télécharge, compile et install boost
./boost.sh
export BOOSTROOT=$boostDir
#télécharge, compile et install gmsh
export GMSH_DIR=$gmshDir
./gmsh.sh
export PATH=$GMSH_DIR/bin:$PATH
#télécharge, compile et install petsc
./petsc.sh
export PETSC_DIR=$petscDir
#télécharge, compile et install feelpp
./feelpp.sh
