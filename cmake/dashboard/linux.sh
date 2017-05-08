#!/bin/bash

#set -x

function run {
  echo $boost
  echo $petsc
  echo $ccomp
  if [[ $ccomp == "gcc"* ]]; then
    cxx=`echo $ccomp | sed 's/gcc/g++/g'`
  else
    cxx=`echo $ccomp | sed 's/clang/clang++/g'`
  fi
  export FEELPP_WORKDIR=/tmp/feel-$ccomp
  rm -rf $FEELPP_WORKDIR
  $COMMON,FEELPP_CXXNAME=$ccomp,FEELPP_CXX=`which $cxx`,FEELPP_C=`which $ccomp`
  rm -rf $FEELPP_WORKDIR
}

function petsc {
  # Set the petsc_dir
  for petsc in ${petsc_dir[@]}; do
    if [ "$petsc" != "default" ]
    then
      if [ -d $petsc ]; then
        echo "\texport PETSC_DIR=$petsc"
        run
      fi
    else
      echo "\tunset PETSC_DIR"
      run
    fi
  done
}


function boost {
  ## Set the boost_root
  for boost in ${boost_dir[@]}; do
    if [ "$boost" == "default" ]
    then
      echo "\tunset BOOST_ROOT"
      petsc 
    else
      if [ -d $boost ] ;
      then
        echo "\texport BOOST_ROOT=$boost"
        petsc
      fi
    fi # boost = default ?
  done # boost
}

# $1 provides the path to the toplevel source Feel++ directory
# e.g. $HOME/Devel/FEEL/feelpp.git
COMMON="ctest -VV -S $1/cmake/dashboard/testsuite.cmake,FEELPP_CTEST_CONFIG=$1/cmake/dashboard/feelpp.site.`hostname -s`.cmake,FEELPP_MODE=$2"

## Where are installed the custom libraries ?
lib_dir=${4:-"/data/software/install"}
## List of compiler (installed in the PATH) to test
compiler=("gcc-4.9 gcc-5.0 gcc-5 gcc-6 clang-3.6 clang-3.7 clang-3.8 clang-3.9");
## List of petsc installation to test (default = system one)
petsc_dir=("$lib_dir/petsc-3.6.3 $lib_dir/petsc-3.6.1 $lib_dir/petsc-3.7.0 default")
## List of boost installation to test (default = system one)
boost_dir=("$lib_dir/boost-1.58 $lib_dir/boost-1.59.0 default")

for ccomp in ${compiler[@]}; do
  echo "**"
  which $ccomp
  if [[ $? == 0 ]]; then
    boost
  fi
done #Compiler
