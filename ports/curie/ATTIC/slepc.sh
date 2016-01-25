#! /bin/bash

. ~/.bash_profile
module load c++/gnu/4.5.1   
prefix=$WORKDIR/local-gcc45
export SLEPC_DIR=`pwd`
./configure --prefix=$prefix 
