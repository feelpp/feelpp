#!/bin/sh

profname=gcc630

./install_cmake.sh $profname
./install_cln.sh $profname
#./boost.sh $profname
./install_petsc.sh $profname
./install_gmsh.sh $profname
