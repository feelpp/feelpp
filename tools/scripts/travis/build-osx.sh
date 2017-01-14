#! /bin/bash

brew update
brew install openmpi boost petsc slepc gmsh
mkdir build
cd build
../configure -r --cxxflags="-O2 -DNDEBUG"
make

