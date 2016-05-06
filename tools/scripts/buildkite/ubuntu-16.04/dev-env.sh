#!/bin/bash

set -e

if [ -f  /opt/config/etc/feelpprc.sh  ]; then
    . /opt/config/etc/feelpprc.sh
fi
module load libs/boost-1.60.feelpp   libs/petsc-3.6.3.feelpp    mpi/openmpi-1.10.2.feelpp
module load libs/hdf5-1.8.16.feelpp  science/gmsh-2.10.1.feelpp


