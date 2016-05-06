#!/bin/bash

set -e

echo '--- setup module'
if [ -f /usr/share/modules/init/sh ]; then
    . /usr/share/modules/init/sh
    echo '--- setup module done'
fi
export FEELPP_HPCNAME=sd-87660
export FEELPP_CONFIG_PATH=/opt/config
export FEELPP_SHARE_PATH=/opt/software/install
export FEELPP_MODULE_PATH=/opt/config/modules

echo '--- setup feelpp module'
if [ -f  /opt/config/etc/feelpprc.sh  ]; then
    . /opt/config/etc/feelpprc.sh
    module load libs/boost-1.60.feelpp   libs/petsc-3.6.3.feelpp    mpi/openmpi-1.10.2.feelpp
    module load libs/hdf5-1.8.16.feelpp  science/gmsh-2.10.1.feelpp
fi



