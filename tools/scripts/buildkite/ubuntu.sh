#!/bin/bash

set -e

echo '+++ install apt-add-repository'
apt-get install software-properties-common python-software-properties

echo '+++ added gcc toolchain and clang repo'
apt-add-repository 'deb http://ppa.launchpad.net/ubuntu-toolchain-r/test/ubuntu trusty main'
apt-add-repository 'deb http://llvm.org/apt/trusty/ llvm-toolchain-trusty-3.7 main'

echo '--- apt-get update'
apt-get -qq update

echo '--- apt-get install'
apt-get -y --force-yes install \
        xauth cmake flex g++-4.9 clang-3.7 gfortran git ipython openmpi-bin pkg-config \
            wget bison sudo \
            libbz2-dev \
            libboost-all-dev libboost-mpi-dev \
            automake autoconf libtool \
            libopenblas-dev libcln-dev libcppunit-dev \
            libeigen3-dev liblapack-dev libmpfr-dev \
            libslepc3.4.2-dev \
            libopenmpi-dev libann-dev libglpk-dev \
            python-dev libhwloc-dev libvtk6-dev libpcre3-dev \
            libhdf5-openmpi-dev libeigen3-dev libcgal-dev \
            python-numpy python-vtk6 python-six python-ply \
            python-h5py python-urllib3 xterm tmux screen

echo '--- apt-get clean'
apt-get clean
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
