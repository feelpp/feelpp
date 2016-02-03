#! /bin/bash

unset SLEPC_DIR

configure \
    --prefix=$WORKDIR/packages-install/slepc-3.6.1 --with-clean=1 \
    --with-arpack=1 --with-arpack-dir=/usr/local/slepc-3.5.3/lib --with-arpack-flags="-lparpack -larpack"

