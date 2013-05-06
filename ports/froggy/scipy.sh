#! /bin/bash

. ~/.bash_profile

# if no loaded
module load python/2.7.3_gcc-4.6.2
module load lapack/3.4.0_gcc-4.6.2
module load swig/2.0.9_gcc-4.6.2

pushd ${HOME}/packages/scipy-0.12.0

export BLAS=/home/chabannes/packages/blas/BLAS/libblas.a
export LAPACK=/applis/ciment/v2/stow/x86_64/gcc_4.6.2/lapack_3.4.0/liblapack.a

export PYTHONPATH=${HOME}/packages/lib/python2.7/site-packages/

unset LDFLAGS
unset CFLAGS
python setup.py build
python setup.py install --prefix= ${HOME}/packages

popd
