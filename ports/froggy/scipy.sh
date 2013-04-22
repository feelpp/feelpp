#! /bin/bash

. ~/.bash_profile

# if no loaded
module load python/2.7.3_gcc-4.6.2
module load lapack/3.4.0_gcc-4.6.2
module load swig/2.0.9_gcc-4.6.2

pushd ${HOME}/packages/scipy-0.12.0

export BLAS=/home/chabannes/packages/blas/BLAS/blas_LINUX.a
export LAPACK=/applis/ciment/v2/stow/x86_64/gcc_4.6.2/lapack_3.4.0/liblapack.a

export PYTHONPATH=${HOME}/python/lib/python2.7/site-packages/

python setup.py config_fc --noarch build
python setup.py install --prefix= ${HOME}/packages

popd
