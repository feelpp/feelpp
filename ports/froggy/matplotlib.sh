#! /bin/bash

. ~/.bash_profile

# if no loaded
#module load python/2.7.3_gcc-4.6.2
module load libpng/1.5.9_gcc-4.6.2
module load freetype/2.4.9_gcc-4.6.2
module load gtk+/2.24.10_gcc-4.6.2

if [ ! -d "${HOME}/packages/build" ]; then
  mkdir -p ${HOME}/packages/build
fi
pushd ${HOME}/packages/matplotlib-1.2.0

export PYTHONPATH=${HOME}/packages/lib/python2.7/site-packages/

python setup.py build
python setup.py install --prefix=${HOME}/packages

popd
