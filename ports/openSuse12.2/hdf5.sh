# Script for building hdf5 library
# We need to build the hdf5 library and install h5cc which is used by cmake to detect the installation parameters
# otherwise Feel++ will automatically revert on system installed hdf5 library

# Does not work with cmake
# Following command line was tested
# cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local/feelpp/install/with-gcc-4.8.1/with-openmpi-1.8.3/hdf5-1.8.13 -DBUILD_SHARED_LIBS=ON -DHDF5_BUILD_TOOLS=ON -DHDF5_ENABLE_PARALLEL=ON -DHDF5_ENABLE_Z_LIB_SUPPORT=ON -DZLIB_DIR=/usr/lib64 -DCMAKE_BUILD_TYPE=Release

# Example with hdf5-1.8.13.tar.gz
tar zxvf hdf5-1.8.13.tar.gz
cd hdf5-1.8.13/

# SET MPICC_PREFIX and INSTALL_PREFIX to the correct values
CC=${MPICC_PREFIX}/mpicc ./configure --enable-parallel --prefix=${INSTALL_PREFIX}/hdf5-1.8.13/ --enable-build-all --enable-production

# build and install the lib
make -j ${NJOBS}
make install

# Important step:
# cmake uses h5cc and h5pcc to detect the hdf5 library and include directories. 
# Only h5pcc is installed by default, we just symlink it to h5cc to ensure that cmake doesn't try to use the
# system installation of h5cc
ln -s ${INSTALL_PREFIX}/hdf5-1.8.13/bin/h5pcc ${INSTALL_PREFIX}/hdf5-1.8.13/bin/h5cc
