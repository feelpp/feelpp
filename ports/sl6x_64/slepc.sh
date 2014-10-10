# download slepc package and untar it
cd $SLEPC_PACKAGE_DIR

export SLEPC_DIR=`pwd`
./configure --prefix=$SLEPC_INSTALL_PATH

# follow the instructions for building and installing using the provided command
