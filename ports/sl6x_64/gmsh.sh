#GMSH
export GMSH_DIR=$workdir/gmsh
mkdir $GMSH_DIR
cd $workdir
wget http://geuz.org/gmsh/src/gmsh-2.6.1-source.tgz
tar xzf gmsh-2.6.1-source.tgz
cd gmsh-2.6.1-source
mkdir BuildDir
cd BuildDir
cmake -DCMAKE_INSTALL_PREFIX=$GMSH_DIR -DCMAKE_BUILD_TYPE=release -DCMAKE_PREFIX_PATH=$workdir/openmpi/include ..
make lib
make shared
make -j8 install
export PATH=$PATH:$GMSH_DIR/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GMSH_DIR/lib
