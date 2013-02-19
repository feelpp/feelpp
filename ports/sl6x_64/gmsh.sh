#GMSH
mkdir $workdir/_gmsh
cd $workdir/_gmsh
wget -c http://geuz.org/gmsh/src/gmsh-2.6.1-source.tgz
if [ ! -d "gmsh-2.6.1-source" ]; then
tar xzf gmsh-2.6.1-source.tgz
fi
cd gmsh-2.6.1-source
mkdir BuildDir
cd BuildDir
cmake28 -DCMAKE_INSTALL_PREFIX=$GMSH_DIR -DCMAKE_BUILD_TYPE=release -DCMAKE_PREFIX_PATH=$openmpiDir/include ..
make -j$nbProc lib
make -j$nbProc shared
make -j$nbProc install
