#GMSH
mkdir -p $gmshDir/src
cd $gmshDir/src
wget -c http://www.geuz.org/gmsh/bin/Linux/gmsh-svn-Linux64.tgz
tar xzf gmsh-svn-Linux64.tgz
cd gmsh-svn-Linux64.tgz
mkdir BuildDir
cd BuildDir
cmake -DCMAKE_INSTALL_PREFIX=$GMSH_DIR -DCMAKE_BUILD_TYPE=release -DENABLE_BUILD_LIB=1 -DENABLE_BUILD_SHARED=1
make -j$nbProc lib
make -j$nbProc shared
make -j$nbProc install
