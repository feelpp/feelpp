#GMSH
mkdir -p $gmshDir/src
cd $gmshDir/src
pwd
wget -c http://www.geuz.org/gmsh/src/gmsh-svn-source.tgz
tar xzf gmsh-svn-source.tgz
cd `tar -tvf gmsh-svn-source.tgz | head -n 1 | sed "s/\// /g" | awk '{print $7}'`
mkdir BuildDir
cd BuildDir
cmake -DCMAKE_INSTALL_PREFIX=$GMSH_DIR -DCMAKE_BUILD_TYPE=release -DENABLE_BUILD_LIB=1 -DENABLE_BUILD_SHARED=1 ..
make -j4 lib
make -j4 shared
make -j4 install
