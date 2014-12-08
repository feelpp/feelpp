#GMSH
mkdir $workdir/_gmsh
cd $workdir/_gmsh

wget -c http://geuz.org/gmsh/src/gmsh-svn-source.tgz
tar xzf gmsh-svn-source.tgz

cd `tar -tvf gmsh-svn-source.tgz | head -n 1 | sed "s/\// /g" | awk '{print $7}'`

mkdir BuildDir
cd BuildDir
cmake -DCMAKE_INSTALL_PREFIX=$GMSH_DIR -DENABLE_MPI=ON -DENABLE_BUILD_LIB=ON -DENABLE_BUILD_SHARED=ON -DCMAKE_BUILD_TYPE=release -DCMAKE_PREFIX_PATH=$openmpiDir/include ..
make -j$nbProc lib
make -j$nbProc shared
make -j$nbProc install
