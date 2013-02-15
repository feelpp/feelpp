# On clone le dépot GIT
cd $workdir
git clone https://github.com/feelpp/feelpp.git
mkdir feelppLib
mkdir feelpp_build_dir 
cd feelpp_build_dir
cmake $workdir/feelpp \
      -DCMAKE_BUILD_TYPE=release \
      -DCMAKE_INSTALL_PREFIX=$workdir/feelppLib
make -j4 install
