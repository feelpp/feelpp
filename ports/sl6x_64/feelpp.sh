# On clone le dépot GIT
#mkdir $workdir/_feel
cd $workdir/_feel
#git clone https://github.com/feelpp/feelpp.git
#mkdir feelpp_build_dir 
cd feelpp_build_dir
cmake28 ../feelpp \
      -DBoost_NO_BOOST_CMAKE=TRUE \
      -DCMAKE_BUILD_TYPE=release \
      -DCMAKE_INSTALL_PREFIX=$feelppDir
make -j$nbProc install
