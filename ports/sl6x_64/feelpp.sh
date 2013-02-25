# On clone le dépot GIT
mkdir $workdir/_feel
mkdir $workdir/_feel/gcc_$gccVersion
cd $workdir/_feel/gcc_$gccVersion
git clone https://github.com/feelpp/feelpp.git
mkdir feelpp_build_dir 
cd feelpp_build_dir
cmake28 ../feelpp \
      -DBoost_NO_BOOST_CMAKE=TRUE \
      -DFEELPP_ENABLE_DOCUMENTATION=off \
      -DCMAKE_BUILD_TYPE=release \
      -DCMAKE_INSTALL_PREFIX=$feelppDir
make -j$nbProc install
