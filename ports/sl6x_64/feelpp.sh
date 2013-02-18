# On clone le dépot GIT
mkdir $workdir/_feel
cd $workdir/_feel
git clone https://github.com/feelpp/feelpp.git
mkdir feelpp_build_dir 
cd feelpp_build_dir
cmake ../feelpp \
      -DCMAKE_BUILD_TYPE=release \
      -DFEELPP_ENABLE_DOCUMENTATION=off \
      -DCMAKE_INSTALL_PREFIX=$feelppDir
make -j4 install
