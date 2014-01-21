# On clone le dépot GIT
cd $baseDir
git clone https://github.com/feelpp/feelpp.git 
mkdir feelpp_build_dir 
cd feelpp_build_dir
cmake ../feelpp -DFEELPP_ENABLE_DOCUMENTATION=off \
make -j$nbProc
