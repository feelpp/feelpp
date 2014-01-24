# On clone le dépot GIT
cd $basedir
git clone https://github.com/feelpp/feelpp.git 
mkdir -p $basedir/feelpp_build_dir 
cd $basedir/feelpp_build_dir
cmake $basedir/feelpp -DFEELPP_ENABLE_DOCUMENTATION=off 
make -j4
