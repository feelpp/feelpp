# On clone le dépot GIT
cd $basedir
if [ -d "$basedir/feelpp" ]; then
cd $basedir/feelpp; git pull
else;
git clone https://github.com/feelpp/feelpp.git 
fi
mkdir -p $basedir/feelpp_build_dir 
cd $basedir/feelpp_build_dir
cmake $basedir/feelpp -DFEELPP_ENABLE_DOCUMENTATION=off 
make -j4
