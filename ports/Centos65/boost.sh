#Installation de BOOST
mkdir -p $boostDir/src
cd $boostDir/src
_ver=54
if [ $_ver -eq 54 ]
  then wget -c http://netcologne.dl.sourceforge.net/project/boost/boost/1.54.0/boost_1_54_0.tar.bz2
  else wget -c http://kent.dl.sourceforge.net/project/boost/boost/1.55.0/boost_1_55_0.tar.bz2
fi

tar xjf boost_1_$ver_0.tar.bz2
cd boost_1_$ver_0
rm user-config.jam
echo "using mpi ;" >> user-config.jam
echo "" >> user-config.jam
./bootstrap.sh
./bjam install \
      --layout=tagged \
      --prefix=$boostDir\
      --user-config=user-config.jam \
      variant=release \
      threading=single,multi \
      link=static,shared
