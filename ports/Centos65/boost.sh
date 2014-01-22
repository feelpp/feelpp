#Installation de BOOST
mkdir -p $boostDir/src
cd $boostDir/src
wget -c http://kent.dl.sourceforge.net/project/boost/boost/1.55.0/boost_1_55_0.tar.bz2
tar xjf boost_1_55_0.tar.bz2
cd boost_1_55_0
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
