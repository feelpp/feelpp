#Installation de BOOST
mkdir $workdir/_boost
cd $workdir/_boost
wget -c http://heanet.dl.sourceforge.net/project/boost/boost/1.55.0/boost_1_55_0.tar.gz
tar xzf boost_1_55_0.tar.gz
cd boost_1_55_0
rm user-config.jam
echo "using mpi ;" >> user-config.jam
echo "" >> user-config.jam
./bootstrap.sh
./bjam -j4 install \
      --layout=tagged \
      --prefix=$boostDir\
      --user-config=user-config.jam \
      variant=release \
      threading=single,multi \
      link=static,shared
