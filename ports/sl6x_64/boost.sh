#Installation de BOOST
mkdir $workdir/_boost
cd $workdir/_boost
wget http://ignum.dl.sourceforge.net/project/boost/boost/1.49.0/boost_1_49_0.tar.bz2
tar xjf boost_1_49_0.tar.bz2
cd boost_1_49_0
rm user-config.jam
echo "using mpi ;" >> user-config.jam
echo "" >> user-config.jam
./bootstrap.sh
./bjam -j4 install \
      --layout=tagged \
      --prefix=$workdir/boost\
      --user-config=user-config.jam \
      variant=release \
      threading=single,multi \
      link=static,shared
export Boost_DIR=$workdir/boost
