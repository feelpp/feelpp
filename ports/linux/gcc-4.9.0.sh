apt-get install libgmp-dev libmpc-dev libmpfr-dev gnat-4.8 libgtk2.0-dev build-essential gawk m4 gcc-multilib libart-2.0-dev libxtst-dev
 ../configure --prefix=/data/software/install/gcc-4.9.0 \
  --enable-bootstrap \
  --enable-shared \
  --enable-threads=posix \
  --enable-checking=release \
  --with-system-zlib \
  --enable-__cxa_atexit \
  --disable-libunwind-exceptions \
  --enable-gnu-unique-object \
  --enable-languages=c,c++,objc,obj-c++,fortran,ada \
  --enable-java-awt=gtk \
  --disable-dssi \
  --enable-libgcj-multifile \
  --with-ppl \
  --with-cloog \
  --with-tune=generic \
  --with-arch_32=i686
make -j
make install
