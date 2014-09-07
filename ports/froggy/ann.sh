#! /bin/bash -x

# Optionnaly set version and nbProc
: ${Version:="1.1.2"}
: ${gccVersion:="4.8.2"}
: ${nbProc:=$(echo `/usr/bin/getconf _NPROCESSORS_ONLN`)}
echo "Installation de ann-$Version"

 
tar zxf ~/Downloads/ann_${Version}+doc.orig.tar.gz 

cd ann-${Version}+doc/

# Eventually apply patches
for file in $(ls ~/Patches/ann-${Version}+doc/*); do
   echo "Applying " $file
   patch -p1 < $file
done

touch NEWS README AUTHORS ChangeLog
#autoupdate --force
#autoreconf --force

rm -f aclocal.m4

libtoolize --automake --force
aclocal --force  
autoheader --force  
autoconf --force  
automake  --force --add-missing

#mkdir -p build
#cd build
find . -name Makefile | xargs  rm

./configure --prefix=/home/trophime/modules/libs/ann-${Version}+doc
make -j ${nbProc}
make -j ${nbProc} install

 
