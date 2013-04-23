#! /bin/bash

#profile should include feelpprc from V.C
. ~/.bash_profile

if [ ! -d "${HOME}/packages/build" ]; then
  mkdir -p ${HOME}/packages/build
fi
pushd ${HOME}/packages/build

export cflags="-g -O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security"
export cxxflags="-g -O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security"
export fflags="-g -O2"
export ldflags="-Wl,-z,relro -Wl,--as-needed"

SRCDIR=${HOME}/packages/bison-2.5.dfsg

touch --date="Jan 01 2000" \
	$SRCDIR/doc/bison.info \
        $SRCDIR/doc/bison.texinfo \
        $SRCDIR/doc/fdl.texi

$SRCDIR/configure \
	--prefix=${HOME}/packages/bison

make || exit
make install  || exit

popd
rm -rf ${HOME}/packages/build

#make all lib shared
