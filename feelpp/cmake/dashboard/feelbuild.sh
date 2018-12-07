#!/bin/sh

export PATH=/usr/bin:/opt/local/bin:$PATH
export FEEL_DIR=$HOME/sources/feel/

ARCH=` uname`
if [ "$ARCH" = "Darwin" ]; then
    export PETSC_DIR=/opt/local/lib/petsc
    export PETSC_ARCH=darwin
fi

# go to source directory
if ! [ -d $HOME/sources/ ]; then mkdir -p $HOME/sources/; fi
cd $HOME/sources/
if ! [ -d $HOME/sources/feel ]; then
    svn checkout svn://scm.forge.imag.fr/var/lib/gforge/chroot/scmrepos/svn/life/trunk/life/trunk feel

    echo "type 'crontab -e' and enter"
    echo " 0  1   *   *   * $HOME/sources/feel/cmake/dashboard/feelbuild.sh Nightly"
    exit
fi

#this is due to a bug, the process extracts information from the svn output which needs to be in english
unset LC_MESSAGES
unset LANG
unset LC_ALL

#do the actual work
if [ -f $HOME/sources/$1.log ]; then  rm -f $HOME/sources/$1.log; fi
if [ -x /usr/bin/ctest ]; then
    /usr/bin/ctest -S $FEEL_DIR/cmake/dashboard/ctest.cmake,$1 -V > $HOME/sources/$1.log 2>&1
elif [ -x /opt/local/bin/ctest ]; then
    /opt/local/bin/ctest -S $FEEL_DIR/cmake/dashboard/ctest.cmake,$1 -V > $HOME/sources/$1.log 2>&1
else
    echo "CTest not found."
fi

