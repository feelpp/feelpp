#!/bin/sh

export FEEL_DIR=$HOME/sources/feel/

# go to source directory
if [! -d $HOME/sources/ ]; then mkdir -p $HOME/sources/; fi
cd $HOME/sources/

#this is due to a bug, the process extracts information from the svn output which needs to be in english
export LC_MESSAGES=en_GB

#do the actual work
if [ -f $HOME/sources/$1.log]; then  rm -f $HOME/sources/$1.log; fi
/usr/bin/ctest -S $FEEL_DIR/cmake/dashboard/ctest.cmake,$1 -V > $HOME/sources/$1.log 2>&1
