#!/bin/bash
if test $# -eq 0 
then
    echo "usage: $0 <toolset>" "[<boost root directory>]" "[<target directory>]"
else
    echo "Running tests for $1"
    if [ ! "$2" = "" ]
    then
        export BOOST_ROOT=$2
        if [ ! "$3" == "" ]
        then
        	export ALL_LOCATE_TARGET=$3
        fi
    else
        if [ "$BOOST_ROOT" = "" ]
        then
        	set BOOST_ROOT `dirname $PWD`
        	set BOOST_ROOT `dirname $BOOST_ROOT`
        	set BOOST_ROOT `dirname $BOOST_ROOT`
        	export BOOST_ROOT
        fi
    fi

    if [ "$ALL_LOCATE_TARGET" == "" ]
    then
        export ALL_LOCATE_TARGET=$BOOST_ROOT
    fi
    bjam --dump-test -sTOOLS=$1 test >bjam.log 2>&1
    process_jam_log <bjam.log $ALL_LOCATE_TARGET
    compiler_status2 --locate-root $ALL_LOCATE_TARGET $BOOST_ROOT compiler_status.html links.html
fi
