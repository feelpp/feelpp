if test $# -eq 0 
then
    echo "usage: $0 <toolset> [<target]>"
else
    if test "$BOOST_ROOT" == ""
    then
        set  +A BOOST_ROOT `dirname $PWD`
        set  +A BOOST_ROOT `dirname $BOOST_ROOT`
        set  +A BOOST_ROOT `dirname $BOOST_ROOT`
    fi
    if test "$2" != ""
    then
        export ALL_LOCATE_TARGET=$2
    fi
    if test "$ALL_LOCATE_TARGET" == ""
    then
        export ALL_LOCATE_TARGET=$BOOST_ROOT
    fi
    export BOOST_ROOT
    echo Running tests for $1 on $BOOST_ROOT to $ALL_LOCATE_TARGET
    bjam --dump-test -sTOOLS=$1 test >bjam.log 2>&1
    process_jam_log <bjam.log $ALL_LOCATE_TARGET
    compiler_status2 --locate-root $ALL_LOCATE_TARGET $BOOST_ROOT compiler_status.html links.html
fi
