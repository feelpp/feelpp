if test $# -eq 0
then
    echo "usage: $0 <test header file> <toolset> [<target>]"
else
    export BOOST_ARCHIVE_LIST=$1
    runtest.ksh $2 $3
fi
 
