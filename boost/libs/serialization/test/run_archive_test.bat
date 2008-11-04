@ echo off
if "%1" == "" goto message
    set BOOST_ARCHIVE_LIST=%1
    runtest.bat %2 %3 %4
    goto end
:message
    echo usage: run_archive_test "<test header file>" "<toolset>" "[<boost root directory>]" "[<target directory>]"
:end
