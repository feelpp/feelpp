@ echo off
if "%1" == "" goto message
    if "%2" == "" goto 1
        set BOOST_ROOT=%2
    if "%3" == "" goto 1
        set ALL_LOCATE_TARGET=%3
        :1
    if "%BOOST_ROOT%" == "" goto message
    if not "%ALL_LOCATE_TARGET%" == "" goto 2
        set ALL_LOCATE_TARGET=%BOOST_ROOT%
    :2
    echo Running tests for %1 on %BOOST_ROOT% to %ALL_LOCATE_TARGET%
    bjam --dump-test -sTOOLS=%1 test >bjam.log 2>&1
    process_jam_log <bjam.log %ALL_LOCATE_TARGET%
    compiler_status2 --locate-root %ALL_LOCATE_TARGET% %BOOST_ROOT%  compiler_status.html links.html
    goto end
:message
    echo usage: runtest "<toolset>" "<boost root directory>" "[<target directory>]"
:end
