SET(MODEL Nightly)
IF(${CTEST_SCRIPT_ARG} MATCHES Experimental)
  SET(MODEL Experimental)
ENDIF()
IF(${CTEST_SCRIPT_ARG} MATCHES Continuous)
  SET(MODEL Continuous)
ENDIF()

SET (CTEST_SOURCE_DIRECTORY "/home/prudhomm/sources/life")
set(CTEST_BINARY_DIRECTORY  "${CTEST_SOURCE_DIRECTORY}/sources/build-${CTEST_BUILD_NAME}")
set (CTEST_COMMAND "ctest -D ${MODEL}" )
SET (CTEST_CMAKE_COMMAND "cmake" )

SET (CTEST_SVN_COMMAND    "svn" )
SET (CTEST_SVN_CHECKOUT   "${CTEST_SVN_COMMAND} co -q svn://scm.ljkforge.imag.fr/svn/life/life/trunk ${CTEST_SOURCE_DIRECTORY}")
set (CTEST_UPDATE_COMMAND "${CTEST_SVN_COMMAND}")



if (${MODEL} MATCHES Nightly )

  # should ctest wipe the binary tree before running
  SET (CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE)

endif()

if ( ${MODEL} MATCHES Continuous )
  while (${CTEST_ELAPSED_TIME} LESS 36000)
    set (START_TIME ${CTEST_ELAPSED_TIME})
    ctest_start (Continuous)

    SET (CTEST_START_WITH_EMPTY_BINARY_DIRECTORY_ONCE 1)
    
    ctest_sleep( ${START_TIME} 300 ${CTEST_ELAPSED_TIME})
  endwhile()
endif()
