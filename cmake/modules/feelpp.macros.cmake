# - Find Feel
if (POLICY CMP0045)
    # error on non-existent target in get_target_property
    cmake_policy(SET CMP0045 OLD)
endif()
INCLUDE(feelpp.precompiled.headers)
INCLUDE(ParseArguments)

# define CMAKE_INSTALL_DOCDIR
include(GNUInstallDirs)


# list the subdicrectories of directory 'curdir'
macro(feelpp_list_subdirs result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
endmacro(feelpp_list_subdirs)

macro(feelpp_add_testcase )
  PARSE_ARGUMENTS(FEELPP_CASE
    "NAME;PREFIX;DEPS;CATEGORY"
    ""
    ${ARGN}
    )
  CAR(FEELPP_CASE_NAME ${FEELPP_CASE_DEFAULT_ARGS})
  if ( FEELPP_CASE_PREFIX )
    set( target ${FEELPP_CASE_PREFIX}_add_testcase_${FEELPP_CASE_NAME})
  else()
    set( target feelpp_add_testcase_${FEELPP_CASE_NAME})
  endif()
  add_custom_target(${target})
  if ( FEELPP_CASE_DEPS )
    foreach(case ${FEELPP_CASE_DEPS})
      add_dependencies(${target} ${FEELPP_CASE_PREFIX}_add_testcase_${case})
    endforeach()
  endif()
  ADD_CUSTOM_COMMAND(
    TARGET ${target}
    POST_BUILD
    COMMAND rsync
    ARGS -aLv
    ${CMAKE_CURRENT_SOURCE_DIR}/${FEELPP_CASE_NAME}
    ${CMAKE_CURRENT_BINARY_DIR}/
    COMMENT "Syncing testcase ${testcase} in ${CMAKE_CURRENT_BINARY_DIR} from ${CMAKE_CURRENT_SOURCE_DIR}/${FEELPP_CASE_NAME}"
    )
  #execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${testcase} ${CMAKE_CURRENT_BINARY_DIR} )
  #file(COPY ${testcase} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  #get_filename_component( relpath ${CMAKE_CURRENT_SOURCE_DIR}/${FEELPP_CASE_NAME} DIRECTORY BASE_DIR ${FEELPP_SOURCE_DIR})
  #message(STATUS "testcase ${FEELPP_CASE_NAME} -> ${relpath} ")
  if ( FEELPP_CASE_CATEGORY )
    #if (NOT EXISTS ${CMAKE_INSTALL_PREFIX}/share/feel/testcases/${FEELPP_CASE_CATEGORY})
    #  execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/share/feel/)
    #  execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/share/feel/testcases/)
    #  execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/share/feel/testcases/${FEELPP_CASE_CATEGORY})
    #endif()

    #INSTALL(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${FEELPP_CASE_NAME} DESTINATION share/feel/applications/${CATEGORY} COMPONENT install-testcase)
    # ADD_CUSTOM_COMMAND(
    # TARGET ${target}
    # POST_BUILD
    # COMMAND rsync
    # ARGS -av
    # ${CMAKE_CURRENT_SOURCE_DIR}/${FEELPP_CASE_NAME}
    # ${CMAKE_INSTALL_PREFIX}/share/feel/testcases/${FEELPP_CASE_CATEGORY}
    # COMMENT "Syncing testcase ${testcase} in ${CMAKE_INSTALL_PREFIX}/share/feel/testcases/${FEELPP_CASE_CATEGORY} from ${CMAKE_CURRENT_SOURCE_DIR}/${FEELPP_CASE_NAME}")
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${FEELPP_CASE_NAME}
      DESTINATION share/feelpp/testcases/${FEELPP_CASE_CATEGORY} COMPONENT testcases)
    if ( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/README.adoc )
      install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/README.adoc 
        DESTINATION share/feelpp/testcases/${FEELPP_CASE_CATEGORY} COMPONENT testcases)
    endif()
    #add_dependencies(install-testcase ${target})
  endif()
endmacro(feelpp_add_testcase)


# add a new application
macro(feelpp_add_application)

  PARSE_ARGUMENTS(FEELPP_APP
    "SRCS;LINK_LIBRARIES;CFG;GEO;MESH;LABELS;DEFS;DEPS;SCRIPTS;TEST;TIMEOUT;PROJECT;EXEC;MAN"
    "NO_TEST;NO_MPI_TEST;NO_SEQ_TEST;EXCLUDE_FROM_ALL;INCLUDE_IN_ALL;ADD_OT;NO_FEELPP_LIBRARY;INSTALL"
    ${ARGN}
    )
  CAR(FEELPP_APP_NAME ${FEELPP_APP_DEFAULT_ARGS})

  if ( FEELPP_APP_PROJECT )
    set(execname feelpp_${FEELPP_APP_PROJECT}_${FEELPP_APP_NAME})
  else( FEELPP_APP_PROJECT )
    if ( PROJECT_NAME AND
        ( NOT PROJECT_NAME STREQUAL "Feel++" )
        )

      if ( PROJECT_SHORTNAME )
        #message(STATUS "project: ${PROJECT_NAME} shortname: ${PROJECT_SHORTNAME}")
        set(execname feelpp_${PROJECT_SHORTNAME}_${FEELPP_APP_NAME})
      else()
        #message(STATUS "project: ${PROJECT_NAME} ")
        set(execname feelpp_${PROJECT_NAME}_${FEELPP_APP_NAME})
      endif()
    else()
      set(execname feelpp_${FEELPP_APP_NAME})
    endif()
  endif( FEELPP_APP_PROJECT )
  if  (FEELPP_APP_EXEC )
    set( ${FEELPP_APP_EXEC} ${execname} )
  endif()

  if ( FEELPP_ENABLE_VERBOSE_CMAKE )
    MESSAGE("*** Arguments for Feel++ application ${FEELPP_APP_NAME}")
    MESSAGE("    Sources: ${FEELPP_APP_SRCS}")
    MESSAGE("    Link libraries: ${FEELPP_APP_LINK_LIBRARIES}")
    MESSAGE("       Cfg file: ${FEELPP_APP_CFG}")
    MESSAGE("      Deps file: ${FEELPP_APP_DEPS}")
    MESSAGE("      Defs file: ${FEELPP_APP_DEFS}")
    MESSAGE("       Geo file: ${FEELPP_APP_GEO}")
    MESSAGE("       Mesh file: ${FEELPP_APP_MESH}")
    MESSAGE("    test timeout: ${FEELPP_APP_TIMEOUT}")
    MESSAGE("       Exec file: ${execname}")
    MESSAGE("exclude from all: ${FEELPP_APP_EXCLUDE_FROM_ALL}")
    MESSAGE("include from all: ${FEELPP_APP_INCLUDE_IN_ALL}")
  endif()

  if ( FEELPP_APP_EXCLUDE_FROM_ALL)
    add_executable(${execname}  EXCLUDE_FROM_ALL  ${FEELPP_APP_SRCS}  )
  elseif( FEELPP_APP_INCLUDE_IN_ALL)
    add_executable(${execname}  ${FEELPP_APP_SRCS}  )
  else()
    add_executable(${execname}  ${FEELPP_APP_SRCS}  )
  endif()
  if ( FEELPP_APP_DEPS )
    add_dependencies(${execname} ${FEELPP_APP_DEPS})
  endif()
  if ( FEELPP_APP_DEFS )
    set_property(TARGET ${execname} PROPERTY COMPILE_DEFINITIONS ${FEELPP_APP_DEFS})
  endif()
  if ( FEELPP_APP_NO_FEELPP_LIBRARY )
      target_link_libraries( ${execname} ${FEELPP_APP_LINK_LIBRARIES} ${FEELPP_LIBRARIES})
  else()
      target_link_libraries( ${execname} ${FEELPP_LIBRARY} ${FEELPP_APP_LINK_LIBRARIES} ${FEELPP_LIBRARIES})
  endif()

  # Use feel++ lib precompiled headers.
  #if( FEELPP_ENABLE_PCH )
  #    add_precompiled_header( ${FEELPP_LIBRARY} )
  #endif()
  # Create application precompiled headers.
  if( FEELPP_ENABLE_PCH_APPLICATIONS )
    add_precompiled_header( ${execname} )
  endif()

  # install rule if INSTALL if target is marked to be installed
  if ( FEELPP_APP_INSTALL )
    install(TARGETS ${execname} RUNTIME DESTINATION bin COMPONENT Bin)
  endif()

  if ( NOT FEELPP_APP_NO_TEST )
    IF(NOT FEELPP_APP_NO_MPI_TEST AND NProcs2 GREATER 1)
      add_test(NAME ${execname}-np-${NProcs2} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NProcs2} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/${execname} ${FEELPP_APP_TEST} ${MPIEXEC_POSTFLAGS} )
    ENDIF()
    IF(NOT FEELPP_APP_NO_SEQ_TEST)
      add_test(NAME ${execname}-np-1 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/${execname} ${FEELPP_APP_TEST} ${MPIEXEC_POSTFLAGS})
    endif()
  endif()

  #add_dependencies(crb ${execname})
  # add TIMEOUT to test
  if ( FEELPP_APP_TIMEOUT )
    if ( NOT FEELPP_APP_NO_TEST )
      IF(NOT FEELPP_APP_NO_MPI_TEST AND NProcs2 GREATER 1)
        set_property(TEST ${execname}-np-${NProcs2}  PROPERTY TIMEOUT ${FEELPP_APP_TIMEOUT})
      endif()
      IF(NOT FEELPP_APP_NO_SEQ_TEST)
        set_property(TEST ${execname}-np-1  PROPERTY TIMEOUT ${FEELPP_APP_TIMEOUT})
      ENDIF()
    endif()
  else()
    if ( NOT FEELPP_APP_NO_TEST )
      IF(NOT FEELPP_APP_NO_MPI_TEST AND NProcs2 GREATER 1)
        set_property(TEST ${execname}-np-${NProcs2}  PROPERTY TIMEOUT  ${FEELPP_DEFAULT_TEST_TIMEOUT})
      endif()
      IF(NOT FEELPP_APP_NO_SEQ_TEST)
        set_property(TEST ${execname}-np-1  PROPERTY TIMEOUT  ${FEELPP_DEFAULT_TEST_TIMEOUT})
      endif()
    endif()
  endif()
  # Add label if provided
  if ( FEELPP_APP_LABELS )
    set_property(TARGET ${execname} PROPERTY LABELS ${FEELPP_APP_LABELS})
    if ( NOT FEELPP_APP_NO_TEST )
      IF(NOT FEELPP_APP_NO_MPI_TEST AND NProcs2 GREATER 1)
        set_property(TEST ${execname}-np-${NProcs2} PROPERTY LABELS ${FEELPP_APP_LABELS})
      ENDIF()
      IF(NOT FEELPP_APP_NO_SEQ_TEST)
        set_property(TEST ${execname}-np-1 PROPERTY LABELS ${FEELPP_APP_LABELS})
      endif()
    endif()
    foreach(l ${FEELPP_APP_LABELS})
      if ( TARGET ${l} )
        add_dependencies( ${l} ${execname} )
      endif()
    endforeach(l)
  endif()

  # add manual page
  if ( FEELPP_APP_MAN )
    feelpp_add_man( ${execname} ${FEELPP_APP_MAN} 1 )
  endif( FEELPP_APP_MAN )
  
  # include schedulers
  include( feelpp.schedulers )


  if ( FEELPP_APP_CFG )
    foreach(  cfg ${FEELPP_APP_CFG} )
      #      if ( EXISTS ${cfg} )
      # extract cfg filename  to be copied in binary dir
      #get_filename_component( CFG_NAME ${cfg} NAME )
      #configure_file( ${cfg} ${CFG_NAME} )

      file(COPY ${cfg} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

      #INSTALL(FILES "${cfg}"  DESTINATION share/feel/config)
      #      else()
      #        message(WARNING "Executable ${FEELPP_APP_NAME}: configuration file ${cfg} does not exist")
      #      endif()
    endforeach()
  endif(FEELPP_APP_CFG)

  if ( FEELPP_APP_GEO )
    foreach(  geo ${FEELPP_APP_GEO} )
      # extract geo filename  to be copied in binary dir
      #get_filename_component( GEO_NAME ${geo} NAME )
      #configure_file( ${geo} ${GEO_NAME} )
      file(COPY ${geo} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
      INSTALL(FILES "${geo}"  DESTINATION share/feel/geo)
    endforeach()
  endif(FEELPP_APP_GEO)

  if ( FEELPP_APP_MESH )
    foreach(  mesh ${FEELPP_APP_MESH} )
      # extract mesh filename  to be copied in binary dir
      #get_filename_component( MESH_NAME ${mesh} NAME )
      #configure_file( ${mesh} ${MESH_NAME} )
      file(COPY ${mesh} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
      INSTALL(FILES "${mesh}"  DESTINATION share/feel/mesh)
    endforeach()
  endif(FEELPP_APP_MESH)

  if ( FEELPP_APP_SCRIPTS )
    foreach(  script ${FEELPP_APP_SCRIPTS} )
      # extract mesh filename  to be copied in binary dir
      get_filename_component( SCRIPT_NAME ${script} NAME )
      configure_file( ${script} ${SCRIPT_NAME} )
    endforeach()
  endif(FEELPP_APP_SCRIPTS)

  if ( FEELPP_APP_ADD_OT AND OPENTURNS_FOUND )
    set(pycpp "${execname}_pywrapper.cpp")
    set(xml "${execname}_ot.xml")
    set(FEELPP_APP_OT_WRAPPER_NAME "${execname}_ot")
    get_filename_component( FEELPP_APP_OUTPUT_WE ${FEELPP_APP_CFG} NAME_WE )
    set(FEELPP_APP_OUTPUT ${FEELPP_APP_OUTPUT_WE}.res )

    configure_file(${FEELPP_SOURCE_DIR}/cmake/templates/ot_python_command_wrapper.cpp ${pycpp})
    if ( NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${FEELPP_APP_OT_WRAPPER_NAME}.xml)
      configure_file(${FEELPP_SOURCE_DIR}/cmake/templates/ot_python.xml.in ${CMAKE_CURRENT_SOURCE_DIR}/${FEELPP_APP_OT_WRAPPER_NAME}.xml)
    endif()
    configure_file(${FEELPP_APP_OT_WRAPPER_NAME}.xml ${xml})
    feelpp_ot_add_python_module(${FEELPP_APP_OT_WRAPPER_NAME} ${pycpp}
      LINK_LIBRARIES ${OpenTURNS_LIBRARIES}
      CFG ${FEELPP_APP_CFG} XML ${xml} TEST)
  endif()
endmacro(feelpp_add_application)


macro(OVERWITE_IF_DIFFERENT thetarget filename var dummy)
  if ( NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${filename} )
    # be careful if file does not exist we use dummy to generate the cpp file which will
    # then be overwritten using the cmake -E copy_if_different command
    #configure_file(${dummy}  ${CMAKE_CURRENT_BINARY_DIR}/${filename})
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${filename} ${var} )
  endif()
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/copy_${filename} ${var})
  add_custom_command(TARGET ${thetarget} COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_BINARY_DIR}/copy_${filename} ${CMAKE_CURRENT_BINARY_DIR}/${filename}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endmacro()

macro(feelpp_add_test)
  PARSE_ARGUMENTS(FEELPP_TEST
    "SRCS;LINK_LIBRARIES;CFG;GEO;MESH;LABEL;DEFS;DEPS;TIMEOUT;CLI;PROJECT;EXEC"
    "NO_TEST;NO_MPI_TEST;EXCLUDE_FROM_ALL;NO_FEELPP_LIBRARY"
    ${ARGN}
    )

  CAR(FEELPP_TEST_NAME ${FEELPP_TEST_DEFAULT_ARGS})
  get_directory_property( FEELPP_TEST_LABEL_DIRECTORY LABEL )
  #set(targetname feelpp_test_${FEELPP_TEST_NAME})

  if ( NOT FEELPP_TEST_SRCS )
    set(filename test_${FEELPP_TEST_NAME}.cpp)
    if ( FEELPP_TEST_NO_FEELPP_LIBRARY )
        feelpp_add_application( ${FEELPP_TEST_NAME} SRCS ${filename} CFG  ${FEELPP_TEST_CFG} GEO ${FEELPP_TEST_GEO} MESH ${FEELPP_TEST_MESH}  DEFS ${FEELPP_TEST_DEFS} PROJECT ${FEELPP_TEST_PROJECT} EXEC targetname LINK_LIBRARIES ${FEELPP_LIBRARIES} ${FEELPP_TEST_LINK_LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}  NO_TEST NO_FEELPP_LIBRARY )
    else()
        feelpp_add_application( ${FEELPP_TEST_NAME} SRCS ${filename} CFG  ${FEELPP_TEST_CFG} GEO ${FEELPP_TEST_GEO}  MESH ${FEELPP_TEST_MESH} DEFS ${FEELPP_TEST_DEFS}  PROJECT ${FEELPP_TEST_PROJECT} EXEC targetname LINK_LIBRARIES ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} ${FEELPP_TEST_LINK_LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}  NO_TEST )
    endif()
    #add_executable(${targetname} ${filename})
  else()
     if ( FEELPP_TEST_NO_FEELPP_LIBRARY )
       feelpp_add_application( ${FEELPP_TEST_NAME} SRCS ${FEELPP_TEST_SRCS}  CFG  ${FEELPP_TEST_CFG} GEO ${FEELPP_TEST_GEO}  MESH ${FEELPP_TEST_MESH} DEFS ${FEELPP_TEST_DEFS}   PROJECT ${FEELPP_TEST_PROJECT} TARGET targetname LINK_LIBRARIES ${FEELPP_LIBRARIES} ${FEELPP_TEST_LINK_LIBRARIES}  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}  NO_TEST NO_FEELPP_LIBRARY )
     else()
       feelpp_add_application( ${FEELPP_TEST_NAME} SRCS ${FEELPP_TEST_SRCS}  CFG  ${FEELPP_TEST_CFG} GEO ${FEELPP_TEST_GEO}  MESH ${FEELPP_TEST_MESH} DEFS ${FEELPP_TEST_DEFS}   PROJECT ${FEELPP_TEST_PROJECT} EXEC targetname LINK_LIBRARIES ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} ${FEELPP_TEST_LINK_LIBRARIES}  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}  NO_TEST )
     endif()
    #add_executable(${targetname} ${FEELPP_TEST_SRCS})
  endif()
  set( FEELPP_TEST_EXEC ${targetname} )
  #target_link_libraries(${targetname} ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} ${FEELPP_TEST_LINK_LIBRARIES}  )
  set_property(TARGET ${targetname} PROPERTY LABELS ${FEELPP_TEST_LABEL} ${FEELPP_TEST_LABEL_DIRECTORY})
  if ( TARGET  ${FEELPP_TEST_LABEL_DIRECTORY})
    add_dependencies(  ${FEELPP_TEST_LABEL_DIRECTORY} ${targetname} )
    add_dependencies( testsuite  ${FEELPP_TEST_LABEL_DIRECTORY} )
  elseif( TARGET testsuite )
    add_dependencies(testsuite ${targetname})
  endif()


  if ( NOT FEELPP_TEST_NO_TEST )
      # split command line options by whitespace into cmake list
      separate_arguments(FEELPP_TEST_CLI)
      set(BOOST_TEST_SEPARATOR "")
      if ( Boost_MAJOR_VERSION GREATER 0 AND Boost_MINOR_VERSION GREATER 59 ) #AND (FEELPP_TEST_CLI OR FEELPP_TEST_CFG) )
        set(BOOST_TEST_SEPARATOR "--")
      endif()
      unset( FEELPP_TEST_CFG_CLI )
      if ( FEELPP_TEST_CFG )
        set( FEELPP_TEST_CFG_CLI --config-file=${CMAKE_CURRENT_BINARY_DIR}/${FEELPP_TEST_CFG} )
      endif()
      IF(NOT FEELPP_TEST_NO_MPI_TEST AND NProcs2 GREATER 1)
        add_test(NAME feelpp_test_${FEELPP_TEST_NAME}-np-${NProcs2} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NProcs2} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/${targetname} --log_level=message ${BOOST_TEST_SEPARATOR} ${MPIEXEC_POSTFLAGS} ${FEELPP_TEST_CFG_CLI} ${FEELPP_TEST_CLI} --directory=testsuite/test_${FEELPP_TEST_NAME} --rm )
        set_property(TEST feelpp_test_${FEELPP_TEST_NAME}-np-${NProcs2}  PROPERTY LABELS ${FEELPP_TEST_LABEL}  ${FEELPP_TEST_LABEL_DIRECTORY} )
        if(CMAKE_BUILD_TYPE MATCHES Debug)
          set_tests_properties(feelpp_test_${FEELPP_TEST_NAME}-np-${NProcs2} PROPERTIES ENVIRONMENT "LSAN_OPTIONS=suppressions=${PROJECT_SOURCE_DIR}/suppressions.txt")
        endif()
      ENDIF()
      add_test(NAME feelpp_test_${FEELPP_TEST_NAME}-np-1 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/${targetname} --log_level=message ${BOOST_TEST_SEPARATOR} ${MPIEXEC_POSTFLAGS} ${FEELPP_TEST_CFG_CLI} ${FEELPP_TEST_CLI} --directory=testsuite/test_${FEELPP_TEST_NAME}  --rm )
      set_property(TEST feelpp_test_${FEELPP_TEST_NAME}-np-1  PROPERTY LABELS ${FEELPP_TEST_LABEL} ${FEELPP_TEST_LABEL_DIRECTORY} )
      if(CMAKE_BUILD_TYPE MATCHES Debug)
        set_tests_properties(feelpp_test_${FEELPP_TEST_NAME}-np-1 PROPERTIES ENVIRONMENT "LSAN_OPTIONS=suppressions=${PROJECT_SOURCE_DIR}/suppressions.txt")
        message(STATUS "option 1: feelpp_test_${FEELPP_TEST_NAME}-np-1")
      endif()
    endif()


    # add TIMEOUT to test
    if ( NOT FEELPP_TEST_NO_TEST )
      if ( FEELPP_TEST_TIMEOUT )
        IF(NOT FEELPP_TEST_NO_MPI_TEST AND NProcs2 GREATER 1)
          set_property(TEST feelpp_test_${FEELPP_TEST_NAME}-np-${NProcs2}  PROPERTY TIMEOUT ${FEELPP_TEST_TIMEOUT})
        endif()
        set_property(TEST feelpp_test_${FEELPP_TEST_NAME}-np-1  PROPERTY TIMEOUT ${FEELPP_TEST_TIMEOUT})
      else()
        IF(NOT FEELPP_TEST_NO_MPI_TEST AND NProcs2 GREATER 1)
          set_property(TEST feelpp_test_${FEELPP_TEST_NAME}-np-${NProcs2}  PROPERTY TIMEOUT ${FEELPP_DEFAULT_TEST_TIMEOUT})
        endif()
        set_property(TEST feelpp_test_${FEELPP_TEST_NAME}-np-1  PROPERTY TIMEOUT  ${FEELPP_DEFAULT_TEST_TIMEOUT})
      endif()
    endif()

    # cfg
    set(cfgname test_${FEELPP_TEST_NAME}.cfg)
    if ( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${cfgname} )
      configure_file(  ${cfgname} ${cfgname} )
    endif()

    # test_geo
    if ( FEELPP_TEST_GEO )
      foreach(  geo ${FEELPP_TEST_GEO} )
        # extract geo filename  to be copied in binary dir
        get_filename_component( GEO_NAME ${geo} NAME )
        configure_file( ${geo} ${GEO_NAME} )
        if(DEFINED ENV{FEELPP_WORKDIR})
          configure_file( ${geo} $ENV{FEELPP_WORKDIR}/geo/${GEO_NAME})
        else(DEFINED ENV{FEELPP_WORKDIR})
          configure_file( ${geo} $ENV{HOME}/feel/geo/${GEO_NAME})
        endif(DEFINED ENV{FEELPP_WORKDIR})
      endforeach()
    endif(FEELPP_TEST_GEO)

endmacro(feelpp_add_test)

#
# feelpp_add_python_module
#
macro(feelpp_ot_add_python_module)
  if ( FEELPP_ENABLE_OPENTURNS AND OPENTURNS_FOUND )
    PARSE_ARGUMENTS(FEELPP_OT_PYTHON
      "LINK_LIBRARIES;SCRIPTS;XML;CFG"
      "TEST"
      ${ARGN}
      )
    CAR(FEELPP_OT_PYTHON_NAME ${FEELPP_OT_PYTHON_DEFAULT_ARGS})
    CDR(FEELPP_OT_PYTHON_SOURCES ${FEELPP_OT_PYTHON_DEFAULT_ARGS})

    add_library( ${FEELPP_OT_PYTHON_NAME} MODULE  ${FEELPP_OT_PYTHON_SOURCES}  )
    target_link_libraries( ${FEELPP_OT_PYTHON_NAME} ${FEELPP_OT_PYTHON_LINK_LIBRARIES}  )
    set_target_properties( ${FEELPP_OT_PYTHON_NAME} PROPERTIES PREFIX "" )
    set_property(TARGET ${FEELPP_OT_PYTHON_NAME} PROPERTY LABELS feelpp)
    #configure_file(${FEELPP_OT_PYTHON_NAME}.xml.in ${FEELPP_OT_PYTHON_NAME}.xml)

    #  add_dependencies(feelpp ${FEELPP_OT_PYTHON_NAME})

    install(TARGETS ${FEELPP_OT_PYTHON_NAME} DESTINATION lib/openturns/wrappers/ COMPONENT Bin)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${FEELPP_OT_PYTHON_NAME}.xml" DESTINATION lib/openturns/wrappers/ COMPONENT Bin)

    if ( FEELPP_OT_PYTHON_SCRIPTS )
      foreach(  script ${FEELPP_OT_PYTHON_SCRIPTS} )
        configure_file( ${script} ${script} )
        if ( FEELPP_OT_PYTHON_TEST )
          add_test(${script} ${PYTHON_EXECUTABLE} ${script})
          set_property(TEST ${script} PROPERTY LABELS feelpp)
        endif()
      endforeach()
    endif()
    if ( FEELPP_OT_PYTHON_CFG )
      foreach(  cfg ${FEELPP_OT_PYTHON_CFG} )
        configure_file( ${cfg} ${cfg} )
        INSTALL(FILES "${cfg}"  DESTINATION share/feel/config)
      endforeach()
    endif()
  endif( FEELPP_ENABLE_OPENTURNS AND OPENTURNS_FOUND )
endmacro(feelpp_ot_add_python_module)

MACRO(find_directories_containing result_list)
  PARSE_ARGUMENTS(FILE
    "FILTER"
    ${ARGN}
    )
  CAR(FILE_NAME ${FILE_DEFAULT_ARGS})

  if ( FEELPP_ENABLE_VERBOSE_CMAKE )
    MESSAGE("*** Arguments for find_directories_containing " ${FILE_DEFAULT_ARGS})
    MESSAGE("  Filters: ${FILE_FILTER}")
  endif()

  FILE(GLOB_RECURSE new_list ${FILE_NAME})
  SET(dir_list "")
  FOREACH(file_path ${new_list})
    GET_FILENAME_COMPONENT(dir_path ${file_path} PATH)
    #MESSAGE("    dir_path: ${dir_path}")
    if ( FILE_FILTER )
      if ( "${dir_path}" MATCHES "(.*)${FILE_FILTER}/(.*)" )
        #MESSAGE(STATUS "${dir_path}")
        LIST(APPEND dir_list ${dir_path})
      endif()
    endif()
  ENDFOREACH()
  LIST(REMOVE_DUPLICATES dir_list)
  #MESSAGE("    LIST: ${dir_list}")
  SET(${result_list} ${dir_list})
ENDMACRO(find_directories_containing)

#
# list of sub-directories of curdir
MACRO(feelpp_list_subdir result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

#
# compute the max of two variables
macro(feelpp_max max var1 var2 )
  if ( ${var1} GREATER ${var2})
    set(${max} ${var1})
  else()
    set(${max} ${var2})
  endif()
endmacro(feelpp_max)

#
# compute the min of two variables
macro(feelpp_min min var1 var2 )
  if ( ${var1} GREATER ${var2})
    set(${min} ${var2})
  else()
    set(${min} ${var1})
  endif()
endmacro(feelpp_min)

# This macros cleans up a variable containing a list of paths
# It:
# - Removes any reference to the original git source directory used for builds (important for instal with tarball)
# - Removes any reference to the original build directory (important for install with tarball)
function(feelpp_clean_variable old_var new_var)
    set(tmp_var "")
    foreach(_entry ${old_var})
        # Try to find build dir reference
        set(_found_position "-1")
        string(FIND ${_entry} ${CMAKE_BINARY_DIR} _found_position)
        if(NOT (${_found_position} MATCHES "0") )
            # Try to find source dir reference
            set(_found_position "-1")
            string(FIND ${_entry} ${CMAKE_SOURCE_DIR} _found_position)
            if(NOT (${_found_position} MATCHES "0"))
                set(tmp_var "${tmp_var};${_entry}")
            endif()
        endif()
    endforeach()
    set(${new_var} ${tmp_var} PARENT_SCOPE)
endfunction(feelpp_clean_variable)

function(feelpp_split_libs libs libnames libpaths)
    set(_paths "")
    set(_names "")
    foreach(_lib ${libs})
        get_filename_component(_path ${_lib} PATH)
        get_filename_component(_name ${_lib} NAME)
        set(_paths ${_paths} ${_path})
        set(_names ${_names} ${_name})
    endforeach()

    set(${libnames} ${_names} PARENT_SCOPE)
    set(${libpaths} ${_paths} PARENT_SCOPE)
endfunction(feelpp_split_libs)

macro(feel_append_src DIRNAME FILES)
  foreach(FILE ${FILES})
    list(APPEND LIST ${DIRNAME}/${FILE})
  endforeach(FILE)
  set(FEELPP_SRCS ${FEELPP_SRCS};${LIST} PARENT_SCOPE)
  set(FEELPP_DIRS ${FEELPP_DIRS};${DIRNAME} PARENT_SCOPE)
endmacro(feel_append_src)

# This function set two variables
# <prefix_name>_LIBRARIES
# <prefix_name>_LIBRARY_DIRS
# from a target
# Usage example:
#     feelpp_expand_target_libraries( FEELPP_VTK ${VTK_LIBRARIES})
#
macro( feelpp_expand_target_libraries prefix_name)
    set( ${prefix_name}_LIBRARIES )
    set( ${prefix_name}_LIBRARY_DIRS )
    # Search VTK_LIBRARIES full path (for cling).
    foreach(LIB ${ARGN})
        # We check if it is a path to the shared lib.
        string( FIND "${LIB}" "/" ISPATH )
        if( "${ISPATH}" STRGREATER "-1") # is a path
            list( APPEND ${prefix_name}_LIBRARIES ${LIB} )
            get_filename_component( LIBDIR ${LIB} DIRECTORY )
            list( APPEND ${prefix_name}_LIBRARY_DIRS ${LIBDIR} )
        else() # not a path
            get_target_property( LIBSO ${LIB} LOCATION )
            if(NOT LIBSO)
                list( APPEND ${prefix_name}_LIBRARIES ${LIB} )
                message( STATUS "Shared library not found: -${LIB}")
                if( $ENV{VERBOSE} )
                    message( STATUS "-${LIB} => ${LIBSO}")
                endif()
            else()
                list( APPEND ${prefix_name}_LIBRARIES ${LIBSO} )
                get_filename_component( LIBDIR ${LIBSO} DIRECTORY )
                list( APPEND ${prefix_name}_LIBRARY_DIRS ${LIBDIR} )
            endif()
        endif()
    endforeach()
    list( REMOVE_DUPLICATES ${prefix_name}_LIBRARIES )
    list( REMOVE_DUPLICATES ${prefix_name}_LIBRARY_DIRS )
    if( $ENV{VERBOSE} )
        message( "-----------------------------------------------------------" )
        message( "FEELPP_VTK_LIBRARY_DIRS=${${prefix_name}_LIBRARY_DIRS}" )
        message( "FEELPP_VTK_LIBRARIES=${${prefix_name}_LIBRARIES}" )
        message( "-----------------------------------------------------------" )
    endif()
endmacro()

if (NOT TARGET man )
  add_custom_target (man)
endif()
if (NOT TARGET html)
  add_custom_target (html)
endif()
macro (feelpp_add_man NAME MAN SECT)
   if (NOT FEELPP_HAS_ASCIIDOCTOR )
     return()
   endif()
    message(STATUS "building manual page ${NAME}.${SECT}")

    if ( FEELPP_HAS_ASCIIDOCTOR_MANPAGE )
      add_custom_target(${NAME}.${SECT})

      add_custom_command (
        TARGET ${NAME}.${SECT}
        #OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.${SECT}
        COMMAND ${FEELPP_A2M} -o ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.${SECT} ${CMAKE_CURRENT_SOURCE_DIR}/${MAN}.adoc
        MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/${MAN}.adoc
        )
      #add_custom_target(${NAME}.${SECT} DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.${SECT})
      if (TARGET man)
        add_dependencies(man ${NAME}.${SECT})
      endif()
      if ( TARGET ${NAME} )
        add_dependencies(${NAME} ${NAME}.${SECT})
      endif()
      install(CODE "execute_process(COMMAND \"bash\" \"-c\" \"${FEELPP_A2M_STR} -o ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.${SECT} ${CMAKE_CURRENT_SOURCE_DIR}/${MAN}.adoc\")" COMPONENT Bin)
      

      install (
        FILES ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.${SECT}
        DESTINATION ${CMAKE_INSTALL_MANDIR}/man${SECT}
        COMPONENT Bin
        )
    endif()
    if ( FEELPP_HAS_ASCIIDOCTOR_HTML5 )
      add_custom_target(${NAME}.${SECT}.html)
      add_custom_command (
        TARGET ${NAME}.${SECT}.html
        #OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.${SECT}.html
        COMMAND ${FEELPP_A2H} -o ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.${SECT}.html ${CMAKE_CURRENT_SOURCE_DIR}/${MAN}.adoc
        DEPENDS ${FEELPP_STYLESHEET}
        MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/${MAN}.adoc
        )
      #add_custom_target(${NAME}.${SECT}.html DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.${SECT}.html)
      if (TARGET html)
        add_dependencies(html ${NAME}.${SECT}.html)
      endif()
      if ( TARGET ${NAME} )
        add_dependencies(${NAME} ${NAME}.${SECT}.html)

      endif()
      install(CODE "execute_process(COMMAND bash \"-c\"  \"${FEELPP_A2H_STR} -o ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.${SECT}.html ${CMAKE_CURRENT_SOURCE_DIR}/${MAN}.adoc\" )" COMPONENT Bin)      
      install (
        FILES ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.${SECT}.html
        DESTINATION ${CMAKE_INSTALL_DOCDIR}
        COMPONENT Bin
        )
      endif()

  endmacro (feelpp_add_man)

  # CRB cmake macros
  include(feelpp.macros.crb)
