# - Find Feel

INCLUDE(ParseArguments)

macro(opus_add_octave_module)

  PARSE_ARGUMENTS(OCTAVE_MODULE
    "LINK_LIBRARIES;SCRIPTS;CFG"
    ""
    ${ARGN}
    )
  CAR(OCTAVE_MODULE_NAME ${OCTAVE_MODULE_DEFAULT_ARGS})
  CDR(OCTAVE_MODULE_SOURCES ${OCTAVE_MODULE_DEFAULT_ARGS})

#  MESSAGE("*** Arguments for Octave module ${OCTAVE_MODULE_NAME}")
#  MESSAGE("    Sources: ${OCTAVE_MODULE_SOURCES}")
#  MESSAGE("    Link libraries: ${OCTAVE_MODULE_LINK_LIBRARIES}")
#  MESSAGE("    Scripts: ${OCTAVE_MODULE_SCRIPTS}")
#  MESSAGE("    Cfg file: ${OCTAVE_MODULE_CFG}")
  set(octname ${OCTAVE_MODULE_NAME}.oct)
  add_library(${octname}  MODULE  ${OCTAVE_MODULE_SOURCES}  )
  target_link_libraries( ${octname} ${OCTAVE_MODULE_LINK_LIBRARIES} )
  set_target_properties( ${octname} PROPERTIES PREFIX "" )
  set_target_properties( ${octname} PROPERTIES SUFFIX "" )
  set_property(TARGET ${octname} PROPERTY LABELS opus)
  INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${octname} DESTINATION ${FEELPP_OCT_DIR})

  add_dependencies(opus ${octname})
  if ( OCTAVE_MODULE_SCRIPTS )
    foreach(  script ${OCTAVE_MODULE_SCRIPTS} )
      configure_file( ${script} ${script} )
      #add_test(${script} ${OCTAVE} ${scipt})
      #set_property(TEST ${script} PROPERTY LABELS opus)

    endforeach()
    INSTALL(FILES ${OCTAVE_MODULE_SCRIPTS}  DESTINATION ${FEELPP_M_DIR})

  endif()
  if ( OCTAVE_MODULE_CFG )
    foreach(  cfg ${OCTAVE_MODULE_CFG} )
      configure_file( ${cfg} ${cfg} )
      INSTALL(FILES "${cfg}"  DESTINATION share/feel/config)
    endforeach()
  endif()
endmacro(opus_add_octave_module)

macro(opus_add_executable)

  PARSE_ARGUMENTS(OPUS_EXEC
    "LINK_LIBRARIES;CFG"
    "NO_TEST"
    ${ARGN}
    )
  CAR(OPUS_EXEC_NAME ${OPUS_EXEC_DEFAULT_ARGS})
  CDR(OPUS_EXEC_SOURCES ${OPUS_EXEC_DEFAULT_ARGS})

#  MESSAGE("*** Arguments for Opus application ${OPUS_EXEC_NAME}")
#  MESSAGE("    Sources: ${OPUS_EXEC_SOURCES}")
#  MESSAGE("    Link libraries: ${OPUS_EXEC_LINK_LIBRARIES}")
#  MESSAGE("    Scripts: ${OPUS_EXEC_SCRIPTS}")
#  MESSAGE("    Cfg file: ${OPUS_EXEC_CFG}")

  set(execname opus_${OPUS_EXEC_NAME})
  add_executable(${execname}    ${OPUS_EXEC_SOURCES}  )
  target_link_libraries( ${execname} ${OPUS_EXEC_LINK_LIBRARIES} )
  set_property(TARGET ${execname} PROPERTY LABELS opus)
  INSTALL(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/${execname}"  DESTINATION bin COMPONENT Bin)
  add_test(${execname} ${CMAKE_CURRENT_BINARY_DIR}/${execname})
  set_property(TEST ${execname} PROPERTY LABELS opus)
  add_dependencies(opus ${execname})
  if ( OPUS_EXEC_CFG )
    foreach(  cfg ${OPUS_EXEC_CFG} )
#      if ( EXISTS ${cfg} )
        configure_file( ${cfg} ${cfg} )
          INSTALL(FILES "${cfg}"  DESTINATION share/feel/config)
#      else()
#        message(WARNING "Executable ${OPUS_EXEC_NAME}: configuration file ${cfg} does not exist")
#      endif()
    endforeach()
  endif()
endmacro(opus_add_executable)

#
# opus_add_python_module
#
macro(opus_add_python_module)

  PARSE_ARGUMENTS(PYTHON
    "LINK_LIBRARIES;SCRIPTS;XML;CFG"
    ""
    ${ARGN}
    )
  CAR(PYTHON_NAME ${PYTHON_DEFAULT_ARGS})
  CDR(PYTHON_SOURCES ${PYTHON_DEFAULT_ARGS})

  add_library( ${PYTHON_NAME} MODULE  ${PYTHON_SOURCES}  )
  target_link_libraries( ${PYTHON_NAME} feel++_opus_models  )
  set_target_properties( ${PYTHON_NAME} PROPERTIES PREFIX "" )
  set_property(TARGET ${PYTHON_NAME} PROPERTY LABELS opus)
  #configure_file(${PYTHON_NAME}.xml.in ${PYTHON_NAME}.xml)

  add_dependencies(opus ${PYTHON_NAME})

  install(TARGETS ${PYTHON_NAME} DESTINATION lib/openturns/wrappers/ COMPONENT Bin)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_NAME}.xml" DESTINATION lib/openturns/wrappers/ COMPONENT Bin)

  if ( PYTHON_SCRIPTS )
    foreach(  script ${PYTHON_SCRIPTS} )
      configure_file( ${script} ${script} )
      add_test(${script} ${PYTHON_EXECUTABLE} ${script})
      set_property(TEST ${script} PROPERTY LABELS opus)
    endforeach()
  endif()
  if ( PYTHON_CFG )
    foreach(  cfg ${PYTHON_CFG} )
      configure_file( ${cfg} ${cfg} )
      INSTALL(FILES "${cfg}"  DESTINATION share/feel/config)
    endforeach()
  endif()

endmacro(opus_add_python_module)

#
# opus_add_model -
#
# generate all C++ files for the various flavors: pfem, crb, scm
# generate all wrappers for python and octave for each flavors
#
macro(opus_add_model)

  PARSE_ARGUMENTS(OPUS_MODEL
    "HDRS;SRCS;LINK_LIBRARIES;CFG;XML;SCRIPTS"
    ""
    ${ARGN}
    )
  CAR(OPUS_MODEL_SHORT_NAME ${OPUS_MODEL_DEFAULT_ARGS})
  CDR(OPUS_MODEL_LONG_NAME ${OPUS_MODEL_DEFAULT_ARGS})

  MESSAGE("*** Arguments for Opus models ${OPUS_MODEL_SHORT_NAME}(${OPUS_MODEL_LONG_NAME})")
  MESSAGE("    Headers: ${OPUS_MODEL_HDRS}")
  MESSAGE("    Sources: ${OPUS_MODEL_SRCS}")
  #MESSAGE("    Link libraries: ${OPUS_MODEL_LINK_LIBRARIES}")
  MESSAGE("    Cfg file: ${OPUS_MODEL_CFG}")
  MESSAGE("    Xml file: ${OPUS_MODEL_XML}")
  MESSAGE("Scripts file: ${OPUS_MODEL_SCRIPTS}")

  include_directories( ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} )
  # generate pfem
  set(CODE "/* this file is generated automatically */
#include <${OPUS_MODEL_SHORT_NAME}.hpp>
#include <feel/feelcrb/opusapp.hpp>

int main( int argc, char** argv )
{
    Feel::OpusApp<Feel::${OPUS_MODEL_LONG_NAME}> app( argc, argv,
                                                      Feel::make${OPUS_MODEL_LONG_NAME}About( \"${OPUS_MODEL_SHORT_NAME}\" ),
                                                      Feel::make${OPUS_MODEL_LONG_NAME}Options()  )\;
    app.run()\;
}
" )
  IF ( EXISTS ${OPUS_MODEL_SHORT_NAME}app.cpp)
    file(READ ${OPUS_MODEL_SHORT_NAME}app.cpp APPCODE)
    IF (NOT ${CODE} STREQUAL ${APPCODE} )
      file(WRITE ${OPUS_MODEL_SHORT_NAME}app.cpp ${CODE})
    ENDIF()
  ELSE()
    file(WRITE ${OPUS_MODEL_SHORT_NAME}app.cpp ${CODE})
  ENDIF()

  opus_add_executable(${OPUS_MODEL_SHORT_NAME}app ${OPUS_MODEL_SHORT_NAME}app.cpp
    LINK_LIBRARIES ${OPUS_MODEL_LINK_LIBRARIES}
    CFG ${OPUS_MODEL_CFG})


  foreach( wrapper pfem scm crb )
    set(pycpp "${OPUS_MODEL_SHORT_NAME}${wrapper}_pywrapper.cpp")
    set(octcpp "${OPUS_MODEL_SHORT_NAME}${wrapper}_octwrapper.cpp")
    set(xml "${OPUS_MODEL_SHORT_NAME}${wrapper}.xml")
    set(OPUS_MODEL_WRAPPER_NAME "opus${OPUS_MODEL_SHORT_NAME}${wrapper}")
    set(OPUS_MODEL_WRAPPER_TYPE "\"${wrapper}\"")
    configure_file(${FEELPP_SOURCE_DIR}/applications/crb/templates/python_wrapper.cpp ${pycpp})
    configure_file(${FEELPP_SOURCE_DIR}/applications/crb/templates/octave_wrapper.cpp ${octcpp})
    configure_file(${OPUS_MODEL_SHORT_NAME}.xml.in ${xml})

    opus_add_python_module(opus${OPUS_MODEL_SHORT_NAME}${wrapper} ${pycpp}
      LINK_LIBRARIES ${OPUS_MODEL_LINK_LIBRARIES}
      CFG ${OPUS_MODEL_CFG} XML ${xml})
    opus_add_octave_module(opus${OPUS_MODEL_SHORT_NAME}${wrapper} ${octcpp}
      LINK_LIBRARIES ${OPUS_MODEL_LINK_LIBRARIES} ${Octave_LIBRARIES}
      CFG ${OPUS_MODEL_CFG} SCRIPTS ${OPUS_MODEL_SCRIPTS} )
  endforeach()



endmacro()
