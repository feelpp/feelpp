# - Find Feel

INCLUDE(ParseArguments)

macro(crb_add_octave_module)
if ( FEELPP_HAS_OCTAVE )
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
  target_link_libraries( ${octname} Feelpp::feelpp  ${OCTAVE_MODULE_LINK_LIBRARIES} )
  set_target_properties( ${octname} PROPERTIES PREFIX "" )
  set_target_properties( ${octname} PROPERTIES SUFFIX "" )
  set_property(TARGET ${octname} PROPERTY LABELS crb)
  INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${octname} DESTINATION ${FEELPP_OCT_DIR})

  #add_dependencies(crb ${octname})
  if ( OCTAVE_MODULE_SCRIPTS )
    foreach(  script ${OCTAVE_MODULE_SCRIPTS} )
      configure_file( ${script} ${script} )
      #add_test(${script} ${OCTAVE} ${scipt})
      #set_property(TEST ${script} PROPERTY LABELS crb)

    endforeach()
    INSTALL(FILES ${OCTAVE_MODULE_SCRIPTS}  DESTINATION ${FEELPP_M_DIR})

  endif()
  if ( OCTAVE_MODULE_CFG )
    foreach(  cfg ${OCTAVE_MODULE_CFG} )
      configure_file( ${cfg} ${cfg} )
      INSTALL(FILES "${cfg}"  DESTINATION share/feel/config)
    endforeach()
  endif()
endif( FEELPP_HAS_OCTAVE )
endmacro(crb_add_octave_module)

macro(crb_add_executable)

  PARSE_ARGUMENTS(CRB_EXEC
    "SOURCES;LINK_LIBRARIES;CFG;GEO;PROJECT;EXEC;MAN"
    "TEST"
    ${ARGN}
    )
  CAR(CRB_EXEC_NAME ${CRB_EXEC_DEFAULT_ARGS})
  CDR(CRB_EXEC_SOURCES ${CRB_EXEC_DEFAULT_ARGS})

  if ( FEELPP_ENABLE_VERBOSE_CMAKE )
    MESSAGE("*** Arguments for Crb application ${CRB_EXEC_NAME}")
    MESSAGE("    Sources: ${CRB_EXEC_SOURCES}")
    MESSAGE("    Link libraries: ${CRB_EXEC_LINK_LIBRARIES}")
    MESSAGE("    Scripts: ${CRB_EXEC_SCRIPTS}")
    MESSAGE("    Cfg file: ${CRB_EXEC_CFG}")
    MESSAGE("    Geo file: ${CRB_EXEC_GEO}")
  endif()

  if ( CRB_EXEC_PROJECT )
    set(execname feelpp_crb_${CRB_EXEC_PROJECT}_${CRB_EXEC_NAME})
  else()
    set(execname feelpp_crb_${CRB_EXEC_NAME})
  endif()
  if  (CRB_EXEC_EXEC )
    set( ${CRB_EXEC_EXEC} ${execname} )
  endif()
  add_executable(${execname}    ${CRB_EXEC_SOURCES} )
  target_link_libraries( ${execname} ${CRB_EXEC_LINK_LIBRARIES} Feelpp::feelpp   )
  set_property(TARGET ${execname} PROPERTY LABELS crb)
  install(TARGETS ${execname} RUNTIME DESTINATION bin COMPONENT Bin)

  # add manual page
  if ( CRB_EXEC_MAN )
    feelpp_add_man( ${execname} ${CRB_EXEC_MAN} 1 )
  endif( CRB_EXEC_MAN )

  if ( CRB_EXEC_TEST )
    add_test(${execname} ${CMAKE_CURRENT_BINARY_DIR}/${execname})
    set_property(TEST ${execname} PROPERTY LABELS crb)
  endif()
  #add_dependencies(crb ${execname})
  if ( CRB_EXEC_CFG )
    foreach(  cfg ${CRB_EXEC_CFG} )
#      if ( EXISTS ${cfg} )
        configure_file( ${cfg} ${cfg} )
          INSTALL(FILES "${cfg}"  DESTINATION share/feel/config)
#      else()
#        message(WARNING "Executable ${CRB_EXEC_NAME}: configuration file ${cfg} does not exist")
#      endif()
    endforeach()
  endif()
  # geo and mesh
  if ( CRB_EXEC_GEO )
    foreach(  geo ${CRB_EXEC_GEO} )
      get_filename_component( GEO_NAME ${geo} NAME )
      configure_file( ${geo} ${GEO_NAME} )
    endforeach()
  endif()

endmacro(crb_add_executable)

macro(crb_add_library)

  PARSE_ARGUMENTS(CRB_LIB
    "SRCS;LINK_LIBRARIES;PROJECT;EXEC;MAN;EXPORT"
    "TEST;NOHEADER"
    ${ARGN}
    )
  CAR(CRB_LIB_NAME ${CRB_LIB_DEFAULT_ARGS})

  if ( FEELPP_ENABLE_VERBOSE_CMAKE )
    MESSAGE("*** Arguments for Crb application ${CRB_LIB_NAME}")
    MESSAGE("    Sources: ${CRB_LIB_SRCS}")
    MESSAGE("    Link libraries: ${CRB_LIB_LINK_LIBRARIES}")
  endif()

  if ( CRB_LIB_PROJECT )
    set(execname feelpp_crb_${CRB_LIB_PROJECT}_${CRB_LIB_NAME})
  else()
    set(execname feelpp_crb_${CRB_LIB_NAME})
  endif()
  if  (CRB_LIB_EXEC )
    set( ${CRB_LIB_EXEC} ${execname} )
  endif()
  add_library(${execname}  SHARED  ${CRB_LIB_SRCS} )
  set_target_properties(${execname} PROPERTIES VERSION 1 SOVERSION 1)
  target_compile_options(${execname} PRIVATE -fvisibility=hidden)
  
  target_include_directories( ${execname} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/feelpp/mor/${CRB_LIB_NAME}>  )

  target_link_libraries( ${execname} PUBLIC ${CRB_LIB_LINK_LIBRARIES} Feelpp::feelpp   )
  set_property(TARGET ${execname} PROPERTY LABELS crb)
  INSTALL(TARGETS ${execname} EXPORT ${CRB_LIB_EXPORT}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/feelpp/mor/${CRB_LIB_NAME}
    )

endmacro(crb_add_library)
#
# crb_add_python_module
#
macro(crb_add_python_module)
if ( FEELPP_HAS_OPENTURNS )
  PARSE_ARGUMENTS(CRB_PYTHON
    "LINK_LIBRARIES;SCRIPTS;XML;CFG;CLASS"
    "TEST"
    ${ARGN}
    )
  CAR(CRB_PYTHON_NAME ${CRB_PYTHON_DEFAULT_ARGS})
  CDR(CRB_PYTHON_SOURCES ${CRB_PYTHON_DEFAULT_ARGS})

  add_library( ${CRB_PYTHON_NAME} MODULE  ${CRB_PYTHON_SOURCES}  )
  target_link_libraries( ${CRB_PYTHON_NAME} Feelpp::feelpp  ${CRB_PYTHON_LINK_LIBRARIES}  ${OpenTURNS_LIBRARIES} )
  set_target_properties( ${CRB_PYTHON_NAME} PROPERTIES PREFIX "" )
  set_property(TARGET ${CRB_PYTHON_NAME} PROPERTY LABELS crb)
  #configure_file(${CRB_PYTHON_NAME}.xml.in ${CRB_PYTHON_NAME}.xml)

  #add_dependencies(crb ${CRB_PYTHON_NAME})

  if ( 0 ) #TODO
  install(TARGETS ${CRB_PYTHON_NAME} DESTINATION lib/openturns/wrappers/ COMPONENT Bin)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${CRB_PYTHON_NAME}.xml" DESTINATION lib/openturns/wrappers/ COMPONENT Bin)
  endif()

  if ( CRB_PYTHON_SCRIPTS )
    foreach(  script ${CRB_PYTHON_SCRIPTS} )
      configure_file( ${script} ${script} )
      if ( CRB_PYTHON_TEST )
        add_test(${script} ${PYTHON_EXECUTABLE} ${script})
        set_property(TEST ${script} PROPERTY LABELS crb)
      endif()
    endforeach()
  endif()
  if ( CRB_PYTHON_CFG )
    foreach(  cfg ${CRB_PYTHON_CFG} )
      configure_file( ${cfg} ${cfg} )
      INSTALL(FILES "${cfg}"  DESTINATION share/feel/config)
    endforeach()
  endif()
endif( FEELPP_HAS_OPENTURNS )
endmacro(crb_add_python_module)

#
# crb_add_model -
#
# generate all C++ files for the various flavors: pfem, crb, scm
# generate all wrappers for python and octave for each flavors
#
macro(crb_add_model)

  PARSE_ARGUMENTS(CRB_MODEL
    "HDRS;SRCS;LINK_LIBRARIES;CFG;XML;SCRIPTS;CLASS;DEFS;GEO;MSH"
    "TEST;ADD_OT"
    ${ARGN}
    )
  CAR(CRB_MODEL_SHORT_NAME ${CRB_MODEL_DEFAULT_ARGS})
  CDR(CRB_MODEL_LONG_NAME ${CRB_MODEL_DEFAULT_ARGS})

  if ( FEELPP_ENABLE_VERBOSE_CMAKE )
    MESSAGE("*** Arguments for Crb models ${CRB_MODEL_SHORT_NAME}(${CRB_MODEL_LONG_NAME})")
    MESSAGE("  q  Headers: ${CRB_MODEL_HDRS}")
    MESSAGE("    Sources: ${CRB_MODEL_SRCS}")
    MESSAGE("    Defs file: ${CRB_MODEL_DEFS}")
    #MESSAGE("    Link libraries: ${CRB_MODEL_LINK_LIBRARIES}")
    MESSAGE("    Cfg file: ${CRB_MODEL_CFG}")
    MESSAGE("    Xml file: ${CRB_MODEL_XML}")
    MESSAGE("    Geo file: ${CRB_MODEL_GEO}")
    MESSAGE("    Msh file: ${CRB_MODEL_MSH}")
    MESSAGE("Scripts file: ${CRB_MODEL_SCRIPTS}")
  endif()

  include_directories( ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} )
  if ( NOT CRB_MODEL_CLASS )
    set(CRB_MODEL_CLASS ${CRB_MODEL_LONG_NAME})
  endif()
    # generate pfem
  set(CODE "/* this file is generated automatically */
#include <${CRB_MODEL_HDRS}>
#include <feel/feelcrb/opusapp.hpp>

int main( int argc, char** argv )
{
    using namespace Feel\;
    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=opusapp_options(\"${CRB_MODEL_SHORT_NAME}\")
                           .add(crbOptions())
                           .add(crbSEROptions())
                           .add(make${CRB_MODEL_LONG_NAME}Options())
                           .add(eimOptions())
                           .add(podOptions())
                           .add(backend_options(\"backend-primal\"))
                           .add(backend_options(\"backend-dual\"))
                           .add(backend_options(\"backend-l2\"))
                           .add(bdf_options(\"${CRB_MODEL_LONG_NAME}\")),
                           _about=make${CRB_MODEL_LONG_NAME}About( \"${CRB_MODEL_SHORT_NAME}\" ) )\;

    Feel::OpusApp<Feel::${CRB_MODEL_CLASS} > app\;
    app.run()\;
}
" )

  IF ( EXISTS ${CRB_MODEL_SHORT_NAME}app.cpp)
    file(READ ${CRB_MODEL_SHORT_NAME}app.cpp APPCODE)
    IF (NOT ${CODE} STREQUAL ${APPCODE} )
      file(WRITE ${CRB_MODEL_SHORT_NAME}app.cpp ${CODE})
    ENDIF()
  ELSE()
    file(WRITE ${CRB_MODEL_SHORT_NAME}app.cpp ${CODE})
  ENDIF()

  if ( CRB_MODEL_TEST )
    crb_add_executable(${CRB_MODEL_SHORT_NAME}app
      ${CRB_MODEL_SHORT_NAME}app.cpp ${CRB_MODEL_SRCS}
      GEO ${CRB_MODEL_GEO}
      LINK_LIBRARIES ${CRB_MODEL_LINK_LIBRARIES}
      CFG ${CRB_MODEL_CFG} TEST )
  else()
    crb_add_executable(${CRB_MODEL_SHORT_NAME}app
      ${CRB_MODEL_SHORT_NAME}app.cpp ${CRB_MODEL_SRCS}
      GEO ${CRB_MODEL_GEO}
      LINK_LIBRARIES ${CRB_MODEL_LINK_LIBRARIES}
      CFG ${CRB_MODEL_CFG} )
  endif()

  # include schedulers
  include( feelpp.schedulers )

  foreach( wrapper pfem scm crb )
    set(pycpp "${CRB_MODEL_SHORT_NAME}${wrapper}_pywrapper.cpp")
    set(octcpp "${CRB_MODEL_SHORT_NAME}${wrapper}_octwrapper.cpp")
    set(xml "${CRB_MODEL_SHORT_NAME}${wrapper}.xml")
    set(CRB_MODEL_WRAPPER_NAME "crb${CRB_MODEL_SHORT_NAME}${wrapper}")
    set(CRB_MODEL_WRAPPER_TYPE "\"${wrapper}\"")
    #configure_file(${FEELPP_SOURCE_DIR}/mor/templates/python_wrapper.cpp ${pycpp})
    if ( EXISTS ${CMAKE_SOURCE_DIR}/mor/templates/ )
      configure_file(${CMAKE_SOURCE_DIR}/mor/templates/ot_python_command_wrapper.cpp ${pycpp})
      configure_file(${CMAKE_SOURCE_DIR}/mor/templates/octave_wrapper.cpp ${octcpp})
    elseif( EXISTS ${FEELPP_DATADIR}/crb/templates )
      configure_file(${FEELPP_DATADIR}/crb/templates/ot_python_command_wrapper.cpp ${pycpp})
      configure_file(${FEELPP_DATADIR}/crb/templates/octave_wrapper.cpp ${octcpp})
    else()
      continue()
    endif()

    if ( CRB_MODEL_DEFS )
      set_property(TARGET ${execname} PROPERTY COMPILE_DEFINITIONS ${CRB_MODEL_DEFS})
    endif()

    if (EXISTS ${CRB_MODEL_SHORT_NAME}.xml.in )
      configure_file(${CRB_MODEL_SHORT_NAME}.xml.in ${xml})

      if ( CRB_MODEL_TEST )
        crb_add_python_module(crb${CRB_MODEL_SHORT_NAME}${wrapper} ${pycpp}
          LINK_LIBRARIES ${CRB_MODEL_LINK_LIBRARIES} CLASS ${CRB_MODEL_CLASS}
          CFG ${CRB_MODEL_CFG} XML ${xml} TEST)
        crb_add_octave_module(crb${CRB_MODEL_SHORT_NAME}${wrapper} ${octcpp}
          LINK_LIBRARIES ${CRB_MODEL_LINK_LIBRARIES} ${Octave_LIBRARIES}
          CFG ${CRB_MODEL_CFG} SCRIPTS ${CRB_MODEL_SCRIPTS} TEST )
      else()
        crb_add_python_module(crb${CRB_MODEL_SHORT_NAME}${wrapper} ${pycpp}
          LINK_LIBRARIES ${CRB_MODEL_LINK_LIBRARIES} CLASS ${CRB_MODEL_CLASS}
          CFG ${CRB_MODEL_CFG} XML ${xml})
        crb_add_octave_module(crb${CRB_MODEL_SHORT_NAME}${wrapper} ${octcpp}
          LINK_LIBRARIES ${CRB_MODEL_LINK_LIBRARIES} ${Octave_LIBRARIES}
          CFG ${CRB_MODEL_CFG} SCRIPTS ${CRB_MODEL_SCRIPTS} )
      endif()
    endif()
  endforeach()

  # Install OpenCL source files
  if ( HARTS_LIBRARIES AND ENABLE_OPENCL )
      set(CRB_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/feel/feelcrb")
      if ( EXISTS ${CRB_INCLUDE_DIR} )
          file(GLOB OPENCL_SOURCE_FILES "${CRB_INCLUDE_DIR}/*.cl")
          if(OPENCL_SOURCE_FILES)
              file(COPY ${OPENCL_SOURCE_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
          endif()
      endif()
  endif()

endmacro()
