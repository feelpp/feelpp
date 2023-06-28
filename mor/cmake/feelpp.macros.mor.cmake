# - Find Feel

INCLUDE(ParseArguments)

macro(mor_add_octave_module)
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
  target_link_libraries( ${octname} Feelpp::feelpp_mor  ${OCTAVE_MODULE_LINK_LIBRARIES} )
  set_target_properties( ${octname} PROPERTIES PREFIX "" )
  set_target_properties( ${octname} PROPERTIES SUFFIX "" )
  set_property(TARGET ${octname} PROPERTY LABELS mor)
  INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${octname} DESTINATION ${FEELPP_OCT_DIR})

  #add_dependencies(crb ${octname})
  if ( OCTAVE_MODULE_SCRIPTS )
    foreach(  script ${OCTAVE_MODULE_SCRIPTS} )
      configure_file( ${script} ${script} )
      #add_test(${script} ${OCTAVE} ${scipt})
      #set_property(TEST ${script} PROPERTY LABELS mor)

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
endmacro(mor_add_octave_module)

macro(mor_add_executable)

  PARSE_ARGUMENTS(mor_EXEC
    "SOURCES;LINK_LIBRARIES;CFG;GEO;PROJECT;EXEC;MAN"
    "TEST"
    ${ARGN}
    )
  CAR(mor_EXEC_NAME ${mor_EXEC_DEFAULT_ARGS})
  CDR(mor_EXEC_SOURCES ${mor_EXEC_DEFAULT_ARGS})

  if ( FEELPP_ENABLE_VERBOSE_CMAKE )
    MESSAGE("*** Arguments for Crb application ${mor_EXEC_NAME}")
    MESSAGE("    Sources: ${mor_EXEC_SOURCES}")
    MESSAGE("    Link libraries: ${mor_EXEC_LINK_LIBRARIES}")
    MESSAGE("    Scripts: ${mor_EXEC_SCRIPTS}")
    MESSAGE("    Cfg file: ${mor_EXEC_CFG}")
    MESSAGE("    Geo file: ${mor_EXEC_GEO}")
  endif()

  if ( mor_EXEC_PROJECT )
    set(execname feelpp_mor_${mor_EXEC_PROJECT}_${mor_EXEC_NAME})
  else()
    set(execname feelpp_mor_${mor_EXEC_NAME})
  endif()
  if  (mor_EXEC_EXEC )
    set( ${mor_EXEC_EXEC} ${execname} )
  endif()
  add_executable(${execname}    ${mor_EXEC_SOURCES} )
  target_link_libraries( ${execname} ${mor_EXEC_LINK_LIBRARIES} Feelpp::feelpp_mor   )
  set_property(TARGET ${execname} PROPERTY LABELS mor)
  install(TARGETS ${execname} RUNTIME DESTINATION bin COMPONENT Bin)

  # add manual page
  if ( mor_EXEC_MAN )
    feelpp_add_man( ${execname} ${mor_EXEC_MAN} 1 )
  endif( mor_EXEC_MAN )

  if ( mor_EXEC_TEST )
    add_test(${execname} ${CMAKE_CURRENT_BINARY_DIR}/${execname})
    set_property(TEST ${execname} PROPERTY LABELS mor)
  endif()
  #add_dependencies(crb ${execname})
  if ( mor_EXEC_CFG )
    foreach(  cfg ${mor_EXEC_CFG} )
#      if ( EXISTS ${cfg} )
        configure_file( ${cfg} ${cfg} )
          INSTALL(FILES "${cfg}"  DESTINATION share/feel/config)
#      else()
#        message(WARNING "Executable ${mor_EXEC_NAME}: configuration file ${cfg} does not exist")
#      endif()
    endforeach()
  endif()
  # geo and mesh
  if ( mor_EXEC_GEO )
    foreach(  geo ${mor_EXEC_GEO} )
      get_filename_component( GEO_NAME ${geo} NAME )
      configure_file( ${geo} ${GEO_NAME} )
    endforeach()
  endif()

endmacro(mor_add_executable)

macro(mor_add_library)

  PARSE_ARGUMENTS(mor_LIB
    "SRCS;LINK_LIBRARIES;PROJECT;EXEC;MAN;EXPORT"
    "TEST;NOHEADER;PLUGIN"
    ${ARGN}
    )
  CAR(mor_LIB_NAME ${mor_LIB_DEFAULT_ARGS})

  if ( FEELPP_ENABLE_VERBOSE_CMAKE )
    MESSAGE("*** Arguments for Crb application ${mor_LIB_NAME}")
    MESSAGE("    Sources: ${mor_LIB_SRCS}")
    MESSAGE("    Link libraries: ${mor_LIB_LINK_LIBRARIES}")
  endif()

  if ( mor_LIB_PLUGIN )
    set(MOR_PREFIX "feelpp_mor_plugin")
  else()
    set(MOR_PLUGIN "feelpp_mor")
  endif()
  if ( mor_LIB_PROJECT )
    set(execname ${MOR_PREFIX}_${mor_LIB_PROJECT}_${mor_LIB_NAME})
  else()
    set(execname ${MOR_PREFIX}_${mor_LIB_NAME})
  endif()
  if  (mor_LIB_EXEC )
    set( ${mor_LIB_EXEC} ${execname} )
  endif()
  add_library(${execname}  SHARED  ${mor_LIB_SRCS} )
  if ( NOT mor_LIB_PLUGIN )
    set_target_properties(${execname} PROPERTIES VERSION 1 SOVERSION 1)
  endif()
  # -fvisibility=hidden not always work (macos, toolbox_heat...)
  #target_compile_options(${execname} PRIVATE -fvisibility=hidden)

  target_include_directories( ${execname} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/feelpp/mor/${mor_LIB_NAME}>  )

  target_link_libraries( ${execname} PUBLIC ${mor_LIB_LINK_LIBRARIES} Feelpp::feelpp_mor   )
  set_property(TARGET ${execname} PROPERTY LABELS mor)

  if(mor_LIB_PLUGIN)
    INSTALL(TARGETS ${execname} EXPORT ${mor_LIB_EXPORT}
      LIBRARY DESTINATION ${FEELPP_LIBDIR}
      INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/feelpp/mor/${mor_LIB_NAME}
      )
  else()
    INSTALL(TARGETS ${execname} EXPORT ${mor_LIB_EXPORT}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/feelpp/mor/${mor_LIB_NAME}
      )
  endif()
endmacro(mor_add_library)
#
# mor_add_python_module
#
macro(mor_add_python_module)
if ( FEELPP_HAS_OPENTURNS )
  PARSE_ARGUMENTS(mor_PYTHON
    "LINK_LIBRARIES;SCRIPTS;XML;CFG;CLASS"
    "TEST"
    ${ARGN}
    )
  CAR(mor_PYTHON_NAME ${mor_PYTHON_DEFAULT_ARGS})
  CDR(mor_PYTHON_SOURCES ${mor_PYTHON_DEFAULT_ARGS})

  add_library( ${mor_PYTHON_NAME} MODULE  ${mor_PYTHON_SOURCES}  )
  target_link_libraries( ${mor_PYTHON_NAME} Feelpp::feelpp_mor  ${mor_PYTHON_LINK_LIBRARIES}  ${OpenTURNS_LIBRARIES} )
  set_target_properties( ${mor_PYTHON_NAME} PROPERTIES PREFIX "" )
  set_property(TARGET ${mor_PYTHON_NAME} PROPERTY LABELS mor)
  #configure_file(${mor_PYTHON_NAME}.xml.in ${mor_PYTHON_NAME}.xml)

  #add_dependencies(crb ${mor_PYTHON_NAME})

  if ( 0 ) #TODO
  install(TARGETS ${mor_PYTHON_NAME} DESTINATION lib/openturns/wrappers/ COMPONENT Bin)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${mor_PYTHON_NAME}.xml" DESTINATION lib/openturns/wrappers/ COMPONENT Bin)
  endif()

  if ( mor_PYTHON_SCRIPTS )
    foreach(  script ${mor_PYTHON_SCRIPTS} )
      configure_file( ${script} ${script} )
      if ( mor_PYTHON_TEST )
        add_test(${script} ${PYTHON_EXECUTABLE} ${script})
        set_property(TEST ${script} PROPERTY LABELS mor)
      endif()
    endforeach()
  endif()
  if ( mor_PYTHON_CFG )
    foreach(  cfg ${mor_PYTHON_CFG} )
      configure_file( ${cfg} ${cfg} )
      INSTALL(FILES "${cfg}"  DESTINATION share/feel/config)
    endforeach()
  endif()
endif( FEELPP_HAS_OPENTURNS )
endmacro(mor_add_python_module)

#
# mor_add_model -
#
# generate all C++ files for the various flavors: pfem, crb, scm
# generate all wrappers for python and octave for each flavors
#
macro(mor_add_model)

  PARSE_ARGUMENTS(mor_MODEL
    "HDRS;SRCS;LINK_LIBRARIES;CFG;XML;SCRIPTS;CLASS;DEFS;GEO;MSH;OPT"
    "TEST;ADD_OT"
    ${ARGN}
    )
  CAR(mor_MODEL_SHORT_NAME ${mor_MODEL_DEFAULT_ARGS})
  CDR(mor_MODEL_LONG_NAME ${mor_MODEL_DEFAULT_ARGS})

  if (NOT DEFINED mor_MODEL_OPT)
    set(mor_MODEL_OPT ${mor_MODEL_SHORT_NAME})
  endif()

  if ( FEELPP_ENABLE_VERBOSE_CMAKE )
    MESSAGE("*** Arguments for Crb models ${mor_MODEL_SHORT_NAME}(${mor_MODEL_LONG_NAME})")
    MESSAGE("  q  Headers: ${mor_MODEL_HDRS}")
    MESSAGE("    Sources: ${mor_MODEL_SRCS}")
    MESSAGE("    Defs file: ${mor_MODEL_DEFS}")
    #MESSAGE("    Link libraries: ${mor_MODEL_LINK_LIBRARIES}")
    MESSAGE("    Cfg file: ${mor_MODEL_CFG}")
    MESSAGE("    Xml file: ${mor_MODEL_XML}")
    MESSAGE("    Geo file: ${mor_MODEL_GEO}")
    MESSAGE("    Msh file: ${mor_MODEL_MSH}")
    MESSAGE("Scripts file: ${mor_MODEL_SCRIPTS}")
  endif()

  include_directories( ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} )
  if ( NOT mor_MODEL_CLASS )
    set(mor_MODEL_CLASS ${mor_MODEL_LONG_NAME})
  endif()
    # generate pfem
  set(CODE "/* this file is generated automatically */
#include <${mor_MODEL_HDRS}>
#include <feel/feelmor/opusapp.hpp>

int main( int argc, char** argv )
{
    using namespace Feel\;
    Feel::Environment env( _argc = argc, _argv = argv,
                           _desc = opusapp_options(\"${mor_MODEL_OPT}\")
                           .add(crbOptions())
                           .add(crbSEROptions())
                           .add(make${mor_MODEL_LONG_NAME}Options())
                           .add(eimOptions())
                           .add(podOptions())
                           .add(backend_options(\"backend-primal\"))
                           .add(backend_options(\"backend-dual\"))
                           .add(backend_options(\"backend-l2\"))
                           .add(bdf_options(\"${mor_MODEL_LONG_NAME}\")),
                           _about = make${mor_MODEL_LONG_NAME}About( \"${mor_MODEL_OPT}\" ) )\;

    Feel::OpusApp<Feel::${mor_MODEL_CLASS} > app\;
    app.run()\;
}
" )

  IF ( EXISTS ${mor_MODEL_SHORT_NAME}app.cpp)
    file(READ ${mor_MODEL_SHORT_NAME}app.cpp APPCODE)
    IF (NOT ${CODE} STREQUAL ${APPCODE} )
      file(WRITE ${mor_MODEL_SHORT_NAME}app.cpp ${CODE})
    ENDIF()
  ELSE()
    file(WRITE ${mor_MODEL_SHORT_NAME}app.cpp ${CODE})
  ENDIF()

  if ( mor_MODEL_TEST )
    mor_add_executable(${mor_MODEL_SHORT_NAME}app
      ${mor_MODEL_SHORT_NAME}app.cpp ${mor_MODEL_SRCS}
      GEO ${mor_MODEL_GEO}
      LINK_LIBRARIES ${mor_MODEL_LINK_LIBRARIES}
      CFG ${mor_MODEL_CFG} TEST )
  else()
    mor_add_executable(${mor_MODEL_SHORT_NAME}app
      ${mor_MODEL_SHORT_NAME}app.cpp ${mor_MODEL_SRCS}
      GEO ${mor_MODEL_GEO}
      LINK_LIBRARIES ${mor_MODEL_LINK_LIBRARIES}
      CFG ${mor_MODEL_CFG} )
  endif()

  # include schedulers
  include( feelpp.schedulers )

  foreach( wrapper pfem scm crb )
    set(pycpp "${mor_MODEL_SHORT_NAME}${wrapper}_pywrapper.cpp")
    set(octcpp "${mor_MODEL_SHORT_NAME}${wrapper}_octwrapper.cpp")
    set(xml "${mor_MODEL_SHORT_NAME}${wrapper}.xml")
    set(mor_MODEL_WRAPPER_NAME "crb${mor_MODEL_SHORT_NAME}${wrapper}")
    set(mor_MODEL_WRAPPER_TYPE "\"${wrapper}\"")
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

    if ( mor_MODEL_DEFS )
      set_property(TARGET ${execname} PROPERTY COMPILE_DEFINITIONS ${mor_MODEL_DEFS})
    endif()

    if (EXISTS ${mor_MODEL_SHORT_NAME}.xml.in )
      configure_file(${mor_MODEL_SHORT_NAME}.xml.in ${xml})

      if ( mor_MODEL_TEST )
        mor_add_python_module(crb${mor_MODEL_SHORT_NAME}${wrapper} ${pycpp}
          LINK_LIBRARIES ${mor_MODEL_LINK_LIBRARIES} CLASS ${mor_MODEL_CLASS}
          CFG ${mor_MODEL_CFG} XML ${xml} TEST)
        mor_add_octave_module(crb${mor_MODEL_SHORT_NAME}${wrapper} ${octcpp}
          LINK_LIBRARIES ${mor_MODEL_LINK_LIBRARIES} ${Octave_LIBRARIES}
          CFG ${mor_MODEL_CFG} SCRIPTS ${mor_MODEL_SCRIPTS} TEST )
      else()
        mor_add_python_module(crb${mor_MODEL_SHORT_NAME}${wrapper} ${pycpp}
          LINK_LIBRARIES ${mor_MODEL_LINK_LIBRARIES} CLASS ${mor_MODEL_CLASS}
          CFG ${mor_MODEL_CFG} XML ${xml})
        mor_add_octave_module(crb${mor_MODEL_SHORT_NAME}${wrapper} ${octcpp}
          LINK_LIBRARIES ${mor_MODEL_LINK_LIBRARIES} ${Octave_LIBRARIES}
          CFG ${mor_MODEL_CFG} SCRIPTS ${mor_MODEL_SCRIPTS} )
      endif()
    endif()
  endforeach()

  # Install OpenCL source files
  if ( HARTS_LIBRARIES AND ENABLE_OPENCL )
      set(mor_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/feel/feelmor")
      if ( EXISTS ${mor_INCLUDE_DIR} )
          file(GLOB OPENCL_SOURCE_FILES "${mor_INCLUDE_DIR}/*.cl")
          if(OPENCL_SOURCE_FILES)
              file(COPY ${OPENCL_SOURCE_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
          endif()
      endif()
  endif()

endmacro()
