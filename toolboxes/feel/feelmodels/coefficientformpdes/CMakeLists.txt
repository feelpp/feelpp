

feelpp_toolboxes_add_library( coefficientformpdebase SRCS coefficientformpdebase.cpp coefficientformpdeboundaryconditions.cpp LINK_LIBRARIES feelpp_modelcore )


if (1)
set( _UNKNOWN_BASIS_TYPE_LIST
  #"Lagrange<0,Scalar,Continuous,PointSetFekete>"
  "Lagrange<1,Scalar,Continuous,PointSetFekete>" "Lagrange<2,Scalar,Continuous,PointSetFekete>"
  "Lagrange<1,Vectorial,Continuous,PointSetFekete>" "Lagrange<2,Vectorial,Continuous,PointSetFekete>"
  "Lagrange<1,Scalar,Discontinuous,PointSetFekete>"
  "Nedelec<0,NedelecKind::NED1>"
  )
set( _UNKNOWN_BASIS_TAG_LIST
  #"Pch0"
  "Pch1" "Pch2"
  "Pchv1" "Pchv2"
  "Pdh1"
  "Ned1h0"
  )
else()
set( _UNKNOWN_BASIS_TYPE_LIST "Lagrange<1,Scalar,Continuous,PointSetFekete>")
set( _UNKNOWN_BASIS_TAG_LIST "Pch1")
endif()

unset( FEELPP_TOOLBOX_COEFFICIENTFORMPDES_REGISTER_ALL_ENTRIES_CLASS_TYPE )
# geo order 1
foreach(_DIM 2 3)
  math(EXPR _T_ORDER "${FEELPP_INSTANTIATION_ORDER_MAX}")
#  feelpp_max( MAX_ORDER 1 ${_T_ORDER})
  foreach(_GEO_ORDER 1 2)
    # message(STATUS "[feelpp] COEFFICIENTFORMPDES ${_DIM}D PNG1 N=1..${MAX_ORDER}")
    message(STATUS "[feelpp] COEFFICIENTFORMPDES ${_DIM}D G${_GEO_ORDER}")
    genLibCoefficientFormPDEs(
      DIM     ${_DIM}
      UNKNOWN_BASIS_TYPE ${_UNKNOWN_BASIS_TYPE_LIST}
      UNKNOWN_BASIS_TAG ${_UNKNOWN_BASIS_TAG_LIST}
      GEO_ORDER ${_GEO_ORDER}
      )
    set(FEELPP_MODELS_LIBRARIES_COEFFICIENTFORMPDES  ${FEELPP_MODELS_LIBRARIES_COEFFICIENTFORMPDES} ${COEFFICIENTFORMPDES_LIB_NAME}  )
    if ( NOT DEFINED FEELPP_TOOLBOX_COEFFICIENTFORMPDES_REGISTER_ALL_ENTRIES_CLASS_TYPE )
      set(FEELPP_TOOLBOX_COEFFICIENTFORMPDES_REGISTER_ALL_ENTRIES_CLASS_TYPE  "${FEELPP_TOOLBOX_COEFFICIENTFORMPDES_REGISTER_ENTRY_CLASS_TYPE}" )
    else()
      set(FEELPP_TOOLBOX_COEFFICIENTFORMPDES_REGISTER_ALL_ENTRIES_CLASS_TYPE "${FEELPP_TOOLBOX_COEFFICIENTFORMPDES_REGISTER_ALL_ENTRIES_CLASS_TYPE} ,\n${FEELPP_TOOLBOX_COEFFICIENTFORMPDES_REGISTER_ENTRY_CLASS_TYPE}" )
    endif()
  endforeach()
endforeach()

#########################################################################
#message( "FEELPP_TOOLBOX_COEFFICIENTFORMPDES_REGISTER_ALL_ENTRIES_CLASS_TYPE=${FEELPP_TOOLBOX_COEFFICIENTFORMPDES_REGISTER_ALL_ENTRIES_CLASS_TYPE}")
set(FEELPP_TOOLBOX_COEFFICIENTFORMPDES_REGISTERD_TYPE_HPP ${CMAKE_CURRENT_BINARY_DIR}/coefficientformpdes_registered_type.hpp )
set(FEELPP_TOOLBOX_COEFFICIENTFORMPDES_REGISTERD_TYPE_HPP_TMP  ${CMAKE_CURRENT_BINARY_DIR}/coefficientformpdes_registered_type.hpp.tmp)
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/coefficientformpdes_registered_type.hpp.in ${FEELPP_TOOLBOX_COEFFICIENTFORMPDES_REGISTERD_TYPE_HPP_TMP} )
if ( NOT EXISTS ${FEELPP_TOOLBOX_COEFFICIENTFORMPDES_REGISTERD_TYPE_HPP} )
  file(WRITE ${FEELPP_TOOLBOX_COEFFICIENTFORMPDES_REGISTERD_TYPE_HPP} "") #write empty file
endif()
# check if the file has changed
execute_process(
  COMMAND ${CMAKE_COMMAND} -E compare_files
  ${FEELPP_TOOLBOX_COEFFICIENTFORMPDES_REGISTERD_TYPE_HPP_TMP}
  ${FEELPP_TOOLBOX_COEFFICIENTFORMPDES_REGISTERD_TYPE_HPP}
  RESULT_VARIABLE FEELPP_TOOLBOX_COEFFICIENTFORMPDES_REGISTERD_TYPE_HPP_NEEDS_UPDATING
  OUTPUT_QUIET ERROR_QUIET
  )
#if change then copy the file
if ( FEELPP_TOOLBOX_COEFFICIENTFORMPDES_REGISTERD_TYPE_HPP_NEEDS_UPDATING )
  execute_process( COMMAND ${CMAKE_COMMAND} -E copy ${FEELPP_TOOLBOX_COEFFICIENTFORMPDES_REGISTERD_TYPE_HPP_TMP} ${FEELPP_TOOLBOX_COEFFICIENTFORMPDES_REGISTERD_TYPE_HPP} )
endif()
INSTALL(FILES ${FEELPP_TOOLBOX_COEFFICIENTFORMPDES_REGISTERD_TYPE_HPP} DESTINATION include/feelpp/toolboxes/feel/feelmodels/coefficientformpdes COMPONENT Devel)
#########################################################################

feelpp_toolboxes_add_library(coefficientformpdes_lib SRCS libcoefficientformpdes.cpp LINK_LIBRARIES ${FEELPP_MODELS_LIBRARIES_COEFFICIENTFORMPDES})

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/coefficientformpdes COMPONENT Devel)
