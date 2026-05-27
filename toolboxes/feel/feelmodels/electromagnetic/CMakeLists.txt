
set( _MAGNETIC_FE_BASIS_TYPE_LIST "Lagrange<1,Vectorial,Continuous,PointSetFekete>" "Nedelec<0,NedelecKind::NED1>")
set( _MAGNETIC_FE_BASIS_TAG_LIST "Pchv1" "Ned1h0")
set( _ELECTRIC_FE_BASIS_TYPE_LIST "Lagrange<1,Scalar,Continuous,PointSetFekete>")
set( _ELECTRIC_FE_BASIS_TAG_LIST "Pch1")

foreach(_DIM 2 3)
  foreach( _FE_BASIS IN ZIP_LISTS
      _ELECTRIC_FE_BASIS_TYPE_LIST _ELECTRIC_FE_BASIS_TAG_LIST
      _MAGNETIC_FE_BASIS_TYPE_LIST _MAGNETIC_FE_BASIS_TAG_LIST )
    message(STATUS "[feelpp] ELECTROMAGNETIC ${_DIM}D Electric:${_FE_BASIS_1} Magnetic:${_FE_BASIS_3}")
    genLibElectromagnetic(
      DIM     ${_DIM}
      ELECTRIC_P_FE_BASIS_TYPE ${_FE_BASIS_0}
      ELECTRIC_P_FE_BASIS_TAG ${_FE_BASIS_1}
      MAGNETIC_P_FE_BASIS_TYPE ${_FE_BASIS_2}
      MAGNETIC_P_FE_BASIS_TAG ${_FE_BASIS_3}
      GEO_ORDER 1
    )
    set(FEELPP_MODELS_LIBRARIES_ELECTROMAGNETIC ${FEELPP_MODELS_LIBRARIES_ELECTROMAGNETIC} ${ELECTROMAGNETIC_LIB_NAME}  )
  endforeach()
endforeach()

add_library(feelpp_toolbox_electromagnetic_lib SHARED libelectromagnetic.cpp )
add_library(Feelpp::feelpp_toolbox_electromagnetic_lib ALIAS feelpp_toolbox_electromagnetic_lib )  # to match exported target
set_target_properties(feelpp_toolbox_electromagnetic_lib PROPERTIES OUTPUT_NAME "feelpp_toolbox_electromagnetic_lib")
target_link_libraries(feelpp_toolbox_electromagnetic_lib ${FEELPP_MODELS_LIBRARIES_ELECTROMAGNETIC})
set_property(TARGET feelpp_toolbox_electromagnetic_lib PROPERTY MACOSX_RPATH ON)
set_target_properties(feelpp_toolbox_electromagnetic_lib PROPERTIES VERSION ${FEELPP_TOOLBOXES_SHARED_VERSION} SOVERSION ${FEELPP_TOOLBOXES_SHARED_SOVERSION})
install(TARGETS feelpp_toolbox_electromagnetic_lib EXPORT feelpp-toolboxes-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/electromagnetic COMPONENT Devel)
