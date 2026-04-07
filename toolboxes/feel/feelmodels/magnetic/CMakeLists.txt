
feelpp_toolboxes_add_library( magneticbase SRCS magneticboundaryconditions.cpp LINK_LIBRARIES feelpp_modelcore )

# set( _UNKNOWN_BASIS_TYPE_LIST
#   #"Lagrange<0,Scalar,Continuous,PointSetFekete>"
#   "Lagrange<1,Scalar,Continuous,PointSetFekete>" "Lagrange<2,Scalar,Continuous,PointSetFekete>"
#   "Lagrange<1,Vectorial,Continuous,PointSetFekete>" "Lagrange<2,Vectorial,Continuous,PointSetFekete>"
#   "Lagrange<1,Scalar,Discontinuous,PointSetFekete>"
#   "Nedelec<0,NedelecKind::NED1>"
# )
# set( _UNKNOWN_BASIS_TAG_LIST
#   #"Pch0"
#   "Pch1" "Pch2"
#   "Pchv1" "Pchv2"
#   "Pdh1"
#   "Ned1h0"
# )

set( _FE_BASIS_TYPE_LIST "Nedelec<0,NedelecKind::NED1>")
set( _FE_BASIS_TAG_LIST "Ned1h0")
foreach( _JJ IN ZIP_LISTS _FE_BASIS_TYPE_LIST _FE_BASIS_TAG_LIST)
  message("JJ=${_JJ_0} and ${_JJ_1}")
endforeach()

# geo order 1
foreach(_DIM 2 3)
  foreach( _FE_BASIS IN ZIP_LISTS _FE_BASIS_TYPE_LIST _FE_BASIS_TAG_LIST)

    #math(EXPR _T_ORDER "${FEELPP_INSTANTIATION_ORDER_MAX}")
    #feelpp_max( MAX_ORDER 1 ${_T_ORDER})
    message(STATUS "[feelpp] MAGNETIC ${_DIM}D ${_FE_BASIS_1}")
    #foreach( _T_ORDER RANGE 1 ${MAX_ORDER} )
    genLibMagnetic(
      DIM     ${_DIM}
      P_FE_BASIS_TYPE ${_FE_BASIS_0}
      P_FE_BASIS_TAG ${_FE_BASIS_1}
      GEO_ORDER 1
    )
    set(FEELPP_MODELS_LIBRARIES_MAGNETIC  ${FEELPP_MODELS_LIBRARIES_MAGNETIC} ${MAGNETIC_LIB_NAME}  )
  endforeach()
endforeach()

feelpp_toolboxes_add_library(magnetic_lib SRCS libmagnetic.cpp LINK_LIBRARIES ${FEELPP_MODELS_LIBRARIES_MAGNETIC})

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/magnetic COMPONENT Devel)
