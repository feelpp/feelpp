
feelpp_toolboxes_add_library( solidbase SRCS
  solidmechanicsboundaryconditions.cpp
  solidmechanics1dreducedboundaryconditions.cpp
  LINK_LIBRARIES feelpp_modelcore )

# geo order 1
foreach( _DIM 2 3)
  math(EXPR _D_ORDER "${FEELPP_INSTANTIATION_ORDER_MAX}-1")
  feelpp_max( MAX_ORDER 1 ${_D_ORDER})
  message(STATUS "[feelpp] SOLID ${_DIM}D PNG1 N=1..${MAX_ORDER}")
  foreach( _DISP_ORDER RANGE 1 ${MAX_ORDER} )
    genLibSolidMechanics(
      DIM ${_DIM}
      DISP_ORDER ${_DISP_ORDER}
      GEO_ORDER 1
      )
    set(FEELPP_MODELS_LIBRARIES_SOLIDMECHANICS  ${FEELPP_MODELS_LIBRARIES_SOLIDMECHANICS} ${SOLIDMECHANICS_LIB_NAME}  )
  endforeach()
endforeach()

# iso-parametric
foreach( _DIM 2 3)
  math(EXPR _D_ORDER "${FEELPP_INSTANTIATION_ORDER_MAX}-1")
  feelpp_min( MAX_ORDER ${FEELPP_MESH_MAX_ORDER} ${_D_ORDER})
  message(STATUS "[feelpp] SOLID ${_DIM}D PNGN N=2..${MAX_ORDER}")
  foreach( _DISP_ORDER RANGE 2 ${MAX_ORDER} )
    genLibSolidMechanics(
      DIM ${_DIM}
      DISP_ORDER ${_DISP_ORDER}
      GEO_ORDER ${_DISP_ORDER}
      )
    set(FEELPP_MODELS_LIBRARIES_SOLIDMECHANICS  ${FEELPP_MODELS_LIBRARIES_SOLIDMECHANICS} ${SOLIDMECHANICS_LIB_NAME}  )
  endforeach()
endforeach()

feelpp_toolboxes_add_library(solid_lib SRCS libsolidmechanics.cpp LINK_LIBRARIES ${FEELPP_MODELS_LIBRARIES_SOLIDMECHANICS})

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/solid COMPONENT Devel)
