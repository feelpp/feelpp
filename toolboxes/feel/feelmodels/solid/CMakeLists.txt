
feelpp_toolboxes_add_library( solidbase SRCS
  solidmechanicsboundaryconditions.cpp
  solidmechanics1dreducedboundaryconditions.cpp
  LINK_LIBRARIES feelpp_modelcore )

# geo order 1
foreach( _DIM 2 3)
  foreach(_DISP_ORDER 1 2 Dynamic)
    foreach( _CONVEX_TYPE Simplex Hypercube )
      message(STATUS "[feelpp] SOLID ${_DIM}D P${_ORDER}G1 ${_CONVEX_TYPE}")
      genLibSolidMechanics(
        DIM ${_DIM}
        DISP_ORDER ${_DISP_ORDER}
        GEO_ORDER 1
        CONVEX ${_CONVEX_TYPE}
      )
      set(FEELPP_MODELS_LIBRARIES_SOLIDMECHANICS  ${FEELPP_MODELS_LIBRARIES_SOLIDMECHANICS} ${SOLIDMECHANICS_LIB_NAME}  )
    endforeach()
  endforeach()
endforeach()

# geo and disp order Dynamic
foreach( _DIM 2 3)
  foreach( _CONVEX_TYPE Simplex Hypercube )
    message(STATUS "[feelpp] SOLID ${_DIM}D PDynamicGDynamic ${_CONVEX_TYPE}")
    genLibSolidMechanics(
      DIM ${_DIM}
      DISP_ORDER Dynamic
      GEO_ORDER Dynamic
      CONVEX ${_CONVEX_TYPE}
    )
    set(FEELPP_MODELS_LIBRARIES_SOLIDMECHANICS  ${FEELPP_MODELS_LIBRARIES_SOLIDMECHANICS} ${SOLIDMECHANICS_LIB_NAME}  )
  endforeach()
endforeach()

feelpp_toolboxes_add_library(solid_lib SRCS libsolidmechanics.cpp LINK_LIBRARIES ${FEELPP_MODELS_LIBRARIES_SOLIDMECHANICS})

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/solid COMPONENT Devel)
