
# geo order 1
foreach( _DIM 2 3)
  foreach( _DISP_ORDER 1 2 3 4 )
    genLibSolidMechanics(
      DIM ${_DIM}
      DISP_ORDER ${_DISP_ORDER}
      GEO_ORDER 1
      DENSITY_COEFFLAME_TYPE P0d
      )
    set(FEELPP_MODELS_LIBRARIES_SOLIDMECHANICS  ${FEELPP_MODELS_LIBRARIES_SOLIDMECHANICS} ${LIBBASE_NAME}  )
  endforeach()
endforeach()

# iso-parametric
foreach( _DIM 2 3)
  foreach( _DISP_ORDER 2 3 4 )
    genLibSolidMechanics(
      DIM ${_DIM}
      DISP_ORDER ${_DISP_ORDER}
      GEO_ORDER ${_DISP_ORDER}
      DENSITY_COEFFLAME_TYPE P0d
      )
    set(FEELPP_MODELS_LIBRARIES_SOLIDMECHANICS  ${FEELPP_MODELS_LIBRARIES_SOLIDMECHANICS} ${LIBBASE_NAME}  )
  endforeach()
endforeach()

add_library(feelpp_model_solidmechanics SHARED libsolidmechanics.cpp )
target_link_libraries(feelpp_model_solidmechanics ${FEELPP_MODELS_LIBRARIES_SOLIDMECHANICS})
