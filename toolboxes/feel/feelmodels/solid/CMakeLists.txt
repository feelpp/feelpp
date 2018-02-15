
# geo order 1
foreach( _DIM 2 3)
  math(EXPR _D_ORDER '${FEELPP_INSTANTIATION_ORDER_MAX}-1')
  feelpp_max( MAX_ORDER 1 ${_D_ORDER})
  message(STATUS "[feelpp] SOLID ${_DIM}D PNG1 N=1..${MAX_ORDER}")
  foreach( _DISP_ORDER RANGE 1 ${MAX_ORDER} )
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
  math(EXPR _D_ORDER '${FEELPP_INSTANTIATION_ORDER_MAX}-1')
  feelpp_min( MAX_ORDER ${FEELPP_MESH_MAX_ORDER} ${_D_ORDER})
  message(STATUS "[feelpp] SOLID ${_DIM}D PNGN N=2..${MAX_ORDER}")
  foreach( _DISP_ORDER RANGE 2 ${MAX_ORDER} )
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
set_property(TARGET feelpp_model_solidmechanics PROPERTY MACOSX_RPATH ON)
install(TARGETS feelpp_model_solidmechanics DESTINATION lib/ COMPONENT Libs)

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/solid COMPONENT Devel)
