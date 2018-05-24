# Fluid: PN+1PNG1, levelset : PKG1
foreach( _DIM 2 3)
  feelpp_max( MAX_ORDER 2 ${FEELPP_INSTANTIATION_ORDER_MAX})
  message(STATUS "[feelpp] MULTIFLUID ${_DIM}D PN+1PNG1_PKG1 N=1..${MAX_ORDER} K=1..${MAX_ORDER}")
  foreach( _FLUID_U_ORDER RANGE 2 ${MAX_ORDER} )
    math(EXPR _FLUID_P_ORDER '${_FLUID_U_ORDER}-1')
    foreach( _LS_ORDER RANGE 1 ${MAX_ORDER} )
        genLibMultifluid(
            DIM ${_DIM}
            FLUID_U_ORDER ${_FLUID_U_ORDER}
            FLUID_P_ORDER ${_FLUID_P_ORDER}
            FLUID_P_CONTINUITY Continuous
            FLUID_DENSITY_VISCOSITY_ORDER 1
            FLUID_DENSITY_VISCOSITY_CONTINUITY Continuous
            LEVELSET_ORDER ${_LS_ORDER}
            LEVELSET_PN_ORDER 2
            GEO_ORDER  1
            )
        set(FEELPP_MODELS_LIBRARIES_MULTIFLUID  ${FEELPP_MODELS_LIBRARIES_MULTIFLUID} ${LIBBASE_NAME}  )
    endforeach()
  endforeach()
endforeach()

add_library(feelpp_model_multifluid SHARED libmultifluid.cpp)
target_link_libraries(feelpp_model_multifluid ${FEELPP_MODELS_LIBRARIES_MULTIFLUID})

install(TARGETS feelpp_model_multifluid DESTINATION lib/ COMPONENT Libs)

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/multifluid COMPONENT Devel)
