
feelpp_toolboxes_add_library( electricbase SRCS electricboundaryconditions.cpp LINK_LIBRARIES feelpp_modelcore )

# geo order 1
foreach(_DIM 2 3)
  math(EXPR _T_ORDER "${FEELPP_INSTANTIATION_ORDER_MAX}-1")
  feelpp_max( MAX_ORDER 1 ${_T_ORDER})
  foreach( _T_ORDER RANGE 1 ${MAX_ORDER} )
    set(_FE_BASIS_TAG "Pch${_T_ORDER}")
    set(_FE_BASIS_TYPE "Lagrange<${_T_ORDER},Scalar,Continuous,PointSetFekete>")
    message(STATUS "[feelpp] ELECTRIC ${_DIM}DG1 ${_FE_BASIS_TAG}")
    genLibElectric(
      DIM     ${_DIM}
      P_FE_BASIS_TYPE ${_FE_BASIS_TYPE}
      P_FE_BASIS_TAG ${_FE_BASIS_TAG}
      GEO_ORDER 1
    )
    set(FEELPP_MODELS_LIBRARIES_ELECTRIC  ${FEELPP_MODELS_LIBRARIES_ELECTRIC} ${ELECTRIC_LIB_NAME}  )

    if ( ${_T_ORDER} GREATER_EQUAL 2 AND ${FEELPP_MESH_MAX_ORDER} GREATER_EQUAL ${_T_ORDER} )
      message(STATUS "[feelpp] ELECTRIC ${_DIM}DG${_T_ORDER} ${_FE_BASIS_TAG}")
      genLibElectric(
        DIM     ${_DIM}
        P_FE_BASIS_TYPE ${_FE_BASIS_TYPE}
        P_FE_BASIS_TAG ${_FE_BASIS_TAG}
        GEO_ORDER ${_T_ORDER}
      )
      set(FEELPP_MODELS_LIBRARIES_ELECTRIC  ${FEELPP_MODELS_LIBRARIES_ELECTRIC} ${ELECTRIC_LIB_NAME}  )
    endif()


    # genLibElectric(
    #   DIM     ${_DIM}
    #   P_ORDER ${_T_ORDER}
    #   GEO_ORDER 1
    #   )
  endforeach()
endforeach()

# # isoparametric
# foreach(_DIM 2 3)
#   math(EXPR _T_ORDER "${FEELPP_INSTANTIATION_ORDER_MAX}-1")
#   feelpp_min( MAX_ORDER ${FEELPP_MESH_MAX_ORDER} ${_T_ORDER})

#   # check that start order is greater than stop order for RANGE below otherwise cmake may have an undefined behavior
#   if(${MAX_ORDER} GREATER_EQUAL 2)
#     message(STATUS "[feelpp] ELECTRIC ${_DIM}D PNGN N=2..${MAX_ORDER}")
#     foreach( _T_ORDER RANGE 2 ${MAX_ORDER} )
#       genLibElectric(
#         DIM     ${_DIM}
#         P_ORDER ${_T_ORDER}
#         GEO_ORDER ${_T_ORDER}
#         )
#       set(FEELPP_MODELS_LIBRARIES_ELECTRIC  ${FEELPP_MODELS_LIBRARIES_ELECTRIC} ${ELECTRIC_LIB_NAME} )
#     endforeach()
#   endif()
# endforeach()

feelpp_toolboxes_add_library(electric_lib SRCS libelectric.cpp LINK_LIBRARIES ${FEELPP_MODELS_LIBRARIES_ELECTRIC})

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/electric COMPONENT Devel)
