
set(FEELPP_MODELS_MESH_SRC fsimesh.cpp)

if( FEELPP_MODELS_HAS_MESHALE )

if( FEELPP_MODELS_HAS_MESHALE_HARMONICEXTENSION )
  set(HARMONICEXTENSION_SRC
    harmonicextension.cpp
    )
endif()


if( FEELPP_MODELS_HAS_MESHALE_WINSLOW )
  set(WINSLOW_SRC 
    winslowgen2d.cpp
    winslowgen3d.cpp
    winslowptfixe.cpp
    winslownewton.cpp
    )
endif()


set(FEELPP_MODELS_MESH_SRC
  ${FEELPP_MODELS_MESH_SRC}
  ${HARMONICEXTENSION_SRC}
  ${WINSLOW_SRC}
  ale2d.cpp 
  ale3d.cpp
  )
set( FEELPP_MODELS_MESH_LIB_DEPENDS feelpp_modelalg feelpp_modelcore )
# add_library(
#   feelpp_modelmesh
#   SHARED
#   ${HARMONICEXTENSION_SRC}
#   ${WINSLOW_SRC}
#   ale2d.cpp 
#   ale3d.cpp
#   )
# target_link_libraries(feelpp_modelmesh feelpp_modelalg feelpp_modelcore ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES})

endif( FEELPP_MODELS_HAS_MESHALE )

add_library(
  feelpp_modelmesh
  SHARED
  ${FEELPP_MODELS_MESH_SRC}
)
target_link_libraries(feelpp_modelmesh ${FEELPP_MODELS_MESH_LIB_DEPENDS} ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES})
