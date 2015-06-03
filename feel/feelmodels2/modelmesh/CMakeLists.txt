

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

add_library(
  feelpp_modelmesh
  SHARED
  ${HARMONICEXTENSION_SRC}
  ${WINSLOW_SRC}
  ale2d.cpp 
  ale3d.cpp
  )
target_link_libraries(feelpp_modelmesh feelpp_modelalg feelpp_modelcore ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES})

endif( FEELPP_MODELS_HAS_MESHALE )
