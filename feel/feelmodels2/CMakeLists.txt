
SET(FEELPP_MODELS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
SET(FEELPP_MODELS_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE INTERNAL "")

option(FEELMODELS_ENABLE_MESHALE "Enable/Disable mesh ALE" ON)
if ( FEELMODELS_ENABLE_MESHALE )
  set( FEELPP_MODELS_HAS_MESHALE 1)
endif()
set(FEELPP_MODELS_HAS_MESHALE_HARMONICEXTENSION 1)
#set(FEELPP_MODELS_HAS_MESHALE_WINSLOW 1)


add_subdirectory(modelcore)
add_subdirectory(modelalg)
add_subdirectory(modelmesh)


# generate models lib
add_subdirectory(thermodyn)
add_subdirectory(solid)
add_subdirectory(fluid)
add_subdirectory(fsi)
