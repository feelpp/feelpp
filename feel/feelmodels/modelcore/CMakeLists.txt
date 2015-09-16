#############################################################################

option(FEELMODELS_ENABLE_FLUID_MECHANICS "Enable/Disable FluidMechanics" ON)

option(FEELMODELS_FLUID_BUILD_LINEAR_CODE "Enable/Disable" ON)
option(FEELMODELS_FLUID_BUILD_JACOBIAN_CODE "Enable/Disable" ON)
option(FEELMODELS_FLUID_BUILD_RESIDUAL_CODE "Enable/Disable" ON)
option(FEELMODELS_FLUID_BUILD_PTFIXED_CODE "Enable/Disable" OFF)
option(FEELMODELS_FLUID_BUILD_STABILISATION_CODE "Enable/Disable" ON)

option(FEELMODELS_SOLID_BUILD_LINEAR_CODE "Enable/Disable" ON)
option(FEELMODELS_SOLID_BUILD_JACOBIAN_CODE "Enable/Disable" ON)
option(FEELMODELS_SOLID_BUILD_RESIDUAL_CODE "Enable/Disable" ON)

SET(FEELMODELS_FSIBASE_NUMBER_OF_PARAMETERS "10" CACHE STRING "number of parameters" )
SET(FEELMODELS_FSIBASE_NUMBER_OF_GEOPARAMETERS "10" CACHE STRING "number of geometric parameters" )
SET(FEELMODELS_FSIBASE_NUMBER_OF_GINACEXPR "10" CACHE STRING "number of ginac expr" )

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/feelmodelscoreconstconfig.h.in ${CMAKE_CURRENT_BINARY_DIR}/feelmodelscoreconstconfig.h )

set(FEELMODELS_LIBBASE_CHECK_DIR "${CMAKE_CURRENT_BINARY_DIR}/checkmodelsbuilt" CACHE STRING "FEELMODELS_LIBBASE_CHECK_DIR")
file(REMOVE_RECURSE ${FEELMODELS_LIBBASE_CHECK_DIR})


#############################################################################
#############################################################################

add_library(feelpp_modelcore
  SHARED
  options.cpp log.cpp timertool.cpp
  modelbase.cpp modelalgebraic.cpp modelnumerical.cpp
  markermanagement.cpp
  )
target_link_libraries(feelpp_modelcore ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES})
INSTALL(TARGETS feelpp_modelcore DESTINATION lib/ COMPONENT LibsFeelppModelsCommon)

include(feelpp.models.macros.cmake)