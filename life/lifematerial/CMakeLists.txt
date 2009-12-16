add_library(lifematerial
  SHARED
  material.cpp
  materiallib.cpp )
set_target_properties(lifematerial PROPERTIES VERSION ${LIFE_SHARED_VERSION} SOVERSION ${LIFE_SHARED_SOVERSION})
