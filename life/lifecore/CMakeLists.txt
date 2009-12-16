add_library(lifecore SHARED
  info.cpp
  application.cpp
  xmlparser.hpp
  applicationxml.cpp
  about.cpp
  debug.cpp
  typeinfo.cpp
  smartassert.cpp
  pslogger.cpp )
set_target_properties( lifecore PROPERTIES VERSION ${LIFE_SHARED_VERSION} SOVERSION ${LIFE_SHARED_SOVERSION})


